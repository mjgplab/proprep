"""Mean-field self-consistent state-map iteration (multi-state).

For each titratable site, run PB on every protonation state (with all
other sites held at their current state-map values), build the
chemistry-group partition function at the system pH, and assign the
site to the dominant chemistry group's most-stable representative
state. Repeat until no states change, or max_iter reached.

This is a generalization of the bundle's `lifted/iterate.py` logic to
arbitrary residue mixtures: AS4/GL4 with proper tautomer averaging,
HIS as a true 3-state titration (HIP/HID/HIE), LYS/TYR/CYS in the
2-state-equivalent picture.

Use this when:
  - Sites are weakly coupled (max |W_ij| < kT)
  - Or coupling is sparse (most pairs > 15 Å apart)
  - Speed matters more than capturing rare correlated flips

Cost per iteration: N_sites × N_states_per_site × 2 (cluster+model) PB
calls. Scales as O(N) for typical residue mixes.

For strongly-coupled systems, prefer `solvers.monte_carlo`.
"""
from __future__ import annotations

import csv
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from ..intrinsic import (
    BackgroundStateMap, SiteKey,
    compute_state_energies,
    populations_by_chemistry_group,
    state_free_energies_at_pH,
    group_states_by_chemistry,
)
from ..residues import TitratableResidue, get_residue, ProtonationState


# ---------------------------------------------------------------------------
# Site abstraction
# ---------------------------------------------------------------------------

@dataclass
class TitratableSite:
    """One titratable site under consideration.

    `state` carries the *currently assigned* ProtonationState — used as
    the background charge distribution for other sites' PB calls in the
    next iteration. For multi-state residues this is the most-stable
    state within the dominant chemistry group at the last iteration.
    """
    resname: str
    resnum:  int                       # 1-based
    residue: TitratableResidue
    state:   ProtonationState

    @property
    def key(self) -> SiteKey:
        return (self.resname, self.resnum)

    def chemistry_label(self) -> str:
        """Short human-readable label for the current chemistry."""
        return self.residue.state_label(self.state)


def make_site(resname: str, resnum: int,
               initial: str = "neutral",
               pH: float = 7.0) -> TitratableSite:
    """Construct a TitratableSite.

    `initial` ∈ {'PROT', 'DEPROT', 'neutral', 'auto'}:
        * 'neutral' (default) — start at the residue's neutral state, so
          the very first iteration's background matches the all-neutral
          reference used by intrinsic.compute_pka_multistate and by the
          coupling layer's ref_chem. This keeps every layer of the
          pipeline (single-site, coupling, mean-field) on the same
          reference.
        * 'PROT' / 'DEPROT' — explicit fully-protonated / fully-
          deprotonated initial state.
        * 'auto' — pH-derived initial state (PROT if pKa_model > pH,
          else DEPROT). Legacy behaviour; equivalent to running constph
          with the model pKa as the only input.
    """
    res = get_residue(resname)
    if initial == "auto":
        initial = "PROT" if res.pka_model > pH else "DEPROT"
    if initial == "neutral":
        state = res.neutral_state
    elif initial == "PROT":
        state = res.prot_state
    else:
        state = res.deprot_state
    return TitratableSite(resname=resname, resnum=resnum,
                           residue=res, state=state)


# ---------------------------------------------------------------------------
# Worker (one site per call)
# ---------------------------------------------------------------------------

_W_CLUSTER_PRM:   Optional[Path] = None
_W_CLUSTER_RST:   Optional[Path] = None
_W_PBSA_PARAMS:   Optional[Dict[str, Any]] = None
_W_WORK_ROOT:     Optional[Path] = None
_W_PH:            float = 7.0
# Cluster cutout + fixed-charge background, mirroring coupling._self_init.
# Without these the worker ran PB on the FULL structure (NATYP huge ->
# "rdparm: a parameter array overflowed"); the cutout keeps each PB call to a
# ~single-site cluster exactly as the pbt-3 single-site path does.
_W_CLUSTER_RADIUS: Optional[float] = None
_W_BOND_DEPTH:     int = 3
_W_SITE_ENVELOPES: Dict[Tuple[str, int],
                        Tuple[List[int], List[List[int]]]] = {}
_W_EXTRA_BACKGROUND: Dict[SiteKey, Any] = {}
_W_MODEL_CACHE_DIR: Optional[Path] = None


def _worker_init(cluster_prm: str, cluster_rst: str,
                  pbsa_params: Dict[str, Any], work_root: str,
                  pH: float, cluster_radius: Optional[float] = None,
                  bond_depth: int = 3,
                  site_envelopes_pkl: Optional[bytes] = None,
                  extra_background_pkl: Optional[bytes] = None,
                  model_cache_dir: Optional[str] = None):
    import pickle
    import warnings
    from parmed.exceptions import AmberWarning
    # Suppress parmed's contiguity warning in worker processes — see
    # coupling._self_init for rationale.
    warnings.filterwarnings(
        "ignore",
        message=r".*Molecule atoms are not contiguous.*",
        category=AmberWarning)
    global _W_CLUSTER_PRM, _W_CLUSTER_RST, _W_PBSA_PARAMS, _W_WORK_ROOT, _W_PH
    global _W_CLUSTER_RADIUS, _W_BOND_DEPTH, _W_SITE_ENVELOPES
    global _W_EXTRA_BACKGROUND, _W_MODEL_CACHE_DIR
    _W_CLUSTER_PRM = Path(cluster_prm)
    _W_CLUSTER_RST = Path(cluster_rst)
    _W_PBSA_PARAMS = pbsa_params
    _W_WORK_ROOT   = Path(work_root)
    _W_PH          = pH
    _W_CLUSTER_RADIUS = cluster_radius
    _W_BOND_DEPTH = bond_depth
    _W_SITE_ENVELOPES = (pickle.loads(site_envelopes_pkl)
                          if site_envelopes_pkl else {})
    _W_EXTRA_BACKGROUND = (pickle.loads(extra_background_pkl)
                            if extra_background_pkl else {})
    _W_MODEL_CACHE_DIR = Path(model_cache_dir) if model_cache_dir else None


def _worker_run(args: Tuple[str, int, BackgroundStateMap, str]):
    resname, resnum, background, iter_tag = args
    work_dir = _W_WORK_ROOT / iter_tag / f"{resname}{resnum:04d}"
    try:
        residue = get_residue(resname)
        # Reconstruct the Site against a freshly loaded structure if envelope
        # info was provided (envelope-aware cutout for multi-residue cofactors;
        # mirrors coupling._self_task).
        site = None
        if (resname, resnum) in _W_SITE_ENVELOPES:
            import parmed as pmd
            from ..sites import Site
            env_idxs, other_envs = _W_SITE_ENVELOPES[(resname, resnum)]
            s = pmd.load_file(str(_W_CLUSTER_PRM), str(_W_CLUSTER_RST))
            envelope = [s.residues[i] for i in env_idxs]
            target = next(r for r in envelope
                          if r.name == resname and r.number + 1 == resnum)
            site = Site(titrating_residue=target,
                        envelope_residues=envelope,
                        other_redox_envelopes=other_envs)
        # Fixed-charge background (metal-coordinating residues etc.) merged on
        # top of the per-iteration mean-field background.
        full_background = dict(background)
        if _W_EXTRA_BACKGROUND:
            full_background.update(_W_EXTRA_BACKGROUND)
        se = compute_state_energies(
            cluster_prmtop=_W_CLUSTER_PRM,
            cluster_rst7=  _W_CLUSTER_RST,
            target_resname=resname,
            target_resnum_1based=resnum,
            residue=residue,
            model_cache_dir=_W_MODEL_CACHE_DIR,
            work_dir=work_dir,
            background=full_background,
            cluster_radius=_W_CLUSTER_RADIUS,
            bond_depth=_W_BOND_DEPTH,
            site=site,
            **_W_PBSA_PARAMS)
        chem_pops = populations_by_chemistry_group(se, _W_PH)
        # Most-populated chemistry group
        dominant_chem = max(chem_pops, key=chem_pops.get)
        # Within that group, pick the lowest-G representative state
        groups = group_states_by_chemistry(se)
        member_idxs = groups[dominant_chem]
        G_per_state = state_free_energies_at_pH(se, _W_PH)
        best_state_idx = min(member_idxs, key=lambda i: G_per_state[i])
        return {
            "resname": resname, "resnum": resnum,
            "dominant_chem": dominant_chem,
            "best_state_idx": best_state_idx,
            "pop_by_chem": chem_pops,
            "error": "",
        }
    except Exception as exc:
        return {"resname": resname, "resnum": resnum, "error": str(exc)}


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def solve(sites: List[TitratableSite],
          cluster_prmtop: Path,
          cluster_rst7:   Path,
          *,
          pH: float = 7.0,
          max_iter: int = 5,
          workers:  int = 4,
          work_root: Path = Path("mean_field_work"),
          pbsa_params: Optional[Dict[str, Any]] = None,
          cluster_radius: Optional[float] = None,
          bond_depth: int = 3,
          site_envelopes: Optional[Dict[Tuple[str, int], Any]] = None,
          extra_background: Optional[Dict[SiteKey, Any]] = None,
          model_cache_dir: Optional[Path] = None,
          flipper_window: int = 3,
          ) -> Dict[str, Any]:
    """Run mean-field iteration. Returns final state map + per-iter results.

    State assignment per iteration:
      1. Run PB for every state of every site (in parallel).
      2. For each site, compute populations by chemistry group at pH.
      3. Pick the dominant chemistry; pick the lowest-G state within it.
      4. If the chosen state changed for any site, iterate again.
    """
    import pickle
    pbsa_params = pbsa_params or {}
    work_root = Path(work_root)
    work_root.mkdir(parents=True, exist_ok=True)

    # Route non-fatal pbsa exit-time-crash salvages to a per-run log (workers
    # inherit this env var) instead of one RuntimeWarning per PB call, then
    # tally it per iteration for a single summary line. See pb_backend.
    salvage_log = work_root / "pbsa_salvage.log"
    os.environ["PROPREP_PBSA_SALVAGE_LOG"] = str(salvage_log)

    def _salvage_count() -> int:
        try:
            with open(salvage_log) as f:
                return sum(1 for _ in f)
        except OSError:
            return 0

    # Distill site_envelopes to a worker-friendly form (parmed residues don't
    # pickle; their indices are stable across reloads of the same prmtop).
    # Mirrors coupling.compute_self_energies.
    envelopes_for_workers: Dict[Tuple[str, int],
                                Tuple[List[int], List[List[int]]]] = {}
    if site_envelopes:
        for key, site_obj in site_envelopes.items():
            envelopes_for_workers[key] = (
                sorted(site_obj.envelope_idxs),
                [list(env) for env in site_obj.other_redox_envelopes],
            )
    envelopes_pkl = (pickle.dumps(envelopes_for_workers)
                      if envelopes_for_workers else b"")
    extra_bg_pkl = pickle.dumps(extra_background) if extra_background else b""
    model_cache_arg = str(model_cache_dir) if model_cache_dir else None

    by_key = {s.key: s for s in sites}
    history: List[List[Dict]] = []
    converged = False
    last_n_changes = 0
    last_changed: List[Tuple[str, int]] = []
    # Per-sweep changed-site lists. A frustrated subset oscillates with a
    # period > 1 (period-2 is common), so any SINGLE sweep undercounts the
    # true flipper set. The union over the last `flipper_window` sweeps
    # captures every site participating in the oscillation.
    changed_history: List[List[Tuple[str, int]]] = []

    print(f"=== mean_field.solve (multi-state) ===")
    print(f"  {len(sites)} sites, pH={pH}, max_iter={max_iter}, "
          f"workers={workers}")
    chem_summary = _summarize_chemistry(sites)
    print(f"  Initial chemistry: {chem_summary}")

    for it in range(1, max_iter + 1):
        t0 = time.time()
        background: BackgroundStateMap = {s.key: s.state for s in sites}
        iter_tag = f"iter{it}"
        print(f"\n  --- iteration {it} ---")

        # Baseline salvage-log line count before this sweep's PB calls, so
        # n_salvaged below reports per-iteration (not cumulative) salvages.
        salvage_before = _salvage_count()
        items = [(s.resname, s.resnum, background, iter_tag) for s in sites]
        results: List[Dict] = []
        with ProcessPoolExecutor(
                max_workers=workers,
                initializer=_worker_init,
                initargs=(str(cluster_prmtop), str(cluster_rst7),
                          pbsa_params, str(work_root), pH,
                          cluster_radius, bond_depth, envelopes_pkl,
                          extra_bg_pkl, model_cache_arg)) as pool:
            futures = {pool.submit(_worker_run, item): item for item in items}
            for fut in as_completed(futures):
                results.append(fut.result())

        n_changes = 0
        n_errors  = 0
        first_error = None
        changed_sites: List[Tuple[str, int]] = []
        for r in sorted(results, key=lambda x: (x["resname"], x["resnum"])):
            if r.get("error"):
                n_errors += 1
                if first_error is None:
                    first_error = r["error"]
                continue
            site = by_key[(r["resname"], r["resnum"])]
            new_state = site.residue.states[r["best_state_idx"]]
            old_chem  = (site.state.prot_count, site.state.pka_corr)
            new_chem  = r["dominant_chem"]
            if old_chem != new_chem:
                n_changes += 1
                changed_sites.append((r["resname"], r["resnum"]))
            site.state = new_state

        last_n_changes = n_changes
        last_changed = changed_sites
        changed_history.append(changed_sites)
        history.append(results)
        elapsed = time.time() - t0
        n_salvaged = _salvage_count() - salvage_before
        salvage_note = (f", pbsa exit-crashes salvaged: {n_salvaged}"
                        if n_salvaged > 0 else "")
        print(f"    {len(results)} sites done in {elapsed:.0f}s. "
              f"Changes: {n_changes}, errors: {n_errors}{salvage_note}")

        # A sweep where every site errored leaves n_changes == 0, which would
        # otherwise be misread as convergence — silently returning the
        # all-neutral initial guess. Treat a fully-failed sweep as a hard error
        # and surface the first failure (commonly the cpinutil/parmed PATH-
        # shadow that only bites cold-cache fork workers).
        if n_errors and n_errors == len(results):
            raise RuntimeError(
                f"Mean-field sweep {it}: all {n_errors} site(s) "
                f"failed — not converged. First error:\n{first_error}")
        if n_errors:
            print(f"    [warning] {n_errors} site(s) errored this sweep; "
                  f"first: {str(first_error)[:200]}")

        if n_changes == 0:
            converged = True
            print(f"    CONVERGED.")
            break
    else:
        print(f"  reached max_iter without convergence")

    chem_summary = _summarize_chemistry(sites)
    print(f"\n  Final chemistry: {chem_summary}")
    return {
        "sites": sites,
        "history": history,
        "converged": converged,
        "iterations": len(history),
        # n_changes / last-sweep changed sites — for the caller's "not
        # converged; add iterations?" prompt.
        "last_n_changes": (0 if converged else last_n_changes),
        "last_changed_sites": ([] if converged else last_changed),
        # Union of changed sites over the last `flipper_window` sweeps. A
        # frustrated subset oscillates (often period-2), so a single sweep
        # undercounts; this is the set targeted coupling should seed from.
        "flipper_sites": ([] if converged else
                          _union_recent(changed_history, flipper_window)),
    }


def _union_recent(changed_history: List[List[Tuple[str, int]]],
                  window: int) -> List[Tuple[str, int]]:
    """Union of changed-site lists over the last `window` sweeps, ordered by
    residue number then name for stable output. `window <= 0` means all sweeps."""
    recent = changed_history if window <= 0 else changed_history[-window:]
    seen = set()
    for sweep in recent:
        seen.update(sweep)
    return sorted(seen, key=lambda k: (k[1], k[0]))


def _summarize_chemistry(sites: List[TitratableSite]) -> str:
    counts: Dict[str, int] = {}
    for s in sites:
        label = f"{s.resname}.{s.chemistry_label()}"
        counts[label] = counts.get(label, 0) + 1
    return ", ".join(f"{k}={v}" for k, v in sorted(counts.items()))


def write_state_map(sites: List[TitratableSite], path: Path) -> None:
    """Write a state-map CSV (resname, resnum, chemistry, state_name)."""
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["resname", "resnum", "chemistry", "state_name",
                     "prot_count", "pka_corr", "pKa_model"])
        for s in sorted(sites, key=lambda x: (x.resname, x.resnum)):
            w.writerow([s.resname, s.resnum, s.chemistry_label(),
                         s.state.name, s.state.prot_count, s.state.pka_corr,
                         s.residue.pka_model])
