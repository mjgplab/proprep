"""Pairwise coupling matrix W_ij(s,t) for the precompute → solver split.

For an N-site titratable system, the total energy under linearized PB is

    G({s_i}, pH) = Σ_i G_intr_i(s_i, pH) + Σ_{i<j} W_ij(s_i, s_j)

where:
  - G_intr_i(s_i, pH) is the intrinsic energy of site i in chemistry s_i
    at the system pH, with all OTHER sites at their reference state.
    This is what `intrinsic.compute_state_energies` already computes.
  - W_ij(s_i, s_j) is the pairwise coupling — how much site j's energy
    changes when site i flips from its reference state to s_i (and j to t_j).

Once (G_intr, W_ij) is precomputed, every solver (mean-field, MC, exact
enumeration) becomes cheap integer arithmetic on these matrices.

Implementation cost (O(N²)):
  - N PB calls for self-energies (one per site, all states; folded into
    intrinsic.compute_state_energies)
  - O(N(N-1)/2 × G²) PB calls for pair couplings, where G is the number
    of chemistry groups per residue (typically 2-3)

Future optimization (O(N) via potential maps): use pbsa's `outphi`
output to extract φ_i(r) maps, then W_ij = Σ over atoms of j of
δq_j · φ_i(r_atom). Deferred — O(N²) is fine for systems up to ~100
sites; the bundle's 122 PRNs would benefit.
"""
from __future__ import annotations

import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import parmed as pmd

from .intrinsic import (
    SiteKey, apply_state_charges,
    compute_state_energies, group_states_by_chemistry,
    state_free_energies_at_pH,
)
from .residues import ProtonationState, TitratableResidue, get_residue
from .pb_backend import RT_LN10, run_pbsa, total_electrostatic
from .cluster import (
    cluster_for_pair, cluster_for_site_pair,
    slice_structure, recommended_pair_cutoff,
)
from .sites import Site


ChemKey = Tuple[int, float]  # (prot_count, pka_corr) — same as intrinsic.py


# ---------------------------------------------------------------------------
# Site descriptor
# ---------------------------------------------------------------------------

@dataclass
class CouplingSite:
    """One titratable site for the coupling/MC layer."""
    resname: str
    resnum:  int                       # 1-based
    residue: TitratableResidue
    chemistries: List[ChemKey]         # ordered list of (n, pka_corr)
    chem_to_state: Dict[ChemKey, ProtonationState]  # representative state
    ref_chem:    ChemKey               # reference chemistry (neutral form)

    @property
    def key(self) -> SiteKey:
        return (self.resname, self.resnum)

    def n_chemistries(self) -> int:
        return len(self.chemistries)


def make_coupling_site(resname: str, resnum: int) -> CouplingSite:
    """Build a CouplingSite from a residue name and number.

    The reference chemistry is set to the residue's neutral form (the
    chemistry whose representative state has charge closest to zero).
    This matches the all-neutral background convention used by the
    single-site step (intrinsic.compute_pka_multistate); keeping the
    two layers on the same reference is what makes the W matrix
    composable with the intrinsic ΔΔGs in MC / exact-enumeration.
    """
    res = get_residue(resname)
    # Group residue states by chemistry; pick first member of each as rep.
    chem_to_state: Dict[ChemKey, ProtonationState] = {}
    for s in res.states:
        key: ChemKey = (s.prot_count, s.pka_corr)
        if key not in chem_to_state:
            chem_to_state[key] = s
    chemistries = list(chem_to_state.keys())
    # Reference = the chemistry whose representative state is the residue's
    # neutral form. This keeps coupling consistent with single-site (which
    # holds non-target titratables at residue.neutral_state).
    neutral = res.neutral_state
    ref_chem: ChemKey = (neutral.prot_count, neutral.pka_corr)
    return CouplingSite(resname=resname, resnum=resnum, residue=res,
                          chemistries=chemistries,
                          chem_to_state=chem_to_state,
                          ref_chem=ref_chem)


# ---------------------------------------------------------------------------
# Self-energies (intrinsic ΔΔG per chemistry per site)
# ---------------------------------------------------------------------------

@dataclass
class SelfEnergies:
    """Intrinsic ΔΔG_PB and pH-dependent terms per (site, chemistry).

    `intrinsic_ddG[(site_key, chem)]` is (E_PB_protein_s − E_PB_model_s)
    for the chemistry-group representative state s. Since within a
    chemistry group all members share (n, pka_corr), the pH term in
    G_eff is (n · pH − pka_corr) — the chemistry key itself.

    The chemistry-group-averaged intrinsic energy (with degeneracy
    correction) is also stored in `chem_group_ddG`.
    """
    sites: List[CouplingSite]
    intrinsic_ddG: Dict[Tuple[SiteKey, ChemKey], float] = field(default_factory=dict)
    state_energies_full: Dict[SiteKey, Dict[int, Dict[str, Any]]] = field(default_factory=dict)


# Globals used by ProcessPoolExecutor workers. Set in `_self_init`.
_W_PRM: Optional[Path] = None
_W_RST: Optional[Path] = None
_W_BACKGROUND: Optional[Dict[SiteKey, ProtonationState]] = None
_W_WORK_ROOT: Optional[Path] = None
_W_PBSA_PARAMS: Optional[Dict[str, Any]] = None
# Writable cache for on-demand model compounds (never the read-only package).
_W_MODEL_CACHE_DIR: Optional[Path] = None


# Globals for self-energy worker
_W_CLUSTER_RADIUS: Optional[float] = None
_W_BOND_DEPTH: int = 3
# Map of (resname, resnum) -> (envelope_idxs, other_redox_envelopes). Each
# worker rebuilds Site objects against its freshly loaded structure so the
# envelope-aware cluster cutout is honored in the self-energy step.
_W_SITE_ENVELOPES: Dict[Tuple[str, int],
                          Tuple[List[int], List[List[int]]]] = {}


def _self_init(prm: str, rst: str, background_pkl: bytes,
                work_root: str, pbsa_params: Dict[str, Any],
                cluster_radius: Optional[float], bond_depth: int,
                site_envelopes_pkl: Optional[bytes] = None,
                model_cache_dir: Optional[str] = None):
    import pickle, warnings
    from parmed.exceptions import AmberWarning
    # Suppress parmed's contiguity warning in worker processes (workflow-
    # level catch_warnings doesn't cross the subprocess boundary).
    warnings.filterwarnings(
        "ignore",
        message=r".*Molecule atoms are not contiguous.*",
        category=AmberWarning)
    global _W_PRM, _W_RST, _W_BACKGROUND, _W_WORK_ROOT, _W_PBSA_PARAMS
    global _W_CLUSTER_RADIUS, _W_BOND_DEPTH, _W_SITE_ENVELOPES
    global _W_MODEL_CACHE_DIR
    _W_PRM = Path(prm)
    _W_RST = Path(rst)
    _W_BACKGROUND = pickle.loads(background_pkl)
    _W_WORK_ROOT = Path(work_root)
    _W_PBSA_PARAMS = pbsa_params
    _W_CLUSTER_RADIUS = cluster_radius
    _W_BOND_DEPTH = bond_depth
    _W_SITE_ENVELOPES = (pickle.loads(site_envelopes_pkl)
                            if site_envelopes_pkl else {})
    _W_MODEL_CACHE_DIR = Path(model_cache_dir) if model_cache_dir else None


def _self_task(args: Tuple[str, int, str]):
    resname, resnum, work_subdir = args
    residue = get_residue(resname)
    # Reconstruct Site against freshly loaded structure if envelope info given.
    site = None
    if (resname, resnum) in _W_SITE_ENVELOPES:
        import parmed as pmd
        from .sites import Site
        env_idxs, other_envs = _W_SITE_ENVELOPES[(resname, resnum)]
        s = pmd.load_file(str(_W_PRM), str(_W_RST))
        envelope = [s.residues[i] for i in env_idxs]
        target = next(r for r in envelope
                      if r.name == resname and r.number + 1 == resnum)
        site = Site(titrating_residue=target,
                    envelope_residues=envelope,
                    other_redox_envelopes=other_envs)
    se = compute_state_energies(
        cluster_prmtop=_W_PRM,
        cluster_rst7=_W_RST,
        target_resname=resname,
        target_resnum_1based=resnum,
        residue=residue,
        model_cache_dir=_W_MODEL_CACHE_DIR,
        work_dir=_W_WORK_ROOT / "self" / work_subdir,
        background=_W_BACKGROUND,
        cluster_radius=_W_CLUSTER_RADIUS,
        bond_depth=_W_BOND_DEPTH,
        site=site,
        **_W_PBSA_PARAMS,
    )
    return (resname, resnum), se


# ---------------------------------------------------------------------------
# Pair-coupling worker pool (PATH B). Mirrors the self-energy pool above: each
# worker holds the structure paths + shared params as globals and computes one
# pair's full set of PB calls (ref + per-chemistry singles + both) in its own
# cluster cutout, returning the W contributions for that pair. Parallelizing
# over pairs is what keeps Phase 2 from running serially.
# ---------------------------------------------------------------------------
_WP_PRM: Path = None
_WP_RST: Path = None
_WP_BACKGROUND: Dict[SiteKey, "ProtonationState"] = {}
_WP_WORK_ROOT: Path = None
_WP_PBSA_PARAMS: Dict[str, Any] = {}
_WP_CLUSTER_RADIUS: Optional[float] = None
_WP_BOND_DEPTH: int = 3
_WP_SITE_ENVELOPES: Dict[Tuple[str, int], Tuple[List[int], List[List[int]]]] = {}
_WP_FULL_STRUCT = None   # lazily loaded once per worker


def _pair_init(prm: str, rst: str, background_pkl: bytes,
               work_root: str, pbsa_params: Dict[str, Any],
               cluster_radius: Optional[float], bond_depth: int,
               site_envelopes_pkl: Optional[bytes] = None):
    import pickle, warnings
    from parmed.exceptions import AmberWarning
    warnings.filterwarnings(
        "ignore",
        message=r".*Molecule atoms are not contiguous.*",
        category=AmberWarning)
    global _WP_PRM, _WP_RST, _WP_BACKGROUND, _WP_WORK_ROOT, _WP_PBSA_PARAMS
    global _WP_CLUSTER_RADIUS, _WP_BOND_DEPTH, _WP_SITE_ENVELOPES, _WP_FULL_STRUCT
    _WP_PRM = Path(prm)
    _WP_RST = Path(rst)
    _WP_BACKGROUND = pickle.loads(background_pkl)
    _WP_WORK_ROOT = Path(work_root)
    _WP_PBSA_PARAMS = pbsa_params
    _WP_CLUSTER_RADIUS = cluster_radius
    _WP_BOND_DEPTH = bond_depth
    _WP_SITE_ENVELOPES = (pickle.loads(site_envelopes_pkl)
                          if site_envelopes_pkl else {})
    _WP_FULL_STRUCT = None


def _pair_task(args: Tuple[str, int, str, int]):
    """Compute the coupling block for ONE site pair in its own cluster cutout.

    args = (resname_i, resnum_i, resname_j, resnum_j). Returns
    ((i_key, ci), (j_key, cj)) -> W for every non-ref chemistry combination,
    so the parent can fill out.W symmetrically. Mirrors the serial PATH B body.
    """
    global _WP_FULL_STRUCT
    rn_i, num_i, rn_j, num_j = args
    site_i = make_coupling_site(rn_i, num_i)
    site_j = make_coupling_site(rn_j, num_j)

    if _WP_FULL_STRUCT is None:
        _WP_FULL_STRUCT = pmd.load_file(str(_WP_PRM), str(_WP_RST))
    full_struct = _WP_FULL_STRUCT

    # Build the per-pair cluster cutout (envelope-aware when available).
    env_i = env_j = None
    if (rn_i, num_i) in _WP_SITE_ENVELOPES and (rn_j, num_j) in _WP_SITE_ENVELOPES:
        from .sites import Site
        ei_idxs, ei_other = _WP_SITE_ENVELOPES[(rn_i, num_i)]
        ej_idxs, ej_other = _WP_SITE_ENVELOPES[(rn_j, num_j)]
        env_i = Site(
            titrating_residue=next(
                r for r in (full_struct.residues[k] for k in ei_idxs)
                if r.name == rn_i and r.number + 1 == num_i),
            envelope_residues=[full_struct.residues[k] for k in ei_idxs],
            other_redox_envelopes=ei_other)
        env_j = Site(
            titrating_residue=next(
                r for r in (full_struct.residues[k] for k in ej_idxs)
                if r.name == rn_j and r.number + 1 == num_j),
            envelope_residues=[full_struct.residues[k] for k in ej_idxs],
            other_redox_envelopes=ej_other)
        keep = cluster_for_site_pair(full_struct, env_i, env_j,
                                     R=_WP_CLUSTER_RADIUS, bond_depth=_WP_BOND_DEPTH)
    else:
        ti = next(r for r in full_struct.residues
                  if r.name == rn_i and (r.number + 1) == num_i)
        tj = next(r for r in full_struct.residues
                  if r.name == rn_j and (r.number + 1) == num_j)
        keep = cluster_for_pair(full_struct, ti, tj,
                                R=_WP_CLUSTER_RADIUS, bond_depth=_WP_BOND_DEPTH)

    label_base = f"{rn_i}{num_i:04d}_{rn_j}{num_j:04d}"
    wd = _WP_WORK_ROOT / "pair_cluster" / label_base

    e_ref_ij = _pb_with_overrides(
        _WP_PRM, _WP_RST, overrides={}, background=_WP_BACKGROUND,
        work_dir=wd, label="ref", residue_idxs_keep=keep, **_WP_PBSA_PARAMS)

    w_block: Dict[Tuple[Tuple[SiteKey, ChemKey], Tuple[SiteKey, ChemKey]], float] = {}
    for ci in site_i.chemistries:
        if ci == site_i.ref_chem:
            continue
        for cj in site_j.chemistries:
            if cj == site_j.ref_chem:
                continue
            e_si = _pb_with_overrides(
                _WP_PRM, _WP_RST,
                overrides={site_i.key: site_i.chem_to_state[ci]},
                background=_WP_BACKGROUND, work_dir=wd,
                label=f"si_n{ci[0]}_p{ci[1]:.1f}",
                residue_idxs_keep=keep, **_WP_PBSA_PARAMS)
            e_sj = _pb_with_overrides(
                _WP_PRM, _WP_RST,
                overrides={site_j.key: site_j.chem_to_state[cj]},
                background=_WP_BACKGROUND, work_dir=wd,
                label=f"sj_n{cj[0]}_p{cj[1]:.1f}",
                residue_idxs_keep=keep, **_WP_PBSA_PARAMS)
            e_both = _pb_with_overrides(
                _WP_PRM, _WP_RST,
                overrides={site_i.key: site_i.chem_to_state[ci],
                           site_j.key: site_j.chem_to_state[cj]},
                background=_WP_BACKGROUND, work_dir=wd,
                label=f"both_n{ci[0]}_n{cj[0]}",
                residue_idxs_keep=keep, **_WP_PBSA_PARAMS)
            w_block[((site_i.key, ci), (site_j.key, cj))] = (
                e_both - e_si - e_sj + e_ref_ij)
    return (site_i.key, site_j.key), w_block


def estimate_self_energy_pb_calls(
        sites: List[CouplingSite],
        precomputed_keys: Optional[set] = None) -> int:
    """Total PB calls compute_self_energies will issue.

    Each non-cached site contributes 2 PB calls per state (cluster + model).
    Sites whose key is in `precomputed_keys` are reused from a prior step's
    pickle and contribute zero calls.
    """
    precomputed_keys = precomputed_keys or set()
    total = 0
    for s in sites:
        if s.key in precomputed_keys:
            continue
        total += 2 * len(s.residue.states)
    return total


def compute_self_energies(sites: List[CouplingSite],
                            cluster_prmtop: Path,
                            cluster_rst7:   Path,
                            *,
                            workers: int = 4,
                            work_root: Path = Path("coupling_work"),
                            pbsa_params: Optional[Dict[str, Any]] = None,
                            cluster_radius: Optional[float] = None,
                            bond_depth: int = 3,
                            site_envelopes: Optional[Dict[Tuple[str, int], Site]] = None,
                            model_cache_dir: Optional[Path] = None,
                            extra_background: Optional[Dict[SiteKey, Any]] = None,
                            on_progress: Optional[Any] = None,
                            precomputed_state_energies: Optional[
                                Dict[SiteKey, Dict[int, Dict[str, Any]]]] = None,
                            ) -> SelfEnergies:
    """Compute intrinsic ΔΔG for every (site, chemistry).

    All other sites are held at their reference chemistry's representative
    state during each PB call. Parallelizes over sites.

    `cluster_radius`: if not None, each PB call uses a single-site cutout
    of that radius around the target. Required for any system >~3000 atoms.

    `site_envelopes`: optional map (resname, resnum) -> Site. When provided,
    each worker rebuilds the Site against its freshly loaded structure
    (using stable parmed residue indices) so the envelope-aware cluster
    cutout is honored. Required to keep integer charge across PB clusters
    that include MCPB-parameterized multi-residue groups (e.g., hemes).

    `precomputed_state_energies`: optional map (resname, resnum) -> the
    `state_energies` dict from a prior `compute_state_energies` /
    `compute_pka_multistate` run. When provided, sites whose key is in the
    map skip the PB calls entirely and reuse the cached state_energies.
    Step 3 (single-site pKa) and step 4 Phase 1 (self-energies) issue
    byte-identical PB calls — same background (all-neutral), same cluster
    radius, same site envelopes, same pbsa params — so the cache hit is
    safe whenever step 3's pickle is still valid. Sites absent from the
    map fall back to fresh PB calls.
    """
    import pickle
    pbsa_params = pbsa_params or {}
    work_root = Path(work_root)
    work_root.mkdir(parents=True, exist_ok=True)
    # Route non-fatal pbsa exit-time-crash salvages to a per-run log (workers
    # inherit this env var) instead of flooding the console with one warning
    # per PB call. See pb_backend._record_pbsa_salvage.
    os.environ["PROPREP_PBSA_SALVAGE_LOG"] = str(work_root / "pbsa_salvage.log")
    precomputed_state_energies = precomputed_state_energies or {}

    background = {s.key: s.chem_to_state[s.ref_chem] for s in sites}
    # Fold in residues that are excluded from titration but must be held at a
    # fixed (non-reference) charge in every PB call — e.g. metal-coordinating
    # acids fixed deprotonated. These are not CouplingSites (never flipped), so
    # they only enter as background. Mirrors the pbt-3 single-site background.
    if extra_background:
        background.update(extra_background)
    bg_pkl = pickle.dumps(background)

    # Distill site_envelopes to a worker-friendly form: parmed residues
    # don't survive pickling, but their indices are stable across reloads
    # of the same prmtop.
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
    model_cache_arg = str(model_cache_dir) if model_cache_dir else None

    out = SelfEnergies(sites=sites)

    # Partition sites into cached (skip PB) and fresh (run PB).
    fresh_sites: List[CouplingSite] = []
    n_cached = 0
    for s in sites:
        if s.key in precomputed_state_energies:
            out.state_energies_full[s.key] = precomputed_state_energies[s.key]
            n_cached += 1
        else:
            fresh_sites.append(s)

    items = [(s.resname, s.resnum, f"{s.resname}{s.resnum:04d}")
             for s in fresh_sites]
    n_fresh = len(items)
    n_sites = len(sites)
    if on_progress is not None and n_cached > 0:
        # Surface cache hits to the caller's progress tracker so the
        # totals line up with `estimate_self_energy_pb_calls`.
        on_progress(n_cached, n_sites, f"{n_cached} cached from step 3")

    if n_fresh == 0:
        # Everything came from the cache; still need the post-extraction.
        pass
    elif workers > 1:
        with ProcessPoolExecutor(
                max_workers=workers,
                initializer=_self_init,
                initargs=(str(cluster_prmtop), str(cluster_rst7),
                          bg_pkl, str(work_root), pbsa_params,
                          cluster_radius, bond_depth, envelopes_pkl,
                          model_cache_arg)) as pool:
            futures = {pool.submit(_self_task, it): it for it in items}
            done = n_cached
            for fut in as_completed(futures):
                key, se = fut.result()
                out.state_energies_full[key] = se
                done += 1
                if on_progress is not None:
                    on_progress(done, n_sites, f"{key[0]}-{key[1]}")
    else:
        _self_init(str(cluster_prmtop), str(cluster_rst7),
                    bg_pkl, str(work_root), pbsa_params,
                    cluster_radius, bond_depth, envelopes_pkl,
                    model_cache_arg)
        done = n_cached
        for it in items:
            key, se = _self_task(it)
            out.state_energies_full[key] = se
            done += 1
            if on_progress is not None:
                on_progress(done, n_sites, f"{key[0]}-{key[1]}")

    # Extract chemistry-group representative ΔΔG
    for site in sites:
        se = out.state_energies_full[site.key]
        for chem in site.chemistries:
            # Find the lowest-G representative within this chemistry at pH=7
            # (any pH; since members share (n, pka_corr), G ordering is pH-independent)
            G_at_pH7 = state_free_energies_at_pH(se, 7.0)
            members = [i for i, info in se.items()
                       if (info["prot_count"], info["pka_corr"]) == chem]
            best_idx = min(members, key=lambda i: G_at_pH7[i])
            out.intrinsic_ddG[(site.key, chem)] = se[best_idx]["ddG_PB_minus_model"]
    return out


# ---------------------------------------------------------------------------
# Pair coupling matrix W_ij(s, t)
# ---------------------------------------------------------------------------

@dataclass
class CouplingMatrix:
    """Full pair coupling tensor for an N-site titratable system.

    `W[(site_i_key, chem_i)][(site_j_key, chem_j)]` is the pair coupling
    energy in kcal/mol. By construction:
      - W[i][j] = W[j][i] (symmetric)
      - W[i, ref][j, *] = 0 and W[i, *][j, ref] = 0 (reference subtracted)
      - W[i, ref][j, ref] = 0 by definition of reference

    Stored as a nested dict for human-debuggable access; for hot loops in
    the MC, you'd want a 2D numpy array.
    """
    sites: List[CouplingSite]
    W: Dict[Tuple[SiteKey, ChemKey], Dict[Tuple[SiteKey, ChemKey], float]] = field(default_factory=dict)
    e_reference: float = 0.0   # PB energy of all-reference background


def _pb_with_overrides(cluster_prmtop: Path, cluster_rst7: Path,
                        overrides: Dict[SiteKey, ProtonationState],
                        background: Dict[SiteKey, ProtonationState],
                        work_dir: Path, label: str,
                        residue_idxs_keep: Optional[List[int]] = None,
                        **pbsa_params) -> float:
    """Run PB with specified site states; return total electrostatic E.

    If `residue_idxs_keep` is provided, slice the structure to those
    0-based residue indices first (cluster cutout). Background sites that
    fall outside the slice are silently skipped (they don't exist in the
    cluster, so there's nothing to apply); overrides must always be in the
    cluster and raise loudly if missing.
    """
    s = pmd.load_file(str(cluster_prmtop), str(cluster_rst7))
    if residue_idxs_keep is not None:
        s = slice_structure(s, residue_idxs_keep)

    # Background: silent-skip when the residue isn't in the (possibly sliced)
    # cluster. Equivalent to intrinsic.apply_background's semantics.
    for (resname, resnum), state in background.items():
        if (resname, resnum) in overrides:
            continue
        for r in s.residues:
            if r.name == resname and (r.number + 1) == resnum:
                for a in r.atoms:
                    if a.name in state.charges:
                        a.charge = state.charges[a.name]
                break

    # Overrides (target sites for the pair / single perturbation): must
    # exist in the cluster; raise loudly if not.
    for (resname, resnum), state in overrides.items():
        apply_state_charges(s, resname, resnum, state)
    e = run_pbsa(s, work_dir, label, fillratio=2.0, **pbsa_params)
    return total_electrostatic(e)


def _site_centroid(structure: pmd.Structure, resname: str, resnum: int):
    """Centroid of a residue's heavy atoms (used for distance cutoff)."""
    for r in structure.residues:
        if r.name == resname and (r.number + 1) == resnum:
            xs, ys, zs, n = 0.0, 0.0, 0.0, 0
            for a in r.atoms:
                if a.atomic_number > 1:
                    xs += a.xx; ys += a.xy; zs += a.xz; n += 1
            return (xs/n, ys/n, zs/n)
    raise ValueError(f"{resname}-{resnum} not found")


def estimate_pair_pb_calls(sites: List[CouplingSite],
                            cluster_prmtop: Path,
                            cluster_rst7:   Path,
                            *,
                            cluster_radius: Optional[float] = None,
                            distance_cutoff: Optional[float] = None,
                            pbsa_params: Optional[Dict[str, Any]] = None,
                            ) -> Dict[str, int]:
    """Pre-compute the pair-coupling work that compute_pair_couplings will do.

    Returns a dict with:
        n_pairs_total      — N·(N-1)/2 site pairs
        n_pairs_kept       — pairs within the distance cutoff
        n_pairs_skipped    — pairs beyond the cutoff (W set to 0)
        distance_cutoff_A  — the cutoff in Å (auto-derived from istrng if needed)
        n_pb_calls         — total PB calls the function will issue
    """
    pbsa_params = pbsa_params or {}
    if distance_cutoff is None:
        istrng = pbsa_params.get("istrng", 150.0)
        distance_cutoff = recommended_pair_cutoff(istrng) if cluster_radius else 1e6

    # Per-site non-ref chemistry counts (everything else is constant per site)
    n_nonref = {s.key: sum(1 for c in s.chemistries if c != s.ref_chem) for s in sites}

    if cluster_radius is None:
        # PATH A: shared full cluster.
        # Calls: 1 reference + Σ_i n_nonref(i) singles + Σ_{i<j} n_nonref(i)·n_nonref(j) pairs.
        n_singles = sum(n_nonref.values())
        n_pair_calls = 0
        for i, si in enumerate(sites):
            for sj in sites[i+1:]:
                n_pair_calls += n_nonref[si.key] * n_nonref[sj.key]
        n_pb = 1 + n_singles + n_pair_calls
        return {
            "n_pairs_total":   len(sites) * (len(sites) - 1) // 2,
            "n_pairs_kept":    len(sites) * (len(sites) - 1) // 2,
            "n_pairs_skipped": 0,
            "distance_cutoff_A": float("inf"),
            "n_pb_calls":      n_pb,
        }

    # PATH B: per-pair cluster cutouts. Need centroids to count kept pairs.
    full = pmd.load_file(str(cluster_prmtop), str(cluster_rst7))
    centroids = {s.key: _site_centroid(full, s.resname, s.resnum) for s in sites}
    n_kept = 0
    n_skipped = 0
    n_pb = 0
    for i, si in enumerate(sites):
        for sj in sites[i+1:]:
            ci_xyz = centroids[si.key]; cj_xyz = centroids[sj.key]
            d = sum((ci_xyz[k]-cj_xyz[k])**2 for k in range(3)) ** 0.5
            if d > distance_cutoff:
                n_skipped += 1
                continue
            n_kept += 1
            # Per pair: 1 ref + 3 PB calls per (ci, cj) non-ref combination
            n_pb += 1 + 3 * n_nonref[si.key] * n_nonref[sj.key]
    return {
        "n_pairs_total":   n_kept + n_skipped,
        "n_pairs_kept":    n_kept,
        "n_pairs_skipped": n_skipped,
        "distance_cutoff_A": distance_cutoff,
        "n_pb_calls":      n_pb,
    }


def compute_pair_couplings(sites: List[CouplingSite],
                            cluster_prmtop: Path,
                            cluster_rst7:   Path,
                            *,
                            workers: int = 4,
                            work_root: Path = Path("coupling_work"),
                            pbsa_params: Optional[Dict[str, Any]] = None,
                            cluster_radius: Optional[float] = None,
                            bond_depth: int = 3,
                            distance_cutoff: Optional[float] = None,
                            site_envelopes: Optional[Dict[Tuple[str, int], Site]] = None,
                            on_progress: Optional[Any] = None,
                            ) -> CouplingMatrix:
    """Compute W_ij(s,t) for every site pair and every (non-ref) chemistry pair.

    `cluster_radius`:
      None (default) → use full structure for all PB calls (one shared
        cluster). Cost: 1 + Σ(n_chem-1) + Σ_pair (n_chem_i-1)(n_chem_j-1).
      float → per-pair cluster cutout (union of two single-site spheres).
        Required for systems too big to PB whole. Cost: 4 PB calls per
        kept pair (each within its small cluster).

    `distance_cutoff`:
      None → auto-compute from istrng in pbsa_params via Debye length.
      float → explicit cutoff in Å (skip pairs with site-site distance >).
      Pairs with W_ij set to 0 in the returned matrix.

    For the bundle (122 PRNs at 650 mM) the auto cutoff is ~11 Å,
    keeping mostly intra-heme pairs. For BPTI/RNase (small), cluster_radius
    None preserves the original behavior.
    """
    pbsa_params = pbsa_params or {}
    work_root = Path(work_root)
    work_root.mkdir(parents=True, exist_ok=True)
    # Route non-fatal pbsa exit-time-crash salvages to a per-run log (workers
    # inherit this env var) instead of flooding the console with one warning
    # per PB call. See pb_backend._record_pbsa_salvage.
    os.environ["PROPREP_PBSA_SALVAGE_LOG"] = str(work_root / "pbsa_salvage.log")

    background = {s.key: s.chem_to_state[s.ref_chem] for s in sites}

    # Auto distance cutoff from istrng if not given
    if distance_cutoff is None:
        istrng = pbsa_params.get("istrng", 150.0)
        distance_cutoff = recommended_pair_cutoff(istrng) if cluster_radius else 1e6

    # Initialize coupling matrix (with all pairs at 0 for cutoff/below-threshold)
    out = CouplingMatrix(sites=sites, e_reference=0.0)
    for site in sites:
        for chem in site.chemistries:
            out.W[(site.key, chem)] = {}

    # ====================================================================
    # PATH A: full structure (cluster_radius is None) — original O(N²) but
    # all PB calls share one cluster, so reference + singles are amortized.
    # ====================================================================
    if cluster_radius is None:
        n_nonref = {s.key: sum(1 for c in s.chemistries if c != s.ref_chem) for s in sites}
        n_singles = sum(n_nonref.values())
        n_pair_calls = sum(n_nonref[si.key] * n_nonref[sj.key]
                            for i, si in enumerate(sites) for sj in sites[i+1:])
        total_pb = 1 + n_singles + n_pair_calls
        done = 0

        e_ref = _pb_with_overrides(
            cluster_prmtop, cluster_rst7,
            overrides={}, background=background,
            work_dir=work_root / "pair", label="reference",
            **pbsa_params)
        out.e_reference = e_ref
        done += 1
        if on_progress is not None:
            on_progress(done, total_pb, "reference")

        e_single: Dict[Tuple[SiteKey, ChemKey], float] = {}
        for site in sites:
            for chem in site.chemistries:
                if chem == site.ref_chem: continue
                label = f"single_{site.resname}{site.resnum:04d}_n{chem[0]}_p{chem[1]:.1f}"
                e = _pb_with_overrides(
                    cluster_prmtop, cluster_rst7,
                    overrides={site.key: site.chem_to_state[chem]},
                    background=background,
                    work_dir=work_root / "pair", label=label,
                    **pbsa_params)
                e_single[(site.key, chem)] = e
                done += 1
                if on_progress is not None:
                    on_progress(done, total_pb,
                                  f"single {site.resname}-{site.resnum}")

        e_pair: Dict[Tuple[Tuple[SiteKey,ChemKey], Tuple[SiteKey,ChemKey]], float] = {}
        for i, site_i in enumerate(sites):
            for site_j in sites[i+1:]:
                for ci in site_i.chemistries:
                    if ci == site_i.ref_chem: continue
                    for cj in site_j.chemistries:
                        if cj == site_j.ref_chem: continue
                        label = (f"pair_{site_i.resname}{site_i.resnum:04d}_n{ci[0]}_"
                                 f"x{site_j.resname}{site_j.resnum:04d}_n{cj[0]}")
                        e = _pb_with_overrides(
                            cluster_prmtop, cluster_rst7,
                            overrides={site_i.key: site_i.chem_to_state[ci],
                                       site_j.key: site_j.chem_to_state[cj]},
                            background=background,
                            work_dir=work_root / "pair", label=label,
                            **pbsa_params)
                        e_pair[((site_i.key, ci), (site_j.key, cj))] = e
                        done += 1
                        if on_progress is not None:
                            on_progress(done, total_pb,
                                          f"pair {site_i.resname}-{site_i.resnum} × "
                                          f"{site_j.resname}-{site_j.resnum}")

        for (ki_ci, kj_cj), e_p in e_pair.items():
            e_si = e_single[ki_ci]
            e_sj = e_single[kj_cj]
            W = e_p - e_si - e_sj + e_ref
            out.W[ki_ci][kj_cj] = W
            out.W[kj_cj][ki_ci] = W   # symmetric
        return out

    # ====================================================================
    # PATH B: cluster_radius set — per-pair clusters with distance cutoff.
    # 4 PB calls per kept pair within its own cluster_ij.
    # ====================================================================
    full_struct = pmd.load_file(str(cluster_prmtop), str(cluster_rst7))
    centroids = {s.key: _site_centroid(full_struct, s.resname, s.resnum)
                 for s in sites}

    # Pre-pass: enumerate kept pairs and total PB-call count for progress.
    n_nonref = {s.key: sum(1 for c in s.chemistries if c != s.ref_chem) for s in sites}
    pairs_to_run = []
    n_skipped = 0
    for i, site_i in enumerate(sites):
        for site_j in sites[i+1:]:
            ci_xyz = centroids[site_i.key]
            cj_xyz = centroids[site_j.key]
            d = sum((ci_xyz[k]-cj_xyz[k])**2 for k in range(3)) ** 0.5
            if d > distance_cutoff:
                n_skipped += 1
                continue
            pairs_to_run.append((site_i, site_j))
    n_kept = len(pairs_to_run)
    total_pb = sum(1 + 3 * n_nonref[si.key] * n_nonref[sj.key]
                    for si, sj in pairs_to_run)
    done = 0

    pair_label = lambda si, sj: (f"{si.resname}-{si.resnum} × "
                                   f"{sj.resname}-{sj.resnum}")

    # Each kept pair is independent (its own cluster cutout + 4 PB calls), so
    # dispatch one _pair_task per pair across the worker pool. This is what
    # parallelizes Phase 2; previously this loop ran every PB call serially.
    import pickle
    bg_pkl = pickle.dumps(background)
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

    items = [(si.resname, si.resnum, sj.resname, sj.resnum)
             for (si, sj) in pairs_to_run]
    label_for = {(si.resname, si.resnum, sj.resname, sj.resnum):
                 pair_label(si, sj) for (si, sj) in pairs_to_run}

    def _store(result):
        (i_key, j_key), w_block = result
        for (ki_ci, kj_cj), W in w_block.items():
            out.W[ki_ci][kj_cj] = W
            out.W[kj_cj][ki_ci] = W   # symmetric

    initargs = (str(cluster_prmtop), str(cluster_rst7), bg_pkl,
                str(work_root), pbsa_params, cluster_radius, bond_depth,
                envelopes_pkl)
    if workers > 1:
        with ProcessPoolExecutor(
                max_workers=workers,
                initializer=_pair_init,
                initargs=initargs) as pool:
            futures = {pool.submit(_pair_task, it): it for it in items}
            for fut in as_completed(futures):
                _store(fut.result())
                done += 1
                if on_progress is not None:
                    it = futures[fut]
                    on_progress(done, n_kept,
                                  f"[pair {done}/{n_kept}] "
                                  f"{label_for[it]}")
    else:
        _pair_init(*initargs)
        for it in items:
            _store(_pair_task(it))
            done += 1
            if on_progress is not None:
                on_progress(done, n_kept,
                              f"[pair {done}/{n_kept}] {label_for[it]}")

    print(f"  pair couplings: {n_kept} pairs computed, {n_skipped} pairs "
          f"skipped (d > {distance_cutoff:.2f} Å)")
    return out


# ===========================================================================
# O(N) coupling via PB potential maps (Beroza/Bashford 1991 trick)
# ===========================================================================

def _zero_all_charges(s: pmd.Structure) -> None:
    """Set every atom's charge to 0 (in place)."""
    for a in s.atoms:
        a.charge = 0.0


def _delta_charges_for(site: CouplingSite, chem: ChemKey) -> Dict[str, float]:
    """Return {atom_name: q_perturbed - q_ref} for one chemistry of one site.

    The "ref" charges are from the site's reference chemistry's
    representative state. The "perturbed" charges are from the
    representative state of `chem`.
    """
    ref_state = site.chem_to_state[site.ref_chem]
    pert_state = site.chem_to_state[chem]
    out: Dict[str, float] = {}
    for atom_name in pert_state.charges:
        ref_q = ref_state.charges.get(atom_name, 0.0)
        out[atom_name] = pert_state.charges[atom_name] - ref_q
    return out


def compute_pair_couplings_o_n(
        sites: List[CouplingSite],
        cluster_prmtop: Path,
        cluster_rst7:   Path,
        *,
        work_root: Path = Path("coupling_o_n_work"),
        pbsa_params: Optional[Dict[str, Any]] = None,
        cluster_radius: Optional[float] = None,
        bond_depth: int = 3,
        distance_cutoff: Optional[float] = None,
        ) -> CouplingMatrix:
    """O(N) computation of W_ij(s,t) via PB potential maps.

    KNOWN ACCURACY BUG (2026-05-07): cross-validation against the explicit
    method on BPTI 3-site shows W_ij values ~5–10× too small. Root cause
    not yet localized — likely related to how pbsa writes phi maps on its
    focused grid (smoothing of singular field components near the
    perturbed atoms, or some convention mismatch in the dumped phi). DO
    NOT USE for production pKa calculations until cross-validation passes.
    The explicit method (`compute_pair_couplings`) with cluster cutout +
    distance cutoff is bundle-feasible and produces correct numbers.

    Algorithm (subtraction-based — most robust against pbsa numerics):
      1. Run PB once with all sites at REFERENCE chemistry; save phi_ref map.
      2. For each (site i, non-ref chemistry s_i):
           - Run PB with site i at chem s_i, all others at REF; save phi_pert map.
           - Compute Δφ_i = phi_pert − phi_ref (subtract grid values).
      3. For each pair (i, j) within distance_cutoff and each (s_i, s_j):
           W_ij(s_i, s_j) = Σ over j's atoms: Δq_j(atom, s_j) · Δφ_i(r_atom)

    Cost: 1 + Σ(n_chem−1) PB calls (linear in N), vs O(N²) for the
    explicit method's pair calls. Mathematically equivalent in linearized
    PB; cross-validate with `compute_pair_couplings` on a small system.
    """
    import numpy as np
    from .potentials import read_pbsa_phi, interpolate, PhiGrid

    pbsa_params = pbsa_params or {}
    work_root = Path(work_root)
    work_root.mkdir(parents=True, exist_ok=True)
    # Route non-fatal pbsa exit-time-crash salvages to a per-run log (workers
    # inherit this env var) instead of flooding the console with one warning
    # per PB call. See pb_backend._record_pbsa_salvage.
    os.environ["PROPREP_PBSA_SALVAGE_LOG"] = str(work_root / "pbsa_salvage.log")

    if distance_cutoff is None:
        istrng = pbsa_params.get("istrng", 150.0)
        distance_cutoff = recommended_pair_cutoff(istrng) if cluster_radius else 1e6

    # Load full structure once (for centroids + slicing per site)
    full_struct = pmd.load_file(str(cluster_prmtop), str(cluster_rst7))
    centroids = {s.key: _site_centroid(full_struct, s.resname, s.resnum)
                 for s in sites}
    background = {s.key: s.chem_to_state[s.ref_chem] for s in sites}

    # Initialize coupling matrix
    out = CouplingMatrix(sites=sites, e_reference=0.0)
    for site in sites:
        for chem in site.chemistries:
            out.W[(site.key, chem)] = {}

    # ====================================================================
    # Phase 1a: reference PB (all sites at ref chemistry)
    # ====================================================================
    phi_dir = work_root / "phi_maps"
    phi_dir.mkdir(parents=True, exist_ok=True)
    s_ref = pmd.load_file(str(cluster_prmtop), str(cluster_rst7))
    if cluster_radius is None:
        # Apply background ref states for the reference call
        for (rn, ru), state in background.items():
            apply_state_charges(s_ref, rn, ru, state)
        run_pbsa(s_ref, phi_dir, "ref",
                  phiout=1, phiform=1, fillratio=2.0, **pbsa_params)
        phi_ref_global = read_pbsa_phi(phi_dir / "ref.phi")
    else:
        phi_ref_global = None  # per-site clusters use per-site ref maps

    # ====================================================================
    # Phase 1b: per-(site, chem) PB call → phi_perturbed → Δφ_i
    # ====================================================================
    phi_delta_paths: Dict[Tuple[SiteKey, ChemKey], PhiGrid] = {}
    site_clusters: Dict[SiteKey, List[int]] = {}
    for site in sites:
        if cluster_radius is not None:
            from .cluster import cluster_for_target
            target_in_full = next(r for r in full_struct.residues
                                   if r.name == site.resname
                                   and (r.number+1) == site.resnum)
            keep = cluster_for_target(full_struct, target_in_full,
                                        R=cluster_radius, bond_depth=bond_depth)
            site_clusters[site.key] = keep
            # Build per-site reference map
            s_ref_local = slice_structure(full_struct, keep)
            for (rn, ru), state in background.items():
                apply_state_charges(s_ref_local, rn, ru, state)
            wd = phi_dir / f"{site.resname}{site.resnum:04d}"
            wd.mkdir(parents=True, exist_ok=True)
            run_pbsa(s_ref_local, wd, "ref",
                      phiout=1, phiform=1, fillratio=2.0, **pbsa_params)
            phi_ref_local = read_pbsa_phi(wd / "ref.phi")
        else:
            phi_ref_local = phi_ref_global

        for chem in site.chemistries:
            if chem == site.ref_chem:
                continue
            label = f"phi_{site.resname}{site.resnum:04d}_n{chem[0]}_p{chem[1]:.1f}"
            wd = (phi_dir / f"{site.resname}{site.resnum:04d}"
                   if cluster_radius is not None else phi_dir)
            wd.mkdir(parents=True, exist_ok=True)

            # Build perturbed structure
            if cluster_radius is not None:
                s_pert = slice_structure(full_struct, site_clusters[site.key])
            else:
                s_pert = pmd.load_file(str(cluster_prmtop), str(cluster_rst7))
            for (rn, ru), state in background.items():
                apply_state_charges(s_pert, rn, ru, state)
            apply_state_charges(s_pert, site.resname, site.resnum,
                                  site.chem_to_state[chem])
            run_pbsa(s_pert, wd, label,
                      phiout=1, phiform=1, fillratio=2.0, **pbsa_params)
            phi_pert = read_pbsa_phi(wd / f"{label}.phi")
            # Δφ = phi_pert − phi_ref (must share the grid; pbsa fixes grid
            # by atom positions which are unchanged across charge swaps)
            if (phi_pert.dims != phi_ref_local.dims or
                    abs(phi_pert.h - phi_ref_local.h) > 1e-6 or
                    any(abs(phi_pert.origin[k] - phi_ref_local.origin[k]) > 1e-6
                        for k in range(3))):
                raise RuntimeError(
                    f"phi grid mismatch between perturbed and reference: "
                    f"dims={phi_pert.dims} vs {phi_ref_local.dims}, "
                    f"h={phi_pert.h} vs {phi_ref_local.h}, "
                    f"origin={phi_pert.origin} vs {phi_ref_local.origin}")
            from .potentials import PhiGrid as PG
            delta_phi = PG(h=phi_pert.h, origin=phi_pert.origin,
                            dims=phi_pert.dims,
                            phi=phi_pert.phi - phi_ref_local.phi)
            phi_delta_paths[(site.key, chem)] = delta_phi

    # ====================================================================
    # Phase 2: O(N²) interpolations — for each pair, evaluate Δφ_i at j's atoms
    # ====================================================================
    n_kept = 0
    n_skipped = 0
    # Cache per-site atom coords + per-(site,chem) Δq arrays
    site_atom_xyz: Dict[SiteKey, np.ndarray] = {}
    site_atom_names: Dict[SiteKey, List[str]] = {}
    for site in sites:
        target = next(r for r in full_struct.residues
                       if r.name == site.resname
                       and (r.number+1) == site.resnum)
        coords = []
        names = []
        for a in target.atoms:
            coords.append((a.xx, a.xy, a.xz))
            names.append(a.name)
        site_atom_xyz[site.key] = np.array(coords, dtype=np.float64)
        site_atom_names[site.key] = names

    for i, site_i in enumerate(sites):
        for site_j in sites[i+1:]:
            ci_xyz = centroids[site_i.key]
            cj_xyz = centroids[site_j.key]
            d = sum((ci_xyz[k]-cj_xyz[k])**2 for k in range(3)) ** 0.5
            if d > distance_cutoff:
                n_skipped += 1
                continue
            n_kept += 1
            for chem_i in site_i.chemistries:
                if chem_i == site_i.ref_chem: continue
                for chem_j in site_j.chemistries:
                    if chem_j == site_j.ref_chem: continue
                    delta_phi_i = phi_delta_paths[(site_i.key, chem_i)]
                    dq_j = _delta_charges_for(site_j, chem_j)
                    j_xyz = site_atom_xyz[site_j.key]
                    j_names = site_atom_names[site_j.key]
                    dq_array = np.array([dq_j.get(n, 0.0) for n in j_names])
                    phi_at_j = interpolate(delta_phi_i, j_xyz)
                    W = float(np.dot(dq_array, phi_at_j))
                    out.W[(site_i.key, chem_i)][(site_j.key, chem_j)] = W
                    out.W[(site_j.key, chem_j)][(site_i.key, chem_i)] = W

    n_pb_calls = (1 if cluster_radius is None else len(sites)) + \
                 sum(s.n_chemistries() - 1 for s in sites)
    print(f"  pair couplings (O(N) method): {n_kept} pairs evaluated, "
          f"{n_skipped} pairs skipped (d > {distance_cutoff:.2f} Å)")
    print(f"  PB calls used: {n_pb_calls}")
    return out


def coupling_summary(W: CouplingMatrix) -> str:
    """Human-readable table of |W_ij| max per site pair (for diagnostics).

    Only pairs with a non-negligible coupling are listed: a pair whose max
    |W| rounds to 0.000 kcal/mol at the displayed precision is omitted (most
    in-cutoff pairs are negligibly coupled and just bury the signal). The
    count of omitted near-zero pairs is reported at the end.
    """
    kT = RT_LN10 / 2.302585  # 0.5925 kcal/mol
    DISPLAY_EPS = 5e-4        # below this, max |W| prints as +0.000 → hide it
    lines = [f"Coupling matrix summary ({len(W.sites)} sites):"]
    n_hidden = 0
    for i, site_i in enumerate(W.sites):
        for site_j in W.sites[i+1:]:
            max_W = 0.0
            for chem_i in site_i.chemistries:
                if chem_i == site_i.ref_chem: continue
                for chem_j in site_j.chemistries:
                    if chem_j == site_j.ref_chem: continue
                    w = W.W.get((site_i.key, chem_i), {}).get((site_j.key, chem_j), 0.0)
                    if abs(w) > abs(max_W):
                        max_W = w
            if abs(max_W) < DISPLAY_EPS:
                n_hidden += 1
                continue
            tag = "  *" if abs(max_W) > 2 * kT else ""
            lines.append(f"  {site_i.resname}-{site_i.resnum:<3} ↔ "
                          f"{site_j.resname}-{site_j.resnum:<3}: "
                          f"max |W| = {max_W:+.3f} kcal/mol "
                          f"({max_W/kT:+.2f} kT){tag}")
    if len(lines) == 1:
        lines.append("  (no pair exceeds 0.000 kcal/mol — all couplings negligible)")
    if n_hidden:
        lines.append(f"  [{n_hidden} pair(s) with negligible coupling "
                     f"(|W| < 0.001 kcal/mol) omitted]")
    return "\n".join(lines)
