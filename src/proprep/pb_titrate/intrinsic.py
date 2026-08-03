"""Intrinsic ΔpKa per site via the Bashford 4-state cycle.

For one titratable site in a protein, the standard cycle is:

    ΔpKa = [(E_DEPROT − E_PROT)_protein − (E_DEPROT − E_PROT)_model]
            / (RT·ln10)
    pKa  = pKa_model + ΔpKa

The four PB calls are:
    cluster_PROT, cluster_DEPROT  (the protein with site i in each state)
    model_PROT,   model_DEPROT    (the ACE-X-NME tripeptide)

"Intrinsic" here means the pKa of site i with all other titratable sites
held in some background state — typically every other site at its
solution-pKa-determined state. This is the input to the precompute
layer of the unified solver; mean-field iteration or MC then folds in
the coupling matrix W_ij.

Generalized over residue type via `titrate.residues` (charge tables) and
`titrate.model_compounds` (on-demand ACE-X-NME builds for amino acids,
methyl-capped singletons for heme titratables like PRN).
"""
from __future__ import annotations

from pathlib import Path
from collections import defaultdict
from math import exp, log, log10
from typing import Any, Dict, List, Optional, Tuple

import parmed as pmd

from .residues import TitratableResidue, get_residue, ProtonationState
from .model_compounds import build_model
from .pb_backend import (
    RT_LN10,
    run_pbsa,
    total_electrostatic,
)
from .cluster import cluster_for_site, cluster_for_target, slice_structure
from .sites import Site


# A site key is (resname, resnum_1based). The background state map is
# {site_key: ProtonationState} — one ProtonationState per non-target site,
# applied to its atoms before each PB call. Use a state object (not just
# 'PROT'/'DEPROT') so multi-state residues (HIP) work uniformly.
SiteKey = Tuple[str, int]
BackgroundStateMap = Dict[SiteKey, ProtonationState]


def _state_delta(resname: str,
                  target_state: ProtonationState,
                  from_state: Optional[ProtonationState] = None,
                  ) -> Dict[str, float]:
    """Per-atom-name charge change from `from_state` to `target_state`.

    `from_state` defaults to the residue's STATE 0 — the form the prmtop
    is in when leaprc.constph builds the topology. Atoms not present in
    BOTH states are skipped (no-op).

    Computing a *delta* (instead of overwriting absolute charges) lets
    us correctly handle:
      * any prmtop whose backbone parameterization differs slightly from
        cpinutil's emitted table (cpinutil and leaprc.constph have
        historically tracked different versions of the AA charges);
      * residues whose charges were set by a non-constph build but agree
        with cpinutil per-residue.
    Step-2 already excludes terminal residues from the site list, so we
    don't have to worry about C/N-terminal-specific atoms here.
    """
    if from_state is None:
        from_state = get_residue(resname).states[0]
    delta: Dict[str, float] = {}
    for atom_name, q_target in target_state.charges.items():
        q_from = from_state.charges.get(atom_name, 0.0)
        d = q_target - q_from
        if abs(d) > 1.0e-9:
            delta[atom_name] = d
    return delta


def apply_background(structure: pmd.Structure,
                      background: BackgroundStateMap,
                      target: SiteKey) -> None:
    """Apply each background state's delta to its residue, except the target.

    Silent no-op if a key isn't found in the structure (legitimate when
    the cluster cutout excluded that site).
    """
    for (resname, resnum), state in background.items():
        if (resname, resnum) == target:
            continue
        delta = _state_delta(resname, state)
        if not delta:
            continue
        for r in structure.residues:
            if r.name == resname and (r.number + 1) == resnum:
                for a in r.atoms:
                    if a.name in delta:
                        a.charge += delta[a.name]
                break


def apply_state_charges(structure: pmd.Structure,
                         resname: str,
                         resnum_1based: int,
                         state: ProtonationState) -> float:
    """Apply the protonation-state delta to a residue's atom charges.

    Computes the per-atom Δq from STATE 0 (the constph build default) to
    `state` and *adds* the delta to the residue's existing charges. This
    preserves any backbone or non-listed atom charges already in the
    prmtop, so applying STATE 0 to a STATE-0 residue is a true no-op.

    Returns the new system net charge. Raises if the residue isn't found.
    """
    target = None
    for r in structure.residues:
        if r.name == resname and (r.number + 1) == resnum_1based:
            target = r
            break
    if target is None:
        raise ValueError(
            f"Residue {resname}-{resnum_1based} not found in structure")
    delta = _state_delta(resname, state)
    for a in target.atoms:
        if a.name in delta:
            a.charge += delta[a.name]
    return sum(a.charge for a in structure.atoms)


def _find_target_in_model(model: pmd.Structure, resname: str) -> int:
    """Return the 1-based residue number of `resname` in an ACE-X-NME model."""
    for r in model.residues:
        if r.name == resname:
            return r.number + 1
    raise ValueError(f"{resname} not found in model compound")


def compute_intrinsic_pka(
        cluster_prmtop: Path,
        cluster_rst7:   Path,
        target_resname: str,
        target_resnum_1based: int,
        *,
        residue: Optional[TitratableResidue] = None,
        model_prmtop: Optional[Path] = None,
        model_rst7:   Optional[Path] = None,
        model_cache_dir: Optional[Path] = None,
        work_dir: Path,
        prot_state: Optional[ProtonationState] = None,
        deprot_state: Optional[ProtonationState] = None,
        background: Optional[BackgroundStateMap] = None,
        cluster_radius: Optional[float] = None,
        bond_depth: int = 3,
        cluster_fillratio: float = 2.0,
        model_fillratio:   float = 4.0,
        site: Optional[Site] = None,
        **pbsa_params: Any,
) -> Dict[str, Any]:
    """4-state Bashford cycle for one titratable site.

    Parameters
    ----------
    cluster_prmtop, cluster_rst7
        The protein (or cluster cutout). Must have the target residue
        present with the constph-style atom set (all titrating Hs present).
    target_resname, target_resnum_1based
        The residue to titrate. Resname matches the topology after any
        constph rename (ASP→AS4, etc.).
    residue
        TitratableResidue object. If omitted, fetched via
        `get_residue(target_resname)`.
    model_prmtop, model_rst7
        Model compound (ACE-X-NME for amino acids; methyl-capped singleton
        for heme titratables like PRN). If omitted, built on demand via
        `model_compounds.build_model(target_resname)`.
    prot_state, deprot_state
        Specific ProtonationStates from `residue.states` to use as the
        two ends of the titration. Defaults to `residue.prot_state` and
        `residue.deprot_state` (most/least protonated). Override for
        multi-state residues like HIP where you want a specific tautomer.
    work_dir
        Directory to drop pbsa input/output files. Created if absent.
    cluster_fillratio, model_fillratio
        PB grid box-size multipliers. Model is much smaller, so we use
        a larger fillratio to keep the box-edge artifacts small.
    **pbsa_params
        Forwarded to `pb_backend.run_pbsa`. Includes epsin, istrng,
        space, nfocus, bcopt, etc. See `pb_backend.run_pbsa` for defaults.

    Returns
    -------
    dict with at least:
        pKa, dpKa, ddG_kcal,
        ddG_protein, ddG_model,
        e_cluster_PROT, e_cluster_DEPROT,
        e_model_PROT,   e_model_DEPROT,
        n_atoms_cluster, n_atoms_model
    """
    if residue is None:
        residue = get_residue(target_resname)
    if prot_state is None:
        prot_state = residue.prot_state
    if deprot_state is None:
        deprot_state = residue.deprot_state
    if model_prmtop is None or model_rst7 is None:
        model_prmtop, model_rst7 = build_model(
            target_resname, cache_dir=model_cache_dir)

    work_dir.mkdir(parents=True, exist_ok=True)
    energies: Dict[str, Dict[str, float]] = {}
    net_qs:   Dict[str, float] = {}
    target_resnum_in_model = _find_target_in_model(
        pmd.load_file(str(model_prmtop)), target_resname)

    target_key: SiteKey = (target_resname, target_resnum_1based)
    for state_label, state in [("PROT", prot_state), ("DEPROT", deprot_state)]:
        # Cluster: load full structure, optionally cut out, then apply
        # background and target charges. If a Site is provided, the cutout
        # uses the envelope-aware path (preserves multi-residue MCPB groups
        # like hemes; enforces all-or-nothing on overlapping redox sites).
        cluster = pmd.load_file(str(cluster_prmtop), str(cluster_rst7))
        if cluster_radius is not None:
            if site is not None:
                keep = cluster_for_site(cluster, site,
                                          R=cluster_radius, bond_depth=bond_depth)
            else:
                target_in_full = next(
                    r for r in cluster.residues
                    if r.name == target_resname
                    and (r.number + 1) == target_resnum_1based)
                keep = cluster_for_target(cluster, target_in_full,
                                            R=cluster_radius, bond_depth=bond_depth)
            cluster = slice_structure(cluster, keep)
        if background:
            apply_background(cluster, background, target_key)
        net_qs[f"cluster_{state_label}"] = apply_state_charges(
            cluster, target_resname, target_resnum_1based, state)
        energies[f"cluster_{state_label}"] = run_pbsa(
            cluster, work_dir, f"cluster_{state_label}",
            fillratio=cluster_fillratio, **pbsa_params)
        # Model: a single capped tripeptide; background not applicable
        model = pmd.load_file(str(model_prmtop), str(model_rst7))
        net_qs[f"model_{state_label}"] = apply_state_charges(
            model, target_resname, target_resnum_in_model, state)
        energies[f"model_{state_label}"] = run_pbsa(
            model, work_dir, f"model_{state_label}",
            fillratio=model_fillratio, **pbsa_params)

    d_protein = (total_electrostatic(energies["cluster_DEPROT"])
                  - total_electrostatic(energies["cluster_PROT"]))
    d_model   = (total_electrostatic(energies["model_DEPROT"])
                  - total_electrostatic(energies["model_PROT"]))
    ddg = d_protein - d_model
    dpka = ddg / RT_LN10
    pka = residue.pka_model + dpka

    return {
        "resname": target_resname,
        "resnum":  target_resnum_1based,
        "pKa_model": residue.pka_model,
        "pKa": pka,
        "dpKa": dpka,
        "ddG_kcal": ddg,
        "ddG_protein": d_protein,
        "ddG_model": d_model,
        "e_cluster_PROT": energies["cluster_PROT"],
        "e_cluster_DEPROT": energies["cluster_DEPROT"],
        "e_model_PROT": energies["model_PROT"],
        "e_model_DEPROT": energies["model_DEPROT"],
        "net_q_cluster_PROT": net_qs["cluster_PROT"],
        "net_q_cluster_DEPROT": net_qs["cluster_DEPROT"],
        "net_q_model_PROT": net_qs["model_PROT"],
        "net_q_model_DEPROT": net_qs["model_DEPROT"],
    }


# ===========================================================================
# Multi-state machinery (HIS 3-state, AS4/GL4 5-state with tautomer averaging)
# ===========================================================================

def compute_state_energies(
        cluster_prmtop: Path,
        cluster_rst7:   Path,
        target_resname: str,
        target_resnum_1based: int,
        *,
        residue: Optional[TitratableResidue] = None,
        model_prmtop: Optional[Path] = None,
        model_rst7:   Optional[Path] = None,
        model_cache_dir: Optional[Path] = None,
        work_dir: Path,
        background: Optional[BackgroundStateMap] = None,
        cluster_radius: Optional[float] = None,
        bond_depth: int = 3,
        cluster_fillratio: float = 2.0,
        model_fillratio:   float = 4.0,
        site: Optional[Site] = None,
        **pbsa_params: Any,
) -> Dict[int, Dict[str, Any]]:
    """Run PB on every state of a multi-state titratable residue.

    Returns a dict keyed by state index, with values:
        {prot_count, pka_corr, e_protein, e_model, ddG_PB_minus_model}

    The ddG_PB_minus_model value is (E_PB_protein_s - E_PB_model_s) in
    kcal/mol — the protein-environment shift relative to the model
    compound for state s. This is the only state-dependent piece in the
    Bashford framework; combine with `state_free_energies_at_pH` to get
    populations at a target pH.

    For 2-state residues this gives the same information as
    `compute_intrinsic_pka` (in a more general form).
    """
    if residue is None:
        residue = get_residue(target_resname)
    if model_prmtop is None or model_rst7 is None:
        model_prmtop, model_rst7 = build_model(
            target_resname, cache_dir=model_cache_dir)

    work_dir.mkdir(parents=True, exist_ok=True)
    target_resnum_in_model = _find_target_in_model(
        pmd.load_file(str(model_prmtop)), target_resname)
    target_key: SiteKey = (target_resname, target_resnum_1based)

    out: Dict[int, Dict[str, Any]] = {}
    for s_idx, state in enumerate(residue.states):
        # Cluster PB (with optional cutout). Envelope-aware path used when
        # `site` is provided (multi-residue MCPB groups, partial-overlap
        # protection); single-residue path otherwise.
        cluster = pmd.load_file(str(cluster_prmtop), str(cluster_rst7))
        if cluster_radius is not None:
            if site is not None:
                keep = cluster_for_site(cluster, site,
                                          R=cluster_radius, bond_depth=bond_depth)
            else:
                target_in_full = next(
                    r for r in cluster.residues
                    if r.name == target_resname
                    and (r.number + 1) == target_resnum_1based)
                keep = cluster_for_target(cluster, target_in_full,
                                            R=cluster_radius, bond_depth=bond_depth)
            cluster = slice_structure(cluster, keep)
        if background:
            apply_background(cluster, background, target_key)
        apply_state_charges(cluster, target_resname, target_resnum_1based, state)
        e_p = run_pbsa(cluster, work_dir, f"cluster_state{s_idx}",
                        fillratio=cluster_fillratio, **pbsa_params)
        # Model PB
        model = pmd.load_file(str(model_prmtop), str(model_rst7))
        apply_state_charges(model, target_resname, target_resnum_in_model, state)
        e_m = run_pbsa(model, work_dir, f"model_state{s_idx}",
                        fillratio=model_fillratio, **pbsa_params)

        out[s_idx] = {
            "state_name":  state.name,
            "prot_count":  state.prot_count,
            "pka_corr":    state.pka_corr,
            "e_protein":   total_electrostatic(e_p),
            "e_model":     total_electrostatic(e_m),
            "ddG_PB_minus_model":
                total_electrostatic(e_p) - total_electrostatic(e_m),
            "raw_protein": e_p,
            "raw_model":   e_m,
        }
    return out


def state_free_energies_at_pH(
        state_energies: Dict[int, Dict[str, Any]],
        pH: float) -> Dict[int, float]:
    """Effective free energy of each state at given pH (kcal/mol).

    G_eff_s(pH) = (E_PB_protein_s - E_PB_model_s)
                  + RT·ln10 · (n_s · pH − pka_corr[s])

    The constant offset (same for all states) cancels in any
    population calculation, so we report G_eff in absolute terms; only
    differences matter.

    Verified consistency with the 2-state code: for AS4 with
    states={DEPROT(n=0,pka=0), PROT(n=1,pka=4)}, the formula reduces to
    G_PROT − G_DEPROT = RT·ln10·(pH − pKa_predicted).
    """
    out: Dict[int, float] = {}
    for s_idx, info in state_energies.items():
        ph_term = RT_LN10 * (info["prot_count"] * pH - info["pka_corr"])
        out[s_idx] = info["ddG_PB_minus_model"] + ph_term
    return out


def state_populations_at_pH(
        state_energies: Dict[int, Dict[str, Any]],
        pH: float,
        T: float = 298.15) -> Dict[int, float]:
    """Boltzmann populations of each state at given pH. Sums to 1.0."""
    G = state_free_energies_at_pH(state_energies, pH)
    # Subtract the minimum to keep exponentials numerically tame.
    g_min = min(G.values())
    # kT in kcal/mol at given T:
    kT = RT_LN10 * T / (298.15 * 2.302585)
    weights = {s: exp(-(g - g_min) / kT) for s, g in G.items()}
    Z = sum(weights.values())
    return {s: w / Z for s, w in weights.items()}


def populations_by_prot_count(
        state_energies: Dict[int, Dict[str, Any]],
        pH: float,
        T: float = 298.15) -> Dict[int, float]:
    """Sum chemistry-group populations across distinct prot_counts.

    Built on top of `populations_by_chemistry_group` (the right primitive)
    so symmetry-equivalent microstates aren't double-counted.
    """
    chem_pops = populations_by_chemistry_group(state_energies, pH, T)
    by_n: Dict[int, float] = defaultdict(float)
    for (n, _pka), p in chem_pops.items():
        by_n[n] += p
    return dict(by_n)


# ---------------------------------------------------------------------------
# Chemistry grouping — collapses symmetry-equivalent tautomers (AS4's 4
# PROT states, GL4's 4 PROT states) into a single chemical group, while
# keeping genuinely distinct chemistries (HIS HID vs HIE) separate.
# ---------------------------------------------------------------------------

ChemKey = Tuple[int, float]  # (prot_count, pka_corr)


def group_states_by_chemistry(
        state_energies: Dict[int, Dict[str, Any]]
        ) -> Dict[ChemKey, List[int]]:
    """Partition state indices by chemical equivalence.

    Two states are "chemically equivalent" when they share the same
    `prot_count` AND `pka_corr`. cpinutil uses this convention: the 4
    PROT tautomers of AS4 all have (n=1, pka_corr=4.0); HIS HID(n=1,
    pka_corr=0.6) and HIE(n=1, pka_corr=0.0) are kept separate.
    """
    groups: Dict[ChemKey, List[int]] = defaultdict(list)
    for s_idx, info in state_energies.items():
        key: ChemKey = (info["prot_count"], info["pka_corr"])
        groups[key].append(s_idx)
    return dict(groups)


def chemistry_group_free_energy(
        state_energies: Dict[int, Dict[str, Any]],
        member_idxs: List[int],
        pH: float,
        T: float = 298.15) -> float:
    """Effective G of one chemistry group at given pH (kcal/mol).

    G_group = -kT · ln( (1/N_group) · Σ_{i ∈ group} exp(-G_i/kT) )

    The 1/N factor removes the degeneracy: cpinutil's pka_corr is
    calibrated per-chemistry, so summing N tautomers naively would shift
    the effective pKa by +log10(N). Dividing by N restores the calibration
    while still allowing the in-protein favored tautomer to dominate
    when the ΔΔG values across tautomers differ.
    """
    G_per_state = state_free_energies_at_pH(state_energies, pH)
    Gs = [G_per_state[i] for i in member_idxs]
    g_min = min(Gs)
    kT = (RT_LN10 / 2.302585) * (T / 298.15)  # 0.5925 kcal/mol at 298.15 K
    Z_group = sum(exp(-(g - g_min) / kT) for g in Gs) / len(member_idxs)
    return g_min - kT * log(Z_group)


def populations_by_chemistry_group(
        state_energies: Dict[int, Dict[str, Any]],
        pH: float,
        T: float = 298.15) -> Dict[ChemKey, float]:
    """Boltzmann populations of each chemistry group at given pH.

    Properly handles symmetry-equivalent tautomers (AS4 PROT, GL4 PROT)
    so the resulting populations match experimental pKa conventions.
    """
    groups = group_states_by_chemistry(state_energies)
    G_groups = {key: chemistry_group_free_energy(state_energies, members, pH, T)
                for key, members in groups.items()}
    g_min = min(G_groups.values())
    kT = (RT_LN10 / 2.302585) * (T / 298.15)
    weights = {k: exp(-(g - g_min) / kT) for k, g in G_groups.items()}
    Z = sum(weights.values())
    return {k: w / Z for k, w in weights.items()}


def effective_pKa_two_state(
        state_energies: Dict[int, Dict[str, Any]],
        pH_anchor: float = 7.0,
        T: float = 298.15) -> Optional[float]:
    """For a 2-prot-count residue (e.g., AS4 with multi-PROT tautomers
    all at n=1), return the effective pKa from the population ratio.

    Returns None if the residue has more than 2 distinct prot_count
    values (e.g., HIS where {0,1,2} are all distinct → no single pKa).
    """
    pops = populations_by_prot_count(state_energies, pH_anchor, T)
    if len(pops) != 2:
        return None
    n_low, n_high = sorted(pops.keys())
    if pops[n_high] <= 0 or pops[n_low] <= 0:
        return None
    # P(high)/P(low) = 10^((pH_anchor - pKa_eff) * (n_high - n_low))
    # → pKa_eff = pH_anchor - log10(P(high)/P(low)) / (n_high - n_low)
    # We want the standard convention: pKa for the protonated form,
    # i.e., (n_high − n_low) deprotonation events.
    # P(deprot)/P(prot) = 10^((pH - pKa) * Δn)  with Δn = (n_high - n_low) > 0
    delta_n = n_high - n_low
    return pH_anchor - log10(pops[n_low] / pops[n_high]) / delta_n


def compute_pka_multistate(
        cluster_prmtop: Path,
        cluster_rst7:   Path,
        target_resname: str,
        target_resnum_1based: int,
        *,
        pH: float = 7.0,
        residue: Optional[TitratableResidue] = None,
        model_prmtop: Optional[Path] = None,
        model_rst7:   Optional[Path] = None,
        model_cache_dir: Optional[Path] = None,
        work_dir: Path,
        background: Optional[BackgroundStateMap] = None,
        cluster_fillratio: float = 2.0,
        model_fillratio:   float = 4.0,
        site: Optional[Site] = None,
        **pbsa_params: Any,
) -> Dict[str, Any]:
    """Convenience wrapper: run multi-state PB and report populations + pKa.

    Returns a dict with:
        state_energies          (per-state dict)
        populations_per_state   (Boltzmann at given pH)
        populations_by_prot_count (aggregated)
        dominant_state_idx      (argmax over populations)
        dominant_prot_count     (prot_count of the most-populated state)
        pKa_eff                 (None for HIS-style 3-prot-count residues)
    """
    if residue is None:
        residue = get_residue(target_resname)
    se = compute_state_energies(
        cluster_prmtop, cluster_rst7,
        target_resname, target_resnum_1based,
        residue=residue,
        model_prmtop=model_prmtop, model_rst7=model_rst7,
        model_cache_dir=model_cache_dir,
        work_dir=work_dir, background=background,
        cluster_fillratio=cluster_fillratio,
        model_fillratio=model_fillratio,
        site=site,
        **pbsa_params,
    )
    pops = state_populations_at_pH(se, pH)
    pops_n = populations_by_prot_count(se, pH)
    dom_idx = max(pops, key=pops.get)
    return {
        "resname": target_resname,
        "resnum":  target_resnum_1based,
        "pH":      pH,
        "pKa_model": residue.pka_model,
        "state_energies":          se,
        "populations_per_state":   pops,
        "populations_by_prot_count": pops_n,
        "dominant_state_idx":      dom_idx,
        "dominant_prot_count":     se[dom_idx]["prot_count"],
        "pKa_eff":                 effective_pKa_two_state(se, pH),
    }
