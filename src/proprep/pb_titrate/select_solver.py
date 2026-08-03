"""Auto-pick mean-field / MC / exact-enum based on coupling structure.

After the precompute step (`coupling.compute_self_energies` +
`coupling.compute_pair_couplings`), inspect the matrix and the per-site
intrinsic energies to recommend a solver.

Returns a Recommendation with the chosen solver and a human-readable
explanation. Power users override with `--solver mf|mc|enum`.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple

from .coupling import CouplingMatrix, CouplingSite, SelfEnergies
from .pb_backend import RT_LN10


# Maximum n_total (product of chemistries across all sites) for which
# exact enumeration is considered feasible. Pure-Python enumerate_all
# evaluates ~10^5 states in seconds and ~10^6 in a minute on commodity
# hardware. Above ~10^6 it becomes impractical and MC takes over.
ENUMERATE_MAX_STATES = 100_000

# Largest per-cluster joint-state count cluster_mean_field can enumerate.
# Mirrors solvers.cluster_mean_field.solve's max_cluster_states default; if
# the strongly-coupled sites form a connected cluster bigger than this, CMF
# cannot run and Monte Carlo is the right tool.
CMF_MAX_CLUSTER_STATES = 65_536


@dataclass
class Recommendation:
    solver: str          # 'mean_field' | 'monte_carlo' | 'enumerate'
    reason: str
    diagnostics: dict


def _kT(T: float = 298.15) -> float:
    return (RT_LN10 / 2.302585) * (T / 298.15)


def near_pka_sites(self_e: SelfEnergies, pH: float = 7.0,
                    tol_kT: float = 1.0) -> List[Tuple[str, int]]:
    """Sites whose intrinsic energy difference between any two chemistries is
    within `tol_kT × kT` at the system pH (i.e., ill-determined state)."""
    kT = _kT()
    out = []
    for site in self_e.sites:
        # Compare chemistries pairwise
        ddgs = []
        for chem in site.chemistries:
            n, pka = chem
            ddg = self_e.intrinsic_ddG[(site.key, chem)] + RT_LN10 * (n * pH - pka)
            ddgs.append(ddg)
        ddgs.sort()
        # If any two adjacent chemistries are within tol_kT
        for i in range(len(ddgs) - 1):
            if abs(ddgs[i+1] - ddgs[i]) < tol_kT * kT:
                out.append((site.resname, site.resnum))
                break
    return out


def max_coupling(W: CouplingMatrix) -> float:
    """Return the largest |W_ij| in the coupling matrix (kcal/mol)."""
    max_w = 0.0
    for site_i in W.sites:
        for chem_i in site_i.chemistries:
            if chem_i == site_i.ref_chem: continue
            for site_j in W.sites:
                if site_j is site_i: continue
                for chem_j in site_j.chemistries:
                    if chem_j == site_j.ref_chem: continue
                    w = W.W.get((site_i.key, chem_i), {}).get(
                        (site_j.key, chem_j), 0.0)
                    if abs(w) > abs(max_w):
                        max_w = w
    return max_w


def total_states(sites: List[CouplingSite]) -> int:
    """Total chemistry-state combinations (product over sites)."""
    n = 1
    for s in sites:
        n *= s.n_chemistries()
    return n


def frustrated_pairs(self_e: SelfEnergies, W: CouplingMatrix,
                       pH: float = 7.0) -> List[Tuple[str, int, str, int, float]]:
    """Identify (i, j) pairs where |W_ij| exceeds the smaller of the two
    sites' individual preferences for their lowest-G chemistry.

    Such pairs are 'frustrated': individually each site has a clear
    preferred chemistry, but the coupling penalty for both being there
    simultaneously is large enough to flip the weaker-preference site.
    Mean-field cannot represent the resulting anti-correlated joint
    distribution (it either oscillates or converges to a wrong
    asymmetric fixed point depending on the iteration order).

    Returns a list of (resname_i, resnum_i, resname_j, resnum_j, W) tuples,
    sorted by descending |W|.
    """
    # Per-site G_eff for each chemistry at the target pH.
    G_eff: Dict[Tuple[Tuple[str, int], Tuple[int, float]], float] = {}
    for site in self_e.sites:
        for chem in site.chemistries:
            ddg = self_e.intrinsic_ddG.get((site.key, chem), 0.0)
            ph_term = RT_LN10 * (chem[0] * pH - chem[1])
            G_eff[(site.key, chem)] = ddg + ph_term

    # Per-site best (lowest-G) chemistry and its margin over the runner-up.
    best_chem: Dict[Tuple[str, int], Tuple[int, float]] = {}
    pref_margin: Dict[Tuple[str, int], float] = {}
    for site in self_e.sites:
        gs = sorted([(G_eff[(site.key, c)], c) for c in site.chemistries])
        best_chem[site.key] = gs[0][1]
        pref_margin[site.key] = (gs[1][0] - gs[0][0]) if len(gs) > 1 else float("inf")

    out: List[Tuple[str, int, str, int, float]] = []
    sites = self_e.sites
    for i, si in enumerate(sites):
        for sj in sites[i+1:]:
            ci = best_chem[si.key]
            cj = best_chem[sj.key]
            w = W.W.get((si.key, ci), {}).get((sj.key, cj), 0.0)
            margin_i = pref_margin[si.key]
            margin_j = pref_margin[sj.key]
            # Frustrated if |W| can flip the weaker site away from its preference
            if abs(w) > min(margin_i, margin_j):
                out.append((si.resname, si.resnum,
                            sj.resname, sj.resnum, w))
    out.sort(key=lambda r: -abs(r[4]))
    return out


def largest_cluster_at_threshold(sites: List[CouplingSite], W: CouplingMatrix,
                                   threshold_kT: float, T: float = 298.15
                                   ) -> Tuple[int, int]:
    """Threshold the coupling graph at `threshold_kT` and return the
    (n_sites, n_joint_states) of the cluster with the MOST joint states —
    the one that governs whether cluster mean-field can enumerate.

    This is the cheap pre-flight that lets the recommender avoid suggesting
    cluster_mean_field for a system whose strongly-coupled sites percolate
    into one giant component (the failure the user hit: a 43-site cluster
    with ~1e14 states).
    """
    from .solvers.cluster_mean_field import build_clusters, cluster_state_count
    kT = _kT(T)
    site_by_key = {s.key: s for s in sites}
    clusters = build_clusters(sites, W, threshold_kT * kT)
    best_sites, best_states = 0, 1
    for keys in clusters:
        clu = [site_by_key[k] for k in keys]
        st = cluster_state_count(clu)
        if st > best_states:
            best_states, best_sites = st, len(clu)
    return best_sites, best_states


def cluster_feasibility_table(sites: List[CouplingSite], W: CouplingMatrix,
                               thresholds_kT: Tuple[float, ...] = (1.0, 2.0, 3.0, 4.0),
                               max_cluster_states: int = CMF_MAX_CLUSTER_STATES,
                               T: float = 298.15) -> List[dict]:
    """For each candidate threshold, report the largest cluster size and
    whether cluster mean-field could enumerate it. Used to show the user the
    threshold↔feasibility tradeoff BEFORE running (or after a CMF failure)."""
    rows = []
    for thr in thresholds_kT:
        n_sites, n_states = largest_cluster_at_threshold(sites, W, thr, T)
        rows.append({
            "threshold_kT": thr,
            "largest_cluster_sites": n_sites,
            "largest_joint_states": n_states,
            "feasible": n_states <= max_cluster_states,
        })
    return rows


def recommend(self_e: SelfEnergies, W: CouplingMatrix,
               pH: float = 7.0) -> Recommendation:
    """Return solver recommendation.

    Decision tree:
      1. n_total ≤ ENUMERATE_MAX_STATES → enumerate (gold-standard, no
                                   approximation; strictly better than MC
                                   when feasible — no acceptance ratio,
                                   no convergence question, no rare-state
                                   sampling noise)
      2. frustrated pair exists → cluster_mean_field (exact within
                                   strongly-coupled clusters, MF between)
      3. |max W| > 2 kT AND near-pKa site ≥ 1 → monte_carlo (coupling can
                                   flip an ambiguous site)
      4. ≥ 5 near-pKa sites    → monte_carlo (many ambiguous sites)
      5. otherwise             → mean_field (MF marginals are reliable)

    The reason string explains which branch was taken. mean_field can be
    safe even with strong |W| as long as no site is ambiguous (near-pKa)
    and no pair is frustrated — coupling that doesn't flip a preference
    is correctly absorbed into MF's effective field.
    """
    kT = _kT()
    max_w = max_coupling(W)
    near_pka = near_pka_sites(self_e, pH)
    n_total = total_states(self_e.sites)
    frustrated = frustrated_pairs(self_e, W, pH)
    diag = {
        "max_|W_ij|_kcal":  max_w,
        "max_|W_ij|_kT":    max_w / kT,
        "n_near_pka":       len(near_pka),
        "near_pka_sites":   near_pka,
        "n_frustrated":     len(frustrated),
        "frustrated_pairs": frustrated[:5],   # cap for readability
        "n_total_states":   n_total,
        "n_sites":          len(self_e.sites),
        "enumerate_max_states": ENUMERATE_MAX_STATES,
    }

    if n_total <= ENUMERATE_MAX_STATES:
        return Recommendation(
            "enumerate",
            f"{n_total:,} total joint states ≤ {ENUMERATE_MAX_STATES:,} "
            f"feasibility limit — exact enumeration gives the gold-standard "
            f"Boltzmann distribution with no sampling noise or convergence "
            f"question. Strictly better than MC when feasible.",
            diag)

    if frustrated:
        top = frustrated[0]
        # Feasibility pre-flight: cluster mean-field can only run if the
        # strongly-coupled sites don't percolate into a cluster too big to
        # enumerate. Check the default 1 kT threshold; if that's infeasible,
        # probe higher thresholds for one that works — but only recommend CMF
        # if a feasible threshold exists WITHOUT demoting so much coupling
        # that we're back to plain mean-field. Otherwise recommend MC.
        clu_sites, clu_states = largest_cluster_at_threshold(
            self_e.sites, W, 1.0)
        diag["cmf_largest_cluster_sites_at_1kT"] = clu_sites
        diag["cmf_largest_cluster_states_at_1kT"] = clu_states
        diag["cmf_feasibility_table"] = cluster_feasibility_table(
            self_e.sites, W)
        if clu_states <= CMF_MAX_CLUSTER_STATES:
            return Recommendation(
                "cluster_mean_field",
                f"{len(frustrated)} frustrated pair(s) detected — coupling "
                f"strong enough to flip a site away from its individual "
                f"preference, e.g., {top[0]}-{top[1]} × {top[2]}-{top[3]} "
                f"(W={top[4]:+.2f} kcal/mol = {top[4]/kT:+.1f} kT). The "
                f"strongly-coupled sites form clusters small enough to solve "
                f"exactly (largest = {clu_sites} sites, {clu_states:,} joint "
                f"states at 1 kT). Cluster mean-field handles these exactly "
                f"within each connected component while keeping inter-cluster "
                f"terms at MF level — reproducible, no sampling noise.",
                diag)
        # CMF infeasible at the natural threshold → Monte Carlo.
        feas = [r for r in diag["cmf_feasibility_table"] if r["feasible"]]
        feas_note = (
            f"A feasible threshold exists ({feas[0]['threshold_kT']:.0f} kT → "
            f"{feas[0]['largest_cluster_sites']} sites), but raising the "
            f"threshold that far demotes real couplings back to mean-field "
            f"level — defeating the purpose. "
            if feas else
            "No threshold breaks it into enumerable clusters without "
            "discarding the very couplings that matter. ")
        return Recommendation(
            "monte_carlo",
            f"{len(frustrated)} frustrated pair(s) detected (e.g. "
            f"{top[0]}-{top[1]} × {top[2]}-{top[3]}, {top[4]/kT:+.1f} kT), but "
            f"the strongly-coupled sites percolate into a {clu_sites}-site "
            f"cluster (~{clu_states:.1e} joint states at 1 kT) — far beyond "
            f"cluster mean-field's enumeration limit ({CMF_MAX_CLUSTER_STATES:,}). "
            f"{feas_note}Monte Carlo samples this joint distribution directly. "
            f"Use the built-in convergence check to confirm the result is "
            f"stable.",
            diag)

    if abs(max_w) > 2 * kT and len(near_pka) >= 1:
        return Recommendation(
            "monte_carlo",
            f"Strong coupling (max |W|={max_w:+.2f} kcal/mol = "
            f"{max_w/kT:+.1f} kT > 2 kT) with {len(near_pka)} near-pKa "
            f"site(s) — coupling could flip an ambiguous site; MF "
            f"marginals would miss the correlation.",
            diag)

    if len(near_pka) >= 5:
        return Recommendation(
            "monte_carlo",
            f"{len(near_pka)} sites within 1 kT of pKa — many ambiguous "
            f"protonation states; MC samples them properly.",
            diag)

    if abs(max_w) > 2 * kT:
        reason = (f"max |W|={max_w/kT:+.1f} kT is large, but no site is "
                  f"near-pKa and no pair is frustrated — every site has a "
                  f"clear intrinsic preference that coupling cannot flip, "
                  f"so MF will assign each site's dominant chemistry "
                  f"correctly. (Joint correlations between strongly-coupled "
                  f"sites are not captured; if you need pair populations, "
                  f"run MC.)")
    else:
        reason = (f"max |W|={max_w/kT:+.1f} kT (weak), no near-pKa sites, "
                  f"no frustrated pairs — MF marginals are accurate and "
                  f"this is the cheapest path.")
    return Recommendation("mean_field", reason, diag)
