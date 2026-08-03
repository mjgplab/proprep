"""Cluster mean-field solver.

Hybrid of mean-field and exact enumeration: build a coupling graph by
thresholding |W_ij|, take connected components as 'clusters', and treat:

  * intra-cluster terms      → exact enumeration (small clusters → fast)
  * inter-cluster terms      → mean-field over the marginals of others

Self-consistent iteration over the per-cluster marginals converges to a
fixed point. This handles frustrated pairs correctly (the anti-correlated
joint distribution lives entirely inside one cluster) while keeping the
overall cost cheap when most pairs are weakly coupled (which is typical
for protein systems).

Inputs:
  * SelfEnergies + CouplingMatrix from coupling.py (same as MC / enumerate)
  * pH
  * threshold_kT — edge threshold for cluster construction
  * max_iter, tol — iteration controls
  * max_cluster_states — safety bound on per-cluster enumeration cost

Returns:
  * marginals[site_key][chem] — Boltzmann marginal at convergence
  * dominant_chem[site_key]   — argmax of marginals (used for the state map)
  * clusters                   — the connected components found
  * n_iterations, converged    — convergence diagnostics
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from math import exp
from typing import Dict, List, Tuple

from ..coupling import CouplingMatrix, CouplingSite, SelfEnergies
from ..pb_backend import RT_LN10

SiteKey = Tuple[str, int]
ChemKey = Tuple[int, float]


@dataclass
class CMFResult:
    sites:                List[CouplingSite]
    pH:                   float
    threshold_kT:         float
    clusters:             List[List[SiteKey]]
    cluster_sizes:        List[int]
    n_iterations:         int
    converged:            bool
    marginals:            Dict[SiteKey, Dict[ChemKey, float]]
    dominant_chem:        Dict[SiteKey, ChemKey]


def _max_pair_W(W: CouplingMatrix, si: CouplingSite, sj: CouplingSite) -> float:
    """Largest |W_ij(c_i, c_j)| over the pair's non-ref chemistries."""
    max_w = 0.0
    for c_i in si.chemistries:
        if c_i == si.ref_chem:
            continue
        for c_j in sj.chemistries:
            if c_j == sj.ref_chem:
                continue
            w = W.W.get((si.key, c_i), {}).get((sj.key, c_j), 0.0)
            if abs(w) > abs(max_w):
                max_w = w
    return max_w


def build_clusters(sites: List[CouplingSite],
                    W: CouplingMatrix,
                    threshold_kcal: float) -> List[List[SiteKey]]:
    """Connected components of the coupling graph (edges where |W|>threshold)."""
    n = len(sites)
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x: int, y: int) -> None:
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for i in range(n):
        for j in range(i + 1, n):
            if abs(_max_pair_W(W, sites[i], sites[j])) > threshold_kcal:
                union(i, j)

    groups: Dict[int, List[SiteKey]] = {}
    for i, s in enumerate(sites):
        root = find(i)
        groups.setdefault(root, []).append(s.key)
    # Sort clusters by descending size for predictable output / better warm-up
    return sorted(groups.values(), key=lambda g: -len(g))


def cluster_state_count(cluster: List[CouplingSite]) -> int:
    """Total joint chemistry-state combinations for one cluster."""
    n = 1
    for s in cluster:
        n *= s.n_chemistries()
    return n


def solve(self_e: SelfEnergies,
           W:      CouplingMatrix,
           *,
           pH: float = 7.0,
           T:  float = 298.15,
           threshold_kT: float = 1.0,
           max_iter: int = 100,
           tol: float = 1.0e-5,
           max_cluster_states: int = 65536,
           ) -> CMFResult:
    """Run cluster mean-field iteration to self-consistency.

    `threshold_kT` controls the cluster-graph edge threshold. Lowering
    it merges more sites into bigger clusters (more accurate, more
    expensive). 1 kT is a sensible default — pairs whose interaction
    can flip a site's preference all stay in the same cluster.

    Raises ValueError if any cluster's joint state space exceeds
    `max_cluster_states` — that's a sign your threshold is too low for
    your system, or you have a true many-body coupling that needs MC.
    """
    sites = self_e.sites
    kT = (RT_LN10 / 2.302585) * (T / 298.15)
    threshold_kcal = threshold_kT * kT

    cluster_keys = build_clusters(sites, W, threshold_kcal)
    site_by_key = {s.key: s for s in sites}
    clusters = [[site_by_key[k] for k in keys] for keys in cluster_keys]

    cluster_sizes = [len(c) for c in clusters]
    cluster_state_counts = [cluster_state_count(c) for c in clusters]
    biggest = max(cluster_state_counts) if cluster_state_counts else 1
    if biggest > max_cluster_states:
        # Identify offending cluster
        idx = max(range(len(cluster_state_counts)),
                   key=lambda i: cluster_state_counts[i])
        raise ValueError(
            f"Cluster {idx} has {cluster_state_counts[idx]} joint states "
            f"({cluster_sizes[idx]} sites), exceeding max_cluster_states="
            f"{max_cluster_states}. Either raise the threshold (current "
            f"{threshold_kT} kT) to break it up, raise max_cluster_states, "
            f"or use Monte Carlo for this system.")

    # Per-site intrinsic free energy at pH (independent of the iteration).
    G_intr: Dict[Tuple[SiteKey, ChemKey], float] = {}
    for s in sites:
        for c in s.chemistries:
            ddg = self_e.intrinsic_ddG.get((s.key, c), 0.0)
            ph_term = RT_LN10 * (c[0] * pH - c[1])
            G_intr[(s.key, c)] = ddg + ph_term

    # Initial marginals: every site at ref_chem with probability 1.
    marginals: Dict[SiteKey, Dict[ChemKey, float]] = {
        s.key: {c: (1.0 if c == s.ref_chem else 0.0) for c in s.chemistries}
        for s in sites
    }

    converged = False
    iter_count = 0
    for it in range(max_iter):
        iter_count = it + 1
        new_marginals: Dict[SiteKey, Dict[ChemKey, float]] = {}
        max_change = 0.0

        for cluster in clusters:
            cluster_set = {s.key for s in cluster}

            # External mean field on each (site_in_cluster, chem) from sites
            # outside the cluster, weighted by their current marginals.
            ext_field: Dict[Tuple[SiteKey, ChemKey], float] = {}
            for s_in in cluster:
                for c_i in s_in.chemistries:
                    field = 0.0
                    for s_out in sites:
                        if s_out.key in cluster_set:
                            continue
                        p_out = marginals[s_out.key]
                        for c_j, p_j in p_out.items():
                            if p_j == 0.0:
                                continue
                            w = W.W.get((s_in.key, c_i), {}).get(
                                (s_out.key, c_j), 0.0)
                            field += p_j * w
                    ext_field[(s_in.key, c_i)] = field

            # Enumerate cluster's joint states.
            chem_lists = [s.chemistries for s in cluster]
            energies: List[Tuple[Tuple[ChemKey, ...], float]] = []
            for state_tuple in product(*chem_lists):
                G = 0.0
                for s_in, c_i in zip(cluster, state_tuple):
                    G += G_intr[(s_in.key, c_i)] + ext_field[(s_in.key, c_i)]
                # Intra-cluster pair coupling (exact)
                k = len(cluster)
                for ii in range(k):
                    s_i = cluster[ii]
                    c_i = state_tuple[ii]
                    for jj in range(ii + 1, k):
                        s_j = cluster[jj]
                        c_j = state_tuple[jj]
                        G += W.W.get((s_i.key, c_i), {}).get(
                            (s_j.key, c_j), 0.0)
                energies.append((state_tuple, G))

            # Boltzmann marginalization (subtract g_min for numerical stability).
            g_min = min(g for _, g in energies)
            ws = [(st, exp(-(g - g_min) / kT)) for st, g in energies]
            Z = sum(w for _, w in ws)

            for ii, s_in in enumerate(cluster):
                p_new = {c: 0.0 for c in s_in.chemistries}
                for st, w in ws:
                    p_new[st[ii]] += w / Z
                new_marginals[s_in.key] = p_new
                for c in s_in.chemistries:
                    d = abs(p_new[c] - marginals[s_in.key][c])
                    if d > max_change:
                        max_change = d

        marginals = new_marginals
        if max_change < tol:
            converged = True
            break

    dominant = {s.key: max(marginals[s.key], key=marginals[s.key].get)
                for s in sites}

    return CMFResult(
        sites=sites, pH=pH, threshold_kT=threshold_kT,
        clusters=cluster_keys, cluster_sizes=cluster_sizes,
        n_iterations=iter_count, converged=converged,
        marginals=marginals, dominant_chem=dominant,
    )
