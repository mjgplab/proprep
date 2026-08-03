"""Metropolis Monte Carlo over chemistry-group state assignments.

Operates on the precomputed (SelfEnergies, CouplingMatrix) tuple from
`titrate.coupling`. No PB calls during MC — this is integer arithmetic
on the precomputed matrices, so 10⁵–10⁶ steps run in seconds.

The total energy of a state assignment {s_i} at pH is

    G_total({s_i}, pH) = Σ_i [intrinsic_ΔΔG_i(s_i)
                              + RT·ln10 · (n_{s_i} · pH − pka_corr_{s_i})]
                        + Σ_{i<j} W_ij(s_i, s_j)

Returns marginal populations per site per chemistry, plus joint
distributions for any subset of sites.

For weakly-coupled systems mean-field will give the same answer
faster — use this when:
  - Coupling-matrix max |W_ij| > 2 kT for any pair near pKa
  - Mean-field doesn't converge (oscillates between states)
  - You want explicit joint distributions, not just marginals
"""
from __future__ import annotations

import random
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from math import exp, log
from typing import Any, Dict, List, Optional, Tuple

from ..coupling import (
    CouplingMatrix, CouplingSite, SelfEnergies, ChemKey,
)
from ..intrinsic import SiteKey
from ..pb_backend import RT_LN10


# State vector: list of chemistries, indexed by site position in the SelfEnergies.sites list
StateVector = Tuple[ChemKey, ...]


def site_intrinsic_G(self_e: SelfEnergies, site_idx: int,
                      chem: ChemKey, pH: float) -> float:
    """Intrinsic G of one site in chemistry `chem` at pH (kcal/mol).

        G_intr_i(s_i, pH) = ΔΔG_PB[s_i] + RT·ln10·(n_{s_i}·pH − pka_corr_{s_i})

    The chemistry tuple itself carries (n, pka_corr).
    """
    site = self_e.sites[site_idx]
    n_si, pka_si = chem
    ddg = self_e.intrinsic_ddG[(site.key, chem)]
    return ddg + RT_LN10 * (n_si * pH - pka_si)


def total_energy(self_e: SelfEnergies, W: CouplingMatrix,
                  state: StateVector, pH: float) -> float:
    """Sum of intrinsic + pairwise terms for one state assignment."""
    g = 0.0
    for i, chem_i in enumerate(state):
        g += site_intrinsic_G(self_e, i, chem_i, pH)
    for i in range(len(state)):
        ki = self_e.sites[i].key
        for j in range(i + 1, len(state)):
            kj = self_e.sites[j].key
            g += W.W.get((ki, state[i]), {}).get((kj, state[j]), 0.0)
    return g


def delta_energy_flip(self_e: SelfEnergies, W: CouplingMatrix,
                       state: List[ChemKey], pH: float,
                       site_idx: int, new_chem: ChemKey) -> float:
    """ΔG for flipping site `site_idx` from current chem to `new_chem`.

    Computed by direct difference (only terms involving site i change).
    O(N) per move, vs O(N²) for recomputing total_energy.
    """
    old_chem = state[site_idx]
    if old_chem == new_chem:
        return 0.0
    site = self_e.sites[site_idx]
    ki = site.key

    # Intrinsic difference at site i
    dG = (site_intrinsic_G(self_e, site_idx, new_chem, pH)
          - site_intrinsic_G(self_e, site_idx, old_chem, pH))

    # Coupling to all other sites
    for j, chem_j in enumerate(state):
        if j == site_idx:
            continue
        kj = self_e.sites[j].key
        w_old = W.W.get((ki, old_chem), {}).get((kj, chem_j), 0.0)
        w_new = W.W.get((ki, new_chem), {}).get((kj, chem_j), 0.0)
        dG += (w_new - w_old)
    return dG


# ---------------------------------------------------------------------------
# Sampler
# ---------------------------------------------------------------------------

@dataclass
class MCResult:
    sites:                    List[CouplingSite]
    pH:                       float
    n_steps_total:            int
    n_steps_equilibration:    int
    acceptance_ratio:         float
    energy_trace:             List[float]
    marginals:                Dict[SiteKey, Dict[ChemKey, float]]
    final_state:              List[ChemKey]
    # pair_counts[(site_i_key, site_j_key)] = number of post-equilibration
    # steps where BOTH sites were in their most-protonated chemistry.
    # Populated only when run_mcmc(..., record_pair_counts=True). Used by
    # `analysis.compute_susceptibility`.
    pair_counts:              Dict[Tuple[SiteKey, SiteKey], int] = field(default_factory=dict)


def run_mcmc(self_e: SelfEnergies,
              W:      CouplingMatrix,
              *,
              pH:                    float = 7.0,
              n_steps:               int   = 100_000,
              n_equil:               int   = 10_000,
              T:                     float = 298.15,
              initial_state:         Optional[List[ChemKey]] = None,
              record_energy_stride:  int   = 100,
              record_pair_counts:    bool  = False,
              seed:                  int   = 42,
              ) -> MCResult:
    """Run Metropolis MC and return marginals + trace.

    `record_pair_counts=True` accumulates joint occurrence counts
    (both sites in their most-protonated chemistry) for susceptibility
    analysis. Adds O(N²) per recorded step but enables `χ_ij`.
    """
    rng = random.Random(seed)
    sites = self_e.sites
    n_sites = len(sites)
    # Precompute most-protonated chemistry per site (for pair-count recording)
    most_prot_chems = {s.key: max(s.chemistries, key=lambda c: c[0])
                       for s in sites}

    if initial_state is None:
        state: List[ChemKey] = [s.ref_chem for s in sites]
    else:
        state = list(initial_state)

    kT_kcal = (RT_LN10 / 2.302585) * (T / 298.15)
    n_accept = 0
    n_propose = 0
    energy_trace: List[float] = []
    chem_counts: Dict[SiteKey, Counter] = {s.key: Counter() for s in sites}
    pair_counts: Dict[Tuple[SiteKey, SiteKey], int] = (
        defaultdict(int) if record_pair_counts else None)

    g_current = total_energy(self_e, W, tuple(state), pH)

    for step in range(n_steps):
        i = rng.randrange(n_sites)
        site_i = sites[i]
        if site_i.n_chemistries() <= 1:
            continue
        old_chem = state[i]
        new_chem = old_chem
        while new_chem == old_chem:
            new_chem = rng.choice(site_i.chemistries)
        n_propose += 1

        dG = delta_energy_flip(self_e, W, state, pH, i, new_chem)
        if dG <= 0 or rng.random() < exp(-dG / kT_kcal):
            state[i] = new_chem
            g_current += dG
            n_accept += 1

        if step >= n_equil:
            for k, chem in enumerate(state):
                chem_counts[sites[k].key][chem] += 1
            if record_pair_counts:
                # Record pairs where BOTH sites are in their most-protonated chemistry
                for k1 in range(n_sites):
                    if state[k1] != most_prot_chems[sites[k1].key]:
                        continue
                    for k2 in range(n_sites):
                        if state[k2] == most_prot_chems[sites[k2].key]:
                            pair_counts[(sites[k1].key, sites[k2].key)] += 1
        if step % record_energy_stride == 0:
            energy_trace.append(g_current)

    n_post_equil = max(1, n_steps - n_equil)
    marginals = {
        sk: {chem: cnt / n_post_equil for chem, cnt in counts.items()}
        for sk, counts in chem_counts.items()
    }
    for site in sites:
        for chem in site.chemistries:
            marginals[site.key].setdefault(chem, 0.0)

    return MCResult(
        sites=sites, pH=pH,
        n_steps_total=n_steps,
        n_steps_equilibration=n_equil,
        acceptance_ratio=n_accept / max(1, n_propose),
        energy_trace=energy_trace,
        marginals=marginals,
        final_state=state,
        pair_counts=dict(pair_counts) if pair_counts else {},
    )


def assign_states_from_mc(result: MCResult) -> Dict[SiteKey, ChemKey]:
    """Pick the most-populated chemistry per site (for production prmtop)."""
    out = {}
    for sk, marg in result.marginals.items():
        out[sk] = max(marg, key=marg.get)
    return out


@dataclass
class MCConvergence:
    """Result of a multi-replicate MC convergence check."""
    n_replicates:     int
    seeds:            List[int]
    n_steps:          int
    acceptance_ratios: List[float]
    # Per-site: list of dominant chemistries across replicates (one per seed).
    per_site_dominant: Dict[SiteKey, List[ChemKey]]
    # Per-site fraction of replicates that agree with the modal dominant state.
    per_site_agreement: Dict[SiteKey, float]
    # Sites where replicates DISAGREE on the dominant state (unconverged).
    unstable_sites:    List[SiteKey]
    # Consensus dominant chemistry per site (modal across replicates).
    consensus:         Dict[SiteKey, ChemKey]
    # Mean dominant-state marginal across replicates, per site (confidence).
    mean_marginal:     Dict[SiteKey, float]
    all_agree:         bool


def check_convergence(self_e: SelfEnergies,
                       W:      CouplingMatrix,
                       *,
                       pH:        float = 7.0,
                       n_steps:   int   = 1_000_000,
                       n_equil:   Optional[int] = None,
                       n_replicates: int = 4,
                       T:         float = 298.15,
                       base_seed: int   = 1,
                       on_progress: Optional[Any] = None,
                       ) -> MCConvergence:
    """Run `n_replicates` independent MC chains (different seeds) and report
    whether they agree on the dominant chemistry per site.

    This is the automated version of the manual "run it again longer and
    compare" test: if independent chains land on the same dominant state for
    every site, the result is converged; sites where they disagree are
    genuinely undersampled or near-degenerate and should be flagged (and
    typically hand-resolved in the review step).

    MC here is pure arithmetic on the precomputed matrices (no PB calls), so
    several replicates of 1e6 steps run in seconds–minutes.
    """
    if n_equil is None:
        n_equil = n_steps // 10
    seeds = [base_seed + r for r in range(n_replicates)]
    sites = self_e.sites

    per_site_dominant: Dict[SiteKey, List[ChemKey]] = {s.key: [] for s in sites}
    per_site_marginals: Dict[SiteKey, List[float]] = {s.key: [] for s in sites}
    acceptance_ratios: List[float] = []

    for idx, seed in enumerate(seeds):
        res = run_mcmc(self_e, W, pH=pH, n_steps=n_steps, n_equil=n_equil,
                        T=T, seed=seed)
        acceptance_ratios.append(res.acceptance_ratio)
        for s in sites:
            marg = res.marginals[s.key]
            dom = max(marg, key=marg.get)
            per_site_dominant[s.key].append(dom)
            per_site_marginals[s.key].append(marg[dom])
        if on_progress is not None:
            on_progress(idx + 1, len(seeds), seed, res.acceptance_ratio)

    per_site_agreement: Dict[SiteKey, float] = {}
    consensus: Dict[SiteKey, ChemKey] = {}
    mean_marginal: Dict[SiteKey, float] = {}
    unstable: List[SiteKey] = []
    for s in sites:
        doms = per_site_dominant[s.key]
        modal = Counter(doms).most_common(1)[0][0]
        agree = sum(1 for d in doms if d == modal) / len(doms)
        per_site_agreement[s.key] = agree
        consensus[s.key] = modal
        mm = per_site_marginals[s.key]
        mean_marginal[s.key] = sum(mm) / len(mm)
        if agree < 1.0:
            unstable.append(s.key)

    return MCConvergence(
        n_replicates=n_replicates, seeds=seeds, n_steps=n_steps,
        acceptance_ratios=acceptance_ratios,
        per_site_dominant=per_site_dominant,
        per_site_agreement=per_site_agreement,
        unstable_sites=unstable,
        consensus=consensus,
        mean_marginal=mean_marginal,
        all_agree=(len(unstable) == 0),
    )


def chemistry_label(site: CouplingSite, chem: ChemKey) -> str:
    """Short label for a chemistry (HIP/HID/HIE for HIS, PROT/DEPROT otherwise)."""
    n, pka = chem
    if site.residue.name == "HIP":
        if n == 2: return "HIP"
        if pka > 0.3: return "HID"
        return "HIE"
    n_max = max(c[0] for c in site.chemistries)
    n_min = min(c[0] for c in site.chemistries)
    if n == n_max: return "PROT"
    if n == n_min: return "DEPROT"
    return f"N={n}"
