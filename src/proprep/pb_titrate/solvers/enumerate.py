"""Exact enumeration over all chemistry-state combinations.

Gold standard when feasible (≤ ~16 total states). Computes the full
partition function exactly: no MC sampling error, no mean-field
approximation. Useful for validating MC and mean-field on small systems.

Cost: O(product of n_chemistries per site). 16 states is the practical
ceiling — beyond that, MC.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import product
from math import exp
from typing import Any, Dict, List, Tuple

from ..coupling import CouplingMatrix, CouplingSite, SelfEnergies, ChemKey
from ..intrinsic import SiteKey
from ..pb_backend import RT_LN10


@dataclass
class EnumResult:
    sites:                List[CouplingSite]
    pH:                   float
    n_states_enumerated:  int
    marginals:            Dict[SiteKey, Dict[ChemKey, float]]
    dominant_state:       List[ChemKey]
    state_probabilities:  Dict[Tuple[ChemKey, ...], float] = field(default_factory=dict)


def enumerate_all(self_e: SelfEnergies,
                   W:      CouplingMatrix,
                   *,
                   pH: float = 7.0,
                   T:  float = 298.15,
                   ) -> EnumResult:
    """Enumerate all state combinations, compute exact populations."""
    sites = self_e.sites
    kT_kcal = (RT_LN10 / 2.302585) * (T / 298.15)

    # Enumerate
    chem_lists = [s.chemistries for s in sites]
    weights: Dict[Tuple[ChemKey, ...], float] = {}
    Z = 0.0
    g_min = float("inf")
    G_per_state: Dict[Tuple[ChemKey, ...], float] = {}

    # Two-pass: first compute all G to find G_min, then weights
    for assignment in product(*chem_lists):
        g = 0.0
        for i, chem_i in enumerate(assignment):
            n_si, pka_si = chem_i
            g += (self_e.intrinsic_ddG[(sites[i].key, chem_i)]
                  + RT_LN10 * (n_si * pH - pka_si))
        for i in range(len(assignment)):
            ki = sites[i].key
            for j in range(i + 1, len(assignment)):
                kj = sites[j].key
                g += W.W.get((ki, assignment[i]), {}).get(
                    (kj, assignment[j]), 0.0)
        G_per_state[assignment] = g
        if g < g_min:
            g_min = g

    # Compute weights and partition function
    for assignment, g in G_per_state.items():
        w = exp(-(g - g_min) / kT_kcal)
        weights[assignment] = w
        Z += w

    state_probabilities = {a: w / Z for a, w in weights.items()}

    # Marginal populations per site per chemistry
    marginals: Dict[SiteKey, Dict[ChemKey, float]] = {
        s.key: {c: 0.0 for c in s.chemistries} for s in sites
    }
    for assignment, p in state_probabilities.items():
        for i, chem_i in enumerate(assignment):
            marginals[sites[i].key][chem_i] += p

    dominant_state = list(max(state_probabilities, key=state_probabilities.get))

    return EnumResult(
        sites=sites, pH=pH,
        n_states_enumerated=len(weights),
        marginals=marginals,
        dominant_state=dominant_state,
        state_probabilities=state_probabilities,
    )
