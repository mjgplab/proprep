"""Coupling-aware ("effective") pKa via a coupled pH titration.

The single-site PB pKa (intrinsic) is computed with every other titratable
site held at its reference state. The *effective* pKa accounts for the
electrostatic coupling between sites: a charged neighbour shifts a site's real
titration midpoint. We obtain it by titrating the whole coupled system — running
the Monte-Carlo solver on the precomputed (SelfEnergies, CouplingMatrix) across
a pH grid and finding, per site, the pH where its protonation fraction crosses
0.5.

This is pure arithmetic on the precomputed coupling matrix — NO PB calls — so a
~20-point pH sweep runs in seconds to a couple of minutes.

Protonation fraction f_prot(pH) of a site = the summed marginal probability of
the chemistries at that site's MAXIMUM proton count (the "protonated extreme").
For two-state acids/bases that's P(protonated); for HIS (HIP n=2 vs HID/HIE n=1)
it's P(HIP), i.e. the His pKa. f_prot decreases as pH rises; the effective pKa is
the first downward 0.5 crossing, linearly interpolated.

Sites whose fraction never crosses 0.5 in the grid are "locked" (their state is
pH-independent over the sampled range); they return a string sentinel ('> phmax'
or '< phmin') so callers can exclude them from numeric statistics.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple, Union

from .coupling import CouplingMatrix, CouplingSite, SelfEnergies, ChemKey, SiteKey
from .solvers.monte_carlo import run_mcmc

# A per-site effective pKa is a float when the titration crosses 0.5 in-range,
# else a string sentinel describing where it lies relative to the grid.
EffPka = Union[float, str]


def _protonated_chems(site: CouplingSite) -> List[ChemKey]:
    """Chemistries at the site's maximum proton count (the protonated extreme).

    f_prot is the summed marginal over these. For AS4/GL4/LYS/CYS/TYR this is
    the single protonated chemistry; for HIS it is HIP (n=2) vs the neutral
    tautomers (n=1).
    """
    if not site.chemistries:
        return []
    max_n = max(n for (n, _pka) in site.chemistries)
    return [c for c in site.chemistries if c[0] == max_n]


def _f_prot(site: CouplingSite, marginals: Dict[ChemKey, float]) -> float:
    """Protonation fraction for ``site`` from its marginal distribution."""
    prot = set(_protonated_chems(site))
    return sum(p for chem, p in marginals.items() if chem in prot)


def _half_crossing(ph_grid: List[float], fracs: List[float]) -> EffPka:
    """First downward crossing of f_prot = 0.5, linearly interpolated.

    ``fracs`` is f_prot at each pH in ``ph_grid`` (ascending pH). f_prot should
    decrease with pH; we scan from low pH for the first point at/below 0.5 and
    interpolate against the previous point. Robust to mild MC non-monotonicity
    because we take the FIRST crossing from the protonated (low-pH) side.
    """
    if not ph_grid:
        return "n/a"
    # Already below 0.5 at the lowest pH → midpoint is below the grid.
    if fracs[0] < 0.5:
        return f"< {ph_grid[0]:.1f}"
    for i in range(1, len(ph_grid)):
        if fracs[i] < 0.5 <= fracs[i - 1]:
            # linear interpolation between (i-1) and i for f_prot == 0.5
            f0, f1 = fracs[i - 1], fracs[i]
            p0, p1 = ph_grid[i - 1], ph_grid[i]
            if f0 == f1:                      # flat (shouldn't happen given the
                return p0                     # bracket test, but guard /0)
            frac = (f0 - 0.5) / (f0 - f1)
            return p0 + frac * (p1 - p0)
    # Never dropped below 0.5 → midpoint is above the grid.
    return f"> {ph_grid[-1]:.1f}"


def coupled_pka_titration(
        self_e: SelfEnergies,
        W: CouplingMatrix,
        *,
        ph_min: float = 2.0,
        ph_max: float = 12.0,
        ph_step: float = 0.5,
        n_steps: int = 200_000,
        n_equil: int = 20_000,
        seed: int = 1,
        on_progress: Optional[Any] = None,
) -> Tuple[Dict[SiteKey, EffPka], Dict[SiteKey, List[Tuple[float, float]]]]:
    """Titrate the coupled system across a pH grid; return effective pKa per site.

    Returns ``(pka, curves)``:
      - ``pka[site.key]`` is the effective pKa (float) or a sentinel string
        ('< {ph_min}' / '> {ph_max}' / 'n/a') for sites that don't cross 0.5.
      - ``curves[site.key]`` is the list of ``(pH, f_prot)`` points (for CSV /
        diagnostics).

    No PB calls — MC runs on the precomputed matrices. ``on_progress(done,
    total, pH)`` is called per pH point if supplied.
    """
    # Build the ascending pH grid (inclusive of ph_max within float tolerance).
    ph_grid: List[float] = []
    pH = ph_min
    while pH <= ph_max + 1e-9:
        ph_grid.append(round(pH, 3))
        pH += ph_step
    n_pts = len(ph_grid)

    sites = self_e.sites
    # frac_by_site[key] accumulates f_prot at each pH (parallel to ph_grid).
    frac_by_site: Dict[SiteKey, List[float]] = {s.key: [] for s in sites}

    for idx, pH in enumerate(ph_grid, 1):
        res = run_mcmc(self_e, W, pH=pH, n_steps=n_steps, n_equil=n_equil,
                       seed=seed)
        for s in sites:
            frac_by_site[s.key].append(_f_prot(s, res.marginals[s.key]))
        if on_progress is not None:
            on_progress(idx, n_pts, pH)

    pka: Dict[SiteKey, EffPka] = {}
    curves: Dict[SiteKey, List[Tuple[float, float]]] = {}
    for s in sites:
        fracs = frac_by_site[s.key]
        pka[s.key] = _half_crossing(ph_grid, fracs)
        curves[s.key] = list(zip(ph_grid, fracs))
    return pka, curves
