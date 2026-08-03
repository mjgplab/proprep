"""Post-processing analyses for the titrate pipeline.

CSV outputs only — no plotting dependency. The user can plot any of the
CSVs in their tool of choice (matplotlib, gnuplot, R, Excel).

Adapted from `pka_monte_carlo/src/analysis.py`. Provides:

- write_marginals_csv          MC marginals → CSV
- write_effective_pkas_csv     per-site effective pKa from any solver
- titration_scan               multi-pH MC sweep
- henderson_hasselbalch_fit    H-H fit gives (pKa, Hill coefficient)
- write_titration_csv          titration curves + per-site pKa+Hill
- write_coupling_matrix_csv    W_ij(s,t) as a flat CSV
- write_significant_pairs_csv  filtered to |W| > threshold (network-style)
- compute_susceptibility       χ_ij(pH) from MC marginals + joint counts
- write_susceptibility_csv     per-pair susceptibility CSV
- compare_state_maps           diff two state-map CSVs

The susceptibility analysis requires MC to record pair counts; see the
extension in `solvers.monte_carlo.run_mcmc(record_pair_counts=True)`.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from math import log10
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .coupling import CouplingMatrix, CouplingSite, ChemKey
from .intrinsic import SiteKey
from .pb_backend import RT_LN10
from .solvers.monte_carlo import MCResult, chemistry_label


# ---------------------------------------------------------------------------
# Tier 2: marginals + effective pKa
# ---------------------------------------------------------------------------

def write_marginals_csv(result: MCResult, path: Path) -> None:
    """Per-site, per-chemistry MC marginal occupation probabilities."""
    with Path(path).open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["resname", "resnum", "chemistry_label",
                     "prot_count", "pka_corr", "probability"])
        for site in result.sites:
            for chem, p in sorted(result.marginals[site.key].items()):
                w.writerow([site.resname, site.resnum,
                             chemistry_label(site, chem),
                             chem[0], chem[1], f"{p:.6f}"])


def effective_pka_from_marginal(site: CouplingSite,
                                  marginal: Dict[ChemKey, float],
                                  pH: float) -> Optional[float]:
    """Per-site effective pKa from MC populations at one pH.

    For 2-chemistry-per-prot_count residues:
        pKa_eff = pH - log10(P(deprot) / P(prot))

    Returns None if the residue has more than 2 distinct prot_counts
    (HIS — use `histidine_state_fractions` for those).
    """
    by_n: Dict[int, float] = defaultdict(float)
    for chem, p in marginal.items():
        by_n[chem[0]] += p
    if len(by_n) != 2:
        return None
    n_low, n_high = sorted(by_n.keys())
    p_low, p_high = by_n[n_low], by_n[n_high]
    if p_low <= 0 or p_high <= 0:
        return None
    delta_n = n_high - n_low
    return pH - log10(p_low / p_high) / delta_n


def write_effective_pkas_csv(result: MCResult, path: Path) -> None:
    """One row per site: effective pKa at the MC's pH (where defined)."""
    with Path(path).open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["resname", "resnum", "n_chemistries",
                     "dominant_chemistry_label",
                     "P_dominant", "pKa_effective"])
        for site in result.sites:
            marg = result.marginals[site.key]
            dom = max(marg, key=marg.get)
            pka = effective_pka_from_marginal(site, marg, result.pH)
            w.writerow([site.resname, site.resnum, len(marg),
                         chemistry_label(site, dom),
                         f"{marg[dom]:.4f}",
                         "" if pka is None else f"{pka:.3f}"])


# ---------------------------------------------------------------------------
# Titration curves: multi-pH MC sweep + H-H fit
# ---------------------------------------------------------------------------

def titration_scan(self_e, W: CouplingMatrix,
                    pH_grid: List[float],
                    *,
                    n_steps: int = 50_000,
                    n_equil: int = 5_000,
                    seed: int = 42,
                    record_pair_counts: bool = False,
                    ) -> Dict[float, MCResult]:
    """Run MC at every pH in `pH_grid` (sequentially); return per-pH MCResult."""
    from .solvers.monte_carlo import run_mcmc  # avoid circular at module load
    out: Dict[float, MCResult] = {}
    print(f"  titration_scan over {len(pH_grid)} pH values...")
    for ph in pH_grid:
        out[ph] = run_mcmc(self_e, W, pH=ph, n_steps=n_steps,
                            n_equil=n_equil, seed=seed,
                            record_pair_counts=record_pair_counts)
        print(f"    pH={ph:>4.1f} done (acceptance {out[ph].acceptance_ratio*100:.1f}%)")
    return out


def henderson_hasselbalch_fit(pH_array: List[float],
                                fraction_protonated: List[float],
                                ) -> Tuple[Optional[float], Optional[float]]:
    """Fit f(pH) = 1 / (1 + 10**(n*(pH − pKa))) for (pKa, Hill n).

    Hill coefficient n = 1 means independent (no cooperativity). n > 1
    means cooperative (sharper transition than H-H). n < 1 means
    anti-cooperative (broader transition).

    Returns (pKa, n). If scipy isn't available or fit fails, returns
    (None, None) and the caller should fall back to bisection or pick
    the data point closest to f=0.5.
    """
    try:
        from scipy.optimize import curve_fit
    except ImportError:
        return None, None
    if len(pH_array) < 4:
        return None, None
    def model(pH, pKa, n):
        # pH and pKa as numpy arrays; n is scalar
        from math import e
        # Use 1/(1+10^(n(pH-pKa)))
        return 1.0 / (1.0 + 10.0 ** (n * (pH - pKa)))
    try:
        popt, _ = curve_fit(model, pH_array, fraction_protonated,
                             p0=[7.0, 1.0], maxfev=2000)
        return float(popt[0]), float(popt[1])
    except Exception:
        return None, None


def write_titration_csv(scan: Dict[float, MCResult], path: Path,
                          fit_path: Optional[Path] = None) -> None:
    """Wide CSV: rows = pH, columns = per-site fraction-protonated.

    Optionally also write a per-site H-H fit summary CSV.
    """
    pH_list = sorted(scan.keys())
    sites = scan[pH_list[0]].sites
    # Build per-site fraction-protonated time series
    series: Dict[SiteKey, List[float]] = {}
    for site in sites:
        f_protonated = []
        n_max = max(c[0] for c in site.chemistries)
        for ph in pH_list:
            marg = scan[ph].marginals[site.key]
            p_top = sum(p for c, p in marg.items() if c[0] == n_max)
            f_protonated.append(p_top)
        series[site.key] = f_protonated

    with Path(path).open("w", newline="") as f:
        w = csv.writer(f)
        header = ["pH"] + [f"{s.resname}_{s.resnum}" for s in sites]
        w.writerow(header)
        for k, ph in enumerate(pH_list):
            row = [f"{ph:.2f}"] + [f"{series[s.key][k]:.4f}" for s in sites]
            w.writerow(row)

    if fit_path is not None:
        with Path(fit_path).open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["resname", "resnum", "pKa_fit", "Hill_n",
                         "pKa_model", "interpretation"])
            for site in sites:
                pka, hill = henderson_hasselbalch_fit(pH_list, series[site.key])
                if pka is None:
                    interp = "(fit failed; check titration_curves.csv)"
                elif hill is None:
                    interp = ""
                elif hill > 1.2:
                    interp = "cooperative (n>1.2)"
                elif hill < 0.8:
                    interp = "anti-cooperative (n<0.8)"
                else:
                    interp = "near-ideal H-H (n≈1)"
                w.writerow([site.resname, site.resnum,
                             "" if pka is None else f"{pka:.3f}",
                             "" if hill is None else f"{hill:.3f}",
                             site.residue.pka_model,
                             interp])


# ---------------------------------------------------------------------------
# Coupling matrix outputs
# ---------------------------------------------------------------------------

def write_coupling_matrix_csv(W: CouplingMatrix, path: Path) -> None:
    """All non-zero W_ij(s,t) as a long-form CSV.

    Columns: resname_i, resnum_i, chem_i_label, resname_j, resnum_j,
             chem_j_label, W_kcal, W_kT
    """
    kT = RT_LN10 / 2.302585
    with Path(path).open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["resname_i","resnum_i","chem_i",
                     "resname_j","resnum_j","chem_j",
                     "W_kcal_per_mol","W_kT"])
        for i, site_i in enumerate(W.sites):
            for site_j in W.sites[i+1:]:
                for chem_i in site_i.chemistries:
                    if chem_i == site_i.ref_chem: continue
                    for chem_j in site_j.chemistries:
                        if chem_j == site_j.ref_chem: continue
                        v = W.W.get((site_i.key, chem_i), {}).get(
                            (site_j.key, chem_j), 0.0)
                        if v == 0.0: continue
                        w.writerow([
                            site_i.resname, site_i.resnum,
                            chemistry_label(site_i, chem_i),
                            site_j.resname, site_j.resnum,
                            chemistry_label(site_j, chem_j),
                            f"{v:.4f}", f"{v/kT:.4f}",
                        ])


def write_significant_pairs_csv(W: CouplingMatrix, path: Path,
                                  threshold_kT: float = 1.0) -> None:
    """Just the pairs with |W| > threshold, for network/visualization use.

    One row per site pair (the maximum |W| over chemistry combinations
    is reported). Suitable input for a coupling-network diagram.
    """
    kT = RT_LN10 / 2.302585
    rows = []
    for i, site_i in enumerate(W.sites):
        for site_j in W.sites[i+1:]:
            max_W = 0.0
            for chem_i in site_i.chemistries:
                if chem_i == site_i.ref_chem: continue
                for chem_j in site_j.chemistries:
                    if chem_j == site_j.ref_chem: continue
                    v = W.W.get((site_i.key, chem_i), {}).get(
                        (site_j.key, chem_j), 0.0)
                    if abs(v) > abs(max_W):
                        max_W = v
            if abs(max_W) >= threshold_kT * kT:
                rows.append((site_i.resname, site_i.resnum,
                              site_j.resname, site_j.resnum,
                              max_W, max_W / kT))
    rows.sort(key=lambda r: -abs(r[4]))  # strongest first
    with Path(path).open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["resname_i","resnum_i","resname_j","resnum_j",
                     "max_W_kcal","max_W_kT"])
        for r in rows:
            w.writerow([r[0],r[1],r[2],r[3], f"{r[4]:+.4f}", f"{r[5]:+.3f}"])


# ---------------------------------------------------------------------------
# Susceptibility from MC pair counts
# ---------------------------------------------------------------------------

def compute_susceptibility(result: MCResult) -> Dict[Tuple[SiteKey, SiteKey], float]:
    """χ_ij = ⟨P_i_protonated · P_j_protonated⟩ − ⟨P_i⟩⟨P_j⟩

    Off-diagonal: cooperative coupling (positive = sites titrate together).
    Diagonal (i=j): protonation capacitance, peaks at the site's pKa.

    Requires MC to have been run with `record_pair_counts=True`. If the
    pair-count store is empty, returns an empty dict.
    """
    if not getattr(result, "pair_counts", None):
        return {}
    n_post_equil = result.n_steps_total - result.n_steps_equilibration
    out: Dict[Tuple[SiteKey, SiteKey], float] = {}
    # Marginals: P(site i is in dominant-chemistry's most-protonated state)
    p_protonated: Dict[SiteKey, float] = {}
    for site in result.sites:
        n_max = max(c[0] for c in site.chemistries)
        p_protonated[site.key] = sum(
            p for c, p in result.marginals[site.key].items() if c[0] == n_max)
    # Joint and cross-correlation
    for i, site_i in enumerate(result.sites):
        for j, site_j in enumerate(result.sites):
            joint_count = result.pair_counts.get((site_i.key, site_j.key), 0)
            joint = joint_count / max(1, n_post_equil)
            chi = joint - p_protonated[site_i.key] * p_protonated[site_j.key]
            out[(site_i.key, site_j.key)] = chi
    return out


def write_susceptibility_csv(chi: Dict[Tuple[SiteKey, SiteKey], float],
                               sites: List[CouplingSite], path: Path) -> None:
    """Long-form: resname_i, resnum_i, resname_j, resnum_j, chi_ij."""
    with Path(path).open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["resname_i","resnum_i","resname_j","resnum_j",
                     "chi","interpretation"])
        for (ki, kj), v in chi.items():
            si = next(s for s in sites if s.key == ki)
            sj = next(s for s in sites if s.key == kj)
            if ki == kj:
                interp = "diagonal (capacitance)"
            elif v > 0.05:
                interp = "cooperative (titrate together)"
            elif v < -0.05:
                interp = "anti-cooperative (titrate oppositely)"
            else:
                interp = ""
            w.writerow([si.resname, si.resnum, sj.resname, sj.resnum,
                         f"{v:+.6f}", interp])


# ---------------------------------------------------------------------------
# State-map comparison
# ---------------------------------------------------------------------------

def compare_state_maps(csv_a: Path, csv_b: Path,
                        out_path: Path) -> Dict[str, int]:
    """Diff two state-map CSVs. Returns a count summary."""
    def read_map(p: Path):
        out = {}
        with Path(p).open() as f:
            for row in csv.DictReader(f):
                key = (row["resname"], int(row["resnum"]))
                out[key] = (row.get("chemistry") or row.get("state", ""),
                            row.get("state_name", ""))
        return out
    A = read_map(csv_a)
    B = read_map(csv_b)
    keys = sorted(set(A) | set(B))
    n_match = n_diff = n_only_a = n_only_b = 0
    rows = []
    for k in keys:
        a = A.get(k); b = B.get(k)
        if a is None:
            n_only_b += 1; status = "only_in_B"
        elif b is None:
            n_only_a += 1; status = "only_in_A"
        elif a[0] == b[0]:
            n_match += 1; status = "match"
        else:
            n_diff += 1; status = "DIFFER"
        rows.append((k[0], k[1],
                      a[0] if a else "",
                      b[0] if b else "",
                      status))
    with Path(out_path).open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["resname","resnum","chemistry_A","chemistry_B","status"])
        for r in rows:
            w.writerow(r)
    return {"match": n_match, "differ": n_diff,
            "only_in_A": n_only_a, "only_in_B": n_only_b}
