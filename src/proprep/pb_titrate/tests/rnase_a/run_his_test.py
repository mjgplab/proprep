#!/usr/bin/env python3
"""Full pipeline on RNase A's 4 HIS residues.

Expected NMR pKas (from various NMR studies, ~150 mM, 25°C):
  HIS-12  ≈ 5.8   (active site)
  HIS-48  ≈ 6.4   (buried, conserved)
  HIS-105 ≈ 6.7
  HIS-119 ≈ 6.2   (active site, partner to HIS-12)

The active-site pair HIS-12/HIS-119 is ~5 Å apart and shows cooperative
behavior. This is the textbook test for whether a pKa method captures
HIS coupling.
"""
import sys
import pickle
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent))

from proprep.pb_titrate.coupling import (
    make_coupling_site, compute_self_energies, compute_pair_couplings,
    coupling_summary,
)
from proprep.pb_titrate.solvers.monte_carlo import (
    run_mcmc, assign_states_from_mc, chemistry_label,
)
from proprep.pb_titrate.solvers.enumerate import enumerate_all
from proprep.pb_titrate.select_solver import recommend
from proprep.pb_titrate.intrinsic import (
    populations_by_chemistry_group, state_populations_at_pH,
    state_free_energies_at_pH,
)
from proprep.pb_titrate.analysis import (
    write_marginals_csv, write_effective_pkas_csv,
    titration_scan, write_titration_csv,
    write_coupling_matrix_csv, write_significant_pairs_csv,
    compute_susceptibility, write_susceptibility_csv,
)

PRM = HERE / "rnase_a.prmtop"
RST = HERE / "rnase_a.rst7"
OUT = HERE / "his_test_out"
OUT.mkdir(exist_ok=True)
PBSA = {"epsin": 4.0, "space": 0.5, "nfocus": 2, "bcopt": 5, "istrng": 150}


# Expected NMR values (mid-range, from multiple sources)
NMR = {12: 5.8, 48: 6.4, 105: 6.7, 119: 6.2}


def main():
    t0 = time.time()
    print("=== RNase A: 4 HIS multi-state pipeline ===")
    print(f"  Reference NMR pKas: {NMR}")

    sites = [make_coupling_site("HIP", n) for n in [12, 48, 105, 119]]
    for s in sites:
        labels = [chemistry_label(s, c) for c in s.chemistries]
        print(f"  {s.resname}-{s.resnum}: chemistries={labels}, "
              f"ref={chemistry_label(s, s.ref_chem)}")

    # --- Step 1: self-energies ---
    print(f"\n--- Step 1: self-energies ({len(sites)} sites × 3 states × 2 PB calls) ---")
    se = compute_self_energies(
        sites, PRM, RST, workers=4,
        work_root=OUT / "coupling_work", pbsa_params=PBSA,
    )
    for site in sites:
        print(f"  {site.resname}-{site.resnum:<3}:")
        for chem in site.chemistries:
            ddg = se.intrinsic_ddG[(site.key, chem)]
            print(f"    {chemistry_label(site, chem):<5} ΔΔG_PB={ddg:+.3f} kcal/mol")

    # --- Step 2: pair couplings ---
    print(f"\n--- Step 2: pair couplings ---")
    W = compute_pair_couplings(
        sites, PRM, RST, workers=1,
        work_root=OUT / "coupling_work", pbsa_params=PBSA,
    )
    print(coupling_summary(W))

    # Save precompute for re-use
    pkl = OUT / "coupling.pkl"
    pkl.write_bytes(pickle.dumps({
        "sites": sites, "self_energies": se, "coupling": W}))

    # --- Step 3: Solver auto-pick + run ---
    print(f"\n--- Step 3: auto-pick solver ---")
    rec = recommend(se, W, pH=7.0)
    print(f"  Recommended: {rec.solver}")
    print(f"  Reason: {rec.reason}")

    # 4 HIS × 3 chemistries each = 81 total combos → enumerate
    print(f"\n--- Step 4: exact enumeration (81 states, gold standard) ---")
    enum = enumerate_all(se, W, pH=7.0)
    print(f"  Marginals at pH 7:")
    for site in sites:
        marg = enum.marginals[site.key]
        labeled = ", ".join(f"{chemistry_label(site, c)}={p:.3f}"
                              for c, p in sorted(marg.items()))
        print(f"    {site.resname}-{site.resnum:<3}  {labeled}")

    # Top 3 joint states
    print(f"\n  Top 5 joint states (out of {enum.n_states_enumerated}):")
    for state, p in sorted(enum.state_probabilities.items(),
                             key=lambda x: -x[1])[:5]:
        labels = ' '.join(f"{site.resname[-3:]}{site.resnum}={chemistry_label(site, state[i])}"
                           for i, site in enumerate(sites))
        print(f"    P={p:.4f}  [{labels}]")

    # --- Step 5: MC for analyses ---
    print(f"\n--- Step 5: MC at pH 7 with pair counts ---")
    mc = run_mcmc(se, W, pH=7.0, n_steps=300_000, n_equil=30_000,
                    record_pair_counts=True, seed=42)
    print(f"  Acceptance: {mc.acceptance_ratio*100:.1f}%")
    write_marginals_csv(mc, OUT / "mc_marginals_pH7.csv")
    write_effective_pkas_csv(mc, OUT / "mc_effective_pkas_pH7.csv")
    chi = compute_susceptibility(mc)
    write_susceptibility_csv(chi, sites, OUT / "susceptibility_pH7.csv")

    # --- Step 6: Titration scan (the crown jewel) ---
    print(f"\n--- Step 6: titration scan over pH 3-10 ---")
    pH_grid = [3, 4, 5, 5.5, 6, 6.25, 6.5, 6.75, 7, 7.5, 8, 9, 10]
    scan = titration_scan(se, W, pH_grid,
                            n_steps=80_000, n_equil=8_000)
    write_titration_csv(scan, OUT / "titration_curves.csv",
                          fit_path=OUT / "titration_fit.csv")

    # --- Comparison vs NMR ---
    print(f"\n=== Comparison: predicted vs NMR ===")
    print(f"  {'site':<10}  {'NMR':<6}  {'pred':<6}  {'Hill':<6}  {'Δ':<6}")
    import csv
    with (OUT / "titration_fit.csv").open() as f:
        for row in csv.DictReader(f):
            n = int(row["resnum"])
            pka = float(row["pKa_fit"]) if row["pKa_fit"] else None
            hill = float(row["Hill_n"]) if row["Hill_n"] else None
            nmr = NMR[n]
            d = pka - nmr if pka is not None else None
            print(f"  HIS-{n:<5}  {nmr:<6.2f}  "
                  f"{'-' if pka is None else f'{pka:<6.2f}'}  "
                  f"{'-' if hill is None else f'{hill:<6.2f}'}  "
                  f"{'-' if d is None else f'{d:+5.2f}'}")

    # Coupling matrix CSVs
    write_coupling_matrix_csv(W, OUT / "coupling_matrix.csv")
    write_significant_pairs_csv(W, OUT / "significant_pairs.csv",
                                   threshold_kT=0.5)
    print(f"\n=== All outputs in {OUT} ===")
    for p in sorted(OUT.iterdir()):
        if p.is_file():
            print(f"  {p.name}")
    print(f"\n  Total wall time: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
