#!/usr/bin/env python3
"""Active-site cluster test: 4 HIS + LYS-41 + ASP-121 + GLU-111.

Diagnostic for whether the HIS-12 over-acidification (predicted 3.35 vs
NMR 5.8) is fixed by letting LYS-41 also titrate (it sits ~5 Å from
HIS-12 and was previously locked at +1).
"""
import sys, time, pickle
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent))

from proprep.pb_titrate.coupling import (
    make_coupling_site, compute_self_energies, compute_pair_couplings,
    coupling_summary)
from proprep.pb_titrate.solvers.monte_carlo import run_mcmc, chemistry_label
from proprep.pb_titrate.analysis import (
    titration_scan, write_titration_csv, write_significant_pairs_csv,
    write_coupling_matrix_csv)

PRM = HERE / "rnase_a.prmtop"
RST = HERE / "rnase_a.rst7"
OUT = HERE / "active_site_out"
OUT.mkdir(exist_ok=True)
PBSA = {"epsin": 4.0, "space": 0.5, "nfocus": 2, "bcopt": 5, "istrng": 150}

NMR = {("HIP",12): 5.8, ("HIP",48): 6.4, ("HIP",105): 6.7, ("HIP",119): 6.2,
       ("LYS",41): 9.0,  # downshifted active-site Lys, exact value debated
       ("AS4",121): 4.5, # active-site Asp, slightly upshifted from model
       ("GL4",111): 3.5}  # surface, model 4.4


def main():
    t0 = time.time()
    site_specs = [("HIP",12),("HIP",48),("HIP",105),("HIP",119),
                   ("LYS",41),("AS4",121),("GL4",111)]
    sites = [make_coupling_site(rn, n) for rn, n in site_specs]

    print("=== RNase A active-site cluster (7 sites including LYS-41, ASP-121) ===")
    for s in sites:
        labels = [chemistry_label(s, c) for c in s.chemistries]
        print(f"  {s.resname}-{s.resnum}: chemistries={labels}")

    se = compute_self_energies(sites, PRM, RST, workers=4,
        work_root=OUT / "coupling_work", pbsa_params=PBSA)
    W = compute_pair_couplings(sites, PRM, RST, workers=1,
        work_root=OUT / "coupling_work", pbsa_params=PBSA)
    pickle.dumps  # noqa
    (OUT / "coupling.pkl").write_bytes(pickle.dumps({
        "sites": sites, "self_energies": se, "coupling": W}))
    print(coupling_summary(W))

    pH_grid = [3, 4, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]
    print(f"\n--- titration scan ({len(pH_grid)} pH values) ---")
    scan = titration_scan(se, W, pH_grid, n_steps=80_000, n_equil=8_000)
    write_titration_csv(scan, OUT / "titration_curves.csv",
                          fit_path=OUT / "titration_fit.csv")
    write_coupling_matrix_csv(W, OUT / "coupling_matrix.csv")
    write_significant_pairs_csv(W, OUT / "significant_pairs.csv",
                                   threshold_kT=1.0)

    print(f"\n=== Comparison vs NMR ===")
    print(f"  {'site':<14}  {'NMR':<5}  {'pred':<6}  {'Hill':<6}  {'Δ':<6}")
    import csv
    with (OUT / "titration_fit.csv").open() as f:
        for row in csv.DictReader(f):
            key = (row["resname"], int(row["resnum"]))
            nmr = NMR.get(key, "?")
            pka = float(row["pKa_fit"]) if row["pKa_fit"] else None
            hill = float(row["Hill_n"]) if row["Hill_n"] else None
            d = (pka - nmr) if (pka is not None and nmr != "?") else None
            label = f"{row['resname']}-{row['resnum']}"
            print(f"  {label:<14}  "
                  f"{('?' if nmr=='?' else f'{nmr:<5.2f}'):<5}  "
                  f"{'-' if pka is None else f'{pka:<6.2f}'}  "
                  f"{'-' if hill is None else f'{hill:<6.2f}'}  "
                  f"{'-' if d is None else f'{d:+5.2f}'}")
    print(f"\n  Wall: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
