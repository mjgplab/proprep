#!/usr/bin/env python3
"""End-to-end test of all titrate.analysis outputs on BPTI 3-site.

Reuses the precomputed coupling.pkl from mc_test.py.
"""
import sys
import pickle
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent))

from proprep.pb_titrate.analysis import (
    write_marginals_csv, write_effective_pkas_csv,
    titration_scan, write_titration_csv,
    write_coupling_matrix_csv, write_significant_pairs_csv,
    compute_susceptibility, write_susceptibility_csv,
    compare_state_maps,
)
from proprep.pb_titrate.solvers.monte_carlo import run_mcmc
from proprep.pb_titrate.coupling import make_coupling_site, compute_self_energies, compute_pair_couplings

PRM = HERE.parent / "bpti_asp3" / "bpti.prmtop"
RST = HERE.parent / "bpti_asp3" / "bpti.rst7"


def main():
    out = HERE / "analysis_out"
    out.mkdir(exist_ok=True)

    # Load or compute precompute
    pkl = HERE / "coupling.pkl"
    if pkl.exists():
        print(f"Loading cached precompute from {pkl}")
        data = pickle.loads(pkl.read_bytes())
    else:
        sites = [make_coupling_site("AS4",3), make_coupling_site("TYR",10),
                  make_coupling_site("LYS",15)]
        pbsa_params = {"epsin":4.0,"space":0.5,"nfocus":2,"bcopt":5,"istrng":150}
        se = compute_self_energies(sites, PRM, RST, workers=1,
            work_root=HERE / "coupling_work", pbsa_params=pbsa_params)
        W  = compute_pair_couplings(sites, PRM, RST, workers=1,
            work_root=HERE / "coupling_work", pbsa_params=pbsa_params)
        data = {"sites": sites, "self_energies": se, "coupling": W}
        pkl.write_bytes(pickle.dumps(data))
    se = data["self_energies"]; W = data["coupling"]; sites = data["sites"]

    # --- Test 1: coupling matrix CSVs (no MC needed) ---
    print("\n--- Test 1: coupling matrix CSVs ---")
    write_coupling_matrix_csv(W, out / "coupling_matrix.csv")
    write_significant_pairs_csv(W, out / "significant_pairs.csv", threshold_kT=0.5)
    print(f"  wrote {out / 'coupling_matrix.csv'}")
    print(f"  wrote {out / 'significant_pairs.csv'}")

    # --- Test 2: single-pH MC marginals + effective pKa ---
    print("\n--- Test 2: MC at pH 7 ---")
    mc7 = run_mcmc(se, W, pH=7.0, n_steps=200_000, n_equil=20_000,
                    record_pair_counts=True, seed=42)
    write_marginals_csv(mc7, out / "marginals_pH7.csv")
    write_effective_pkas_csv(mc7, out / "effective_pkas_pH7.csv")
    print(f"  wrote {out / 'marginals_pH7.csv'}")
    print(f"  wrote {out / 'effective_pkas_pH7.csv'}")

    # --- Test 3: susceptibility ---
    print("\n--- Test 3: susceptibility χ_ij from MC pair counts ---")
    chi = compute_susceptibility(mc7)
    write_susceptibility_csv(chi, sites, out / "susceptibility_pH7.csv")
    print(f"  wrote {out / 'susceptibility_pH7.csv'}")
    # Print salient values
    for site in sites:
        chi_diag = chi[(site.key, site.key)]
        print(f"    χ[{site.resname}-{site.resnum}, self] = {chi_diag:+.4f} "
              f"(capacitance peaks near pKa)")

    # --- Test 4: titration curves ---
    print("\n--- Test 4: titration scan over pH 3-12 ---")
    pH_grid = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    scan = titration_scan(se, W, pH_grid,
                            n_steps=50_000, n_equil=5_000)
    write_titration_csv(scan, out / "titration_curves.csv",
                          fit_path=out / "titration_fit.csv")
    print(f"  wrote {out / 'titration_curves.csv'}")
    print(f"  wrote {out / 'titration_fit.csv'}")

    # Print fit summary
    print(f"\n  H-H fit results:")
    import csv
    with (out / "titration_fit.csv").open() as f:
        for row in csv.DictReader(f):
            print(f"    {row['resname']}-{row['resnum']}  "
                  f"pKa_fit={row['pKa_fit']}  Hill={row['Hill_n']}  "
                  f"({row['interpretation']})")

    # --- Test 5: state-map compare ---
    print("\n--- Test 5: compare state maps ---")
    # Make two: one from MC at pH 5, one at pH 9 (very different states)
    mc5 = run_mcmc(se, W, pH=5.0, n_steps=50_000, n_equil=5_000, seed=43)
    mc9 = run_mcmc(se, W, pH=9.0, n_steps=50_000, n_equil=5_000, seed=43)
    from proprep.pb_titrate.solvers.monte_carlo import assign_states_from_mc
    def write_assign(mc, p):
        with p.open("w") as f:
            f.write("resname,resnum,chemistry,prot_count,pka_corr,state_name\n")
            assign = assign_states_from_mc(mc)
            for site in sites:
                chem = assign[site.key]
                state = site.chem_to_state[chem]
                from proprep.pb_titrate.solvers.monte_carlo import chemistry_label
                f.write(f"{site.resname},{site.resnum},"
                        f"{chemistry_label(site, chem)},{chem[0]},{chem[1]},{state.name}\n")
    write_assign(mc5, out / "state_map_pH5.csv")
    write_assign(mc9, out / "state_map_pH9.csv")
    summary = compare_state_maps(out / "state_map_pH5.csv",
                                   out / "state_map_pH9.csv",
                                   out / "state_map_diff.csv")
    print(f"  Compare pH 5 vs pH 9: {summary}")
    print(f"  wrote {out / 'state_map_diff.csv'}")

    print(f"\n=== All outputs in {out} ===")
    for p in sorted(out.iterdir()):
        print(f"  {p.name}")


if __name__ == "__main__":
    main()
