#!/usr/bin/env python3
"""End-to-end pipeline test on BPTI 3-site, with cluster cutout enabled.

Same 3 sites as bpti_multi/run_test.py (AS4-3, TYR-10, LYS-15) but with
`cluster_radius=15` to exercise the cluster-cutout code path through the
full pipeline:

  load -> per-site precompute (compute_state_energies with cutout)
       -> mean-field iteration (background propagation across sites)
       -> chemistry-group populations -> state assignment -> CSV.

Pass criteria at pH 7:
  AS4-3   dominant chemistry = DEPROT  (n=0)        NMR pKa ≈ 3.0–3.4
  TYR-10  dominant chemistry = PROT    (n=1)        NMR pKa ≈ 10.5
  LYS-15  dominant chemistry = PROT    (n=1)        pKa_model 10.4

The test prints PASS/FAIL per site and exits nonzero if any site is wrong.
This is the gate to clear before handing the package to another agent.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent))

from proprep.pb_titrate.solvers.mean_field import make_site, solve, write_state_map
from proprep.pb_titrate.apply_state import build_production_prmtop, state_map_from_sites

BPTI_PRMTOP = HERE.parent / "bpti_asp3" / "bpti.prmtop"
BPTI_RST7   = HERE.parent / "bpti_asp3" / "bpti.rst7"

# cpinutil prot_count convention is "absolute # of titrating H on the residue":
#   AS4 DEPROT = 0,   AS4 PROT = 1
#   TYR PROT   = 1,   TYR DEPROT = 0
#   LYS PROT   = 3 (NZ-H3+),   LYN (DEPROT) = 2 (NZ-H2)
EXPECTED_PROT_COUNT = {
    ("AS4",  3): 0,   # DEPROT
    ("TYR", 10): 1,   # PROT
    ("LYS", 15): 3,   # PROT (NZ-H3+, charged)
}


def main():
    sites = [
        make_site("AS4",  3, initial="auto", pH=7.0),
        make_site("TYR", 10, initial="auto", pH=7.0),
        make_site("LYS", 15, initial="auto", pH=7.0),
    ]
    print("=== BPTI end-to-end (with cluster_radius=15) ===")
    for s in sites:
        print(f"  init  {s.resname}-{s.resnum:<3}  pKa_model={s.residue.pka_model:>5.2f} "
              f" -> {s.chemistry_label()}")

    result = solve(
        sites=sites,
        cluster_prmtop=BPTI_PRMTOP,
        cluster_rst7=  BPTI_RST7,
        pH=7.0,
        max_iter=3,
        workers=3,
        work_root=HERE / "mf_work",
        # cluster_radius is consumed by compute_state_energies (named kwarg);
        # remaining keys forward to run_pbsa.
        pbsa_params={"cluster_radius": 15.0, "bond_depth": 3,
                     "epsin": 4.0, "space": 0.5, "nfocus": 2,
                     "bcopt": 5, "istrng": 150},
    )
    write_state_map(result["sites"], HERE / "state_final.csv")

    print("\n=== Final state map ===")
    last = {(r["resname"], r["resnum"]): r for r in result["history"][-1]}
    n_fail = 0
    for s in result["sites"]:
        r = last.get(s.key)
        if r is None or r.get("error"):
            print(f"  FAIL  {s.resname}-{s.resnum}  error: {r.get('error') if r else 'no result'}")
            n_fail += 1
            continue
        # Aggregate populations by prot_count (chemistry key is (n, pka_corr))
        pops_by_n: dict = {}
        for (n, _pka), p in r["pop_by_chem"].items():
            pops_by_n[n] = pops_by_n.get(n, 0.0) + p
        dominant_n = max(pops_by_n, key=pops_by_n.get)
        expected_n = EXPECTED_PROT_COUNT[s.key]
        ok = dominant_n == expected_n
        tag = "PASS" if ok else "FAIL"
        if not ok: n_fail += 1
        pop_str = "  ".join(f"n={n}: {p:.3f}" for n, p in sorted(pops_by_n.items()))
        print(f"  {tag}  {s.resname}-{s.resnum:<3}  dominant n={dominant_n} "
              f"(expected {expected_n})  {pop_str}")

    print(f"\nWrote {HERE / 'state_final.csv'}")
    if n_fail:
        print(f"\n*** {n_fail} site(s) FAILED at solver step ***")
        sys.exit(1)

    # ---- Apply state map and write production prmtop ----
    # Note: BPTI is a dry protein with no ions, so we cannot exercise
    # rebalance_ions here. The rebalance path is tested separately in
    # tests/rebalance_unit/.
    print("\n=== build_production_prmtop (no rebalance — BPTI has no ions) ===")
    sm = state_map_from_sites(result["sites"])
    out_info = build_production_prmtop(
        input_prmtop=BPTI_PRMTOP,
        input_rst7=  BPTI_RST7,
        state_map=   sm,
        output_prefix=HERE / "bpti_production",
        rebalance=False,
    )
    print(f"  net q  in: {out_info['net_q_in']:+.3f}  "
          f"after_apply: {out_info['net_q_after_apply']:+.3f}  "
          f"out: {out_info['net_q_out']:+.3f}")
    print(f"  atoms  in: {out_info['n_atoms_in']}   "
          f"out: {out_info['n_atoms_out']}")
    print(f"  per-site charge check:")
    n_charge_fail = 0
    for site_id, info in out_info["site_check"].items():
        ok = abs(info["delta"]) < 1e-3
        tag = "ok" if ok else "FAIL"
        if not ok: n_charge_fail += 1
        print(f"    {tag}  {site_id[0]}-{site_id[1]:<3}  "
              f"actual={info['actual']:+.4f}  expected={info['expected']:+.4f}")
    print(f"  output: {out_info['output_prmtop'].name} + .rst7")

    if n_charge_fail:
        print(f"\n*** {n_charge_fail} per-site charge mismatch ***")
        sys.exit(1)
    print("\n*** ALL CHECKS PASSED ***")


if __name__ == "__main__":
    main()
