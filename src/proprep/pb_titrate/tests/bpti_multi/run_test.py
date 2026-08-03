#!/usr/bin/env python3
"""Multi-site BPTI mean-field smoke test.

Three sites of three different residue types: AS4-3, TYR-10, LYS-15.
At pH 7:
  - AS4-3  expected DEPROT (NMR pKa ≈ 3.0–3.4)
  - TYR-10 expected PROT   (NMR pKa ≈ 10.5)
  - LYS-15 expected PROT   (LYS pKa_model = 10.4, exposed)

This validates:
  (a) mean_field.solve runs over a heterogeneous site list
  (b) background_state_map propagation works (each site's PB calc
      sees the others at their current state)
  (c) state assignments converge to chemically reasonable answers
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent))  # titrate/

from proprep.pb_titrate.solvers.mean_field import make_site, solve, write_state_map

BPTI_PRMTOP = HERE.parent / "bpti_asp3" / "bpti.prmtop"
BPTI_RST7   = HERE.parent / "bpti_asp3" / "bpti.rst7"


def main():
    sites = [
        make_site("AS4",  3, initial="auto", pH=7.0),
        make_site("TYR", 10, initial="auto", pH=7.0),
        make_site("LYS", 15, initial="auto", pH=7.0),
    ]
    print("Initial state assignments (from pKa_model vs pH=7):")
    for s in sites:
        print(f"  {s.resname}-{s.resnum:<3} pKa_model={s.residue.pka_model:>5.2f}  "
              f"-> {s.chemistry_label()}")

    result = solve(
        sites=sites,
        cluster_prmtop=BPTI_PRMTOP,
        cluster_rst7=  BPTI_RST7,
        pH=7.0,
        max_iter=3,
        workers=3,  # 3 sites; use one worker each
        work_root=HERE / "mf_work",
        pbsa_params={"epsin": 4.0, "space": 0.5, "nfocus": 2,
                     "bcopt": 5, "istrng": 150},
    )
    write_state_map(result["sites"], HERE / "state_final.csv")

    print("\nFinal state map:")
    for s in result["sites"]:
        for r in result["history"][-1]:
            if (r["resname"], r["resnum"]) == s.key and not r.get("error"):
                pops_str = "  ".join(
                    f"{n}={p:.3f}"
                    for (n, _pka), p in sorted(r["pop_by_chem"].items()))
                print(f"  {s.resname}-{s.resnum:<3}  chemistry="
                      f"{s.chemistry_label():<6}  pops_by_chem: {pops_str}")
                break

    print(f"\nWrote {HERE / 'state_final.csv'}")


if __name__ == "__main__":
    main()
