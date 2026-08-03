#!/usr/bin/env python3
"""BPTI ASP-3 test (v2) — using titrate.intrinsic.

Should reproduce the v1 result exactly: pKa(AS4-3) = 3.27.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent))  # titrate/

from proprep.pb_titrate.intrinsic import compute_intrinsic_pka


def main():
    r = compute_intrinsic_pka(
        cluster_prmtop=HERE / "bpti.prmtop",
        cluster_rst7=  HERE / "bpti.rst7",
        target_resname="AS4",
        target_resnum_1based=3,
        work_dir=HERE / "pbsa_work_v2",
        # PB params (forwarded to run_pbsa):
        epsin=4.0, space=0.5, nfocus=2, bcopt=5, istrng=150,
    )
    print(f"=== BPTI AS4-3 via titrate.intrinsic ===")
    print(f"  ΔG(protein) = {r['ddG_protein']:+8.3f}  "
          f"ΔG(model) = {r['ddG_model']:+8.3f}")
    print(f"  ΔΔG = {r['ddG_kcal']:+8.3f}  ΔpKa = {r['dpKa']:+6.3f}  "
          f"pKa = {r['pKa']:6.3f}")
    print(f"  v1 baseline: pKa = 3.271 (NMR ~3.0–3.4)")
    print(f"  Δ(v2 − v1)  = {r['pKa'] - 3.271:+.4f}")


if __name__ == "__main__":
    main()
