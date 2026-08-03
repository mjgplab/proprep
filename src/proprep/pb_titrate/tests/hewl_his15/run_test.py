#!/usr/bin/env python3
"""HEWL HIS-15 — first real-protein 3-state HIS test.

HEWL HIS-15 NMR pKa is around 5.5-5.8 (down-shifted from cpinutil
pKa_model 6.6 due to its location near the active-site cleft).

Validates:
  - Multi-state machinery on a real HIS-containing protein
  - Tautomer preference (HID vs HIE) in a real environment
  - Effective HIS pKa calculation via the partition function

Output: HIP/HID/HIE populations across a pH range, effective pKa.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent))

from proprep.pb_titrate.intrinsic import (
    compute_state_energies, populations_by_chemistry_group,
    populations_by_prot_count, state_populations_at_pH,
    state_free_energies_at_pH,
)
from proprep.pb_titrate.residues import get_residue


def main():
    res = get_residue("HIP")
    print(f"=== HEWL HIS-15 multi-state PB ===")
    print(f"  pKa_model (cpinutil) = {res.pka_model}")
    print(f"  Reference: NMR pKa(HIS-15) ≈ 5.5-5.8")
    print(f"  εin=4, εout=80, I=150 mM, bcopt=5\n")

    se = compute_state_energies(
        cluster_prmtop=HERE / "hewl.prmtop",
        cluster_rst7=  HERE / "hewl.rst7",
        target_resname="HIP",
        target_resnum_1based=15,
        residue=res,
        work_dir=HERE / "pbsa_work",
        epsin=4.0, space=0.5, nfocus=2, bcopt=5, istrng=150,
    )

    state_labels = {0: "HIP", 1: "HID", 2: "HIE"}
    print(f"  Per-state ΔΔG (vs ACE-HIP-NME model, kcal/mol):")
    for s_idx, info in se.items():
        print(f"    {state_labels[s_idx]} (n={info['prot_count']}, "
              f"pka={info['pka_corr']}): "
              f"E_p={info['e_protein']:+9.3f}  E_m={info['e_model']:+8.3f}  "
              f"ΔΔG={info['ddG_PB_minus_model']:+7.3f}")

    print(f"\n  Populations across pH:")
    print(f"    pH    HIP      HID      HIE      total_DEPROT (HID+HIE)")
    for pH in [4.0, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0]:
        pops = state_populations_at_pH(se, pH)
        labeled = {state_labels[s]: p for s, p in pops.items()}
        deprot_total = labeled["HID"] + labeled["HIE"]
        print(f"    {pH:>4.1f}  {labeled['HIP']:.3f}    "
              f"{labeled['HID']:.3f}    {labeled['HIE']:.3f}    "
              f"{deprot_total:.3f}")

    # Find effective pKa by bisection — pH where P(HIP) = 0.5
    print(f"\n  Bisection for effective pKa (P(HIP)=0.5):")
    lo, hi = 0.0, 14.0
    for _ in range(40):
        mid = (lo + hi) / 2
        pops = state_populations_at_pH(se, mid)
        p_hip = pops[0]
        if p_hip > 0.5:
            lo = mid
        else:
            hi = mid
    eff_pKa = (lo + hi) / 2
    print(f"    Effective HIS pKa ≈ {eff_pKa:.3f}")
    print(f"    NMR target: 5.5-5.8")
    print(f"    Δ(predicted − NMR_mid) = {eff_pKa - 5.65:+.3f}")

    # Tautomer preference at high pH
    pops_high = state_populations_at_pH(se, 10.0)
    hid_hie_ratio = pops_high[1] / max(pops_high[2], 1e-30)
    print(f"\n  At pH=10 (deprotonated regime):")
    print(f"    HID/HIE ratio = {hid_hie_ratio:.2f}")
    print(f"    (cpinutil pure-solution: 4:1 HID-favored)")


if __name__ == "__main__":
    main()
