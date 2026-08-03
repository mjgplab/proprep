#!/usr/bin/env python3
"""End-to-end multi-state PB tests.

(1) ACE-AS4-NME standalone (model = "protein" = same compound):
    Should give populations matching pure-solution expectations.
    All 4 PROT tautomers averaged; effective pKa should = 4.0.

(2) ACE-HIP-NME standalone:
    HIP/HID/HIE populations should match cpinutil pure-solution conventions.
    At pH 7: HIE≈16%, HID≈64%, HIP≈20%; effective HIS pKa ≈ 6.5.

(3) BPTI AS4-3 with multi-state averaging:
    Replaces the arbitrary "state 1" tautomer choice with proper Boltzmann
    average over all 4. Should be very close to the single-tautomer result
    (3.27) since vacuum tautomers are equivalent and the protein
    environment shouldn't strongly prefer one O over the other.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent))  # titrate/

from proprep.pb_titrate.intrinsic import (
    compute_state_energies, populations_by_chemistry_group,
    populations_by_prot_count, effective_pKa_two_state,
    state_populations_at_pH,
)
from proprep.pb_titrate.residues import get_residue
from proprep.pb_titrate.model_compounds import build_model

PBSA_PARAMS = {"epsin": 4.0, "space": 0.5, "nfocus": 2,
                "bcopt": 5, "istrng": 150}


def test_as4_standalone():
    print("\n=== Test 1: ACE-AS4-NME standalone (5 states, 4 PROT tautomers) ===")
    res = get_residue("AS4")
    prm, rst = build_model("AS4")

    se = compute_state_energies(
        cluster_prmtop=prm, cluster_rst7=rst,
        target_resname="AS4", target_resnum_1based=2,
        residue=res, model_prmtop=prm, model_rst7=rst,
        work_dir=HERE / "as4_standalone",
        cluster_fillratio=4.0, model_fillratio=4.0,
        **PBSA_PARAMS,
    )
    print(f"  ΔΔG_PB_minus_model per state (should be ~0 since model=protein):")
    for s_idx, info in se.items():
        print(f"    state {s_idx} (n={info['prot_count']}, "
              f"pka={info['pka_corr']}): ΔΔG={info['ddG_PB_minus_model']:+.4f} kcal/mol")

    print("\n  Populations and effective pKa:")
    for pH in [3.0, 4.0, 5.0]:
        chem = populations_by_chemistry_group(se, pH)
        by_n = populations_by_prot_count(se, pH)
        eff_pKa = effective_pKa_two_state(se, pH)
        print(f"    pH={pH:>3.1f}: P(DEPROT)={by_n[0]:.4f}  "
              f"P(PROT)={by_n[1]:.4f}  effective_pKa={eff_pKa:.4f}")
    print(f"    Expected at pH=4: 50/50, effective_pKa=4.00")


def test_hip_standalone():
    print("\n=== Test 2: ACE-HIP-NME standalone (3 states: HIP, HID, HIE) ===")
    res = get_residue("HIP")
    prm, rst = build_model("HIP")

    se = compute_state_energies(
        cluster_prmtop=prm, cluster_rst7=rst,
        target_resname="HIP", target_resnum_1based=2,
        residue=res, model_prmtop=prm, model_rst7=rst,
        work_dir=HERE / "hip_standalone",
        cluster_fillratio=4.0, model_fillratio=4.0,
        **PBSA_PARAMS,
    )
    print(f"  ΔΔG_PB_minus_model per state (should be ~0):")
    state_labels = {0: "HIP", 1: "HID", 2: "HIE"}
    for s_idx, info in se.items():
        print(f"    {state_labels[s_idx]} (n={info['prot_count']}, "
              f"pka={info['pka_corr']}): ΔΔG={info['ddG_PB_minus_model']:+.4f} kcal/mol")

    print("\n  Populations:")
    for pH in [4.0, 6.0, 6.5, 7.0, 8.0, 10.0]:
        pops = state_populations_at_pH(se, pH)
        chem = populations_by_chemistry_group(se, pH)
        labeled = {state_labels[s_idx]: p for s_idx, p in pops.items()}
        print(f"    pH={pH:>3.1f}: HIP={labeled['HIP']:.3f}  "
              f"HID={labeled['HID']:.3f}  HIE={labeled['HIE']:.3f}")
    print(f"    Expected at pH=7: HID/HIE ≈ 4:1, HIP ≈ 20%")


def test_bpti_as4_multistate():
    print("\n=== Test 3: BPTI AS4-3 with proper multi-state averaging ===")
    res = get_residue("AS4")
    bpti_prm = HERE.parent / "bpti_asp3" / "bpti.prmtop"
    bpti_rst = HERE.parent / "bpti_asp3" / "bpti.rst7"

    se = compute_state_energies(
        cluster_prmtop=bpti_prm, cluster_rst7=bpti_rst,
        target_resname="AS4", target_resnum_1based=3,
        residue=res, work_dir=HERE / "bpti_as4_multistate",
        **PBSA_PARAMS,
    )
    print(f"  Per-state ΔΔG (vs ACE-AS4-NME model):")
    for s_idx, info in se.items():
        print(f"    state {s_idx} (n={info['prot_count']}, "
              f"pka={info['pka_corr']}): "
              f"E_p={info['e_protein']:+9.3f}  E_m={info['e_model']:+8.3f}  "
              f"ΔΔG={info['ddG_PB_minus_model']:+7.3f}")

    by_n = populations_by_prot_count(se, pH=7.0)
    eff_pKa = effective_pKa_two_state(se, pH_anchor=7.0)
    print(f"\n  At pH=7:  P(DEPROT)={by_n[0]:.4f}  P(PROT)={by_n[1]:.4f}")
    print(f"  Effective pKa = {eff_pKa:.4f}")
    print(f"  Single-tautomer baseline (state 1 only): pKa = 3.271")
    print(f"  Δ(multistate − single)  = {eff_pKa - 3.271:+.4f}")


def main():
    test_as4_standalone()
    test_hip_standalone()
    test_bpti_as4_multistate()


if __name__ == "__main__":
    main()
