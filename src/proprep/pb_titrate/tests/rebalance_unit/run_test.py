#!/usr/bin/env python3
"""Unit test for titrate.rebalance.rebalance_ions.

Reloads the bundle prmtop (484 Na+, 343 Cl-, 32 Ca2+, 152k waters) for
each case so cases are independent. Reload is ~22s.

NO strip is applied: the realistic production flow is
  full solvated prmtop -> apply state map -> rebalance bulk ions
not strip-then-rebalance. The 32 CryoEM Ca2+ are the only structural
ions; they are passed via keep_ions to be protected from removal.

  Case A: target = current → no-op.
  Case B: protein +5 → 5 Na+ removed.
  Case C: protein -3 → 3 Cl- removed.
  Case D: cation_preference forces Ca2+ first; structural Ca2+ skipped,
          falls back to Na+. Need ΔQ=-2 → 2 Na+ removed.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent))

import parmed as pmd

from proprep.pb_titrate.rebalance import rebalance_ions

BUNDLE_PRM = Path("/workhorse/9YUQ/SecondPass/MS001/transformed_microstate_001_fixed.prmtop")
BUNDLE_RST = Path("/workhorse/9YUQ/SecondPass/MS001/transformed_microstate_001.rst7")

# The 32 structural Ca2+ from CryoEM (residues 3352-3383, name "CA").
KEEP_IONS = {("CA", n) for n in range(3352, 3384)}


def perturb_net_q(structure: pmd.Structure, delta: float):
    """Bump the first protein residue's first-atom charge by `delta` to
    shift the system net q by exactly `delta`."""
    for r in structure.residues:
        if len(r.atoms) > 1 and r.name.strip() not in {"WAT", "HOH"}:
            r.atoms[0].charge += delta
            return
    raise RuntimeError("No protein residue found to perturb")


def assert_close(actual, expected, tol, label):
    if abs(actual - expected) > tol:
        print(f"  FAIL  {label}: got {actual:.4f}, expected {expected:.4f} (tol {tol})")
        return False
    print(f"  ok    {label}: {actual:.4f} ≈ {expected:.4f}")
    return True


def case_A_noop():
    print("\n--- Case A: target = current → no-op ---")
    s = pmd.load_file(str(BUNDLE_PRM), str(BUNDLE_RST))
    q0 = sum(a.charge for a in s.atoms)
    print(f"  initial net q = {q0:+.4f}")
    info = rebalance_ions(s, target_net_q=q0, keep_ions=KEEP_IONS)
    n_removed = sum(info["removed_per_species"].values())
    ok = assert_close(info["net_q_after"], q0, 0.01, "no-op net q unchanged")
    if n_removed != 0:
        print(f"  FAIL  expected 0 ions removed, got {n_removed}")
        ok = False
    else:
        print(f"  ok    removed: {info['removed_per_species']}")
    return ok


def case_B_remove_cations():
    print("\n--- Case B: ΔQ_protein = +5 → remove 5 Na+ ---")
    s = pmd.load_file(str(BUNDLE_PRM), str(BUNDLE_RST))
    q_target = sum(a.charge for a in s.atoms)
    perturb_net_q(s, +5.0)
    print(f"  net q target={q_target:+.4f}, after +5 perturb="
          f"{sum(a.charge for a in s.atoms):+.4f}")
    info = rebalance_ions(s, target_net_q=q_target,
                           keep_ions=KEEP_IONS, verbose=True)
    ok = assert_close(info["net_q_after"], q_target, 0.01, "rebalanced net q")
    n_na = info["removed_per_species"].get("Na+", 0)
    if n_na != 5:
        print(f"  FAIL  expected 5 Na+ removed, got {info['removed_per_species']}")
        ok = False
    else:
        print(f"  ok    removed exactly 5 Na+")
    n_ca = sum(1 for r in s.residues
                if r.name.strip() == "CA" and len(r.atoms) == 1)
    if n_ca != 32:
        print(f"  FAIL  expected 32 structural Ca2+ remain, got {n_ca}")
        ok = False
    else:
        print(f"  ok    32 structural Ca2+ untouched")
    return ok


def case_C_remove_anions():
    print("\n--- Case C: ΔQ_protein = -3 → remove 3 Cl- ---")
    s = pmd.load_file(str(BUNDLE_PRM), str(BUNDLE_RST))
    q_target = sum(a.charge for a in s.atoms)
    perturb_net_q(s, -3.0)
    print(f"  net q target={q_target:+.4f}, after -3 perturb="
          f"{sum(a.charge for a in s.atoms):+.4f}")
    info = rebalance_ions(s, target_net_q=q_target,
                           keep_ions=KEEP_IONS, verbose=True)
    ok = assert_close(info["net_q_after"], q_target, 0.01, "rebalanced net q")
    n_cl = info["removed_per_species"].get("Cl-", 0)
    if n_cl != 3:
        print(f"  FAIL  expected 3 Cl- removed, got {info['removed_per_species']}")
        ok = False
    else:
        print(f"  ok    removed exactly 3 Cl-")
    return ok


def case_D_structural_skipped():
    print("\n--- Case D: Ca2+ first in pref, but structural → fall back to Na+ ---")
    s = pmd.load_file(str(BUNDLE_PRM), str(BUNDLE_RST))
    q_target = sum(a.charge for a in s.atoms)
    perturb_net_q(s, +2.0)
    info = rebalance_ions(s, target_net_q=q_target,
                           keep_ions=KEEP_IONS,
                           cation_preference=["CA", "Mg2+", "Na+", "K+"],
                           verbose=True)
    ok = assert_close(info["net_q_after"], q_target, 0.01, "rebalanced net q")
    n_ca_removed = info["removed_per_species"].get("CA", 0)
    if n_ca_removed != 0:
        print(f"  FAIL  expected 0 Ca2+ removed (structural), got {n_ca_removed}")
        ok = False
    else:
        print(f"  ok    no structural Ca2+ removed")
    n_na = info["removed_per_species"].get("Na+", 0)
    if n_na != 2:
        print(f"  FAIL  expected 2 Na+ removed (fallback), got {info['removed_per_species']}")
        ok = False
    else:
        print(f"  ok    2 Na+ removed as fallback")
    n_ca = sum(1 for r in s.residues
                if r.name.strip() == "CA" and len(r.atoms) == 1)
    if n_ca != 32:
        print(f"  FAIL  expected 32 structural Ca2+ remain, got {n_ca}")
        ok = False
    else:
        print(f"  ok    all 32 Ca2+ still in structure")
    return ok


def main():
    print("=== rebalance_ions unit test ===")
    print(f"Bundle: {BUNDLE_PRM.name}  (reload per case for isolation)")
    cases = [
        ("A", case_A_noop),
        ("B", case_B_remove_cations),
        ("C", case_C_remove_anions),
        ("D", case_D_structural_skipped),
    ]
    results = []
    for name, fn in cases:
        try:
            results.append(fn())
        except Exception as e:
            print(f"  FAIL  Case {name} raised: {type(e).__name__}: {e}")
            results.append(False)
    n_pass = sum(results)
    n_fail = len(results) - n_pass
    print(f"\n=== summary: {n_pass}/{len(results)} cases passed ===")
    if n_fail:
        sys.exit(1)
    print("*** ALL REBALANCE CASES PASSED ***")


if __name__ == "__main__":
    main()
