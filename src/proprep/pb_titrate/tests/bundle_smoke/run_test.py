#!/usr/bin/env python3
"""Bundle smoke test for cluster cutout machinery.

Goals:
  1. Bond walk on PRN-1011 reaches its parent heme (HCO).
  2. Cluster cutout (R=25 Å) gives reasonable atom count.
  3. strip_bulk_ions runs cleanly on the bundle and reports per-ion counts.
  4. Single PB call on the cutout runs without error.
"""
import sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent))

import parmed as pmd
from proprep.pb_titrate.ions import (strip_bulk_ions, within_distance_of_protein,
                            keep_none, keep_all_ions)
from proprep.pb_titrate.cluster import (bond_walk, residues_within, cluster_for_target,
                                slice_structure)
from proprep.pb_titrate.pb_backend import run_pbsa, total_electrostatic

BUNDLE_PRM = Path("/workhorse/9YUQ/SecondPass/MS001/transformed_microstate_001_fixed.prmtop")
BUNDLE_RST = Path("/workhorse/9YUQ/SecondPass/MS001/transformed_microstate_001.rst7")


def main():
    t0 = time.time()
    print(f"=== Bundle smoke test ===")

    print(f"\n--- Step 1: load bundle ({BUNDLE_PRM.name}) ---")
    s = pmd.load_file(str(BUNDLE_PRM), str(BUNDLE_RST))
    print(f"  raw: {len(s.residues)} residues, {len(s.atoms)} atoms, "
          f"net q={sum(a.charge for a in s.atoms):+.3f}")
    print(f"  load took {time.time()-t0:.0f}s")

    print(f"\n--- Step 2: strip_bulk_ions (default rule = within 4 Å of protein) ---")
    t1 = time.time()
    info = strip_bulk_ions(s, verbose=True)
    print(f"  post-strip: {len(s.residues)} residues, {len(s.atoms)} atoms, "
          f"net q={sum(a.charge for a in s.atoms):+.3f}")
    print(f"  strip took {time.time()-t1:.0f}s")
    print(f"  pre_counts (top 6 by count):")
    for k, v in sorted(info["pre_counts"].items(), key=lambda x: -x[1])[:6]:
        print(f"    {k}: {v}")
    print(f"  post_counts (top 6):")
    for k, v in sorted(info["post_counts"].items(), key=lambda x: -x[1])[:6]:
        print(f"    {k}: {v}")
    n_kept = len(info["kept_structural_ions"])
    print(f"  kept {n_kept} structural ions")
    if n_kept > 0:
        print(f"    examples: {info['kept_structural_ions'][:5]}")

    print(f"\n--- Step 3: bond walk PRN-1011 ---")
    target = next(r for r in s.residues
                   if r.name == "PRN" and (r.number + 1) == 1011)
    print(f"  target: {target.name}-{target.number+1} ({len(target.atoms)} atoms)")
    bw3 = bond_walk(target, depth=3)
    bw5 = bond_walk(target, depth=5)
    print(f"  bond_walk(depth=3): {len(bw3)} residues -> "
          f"{sorted({s.residues[i].name+'-'+str(s.residues[i].number+1) for i in bw3})}")
    print(f"  bond_walk(depth=5): {len(bw5)} residues -> "
          f"{sorted({s.residues[i].name+'-'+str(s.residues[i].number+1) for i in bw5})}")
    # Verify heme is reached
    bw3_names = {s.residues[i].name for i in bw3}
    if any(n in ("HCO", "HCR") for n in bw3_names):
        print(f"  ✓ parent heme reached at depth 3")
    else:
        print(f"  ✗ parent heme NOT reached at depth 3 (got names: {bw3_names})")

    print(f"\n--- Step 4: full cluster cutout R=25 Å ---")
    cl = cluster_for_target(s, target, R=25.0, bond_depth=3)
    print(f"  cluster_for_target(R=25): {len(cl)} residues")
    sub = slice_structure(s, cl)
    sub_q = sum(a.charge for a in sub.atoms)
    print(f"  sliced structure: {len(sub.atoms)} atoms, "
          f"{len(sub.residues)} residues, "
          f"net q={sub_q:+.3f}")
    if abs(sub_q - round(sub_q)) > 0.05:
        print(f"  NOTE: net charge is fractional. This is EXPECTED pre-integration:")
        print(f"        the heme is MCPB-parameterized as a multi-residue group")
        print(f"        (HCO + 2×PRN + DHO + PHO + (CYO)) with partial charges")
        print(f"        balanced ACROSS the group. A 25 Å cutout around one PRN")
        print(f"        captures the parent HCO but not necessarily the sibling")
        print(f"        PRN or axial residues, breaking the per-heme integer.")
        print(f"        Will be fixed by the Site abstraction.")
    print(f"  PRN-1011 in subset? "
          f"{any(r.name=='PRN' and r.number+1==1011 for r in sub.residues)}")

    print(f"  R=25 atoms in cluster: {len(sub.atoms)}")

    print(f"\n--- Step 5: PB call on the cluster ---")
    wd = HERE / "pbsa_smoke"
    wd.mkdir(parents=True, exist_ok=True)
    t2 = time.time()
    e = run_pbsa(sub, wd, "smoke",
                  epsin=4.0, space=0.5, nfocus=2, bcopt=5, istrng=650,
                  fillratio=2.0)
    print(f"  PB succeeded in {time.time()-t2:.1f}s")
    print(f"    EPB={e['EPB']:+.3f}  EELEC={e['EELEC']:+.3f}  "
          f"E_elec={total_electrostatic(e):+.3f} kcal/mol")

    print(f"\n  Total wall time: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
