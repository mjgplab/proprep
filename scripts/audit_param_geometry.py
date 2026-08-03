#!/usr/bin/env python3
"""Cross-check shipped parameter sets against their own library coordinates.

Every specialized_residues leaf ships a .frcmod saying what geometry the force
field wants and a .lib holding actual coordinates of the same molecule. Nothing
compares them. This does.

It exists because the his_met_axial_c_type heme frcmods restrained the two TRANS
pyrrole N-Fe-N angles to 90 degrees. NP/NO alternate around the macrocycle, so a
same-type triplet sits across the ring; all three N-Fe-N type combinations had
been written with the cis class average. The set's own library coordinates said
those angles were ~176 degrees. A single comparison would have caught it.

Two checks:

  ANGLES   For each frcmod ANGLE term, measure every instance of that atom-type
           triplet in the library and confirm theta0 falls inside the observed
           range. Terms are type-averaged, so a range test is the honest one --
           testing against a fixed "reasonable" value would invent a threshold
           the physics does not supply. Terms with k == 0 are skipped: conste
           writes unrestrained trans angles as k=0 with a placeholder theta0,
           and flagging those is pure noise.

  CONTACTS Atom pairs closer than they should be, classified by topological
           separation. Blanket distance cutoffs do not work here: a geminal
           (1-3) H..H at 1.56 A is perfectly normal when C-H is 0.970 A, so
           1-3 pairs are judged by their ANGLE, not their distance. Expect the
           non-bonded check to surface intramolecular hydrogen bonds (NAD/NADP
           show O..H near 1.67 A); those are short contacts, not clashes.

Not a product feature and not wired into the workflow -- a developer audit tool.
Findings still need judgement: some flagged geometry is inherited verbatim from
upstream AmberTools (conste.lib's thioether hydrogens, for one) and is not
ProPrep's to change.

Usage:
    python scripts/audit_param_geometry.py [root] [--angle-tol DEG] [--verbose]
"""

import argparse
import itertools
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

DEFAULT_ROOT = Path("src/proprep/forcefield_params/specialized_residues")

# theta0 may sit this far outside the observed range before we complain. Absorbs
# crystal-vs-optimised differences and the spread inside a type class. The defect
# that motivated this tool missed by 86 degrees; legitimate spread ran under 20.
ANGLE_TOL = 25.0

# A 1-3 pair is judged by its angle, never its distance -- and only when no
# metal is involved. Metal coordination legitimately spans a huge range: Fe-S-Fe
# in an Fe4S4 cubane is ~75 deg and cis N-Fe-N is ~84 deg, both correct. Applying
# a light-atom floor to those reports valid inorganic geometry as broken. Around
# saturated C/N/O/S, by contrast, anything far below tetrahedral is real damage.
MIN_13_ANGLE = 85.0

# Anything heavier than Ar, plus the common biological light metals. Used only to
# exempt triplets from the 1-3 angle floor.
METALS = {11, 12, 19, 20} | set(range(21, 84))

# True non-bonded contacts (1-4 and beyond).
MIN_NONBONDED = 1.7


def parse_lib(path):
    """OFF library -> {unit: {names, types, pos, bonds}}."""
    units = defaultdict(lambda: {"names": [], "types": [], "pos": [], "bonds": [],
                                 "elem": []})
    unit = kind = None
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("!"):
            m = re.match(r"!entry\.(\w+)\.unit\.(\S+)", line)
            unit, kind = (m.group(1), m.group(2)) if m else (None, None)
            continue
        if unit is None:
            continue
        if kind == "atoms":
            m = re.match(r'\s*"(\S+)"\s+"(\S+)"'
                         r'\s+\S+\s+\S+\s+\S+\s+\S+\s+(-?\d+)', line)
            if m:
                units[unit]["names"].append(m.group(1))
                units[unit]["types"].append(m.group(2))
                units[unit]["elem"].append(int(m.group(3)))
        elif kind == "positions":
            f = line.split()
            if len(f) == 3:
                try:
                    units[unit]["pos"].append(tuple(map(float, f)))
                except ValueError:
                    pass
        elif kind == "connectivity":
            f = line.split()
            if len(f) == 3:
                units[unit]["bonds"].append((int(f[0]) - 1, int(f[1]) - 1))
    return units


def parse_angles(path):
    """frcmod -> {(t1, t2, t3): (k, theta0)}."""
    out, section = {}, None
    for line in path.read_text(errors="replace").splitlines():
        head = line.strip().upper()
        if head in ("MASS", "BOND", "ANGLE", "ANGL", "DIHEDRAL", "DIHE",
                    "IMPROPER", "IMPR", "NONBON", "NONB"):
            section = "ANGLE" if head in ("ANGLE", "ANGL") else head
            continue
        if section != "ANGLE" or not line.strip():
            continue
        m = re.match(r"^(\S{1,2})-\s?(\S{1,2})-\s?(\S{1,2})\s+([\d.]+)\s+([\d.]+)",
                     line.split("!")[0])
        if m:
            out[(m.group(1), m.group(2), m.group(3))] = (float(m.group(4)),
                                                         float(m.group(5)))
    return out


def angle(p, q, r):
    v1 = [p[i] - q[i] for i in range(3)]
    v2 = [r[i] - q[i] for i in range(3)]
    n1 = math.sqrt(sum(x * x for x in v1))
    n2 = math.sqrt(sum(x * x for x in v2))
    if n1 == 0 or n2 == 0:
        return None
    dot = sum(a * b for a, b in zip(v1, v2))
    return math.degrees(math.acos(max(-1.0, min(1.0, dot / (n1 * n2)))))


def dist(p, q):
    return math.sqrt(sum((p[i] - q[i]) ** 2 for i in range(3)))


def check_leaf(frcmod, angle_tol):
    """Return (angle_findings, contact_findings) for one frcmod + sibling libs."""
    libs = sorted(frcmod.parent.glob("*.lib"))
    if not libs:
        return [], []
    terms = parse_angles(frcmod)
    measured = defaultdict(list)
    contacts = []

    for lib in libs[:1]:                    # sibling libs share one geometry
        for uname, u in parse_lib(lib).items():
            if not u["pos"] or len(u["pos"]) != len(u["names"]):
                continue
            adj = defaultdict(set)
            for a, b in u["bonds"]:
                adj[a].add(b)
                adj[b].add(a)

            for centre, nbrs in adj.items():
                for a, b in itertools.combinations(sorted(nbrs), 2):
                    val = angle(u["pos"][a], u["pos"][centre], u["pos"][b])
                    if val is None:
                        continue
                    trip = (u["types"][a], u["types"][centre], u["types"][b])
                    key = trip if trip in terms else trip[::-1]
                    if key in terms:
                        measured[key].append(val)
                    elem = u.get("elem") or []
                    has_metal = any(elem[x] in METALS
                                    for x in (a, centre, b) if x < len(elem))
                    if val < MIN_13_ANGLE and not has_metal:
                        contacts.append((
                            f"{uname}: {u['names'][a]}-{u['names'][centre]}-{u['names'][b]}",
                            f"1-3 angle {val:.1f} deg", val))

            bonded = {frozenset(b) for b in u["bonds"]}
            for i, j in itertools.combinations(range(len(u["pos"])), 2):
                if frozenset((i, j)) in bonded or (adj[i] & adj[j]):
                    continue                # 1-2 and 1-3 handled above
                d = dist(u["pos"][i], u["pos"][j])
                if d < MIN_NONBONDED:
                    contacts.append((f"{uname}: {u['names'][i]}..{u['names'][j]}",
                                     f"non-bonded {d:.3f} A", d))

    findings = []
    for key, (k, theta0) in sorted(terms.items()):
        if k == 0 or key not in measured:
            continue
        lo, hi = min(measured[key]), max(measured[key])
        if lo - angle_tol <= theta0 <= hi + angle_tol:
            continue
        findings.append({
            "term": "-".join(key), "k": k, "theta0": theta0,
            "lo": lo, "hi": hi, "n": len(measured[key]),
            "off": max(abs(theta0 - lo), abs(theta0 - hi)),
        })
    findings.sort(key=lambda f: -f["off"])
    contacts.sort(key=lambda c: c[2])
    return findings, contacts


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", nargs="?", default=str(DEFAULT_ROOT))
    ap.add_argument("--angle-tol", type=float, default=ANGLE_TOL)
    ap.add_argument("--verbose", action="store_true",
                    help="also list contact findings")
    args = ap.parse_args()

    root = Path(args.root)
    if not root.exists():
        sys.exit(f"no such path: {root}")

    n_leaves = n_flagged = 0
    worst = []
    for frcmod in sorted(root.rglob("*.frcmod")):
        n_leaves += 1
        findings, contacts = check_leaf(frcmod, args.angle_tol)
        if not findings and not (args.verbose and contacts):
            continue
        n_flagged += 1
        print(f"\n{frcmod.relative_to(root)}")
        for f in findings:
            print(f"   {f['term']:<14} k={f['k']:8.2f} theta0={f['theta0']:8.3f}"
                  f"   library geometry {f['lo']:.2f}..{f['hi']:.2f} deg"
                  f" over {f['n']} instance(s)   [off by {f['off']:.0f} deg]")
            worst.append((f["off"], str(frcmod.relative_to(root)), f["term"]))
        if args.verbose:
            for where, what, _ in contacts[:8]:
                print(f"   contact  {where:<34} {what}")
            if len(contacts) > 8:
                print(f"   contact  ... and {len(contacts) - 8} more")

    print(f"\n{'-' * 70}")
    print(f"{n_leaves} frcmod file(s) checked, {n_flagged} with findings"
          f"  (angle tolerance {args.angle_tol:.0f} deg)")
    if worst:
        worst.sort(reverse=True)
        print(f"largest discrepancy: {worst[0][0]:.0f} deg  "
              f"{worst[0][2]} in {worst[0][1]}")
    return 1 if worst else 0


if __name__ == "__main__":
    sys.exit(main())
