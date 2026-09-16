"""
Would a bonded term written under shared atom types reach other residues?

A modified residue typed by antechamber (``-at amber``) carries the same atom
types as the standard residues. Anything ProPrep writes under those types is
loaded after the protein leaprc and applies to every residue in the system,
so a refined value is safe only when nothing else in the system uses that
exact type combination. Two facts decide that for a dihedral:

* the force field has a SPECIFIC term for the quad (no ``X`` wildcard):
  a refit would replace it everywhere it is used;
* a standard residue contains the exact quad, within the residue or across
  a peptide bond: a refit would parameterize that residue too.

A quad the force field covers only by a wildcard, and that no standard
residue contains, is local to the residues that carry it, which for a
modified amino acid means the modified residue itself. That is the usual
situation for a covalent linkage torsion (``CT-S-CT-CT`` with ff14SB, where
methionine and cysteine use ``2C`` for their sulfur-bonded carbons).
"""

from __future__ import annotations

import os
import re
import sys
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

Quad = Tuple[str, str, str, str]

# Protein force-field files whose SPECIFIC dihedral terms a modified residue
# could collide with, and the residue libraries tLEaP loads for them.
_PARM_FILES = ("parm10.dat", "parm19.dat", "frcmod.ff14SB", "frcmod.ff19SB", "frcmod.ff99SB")
_PROTEIN_LIBS = ("amino12.lib", "aminont12.lib", "aminoct12.lib", "amino19.lib",
                 "aminont19.lib", "aminoct19.lib")

_DIHE_LINE = re.compile(r"^(.{2})-(.{2})-(.{2})-(.{2})\s+\d+\s+[-+\d.]+\s+[-+\d.]+\s+[-+\d.]+")


def canonical(quad: Sequence[str]) -> Quad:
    q = tuple(t.strip() for t in quad)
    r = tuple(reversed(q))
    return q if q <= r else r  # type: ignore[return-value]


def _leap_dir() -> Optional[Path]:
    for root in (os.environ.get("AMBERHOME"), sys.prefix):
        if root and (Path(root) / "dat" / "leap").is_dir():
            return Path(root) / "dat" / "leap"
    return None


# ---------------------------------------------------------------------------
# specific terms in the parameter files
# ---------------------------------------------------------------------------

def parm_specific_dihedrals(paths: Iterable[str]) -> Set[Quad]:
    """Type quads with a specific (no ``X``) dihedral line in the given files.

    Works on both parm ``.dat`` files and frcmods: a dihedral line is the
    only bonded line shaped ``A -B -C -D  IDIVF  PK  PHASE  PN`` (an improper
    has no IDIVF, bonds and angles have shorter keys).
    """
    found: Set[Quad] = set()
    for path in paths:
        try:
            with open(path, errors="ignore") as fh:
                for line in fh:
                    m = _DIHE_LINE.match(line)
                    if not m:
                        continue
                    types = tuple(t.strip() for t in m.groups())
                    if all(types) and "X" not in types:
                        found.add(canonical(types))
        except OSError:
            continue
    return found


# ---------------------------------------------------------------------------
# type quads present in the standard residues
# ---------------------------------------------------------------------------

def _parse_off_units(lib_path: str) -> Dict[str, dict]:
    """``{unit: {"types": [..], "bonds": [(i, j), ..], "head": i, "tail": j}}`` (1-based)."""
    units: Dict[str, dict] = {}
    try:
        text = Path(lib_path).read_text(errors="ignore")
    except OSError:
        return units
    unit, table = None, None
    for line in text.splitlines():
        if line.startswith("!entry."):
            m = re.match(r"!entry\.(.+?)\.unit\.(\w+)", line)
            unit, table = (m.group(1), m.group(2)) if m else (None, None)
            if unit is not None:
                units.setdefault(unit, {"types": [], "bonds": [], "head": 0, "tail": 0, "_connect": []})
            continue
        if line.startswith("!") or unit is None or table is None:
            continue
        parts = line.split()
        if table == "atoms" and len(parts) >= 2:
            units[unit]["types"].append(parts[1].strip('"'))
        elif table == "connectivity" and len(parts) >= 2:
            try:
                units[unit]["bonds"].append((int(parts[0]), int(parts[1])))
            except ValueError:
                pass
        elif table == "connect" and parts:
            try:
                units[unit]["_connect"].append(int(parts[0]))
            except ValueError:
                pass
    for u in units.values():
        c = u.pop("_connect")
        if len(c) >= 2:
            u["head"], u["tail"] = c[0], c[1]
    return {k: v for k, v in units.items() if v["types"]}


def _adjacency(bonds: Iterable[Tuple[int, int]]) -> Dict[int, Set[int]]:
    adj: Dict[int, Set[int]] = {}
    for i, j in bonds:
        adj.setdefault(i, set()).add(j)
        adj.setdefault(j, set()).add(i)
    return adj


def standard_residue_type_quads(lib_paths: Iterable[str]) -> Set[Quad]:
    """Every dihedral type quad a chain of these residues can contain.

    Within each unit, every path of four bonded atoms. Across a peptide bond,
    every unit's tail joined to every unit's head, so backbone quads that
    span two residues (``2C-XC-C-N``, ``C-N-XC-2C`` ...) are included.
    """
    quads: Set[Quad] = set()
    all_units: Dict[str, dict] = {}
    for path in lib_paths:
        all_units.update(_parse_off_units(path))
    for u in all_units.values():
        types, adj = u["types"], _adjacency(u["bonds"])
        u["_adj"] = adj
        for b, nb in adj.items():
            for c in nb:
                for a in adj.get(b, ()):
                    if a == c:
                        continue
                    for d in adj.get(c, ()):
                        if d == b or d == a:
                            continue
                        try:
                            quads.add(canonical((types[a - 1], types[b - 1], types[c - 1], types[d - 1])))
                        except IndexError:
                            pass
    # cross-residue: tail(U) - head(V)
    tails = [(u["types"], u["_adj"], u["tail"]) for u in all_units.values() if u["tail"]]
    heads = [(u["types"], u["_adj"], u["head"]) for u in all_units.values() if u["head"]]
    seen_pairs: Set[Tuple[Tuple[str, ...], Tuple[str, ...]]] = set()
    for t_types, t_adj, t in tails:
        t_type = t_types[t - 1]
        t_nb = [(t_types[b - 1], [t_types[a - 1] for a in t_adj.get(b, ()) if a != t]) for b in t_adj.get(t, ())]
        for h_types, h_adj, h in heads:
            h_type = h_types[h - 1]
            h_nb = [(h_types[c - 1], [h_types[d - 1] for d in h_adj.get(c, ()) if d != h]) for c in h_adj.get(h, ())]
            key = (tuple(sorted(str(x) for x in t_nb)) + (t_type,), (h_type,) + tuple(sorted(str(x) for x in h_nb)))
            if key in seen_pairs:
                continue
            seen_pairs.add(key)
            for b_type, a_types in t_nb:
                for a_type in a_types:
                    quads.add(canonical((a_type, b_type, t_type, h_type)))
                for c_type, _ in h_nb:
                    quads.add(canonical((b_type, t_type, h_type, c_type)))
            for c_type, d_types in h_nb:
                for d_type in d_types:
                    quads.add(canonical((t_type, h_type, c_type, d_type)))
    return quads


# ---------------------------------------------------------------------------
# the question a refinement path asks
# ---------------------------------------------------------------------------

@lru_cache(maxsize=1)
def _protein_ff_facts() -> Tuple[Set[Quad], Set[Quad]]:
    leap = _leap_dir()
    if leap is None:
        return set(), set()
    parm = [str(leap / "parm" / f) for f in _PARM_FILES if (leap / "parm" / f).is_file()]
    libs = [str(leap / "lib" / f) for f in _PROTEIN_LIBS if (leap / "lib" / f).is_file()]
    return parm_specific_dihedrals(parm), standard_residue_type_quads(libs)


def reset_cache() -> None:
    _protein_ff_facts.cache_clear()


def shared_type_dihedral_conflict(quad: Sequence[str]) -> Optional[str]:
    """Why a refit of this type quad would reach beyond the modified residue, or None.

    None also when the Amber data cannot be found: the caller then has no
    basis to refuse, and says so.
    """
    specific, standard = _protein_ff_facts()
    q = canonical(quad)
    name = "-".join(q)
    if q in specific:
        return (f"the protein force field defines {name} explicitly; a refit written under "
                "these shared atom types would replace that definition everywhere it is used")
    if q in standard:
        return (f"a standard residue contains the type quad {name}; a refit written under "
                "these shared atom types would parameterize that residue too")
    return None


def amber_data_available() -> bool:
    specific, standard = _protein_ff_facts()
    return bool(specific) and bool(standard)
