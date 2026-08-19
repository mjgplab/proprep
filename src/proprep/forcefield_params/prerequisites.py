"""
Infer which leaprcs a parameter set needs, from the atom types it uses.

A deposited set is only loadable if the force fields defining its atom types
are sourced. That requirement is declarable as ``prerequisites.leaprc_groups``
in metadata, and the Topology Generator already enforces what it finds there --
but almost nothing writes it. MCPB output declared nothing, and neither did an
imported set, so the requirement existed only in the depositor's head.

It does not need to. At deposit time the types are in the files: a GAFF set
uses ``c3``/``os``/``p5``, an MCPB set uses ``N``/``CT``/``XC`` alongside its
``M*``/``Y*``. Reading them says what must be sourced, with no guessing and no
prose about what the molecule is.

Groups are AND-of-OR, matching the schema: every group must be satisfied, and
any member satisfies its group.
"""

from __future__ import annotations

import logging
import os
import re
import sys
from pathlib import Path
from typing import Iterable, List, Sequence, Set

logger = logging.getLogger(__name__)

# Any of these satisfies "a protein force field is sourced". Mirrors the group
# zinc/cys4 declares by hand, which is the only hand-written example in tree.
PROTEIN_LEAPRCS = [
    "leaprc.protein.ff19SB",
    "leaprc.protein.ff14SB",
    "leaprc.protein.ff14SBonlysc",
    "leaprc.constph",
    "leaprc.conste",
]
GAFF2_LEAPRC = "leaprc.gaff2"


def _parm_dir() -> Path | None:
    amberhome = os.environ.get("AMBERHOME")
    for candidate in ([Path(amberhome)] if amberhome else []) + [Path(sys.prefix)]:
        parm = candidate / "dat" / "leap" / "parm"
        if parm.is_dir():
            return parm
    return None


def _mass_types(path: Path, labelled: bool) -> Set[str]:
    """Atom types in a file's MASS block.

    A .dat opens with its MASS block; a frcmod labels it. Either way the block
    ends at the first blank line.
    """
    types: Set[str] = set()
    try:
        with open(path, errors="ignore") as handle:
            in_mass = not labelled
            if not labelled:
                next(handle, None)          # title line
            for line in handle:
                stripped = line.strip()
                if labelled and not in_mass:
                    if stripped.upper() == "MASS":
                        in_mass = True
                    continue
                if not stripped:
                    break
                parts = stripped.split()
                if len(parts) >= 2:
                    try:
                        float(parts[1])
                    except ValueError:
                        continue
                    types.add(parts[0])
    except OSError:
        pass
    return types


def type_sources() -> tuple[Set[str], Set[str]]:
    """``(gaff2_only, protein_only)`` atom types.

    Types defined by both are excluded from each: they discriminate nothing, so
    letting them vote would make every set claim to need everything.
    """
    parm = _parm_dir()
    if parm is None:
        return set(), set()

    gaff: Set[str] = set()
    for name in ("gaff2.dat", "gaff.dat"):
        gaff |= _mass_types(parm / name, labelled=False)

    protein: Set[str] = set()
    for name in ("parm19.dat", "parm10.dat"):
        protein |= _mass_types(parm / name, labelled=False)
    # ff19SB's 2C/3C live in frcmod.ff19SB, not parm19.dat. Only frcmods a
    # leaprc actually sources are consulted; scanning all of them lets niche
    # sets vouch for types they merely happen to reuse.
    for pattern in ("frcmod.ff*", "frcmod.tip*", "frcmod.opc*", "frcmod.spce*"):
        for path in sorted(parm.glob(pattern)):
            protein |= _mass_types(path, labelled=True)

    return gaff - protein, protein - gaff


_LIB_ATOM_RE = re.compile(r'\s*"([^"]+)"\s+"([^"]+)"')


def atom_types_used(paths: Iterable[str | Path]) -> Set[str]:
    """Atom types a set's library and frcmod files actually use."""
    used: Set[str] = set()
    for raw in paths:
        path = Path(raw)
        if not path.is_file():
            continue
        try:
            text = path.read_text(errors="ignore")
        except OSError:
            continue

        if path.suffix.lower() in (".lib", ".off"):
            in_atoms = False
            for line in text.splitlines():
                if line.startswith("!entry.") and ".unit.atoms table" in line:
                    in_atoms = True
                    continue
                if in_atoms:
                    if line.startswith("!"):
                        in_atoms = False
                        continue
                    match = _LIB_ATOM_RE.match(line)
                    if match:
                        used.add(match.group(2).strip())
        elif path.suffix.lower() == ".mol2":
            in_atoms = False
            for line in text.splitlines():
                if line.startswith("@<TRIPOS>ATOM"):
                    in_atoms = True
                    continue
                if in_atoms:
                    if line.startswith("@<TRIPOS>"):
                        break
                    parts = line.split()
                    if len(parts) >= 6:
                        used.add(parts[5].strip())
    return used


def infer_leaprc_groups(paths: Sequence[str | Path]) -> List[List[str]]:
    """AND-groups of leaprcs the given parameter files require.

    Returns ``[]`` when nothing can be determined -- an empty list means "no
    prerequisites declared", which is what the schema already means and what
    every consumer already handles.
    """
    used = atom_types_used(paths)
    if not used:
        return []

    gaff_only, protein_only = type_sources()
    if not gaff_only and not protein_only:
        return []

    groups: List[List[str]] = []
    if used & protein_only:
        groups.append(list(PROTEIN_LEAPRCS))
    if used & gaff_only:
        groups.append([GAFF2_LEAPRC])
    return groups
