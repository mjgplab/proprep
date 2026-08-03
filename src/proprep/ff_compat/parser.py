"""
Parse a bundled or user FF set (frcmod + lib(s)) into a normalized
``FFSignature`` capturing every parameter that could collide with
another set's parameters under tleap's global-atom-type-namespace rule.

The signature is intentionally narrow: it records only what tleap uses
to resolve a parameter lookup at load time. Comments, ordering of
entries within a section, and whitespace formatting are all discarded.

Two same-key entries with identical numeric values are NOT flagged as
conflicts by the comparator (very common across sets that both inherit
ff14SB's standard types). Only differing values are conflicts.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union


# Section headers we recognize in a frcmod. AMBER variants:
#   ANGL / ANGLE both used
#   NONB / NONBON both used
_FRCMOD_SECTIONS = {"MASS", "BOND", "ANGL", "ANGLE", "DIHE", "IMPR", "NONB", "NONBON"}

# Number of trailing numeric values per section line
_TRAILING_VALUES = {
    "MASS":   1,   # type  mass  [polarizability]   comment
    "BOND":   2,   # type-type  K  r0  comment
    "ANGL":   2,   # type-type-type  K  theta0  comment
    "DIHE":   4,   # type-type-type-type  IDIVF  Vn  phase  period  comment
    "IMPR":   3,   # type-type-type-type  Vn  phase  period  comment
    "NONB":   2,   # type  R*  epsilon  comment
}


# ---------------------------------------------------------------------------
# Type-name canonicalization helpers
# ---------------------------------------------------------------------------

def _canon_bond(t1: str, t2: str) -> Tuple[str, str]:
    """Bonds: sort the two endpoints alphabetically."""
    return tuple(sorted([t1, t2]))


def _canon_angle(t1: str, t2: str, t3: str) -> Tuple[str, str, str]:
    """Angles: central atom (t2) is fixed; endpoints sorted by alphabetic order.
    If t1 > t3, reverse so t1 <= t3."""
    if t1 > t3:
        t1, t3 = t3, t1
    return (t1, t2, t3)


def _canon_dihedral(t1: str, t2: str, t3: str, t4: str) -> Tuple[str, str, str, str]:
    """Dihedrals: AMBER canonical form has the lexicographically smaller end
    on the left. If (t1, t2) > (t4, t3) lexically (comparing as reversed
    pairs), reverse the whole quad."""
    forward = (t1, t2, t3, t4)
    reverse = (t4, t3, t2, t1)
    return min(forward, reverse)


def _canon_improper(t1: str, t2: str, t3: str, t4: str) -> Tuple[str, str, str, str]:
    """Impropers: AMBER convention has the central atom in position 3.
    The three substituents (t1, t2, t4) are sorted alphabetically; t3 stays.
    If wildcards (X) are present, they sort after letters by Python default."""
    subs = sorted([t1, t2, t4])
    return (subs[0], subs[1], t3, subs[2])


# ---------------------------------------------------------------------------
# Signature container
# ---------------------------------------------------------------------------

# Atomic-number -> element symbol. Covers the elements that actually appear in
# the bundled FF libs (C, H, N, O, P, S, halogens, common metals); the catch-all
# branch returns None and the caller decides how to handle the gap.
_ELEMENT_SYMBOLS = {
    1: "H",   5: "B",   6: "C",   7: "N",   8: "O",   9: "F",
    11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P",  16: "S",
    17: "Cl", 19: "K",  20: "Ca", 24: "Cr", 25: "Mn", 26: "Fe",
    27: "Co", 28: "Ni", 29: "Cu", 30: "Zn", 33: "As", 34: "Se",
    35: "Br", 42: "Mo", 47: "Ag", 48: "Cd", 53: "I",  74: "W",
    78: "Pt", 79: "Au", 80: "Hg",
}


def _hybridization_for(element: int, degree: int) -> str:
    """Map (atomic number, bond degree) to a tleap addAtomTypes hybridization label.

    tleap accepts "sp", "sp2", "sp3" (and ignores "?" / unknown gracefully).
    Connectivity-based rules — the only ones we use, since they don't depend on
    the atom-type name string:

        H            → sp3       (convention; H has no hybridization but tleap accepts sp3)
        C: deg 4 → sp3, 3 → sp2, ≤2 → sp
        N: deg 4 → sp3, 3 → sp3, 2 → sp2, 1 → sp
        O: deg 2 → sp3, 1 → sp2
        S: deg ≥2 → sp3, 1 → sp2
        P: any deg → sp3
        Halogens (F, Cl, Br, I): sp3
        Group-1/2 metals, transition metals, metalloids: sp3 by convention
        Anything else: sp3 (safe fallback — tleap only uses this to choose a default radius)
    """
    if element == 1:
        return "sp3"
    if element == 6:  # C
        if degree >= 4:
            return "sp3"
        if degree == 3:
            return "sp2"
        return "sp"
    if element == 7:  # N
        if degree >= 4:
            return "sp3"
        if degree == 3:
            return "sp3"
        if degree == 2:
            return "sp2"
        return "sp"
    if element == 8:  # O
        if degree >= 2:
            return "sp3"
        return "sp2"
    if element == 16:  # S
        if degree >= 2:
            return "sp3"
        return "sp2"
    # Phosphorus, halogens, metals, anything else: sp3 fallback.
    return "sp3"


@dataclass
class FFSignature:
    """Normalized parameter signature for a single FF set."""
    set_id: str
    source_files: List[str] = field(default_factory=list)

    # Atomic-level parameters
    mass: Dict[str, float] = field(default_factory=dict)
    nonb: Dict[str, Tuple[float, float]] = field(default_factory=dict)

    # Bonded parameters (keys are canonical type tuples)
    bonds: Dict[Tuple[str, str], Tuple[float, float]] = field(default_factory=dict)
    angles: Dict[Tuple[str, str, str], Tuple[float, float]] = field(default_factory=dict)
    # Dihedrals: keyed by (canon_quad, periodicity_int). Multiple periodicities
    # for the same quad each get their own key.
    dihedrals: Dict[Tuple[Tuple[str, str, str, str], int], Tuple[float, float, float]] = field(default_factory=dict)
    # Impropers: keyed by canon_quad + periodicity_int as well.
    impropers: Dict[Tuple[Tuple[str, str, str, str], int], Tuple[float, float]] = field(default_factory=dict)

    # Atom-type declarations from addAtomTypes lines (rarely in frcmod;
    # usually inferred from lib files or tleap input)
    atom_type_declarations: Dict[str, Tuple[str, str]] = field(default_factory=dict)

    # Lib-derived assignments: (unit_name, atom_name) -> type_name
    library_unit_types: Dict[Tuple[str, str], str] = field(default_factory=dict)
    # Lib-derived (unit_name, atom_name) -> atomic number (column 7 of the atoms table).
    library_unit_elements: Dict[Tuple[str, str], int] = field(default_factory=dict)
    # Lib-derived (unit_name, atom_name) -> bond degree (count from the connectivity table).
    library_unit_degrees: Dict[Tuple[str, str], int] = field(default_factory=dict)

    def all_types(self) -> set:
        """Every atom type this signature mentions, anywhere."""
        types = set(self.mass) | set(self.nonb) | set(self.atom_type_declarations)
        for pair in self.bonds:
            types.update(pair)
        for tri in self.angles:
            types.update(tri)
        for (quad, _period) in self.dihedrals:
            types.update(quad)
        for (quad, _period) in self.impropers:
            types.update(quad)
        types.update(self.library_unit_types.values())
        # Wildcards are not "real" types; exclude
        types.discard("X")
        return types

    def infer_type_declaration(self, type_name: str) -> Optional[Tuple[str, str]]:
        """Return ``(element_symbol, hybridization)`` for ``type_name`` if any
        lib atom uses it; otherwise None.

        Strategy: find every (unit, atom_name) bearing this type, look up
        its element + degree, derive hybridization, and pick the most
        common (element_symbol, hybridization) tuple. Returns None when
        the type is unused in any parsed lib or when the element is
        outside the known periodic table table.
        """
        from collections import Counter

        votes: Counter = Counter()
        for key, t in self.library_unit_types.items():
            if t != type_name:
                continue
            element = self.library_unit_elements.get(key)
            degree = self.library_unit_degrees.get(key, 0)
            if element is None:
                continue
            symbol = _ELEMENT_SYMBOLS.get(element)
            if symbol is None:
                continue
            votes[(symbol, _hybridization_for(element, degree))] += 1
        if not votes:
            return None
        # Most common wins; ties broken by insertion order via Counter semantics.
        return votes.most_common(1)[0][0]


# ---------------------------------------------------------------------------
# frcmod parsing
# ---------------------------------------------------------------------------

def _split_type_spec(text: str, n_types: int) -> List[str]:
    """Split a hyphen-separated type spec into N stripped tokens.

    AMBER frcmods use fixed 2-char columns for type names; a 1-char type
    (e.g., ``C``) is padded with a trailing space. Joining by '-' then
    stripping handles both forms. Caller guarantees text contains exactly
    n_types-1 hyphens between meaningful tokens.
    """
    parts = text.split("-")
    if len(parts) != n_types:
        raise ValueError(
            f"Expected {n_types} hyphen-separated types, got {len(parts)} in {text!r}"
        )
    return [p.strip() for p in parts]


def _parse_frcmod_line(line: str, section: str) -> Optional[Tuple]:
    """Parse a single section line. Returns (key, values) or None for blank/comment."""
    stripped = line.rstrip()
    if not stripped or stripped.lstrip().startswith(("#", "REMARK")):
        return None

    n_trailing = _TRAILING_VALUES[section]
    tokens = stripped.split()
    if len(tokens) < 1 + n_trailing:
        return None  # malformed line; skip silently

    # MASS and NONB have a single-token atom type at the head.
    if section in ("MASS", "NONB"):
        atom_type = tokens[0]
        try:
            values = [float(t) for t in tokens[1:1 + n_trailing]]
        except ValueError:
            return None
        return (atom_type, values)

    # BOND/ANGL/DIHE/IMPR: the type spec may contain a space when a 1-char
    # type is padded (e.g., "C -CX-CT-YZ"). Find the type spec by looking
    # for the first token containing a hyphen, then merge with earlier
    # tokens if needed.
    type_spec_tokens = []
    val_start = 0
    for i, tok in enumerate(tokens):
        type_spec_tokens.append(tok)
        # The type spec ends when we have a token that doesn't end with '-'
        # AND there's a numeric token ahead.
        if "-" in tok and not tok.endswith("-"):
            # Could this be the full spec? It is if remaining tokens are numeric.
            try:
                _ = float(tokens[i + 1])
                val_start = i + 1
                break
            except (IndexError, ValueError):
                continue
    if val_start == 0:
        # Didn't find a numeric break; try assuming the first token IS the full type-spec
        type_spec_tokens = [tokens[0]]
        val_start = 1

    type_spec = "".join(type_spec_tokens)
    n_expected_types = {"BOND": 2, "ANGL": 3, "DIHE": 4, "IMPR": 4}[section]

    try:
        types = _split_type_spec(type_spec, n_expected_types)
    except ValueError:
        return None

    try:
        values = [float(t) for t in tokens[val_start:val_start + n_trailing]]
    except (ValueError, IndexError):
        return None

    return (types, values)


def parse_frcmod(path: Union[str, Path]) -> Dict:
    """Parse an AMBER frcmod into raw section dicts.

    Returns a dict with keys:
        mass:      {type_name: mass}
        nonb:      {type_name: (R*, epsilon)}
        bonds:     [(types, [K, r0]), ...]
        angles:    [(types, [K, theta0]), ...]
        dihedrals: [(types, [idivf, Vn, phase, period]), ...]
        impropers: [(types, [Vn, phase, period]), ...]

    Bonded entries are returned as ordered lists (canonicalization is
    applied later in ``parse_set``).
    """
    out: Dict = {
        "mass": {},
        "nonb": {},
        "bonds": [],
        "angles": [],
        "dihedrals": [],
        "impropers": [],
    }

    section = None
    with open(path) as f:
        for line in f:
            raw = line.rstrip("\n")
            head = raw.lstrip().split(maxsplit=1)[0].upper() if raw.strip() else ""
            if head in _FRCMOD_SECTIONS:
                section = head
                if section == "ANGLE":
                    section = "ANGL"
                if section == "NONBON":
                    section = "NONB"
                continue
            if section is None:
                continue
            if not raw.strip():
                # Blank line; AMBER convention says a blank line ends the section.
                # But many frcmod files have stray blanks inside sections too,
                # so we DON'T reset section here. Let the next header reset it.
                continue

            parsed = _parse_frcmod_line(raw, section)
            if parsed is None:
                continue
            key, values = parsed
            if section == "MASS":
                out["mass"][key] = values[0]
            elif section == "NONB":
                out["nonb"][key] = (values[0], values[1])
            elif section == "BOND":
                out["bonds"].append((key, values))
            elif section == "ANGL":
                out["angles"].append((key, values))
            elif section == "DIHE":
                out["dihedrals"].append((key, values))
            elif section == "IMPR":
                out["impropers"].append((key, values))

    return out


# ---------------------------------------------------------------------------
# lib (OFF) parsing
# ---------------------------------------------------------------------------

_LIB_ATOMS_HEADER_RE = re.compile(r"^!entry\.(\S+)\.unit\.atoms\b")
_LIB_CONNECTIVITY_HEADER_RE = re.compile(r"^!entry\.(\S+)\.unit\.connectivity\b")
# `unit.connect` lists 1-indexed atom indices that bond to adjacent residues.
# Distinct from `unit.connectivity` (the intra-unit bond table).
_LIB_CONNECT_HEADER_RE = re.compile(r"^!entry\.(\S+)\.unit\.connect\s+array")

# Unit names to skip when parsing lib files. These are leftover artifacts from
# upstream parameterization tooling (e.g., antechamber's default ``mol``
# placeholder) that ship in some bundled libs without being referenced by any
# residue template. Including them would generate spurious cross-cofactor
# conflicts on the placeholder name.
_LIB_IGNORE_UNITS = {"mol"}


@dataclass
class LibUnitInfo:
    """Per-unit lib content used by the compatibility tools."""
    # (unit_name, atom_name) -> atom_type
    unit_types: Dict[Tuple[str, str], str] = field(default_factory=dict)
    # (unit_name, atom_name) -> atomic number (column 7 of the atoms table)
    unit_elements: Dict[Tuple[str, str], int] = field(default_factory=dict)
    # (unit_name, atom_name) -> bond degree (number of times the atom appears
    # in the unit.connectivity table, summed over both endpoints).
    unit_degrees: Dict[Tuple[str, str], int] = field(default_factory=dict)


def parse_lib(path: Union[str, Path]) -> LibUnitInfo:
    """Parse an OFF library file's atom and connectivity tables.

    Captures, per (unit, atom_name):
        atom_type   — column 2 of the ``atoms`` table
        element     — column 7 of the ``atoms`` table (atomic number)
        degree      — bond count from the ``connectivity`` table

    Skips charges and other columns we don't need.
    """
    out = LibUnitInfo()
    current_unit: Optional[str] = None
    in_atoms_table = False
    in_connectivity_table = False
    in_connect_array = False
    # Per-unit ordered atom names — needed to map connectivity's 1-indexed
    # atom indices back to atom_name. Atom-table rows are written in seq order,
    # which is the index space connectivity references.
    unit_atom_order: Dict[str, List[str]] = {}

    def _close_active_section() -> None:
        nonlocal in_atoms_table, in_connectivity_table, in_connect_array
        in_atoms_table = False
        in_connectivity_table = False
        in_connect_array = False

    with open(path) as f:
        for line in f:
            raw = line.rstrip("\n")
            stripped = raw.strip()
            if not stripped:
                continue
            m_atoms = _LIB_ATOMS_HEADER_RE.match(stripped)
            m_conn = _LIB_CONNECTIVITY_HEADER_RE.match(stripped)
            m_connpt = _LIB_CONNECT_HEADER_RE.match(stripped)
            if m_atoms:
                unit = m_atoms.group(1)
                if unit in _LIB_IGNORE_UNITS:
                    current_unit = None
                    _close_active_section()
                    continue
                current_unit = unit
                unit_atom_order.setdefault(unit, [])
                _close_active_section()
                in_atoms_table = True
                continue
            if m_conn:
                unit = m_conn.group(1)
                if unit in _LIB_IGNORE_UNITS:
                    current_unit = None
                    _close_active_section()
                    continue
                current_unit = unit
                _close_active_section()
                in_connectivity_table = True
                continue
            if m_connpt:
                unit = m_connpt.group(1)
                if unit in _LIB_IGNORE_UNITS:
                    current_unit = None
                    _close_active_section()
                    continue
                current_unit = unit
                _close_active_section()
                in_connect_array = True
                continue
            if stripped.startswith("!"):
                # any other section header ends whichever table we're inside
                _close_active_section()
                current_unit = None
                continue
            if current_unit is None:
                continue

            if in_atoms_table:
                # Atom table row format:
                #   "ATOM_NAME" "ATOM_TYPE" 0 1 131072 1 6 -0.035100
                tokens = stripped.split()
                if len(tokens) < 2:
                    continue
                atom_name = tokens[0].strip('"')
                atom_type = tokens[1].strip('"')
                out.unit_types[(current_unit, atom_name)] = atom_type
                unit_atom_order[current_unit].append(atom_name)
                # Column 7 (token index 6) is the atomic number. Older / unusual
                # libs may omit it — fall back to silently skipping the element.
                if len(tokens) >= 7:
                    try:
                        out.unit_elements[(current_unit, atom_name)] = int(tokens[6])
                    except ValueError:
                        pass
                continue

            if in_connectivity_table:
                # Connectivity row: " atom1x atom2x flags"  (1-indexed)
                tokens = stripped.split()
                if len(tokens) < 2:
                    continue
                try:
                    a1 = int(tokens[0])
                    a2 = int(tokens[1])
                except ValueError:
                    continue
                order = unit_atom_order.get(current_unit, [])
                if not (1 <= a1 <= len(order) and 1 <= a2 <= len(order)):
                    continue
                name1 = order[a1 - 1]
                name2 = order[a2 - 1]
                k1 = (current_unit, name1)
                k2 = (current_unit, name2)
                out.unit_degrees[k1] = out.unit_degrees.get(k1, 0) + 1
                out.unit_degrees[k2] = out.unit_degrees.get(k2, 0) + 1
                continue

            if in_connect_array:
                # unit.connect row: a single 1-indexed atom number per line.
                # These atoms bond to atoms in adjacent residues (not shown in
                # this unit's connectivity table), so we add +1 to their degree.
                tokens = stripped.split()
                if not tokens:
                    continue
                try:
                    idx = int(tokens[0])
                except ValueError:
                    continue
                order = unit_atom_order.get(current_unit, [])
                if not (1 <= idx <= len(order)):
                    continue
                name = order[idx - 1]
                key = (current_unit, name)
                out.unit_degrees[key] = out.unit_degrees.get(key, 0) + 1
                continue

    # Make sure every atom that appeared in the atoms table has a degree
    # entry (zero if it has no bonds listed — rare but possible for lone metals).
    for key in out.unit_types:
        out.unit_degrees.setdefault(key, 0)

    return out


# ---------------------------------------------------------------------------
# Whole-set assembly
# ---------------------------------------------------------------------------

def parse_set(
    set_id: str,
    frcmod_path: Union[str, Path, List[Union[str, Path]]],
    lib_paths: List[Union[str, Path]],
) -> FFSignature:
    """Parse a complete FF set into a canonical FFSignature.

    Args:
        set_id: Stable identifier (typically "<cofactor_path>:<set_name>",
            e.g., "heme/bis_his_c_type:Henriques_Reduced").
        frcmod_path: Path to the set's frcmod file, OR a list of them. A metal
            site carries more than one (the MCPB bonded frcmod plus each organic
            ligand's own GAFF frcmod); all are parsed and merged into one
            signature, later entries winning on duplicate keys (tleap order).
        lib_paths: Paths to the set's OFF library file(s).
    """
    frcmod_paths = (frcmod_path if isinstance(frcmod_path, (list, tuple))
                    else [frcmod_path])
    sig = FFSignature(
        set_id=set_id,
        source_files=[*map(str, frcmod_paths), *map(str, lib_paths)],
    )

    for one_frcmod in frcmod_paths:
        raw = parse_frcmod(one_frcmod)

        sig.mass.update(raw["mass"])
        sig.nonb.update(raw["nonb"])

        # Canonicalize bonded entries before storing
        for types, values in raw["bonds"]:
            key = _canon_bond(*types)
            # If the same canonical key appears twice with different values, the
            # later one wins by tleap convention. Mirror that here.
            sig.bonds[key] = (values[0], values[1])

        for types, values in raw["angles"]:
            key = _canon_angle(*types)
            sig.angles[key] = (values[0], values[1])

        for types, values in raw["dihedrals"]:
            # values: [idivf, Vn, phase, period]
            # AMBER allows negative periodicity (continuation marker); take abs for keying.
            idivf, vn, phase, period = values
            period_key = int(round(abs(period)))
            canon = _canon_dihedral(*types)
            sig.dihedrals[(canon, period_key)] = (idivf, vn, phase)

        for types, values in raw["impropers"]:
            vn, phase, period = values
            period_key = int(round(abs(period)))
            canon = _canon_improper(*types)
            sig.impropers[(canon, period_key)] = (vn, phase)

    for lib_path in lib_paths:
        info = parse_lib(lib_path)
        sig.library_unit_types.update(info.unit_types)
        sig.library_unit_elements.update(info.unit_elements)
        sig.library_unit_degrees.update(info.unit_degrees)

    return sig
