"""
Named, replayable edits to a molecule's chemistry.

Every change to bond orders, formal charges or hydrogens -- whether it
resolves a problem in a CCD entry or sets a protonation state the user chose
in the viewer -- is a ChemistryEdit: an operation plus atom NAMES, never atom
indices, so the same list applied to a rebuilt molecule gives the same result.
The list is what a session records and replays.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

from rdkit import Chem

from .ccd_chemistry import ChemistryProblem, find_valence_problems

_ORDERS = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}

OPERATIONS = ("set_bond_order", "set_formal_charge", "add_hydrogen", "remove_atom")


@dataclass(frozen=True)
class ChemistryEdit:
    """One change. ``atoms`` holds one name (charge, add H, remove) or two (bond order)."""

    op: str
    atoms: Tuple[str, ...]
    value: Optional[int] = None

    def __post_init__(self):
        if self.op not in OPERATIONS:
            raise ValueError(f"Unknown chemistry edit {self.op!r}; expected one of {', '.join(OPERATIONS)}.")
        expected = 2 if self.op == "set_bond_order" else 1
        if len(self.atoms) != expected:
            raise ValueError(f"{self.op} takes {expected} atom name(s), got {list(self.atoms)}.")
        if self.op in ("set_bond_order", "set_formal_charge") and self.value is None:
            raise ValueError(f"{self.op} needs a value.")
        if self.op == "set_bond_order" and self.value not in _ORDERS:
            raise ValueError(f"Bond order must be 1, 2 or 3 (Kekule form), got {self.value}.")

    def describe(self) -> str:
        if self.op == "set_bond_order":
            return f"bond {self.atoms[0]}-{self.atoms[1]} -> order {self.value}"
        if self.op == "set_formal_charge":
            return f"{self.atoms[0]} formal charge -> {self.value:+d}"
        if self.op == "add_hydrogen":
            return f"add H on {self.atoms[0]}"
        return f"remove {self.atoms[0]}"

    def to_dict(self) -> dict:
        data = {"op": self.op, "atoms": list(self.atoms)}
        if self.value is not None:
            data["value"] = self.value
        return data

    @classmethod
    def from_dict(cls, data: dict) -> "ChemistryEdit":
        return cls(op=data["op"], atoms=tuple(data["atoms"]), value=data.get("value"))


def atom_name(atom: Chem.Atom) -> str:
    return atom.GetProp("name")


def _index_by_name(mol: Chem.Mol) -> Dict[str, int]:
    index: Dict[str, int] = {}
    for atom in mol.GetAtoms():
        if not atom.HasProp("name"):
            raise ValueError(f"Atom {atom.GetIdx() + 1} ({atom.GetSymbol()}) has no name; "
                             "edits address atoms by name.")
        name = atom_name(atom)
        if name in index:
            raise ValueError(f"Atom name {name!r} is used twice; edits need unique names.")
        index[name] = atom.GetIdx()
    return index


def new_hydrogen_name(heavy_name: str, taken: Iterable[str]) -> str:
    """A PDB-width (<= 4 character) hydrogen name not already in use, derived from its heavy atom."""
    taken = set(taken)
    stem = heavy_name[1:] if len(heavy_name) > 1 else heavy_name
    candidates = [f"H{stem}"] + [f"H{stem}{n}" for n in range(1, 10)] + [f"H{n}{stem}" for n in range(1, 10)]
    candidates += [f"H{n}" for n in range(1, 1000)]
    for candidate in candidates:
        if len(candidate) <= 4 and candidate not in taken:
            return candidate
    raise ValueError(f"No free hydrogen name for {heavy_name}.")


def apply_edits(mol: Chem.Mol, edits: Iterable[ChemistryEdit]) -> Tuple[Chem.Mol, List[ChemistryProblem]]:
    """Apply ``edits`` in order to a copy of ``mol``; return it and the problems that remain.

    The molecule is sanitized only when no problems remain, so a partial set of
    edits can be applied, inspected and continued.
    """
    rw = Chem.RWMol(mol)
    for edit in edits:
        names = _index_by_name(rw)
        missing = [name for name in edit.atoms if name not in names]
        if missing:
            raise KeyError(f"{edit.describe()}: no atom named {', '.join(missing)}.")
        if edit.op == "set_bond_order":
            i, j = (names[name] for name in edit.atoms)
            bond = rw.GetBondBetweenAtoms(i, j)
            if bond is None:
                raise ValueError(f"{edit.describe()}: {edit.atoms[0]} and {edit.atoms[1]} are not bonded.")
            bond.SetIsAromatic(False)
            bond.SetBondType(_ORDERS[edit.value])
        elif edit.op == "set_formal_charge":
            rw.GetAtomWithIdx(names[edit.atoms[0]]).SetFormalCharge(int(edit.value))
        elif edit.op == "remove_atom":
            rw.RemoveAtom(names[edit.atoms[0]])
        elif edit.op == "add_hydrogen":
            _add_hydrogen(rw, names[edit.atoms[0]], names)
    for atom in rw.GetAtoms():
        atom.SetNoImplicit(True)
    result = rw.GetMol()
    result.UpdatePropertyCache(strict=False)
    problems = find_valence_problems(result)
    if not problems:
        Chem.SanitizeMol(result)
    return result, problems


def _add_hydrogen(rw: Chem.RWMol, heavy_index: int, names: Dict[str, int]) -> None:
    heavy = rw.GetAtomWithIdx(heavy_index)
    hydrogen = Chem.Atom(1)
    hydrogen.SetNoImplicit(True)
    hydrogen.SetProp("name", new_hydrogen_name(atom_name(heavy), names))
    h_index = rw.AddAtom(hydrogen)
    rw.AddBond(heavy_index, h_index, Chem.BondType.SINGLE)
    if rw.GetNumConformers():
        _place_hydrogen(rw, heavy_index, h_index)


def _place_hydrogen(rw: Chem.RWMol, heavy_index: int, h_index: int) -> None:
    """Position a new H with RDKit's own placement, computed on a scratch copy."""
    scratch = Chem.RWMol(rw)
    scratch.RemoveAtom(h_index)
    scratch.GetAtomWithIdx(heavy_index).SetNumExplicitHs(1)
    scratch.GetAtomWithIdx(heavy_index).SetNoImplicit(True)
    scratch.UpdatePropertyCache(strict=False)
    placed = Chem.AddHs(scratch, addCoords=True, onlyOnAtoms=(heavy_index,), explicitOnly=True)
    new_atom = placed.GetNumAtoms() - 1
    position = placed.GetConformer().GetAtomPosition(new_atom)
    for conformer in rw.GetConformers():
        conformer.SetAtomPosition(h_index, position)
