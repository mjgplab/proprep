"""
Chemical Component Dictionary (CCD) entries as RDKit molecules, checked before use.

The CCD is the source of bond orders, formal charges and hydrogens for any
residue in a deposited structure, but an entry can be internally inconsistent.
Every heme entry (HEM, HEC, HEA, HEB) currently gives the propionate oxygens
O2A and O2D a charge of -1 while also bonding a hydrogen to each; RDKit rejects
the valence, and Meeko's own CCD route then fails. Such problems are reported
atom by atom for the user to resolve; nothing here changes the chemistry.

The CCD text comes from ProPrep's existing cached copy of the full dictionary
(``CCDParser.get_component_block``), so no second download path exists.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import gemmi
from rdkit import Chem, RDLogger
from rdkit.Geometry import Point3D

from proprep.structure_prep.pdb_searcher import NONMETAL_ELEMENTS

_BOND_ORDERS = {
    "SING": Chem.BondType.SINGLE,
    "DOUB": Chem.BondType.DOUBLE,
    "TRIP": Chem.BondType.TRIPLE,
}


@dataclass
class ChemistryProblem:
    """One atom whose element, formal charge and bonds are not a valid valence."""

    atom_name: str
    element: str
    formal_charge: int
    bonded_atoms: List[str]
    detail: str

    def describe(self) -> str:
        charge = f"{self.formal_charge:+d}" if self.formal_charge else "0"
        return (f"{self.atom_name} ({self.element}, charge {charge}) bonded to "
                f"{', '.join(self.bonded_atoms)}: {self.detail}")


@dataclass
class CcdComponent:
    """A CCD entry: its molecule with CCD atom names, and any problems found."""

    code: str
    name: str
    mol: Chem.Mol                      # all CCD atoms incl. H; not sanitized if problems exist
    problems: List[ChemistryProblem] = field(default_factory=list)
    coordinates: str = "ideal"         # "ideal", "model", or "none"

    def atom_index(self) -> Dict[str, int]:
        return {atom.GetProp("name"): atom.GetIdx() for atom in self.mol.GetAtoms()}


def _float_or_none(value: str) -> Optional[float]:
    if value in ("?", "."):
        return None
    return float(value)


def _element(type_symbol: str) -> str:
    symbol = type_symbol.strip()
    return symbol[0].upper() + symbol[1:].lower()


def component_from_block(block_text: str, code: str) -> CcdComponent:
    """Build a CcdComponent from one mmCIF data block (text).

    Bonds use the CCD's Kekule ``value_order``; hydrogens are the CCD's own
    explicit atoms, so no implicit hydrogens are added. Coordinates come from
    the ideal set when every atom has one, otherwise the model set.
    """
    block = gemmi.cif.read_string(block_text).sole_block()
    code = code.upper()
    atoms = block.find(
        "_chem_comp_atom.",
        ["atom_id", "type_symbol", "charge",
         "?pdbx_model_Cartn_x_ideal", "?pdbx_model_Cartn_y_ideal", "?pdbx_model_Cartn_z_ideal",
         "?model_Cartn_x", "?model_Cartn_y", "?model_Cartn_z"],
    )
    rw = Chem.RWMol()
    names: Dict[str, int] = {}
    ideal, model = [], []
    for row in atoms:
        name = gemmi.cif.as_string(row[0])
        atom = Chem.Atom(_element(row[1]))
        charge = row[2]
        atom.SetFormalCharge(0 if charge in ("?", ".") else int(charge))
        atom.SetNoImplicit(True)
        atom.SetProp("name", name)
        names[name] = rw.AddAtom(atom)
        ideal.append(tuple(_float_or_none(row[i]) if row.has(i) else None for i in (3, 4, 5)))
        model.append(tuple(_float_or_none(row[i]) if row.has(i) else None for i in (6, 7, 8)))

    for row in block.find("_chem_comp_bond.", ["atom_id_1", "atom_id_2", "value_order"]):
        order = _BOND_ORDERS.get(row[2].upper())
        if order is None:
            raise ValueError(f"CCD {code}: unsupported bond order {row[2]!r} "
                             f"between {row[0]} and {row[1]}")
        rw.AddBond(names[gemmi.cif.as_string(row[0])], names[gemmi.cif.as_string(row[1])], order)

    coordinates = "none"
    for label, coords in (("ideal", ideal), ("model", model)):
        if coords and all(None not in xyz for xyz in coords):
            conformer = Chem.Conformer(rw.GetNumAtoms())
            for index, xyz in enumerate(coords):
                conformer.SetAtomPosition(index, Point3D(*xyz))
            rw.AddConformer(conformer, assignId=True)
            coordinates = label
            break

    make_metal_bonds_dative(rw)
    mol = rw.GetMol()
    mol.UpdatePropertyCache(strict=False)
    problems = find_valence_problems(mol)
    if not problems:
        Chem.SanitizeMol(mol)
    name = block.find_value("_chem_comp.name")
    return CcdComponent(code=code, name=gemmi.cif.as_string(name) if name else "",
                        mol=mol, problems=problems, coordinates=coordinates)


def is_metal_symbol(symbol: str) -> bool:
    """An element symbol outside ProPrep's closed list of nonmetals (which includes D and T)."""
    return symbol.strip() not in NONMETAL_ELEMENTS


def is_metal(atom: Chem.Atom) -> bool:
    return is_metal_symbol(atom.GetSymbol())


def make_metal_bonds_dative(rw: Chem.RWMol) -> int:
    """Turn single bonds between a nonmetal and a metal into dative bonds pointing at the metal.

    The CCD draws coordination (Fe-N of a porphyrin, Mg-N of chlorophyll) as
    ordinary single bonds, which RDKit counts toward the donor's valence: a
    pyrrolide N- with three bonds, or a neutral N with four, then reads as
    over-valent. A dative bond does not count toward the donor. Bonds between
    two metals, or between two nonmetals, are left alone. Returns the number
    of bonds changed.
    """
    changed = 0
    for bond in list(rw.GetBonds()):
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        if bond.GetBondType() != Chem.BondType.SINGLE or is_metal(a) == is_metal(b):
            continue
        donor, metal = (b, a) if is_metal(a) else (a, b)
        rw.RemoveBond(a.GetIdx(), b.GetIdx())
        rw.AddBond(donor.GetIdx(), metal.GetIdx(), Chem.BondType.DATIVE)
        changed += 1
    return changed


def find_valence_problems(mol: Chem.Mol) -> List[ChemistryProblem]:
    """Atoms whose element, formal charge and bonds are not an allowed valence."""
    problems: List[ChemistryProblem] = []
    RDLogger.DisableLog("rdApp.error")      # the structured list below replaces RDKit's log lines
    try:
        issues = Chem.DetectChemistryProblems(mol)
    finally:
        RDLogger.EnableLog("rdApp.error")
    for issue in issues:
        if not hasattr(issue, "GetAtomIdx"):
            problems.append(ChemistryProblem("?", "?", 0, [], issue.Message()))
            continue
        atom = mol.GetAtomWithIdx(issue.GetAtomIdx())
        neighbours = [n.GetProp("name") if n.HasProp("name") else f"{n.GetSymbol()}{n.GetIdx() + 1}"
                      for n in atom.GetNeighbors()]
        problems.append(ChemistryProblem(
            atom_name=atom.GetProp("name") if atom.HasProp("name") else f"{atom.GetSymbol()}{atom.GetIdx() + 1}",
            element=atom.GetSymbol(),
            formal_charge=atom.GetFormalCharge(),
            bonded_atoms=neighbours,
            detail=issue.Message(),
        ))
    return problems


def fetch_component(code: str, parser=None) -> CcdComponent:
    """The CCD entry for ``code`` from ProPrep's cached dictionary.

    Raises LookupError when the CCD has no such code; network errors from the
    first download propagate unchanged.
    """
    if parser is None:
        from proprep.structure_prep.chem_comp_dict_fetcher import CCDParser
        parser = CCDParser()
    block = parser.get_component_block(code)
    if block is None:
        raise LookupError(f"The Chemical Component Dictionary has no entry {code.upper()!r}.")
    return component_from_block(block, code)
