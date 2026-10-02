"""
Ligand PDBQT with Meeko, and exact control over which bonds rotate.

Meeko starts with every single bond rotatable, then removes amides (unless
``flexible_amides``), bonds next to a triple bond, ring bonds and bonds whose
rotation moves nothing. Its own override, ``rigidify_bonds_smarts`` /
``rigidify_bonds_indices``, selects bonds by SMARTS, which cannot tell
symmetry-equivalent bonds apart (one ester tail of a symmetric diester from
the other). The user chooses individual bonds, so the choice is applied by
wrapping Meeko's bond typer: it runs after Meeko's rules and before the torsion
tree is built, and sets the chosen atom pairs rigid or rotatable exactly.

Only two kinds of change are offered, because only these leave Meeko's torsion
tree valid: a rotatable bond may be made rigid, and an amide C-N bond may be
made rotatable. Ring bonds stay rigid (Vina does not sample ring conformations).

Charges are Gasteiger, the model AutoDock4's scoring function was calibrated
with; Vina ignores charges.
"""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, Iterable, List, Optional, Sequence, Set, Tuple

from rdkit import Chem

BondNames = Tuple[str, str]


def _key(names: Sequence[str]) -> FrozenSet[str]:
    if len(names) != 2 or names[0] == names[1]:
        raise ValueError(f"A bond is two different atom names, got {list(names)}.")
    return frozenset(names)


@dataclass
class TorsionBond:
    """A bond the user can choose to rotate or not."""

    atoms: BondNames
    kind: str                  # "single", "conjugated single" or "amide"
    rotatable_by_default: bool

    def label(self) -> str:
        return f"{self.atoms[0]}-{self.atoms[1]}"


_GLUE_TYPE = re.compile(r"^(CG|G)\d+$")


def pdbqt_atom_type(line: str) -> str:
    """The AutoDock type of an ATOM/HETATM line. It starts at column 78 and may run past
    column 79: Meeko's macrocycle glue types (CG0, G0) are three characters."""
    return line[77:].strip()


@dataclass
class PreparedLigand:
    pdbqt: str
    rotatable: List[BondNames]     # bonds that rotate during docking (incl. an opened macrocycle's)
    torsdof: int                   # the file's TORSDOF: torsions of the unopened molecule (AD4 entropy term)
    atom_types: Dict[str, int]
    gasteiger_charge_sum: float
    notes: List[str] = field(default_factory=list)


class _ExactBondChoices:
    """Meeko bond typer wrapper: Meeko's rules first, then the user's exact choices."""

    def __init__(self, inner, rigid: Set[Tuple[int, int]], rotatable: Set[Tuple[int, int]]):
        self.inner = inner
        self.rigid = rigid
        self.rotatable = rotatable

    def __call__(self, setup, flexible_amides, rigidify_bonds_smarts, rigidify_bonds_indices):
        self.inner(setup, flexible_amides, rigidify_bonds_smarts, rigidify_bonds_indices)
        for bond in self.rigid:
            setup.bond_info[bond].rotatable = False
        for bond in self.rotatable:
            setup.bond_info[bond].rotatable = True


def _names_to_index(mol: Chem.Mol) -> Dict[str, int]:
    return {atom.GetProp("name"): atom.GetIdx() for atom in mol.GetAtoms()}


def _bond_id(i: int, j: int) -> Tuple[int, int]:
    return (i, j) if i < j else (j, i)


_CHARGE_PROP = "proprep_partial_charge"


def _with_charges(mol: Chem.Mol, metal_charges: Optional[Dict[str, int]]) -> Tuple[Chem.Mol, Optional[str], List[str]]:
    """A working copy of ``mol`` and, when it holds metals, charges for Meeko to read.

    Without metals Meeko computes Gasteiger charges itself. With metals its own
    route removes each metal and protonates the atoms bound to it, which fails
    for chemistry drawn with dative bonds and charged donors (a heme's
    pyrrolide nitrogens), and Meeko 0.8.0 then crashes on an unbound variable.
    Here the metal-free molecule, with its own formal charges, gets Gasteiger
    charges, and each metal carries the formal charge the user chose, as in the
    receptor. The total equals the sum of formal charges.
    """
    from .ccd_chemistry import is_metal
    metals = [a for a in mol.GetAtoms() if is_metal(a)]
    if not metals:
        if metal_charges:
            raise ValueError("Metal charges given for a ligand without metals.")
        return mol, None, []
    names = [a.GetProp("name") for a in metals]
    missing = [n for n in names if n not in (metal_charges or {})]
    if missing:
        raise ValueError(f"Choose a formal charge for every metal in the ligand: {', '.join(missing)}.")
    work = Chem.Mol(mol)
    for atom in work.GetAtoms():
        if is_metal(atom):
            atom.SetFormalCharge(int(metal_charges[atom.GetProp("name")]))
    work.UpdatePropertyCache(strict=False)
    stripped = Chem.RWMol(work)
    for atom in sorted((a for a in work.GetAtoms() if is_metal(a)), key=lambda a: a.GetIdx(), reverse=True):
        stripped.RemoveAtom(atom.GetIdx())
    stripped = stripped.GetMol()
    Chem.SanitizeMol(stripped)
    from rdkit.Chem import rdPartialCharges
    rdPartialCharges.ComputeGasteigerCharges(stripped)
    organic = iter(float(a.GetDoubleProp("_GasteigerCharge")) for a in stripped.GetAtoms())
    for atom in work.GetAtoms():
        charge = float(atom.GetFormalCharge()) if is_metal(atom) else next(organic)
        if charge != charge:                            # NaN: Gasteiger has no parameters for this atom
            raise ValueError(f"Gasteiger has no parameters for {atom.GetProp('name')} ({atom.GetSymbol()}).")
        atom.SetProp(_CHARGE_PROP, repr(charge))
    notes = [f"Metal(s) {', '.join(f'{n} {int(metal_charges[n]):+d}' for n in names)} carry their chosen formal "
             "charge; Gasteiger charges for the rest were computed without the metal (the receptor convention, "
             "not Meeko's surrogate-proton scheme)."]
    return work, _CHARGE_PROP, notes


def _preparation(flexible_amides: bool, rigid=frozenset(), rotatable=frozenset(), rigid_macrocycles=False,
                 charge_prop: Optional[str] = None):
    from meeko import MoleculePreparation
    charges = dict(charge_model="read", charge_atom_prop=charge_prop) if charge_prop else dict(charge_model="gasteiger")
    mk = MoleculePreparation(flexible_amides=flexible_amides, rigid_macrocycles=rigid_macrocycles, **charges)
    if not hasattr(mk, "_bond_typer"):
        raise RuntimeError("This Meeko version has no _bond_typer; per-bond torsion choices "
                           "were validated with Meeko 0.8.0.")
    mk._bond_typer = _ExactBondChoices(mk._bond_typer, set(rigid), set(rotatable))
    return mk


def _rotatable_ids(mol: Chem.Mol, flexible_amides: bool, rigid_macrocycles: bool,
                   charge_prop: Optional[str] = None) -> Set[Tuple[int, int]]:
    setup = _preparation(flexible_amides, rigid_macrocycles=rigid_macrocycles, charge_prop=charge_prop).prepare(mol)[0]
    n_atoms = mol.GetNumAtoms()
    return {bond_id for bond_id, info in setup.bond_info.items()
            if info.rotatable and max(bond_id) < n_atoms}         # skip pseudo-atom bonds


def torsion_bonds(mol: Chem.Mol, *, rigid_macrocycles: bool,
                  metal_charges: Optional[Dict[str, int]] = None) -> List[TorsionBond]:
    """Every bond the user may switch: Meeko's default rotatable bonds, plus amide C-N bonds.

    Listed in atom order. ``conjugated single`` marks single bonds RDKit
    perceives as conjugated (an ester C(=O)-O, the single bonds of a polyene),
    which the user may want to hold rigid. ``macrocycle ring`` marks ring bonds
    that rotate because Meeko opens the macrocycle; they exist only with
    ``rigid_macrocycles`` False, the same setting prepare_ligand must be given.
    """
    mol, charge_prop, _ = _with_charges(mol, metal_charges)
    default = _rotatable_ids(mol, False, rigid_macrocycles, charge_prop)
    with_amides = _rotatable_ids(mol, True, rigid_macrocycles, charge_prop)
    bonds = []
    for bond_id in sorted(default | with_amides):
        i, j = bond_id
        rd_bond = mol.GetBondBetweenAtoms(i, j)
        if bond_id in with_amides and bond_id not in default:
            kind = "amide"
        elif rd_bond.IsInRing():
            kind = "macrocycle ring"
        elif rd_bond.GetIsConjugated():
            kind = "conjugated single"
        else:
            kind = "single"
        names = (mol.GetAtomWithIdx(i).GetProp("name"), mol.GetAtomWithIdx(j).GetProp("name"))
        bonds.append(TorsionBond(names, kind, bond_id in default))
    return bonds


def prepare_ligand(mol: Chem.Mol, rigid: Iterable[BondNames] = (),
                   rotatable_amides: Iterable[BondNames] = (), *,
                   rigid_macrocycles: bool, metal_charges: Optional[Dict[str, int]] = None) -> PreparedLigand:
    """PDBQT for ``mol`` with Meeko's default torsions, minus ``rigid``, plus ``rotatable_amides``.

    Every named bond must be one ``torsion_bonds`` offers: ``rigid`` a bond
    rotatable by default, ``rotatable_amides`` an amide bond. Anything else is
    refused with the bond named.

    ``rigid_macrocycles`` False lets Meeko open a macrocycle at one bond, keeping
    the ends together with glue pseudo-atoms (types CG/G) so the ring flexes
    during docking; True keeps every ring as built. Which happened is noted.
    """
    from meeko import PDBQTWriterLegacy
    index = _names_to_index(mol)
    offered = {_key(b.atoms): b for b in torsion_bonds(mol, rigid_macrocycles=rigid_macrocycles,
                                                       metal_charges=metal_charges)}
    mol, charge_prop, charge_notes = _with_charges(mol, metal_charges)

    def resolve(pairs, wanted_kind):
        ids = set()
        for pair in pairs:
            key = _key(pair)
            missing = [name for name in pair if name not in index]
            if missing:
                raise KeyError(f"No atom named {', '.join(missing)}.")
            bond = offered.get(key)
            if bond is None:
                raise ValueError(f"{pair[0]}-{pair[1]} is not a bond that can rotate "
                                 "(not bonded, in a ring, double/triple, or moves nothing).")
            if wanted_kind == "rigid" and not bond.rotatable_by_default:
                raise ValueError(f"{bond.label()} does not rotate by default ({bond.kind}); nothing to rigidify.")
            if wanted_kind == "amide" and bond.kind != "amide":
                raise ValueError(f"{bond.label()} is not an amide bond; it rotates by default.")
            ids.add(_bond_id(index[pair[0]], index[pair[1]]))
        return ids

    rigid_ids = resolve(list(rigid), "rigid")
    amide_ids = resolve(list(rotatable_amides), "amide")
    overlap = rigid_ids & amide_ids
    if overlap:
        raise ValueError("A bond cannot be both rigid and rotatable.")

    setups = _preparation(False, rigid_ids, amide_ids, rigid_macrocycles, charge_prop).prepare(mol)
    if len(setups) != 1:
        raise RuntimeError(f"Meeko produced {len(setups)} setups for one ligand; expected 1.")
    setup = setups[0]
    pdbqt, ok, error = PDBQTWriterLegacy.write_string(setup)
    if not ok:
        raise RuntimeError(f"Meeko could not write the ligand PDBQT: {error}")

    rotatable = sorted(
        (mol.GetAtomWithIdx(i).GetProp("name"), mol.GetAtomWithIdx(j).GetProp("name"))
        for (i, j), info in setup.bond_info.items() if info.rotatable
    )
    records = [line for line in pdbqt.splitlines() if line.startswith(("ATOM", "HETATM"))]
    torsdof = [int(line.split()[1]) for line in pdbqt.splitlines() if line.startswith("TORSDOF")]
    notes = [f"{len(rotatable)} rotatable bonds"
             + (f" ({len(rigid_ids)} held rigid by choice)" if rigid_ids else "")
             + (f", {len(amide_ids)} amide bond(s) made rotatable" if amide_ids else "") + ".",
             "Partial charges: Gasteiger (Meeko); non-polar hydrogens merged into their carbons."]
    notes.extend(charge_notes)
    glue = [line for line in records if _GLUE_TYPE.match(pdbqt_atom_type(line))]
    if glue:
        notes.append(f"A macrocycle was opened for flexible docking ({len(glue)} glue atoms, types CG/G).")
    elif rigid_macrocycles:
        notes.append("Rings kept rigid as built (macrocycle flexibility off).")
    return PreparedLigand(
        pdbqt=pdbqt,
        rotatable=rotatable,
        torsdof=torsdof[0] if torsdof else 0,
        atom_types=dict(Counter(pdbqt_atom_type(line) for line in records)),
        gasteiger_charge_sum=round(sum(float(line[70:76]) for line in records), 3),
        notes=notes,
    )
