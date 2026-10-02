"""
Meeko residue templates for what its built-in set does not cover, or covers
with a choice the user must make.

Meeko matches every receptor residue to a template (a SMILES with atom names
and "link labels" on atoms bonded to neighbouring residues) and rebuilds the
hydrogens from it. Its built-in set covers standard residues and a few ions,
but a cofactor it does not know is downloaded from the CCD on every run, and
its ion templates fix oxidation states nobody chose (FE is [Fe+3]). This
module builds those templates from inputs ProPrep has checked:

* metal ions with the charge the user chose. The residue name is the element
  plus the charge (FE2, FE3, CU1, ZN2), so two oxidation states of one
  element can coexist in a receptor.
* cofactors from their checked CCD chemistry (``ccd_chemistry`` +
  ``chemistry_edits``), with metal atoms split out as ions of their own,
  because Meeko does not match a template in which the metal sits as a
  disconnected fragment.
* covalent links: a residue bonded to another (Cys SG to a heme c vinyl
  carbon) needs a template variant with one hydrogen fewer on the link atom
  and a link label there, on both sides of the bond.

A link needs a "padder", a reaction that caps the residue with the neighbour's
atoms when Meeko types it on its own. Meeko applies the padder's reactant
SMARTS to the residue and requires EXACTLY ONE match. Two chemically similar
link atoms in one residue (heme c's CAB and CAC) therefore need patterns large
enough to tell them apart; ``unique_environment_smarts`` grows the atom
environment around the link atom until the pattern matches once.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from rdkit import Chem

from .ccd_chemistry import CcdComponent, is_metal


@dataclass
class TemplateSet:
    """Templates and padders to add to Meeko's defaults."""

    residue_templates: Dict[str, dict] = field(default_factory=dict)
    padders: Dict[str, Tuple[str, str]] = field(default_factory=dict)   # label -> (rxn_smarts, adjacent_smarts)
    notes: List[str] = field(default_factory=list)

    def merge(self, other: "TemplateSet") -> None:
        for key, value in other.residue_templates.items():
            if key in self.residue_templates and self.residue_templates[key] != value:
                raise ValueError(f"Two different templates named {key}.")
            self.residue_templates[key] = value
        for key, value in other.padders.items():
            if key in self.padders and self.padders[key] != value:
                raise ValueError(f"Two different padders labelled {key}.")
            self.padders[key] = value
        self.notes.extend(other.notes)

    def to_meeko(self):
        """Meeko's defaults plus these. Padders are assigned directly: in Meeko 0.8.0,
        ResidueChemTemplates.add_dict builds padders from the wrong variable."""
        from meeko import ResidueChemTemplates, ResiduePadder
        chem = ResidueChemTemplates.create_from_defaults()
        if self.residue_templates:
            chem.add_dict({"residue_templates": self.residue_templates}, overwrite=True)
        for label, (rxn, adjacent) in self.padders.items():
            chem.padders[label] = ResiduePadder(rxn, adjacent)
        return chem


# ------------------------------------------------------------------ ions ---

def ion_residue_name(element: str, charge: int) -> str:
    """FE + 2 -> FE2. At most three characters, as the PDB residue field allows."""
    if charge < 0 or charge > 9:
        raise ValueError(f"Metal charge {charge:+d} cannot be written as a residue name.")
    name = f"{element.upper()}{charge}"
    if len(name) > 3:
        raise ValueError(f"{element} with charge {charge} gives residue name {name!r}, longer than 3 characters.")
    return name


def ion_template(element: str, charge: int, atom_name: Optional[str] = None) -> TemplateSet:
    """Template for a single metal ion. ``atom_name`` is the ion's name in the structure
    (ZN for a zinc ion, FE for a heme iron, FE1 for one iron of a cluster); default: the element."""
    symbol = element[0].upper() + element[1:].lower()
    atom = Chem.Atom(symbol)
    if not is_metal(atom):
        raise ValueError(f"{symbol} is not a metal.")
    atom.SetFormalCharge(charge)
    atom.SetNoImplicit(True)
    ion = Chem.RWMol()
    ion.AddAtom(atom)
    smiles = Chem.MolToSmiles(ion.GetMol())
    name = ion_residue_name(element, charge)
    return TemplateSet(residue_templates={name: {"smiles": smiles, "atom_name": [atom_name or element.upper()],
                                                 "link_labels": {}}})


# ------------------------------------------------------------- templates ---

def template_dict_from_mol(mol: Chem.Mol, link_labels: Optional[Dict[str, str]] = None) -> dict:
    """Meeko template: SMILES, atom names in SMILES order, link labels by name.

    Hydrogens must be atoms. A bracket hydrogen count such as [NH] reads to
    Meeko as an implicit hydrogen and it refuses the residue, so every atom is
    written with its hydrogen atoms only. An atom other than a link atom that
    is short of hydrogens means the chemistry is incomplete and is refused.
    """
    link_labels = dict(link_labels or {})
    names_present = {a.GetProp("name") for a in mol.GetAtoms()}
    for atom_name, label in link_labels.items():
        if atom_name not in names_present:
            raise KeyError(f"No atom {atom_name} for link label {label}.")
    probe = Chem.Mol(mol)
    for atom in probe.GetAtoms():
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    probe.UpdatePropertyCache(strict=False)
    short = [f"{a.GetProp('name')} ({a.GetNumImplicitHs()} H)" for a in probe.GetAtoms()
             if a.GetNumImplicitHs() and a.GetProp("name") not in link_labels and not is_metal(a)]
    if short:
        raise ValueError("Template atoms short of hydrogens that are not link atoms: " + ", ".join(short))
    written = Chem.Mol(mol)
    for atom in written.GetAtoms():
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(True)
    written.UpdatePropertyCache(strict=False)
    smiles = Chem.MolToSmiles(written, allHsExplicit=True)
    order = list(map(int, written.GetProp("_smilesAtomOutputOrder").strip("[],").split(",")))
    names = [written.GetAtomWithIdx(i).GetProp("name") for i in order]
    labels = {str(names.index(atom_name)): label for atom_name, label in link_labels.items()}
    return {"smiles": smiles, "atom_name": names, "link_labels": labels}


def split_metals(mol: Chem.Mol) -> Tuple[Chem.Mol, List[Tuple[str, str]]]:
    """Remove metal atoms from a residue; return the rest and [(atom name, element)] removed.

    Coordination (dative) bonds go with the metal. The formal charges of the
    remaining atoms are kept (a porphyrin's pyrrolide N- stay charged).
    """
    rw = Chem.RWMol(mol)
    removed = []
    for atom in sorted((a for a in mol.GetAtoms() if is_metal(a)), key=lambda a: a.GetIdx(), reverse=True):
        removed.append((atom.GetProp("name"), atom.GetSymbol()))
        rw.RemoveAtom(atom.GetIdx())
    result = rw.GetMol()
    Chem.SanitizeMol(result)
    return result, list(reversed(removed))


def remove_link_hydrogens(mol: Chem.Mol, link_atoms: Dict[str, int]) -> Tuple[Chem.Mol, List[str]]:
    """Remove ``count`` hydrogens from each named link atom (one per bond it forms to a neighbour).

    Returns the molecule and the names of the hydrogens removed. A link atom
    without enough hydrogens is refused: forming that bond changes more than a
    hydrogen (a leaving group, a bond order) and needs an explicit edit.
    """
    rw = Chem.RWMol(mol)
    index = {a.GetProp("name"): a.GetIdx() for a in rw.GetAtoms()}
    to_remove, removed_names = [], []
    for name, count in link_atoms.items():
        if name not in index:
            raise KeyError(f"No link atom named {name}.")
        hydrogens = sorted((n for n in rw.GetAtomWithIdx(index[name]).GetNeighbors() if n.GetAtomicNum() == 1),
                           key=lambda n: n.GetProp("name"))
        if len(hydrogens) < count:
            raise ValueError(f"Link atom {name} has {len(hydrogens)} hydrogen(s) but forms {count} bond(s) "
                             "to other residues; the linkage needs an explicit chemistry edit.")
        for hydrogen in hydrogens[-count:]:          # equivalent H; the last by name goes (reported)
            to_remove.append(hydrogen.GetIdx())
            removed_names.append(hydrogen.GetProp("name"))
    for i in sorted(to_remove, reverse=True):
        rw.RemoveAtom(i)
    result = rw.GetMol()
    Chem.SanitizeMol(result)
    return result, removed_names


def unique_environment_smarts(mol: Chem.Mol, atom_index: int, max_radius: int = 16) -> Tuple[str, int]:
    """The smallest atom environment around ``atom_index`` whose SMARTS matches ``mol`` once.

    The atom itself carries map number 1, the others 2, 3, ... Raises when even
    the largest environment matches more than once, which happens only for
    atoms related by the molecule's own symmetry.
    """
    for radius in range(0, max_radius + 1):
        if radius == 0:
            atoms = [atom_index]
        else:
            env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atom_index)
            atoms = sorted({atom_index} | {mol.GetBondWithIdx(b).GetBeginAtomIdx() for b in env}
                           | {mol.GetBondWithIdx(b).GetEndAtomIdx() for b in env})
        mapped = Chem.RWMol(mol)
        for a in mapped.GetAtoms():
            a.SetAtomMapNum(0)
        for n, i in enumerate(atoms, start=1):
            mapped.GetAtomWithIdx(i).SetAtomMapNum(1 if i == atom_index else n + 1)
        smarts = Chem.MolFragmentToSmarts(mapped, atomsToUse=atoms)
        if len(mol.GetSubstructMatches(Chem.MolFromSmarts(smarts), uniquify=True, useChirality=False)) == 1:
            return smarts, radius
    raise ValueError(f"Atom {atom_index + 1} has no unique environment within radius {max_radius} "
                     "(it is symmetry-equivalent to another atom).")


def link_padder(residue_mol: Chem.Mol, link_atom: str, neighbour_atoms: str) -> Tuple[str, str]:
    """Padder for a link: attach ``neighbour_atoms`` (SMARTS, first atom bonded to the link atom,
    mapped from 101, e.g. '[S:101][C:102]') to ``link_atom`` of ``residue_mol``.

    Returns (reaction SMARTS, adjacent-residue SMARTS). The template molecule
    must carry its hydrogens: Meeko pads the residue as built from its template.
    """
    index = {a.GetProp("name"): a.GetIdx() for a in residue_mol.GetAtoms()}
    reactant, _ = unique_environment_smarts(residue_mol, index[link_atom])
    query = Chem.RWMol(Chem.MolFromSmarts(reactant))
    anchor = next(a.GetIdx() for a in query.GetAtoms() if a.GetAtomMapNum() == 1)
    added = Chem.MolFromSmarts(neighbour_atoms)
    offset = query.GetNumAtoms()
    query.InsertMol(added)
    query.AddBond(anchor, offset, Chem.BondType.SINGLE)
    return f"{reactant}>>{Chem.MolToSmarts(query)}", neighbour_atoms


def neighbour_smarts(neighbour_mol: Chem.Mol, neighbour_link_atom: str) -> str:
    """SMARTS for the atoms that cap the other residue at a link, positions copied from the neighbour.

    Meeko requires this pattern to match the neighbour exactly once among
    matches containing its link atom. With a single heavy neighbour beyond the
    link atom (Cys SG-CB) the cap is the pair, '[#16:101][#6:102]'; with more
    (heme c CAB, bonded to C3B and CBB) a pair would match twice, so the cap is
    the link atom alone.
    """
    index = {a.GetProp("name"): a.GetIdx() for a in neighbour_mol.GetAtoms()}
    atom = neighbour_mol.GetAtomWithIdx(index[neighbour_link_atom])
    heavy = [n for n in atom.GetNeighbors() if n.GetAtomicNum() > 1 and not is_metal(n)]
    if len(heavy) == 1:
        return f"[#{atom.GetAtomicNum()}:101][#{heavy[0].GetAtomicNum()}:102]"
    return f"[#{atom.GetAtomicNum()}:101]"


def meeko_template_mol(template_key: str) -> Chem.Mol:
    """One of Meeko's built-in templates as an RDKit molecule with its atom names."""
    from meeko import ResidueChemTemplates
    template = ResidueChemTemplates.create_from_defaults().residue_templates[template_key]
    mol = Chem.Mol(template.mol)
    for atom, name in zip(mol.GetAtoms(), template.atom_names):
        atom.SetProp("name", name)
    return mol, {int(k): v for k, v in (template.link_labels or {}).items()}


def linked_variant(template_key: str, link_atom: str, label: str, variant_key: str) -> dict:
    """A built-in template with one hydrogen fewer on ``link_atom`` and a link label there,
    keeping the template's own link labels (C-term/N-term for an amino acid)."""
    mol, labels = meeko_template_mol(template_key)
    kept = {mol.GetAtomWithIdx(i).GetProp("name"): lab for i, lab in labels.items()}
    stripped, _ = remove_link_hydrogens(mol, {link_atom: 1})
    kept[link_atom] = label
    return template_dict_from_mol(stripped, kept)
