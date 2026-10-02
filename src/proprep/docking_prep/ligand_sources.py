"""
Ligand intake: every route ends as one RDKit molecule with hydrogens, bond
orders, formal charges, named atoms and 3D coordinates.

Routes:

* ``from_smiles``    -- a SMILES string; 3D coordinates from RDKit ETKDG with an
                        explicit seed, optionally MMFF94-minimised.
* ``from_file``      -- an SDF or mol2 file (one record chosen explicitly when
                        the file holds several).
* ``from_structure_residue`` -- a HETATM residue of a loaded structure: crystal
                        coordinates for the heavy atoms, bond orders, charges and
                        hydrogens from the residue's CCD entry.

Nothing here decides chemistry silently. Whatever a route had to settle on its
own -- hydrogens added from valence, a stereocentre the SMILES left open,
crystal hydrogens replaced by the CCD's -- is written to ``LigandBuild.notes``
for the caller to show, and anything that cannot be settled safely (missing
atoms, an unresolved CCD problem, several records and no choice) raises with
the specifics.
"""

from __future__ import annotations

import os
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

from .ccd_chemistry import CcdComponent, find_valence_problems

OPTIMIZERS = ("mmff94", "none")


@dataclass
class LigandBuild:
    """A ligand ready for docking preparation, with how it was obtained."""

    mol: Chem.Mol
    source: Dict[str, object]                 # everything needed to rebuild it (session replay)
    notes: List[str] = field(default_factory=list)

    @property
    def net_charge(self) -> int:
        return Chem.GetFormalCharge(self.mol)


# ---------------------------------------------------------------- naming ----

def assign_atom_names(mol: Chem.Mol) -> int:
    """Give every atom a unique ``name`` property; keep names it already has.

    Existing names come from a ``name`` property, a mol2 atom name, or PDB
    residue info. Unnamed atoms get element + running number (C1, C2, O1...).
    Returns the number of atoms that were named here.
    """
    taken = set()
    for atom in mol.GetAtoms():
        for source in ("name", "_TriposAtomName"):
            if atom.HasProp(source) and atom.GetProp(source).strip():
                atom.SetProp("name", atom.GetProp(source).strip())
                break
        else:
            info = atom.GetPDBResidueInfo()
            if info is not None and info.GetName().strip():
                atom.SetProp("name", info.GetName().strip())
        if atom.HasProp("name"):
            if atom.GetProp("name") in taken:
                atom.ClearProp("name")          # duplicate: rename below
            else:
                taken.add(atom.GetProp("name"))
    counters: Counter = Counter()
    named = 0
    for atom in mol.GetAtoms():
        if atom.HasProp("name"):
            continue
        symbol = atom.GetSymbol()
        while True:
            counters[symbol] += 1
            candidate = f"{symbol.upper() if len(symbol) == 1 else symbol}{counters[symbol]}"
            if candidate not in taken and len(candidate) <= 4:
                break
        atom.SetProp("name", candidate)
        taken.add(candidate)
        named += 1
    return named


def _stereo_report(before: Chem.Mol, after: Chem.Mol) -> List[str]:
    """Stereo elements the input left unspecified, and what the 3D build chose."""
    notes = []
    unspecified = [info for info in Chem.FindPotentialStereo(before)
                   if info.specified == Chem.StereoSpecified.Unspecified]
    if not unspecified:
        return notes
    Chem.AssignStereochemistryFrom3D(after)
    from rdkit.Chem import rdCIPLabeler
    rdCIPLabeler.AssignCIPLabels(after)
    for info in unspecified:
        if info.type == Chem.StereoType.Atom_Tetrahedral:
            atom = after.GetAtomWithIdx(info.centeredOn)
            label = atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else "?"
            notes.append(f"Stereocentre {atom.GetProp('name')} was not specified; "
                         f"the 3D build made it {label}.")
        elif info.type == Chem.StereoType.Bond_Double:
            bond = after.GetBondWithIdx(info.centeredOn)
            label = bond.GetProp("_CIPCode") if bond.HasProp("_CIPCode") else "?"
            a, b = bond.GetBeginAtom().GetProp("name"), bond.GetEndAtom().GetProp("name")
            notes.append(f"Double bond {a}={b} geometry was not specified; the 3D build made it {label}.")
        else:
            notes.append(f"A stereo element of type {info.type} was not specified.")
    return notes


# ----------------------------------------------------------------- SMILES ---

def from_smiles(smiles: str, embed_seed: int, optimizer: str, residue_name: str = "LIG") -> LigandBuild:
    """Build a 3D ligand from SMILES.

    ``embed_seed`` makes the ETKDG coordinates reproducible and is required.
    ``optimizer`` is "mmff94" (minimise with MMFF94; fails if MMFF has no
    parameters for an atom) or "none". Hydrogens are added from the SMILES'
    valences at the formal charges written.
    """
    if optimizer not in OPTIMIZERS:
        raise ValueError(f"optimizer must be one of {', '.join(OPTIMIZERS)}, got {optimizer!r}.")
    parsed = Chem.MolFromSmiles(smiles)
    if parsed is None:
        raise ValueError(f"RDKit could not read the SMILES {smiles!r}.")
    before = Chem.AddHs(parsed)
    generated = assign_atom_names(before)
    mol = Chem.Mol(before)
    params = AllChem.ETKDGv3()
    params.randomSeed = int(embed_seed)
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise RuntimeError(f"RDKit ETKDG could not generate 3D coordinates for {smiles!r} "
                           f"(seed {embed_seed}).")
    notes = [f"3D coordinates from RDKit ETKDGv3, seed {embed_seed}.",
             f"Atom names generated for all {generated} atoms (element + number)."]
    notes += _optimize(mol, optimizer)
    notes += _stereo_report(before, mol)
    _set_residue_info(mol, residue_name)
    return LigandBuild(mol=mol, notes=notes, source={
        "kind": "smiles", "smiles": smiles, "embed_seed": int(embed_seed), "optimizer": optimizer,
        "residue_name": residue_name})


def _optimize(mol: Chem.Mol, optimizer: str) -> List[str]:
    if optimizer == "none":
        return ["Geometry not minimised."]
    if not AllChem.MMFFHasAllMoleculeParams(mol):
        raise ValueError("MMFF94 has no parameters for at least one atom of this ligand; "
                         "choose no minimisation instead.")
    status = AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    if status != 0:
        return ["MMFF94 minimisation stopped at 2000 iterations before converging."]
    return ["Geometry minimised with MMFF94 (converged)."]


# ------------------------------------------------------------------- files ---

def _mol2_records(text: str) -> List[str]:
    parts = text.split("@<TRIPOS>MOLECULE")
    return ["@<TRIPOS>MOLECULE" + part for part in parts[1:]]


def from_file(path: str, record: Optional[int] = None, residue_name: str = "LIG") -> LigandBuild:
    """Read a ligand from an SDF (.sdf/.mol) or mol2 file.

    ``record`` is the 0-based record to use; it is required when the file holds
    more than one. Coordinates must be 3D. When the file carries no hydrogens
    they are added from valence at the formal charges written, and a note says
    how many.
    """
    extension = os.path.splitext(path)[1].lower()
    if extension in (".sdf", ".mol", ".sd"):
        supplier = Chem.SDMolSupplier(path, removeHs=False, sanitize=True)
        records = [supplier[i] for i in range(len(supplier))]
        titles = [m.GetProp("_Name") if m is not None and m.HasProp("_Name") else "" for m in records]
    elif extension == ".mol2":
        with open(path) as handle:
            blocks = _mol2_records(handle.read())
        records = [Chem.MolFromMol2Block(block, removeHs=False, sanitize=True) for block in blocks]
        titles = [block.splitlines()[1].strip() if len(block.splitlines()) > 1 else "" for block in blocks]
    else:
        raise ValueError(f"Unsupported ligand file type {extension!r}; use .sdf, .mol or .mol2.")
    if not records:
        raise ValueError(f"{os.path.basename(path)} contains no molecules.")
    if record is None:
        if len(records) > 1:
            listing = "; ".join(f"{i}: {t or '(untitled)'}" for i, t in enumerate(titles))
            raise ValueError(f"{os.path.basename(path)} holds {len(records)} molecules; choose one ({listing}).")
        record = 0
    if not 0 <= record < len(records):
        raise IndexError(f"Record {record} requested but {os.path.basename(path)} holds {len(records)}.")
    mol = records[record]
    if mol is None:
        raise ValueError(f"RDKit could not read record {record} of {os.path.basename(path)}.")
    if mol.GetNumConformers() == 0 or not mol.GetConformer().Is3D() or np.allclose(
            mol.GetConformer().GetPositions()[:, 2], 0.0):
        raise ValueError(f"Record {record} of {os.path.basename(path)} has no 3D coordinates.")
    notes = [f"Read record {record} ({titles[record] or 'untitled'}) of {os.path.basename(path)}."]
    n_h = sum(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms())
    if n_h == 0:
        mol = Chem.AddHs(mol, addCoords=True)
        added = sum(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms())
        notes.append(f"The file has no hydrogens; {added} were added from valence at the formal "
                     f"charges written (net charge {Chem.GetFormalCharge(mol):+d}).")
    generated = assign_atom_names(mol)
    if generated:
        notes.append(f"Atom names generated for {generated} of {mol.GetNumAtoms()} atoms (element + number).")
    _set_residue_info(mol, residue_name)
    return LigandBuild(mol=mol, notes=notes, source={
        "kind": "file", "path": os.path.abspath(path), "record": record, "residue_name": residue_name})


# ------------------------------------------------------- structure residue ---

@dataclass
class _AtomRecord:
    name: str
    element: str
    altloc: str
    xyz: Tuple[float, float, float]


def _residue_atoms(pdb_path: str, chain: str, resseq: int, icode: str, resname: str) -> List[_AtomRecord]:
    from proprep.utils.pdb_format import element_from_name_field
    atoms = []
    with open(pdb_path) as handle:
        for line in handle:
            if line.startswith("ENDMDL"):
                break                          # first model only
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if (line[21] != chain or int(line[22:26]) != resseq or line[26].strip() != icode.strip()
                    or line[17:20].strip() != resname):
                continue
            element = line[76:78].strip() if len(line) >= 78 else ""
            if not element:
                element = element_from_name_field(line[12:16])
            atoms.append(_AtomRecord(line[12:16].strip(), element.capitalize(), line[16].strip(),
                                     (float(line[30:38]), float(line[38:46]), float(line[46:54]))))
    return atoms


def from_structure_residue(pdb_path: str, chain: str, resseq: int, icode: str, resname: str,
                           component: CcdComponent, altloc: Optional[str] = None,
                           altloc_fills: Optional[Dict[str, str]] = None) -> LigandBuild:
    """A HETATM residue with its CCD chemistry on the crystal coordinates.

    ``component`` is the residue's CCD entry, with any problems already
    resolved by edits (see ``chemistry_edits``); it defines bond orders,
    formal charges and which heavy atom carries which hydrogens. Heavy atoms
    are matched to the CCD by name. Hydrogens in the file are not used; the
    CCD's hydrogens are placed on the crystal heavy atoms. ``altloc`` is
    required when the residue has alternate locations; ``altloc_fills`` names
    the atoms it does not model and the alternate each is taken from
    (``altloc_picker.fill_plan``).
    """
    if component.problems:
        listing = "; ".join(p.describe() for p in component.problems)
        raise ValueError(f"The CCD entry {component.code} has unresolved problems: {listing}")
    records = _residue_atoms(pdb_path, chain, resseq, icode, resname)
    where = f"{resname} {chain}{resseq}{icode.strip()}"
    if not records:
        raise LookupError(f"No atoms for {where} in {os.path.basename(pdb_path)}.")
    altlocs = sorted({r.altloc for r in records if r.altloc})
    if altlocs and altloc is None:
        raise ValueError(f"{where} has alternate locations {', '.join(altlocs)}; choose one.")
    if altloc is not None and altloc not in altlocs:
        raise ValueError(f"{where} has no alternate location {altloc!r} (has {', '.join(altlocs) or 'none'}).")
    from proprep.structure_prep.altloc_picker import keeps_atom
    chosen = [r for r in records if keeps_atom(r.name, r.altloc, altloc, altloc_fills)]
    crystal_h = [r for r in chosen if r.element in ("H", "D")]
    heavy = {r.name: r for r in chosen if r.element not in ("H", "D")}

    ccd = component.mol
    ccd_heavy = [a for a in ccd.GetAtoms() if a.GetAtomicNum() != 1]
    ccd_names = {a.GetProp("name") for a in ccd_heavy}
    missing = [a.GetProp("name") for a in ccd_heavy if a.GetProp("name") not in heavy]
    extra = sorted(set(heavy) - ccd_names)
    if missing or extra:
        parts = []
        if missing:
            parts.append(f"missing from the structure: {', '.join(missing)}")
        if extra:
            parts.append(f"not in CCD {component.code}: {', '.join(extra)}")
        raise ValueError(f"{where} does not match its CCD entry ({'; '.join(parts)}).")
    for a in ccd_heavy:
        if a.GetSymbol() != heavy[a.GetProp("name")].element:
            raise ValueError(f"{where} atom {a.GetProp('name')} is {heavy[a.GetProp('name')].element} "
                             f"in the structure but {a.GetSymbol()} in CCD {component.code}.")

    # heavy-atom skeleton with the CCD's hydrogen count on each atom
    hydrogen_names: Dict[str, List[str]] = defaultdict(list)
    for a in ccd.GetAtoms():
        if a.GetAtomicNum() == 1:
            parent = a.GetNeighbors()[0].GetProp("name")
            hydrogen_names[parent].append(a.GetProp("name"))
    skeleton = Chem.RWMol(ccd)
    for index in sorted((a.GetIdx() for a in ccd.GetAtoms() if a.GetAtomicNum() == 1), reverse=True):
        skeleton.RemoveAtom(index)
    for atom in skeleton.GetAtoms():
        atom.SetNumExplicitHs(len(hydrogen_names.get(atom.GetProp("name"), [])))
        atom.SetNoImplicit(True)
    skeleton.RemoveAllConformers()
    conformer = Chem.Conformer(skeleton.GetNumAtoms())
    for atom in skeleton.GetAtoms():
        conformer.SetAtomPosition(atom.GetIdx(), Point3D(*heavy[atom.GetProp("name")].xyz))
    skeleton.AddConformer(conformer, assignId=True)
    mol = skeleton.GetMol()
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol, addCoords=True)

    pending = {k: list(v) for k, v in hydrogen_names.items()}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1 and not atom.HasProp("name"):
            parent = atom.GetNeighbors()[0].GetProp("name")
            atom.SetProp("name", pending[parent].pop(0))
    problems = find_valence_problems(mol)
    if problems:
        raise RuntimeError("Building on the crystal coordinates produced invalid chemistry: "
                           + "; ".join(p.describe() for p in problems))

    notes = [f"{len(heavy)} heavy atoms from {where} in {os.path.basename(pdb_path)}"
             + (f", alternate location {altloc}" if altloc else "")
             + (f" ({', '.join(f'{n} from {l}' for n, l in sorted(altloc_fills.items()))}, which {altloc} "
                "does not model)" if altloc_fills else "") + ".",
             f"Bond orders, formal charges and hydrogens from CCD {component.code} "
             f"(net charge {Chem.GetFormalCharge(mol):+d})."]
    if crystal_h:
        notes.append(f"{len(crystal_h)} hydrogens in the structure were replaced by the CCD's.")
    notes += _crystal_vs_ccd_stereo(mol, component)
    if len(resname) > 3:
        notes.append(f"The PDB/PDBQT residue field holds 3 characters; {resname} is written as {resname[:3]}.")
    _set_residue_info(mol, resname, chain=chain, resseq=resseq, icode=icode)
    return LigandBuild(mol=mol, notes=notes, source={
        "kind": "structure_residue", "path": os.path.abspath(pdb_path), "chain": chain,
        "resseq": resseq, "icode": icode, "resname": resname, "altloc": altloc,
        "altloc_fills": dict(altloc_fills or {}),
        "ccd_code": component.code})


def _crystal_vs_ccd_stereo(mol: Chem.Mol, component: CcdComponent) -> List[str]:
    """Stereocentres whose configuration in the structure differs from the CCD's ideal geometry."""
    if component.coordinates == "none":
        return []
    from rdkit.Chem import rdCIPLabeler
    crystal = Chem.Mol(mol)
    ideal = Chem.Mol(component.mol)
    for m in (crystal, ideal):
        Chem.AssignStereochemistryFrom3D(m)
        rdCIPLabeler.AssignCIPLabels(m)
    ideal_labels = {a.GetProp("name"): a.GetProp("_CIPCode") for a in ideal.GetAtoms() if a.HasProp("_CIPCode")}
    notes = []
    for atom in crystal.GetAtoms():
        name = atom.GetProp("name")
        if atom.HasProp("_CIPCode") and name in ideal_labels and atom.GetProp("_CIPCode") != ideal_labels[name]:
            notes.append(f"Stereocentre {name} is {atom.GetProp('_CIPCode')} in the structure but "
                         f"{ideal_labels[name]} in CCD {component.code}; the structure's geometry is used.")
    return notes


def _set_residue_info(mol: Chem.Mol, resname: str, chain: str = "A", resseq: int = 1, icode: str = "") -> None:
    """PDB residue info from the atom names, so names survive into PDBQT and back."""
    from proprep.utils.pdb_format import atom_name_field
    for atom in mol.GetAtoms():
        info = Chem.AtomPDBResidueInfo()
        info.SetName(atom_name_field(atom.GetProp("name"), atom.GetSymbol()))
        info.SetResidueName(resname[:3].rjust(3))
        info.SetResidueNumber(int(resseq))
        info.SetChainId(chain)
        info.SetInsertionCode(icode or " ")
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)
