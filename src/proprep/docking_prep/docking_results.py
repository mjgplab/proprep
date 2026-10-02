"""
Docked poses back to molecules: the ligand with its own atom names, the
flexible side chains, RMSD to a reference, and an SDF with the energies.

Meeko rebuilds each pose from the SMILES stored in the PDBQT, which carries
no atom names. The pose coordinates are therefore mapped onto the ligand that
went in (by substructure match, hydrogens included), so a pose keeps the names
the user saw and edited. Where the molecule has symmetry the match picks one
of the equivalent assignments; the geometry is the same either way.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

from rdkit import Chem
from rdkit.Chem import rdMolAlign


@dataclass
class Pose:
    rank: int
    ligand: Chem.Mol                          # the input ligand, with this pose's coordinates
    flexible_residues: List[Chem.Mol]         # side chains as docked (empty for a rigid receptor)
    energies: Dict[str, float]
    rmsd_to_reference: Optional[float] = None


def poses_from_run(poses_pdbqt: str, ligand: Chem.Mol, energies: Sequence[Dict[str, float]]) -> List[Pose]:
    """One Pose per model in ``poses_pdbqt``, in Vina's order (best first)."""
    from meeko import PDBQTMolecule, RDKitMolCreate
    pdbqt_mol = PDBQTMolecule(poses_pdbqt, skip_typing=True)
    exported = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol, keep_flexres=True)
    if not exported or exported[0] is None:
        raise RuntimeError("Meeko could not rebuild the docked ligand from the poses.")
    docked, flexres = exported[0], [m for m in exported[1:] if m is not None]
    match = _atom_match(ligand, docked)
    if len(match) != docked.GetNumAtoms() or docked.GetNumAtoms() != ligand.GetNumAtoms():
        raise RuntimeError("The docked ligand does not match the ligand that was prepared "
                           f"({docked.GetNumAtoms()} vs {ligand.GetNumAtoms()} atoms).")
    n_conformers = docked.GetNumConformers()
    if n_conformers != len(energies):
        raise RuntimeError(f"{n_conformers} poses but {len(energies)} energy rows.")
    poses = []
    for rank, conformer in enumerate(docked.GetConformers()):
        mol = Chem.Mol(ligand)
        mol.RemoveAllConformers()
        placed = Chem.Conformer(ligand.GetNumAtoms())
        for docked_index, ligand_index in enumerate(match):
            placed.SetAtomPosition(ligand_index, conformer.GetAtomPosition(docked_index))
        placed.Set3D(True)
        mol.AddConformer(placed, assignId=True)
        sides = []
        for side_chain in flexres:
            single = Chem.Mol(side_chain)
            keep = side_chain.GetConformers()[rank] if side_chain.GetNumConformers() > rank else None
            single.RemoveAllConformers()
            if keep is not None:
                single.AddConformer(Chem.Conformer(keep), assignId=True)
            sides.append(single)
        poses.append(Pose(rank=rank + 1, ligand=mol, flexible_residues=sides, energies=dict(energies[rank])))
    return poses


def _topology(mol: Chem.Mol) -> Chem.Mol:
    """The molecule as elements and connections only: every bond single, nothing aromatic, no charges."""
    rw = Chem.RWMol(mol)
    for bond in rw.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    for atom in rw.GetAtoms():
        atom.SetIsAromatic(False)
        atom.SetFormalCharge(0)
        atom.SetNoImplicit(True)
    out = rw.GetMol()
    out.UpdatePropertyCache(strict=False)
    return out


def _atom_match(ligand: Chem.Mol, docked: Chem.Mol) -> tuple:
    """Which ligand atom each docked atom is. Meeko rebuilds the pose from SMILES and may draw a
    fused aromatic system with other bond orders (caffeine), so when the exact match fails the
    atoms are matched by element and connectivity alone, which is all a pose needs."""
    match = ligand.GetSubstructMatch(docked)
    if len(match) == docked.GetNumAtoms():
        return match
    return _topology(ligand).GetSubstructMatch(_topology(docked))


def heavy_atom_rmsd(pose: Chem.Mol, reference: Chem.Mol) -> float:
    """Symmetry-aware heavy-atom RMSD in place (no superposition), in Angstrom."""
    return float(rdMolAlign.CalcRMS(Chem.RemoveHs(pose), Chem.RemoveHs(reference)))


def write_poses_sdf(path: str, poses: Sequence[Pose], properties: Dict[str, object]) -> None:
    """All poses in one SDF, best first; each record carries its rank, energies and ``properties``."""
    writer = Chem.SDWriter(path)
    try:
        for pose in poses:
            mol = Chem.Mol(pose.ligand)
            mol.SetProp("_Name", f"pose {pose.rank}")
            mol.SetIntProp("rank", pose.rank)
            for name, value in pose.energies.items():
                mol.SetProp(f"energy {name} (kcal/mol)", f"{float(value):.3f}")    # Vina reports 3 decimals
            if pose.rmsd_to_reference is not None:
                mol.SetProp("RMSD to reference (A)", f"{float(pose.rmsd_to_reference):.3f}")
            for name, value in properties.items():
                mol.SetProp(str(name), str(value))
            writer.write(mol)
    finally:
        writer.close()
