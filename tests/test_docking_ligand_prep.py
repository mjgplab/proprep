#!/usr/bin/env python3
"""
Ligand PDBQT preparation with Meeko and per-bond torsion choices.

Meeko's own rigidification selects bonds by SMARTS, which cannot separate
symmetry-equivalent bonds; the module applies the user's choice to exact atom
pairs instead. These tests pin that, the amide opt-in, the labelling of
conjugated single bonds (Meeko rotates ester C(=O)-O bonds by default), and
macrocycle handling, whose glue atom types (CG0, G0) run past the two-column
PDBQT type field and whose TORSDOF counts the unopened molecule.

Run with: pytest tests/test_docking_ligand_prep.py
"""

import os

import pytest

pytest.importorskip("rdkit", reason="RDKit is required for docking ligand preparation")
pytest.importorskip("meeko", reason="Meeko is required for docking ligand preparation")

from proprep.docking_prep import ligand_sources as ls  # noqa: E402
from proprep.docking_prep.ccd_chemistry import component_from_block  # noqa: E402
from proprep.docking_prep.ligand_prep import prepare_ligand, torsion_bonds  # noqa: E402

DATA = os.path.join(os.path.dirname(__file__), "data", "docking")

# two identical ester tails on a glycerol-like core, an amide, and a conjugated diene
DIESTER = "CCCCC(=O)OCC(COC(=O)CCCC)NC(=O)C=CC=CC"


@pytest.fixture(scope="module")
def diester():
    return ls.from_smiles(DIESTER, embed_seed=5, optimizer="mmff94").mol


def _by_kind(bonds, kind):
    return [b for b in bonds if b.kind == kind]


def test_one_of_two_symmetric_bonds_can_be_held_rigid(diester):
    bonds = torsion_bonds(diester, rigid_macrocycles=False)
    ester_o_c = [b for b in _by_kind(bonds, "single") if {n[0] for n in b.atoms} == {"C", "O"}]
    assert len(ester_o_c) == 2
    prepared = prepare_ligand(diester, rigid=[ester_o_c[0].atoms], rigid_macrocycles=False)
    rotating = {frozenset(pair) for pair in prepared.rotatable}
    assert frozenset(ester_o_c[0].atoms) not in rotating
    assert frozenset(ester_o_c[1].atoms) in rotating


def test_amide_rotates_only_when_chosen(diester):
    bonds = torsion_bonds(diester, rigid_macrocycles=False)
    amides = _by_kind(bonds, "amide")
    assert len(amides) == 1 and not amides[0].rotatable_by_default
    default = prepare_ligand(diester, rigid_macrocycles=False)
    chosen = prepare_ligand(diester, rotatable_amides=[amides[0].atoms], rigid_macrocycles=False)
    assert len(chosen.rotatable) == len(default.rotatable) + 1
    with pytest.raises(ValueError, match="does not rotate by default"):
        prepare_ligand(diester, rigid=[amides[0].atoms], rigid_macrocycles=False)


def test_conjugated_single_bonds_are_labelled_and_can_be_frozen_together(diester):
    bonds = torsion_bonds(diester, rigid_macrocycles=False)
    conjugated = _by_kind(bonds, "conjugated single")
    assert len(conjugated) == 4          # two ester C(=O)-O, two diene single bonds
    default = prepare_ligand(diester, rigid_macrocycles=False)
    frozen = prepare_ligand(diester, rigid=[b.atoms for b in conjugated], rigid_macrocycles=False)
    assert len(frozen.rotatable) == len(default.rotatable) - 4


def test_bonds_that_cannot_rotate_are_refused_by_name(diester):
    with pytest.raises(ValueError, match="two different atom names"):
        prepare_ligand(diester, rigid=[("C1", "C1")], rigid_macrocycles=False)
    with pytest.raises(KeyError, match="NOPE"):
        prepare_ligand(diester, rigid=[("C1", "NOPE")], rigid_macrocycles=False)
    with pytest.raises(ValueError, match="not a bond that can rotate"):
        prepare_ligand(diester, rigid=[("C1", "C17")], rigid_macrocycles=False)


def test_macrocycle_opening_is_explicit_and_reported():
    lactone = ls.from_smiles("O=C1CCCCCCCCCCCO1", embed_seed=3, optimizer="mmff94").mol
    opened = prepare_ligand(lactone, rigid_macrocycles=False)
    assert opened.atom_types.get("CG0") == 2 and opened.atom_types.get("G0") == 2
    assert len(opened.rotatable) == 13 and opened.torsdof == 0
    assert any("macrocycle was opened" in note for note in opened.notes)
    assert all(b.kind == "macrocycle ring" for b in torsion_bonds(lactone, rigid_macrocycles=False))
    rigid = prepare_ligand(lactone, rigid_macrocycles=True)
    assert "CG0" not in rigid.atom_types and rigid.rotatable == []
    assert torsion_bonds(lactone, rigid_macrocycles=True) == []


def test_crystal_ligand_keeps_its_ccd_names_in_the_pdbqt(tmp_path):
    with open(os.path.join(DATA, "ccd_BTN.cif.txt")) as handle:
        btn = component_from_block(handle.read(), "BTN")
    pdb = tmp_path / "1STP_BTN.pdb"
    with open(os.path.join(DATA, "1STP_BTN.pdb.txt")) as handle:
        pdb.write_text(handle.read())
    build = ls.from_structure_residue(str(pdb), "A", 300, "", "BTN", btn)
    prepared = prepare_ligand(build.mol, rigid_macrocycles=False)
    names = {line[12:16].strip() for line in prepared.pdbqt.splitlines() if line.startswith(("ATOM", "HETATM"))}
    heavy = {a.GetProp("name") for a in btn.mol.GetAtoms() if a.GetAtomicNum() > 1}
    assert heavy <= names
    assert abs(prepared.gasteiger_charge_sum - build.net_charge) < 0.01
    assert prepared.atom_types.get("SA") == 1


def _heme_ligand(tmp_path):
    from proprep.docking_prep.chemistry_edits import ChemistryEdit, apply_edits
    with open(os.path.join(DATA, "ccd_HEM.cif.txt")) as handle:
        hem = component_from_block(handle.read(), "HEM")
    hem.mol, hem.problems = apply_edits(hem.mol, [ChemistryEdit("remove_atom", ("H2A",)),
                                                  ChemistryEdit("remove_atom", ("H2D",))])
    pdb = tmp_path / "1A6M_HEM.pdb"
    with open(os.path.join(DATA, "1A6M_HEM.pdb.txt")) as handle:
        pdb.write_text(handle.read())
    return ls.from_structure_residue(str(pdb), "A", 154, "", "HEM", hem).mol


def test_a_metal_in_the_ligand_needs_its_charge_and_never_rotates(tmp_path):
    """Meeko's own charge route protonates metal-bound atoms and fails on a heme; the metal-free
    part gets Gasteiger charges here and the metal its chosen formal charge."""
    heme = _heme_ligand(tmp_path)
    with pytest.raises(ValueError, match="formal charge for every metal in the ligand: FE"):
        prepare_ligand(heme, rigid_macrocycles=True)
    bonds = torsion_bonds(heme, rigid_macrocycles=True, metal_charges={"FE": 2})
    assert bonds and not any("FE" in b.atoms for b in bonds)
    prepared = prepare_ligand(heme, rigid_macrocycles=True, metal_charges={"FE": 2})
    iron = [l for l in prepared.pdbqt.splitlines() if l[12:16].strip() == "FE"]
    assert len(iron) == 1 and iron[0][70:76].strip() == "+2.000"
    assert abs(prepared.gasteiger_charge_sum - (-2)) < 0.01          # porphyrin -2, propionates -2, Fe +2
    assert any("chosen formal charge" in note for note in prepared.notes)
