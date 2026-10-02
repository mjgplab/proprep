#!/usr/bin/env python3
"""
Docking ligand intake: CCD chemistry checks, named edits, and the three ligand routes.

Fixtures in tests/data/docking are unmodified real data: the BTN and HEM
entries of the Chemical Component Dictionary, the biotin HETATM records of
PDB 1STP and the heme HETATM records of PDB 1A6M.

The CCD's heme entries give O2A and O2D a charge of -1 while also bonding a
hydrogen to each; these tests pin that the problem is reported, that nothing
resolves it silently, and that either explicit resolution (carboxylate or
carboxylic acid) yields valid chemistry.

Run with: pytest tests/test_docking_ligand_intake.py
"""

import os

import numpy as np
import pytest

Chem = pytest.importorskip("rdkit.Chem", reason="RDKit is required for docking ligand intake")
pytest.importorskip("gemmi", reason="gemmi is required to read CCD entries")

from proprep.docking_prep import ligand_sources as ls  # noqa: E402
from proprep.docking_prep.ccd_chemistry import component_from_block  # noqa: E402
from proprep.docking_prep.chemistry_edits import ChemistryEdit, apply_edits  # noqa: E402

DATA = os.path.join(os.path.dirname(__file__), "data", "docking")


def _read(name):
    with open(os.path.join(DATA, name)) as handle:
        return handle.read()


@pytest.fixture
def btn():
    return component_from_block(_read("ccd_BTN.cif.txt"), "BTN")


@pytest.fixture
def hem():
    return component_from_block(_read("ccd_HEM.cif.txt"), "HEM")


def _names(mol):
    return {atom.GetProp("name") for atom in mol.GetAtoms()}


def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


# ------------------------------------------------------------- CCD checks ---

def test_clean_ccd_entry_has_no_problems(btn):
    assert btn.problems == []
    assert btn.mol.GetNumAtoms() == 32
    assert btn.coordinates == "ideal"
    assert "BIOTIN" in btn.name.upper()


def test_heme_entry_reports_only_the_protonated_carboxylate_oxygens(hem):
    assert sorted(p.atom_name for p in hem.problems) == ["O2A", "O2D"]
    for problem in hem.problems:
        assert problem.formal_charge == -1
        assert any(name.startswith("H2") for name in problem.bonded_atoms)


def test_metal_coordination_bonds_are_dative_not_counted_against_nitrogen(hem):
    dative = [b for b in hem.mol.GetBonds() if b.GetBondType() == Chem.BondType.DATIVE]
    assert len(dative) == 4
    assert all(b.GetEndAtom().GetSymbol() == "Fe" for b in dative)
    assert not any(p.atom_name in ("NA", "NB", "NC", "ND") for p in hem.problems)


# ----------------------------------------------------------------- edits ---

def test_either_explicit_resolution_of_the_heme_gives_valid_chemistry(hem):
    carboxylate, left = apply_edits(hem.mol, [ChemistryEdit("remove_atom", ("H2A",)),
                                              ChemistryEdit("remove_atom", ("H2D",))])
    assert left == [] and Chem.GetFormalCharge(carboxylate) == -4
    acid, left = apply_edits(hem.mol, [ChemistryEdit("set_formal_charge", ("O2A",), 0),
                                       ChemistryEdit("set_formal_charge", ("O2D",), 0)])
    assert left == [] and Chem.GetFormalCharge(acid) == -2


def test_partial_edits_leave_the_remaining_problem_visible(hem):
    _, left = apply_edits(hem.mol, [ChemistryEdit("remove_atom", ("H2A",))])
    assert [p.atom_name for p in left] == ["O2D"]


def test_added_hydrogen_is_named_and_placed_at_a_bonding_distance(hem):
    carboxylate, _ = apply_edits(hem.mol, [ChemistryEdit("remove_atom", ("H2A",)),
                                           ChemistryEdit("remove_atom", ("H2D",))])
    acid, left = apply_edits(carboxylate, [ChemistryEdit("set_formal_charge", ("O2A",), 0),
                                           ChemistryEdit("add_hydrogen", ("O2A",))])
    assert left == []
    index = {a.GetProp("name"): a.GetIdx() for a in acid.GetAtoms()}
    hydrogens = [n for n in acid.GetAtomWithIdx(index["O2A"]).GetNeighbors() if n.GetSymbol() == "H"]
    assert len(hydrogens) == 1
    positions = acid.GetConformer().GetPositions()
    distance = np.linalg.norm(positions[index["O2A"]] - positions[hydrogens[0].GetIdx()])
    assert 0.9 < distance < 1.1
    assert len(_names(acid)) == acid.GetNumAtoms()


@pytest.mark.parametrize("edit, error", [
    (ChemistryEdit("set_bond_order", ("C1A", "C1B"), 2), ValueError),     # not bonded
    (ChemistryEdit("remove_atom", ("NOPE",)), KeyError),
])
def test_bad_edits_fail_with_the_atom_names(hem, edit, error):
    with pytest.raises(error, match=edit.atoms[-1]):
        apply_edits(hem.mol, [edit])


def test_edit_validation_and_round_trip():
    with pytest.raises(ValueError):
        ChemistryEdit("set_bond_order", ("C1",), 2)
    with pytest.raises(ValueError):
        ChemistryEdit("set_bond_order", ("C1", "C2"), 4)
    with pytest.raises(ValueError):
        ChemistryEdit("rename", ("C1",))
    edit = ChemistryEdit("set_formal_charge", ("N1",), 1)
    assert ChemistryEdit.from_dict(edit.to_dict()) == edit


# ------------------------------------------------------ structure residue ---

def test_structure_residue_uses_crystal_heavy_atoms_and_ccd_chemistry(tmp_path, btn):
    path = _write(tmp_path, "1STP_BTN.pdb", _read("1STP_BTN.pdb.txt"))
    build = ls.from_structure_residue(path, "A", 300, "", "BTN", btn)
    mol = build.mol
    assert _names(mol) == _names(btn.mol)
    assert build.net_charge == 0
    crystal = {line[12:16].strip(): [float(line[30:38]), float(line[38:46]), float(line[46:54])]
               for line in _read("1STP_BTN.pdb.txt").splitlines() if line.startswith("HETATM")}
    positions = mol.GetConformer().GetPositions()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() > 1:
            assert np.allclose(positions[atom.GetIdx()], crystal[atom.GetProp("name")], atol=1e-3)
    assert build.source["kind"] == "structure_residue" and build.source["ccd_code"] == "BTN"


def test_unresolved_ccd_problem_is_refused_and_resolved_entry_builds(tmp_path, hem):
    path = _write(tmp_path, "1A6M_HEM.pdb", _read("1A6M_HEM.pdb.txt"))
    with pytest.raises(ValueError, match="O2A"):
        ls.from_structure_residue(path, "A", 154, "", "HEM", hem)
    hem.mol, hem.problems = apply_edits(hem.mol, [ChemistryEdit("remove_atom", ("H2A",)),
                                                  ChemistryEdit("remove_atom", ("H2D",))])
    build = ls.from_structure_residue(path, "A", 154, "", "HEM", hem)
    assert build.mol.GetNumAtoms() == 73
    assert build.net_charge == -4


def test_alternate_locations_need_an_explicit_choice(tmp_path, btn):
    lines = [l for l in _read("1STP_BTN.pdb.txt").splitlines() if l.startswith("HETATM")]
    doubled = []
    for line in lines:
        doubled.append(line[:16] + "A" + line[17:])
        shifted = f"{float(line[30:38]) + 1.0:8.3f}"
        doubled.append(line[:16] + "B" + line[17:30] + shifted + line[38:])
    path = _write(tmp_path, "altloc.pdb", "\n".join(doubled) + "\nEND\n")
    with pytest.raises(ValueError, match="alternate locations A, B"):
        ls.from_structure_residue(path, "A", 300, "", "BTN", btn)
    a = ls.from_structure_residue(path, "A", 300, "", "BTN", btn, altloc="A")
    b = ls.from_structure_residue(path, "A", 300, "", "BTN", btn, altloc="B")
    shift = b.mol.GetConformer().GetPositions() - a.mol.GetConformer().GetPositions()
    heavy = [atom.GetIdx() for atom in a.mol.GetAtoms() if atom.GetAtomicNum() > 1]
    assert np.allclose(shift[heavy, 0], 1.0, atol=1e-3)
    assert "alternate location B" in b.notes[0]


def test_missing_heavy_atoms_are_listed_not_guessed(tmp_path, btn):
    lines = [l for l in _read("1STP_BTN.pdb.txt").splitlines()
             if l.startswith("HETATM") and l[12:16].strip() != "O11"]
    path = _write(tmp_path, "missing.pdb", "\n".join(lines) + "\nEND\n")
    with pytest.raises(ValueError, match="missing from the structure: O11"):
        ls.from_structure_residue(path, "A", 300, "", "BTN", btn)


# ------------------------------------------------------------------ SMILES ---

BIOTIN_STEREO = "[O-]C(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12"
BIOTIN_FLAT = "[O-]C(=O)CCCCC1SCC2NC(=O)NC12"


def test_smiles_build_is_reproducible_with_the_same_seed():
    first = ls.from_smiles(BIOTIN_STEREO, embed_seed=11, optimizer="none")
    second = ls.from_smiles(BIOTIN_STEREO, embed_seed=11, optimizer="none")
    assert np.allclose(first.mol.GetConformer().GetPositions(), second.mol.GetConformer().GetPositions())
    assert first.net_charge == -1
    assert first.source == {"kind": "smiles", "smiles": BIOTIN_STEREO, "embed_seed": 11,
                            "optimizer": "none", "residue_name": "LIG"}


def test_unspecified_stereocentres_are_reported_with_what_was_built():
    build = ls.from_smiles(BIOTIN_FLAT, embed_seed=7, optimizer="mmff94")
    stereo_notes = [n for n in build.notes if n.startswith("Stereocentre")]
    assert len(stereo_notes) == 3
    specified = ls.from_smiles(BIOTIN_STEREO, embed_seed=7, optimizer="mmff94")
    assert not any(n.startswith("Stereocentre") for n in specified.notes)


def test_smiles_rejects_unknown_optimizer_and_bad_smiles():
    with pytest.raises(ValueError, match="optimizer"):
        ls.from_smiles("CCO", embed_seed=1, optimizer="uff")
    with pytest.raises(ValueError, match="could not read"):
        ls.from_smiles("C1CC", embed_seed=1, optimizer="none")


# ------------------------------------------------------------------- files ---

def test_file_with_several_records_needs_a_choice(tmp_path):
    path = str(tmp_path / "two.sdf")
    writer = Chem.SDWriter(path)
    for name, smiles in (("biotin", BIOTIN_STEREO), ("ethanol", "CCO")):
        mol = ls.from_smiles(smiles, embed_seed=3, optimizer="none").mol
        mol.SetProp("_Name", name)
        writer.write(mol)
    writer.close()
    with pytest.raises(ValueError, match="0: biotin; 1: ethanol"):
        ls.from_file(path)
    assert ls.from_file(path, record=1).mol.GetNumAtoms() == 9


def test_file_without_hydrogens_reports_how_many_were_added(tmp_path, btn):
    pdb = _write(tmp_path, "1STP_BTN.pdb", _read("1STP_BTN.pdb.txt"))
    heavy = Chem.RemoveHs(ls.from_structure_residue(pdb, "A", 300, "", "BTN", btn).mol)
    path = str(tmp_path / "noH.mol")
    Chem.MolToMolFile(heavy, path)
    build = ls.from_file(path)
    assert any("16 were added from valence" in note for note in build.notes)
    assert any("Atom names generated" in note for note in build.notes)


def test_flat_file_is_refused(tmp_path):
    mol = Chem.MolFromSmiles("CCO")
    Chem.rdDepictor.Compute2DCoords(mol)
    path = str(tmp_path / "flat.mol")
    Chem.MolToMolFile(mol, path)
    with pytest.raises(ValueError, match="no 3D coordinates"):
        ls.from_file(path)


def test_atom_names_are_unique_after_assignment():
    mol = Chem.AddHs(Chem.MolFromSmiles("OCCO"))
    mol.GetAtomWithIdx(0).SetProp("name", "O1")
    mol.GetAtomWithIdx(3).SetProp("name", "O1")          # duplicate given name
    ls.assign_atom_names(mol)
    assert len(_names(mol)) == mol.GetNumAtoms()
