#!/usr/bin/env python3
"""
Docking runs end to end with the real Vina (and autogrid4 for AD4 scoring).

The case is biotin redocked into streptavidin (PDB 1STP, unmodified in
tests/data/docking), with the receptor built from pdb2pqr's protonation
proposals when pdb2pqr is available and biotin as its carboxylate. The
reference is the crystal pose; the top pose must lie within 2 A of it.

Also pinned: a seed repeats exactly, seed 0 (Vina's "random") is refused,
Vina's console warnings end up in the run's notes rather than on the
terminal, and poses come back with the ligand's own atom names.

Run with: pytest tests/test_docking_run.py
"""

import os
import shutil

import numpy as np
import pytest

pytest.importorskip("rdkit", reason="RDKit is required for docking")
pytest.importorskip("meeko", reason="Meeko is required for docking")
pytest.importorskip("vina", reason="AutoDock Vina is required for docking")

from rdkit import Chem  # noqa: E402

from proprep.docking_prep import ligand_sources as ls  # noqa: E402
from proprep.docking_prep.ccd_chemistry import component_from_block  # noqa: E402
from proprep.docking_prep.chemistry_edits import ChemistryEdit, apply_edits  # noqa: E402
from proprep.docking_prep.dependencies import find_executable  # noqa: E402
from proprep.docking_prep.docking_results import heavy_atom_rmsd, poses_from_run, write_poses_sdf  # noqa: E402
from proprep.docking_prep.docking_run import DockingSettings, box_around, run_docking  # noqa: E402
from proprep.docking_prep.ligand_prep import prepare_ligand  # noqa: E402
from proprep.docking_prep.receptor_prep import prepare_receptor  # noqa: E402

DATA = os.path.join(os.path.dirname(__file__), "data", "docking")


def _settings(**changes):
    values = dict(scoring="vina", seed=42, exhaustiveness=8, n_poses=9, energy_range=3.0,
                  min_rmsd=1.0, max_evals=0, spacing=0.375, cpu=0)
    values.update(changes)
    return DockingSettings(**values)


@pytest.fixture(scope="module")
def system(tmp_path_factory):
    folder = tmp_path_factory.mktemp("stp")
    pdb = str(folder / "1STP.pdb")
    shutil.copy(os.path.join(DATA, "1STP.pdb.txt"), pdb)
    with open(os.path.join(DATA, "ccd_BTN.cif.txt")) as handle:
        btn = component_from_block(handle.read(), "BTN")
    # explicit ligand choice: carboxylate (remove the acid H from O12, charge -1)
    acid_h = [a for a in btn.mol.GetAtoms() if a.GetSymbol() == "H" and a.GetNeighbors()[0].GetSymbol() == "O"]
    assert len(acid_h) == 1
    oxygen = acid_h[0].GetNeighbors()[0].GetProp("name")
    btn.mol, btn.problems = apply_edits(btn.mol, [ChemistryEdit("remove_atom", (acid_h[0].GetProp("name"),)),
                                                  ChemistryEdit("set_formal_charge", (oxygen,), -1)])
    ligand = ls.from_structure_residue(pdb, "A", 300, "", "BTN", btn)
    prepared = prepare_ligand(ligand.mol, rigid_macrocycles=False)
    receptor = prepare_receptor(pdb, ["A"], [], protonation={"A:87": "HID", "A:127": "HID"},
                                termini={"A:13": "chain break", "A:133": "chain break"},
                                metal_charges={}, altlocs={}, cofactor_edits={})
    flexible = prepare_receptor(pdb, ["A"], [], protonation={"A:87": "HID", "A:127": "HID"},
                                termini={"A:13": "chain break", "A:133": "chain break"},
                                metal_charges={}, altlocs={}, cofactor_edits={}, flexible=["A:45", "A:88"])
    box = box_around(Chem.RemoveHs(ligand.mol).GetConformer().GetPositions(), 8.0, "crystal BTN A:300")
    return dict(ligand=ligand, prepared=prepared, receptor=receptor, flexible=flexible, box=box, folder=folder)


def _dock(system, settings, receptor="receptor"):
    rec = system[receptor]
    run = run_docking(rec.rigid_pdbqt, rec.flex_pdbqt, system["prepared"].pdbqt, system["box"], settings)
    poses = poses_from_run(run.poses_pdbqt, system["ligand"].mol, run.energies)
    return run, poses


# ------------------------------------------------------------------- units ---

def test_box_encloses_the_atoms_with_padding():
    box = box_around(np.array([[0.0, 0.0, 0.0], [2.0, 4.0, 6.0]]), 5.0, "test")
    assert box.center == (1.0, 2.0, 3.0) and box.size == (12.0, 14.0, 16.0)
    assert "padding 5 A" in box.source
    with pytest.raises(ValueError):
        box_around(np.zeros((0, 3)), 5.0, "empty")


def test_seed_zero_is_refused_because_vina_reads_it_as_random():
    with pytest.raises(ValueError, match="random"):
        _settings(seed=0).validate()
    with pytest.raises(ValueError, match="Scoring"):
        _settings(scoring="glide").validate()


# ---------------------------------------------------------------- docking ---

def test_vina_redocks_biotin_and_repeats_with_the_same_seed(system):
    first_run, first = _dock(system, _settings())
    second_run, second = _dock(system, _settings())
    top = first[0]
    assert heavy_atom_rmsd(top.ligand, system["ligand"].mol) < 2.0
    assert [p.energies for p in first] == [p.energies for p in second]
    assert np.allclose(top.ligand.GetConformer().GetPositions(), second[0].ligand.GetConformer().GetPositions())
    names = {a.GetProp("name") for a in top.ligand.GetAtoms()}
    assert names == {a.GetProp("name") for a in system["ligand"].mol.GetAtoms()}
    assert set(top.energies) == {"total", "inter", "intra", "torsions", "intra best pose"}


def test_vina_console_warnings_are_captured_into_the_notes(system, capfd):
    run, _ = _dock(system, _settings(exhaustiveness=1))
    captured = capfd.readouterr()
    assert "WARNING" not in captured.out and "WARNING" not in captured.err
    assert all(isinstance(note, str) for note in run.notes)


def test_flexible_side_chains_come_back_with_every_pose(system):
    run, poses = _dock(system, _settings(), receptor="flexible")
    assert len(poses[0].flexible_residues) == 2
    assert heavy_atom_rmsd(poses[0].ligand, system["ligand"].mol) < 2.0


@pytest.mark.skipif(find_executable("autogrid4") is None, reason="autogrid4 is required for AD4 scoring")
def test_ad4_scoring_through_autogrid_maps(system):
    run, poses = _dock(system, _settings(scoring="ad4"))
    assert heavy_atom_rmsd(poses[0].ligand, system["ligand"].mol) < 2.0
    assert set(poses[0].energies) == {"total", "inter", "intra", "torsions", "-intra"}
    assert any("AutoGrid maps" in note and "distance-dependent" in note for note in run.notes)


def test_poses_sdf_carries_rank_energies_and_run_properties(system, tmp_path):
    _, poses = _dock(system, _settings())
    for pose in poses:
        pose.rmsd_to_reference = heavy_atom_rmsd(pose.ligand, system["ligand"].mol)
    path = str(tmp_path / "poses.sdf")
    write_poses_sdf(path, poses, {"scoring": "vina", "seed": 42})
    records = [m for m in Chem.SDMolSupplier(path, removeHs=False)]
    assert len(records) == len(poses)
    first = records[0]
    assert first.GetProp("rank") == "1" and first.GetProp("seed") == "42"
    assert first.GetProp("energy total (kcal/mol)") == f"{poses[0].energies['total']:.3f}"
