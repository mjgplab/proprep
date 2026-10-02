#!/usr/bin/env python3
"""
Docking campaigns: a ligand library docked into one receptor and box.

A campaign cannot ask about each ligand, so the single-ligand questions
become recorded policies (unspecified stereo, 2D records, 3D build seed and
minimisation, torsions). These tests pin the policies, that a ligand which
cannot be read or prepared is recorded with its reason while the campaign
goes on, that an interrupted campaign resumes without redoing finished
ligands, that AD4 maps are built once for the whole library, and the menu.

Streptavidin (PDB 1STP, tests/data/docking) with biotin, desthiobiotin,
caffeine and a broken SMILES. A case this caught: Meeko rebuilds caffeine's
fused aromatic ring with other bond orders, so its poses are matched to the
ligand by element and connectivity when the exact match fails.

Run with: pytest tests/test_docking_campaign.py
"""

import csv
import io
import json
import os
import shutil

import pytest

Chem = pytest.importorskip("rdkit.Chem", reason="RDKit is required for docking campaigns")

from proprep.docking_prep import batch  # noqa: E402

DATA = os.path.join(os.path.dirname(__file__), "data", "docking")
BIOTIN = "[O-]C(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12"
LIBRARY = f"""# streptavidin test library
{BIOTIN} biotin
[O-]C(=O)CCCCC1SCC2NC(=O)NC12 biotin no stereo
Cn1cnc2c1c(=O)n(C)c(=O)n2C caffeine
C1CC broken
{BIOTIN} biotin
"""


def _policies(**changes):
    values = dict(unspecified_stereo="skip", flat_records="build", embed_seed=42, optimizer="mmff94",
                  rigid_conjugated=False, rigid_macrocycles=False)
    values.update(changes)
    return batch.BatchPolicies(**values)


@pytest.fixture
def smi(tmp_path):
    path = tmp_path / "library.smi"
    path.write_text(LIBRARY)
    return str(path)


# ---------------------------------------------------------------- library ---

def test_smiles_library_keeps_order_names_and_errors(smi):
    entries = batch.read_library(smi)
    assert [e.name for e in entries] == ["biotin", "biotin_no_stereo", "caffeine", "broken", "biotin_5"]
    assert [e.index for e in entries] == [1, 2, 3, 4, 5]
    assert entries[3].error and "C1CC" in entries[3].error
    assert all(e.error is None for i, e in enumerate(entries) if i != 3)


def test_2d_records_are_built_or_skipped_by_policy(tmp_path):
    from rdkit.Chem import rdDepictor
    path = str(tmp_path / "flat.sdf")
    writer = Chem.SDWriter(path)
    mol = Chem.MolFromSmiles(BIOTIN)
    rdDepictor.Compute2DCoords(mol)
    mol.SetProp("_Name", "biotin")
    writer.write(mol)
    writer.close()
    [entry] = batch.read_library(path)
    with pytest.raises(batch.LigandSkipped, match="2D record"):
        batch.build_entry(entry, _policies(flat_records="skip"))
    build = batch.build_entry(entry, _policies())
    assert build.mol.GetConformer().Is3D() and "2D record" in build.notes[0]
    Chem.AssignStereochemistryFrom3D(build.mol)
    assert Chem.MolToSmiles(Chem.RemoveHs(build.mol)) == Chem.MolToSmiles(Chem.MolFromSmiles(BIOTIN))


def test_unspecified_stereo_is_skipped_or_docked_by_policy(smi):
    entry = batch.read_library(smi)[1]
    with pytest.raises(batch.LigandSkipped, match="3 unspecified stereo"):
        batch.build_entry(entry, _policies())
    assert any(n.startswith("Stereocentre") for n in batch.build_entry(entry, _policies(unspecified_stereo="dock")).notes)


def test_policies_are_validated():
    with pytest.raises(ValueError):
        _policies(unspecified_stereo="maybe").validate()
    with pytest.raises(ValueError):
        _policies(embed_seed=0).validate()


# --------------------------------------------------------------- campaign ---

@pytest.fixture
def state(tmp_path):
    pytest.importorskip("meeko", reason="Meeko is required for docking")
    pytest.importorskip("vina", reason="AutoDock Vina is required for docking")
    from proprep.docking_prep.docking_run import DockingSettings
    from proprep.docking_prep.docking_state import DockingState
    pdb = tmp_path / "1STP.pdb"
    shutil.copy(os.path.join(DATA, "1STP.pdb.txt"), pdb)
    s = DockingState()
    s.receptor.pdb_path, s.receptor.chains, s.receptor.ph = str(pdb), ["A"], 7.0
    s.receptor.hetero_reviewed = True                     # menu 2 seen: crystal biotin left out
    s.receptor.protonation = {"A:87": "HID", "A:127": "HID"}
    s.receptor.termini = {"A:13": "chain break", "A:133": "chain break"}
    s.box = {"center": [11.419, 2.371, -11.373], "size": [22.545, 22.004, 23.519], "source": "crystal BTN A:300"}
    s.settings = DockingSettings("vina", 42, 8, 9, 3.0, 1.0, 0, 0.375, 0).to_dict()
    return s


def _csv(folder):
    with open(os.path.join(folder, "results.csv")) as handle:
        return {row["name"]: row for row in csv.DictReader(handle)}


def test_campaign_records_every_ligand_and_goes_on_past_failures(state, smi, tmp_path):
    result = batch.run_campaign(state, smi, _policies(), str(tmp_path / "docking"))
    rows = _csv(result.folder)
    assert rows["broken"]["status"] == "skipped" and "C1CC" in rows["broken"]["notes"]
    assert rows["biotin_no_stereo"]["status"] == "skipped"
    assert rows["caffeine"]["status"] == "docked"            # matched by topology after Meeko's rebuild
    assert rows["biotin"]["status"] == rows["biotin_5"]["status"] == "docked"
    assert rows["biotin"]["best_score"] == rows["biotin_5"]["best_score"]   # same seed, same result anywhere in the file
    assert rows["biotin"]["net_charge"] == "-1"
    top = [m.GetProp("_Name") for m in Chem.SDMolSupplier(os.path.join(result.folder, "top_poses.sdf"))]
    assert set(top) == {"biotin", "biotin_5", "caffeine"}
    assert top[-1] == "caffeine"                             # ranked best score first
    manifest = json.load(open(os.path.join(result.folder, "manifest.json")))
    assert manifest["finished"] is True and manifest["manifest"]["policies"]["unspecified_stereo"] == "skip"


def test_an_interrupted_campaign_resumes_without_redoing_finished_ligands(state, smi, tmp_path):
    base = str(tmp_path / "docking")
    seen = []

    def stop_after_first(text):
        seen.append(text)
        if text.startswith("[1/"):
            raise KeyboardInterrupt
    first = batch.run_campaign(state, smi, _policies(), base, progress=stop_after_first)
    assert first.interrupted
    from proprep.docking_prep.docking_pipeline import build_receptor
    manifest = batch._manifest(state, smi, _policies(), build_receptor(state))
    folder = batch.find_resumable(base, manifest)
    assert folder == first.folder
    assert batch.find_resumable(base, batch._manifest(state, smi, _policies(embed_seed=7), build_receptor(state))) is None
    log = []
    second = batch.run_campaign(state, smi, _policies(), base, resume_folder=folder, progress=log.append)
    assert not second.interrupted
    assert not any(line.startswith("[") and "biotin:" in line and "biotin_5" not in line for line in log)
    assert len(_csv(folder)) == 5 and batch.find_resumable(base, manifest) is None


def test_ad4_maps_are_built_once_for_the_whole_library(state, smi, tmp_path):
    from proprep.docking_prep.dependencies import find_executable
    if find_executable("autogrid4") is None:
        pytest.skip("autogrid4 is required for AD4 scoring")
    state.settings["scoring"] = "ad4"
    result = batch.run_campaign(state, smi, _policies(), str(tmp_path / "docking"))
    assert os.path.exists(os.path.join(result.folder, "ad4_maps", "receptor.gpf"))
    ligand_dirs = os.listdir(os.path.join(result.folder, "ligands"))
    assert ligand_dirs and not any(os.path.exists(os.path.join(result.folder, "ligands", d, "receptor.gpf"))
                                   for d in ligand_dirs)
    assert _csv(result.folder)["biotin"]["status"] == "docked"


def test_campaign_from_the_menu(state, smi, tmp_path, monkeypatch):
    from tests.test_docking_module import _Processor, _scripted
    from proprep.docking_prep.docking_module import MolecularDockingModule
    module = MolecularDockingModule()
    module.processor = _Processor(tmp_path)
    module.initialize()
    module.processor.workspace["ligand_file"] = smi
    module.processor.workspace["docking_state"] = state.to_dict()
    import proprep.utils.structure_selector as selector
    monkeypatch.setattr(selector.StructureSelector, "get_structure", lambda self, *a, **k: state.receptor.pdb_path)
    answers = ["c", "s", "b", "42", "mmff94", "n", "y", "y",       # policies, then start
               "v", "1",                                          # campaign results (no single run, so not asked which), rank 1
               "b"]
    _scripted(monkeypatch, answers, [])
    assert module.process(module.processor.workspace) is True
    assert answers == []
    shown = module.processor.out.getvalue()
    assert "3 docked, 2 skipped, 0 failed." in shown and "Top 3 by score" in shown
    assert module.processor.workspace["docking_campaign_dir"].startswith(str(tmp_path))
