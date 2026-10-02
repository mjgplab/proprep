"""Structure Loader option 4 loads a ligand file (SDF/mol2) for Molecular Docking.

Input files enter ProPrep through the Structure Loader; Molecular Docking
reads the workspace key ``ligand_file`` and does not browse for files
itself. The option is appended after the existing three so their keys, and
recorded sessions, are unchanged."""

import io
import os

import pytest
from rich.console import Console

import proprep.qmmm_prep.frame_extractor as fe
from proprep.structure_prep.structure_loader import StructureLoaderModule


class _Workspace(dict):
    def has(self, k):
        return k in self

    def set(self, k, v):
        self[k] = v


class _Processor:
    def __init__(self):
        self.output = io.StringIO()
        self.console = Console(file=self.output, width=200)
        self.workspace = _Workspace()

    def _get_workspace(self):
        return self.workspace


def _module():
    m = StructureLoaderModule()
    m.processor = _Processor()
    m.initialize()
    m.console = m.processor.console          # initialize() makes its own Console; record this one
    return m


def test_ligand_option_is_appended_after_the_existing_three():
    m = _module()
    assert list(m.get_menu_options()) == ["load", "md_setup", "metadata", "ligand_file"]
    keys = [o.key for o in m.get_enhanced_menu_options(m.processor.workspace)]
    assert keys == ["1", "2", "3", "4"]


def test_ligand_file_is_stored_and_its_molecules_listed(tmp_path, monkeypatch):
    Chem = pytest.importorskip("rdkit.Chem")
    from rdkit.Chem import AllChem
    path = tmp_path / "two.sdf"
    writer = Chem.SDWriter(str(path))
    for name, smiles in (("ethanol", "CCO"), ("acetate", "CC(=O)[O-]")):
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMolecule(mol, randomSeed=1)
        mol.SetProp("_Name", name)
        writer.write(mol)
    writer.close()
    monkeypatch.setattr(fe.FrameExtractor, "find_or_browse_file", lambda self, *a, **k: str(path))
    m = _module()
    assert m.handle_menu_option("ligand_file")
    assert m.processor.workspace["ligand_file"] == str(path.absolute())
    shown = m.processor.output.getvalue()
    assert "2 molecule(s)" in shown and "0: ethanol" in shown and "1: acetate" in shown
    options = {o.key: o for o in m.get_enhanced_menu_options(m.processor.workspace)}
    assert options["4"].status.name == "COMPLETED"
