"""The viewer tells NGL each file's format. AlphaFold DB downloads are mmCIF;
with ``ext`` hard-wired to "pdb" they parsed as an empty PDB and never
appeared in the viewer."""

import re
from pathlib import Path

from proprep.structure_prep.interactive_structure_viewer import (
    InteractiveStructureViewer,
    structure_display_name,
    structure_file_ext,
)


def test_ext_from_suffix():
    assert structure_file_ext("/x/AF-P00520-F1-model_v4.cif") == "cif"
    assert structure_file_ext("/x/1abc.pdb") == "pdb"
    assert structure_file_ext("/x/1abc.ent") == "pdb"
    assert structure_file_ext("/x/1abc.pdb.gz") == "pdb"
    assert structure_file_ext("/x/lig.mol2") == "mol2"
    assert structure_file_ext("/x/no_suffix") == "pdb"


def test_display_name_strips_any_structure_suffix():
    assert structure_display_name("/x/AF-P00520-F1-model_v4.cif") == "AF-P00520-F1-model_v4"
    assert structure_display_name("/x/1abc.pdb") == "1abc"
    assert structure_display_name("/x/1abc.pdb.gz") == "1abc"
    assert structure_display_name("/x/notes.txt") == "notes.txt"


def test_viewer_config_carries_ext_per_structure():
    v = InteractiveStructureViewer()
    v.selected_structures = ["/tmp/1abc.pdb", "/tmp/AF-P00520-F1-model_v4.cif"]
    v.viewer_config = {}
    structures = v._build_viewer_config()["structures"]
    assert [s["ext"] for s in structures] == ["pdb", "cif"]
    assert [s["name"] for s in structures] == ["1abc", "AF-P00520-F1-model_v4"]


def test_scene_config_carries_ext():
    v = InteractiveStructureViewer()
    v.selected_structures = ["/tmp/AF-P00520-F1-model_v4.cif"]
    v.viewer_config = {}
    v.scene_override = {
        "_for": list(v.selected_structures),
        "representations": {0: [{"id": "x", "type": "cartoon", "sele": "protein"}]},
    }
    structures = v._build_viewer_config()["structures"]
    assert structures[0]["ext"] == "cif"


def test_template_passes_ext_to_ngl():
    html = Path("src/proprep/structure_prep/templates/ngl_viewer.html").read_text()
    assert re.search(r'ext:\s*structInfo\.ext\s*\|\|\s*"pdb"', html)
    assert 'ext: "pdb", firstModelOnly' not in html
