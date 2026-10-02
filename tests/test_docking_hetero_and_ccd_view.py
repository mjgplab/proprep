"""
Molecular Docking: what menu 2 keeps, and the CCD entry shown while it is reviewed.

Line 2 of the dashboard used to read "Kept HETATM: 0 residue(s)" both before
menu 2 was opened and after everything was deliberately left out, and a run
could start without menu 2 ever being seen, silently dropping every cofactor
and metal from the receptor. It now says what is kept, what is left out and
which residue is the ligand, and a run needs menu 2 whenever the file has
HETATM residues besides the ligand and waters.

The CCD review (choosing a ligand from the structure, or a cofactor in menu 4)
showed nothing in the viewer. It now shows the entry on its residue: heavy
atoms on the crystal atoms (matched by name), hydrogens placed by local fits
(a rigid fit of the ideal conformer misses flexible ligands by Angstroms),
and accepts 'p' to pick an atom or bond there.

Streptavidin with biotin (PDB 1STP, tests/data/docking).

Run with: pytest tests/test_docking_hetero_and_ccd_view.py
"""

import os
import shutil
from unittest.mock import MagicMock

import numpy as np
import pytest

Chem = pytest.importorskip("rdkit.Chem", reason="RDKit is required for docking")
pytest.importorskip("gemmi", reason="gemmi is required to read CCD entries")

from tests.test_docking_module import _Processor, _scripted  # noqa: E402

DATA = os.path.join(os.path.dirname(__file__), "data", "docking")


@pytest.fixture
def module(tmp_path):
    from proprep.docking_prep.docking_module import MolecularDockingModule
    pdb = tmp_path / "1STP.pdb"
    shutil.copy(os.path.join(DATA, "1STP.pdb.txt"), pdb)
    m = MolecularDockingModule()
    m.processor = _Processor(tmp_path)
    m.initialize()
    m.state.receptor.pdb_path, m.state.receptor.chains = str(pdb), ["A"]
    return m


@pytest.fixture
def btn():
    from proprep.docking_prep.ccd_chemistry import component_from_block
    with open(os.path.join(DATA, "ccd_BTN.cif.txt")) as handle:
        return component_from_block(handle.read(), "BTN")


# --------------------------------------------------------------- menu 2 ---

def test_dashboard_line_2_says_what_is_kept_and_a_run_needs_menu_2(module):
    s = module.state
    assert module._hetero_summary() == "not reviewed; 1 in the file: BTN A:300"
    assert "which HETATM residues stay (3)" in s.missing()

    s.receptor.hetero_reviewed = True                       # menu 2 seen, biotin left out
    assert module._hetero_summary() == "0 of 1 kept; left out: BTN A:300"
    assert "which HETATM residues stay (3)" not in s.missing()

    s.receptor.keep_hetero = ["A:300", "A:442", "A:443"]
    assert module._hetero_summary() == "1 of 1 kept: BTN A:300; 2 of 84 waters kept"


def test_a_ligand_from_the_structure_is_named_on_line_2_and_needs_no_review(module):
    s = module.state
    s.ligand.source = {"kind": "structure_residue", "chain": "A", "resseq": 300, "icode": "", "resname": "BTN"}
    assert module._hetero_summary() == "no other HETATM residues in the file; BTN A:300 is the ligand (2)"
    assert "which HETATM residues stay (3)" not in s.missing()


def test_the_review_is_saved_with_the_state(module):
    from proprep.docking_prep.docking_state import DockingState
    module.state.receptor.hetero_reviewed = True
    assert DockingState.from_dict(module.state.to_dict()).receptor.hetero_reviewed is True
    assert DockingState.from_dict({"receptor": {"pdb_path": "x"}}).receptor.hetero_reviewed is False


# ------------------------------------------------------ CCD in the viewer ---

def test_the_ccd_entry_is_drawn_on_its_crystal_residue(module, btn, monkeypatch):
    import proprep.structure_prep.viewer_coordinator as vc
    from proprep.docking_prep.receptor_prep import read_atoms
    fake = MagicMock()
    monkeypatch.setattr(vc, "viewer", fake)

    assert module.show_ccd_in_viewer(btn.mol, "BTN", "A:300") is True

    receptor, preview = fake.show_structures.call_args.args[0]
    assert receptor == module.state.receptor.pdb_path
    fake.focus_on.assert_called_with(":A and 300")
    shown = Chem.MolFromMolFile(preview, removeHs=False, sanitize=False)
    assert shown.GetNumAtoms() == btn.mol.GetNumAtoms()           # same atoms, same order: picks map back
    positions = shown.GetConformer().GetPositions()
    crystal = {a.name: np.array(a.xyz) for a in read_atoms(receptor) if a.res_id == "A:300"}
    names = [a.GetProp("name") for a in btn.mol.GetAtoms()]
    heavy = [np.linalg.norm(positions[i] - crystal[n]) for i, n in enumerate(names) if n in crystal]
    assert len(heavy) == 16 and max(heavy) < 1e-3
    to_h = [np.linalg.norm(positions[b.GetBeginAtomIdx()] - positions[b.GetEndAtomIdx()])
            for b in btn.mol.GetBonds() if 1 in (b.GetBeginAtom().GetAtomicNum(), b.GetEndAtom().GetAtomicNum())]
    assert to_h and 0.9 < min(to_h) and max(to_h) < 1.2


def test_without_matching_atoms_the_residue_is_highlighted_instead(module, btn, monkeypatch):
    import proprep.structure_prep.viewer_coordinator as vc
    fake = MagicMock()
    monkeypatch.setattr(vc, "viewer", fake)
    assert module.show_ccd_in_viewer(btn.mol, "BTN", "A:449") is False   # a water: no BTN atoms there
    assert fake.show_structures.call_args.args[0] == [module.state.receptor.pdb_path]
    assert fake.highlight.call_args.args[0] == ":A and 449"


def test_a_bond_order_is_picked_in_the_viewer_during_the_ccd_review(module, btn, monkeypatch):
    import proprep.docking_prep.docking_menus_receptor as menus
    from proprep.docking_prep import docking_ui
    import proprep.structure_prep.viewer_coordinator as vc
    monkeypatch.setattr(vc, "viewer", MagicMock())
    monkeypatch.setattr(menus, "fetch_component", lambda code: btn)
    index = {a.GetProp("name"): a.GetIdx() for a in btn.mol.GetAtoms()}
    monkeypatch.setattr(docking_ui, "pick", lambda m, kind, prompt, structure_index: {
        "atom1": {"index": index["C2"]}, "atom2": {"index": index["S1"]}})
    asked = []
    answers = ["p", "b", "1", ""]                           # pick a bond, keep it single, then done
    _scripted(monkeypatch, answers, asked)
    store = {}
    menus.review_ccd_chemistry(module, "BTN", store, "A:300")
    assert answers == []
    assert asked[0].startswith("Edit for BTN (") and "'p' to pick in the viewer" in asked[0]
    assert asked[2] == "Bond order for C2-S1"
    assert store == {"BTN": [{"op": "set_bond_order", "atoms": ["C2", "S1"], "value": 1}]}


# ------------------------------------------- one edit prompt for a ligand ---

def test_an_inconsistent_ccd_entry_is_reviewed_before_the_ligand_is_built(tmp_path, monkeypatch):
    """HEM's CCD entry protonates O2A and O2D yet charges them -1: it cannot be built until
    fixed, so its review prompt comes first. A sound entry (BTN, in test_docking_module)
    goes straight to the one ligand edit prompt."""
    import proprep.docking_prep.docking_menus_receptor as menus
    import proprep.structure_prep.viewer_coordinator as vc
    from proprep.docking_prep.ccd_chemistry import component_from_block
    from proprep.docking_prep.docking_menus_ligand import choose_ligand
    from proprep.docking_prep.docking_module import MolecularDockingModule
    with open(os.path.join(DATA, "ccd_HEM.cif.txt")) as handle:
        block = handle.read()
    monkeypatch.setattr(menus, "fetch_component", lambda code: component_from_block(block, "HEM"))
    monkeypatch.setattr(vc, "viewer", MagicMock())
    pdb = tmp_path / "hem.pdb"
    shutil.copy(os.path.join(DATA, "1A6M_HEM.pdb.txt"), pdb)
    m = MolecularDockingModule()
    m.processor = _Processor(tmp_path)
    m.initialize()
    m.state.receptor.pdb_path = str(pdb)
    answers = ["s", "1", "remove H2A", "remove H2D", "",     # the CCD entry, fixed
               "",                                           # the built ligand: no further edits
               "2"]                                          # iron's formal charge
    asked = []
    _scripted(monkeypatch, answers, asked)
    choose_ligand(m)
    assert answers == []
    assert [p.split(" (")[0] for p in asked if p.startswith("Edit for")] == ["Edit for HEM"] * 3 + ["Edit for the ligand"]
    assert m.state.ligand.ccd_edits == [{"op": "remove_atom", "atoms": ["H2A"]},
                                        {"op": "remove_atom", "atoms": ["H2D"]}]


# ---------------------------------------------------- viewer reloads ---

def test_the_viewer_reloads_only_when_what_it_shows_has_changed(module, tmp_path, monkeypatch):
    import proprep.structure_prep.viewer_coordinator as vc
    fake = MagicMock()
    shown = []
    fake.current_structures.side_effect = lambda: list(shown)
    fake.show_structures.side_effect = lambda paths, **kwargs: shown.__setitem__(slice(None), paths)
    monkeypatch.setattr(vc, "viewer", fake)
    preview = tmp_path / "ligand_preview.sdf"
    preview.write_text("first\n")
    files = [module.state.receptor.pdb_path, str(preview)]

    module._show_files(files)                 # new files: shown
    module._show_files(files)                 # same files, same contents: nothing
    assert fake.show_structures.call_count == 1 and fake.refresh_structure.call_count == 0
    preview.write_text("after an edit\n")
    module._show_files(files)                 # same files, new contents: re-read
    assert fake.show_structures.call_count == 1 and fake.refresh_structure.call_count == 1


# ------------------------------------------- crystal copy under the preview ---

def test_the_crystal_ligand_is_set_aside_under_its_preview(module, monkeypatch):
    """The crystal ligand in the receptor file has no bond orders and lies exactly under
    the preview: it leaves structure 0's default representations for its own, hidden one.
    Cofactors kept in the receptor stay in the default Ligands representation."""
    import proprep.structure_prep.viewer_coordinator as vc
    fake = MagicMock()
    fake.current_structures.return_value = []
    monkeypatch.setattr(vc, "viewer", fake)
    module.state.ligand.source = {"kind": "structure_residue", "chain": "A", "resseq": 300, "icode": "",
                                  "resname": "BTN"}
    build = MagicMock()
    build.mol = Chem.AddHs(Chem.MolFromSmiles("CCO"), addCoords=False)
    from rdkit.Chem import AllChem
    AllChem.EmbedMolecule(build.mol, randomSeed=1)
    module.show_ligand_in_viewer(build)
    assert fake.show_structures.call_args.kwargs["set_aside"] == {
        0: {"selection": ":A and 300", "label": "Crystal ligand (A:300)"}}

    module.state.ligand.source = {"kind": "smiles", "smiles": "CCO"}       # nothing in the receptor to hide
    module.show_ligand_in_viewer(build)
    assert fake.show_structures.call_args.kwargs["set_aside"] is None


def test_the_viewer_draws_a_set_aside_residue_separately_and_hidden(tmp_path):
    from proprep.structure_prep.interactive_structure_viewer import InteractiveStructureViewer
    v = InteractiveStructureViewer.__new__(InteractiveStructureViewer)
    v.selected_structures = [str(tmp_path / "1SDU.pdb"), str(tmp_path / "ligand_preview.sdf")]
    v.available_annotations, v.annotation_config, v.shape_config = {}, {}, {}
    v.trajectory_files, v.density_by_file, v.scene_override, v.processor = {}, {}, None, None
    v._scene_dir = lambda: str(tmp_path)
    v.viewer_config = {"set_aside": {0: {"selection": ":B and 902", "label": "Crystal ligand (B:902)"}}}
    structures = v._build_viewer_config()["structures"]
    receptor = {r["id"]: r for r in structures[0]["representations"]}
    assert receptor["default_ligands"]["selection"] == "(hetero and not (water or ion)) and not (:B and 902)"
    assert receptor["default_ligands"]["visible"] is True                  # a kept heme still shows
    assert receptor["set_aside"] == {"id": "set_aside", "label": "Crystal ligand (B:902)",
                                     "selection": ":B and 902", "style": "ball+stick", "color": "element",
                                     "visible": False}
    preview = {r["id"]: r for r in structures[1]["representations"]}
    assert "set_aside" not in preview and preview["default_ligands"]["multiple_bond"] is True


# --------------------------------------------------- taking back an edit ---

def test_drop_takes_back_the_last_edit_and_never_collides_with_the_session_rewind(module, btn, monkeypatch):
    """'undo' at any prompt rewinds the whole session (utils.session_rewind) before the
    module sees it, so the edit prompts' own 'undo' could never run. They use 'drop'."""
    import proprep.docking_prep.docking_menus_receptor as menus
    from proprep.docking_prep import docking_ui
    from proprep.utils.session_rewind import is_rewind_keyword
    import proprep.structure_prep.viewer_coordinator as vc
    monkeypatch.setattr(vc, "viewer", MagicMock())
    monkeypatch.setattr(menus, "fetch_component", lambda code: btn)
    assert not is_rewind_keyword(docking_ui.DROP)
    answers = ["remove HO2", "charge O12 -1", "drop", ""]
    asked = []
    _scripted(monkeypatch, answers, asked)
    store = {}
    menus.review_ccd_chemistry(module, "BTN", store, "A:300")
    assert store == {"BTN": [{"op": "remove_atom", "atoms": ["HO2"]}]}
    assert "'drop' to take back the last edit" in asked[0] and "undo" not in asked[0]
    assert "Dropped: " in module.processor.out.getvalue()


# ---------------------------------------------- menu 2: table and waters ---

def _biotin_ligand(module):
    module.state.ligand.source = {"kind": "structure_residue", "chain": "A", "resseq": 300, "icode": "",
                                  "resname": "BTN"}
    module.save_state = lambda: None


def test_menu_2_lists_distances_and_the_waters_near_the_ligand(module, monkeypatch):
    """Waters within a hydrogen-bond distance of the crystal ligand get rows with their polar
    contacts; any other water is switched by residue id; the proposal is said, not hidden."""
    from proprep.docking_prep.docking_menus_receptor import choose_hetero
    import proprep.structure_prep.viewer_coordinator as vc
    fake = MagicMock()
    fake.current_structures.return_value = []
    monkeypatch.setattr(vc, "viewer", fake)
    _biotin_ligand(module)
    answers = ["4", "A:311", ""]                  # keep water row 4 (A:445), and A:311 by its id
    _scripted(monkeypatch, answers, [])
    choose_hetero(module)
    assert answers == []
    shown = module.processor.out.getvalue()
    assert "Proposed: every HETATM residue other than water and the ligand, as found" in shown
    for rid in ("A:315", "A:398", "A:445"):
        assert rid in shown
    assert "BTN O11, BTN O12" in shown            # A:445 bridges biotin's carboxylate oxygens
    assert "BTN O11, A:86 O, A:112 OG" in shown    # A:398 bridges biotin and the protein
    assert sorted(module.state.receptor.keep_hetero) == ["A:311", "A:445"]
    assert "Waters: 84 in the file, 2 kept (A:311, A:445)." in shown
    assert "AD4's hydrogen-bond term does" in shown
    assert module.state.receptor.hetero_reviewed


def test_menu_2_draws_each_row_with_its_number_in_one_viewer_update(module, monkeypatch):
    from proprep.docking_prep.docking_menus_receptor import choose_hetero
    import proprep.structure_prep.viewer_coordinator as vc
    fake = MagicMock()
    fake.current_structures.return_value = []
    monkeypatch.setattr(vc, "viewer", fake)
    _biotin_ligand(module)
    _scripted(monkeypatch, [""], [])
    choose_hetero(module)
    draw, clear = fake.replace_annotations.call_args_list
    prefix, entries = draw.args
    assert prefix == "dock_het_" and clear.args == ("dock_het_", [])        # cleared on leaving
    by_label = {e["label"]: e for e in entries}
    assert by_label["dock_het_1"]["display_label"] == "#1 BTN A:300: ligand to dock"
    assert by_label["dock_het_1"]["text"] == "1"                 # the number shares the row's panel entry
    assert not [e for e in entries if e.get("style") == "tag"]
    waters = [e for e in entries if e["label"].startswith("dock_het_") and "HOH" in e.get("display_label", "")]
    assert waters and all(e["style"] == "line" for e in waters)              # left out: thin lines
    aside = fake.show_structures.call_args.kwargs["set_aside"][0]
    assert ":A and 300" in aside["selection"] and ":A and 445" in aside["selection"]


def test_replace_annotations_swaps_one_prefix_in_a_single_update():
    from proprep.structure_prep import viewer_coordinator as vc_mod
    v = MagicMock()
    v.selected_structures = ["x.pdb"]
    v.annotation_config = {"other": {"style": "halo"}, "dock_het_old": {"style": "line"}}
    c = vc_mod.ViewerCoordinator()
    c._ensure_viewer = lambda: v
    c.is_running = lambda: True
    c._owns_live_server = lambda v: True
    c.replace_annotations("dock_het_", [{"label": "dock_het_tag_1", "selection": ":A and 5", "style": "tag",
                                         "text": "1", "display_label": "#1 number"}])
    assert set(v.annotation_config) == {"other", "dock_het_tag_1"}
    assert v.annotation_config["dock_het_tag_1"]["text"] == "1"
    assert v.annotation_config["dock_het_tag_1"]["label"] == "#1 number"
    assert v.update_annotations.call_count == 1


def test_a_water_is_kept_by_clicking_it_and_every_water_is_drawn(module, monkeypatch):
    """Without a ligand from the structure no water is listed; any one can still be kept, by its
    residue id or by clicking it ('p'). Every water is drawn so it can be seen and clicked, kept
    ones larger; clicking the protein says why nothing changed."""
    from proprep.docking_prep.docking_menus_receptor import choose_hetero
    from proprep.docking_prep import docking_ui
    import proprep.structure_prep.viewer_coordinator as vc
    fake = MagicMock()
    fake.current_structures.return_value = []
    monkeypatch.setattr(vc, "viewer", fake)
    module.save_state = lambda: None
    picks = iter([{"atom": {"index": 1, "name": "N", "resname": "ASP", "resno": 128, "chain": "A", "inscode": ""}},
                  {"atom": {"index": 2, "name": "O", "resname": "HOH", "resno": 445, "chain": "A", "inscode": ""}}])
    monkeypatch.setattr(docking_ui, "pick", lambda m, kind, prompt, index: next(picks))
    _scripted(monkeypatch, ["p", "p", ""], [])
    choose_hetero(module)
    shown = module.processor.out.getvalue()
    assert "A:128 (ASP) is part of the protein, not a HETATM residue." in shown
    assert "To keep one water, type its residue id (e.g. B:312) or 'p' and click it in the viewer; " \
           "'w' keeps or leaves out all 84 at once." in shown.replace("\n", " ").replace("  ", " ")
    assert "A:445" in module.state.receptor.keep_hetero
    last_draw = [c for c in fake.replace_annotations.call_args_list if c.args[1]][-1].args[1]
    reps = {e["label"]: e for e in last_draw}
    assert reps["dock_het_waters_kept"]["selection"] == "(:A and 445)"
    assert reps["dock_het_waters_kept"]["style"] == "spacefill"
    assert reps["dock_het_waters"]["display_label"] == "Waters left out (83): 'p' picks one"


def test_option_4_says_why_a_kept_residue_is_reviewed_and_what_enter_does(module, btn, monkeypatch):
    """Walking the dashboard in order keeps the crystal ligand as a cofactor (option 3's proposal);
    option 4 then reviews it. It says the residue is in the receptor, how to dock it instead, and
    that Enter uses a consistent entry as written."""
    import proprep.docking_prep.docking_menus_receptor as menus
    import proprep.structure_prep.viewer_coordinator as vc
    monkeypatch.setattr(vc, "viewer", MagicMock())
    monkeypatch.setattr(menus, "fetch_component", lambda code: btn)
    module.save_state = lambda: None
    module.state.receptor.keep_hetero = ["A:300"]
    module.state.receptor.hetero_reviewed = True
    _scripted(monkeypatch, [""], [])
    menus.choose_metals_and_cofactors(module)
    shown = module.processor.out.getvalue().replace("\n", " ").replace("  ", " ")
    assert "BTN is kept in the receptor (A:300; option 3). If it is the ligand you mean to dock, " \
           "choose it in option 2 instead, which takes it out of the receptor." in shown
    assert "The entry is consistent. Press Enter to use it as written" in shown
    assert module.state.receptor.cofactor_edits == {"BTN": []}
