"""A batch download fills rcsb_pdb_files and also points rcsb_pdb_file at the last file, the
one ProPrep's other tools use. The selector listed that structure twice: downloading 2RVD and
1UAO gave the viewer a table of three."""

import io

from rich.console import Console

import proprep.utils.prompts as prompts
from proprep.utils.structure_selector import StructureSelector


def _workspace(tmp_path, single, batch, **more):
    paths = {}
    for name in set([single, *batch, *more.values()]):
        path = tmp_path / f"{name}.pdb"
        path.write_text("ATOM\n")
        paths[name] = str(path)
    workspace = {"rcsb_pdb_file": paths[single], "rcsb_pdb_files": [paths[n] for n in batch],
                 "rcsb_download_info": [{"pdb_id": n} for n in batch]}
    workspace.update({key: paths[name] for key, name in more.items()})
    return workspace, paths


def _rows(selector):
    return [(i.structure_type.workspace_key, i.structure_type.display_name, i.is_current)
            for i in selector.get_available_structures()]


def test_the_structure_held_by_both_keys_is_listed_once_and_marked_current(tmp_path):
    workspace, _ = _workspace(tmp_path, "1UAO", ["2RVD", "1UAO"])
    assert _rows(StructureSelector(workspace)) == [
        ("rcsb_pdb_file", "RCSB PDB (1UAO)", True),          # the single key's row, now named
        ("rcsb_pdb_files[0]", "RCSB PDB (2RVD)", False)]


def test_other_structures_and_their_order_are_untouched(tmp_path):
    workspace, _ = _workspace(tmp_path, "1UAO", ["2RVD", "1UAO"], local_pdb_file="mine")
    keys = [row[0] for row in _rows(StructureSelector(workspace))]
    assert keys == ["rcsb_pdb_file", "rcsb_pdb_files[0]", "local_pdb_file"]


def test_nothing_is_folded_when_the_single_key_holds_a_different_structure(tmp_path):
    workspace, _ = _workspace(tmp_path, "4XYZ", ["2RVD", "1UAO"])     # a later single download
    assert _rows(StructureSelector(workspace)) == [
        ("rcsb_pdb_file", "RCSB PDB", False),
        ("rcsb_pdb_files[0]", "RCSB PDB (2RVD)", False),
        ("rcsb_pdb_files[1]", "RCSB PDB (1UAO)", False)]


def test_the_same_file_by_another_path_is_still_one_structure(tmp_path):
    workspace, paths = _workspace(tmp_path, "1UAO", ["2RVD", "1UAO"])
    link = tmp_path / "link.pdb"
    link.symlink_to(paths["1UAO"])
    workspace["rcsb_pdb_file"] = str(link)
    assert len(_rows(StructureSelector(workspace))) == 2


def test_priority_selection_still_answers_with_the_single_key(tmp_path):
    workspace, paths = _workspace(tmp_path, "1UAO", ["2RVD", "1UAO"])
    console = Console(file=io.StringIO())
    selector = StructureSelector(workspace, console=console)
    assert selector._priority_selection(selector.get_available_structures(), silent=True, return_key=True) == \
        (paths["1UAO"], "rcsb_pdb_file")


def test_table_says_current_but_the_recorded_label_does_not(tmp_path, monkeypatch):
    """Replay resolves a choice by its label. The merged row's label is exactly the batch item's,
    so a session that chose the old duplicate row (3) lands on it; the old single-key row was 1
    and still is."""
    import types
    import proprep.application.menu_commands as menu_commands
    from proprep.utils.prompts import resolve_replayed_key

    workspace, paths = _workspace(tmp_path, "1UAO", ["2RVD", "1UAO"])
    console = Console(width=250, record=True, file=io.StringIO(), force_terminal=False)
    selector = StructureSelector(workspace, console=console, processor=object())
    seen = {}

    def fake_prompt(processor=None, prompt=None, **kwargs):
        seen.update(kwargs)
        return "1"
    monkeypatch.setattr(menu_commands, "prompt_with_context", fake_prompt)

    chosen = selector._interactive_multi_selection(selector.get_available_structures())

    assert chosen == [paths["1UAO"]]
    out = console.export_text()
    assert "RCSB PDB (1UAO), current" in out
    assert "current: of the structures downloaded together" in out        # the word is explained where it is used
    assert seen["options_map"] == {"1": f"RCSB PDB (1UAO) ({paths['1UAO']})",
                                   "2": f"RCSB PDB (2RVD) ({paths['2RVD']})",
                                   "all": "All available structures"}

    # a session recorded when the table had three rows: "3" was the duplicate of 1UAO
    def replaying(label):
        replayer = types.SimpleNamespace(replaying=True, last_returned_interaction={"context": {"option_label": label}})
        return types.SimpleNamespace(session_manager=types.SimpleNamespace(replayer=replayer))
    assert resolve_replayed_key(replaying(f"RCSB PDB (1UAO) ({paths['1UAO']})"), "3", seen["options_map"]) == "1"
    assert resolve_replayed_key(replaying(f"RCSB PDB ({paths['1UAO']})"), "1", seen["options_map"]) == "1"
    assert resolve_replayed_key(replaying(f"RCSB PDB (2RVD) ({paths['2RVD']})"), "2", seen["options_map"]) == "2"
