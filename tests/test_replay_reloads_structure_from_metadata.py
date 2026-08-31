"""A session recorded from `proprep --pdbid X` / `--pdbfile F` has no loader
prompts; replay must reload the structure from the session metadata."""
import json
import os

from proprep import main as proprep_main


class _Loader:
    def __init__(self):
        self.calls = []
    def _download_and_load_pdb(self, pdb_id, workspace):
        self.calls.append(("id", pdb_id))
    def _load_local_file_by_path(self, path, workspace):
        self.calls.append(("file", path))


class _Proc:
    def __init__(self, project_dir):
        self.loader = _Loader()
        self._ws = {"project_directory": project_dir}
    def get_module_instance(self, name):
        assert name == "Structure Loader"; return self.loader
    def _get_workspace(self):
        return self._ws


def _session(tmp_path, metadata):
    f = tmp_path / "proprep_session_x.json"
    f.write_text(json.dumps({"metadata": metadata, "interactions": []}))
    return str(f)


def test_pdbid_session_redownloads(tmp_path):
    proc = _Proc(str(tmp_path))
    assert proprep_main.reload_structure_from_session_metadata(proc, _session(tmp_path, {"pdb_id": "1M1Q"}))
    assert proc.loader.calls == [("id", "1M1Q")]


def test_pdbfile_from_another_machine_is_found_by_basename(tmp_path):
    local = tmp_path / "1M1Q.pdb"; local.write_text("END\n")
    proc = _Proc(str(tmp_path))
    meta = {"pdb_file": "/Users/someone/elsewhere/1M1Q.pdb"}      # absolute path from the other Mac
    assert proprep_main.reload_structure_from_session_metadata(proc, _session(tmp_path, meta))
    assert proc.loader.calls == [("file", str(local))]


def test_missing_pdbfile_reports_and_loads_nothing(tmp_path):
    proc = _Proc(str(tmp_path))
    meta = {"pdb_file": "/nowhere/gone.pdb"}
    assert proprep_main.reload_structure_from_session_metadata(proc, _session(tmp_path, meta)) is False
    assert proc.loader.calls == []


def test_interactive_session_has_nothing_to_reload(tmp_path):
    proc = _Proc(str(tmp_path))
    assert proprep_main.reload_structure_from_session_metadata(proc, _session(tmp_path, {"auto_recorded": True})) is False
    assert proc.loader.calls == []
