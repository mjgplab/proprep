"""
Resuming a pending parameterization must not crash.

The Force Field Parameterizer offers "Resume pending workflow" whenever
``pending_parameterizations`` is non-empty, then dispatches on the entry's
type::

    if param_type == "modified_amino_acid":
        self._resume_modified_amino_acid_workflow(...)   # never defined
    elif param_type == "small_molecule":
        self._resume_small_molecule_workflow(...)        # never defined
    elif param_type == "metal_site":
        self._resume_metal_site_workflow(...)            # imported a dead module

All three branches were broken. Two raised AttributeError; the third imported
``proprep.ff_prep``, a path predating the rename to ``forcefield_prep``, and
reported "parameterizer not available" from the caught ImportError.

The small-molecule branch is reachable in ordinary use: that workflow pauses
for Gaussian and records ``"type": "small_molecule"``, so it is the entry most
likely to be resumed.

A misspelled or missing ``self.<method>`` is a clean AttributeError at call
time and invisible before it, which is why these survived.
"""

from types import SimpleNamespace

import pytest
from rich.console import Console

from proprep.forcefield_prep.forcefield_parameterizer import (
    ForcefieldParameterizer,
)


def _parameterizer(tmp_path=None, workspace=None):
    p = ForcefieldParameterizer.__new__(ForcefieldParameterizer)
    p.console = Console(record=True, width=100)
    p.processor = SimpleNamespace()
    store = dict(workspace or {})
    p.get_from_workspace = lambda key, default=None: store.get(key, default)
    p.update_workspace = lambda key, value: store.__setitem__(key, value)
    p._store = store
    return p


# --------------------------------------------------------------------------- #
# the methods exist at all
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("name", [
    "_get_site_attribute",
    "_resume_modified_amino_acid_workflow",
    "_resume_small_molecule_workflow",
    "_resume_metal_site_workflow",
])
def test_the_dispatched_methods_exist(name):
    assert hasattr(ForcefieldParameterizer, name)


# --------------------------------------------------------------------------- #
# _get_site_attribute
# --------------------------------------------------------------------------- #

def test_site_attribute_delegates_to_the_processor():
    """It is defined on PDBProcessor; delegating keeps the two from drifting."""
    p = _parameterizer()
    p.processor = SimpleNamespace(
        _get_site_attribute=lambda site, attr, default=None: "from-processor")

    assert p._get_site_attribute({"x": 1}, "x") == "from-processor"


def test_site_attribute_reads_a_dict_without_a_processor():
    p = _parameterizer()

    assert p._get_site_attribute({"metal_element": "FE"}, "metal_element") == "FE"


def test_site_attribute_reads_an_object_without_a_processor():
    p = _parameterizer()
    site = SimpleNamespace(metal_element="ZN")

    assert p._get_site_attribute(site, "metal_element") == "ZN"


def test_site_attribute_returns_the_default_when_absent():
    p = _parameterizer()

    assert p._get_site_attribute({}, "missing", "?") == "?"


# --------------------------------------------------------------------------- #
# resuming without an output directory
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("method", [
    "_resume_modified_amino_acid_workflow",
    "_resume_small_molecule_workflow",
    "_resume_metal_site_workflow",
])
def test_a_missing_output_dir_reports_rather_than_raises(method):
    p = _parameterizer()

    getattr(p, method)("LIG", {})       # must not raise

    assert "output directory" in p.console.export_text().lower()


# --------------------------------------------------------------------------- #
# the small-molecule branch, the reachable one
# --------------------------------------------------------------------------- #

def test_small_molecule_resume_re_enters_the_workflow(monkeypatch, tmp_path):
    """
    The parameterizer is checklist-driven and keeps its state in the output
    directory, so re-entering run_workflow there resumes rather than restarts.
    """
    calls = {}

    def fake_run_workflow(**kwargs):
        calls.update(kwargs)
        return {"success": True, "status": "completed"}

    monkeypatch.setattr(
        "proprep.forcefield_prep.small_molecule_parameterizer.run_workflow",
        fake_run_workflow)

    p = _parameterizer(workspace={"pending_parameterizations": {"LIG": {}}})
    p._update_parameterization_results = lambda *a, **k: None

    p._resume_small_molecule_workflow("LIG", {"output_dir": str(tmp_path)})

    assert calls["residue_name"] == "LIG"
    assert calls["output_dir"] == str(tmp_path)
    # Not passed, so the parameter's default applies -- which must be False,
    # or "resume" would silently rebuild from scratch.
    assert not calls.get("regenerate"), "resuming must not force a rebuild"


def test_the_regenerate_default_is_what_resuming_relies_on():
    """The resume path passes no `regenerate`, so the default is the contract."""
    import inspect

    from proprep.forcefield_prep.small_molecule_parameterizer import run_workflow

    default = inspect.signature(run_workflow).parameters["regenerate"].default
    assert default is False


def test_a_completed_resume_clears_the_pending_entry(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "proprep.forcefield_prep.small_molecule_parameterizer.run_workflow",
        lambda **k: {"success": True, "status": "completed"})

    p = _parameterizer(workspace={"pending_parameterizations": {"LIG": {"type": "small_molecule"}}})
    p._update_parameterization_results = lambda *a, **k: None

    p._resume_small_molecule_workflow("LIG", {"output_dir": str(tmp_path)})

    assert "LIG" not in p._store["pending_parameterizations"]


def test_an_incomplete_resume_keeps_the_pending_entry(monkeypatch, tmp_path):
    """Still waiting on Gaussian: the entry must survive for the next attempt."""
    monkeypatch.setattr(
        "proprep.forcefield_prep.small_molecule_parameterizer.run_workflow",
        lambda **k: {"success": False, "message": "waiting for Gaussian",
                     "missing_files": ["mol.log"]})

    p = _parameterizer(workspace={"pending_parameterizations": {"LIG": {}}})

    p._resume_small_molecule_workflow("LIG", {"output_dir": str(tmp_path)})
    text = p.console.export_text()

    assert "LIG" in p._store["pending_parameterizations"]
    assert "waiting for Gaussian" in text
    assert "mol.log" in text


def test_a_raising_workflow_does_not_kill_the_menu(monkeypatch, tmp_path):
    def boom(**kwargs):
        raise RuntimeError("antechamber exploded")

    monkeypatch.setattr(
        "proprep.forcefield_prep.small_molecule_parameterizer.run_workflow", boom)
    p = _parameterizer()

    p._resume_small_molecule_workflow("LIG", {"output_dir": str(tmp_path)})

    assert "antechamber exploded" in p.console.export_text()


# --------------------------------------------------------------------------- #
# the modified-amino-acid branch
# --------------------------------------------------------------------------- #

def test_modaa_resume_uses_the_real_entry_point(monkeypatch, tmp_path):
    """resume_paused_workflow checks whether the awaited QM logs have appeared."""
    calls = {}

    def fake_resume(amino_acid, output_dir):
        calls["amino_acid"], calls["output_dir"] = amino_acid, output_dir
        return {"success": True, "status": "completed"}

    monkeypatch.setattr(
        "proprep.forcefield_prep.modified_amino_acid_parameterizer.resume_paused_workflow",
        fake_resume)

    p = _parameterizer(workspace={"pending_parameterizations": {"SEP": {}}})
    p._update_parameterization_results = lambda *a, **k: None

    p._resume_modified_amino_acid_workflow("SEP", {"output_dir": str(tmp_path)})

    assert calls == {"amino_acid": "SEP", "output_dir": str(tmp_path)}


def test_modaa_reports_the_files_it_is_still_waiting_on(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "proprep.forcefield_prep.modified_amino_acid_parameterizer.resume_paused_workflow",
        lambda amino_acid, output_dir: {
            "success": False, "message": "Workflow cannot be resumed yet.",
            "missing_files": ["SEP_opt.log", "SEP_esp.log"]})

    p = _parameterizer()
    p._resume_modified_amino_acid_workflow("SEP", {"output_dir": str(tmp_path)})
    text = p.console.export_text()

    assert "SEP_opt.log" in text and "SEP_esp.log" in text


# --------------------------------------------------------------------------- #
# the metal-site branch
# --------------------------------------------------------------------------- #

def test_metal_site_resume_points_at_the_checklist(tmp_path):
    """
    Its old import target does not exist, so rather than re-implement the
    resume it names where one actually happens.
    """
    (tmp_path / "workflow_state.json").write_text("{}")
    p = _parameterizer()

    p._resume_metal_site_workflow("FES", {"output_dir": str(tmp_path)})
    text = p.console.export_text()

    assert "checklist" in text.lower()
    assert "workflow_state.json" in text


def test_metal_site_resume_says_when_there_is_no_saved_state(tmp_path):
    p = _parameterizer()

    p._resume_metal_site_workflow("FES", {"output_dir": str(tmp_path)})

    assert "No saved state" in p.console.export_text()


def test_the_dead_module_path_is_gone():
    """proprep.ff_prep predates the rename to forcefield_prep."""
    import inspect

    source = inspect.getsource(ForcefieldParameterizer._resume_metal_site_workflow)

    assert "from proprep.ff_prep" not in source
