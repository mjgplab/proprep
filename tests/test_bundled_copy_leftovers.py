"""Files left in the installed package by AmberTools' bundled ProPrep are removed.

AmberTools 26 ships ProPrep 1.0.0 in the same site-packages/proprep folder;
installing ProPrep over it left 35 files only the old copy had (a retired
built-in MD workflow the MD Manager offered, old templates and parameters).
proprep.utils.bundled_copy removes every file under the package folder that
conda's record of the installed ProPrep package does not list, and nothing
when there is no such record.
"""

import json
import os

from proprep.utils import bundled_copy


def _install(tmp_path, owned, extra, record=True):
    prefix = tmp_path / "env"
    package = prefix / "lib" / "python3.12" / "site-packages" / "proprep"
    for rel in owned + extra:
        path = package / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("x")
    (package / "md_prep" / "__pycache__").mkdir(parents=True, exist_ok=True)
    (package / "md_prep" / "__pycache__" / "old.cpython-312.pyc").write_text("x")
    (prefix / "conda-meta").mkdir(parents=True)
    (prefix / "conda-meta" / "ambertools-dac-26.0.0-py312_0.json").write_text(json.dumps({"files": []}))
    if record:
        files = [f"lib/python3.12/site-packages/proprep/{rel}" for rel in owned]
        (prefix / "conda-meta" / "proprep-1.23.0-py_0.json").write_text(json.dumps(
            {"name": "proprep", "version": "1.23.0", "build": "py_0", "files": files}))
    return str(prefix), str(package)


OWNED = ["__init__.py", "md_workflows/builtin/basic_equilibration.json", "md_prep/workflow_loader.py"]
EXTRA = ["md_workflows/builtin/protein_equilibration.json",
         "md_templates/builtin/standard_protein/00_initial_minimization.mdin",
         "md_prep/amber_input_generator_backup.py"]


def test_files_the_package_does_not_own_are_removed(tmp_path):
    prefix, package = _install(tmp_path, OWNED, EXTRA)
    record, removed = bundled_copy.remove_files_not_in_package(prefix, package)
    assert record == "proprep-1.23.0-py_0"
    assert sorted(os.path.relpath(p, package) for p in removed) == sorted(EXTRA)
    for rel in OWNED:
        assert os.path.exists(os.path.join(package, rel))
    for rel in EXTRA:
        assert not os.path.exists(os.path.join(package, rel))
    assert not os.path.exists(os.path.join(package, "md_templates"))     # emptied folders go too
    assert os.path.exists(os.path.join(package, "md_prep", "__pycache__", "old.cpython-312.pyc"))   # caches untouched


def test_a_dry_run_lists_without_removing(tmp_path):
    prefix, package = _install(tmp_path, OWNED, EXTRA)
    _, listed = bundled_copy.remove_files_not_in_package(prefix, package, dry_run=True)
    assert len(listed) == 3 and all(os.path.exists(p) for p in listed)


def test_nothing_is_removed_without_conda_s_record_of_proprep(tmp_path):
    """A pip or editable install, or the AmberTools build: no proprep record, nothing to judge by."""
    prefix, package = _install(tmp_path, OWNED, EXTRA, record=False)
    record, removed = bundled_copy.remove_files_not_in_package(prefix, package)
    assert record is None and removed == []
    assert all(os.path.exists(os.path.join(package, rel)) for rel in OWNED + EXTRA)


def test_the_command_reports_what_it_removed(tmp_path, monkeypatch, capsys):
    prefix, package = _install(tmp_path, OWNED, EXTRA)
    monkeypatch.setattr(bundled_copy.sys, "prefix", prefix)
    monkeypatch.setattr(bundled_copy, "_package_dir", lambda: package)
    assert bundled_copy.main([]) == 0
    out = capsys.readouterr().out
    assert "Removed 3 file(s) left in the ProPrep package by another copy (not part of proprep-1.23.0-py_0)" in out
    assert bundled_copy.main([]) == 0
    assert "Every file under the ProPrep package belongs to proprep-1.23.0-py_0." in capsys.readouterr().out
