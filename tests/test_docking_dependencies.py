#!/usr/bin/env python3
"""
Docking dependency report: a missing package is named with what it is needed
for and how to install it, instead of failing on the first import.

Run with: pytest tests/test_docking_dependencies.py
"""

import importlib

from proprep.docking_prep import dependencies


def _report_with_missing(monkeypatch, missing_module):
    real_import = importlib.import_module

    def fake_import(name, *args, **kwargs):
        if name == missing_module:
            raise ImportError(f"No module named '{missing_module}'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(dependencies.importlib, "import_module", fake_import)
    monkeypatch.setattr(dependencies, "find_executable", lambda name: None)
    return {d.name: d for d in dependencies.check_dependencies()}


def test_missing_module_is_reported_with_purpose_and_install_hint(monkeypatch):
    report = _report_with_missing(monkeypatch, "meeko")
    meeko = report["Meeko"]
    assert not meeko.available
    text = meeko.describe()
    assert "MISSING" in text and "PDBQT" in text and "No module named 'meeko'" in text
    assert "Install:" in text and "prody" in text


def test_autogrid_blocks_only_when_ad4_scoring_is_wanted(monkeypatch):
    report = list(_report_with_missing(monkeypatch, "no_such_module").values())
    names_vina = {d.name for d in dependencies.missing(report, include_ad4=False)}
    names_ad4 = {d.name for d in dependencies.missing(report, include_ad4=True)}
    assert "autogrid4" not in names_vina
    assert "autogrid4" in names_ad4


def test_version_mismatch_is_reported_not_refused():
    dep = dependencies.Dependency("X", "testing", "1.0", available=True, version="2.0")
    assert dep.describe() == "X 2.0 (validated with 1.0)"
    same = dependencies.Dependency("X", "testing", "1.0", available=True, version="1.0")
    assert same.describe() == "X 1.0"
