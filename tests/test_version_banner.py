#!/usr/bin/env python3
"""
The banner version must follow the source pin when ProPrep runs from a
checkout. An editable install records the version at install time, so the
package metadata said 1.16.0 on the 1.17.0 tree until `pip install -e .` was
re-run; pyproject.toml is the pin the lockstep check verifies.

Run with: pytest tests/test_version_banner.py
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from proprep import main as proprep_main  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent


def test_source_checkout_reports_pyproject_pin():
    pin = re.search(r'^version = "([^"]+)"', (ROOT / "pyproject.toml").read_text(), re.M).group(1)
    assert proprep_main._source_checkout_version() == pin
    assert proprep_main.get_version() == pin


def test_installed_package_falls_back_to_metadata(monkeypatch):
    monkeypatch.setattr(proprep_main, "_source_checkout_version", lambda: None)
    monkeypatch.setattr(proprep_main, "version", lambda name: "9.9.9")
    assert proprep_main.get_version() == "9.9.9"
    monkeypatch.setattr(proprep_main, "version", lambda name: (_ for _ in ()).throw(RuntimeError()))
    assert proprep_main.get_version() == "unknown"
