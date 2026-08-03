"""Regression tests for batch-replay web/terminal mode matching.

A template recorded under ``proprep-web`` (PROPREP_WEB_SHELL set) skips some
mode-gated prompts (e.g. the PDB Filter "view structure?" Y/N). Replaying such
a template under plain ``proprep`` would re-introduce those prompts and desync
the recorded inputs. The batch runner therefore recreates the recorded launch
environment (and stays headless via PROPREP_BATCH). ``--web``/``--no-web``
override the recorded mode.
"""

import json
import os
import tempfile
from unittest import mock

import pytest

from proprep.utils import batch_processor as bp


def _make_template(web_shell_mode):
    d = tempfile.mkdtemp()
    tpl = os.path.join(d, "sess_template.json")
    with open(tpl, "w") as f:
        json.dump(
            {
                "template": True,
                "template_variables": {"input_protein": {"description": "x"}},
                "metadata": {"web_shell_mode": web_shell_mode},
                "interactions": [],
            },
            f,
        )
    lst = os.path.join(d, "list.txt")
    with open(lst, "w") as f:
        f.write("1abc\n")
    return tpl, lst


def _run_capturing_env(web_shell_mode, web_override):
    """Run a one-item batch, capturing the env visible during replay."""
    tpl, lst = _make_template(web_shell_mode)
    captured = {}

    def fake_single(self, template_file, variables, base_dir, template_vars, run_index):
        captured["WEB"] = os.environ.get("PROPREP_WEB_SHELL")
        captured["BATCH"] = os.environ.get("PROPREP_BATCH")
        return True, "out", None

    os.environ.pop("PROPREP_WEB_SHELL", None)
    os.environ.pop("PROPREP_BATCH", None)
    with mock.patch.object(bp.BatchProcessor, "_process_single_protein", fake_single):
        rc = bp.run_batch_processing(
            tpl, lst, base_dir=tempfile.mkdtemp(),
            continue_on_error=True, web_override=web_override,
        )
    return rc, captured


def test_web_recorded_template_replays_in_web_mode():
    rc, cap = _run_capturing_env(web_shell_mode=True, web_override=None)
    assert rc == 0
    assert cap["WEB"] == "1"      # PROPREP_WEB_SHELL set to match recording
    assert cap["BATCH"] == "1"    # headless guard always on during batch


def test_terminal_recorded_template_replays_in_terminal_mode():
    rc, cap = _run_capturing_env(web_shell_mode=False, web_override=None)
    assert rc == 0
    assert cap["WEB"] is None     # not in web mode
    assert cap["BATCH"] == "1"


def test_web_override_forces_web_on_terminal_template():
    rc, cap = _run_capturing_env(web_shell_mode=False, web_override=True)
    assert cap["WEB"] == "1"


def test_no_web_override_forces_terminal_on_web_template():
    rc, cap = _run_capturing_env(web_shell_mode=True, web_override=False)
    assert cap["WEB"] is None


def test_missing_metadata_defaults_to_terminal():
    # Older templates predate web_shell_mode; treat as terminal.
    rc, cap = _run_capturing_env(web_shell_mode=None, web_override=None)
    assert cap["WEB"] is None
    assert cap["BATCH"] == "1"


def test_environment_restored_after_batch():
    # Pre-existing PROPREP_WEB_SHELL must be restored verbatim.
    os.environ["PROPREP_WEB_SHELL"] = "preexisting"
    try:
        tpl, lst = _make_template(False)
        with mock.patch.object(
            bp.BatchProcessor, "_process_single_protein",
            lambda *a, **k: (True, "out", None),
        ):
            bp.run_batch_processing(tpl, lst, base_dir=tempfile.mkdtemp())
        assert os.environ.get("PROPREP_WEB_SHELL") == "preexisting"
        assert "PROPREP_BATCH" not in os.environ
    finally:
        os.environ.pop("PROPREP_WEB_SHELL", None)


def test_launch_viewer_is_noop_under_batch():
    from proprep.structure_prep.interactive_structure_viewer import (
        InteractiveStructureViewer,
    )

    v = InteractiveStructureViewer()
    v.selected_structures = ["/nonexistent.pdb"]
    os.environ["PROPREP_BATCH"] = "1"
    try:
        assert v._launch_viewer() is False
    finally:
        os.environ.pop("PROPREP_BATCH", None)
