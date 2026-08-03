"""Tests for plugin-contributed nav-footer shortcuts and the
display-name fallback in the workflow menu.

The renderer + nav-footer logic in ``WorkflowMenuCommand._show_menu``
reads three new pieces of information from PluginMetadata:

* ``display_name`` — branded label rendered in the tool list (falls
  back to the lowercase ``name`` when unset)
* ``stage_shortcut`` — single-letter key bracketed into the nav
  footer alongside the native [l]oad / [d]etect / etc. (silently
  dropped when colliding with an existing shortcut)

These tests don't drive the interactive menu loop — they exercise
the helper that formats the bracketed labels and the registration
flow that feeds the renderer.
"""

from __future__ import annotations

import copy

import pytest

from proprep.application.menu_commands import (
    WorkflowMenuCommand,
    _format_stage_shortcut_label,
)
from proprep.plugins import PluginMetadata
from proprep.utils.workflow_state_manager import WorkflowStateManager


@pytest.fixture(autouse=True)
def _restore_class_state():
    """Same snapshot/restore as test_plugins_menu_splicing — needed
    here too because we mutate the same class attributes."""
    saved = {
        "stages": list(WorkflowStateManager.STAGES),
        "stage_names": dict(WorkflowStateManager.STAGE_NAMES),
        "stage_order": dict(WorkflowStateManager.STAGE_ORDER),
        "stage_tools": copy.deepcopy(WorkflowMenuCommand.STAGE_TOOLS),
        "plugin_instances": dict(WorkflowMenuCommand._plugin_instances),
        "plugin_metadata": dict(WorkflowMenuCommand._plugin_metadata),
    }
    yield
    WorkflowStateManager.STAGES = saved["stages"]
    WorkflowStateManager.STAGE_NAMES = saved["stage_names"]
    WorkflowStateManager.STAGE_ORDER = saved["stage_order"]
    WorkflowMenuCommand.STAGE_TOOLS = saved["stage_tools"]
    WorkflowMenuCommand._plugin_instances = saved["plugin_instances"]
    WorkflowMenuCommand._plugin_metadata = saved["plugin_metadata"]


# ---------------------------------------------------------------------------
# _format_stage_shortcut_label
# ---------------------------------------------------------------------------


def test_format_label_brackets_letter_at_first_occurrence():
    """When the shortcut letter appears in the stage name, bracket
    its first case-insensitive occurrence — matches the native
    style ``\\[l]oad``."""
    assert _format_stage_shortcut_label("Analyze", "z") == "analy\\[z]e"
    assert _format_stage_shortcut_label("Analyze", "a") == "\\[a]nalyze"
    assert _format_stage_shortcut_label("Detect", "d")  == "\\[d]etect"


def test_format_label_appends_when_letter_absent():
    """Falls back to ``<name> \\[<letter>]`` so the user still sees
    which key to press, even when no plain bracketing fits."""
    assert _format_stage_shortcut_label("Cooper", "x") == "cooper \\[x]"


def test_format_label_is_case_insensitive_on_input():
    """Plugin authors may write the shortcut as 'Z' or 'z' — both
    bracket the first lowercase 'z' in the stage name."""
    assert _format_stage_shortcut_label("Analyze", "Z") == "analy\\[z]e"


# ---------------------------------------------------------------------------
# Plugin metadata wires through registration
# ---------------------------------------------------------------------------


def _meta(**overrides):
    base = dict(
        name="x", version="0.0", description="d",
        stage="analyze", stage_name="Analyze", stage_order=60,
    )
    base.update(overrides)
    return PluginMetadata(**base)


def _plugin():
    class _P: pass
    return _P()


def test_register_plugin_tool_records_display_name_via_metadata():
    """The renderer reads display_name off
    ``WorkflowMenuCommand._plugin_metadata[name]`` — so registration
    just has to store the metadata. This pins that contract: a
    plugin with display_name set is retrievable post-registration."""
    meta = _meta(name="etanalyze", display_name="ETAnalyze")
    WorkflowMenuCommand.register_plugin_tool(meta, _plugin())
    stored = WorkflowMenuCommand._plugin_metadata["etanalyze"]
    assert stored.display_name == "ETAnalyze"


def test_register_plugin_tool_default_display_name_is_none():
    """When the plugin omits display_name, the renderer's fallback
    path (lowercase ``name``) kicks in. Validates the stored
    metadata reflects the omission rather than auto-deriving."""
    meta = _meta(name="quickplugin")
    WorkflowMenuCommand.register_plugin_tool(meta, _plugin())
    stored = WorkflowMenuCommand._plugin_metadata["quickplugin"]
    assert stored.display_name is None
