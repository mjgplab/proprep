"""Tests for the menu/stage splicing helpers.

WorkflowStateManager.register_plugin_stage and
WorkflowMenuCommand.register_plugin_tool mutate class-level state
to splice plugin-declared entries into the workflow menu. These
tests validate the splice without touching PDBProcessor (heavy) or
the actual menu loop (interactive).
"""

from __future__ import annotations

import copy

import pytest

from proprep.application.menu_commands import WorkflowMenuCommand
from proprep.plugins import PluginMetadata
from proprep.utils.workflow_state_manager import WorkflowStateManager


@pytest.fixture(autouse=True)
def _restore_class_state():
    """Snapshot the mutable class attributes both helpers touch and
    restore them after every test. Keeps tests independent and
    prevents order-of-execution coupling — class-level mutation is
    intentional in production but a hazard for tests."""
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


def _meta(stage="analyze", stage_name="Analyze", stage_order=60,
          name="x", tool_order=0):
    return PluginMetadata(
        name=name, version="1.0", description=f"{name} description",
        stage=stage, stage_name=stage_name, stage_order=stage_order,
        tool_order=tool_order,
    )


def _plugin():
    class _P: pass
    return _P()


# ---------------------------------------------------------------------------
# register_plugin_stage
# ---------------------------------------------------------------------------


def test_register_plugin_stage_adds_new_stage():
    """A plugin declaring a previously-unknown stage extends the
    class-level lists/dicts."""
    meta = _meta()
    WorkflowStateManager.register_plugin_stage(meta)
    assert "analyze" in WorkflowStateManager.STAGES
    assert WorkflowStateManager.STAGE_NAMES["analyze"] == "Analyze"
    assert WorkflowStateManager.STAGE_ORDER["analyze"] == 60


def test_register_plugin_stage_orders_by_stage_order():
    """Stages render in stage_order, not registration order."""
    early = _meta(stage="zearly", stage_name="ZEarly", stage_order=15)
    late = _meta(stage="zlate", stage_name="ZLate", stage_order=200)
    # Register in opposite of the desired display order.
    WorkflowStateManager.register_plugin_stage(late)
    WorkflowStateManager.register_plugin_stage(early)
    stages = WorkflowStateManager.STAGES
    # Native "load" is stage_order=10, so should still be first;
    # zearly (15) immediately after.
    assert stages.index("zearly") == 1
    assert stages.index("zlate") == len(stages) - 1


def test_register_plugin_stage_idempotent():
    """Calling twice with the same stage doesn't double-add."""
    meta = _meta()
    WorkflowStateManager.register_plugin_stage(meta)
    before = list(WorkflowStateManager.STAGES)
    WorkflowStateManager.register_plugin_stage(meta)
    assert WorkflowStateManager.STAGES == before


def test_register_plugin_stage_first_registration_wins_for_label():
    """Two plugins targeting the same new stage: first one's
    stage_name + stage_order are used; second is ignored."""
    a = _meta(stage="shared", stage_name="A label", stage_order=70, name="a")
    b = _meta(stage="shared", stage_name="B label", stage_order=70, name="b")
    WorkflowStateManager.register_plugin_stage(a)
    WorkflowStateManager.register_plugin_stage(b)
    assert WorkflowStateManager.STAGE_NAMES["shared"] == "A label"


# ---------------------------------------------------------------------------
# register_plugin_tool
# ---------------------------------------------------------------------------


def test_register_plugin_tool_creates_stage_tools_entry_when_absent():
    """If the stage has no native tools yet, register_plugin_tool
    creates the STAGE_TOOLS entry from the plugin's metadata."""
    meta = _meta()
    WorkflowMenuCommand.register_plugin_tool(meta, _plugin())
    assert "analyze" in WorkflowMenuCommand.STAGE_TOOLS
    tools = WorkflowMenuCommand.STAGE_TOOLS["analyze"]["tools"]
    assert "__plugin::x" in tools


def test_register_plugin_tool_appends_under_existing_stage():
    """A plugin can target an existing native stage; its sentinel
    is appended to that stage's tool list."""
    meta = _meta(stage="load", stage_name="Load", name="loadhelper")
    WorkflowMenuCommand.register_plugin_tool(meta, _plugin())
    assert "__plugin::loadhelper" in WorkflowMenuCommand.STAGE_TOOLS["load"]["tools"]
    # Native tools still present.
    assert "Structure Loader" in WorkflowMenuCommand.STAGE_TOOLS["load"]["tools"]


def test_register_plugin_tool_idempotent():
    """Re-registering the same plugin replaces its prior entry
    rather than duplicating."""
    meta = _meta()
    plugin = _plugin()
    WorkflowMenuCommand.register_plugin_tool(meta, plugin)
    WorkflowMenuCommand.register_plugin_tool(meta, plugin)
    tools = WorkflowMenuCommand.STAGE_TOOLS["analyze"]["tools"]
    assert tools.count("__plugin::x") == 1


def test_register_plugin_tool_records_instance_and_metadata():
    """Dispatcher needs to look up the plugin object + its metadata
    by name — these get recorded in the class-level lookup tables."""
    meta = _meta()
    plugin = _plugin()
    WorkflowMenuCommand.register_plugin_tool(meta, plugin)
    assert WorkflowMenuCommand._plugin_instances["x"] is plugin
    assert WorkflowMenuCommand._plugin_metadata["x"] is meta


def test_register_plugin_tool_orders_by_tool_order():
    """Multiple plugins under one stage: tools render in tool_order
    (lower first), with native tools always before plugin tools."""
    a = _meta(stage="analyze", name="a", tool_order=20)
    b = _meta(stage="analyze", name="b", tool_order=10)
    c = _meta(stage="analyze", name="c", tool_order=15)
    WorkflowMenuCommand.register_plugin_tool(a, _plugin())
    WorkflowMenuCommand.register_plugin_tool(b, _plugin())
    WorkflowMenuCommand.register_plugin_tool(c, _plugin())
    tools = WorkflowMenuCommand.STAGE_TOOLS["analyze"]["tools"]
    plugins_only = [t for t in tools if t.startswith("__plugin::")]
    assert plugins_only == ["__plugin::b", "__plugin::c", "__plugin::a"]
