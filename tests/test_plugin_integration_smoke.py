"""End-to-end smoke test for the plugin pipeline.

Exercises the *whole* host→plugin chain in one place rather than
slicing it into protocol / discovery / dispatch parts (those live in
their own files). Drives discover_plugins → register_plugin_stage →
register_plugin_tool → menu sentinel detection → HostContext
construction → plugin.launch() invocation, and asserts the plugin
ran with the expected seed.

Doesn't require etanalyze to be installed — uses an in-memory mock
plugin that returns a known PluginMetadata and records its launch
arguments. The point is to assert ProPrep's wiring is correct,
not to re-test etanalyze's behaviour (which has its own suite).
"""

from __future__ import annotations

import copy
import logging
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Optional
from unittest.mock import MagicMock, patch

import pytest

from proprep.application.menu_commands import WorkflowMenuCommand
from proprep.application.pdbprocessor import PDBProcessor
from proprep.plugins import HostContext, PluginMetadata
from proprep.utils.workflow_state_manager import WorkflowStateManager


# ---------------------------------------------------------------------------
# Test fixtures
# ---------------------------------------------------------------------------


class _SmokePlugin:
    """A plugin object that records what was passed into launch().

    Discovery → menu splice → dispatch all flow through this; the
    HostContext we capture in ``last_host`` is what real plugins
    would receive in production.
    """

    def __init__(self):
        self.last_host: Optional[HostContext] = None
        self.launch_count = 0

    def get_metadata(self) -> PluginMetadata:
        return PluginMetadata(
            name="smoke",
            version="0.1",
            description="Smoke-test plugin",
            stage="analyze",
            stage_name="Analyze (smoke)",
            stage_order=60,
            tool_order=0,
            consumes_workspace_keys=[
                "microstate_metadata_path", "structure_path",
            ],
        )

    def is_available(self) -> bool:
        return True

    def launch(self, host: HostContext) -> None:
        self.launch_count += 1
        self.last_host = host


@pytest.fixture(autouse=True)
def _restore_class_state():
    """Snapshot + restore class-level state mutated by registration
    helpers. Without this, smoke tests leak state into the
    test_plugins_* tests (and vice versa)."""
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


def _fake_processor_with_seed(seed_data: Dict[str, Any]):
    """Stand-in for PDBProcessor that exposes just enough surface for
    _dispatch_plugin: a workspace + console + session_manager. We don't
    instantiate the real PDBProcessor (which imports half of ProPrep's
    heavy modules); we exercise the dispatcher's contract directly."""
    workspace = SimpleNamespace(
        has=lambda key: key in seed_data,
        get=lambda key, default=None: seed_data.get(key, default),
    )
    return SimpleNamespace(
        workspace=workspace,
        console=MagicMock(),
        session_manager=SimpleNamespace(
            recorder=None,
            is_recording=lambda: False,
        ),
        registry=SimpleNamespace(modules={}),
        # build_plugin_seed is the real method bound to our stub.
        build_plugin_seed=(
            lambda meta: PDBProcessor.build_plugin_seed(
                SimpleNamespace(workspace=workspace), meta,
            )
        ),
    )


# ---------------------------------------------------------------------------
# Discovery → registration → menu splicing
# ---------------------------------------------------------------------------


def test_smoke_plugin_registers_into_workflow_menu():
    """When discovery returns a plugin, both register helpers wire
    it into the class-level state the menu reads at render time."""
    plugin = _SmokePlugin()
    meta = plugin.get_metadata()

    WorkflowStateManager.register_plugin_stage(meta)
    WorkflowMenuCommand.register_plugin_tool(meta, plugin)

    # Stage extended.
    assert "analyze" in WorkflowStateManager.STAGES
    assert WorkflowStateManager.STAGE_NAMES["analyze"] == "Analyze (smoke)"
    # Tool present under the right stage.
    tools = WorkflowMenuCommand.STAGE_TOOLS["analyze"]["tools"]
    assert "__plugin::smoke" in tools
    # Dispatcher can find the instance + metadata back out.
    assert WorkflowMenuCommand._plugin_instances["smoke"] is plugin
    assert WorkflowMenuCommand._plugin_metadata["smoke"] is meta


def test_full_pipeline_with_mocked_entry_points():
    """Patch entry_points → discover_plugins → register both
    helpers → assert the menu surfaces the plugin's tool."""
    plugin = _SmokePlugin()
    fake_ep = SimpleNamespace(name="smoke", load=lambda: plugin)

    with patch(
        "proprep.plugins.discovery.entry_points",
        return_value=[fake_ep],
    ):
        from proprep.plugins import discover_plugins
        discovered = discover_plugins()
        assert "smoke" in discovered
        # Replicate what PDBProcessor._discover_and_register_plugins
        # does after discovery.
        for name, p in discovered.items():
            meta = p.get_metadata()
            WorkflowStateManager.register_plugin_stage(meta)
            WorkflowMenuCommand.register_plugin_tool(meta, p)

    assert "__plugin::smoke" in WorkflowMenuCommand.STAGE_TOOLS["analyze"]["tools"]


# ---------------------------------------------------------------------------
# Menu dispatch → HostContext construction → plugin.launch
# ---------------------------------------------------------------------------


def test_dispatch_builds_seed_from_declared_keys_and_invokes_launch():
    """The dispatcher's job: look up the plugin, pull declared keys
    out of the host workspace, build a HostContext, call launch().
    Assert the plugin received a HostContext containing only the
    keys it declared (and no others)."""
    plugin = _SmokePlugin()
    meta = plugin.get_metadata()
    WorkflowStateManager.register_plugin_stage(meta)
    WorkflowMenuCommand.register_plugin_tool(meta, plugin)

    seed_data = {
        "microstate_metadata_path": "/p/meta.json",
        "structure_path": "/p/struct.pdb",
        "irrelevant_key": "ignored",  # not declared → must NOT leak
    }

    proc = _fake_processor_with_seed(seed_data)
    # console and workspace are read-only properties that read from
    # self.processor; just set processor and they pick up.
    cmd = WorkflowMenuCommand.__new__(WorkflowMenuCommand)
    cmd.processor = proc

    cmd._dispatch_plugin("smoke")

    assert plugin.launch_count == 1
    host = plugin.last_host
    assert isinstance(host, HostContext)
    assert host.seed_state == {
        "microstate_metadata_path": "/p/meta.json",
        "structure_path": "/p/struct.pdb",
    }
    # Logger named per the plan's locked convention.
    assert host.logger.name == "proprep.plugin.smoke"


def test_dispatch_isolates_plugin_exceptions():
    """A plugin that raises during launch must not crash ProPrep —
    the dispatcher catches and surfaces via the console."""
    class _BoomPlugin(_SmokePlugin):
        def launch(self, host):
            raise RuntimeError("synthetic plugin crash")

    plugin = _BoomPlugin()
    meta = plugin.get_metadata()
    WorkflowStateManager.register_plugin_stage(meta)
    WorkflowMenuCommand.register_plugin_tool(meta, plugin)

    proc = _fake_processor_with_seed({})
    # console and workspace are read-only properties that read from
    # self.processor; just set processor and they pick up.
    cmd = WorkflowMenuCommand.__new__(WorkflowMenuCommand)
    cmd.processor = proc

    # Must NOT raise — the dispatcher swallows + prints a panel.
    cmd._dispatch_plugin("smoke")

    # Console.print was called with a Panel instance describing the
    # error. We don't pin the exact format, just that something was
    # printed.
    assert proc.console.print.called


def test_dispatch_unknown_plugin_warns_without_crashing():
    """Defensive: if the menu somehow asks to dispatch a name not
    in the lookup tables (stale registration, race), print a warning
    rather than KeyError."""
    proc = _fake_processor_with_seed({})
    # console and workspace are read-only properties that read from
    # self.processor; just set processor and they pick up.
    cmd = WorkflowMenuCommand.__new__(WorkflowMenuCommand)
    cmd.processor = proc

    cmd._dispatch_plugin("nonexistent")

    assert proc.console.print.called


# ---------------------------------------------------------------------------
# HostContext.session_recorder hand-off
# ---------------------------------------------------------------------------


def test_dispatch_passes_host_recorder_when_recording():
    """When the host is recording, its recorder rides through to
    the plugin via HostContext.session_recorder. The plugin (or
    its adapter) decides whether to wrap or use directly."""
    plugin = _SmokePlugin()
    meta = plugin.get_metadata()
    WorkflowStateManager.register_plugin_stage(meta)
    WorkflowMenuCommand.register_plugin_tool(meta, plugin)

    fake_recorder = MagicMock(name="host_recorder")
    proc = _fake_processor_with_seed({})
    proc.session_manager.is_recording = lambda: True
    proc.session_manager.recorder = fake_recorder

    # console and workspace are read-only properties that read from
    # self.processor; just set processor and they pick up.
    cmd = WorkflowMenuCommand.__new__(WorkflowMenuCommand)
    cmd.processor = proc

    cmd._dispatch_plugin("smoke")

    assert plugin.last_host.session_recorder is fake_recorder


def test_dispatch_passes_none_recorder_when_not_recording():
    """No active session → HostContext.session_recorder is None,
    plugin runs without scope annotation."""
    plugin = _SmokePlugin()
    meta = plugin.get_metadata()
    WorkflowStateManager.register_plugin_stage(meta)
    WorkflowMenuCommand.register_plugin_tool(meta, plugin)

    proc = _fake_processor_with_seed({})
    proc.session_manager.is_recording = lambda: False

    # console and workspace are read-only properties that read from
    # self.processor; just set processor and they pick up.
    cmd = WorkflowMenuCommand.__new__(WorkflowMenuCommand)
    cmd.processor = proc

    cmd._dispatch_plugin("smoke")

    assert plugin.last_host.session_recorder is None


# ---------------------------------------------------------------------------
# Full-menu mode (MainMenuCommand) — plugins must work identically
# ---------------------------------------------------------------------------


def test_full_menu_dispatch_invokes_shared_helper():
    """``MainMenuCommand._handle_choice`` recognises the
    ``__plugin::<name>`` sentinel and routes through the shared
    dispatcher, so plugins behave identically across menu modes.
    Regression guard against the original implementation that only
    wired the workflow menu."""
    from proprep.application.menu_commands import MainMenuCommand

    plugin = _SmokePlugin()
    meta = plugin.get_metadata()
    WorkflowStateManager.register_plugin_stage(meta)
    WorkflowMenuCommand.register_plugin_tool(meta, plugin)

    proc = _fake_processor_with_seed({})

    cmd = MainMenuCommand.__new__(MainMenuCommand)
    cmd.processor = proc
    cmd.options_map = {"6a": "__plugin::smoke"}

    result = cmd._handle_choice("6a")

    # Dispatcher returned True (continue menu loop) and the plugin
    # was invoked with a HostContext.
    assert result is True
    assert plugin.launch_count == 1
    assert plugin.last_host is not None


def test_full_menu_render_appends_plugin_section_after_natives():
    """Plugins get their own section in the full-menu module_groups,
    inserted after the last native section (5) and before UTILITIES.
    Native section numbering must stay stable regardless of which
    plugins are installed — users learn 1a/2c/etc. and shouldn't
    see those keys shift."""
    from proprep.application.menu_commands import MainMenuCommand
    from rich.console import Console

    plugin = _SmokePlugin()
    meta = plugin.get_metadata()
    WorkflowStateManager.register_plugin_stage(meta)
    WorkflowMenuCommand.register_plugin_tool(meta, plugin)

    fake_workspace = SimpleNamespace(
        has=lambda k: False, get=lambda k, d=None: d,
    )
    proc = SimpleNamespace(
        workspace=fake_workspace,
        _get_workspace=lambda: fake_workspace,
        registry=SimpleNamespace(modules={}),
        get_module_instance=lambda name: None,
        console=Console(record=True, width=120),
        session_manager=SimpleNamespace(
            recorder=None, is_recording=lambda: False,
        ),
    )

    cmd = MainMenuCommand.__new__(MainMenuCommand)
    cmd.processor = proc
    # _show_menu populates options_map and writes to the recorded
    # console; drive it directly without spinning up the menu loop.
    cmd._show_menu()

    # Native option keys still 1a..5x; plugin lands at 6a (one past
    # the last native section's number).
    assert "6a" in cmd.options_map
    assert cmd.options_map["6a"] == "__plugin::smoke"

    rendered = proc.console.export_text()
    # Plugin section title is the stage name uppercased.
    assert "ANALYZE" in rendered.upper()


def test_full_menu_with_no_plugins_is_unchanged():
    """When no plugins are installed, the full menu's options_map
    has no plugin entries. Native-only workflows must not regress."""
    from proprep.application.menu_commands import MainMenuCommand
    from rich.console import Console

    # _restore_class_state autouse fixture has cleared plugin tables.
    assert WorkflowMenuCommand._plugin_metadata == {}

    fake_workspace = SimpleNamespace(
        has=lambda k: False, get=lambda k, d=None: d,
    )
    proc = SimpleNamespace(
        workspace=fake_workspace,
        _get_workspace=lambda: fake_workspace,
        registry=SimpleNamespace(modules={}),
        get_module_instance=lambda name: None,
        console=Console(record=True, width=120),
        session_manager=SimpleNamespace(
            recorder=None, is_recording=lambda: False,
        ),
    )

    cmd = MainMenuCommand.__new__(MainMenuCommand)
    cmd.processor = proc
    cmd._show_menu()

    # No plugin section means no 6a key (only natives 1a..5x +
    # utility ua + control keys).
    assert "6a" not in cmd.options_map
    # And no sentinel ever lands in options_map values.
    assert not any(
        v.startswith("__plugin::") for v in cmd.options_map.values()
    )


def test_workflow_dispatch_still_works_after_helper_extraction():
    """Regression guard: pulling _dispatch_plugin into a module-
    level helper must not break WorkflowMenuCommand's existing
    dispatch contract."""
    plugin = _SmokePlugin()
    meta = plugin.get_metadata()
    WorkflowStateManager.register_plugin_stage(meta)
    WorkflowMenuCommand.register_plugin_tool(meta, plugin)

    proc = _fake_processor_with_seed({})

    cmd = WorkflowMenuCommand.__new__(WorkflowMenuCommand)
    cmd.processor = proc
    cmd._dispatch_plugin("smoke")

    assert plugin.launch_count == 1
