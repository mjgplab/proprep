"""Structural-typing tests for the plugin Protocols.

The two Protocols are runtime_checkable, so plugins (and host
recorders) can be validated by isinstance(). These tests pin the
shapes both packages have to satisfy. They're cheap, deterministic,
and don't touch ProPrep's heavy modules.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from proprep.plugins import (
    HostContext,
    PLUGIN_API_VERSION,
    PluginMetadata,
    ProPrepPlugin,
    SessionRecorder,
)


# ---------------------------------------------------------------------------
# SessionRecorder Protocol
# ---------------------------------------------------------------------------


class _MinimalRecorder:
    """Bare structural match for SessionRecorder. Mirrors the shape
    both ProPrep's and etanalyze's existing classes already provide."""

    def __init__(self) -> None:
        self.recording: bool = False
        self.record_file: Optional[str] = None
        self.session_data: Dict[str, Any] = {}
        self.calls: List[Dict[str, Any]] = []

    def record_interaction(
        self, interaction_type: str, prompt: str, response: str,
        choices: Optional[List[str]] = None,
        context: Optional[Dict[str, Any]] = None,
    ) -> None:
        self.calls.append({
            "type": interaction_type, "prompt": prompt, "response": response,
            "choices": choices, "context": context,
        })

    def start_recording(self, metadata: Optional[Dict[str, Any]] = None) -> None:
        self.recording = True

    def stop_recording(self) -> None:
        self.recording = False


def test_minimal_recorder_satisfies_protocol():
    """Any class with the right methods/attributes counts as a
    SessionRecorder — no inheritance required. Locks in the
    structural-typing contract."""
    rec = _MinimalRecorder()
    assert isinstance(rec, SessionRecorder)


def test_object_missing_methods_fails_protocol_check():
    """Sanity: an unrelated object isn't accidentally a recorder."""
    assert not isinstance(object(), SessionRecorder)


# ---------------------------------------------------------------------------
# ProPrepPlugin Protocol
# ---------------------------------------------------------------------------


class _MinimalPlugin:
    def get_metadata(self) -> PluginMetadata:
        return PluginMetadata(
            name="x", version="0.0", description="d",
            stage="x", stage_name="X",
        )

    def is_available(self) -> bool:
        return True

    def launch(self, host: HostContext) -> None:
        pass


def test_minimal_plugin_satisfies_protocol():
    assert isinstance(_MinimalPlugin(), ProPrepPlugin)


def test_object_missing_launch_fails_protocol_check():
    class _NoLaunch:
        def get_metadata(self): pass
        def is_available(self): return True

    assert not isinstance(_NoLaunch(), ProPrepPlugin)


# ---------------------------------------------------------------------------
# HostContext + PluginMetadata defaults
# ---------------------------------------------------------------------------


def test_host_context_defaults_match_contract():
    """Working_dir is required; everything else is optional. Empty
    seed_state is the documented "no host data yet" state — plugins
    must tolerate it."""
    from pathlib import Path
    h = HostContext(working_dir=Path("/tmp"))
    assert h.seed_state == {}
    assert h.session_recorder is None
    assert h.console is None
    assert h.logger is None
    assert h.api_version == PLUGIN_API_VERSION


def test_plugin_metadata_consumes_keys_defaults_to_empty_list():
    """A plugin that wants no host data still parses; the seed
    builder loops over an empty list and returns {}."""
    m = PluginMetadata(
        name="x", version="0.0", description="d",
        stage="x", stage_name="X",
    )
    assert m.consumes_workspace_keys == []
    assert m.api_version == PLUGIN_API_VERSION
    assert m.stage_order == 999
    assert m.tool_order == 0
    # Optional UX fields default to None — host renderer falls back
    # to the lowercase ``name`` and emits no nav-shortcut entry.
    assert m.display_name is None
    assert m.stage_shortcut is None


def test_plugin_metadata_accepts_display_name_and_shortcut():
    """Optional fields wire through unchanged — used by the menu
    renderer (display_name) and the nav footer (stage_shortcut)."""
    m = PluginMetadata(
        name="etanalyze", version="1.0", description="d",
        stage="analyze", stage_name="Analyze",
        display_name="ETAnalyze", stage_shortcut="z",
    )
    assert m.display_name == "ETAnalyze"
    assert m.stage_shortcut == "z"
