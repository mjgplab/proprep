"""Tests for plugin discovery via importlib entry-points.

We don't actually install fake packages — instead we monkey-patch
``importlib.metadata.entry_points`` to return mock EntryPoint
objects. This exercises the discovery filter logic (api_version
mismatch, is_available=False, exceptions during get_metadata) in
isolation.
"""

from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import patch

import pytest

from proprep.plugins import PLUGIN_API_VERSION, PluginMetadata
from proprep.plugins.discovery import discover_plugins


def _fake_entry_point(name: str, plugin_obj):
    """An object that quacks like ``importlib.metadata.EntryPoint``."""
    return SimpleNamespace(name=name, load=lambda: plugin_obj)


class _GoodPlugin:
    def get_metadata(self):
        return PluginMetadata(
            name="good", version="1.0", description="d",
            stage="x", stage_name="X",
        )

    def is_available(self):
        return True

    def launch(self, host):
        pass


class _UnavailablePlugin:
    def get_metadata(self):
        return PluginMetadata(
            name="unavail", version="1.0", description="d",
            stage="x", stage_name="X",
        )

    def is_available(self):
        return False

    def launch(self, host):
        pass


class _StaleApiPlugin:
    def get_metadata(self):
        return PluginMetadata(
            name="stale", version="1.0", description="d",
            stage="x", stage_name="X",
            api_version=PLUGIN_API_VERSION + 99,
        )

    def is_available(self):
        return True

    def launch(self, host):
        pass


class _CrashesOnMetadata:
    def get_metadata(self):
        raise RuntimeError("synthetic")

    def is_available(self):
        return True

    def launch(self, host):
        pass


def _patch_entry_points(plugins):
    """Patch entry_points to return the given list of fake EPs."""
    return patch(
        "proprep.plugins.discovery.entry_points",
        return_value=plugins,
    )


def test_discovery_empty_when_no_plugins_registered():
    """Fresh ProPrep install with nothing in the entry-point group →
    discover_plugins returns {} cleanly. Never raises."""
    with _patch_entry_points([]):
        assert discover_plugins() == {}


def test_discovery_loads_good_plugin():
    plugin = _GoodPlugin()
    ep = _fake_entry_point("good", plugin)
    with _patch_entry_points([ep]):
        result = discover_plugins()
    assert "good" in result
    assert result["good"] is plugin


def test_discovery_drops_unavailable_plugin():
    """is_available()=False plugins are silently dropped — no
    confusing menu entries that would crash on launch."""
    ep = _fake_entry_point("unavail", _UnavailablePlugin())
    with _patch_entry_points([ep]):
        result = discover_plugins()
    assert result == {}


def test_discovery_drops_api_version_mismatch():
    """Plugins claiming a different PLUGIN_API_VERSION are dropped
    so they can't crash on contract drift."""
    ep = _fake_entry_point("stale", _StaleApiPlugin())
    with _patch_entry_points([ep]):
        result = discover_plugins()
    assert result == {}


def test_discovery_isolates_metadata_exceptions():
    """A plugin that raises in get_metadata is dropped; other
    plugins still load."""
    bad = _fake_entry_point("bad", _CrashesOnMetadata())
    good = _fake_entry_point("good", _GoodPlugin())
    with _patch_entry_points([bad, good]):
        result = discover_plugins()
    assert "good" in result
    assert "bad" not in result


def test_discovery_isolates_load_failure():
    """Even a broken entry-point (load() raises) doesn't crash
    discovery."""
    broken = SimpleNamespace(
        name="broken",
        load=lambda: (_ for _ in ()).throw(ImportError("synthetic")),
    )
    good = _fake_entry_point("good", _GoodPlugin())
    with _patch_entry_points([broken, good]):
        result = discover_plugins()
    assert "good" in result
    assert "broken" not in result


def test_discovery_handles_entry_points_raising():
    """importlib.metadata can raise on certain mis-installed packages.
    discover_plugins must swallow and return {}."""
    with patch(
        "proprep.plugins.discovery.entry_points",
        side_effect=RuntimeError("packaging metadata broken"),
    ):
        assert discover_plugins() == {}
