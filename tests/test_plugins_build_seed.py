"""Tests for PDBProcessor.build_plugin_seed.

This is the centralised host→plugin seed-builder. It reads
``meta.consumes_workspace_keys`` and copies any matching keys from
the host workspace into the seed dict. Absent keys are silently
omitted (not present as None) so plugins disambiguate "host doesn't
have it" from "host has it and it's null".

We don't instantiate the full PDBProcessor (which imports half of
ProPrep) — instead we exercise the method on a minimal stand-in
that exposes a workspace with .has() / .get(). The method body uses
no other PDBProcessor state.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from proprep.application.pdbprocessor import PDBProcessor
from proprep.plugins import PluginMetadata


def _fake_workspace(data: dict):
    """Minimal workspace stub with the .has/.get surface
    build_plugin_seed depends on."""
    return SimpleNamespace(
        has=lambda key: key in data,
        get=lambda key, default=None: data.get(key, default),
    )


def _bound_build_plugin_seed(ws):
    """Bind the unbound method to a stub processor so we can call
    it without running PDBProcessor.__init__."""
    stub = SimpleNamespace(workspace=ws)
    return lambda meta: PDBProcessor.build_plugin_seed(stub, meta)


def _meta(keys):
    return PluginMetadata(
        name="x", version="1.0", description="d",
        stage="x", stage_name="X",
        consumes_workspace_keys=keys,
    )


def test_seed_copies_present_workspace_keys():
    """The keys a plugin declares + the host has → land in the seed
    with values copied verbatim."""
    ws = _fake_workspace({"foo": [1, 2, 3], "bar": "value"})
    build = _bound_build_plugin_seed(ws)
    seed = build(_meta(["foo", "bar"]))
    assert seed == {"foo": [1, 2, 3], "bar": "value"}


def test_seed_omits_keys_not_in_workspace():
    """Declared-but-absent keys are silently omitted — never present
    as None. Plugins use ``key in seed`` rather than
    ``seed.get(key) is None`` to disambiguate."""
    ws = _fake_workspace({"foo": "value"})
    build = _bound_build_plugin_seed(ws)
    seed = build(_meta(["foo", "missing"]))
    assert seed == {"foo": "value"}
    assert "missing" not in seed


def test_seed_ignores_workspace_keys_not_declared():
    """The host might have many keys; only those the plugin
    explicitly declares come through. Keeps cross-package coupling
    explicit and auditable."""
    ws = _fake_workspace({"declared": 1, "secret": 2, "private": 3})
    build = _bound_build_plugin_seed(ws)
    seed = build(_meta(["declared"]))
    assert seed == {"declared": 1}


def test_seed_empty_when_consumes_list_empty():
    """A plugin with no declared workspace keys gets an empty seed
    — only HostContext fields (working_dir etc.) carry the
    handoff."""
    ws = _fake_workspace({"foo": 1})
    build = _bound_build_plugin_seed(ws)
    seed = build(_meta([]))
    assert seed == {}


def test_seed_preserves_value_identity_not_copies():
    """Pass-through, not deep-copy. The plugin gets a reference to
    the host's value (cheap; documented). Plugins must not mutate
    seed values in place — the contract is read-only consumption."""
    sentinel = [1, 2, 3]
    ws = _fake_workspace({"foo": sentinel})
    build = _bound_build_plugin_seed(ws)
    seed = build(_meta(["foo"]))
    assert seed["foo"] is sentinel
