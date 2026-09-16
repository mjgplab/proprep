"""BLAST runs without a loaded structure; the sequence-source prompt keeps a
constant choice set (session-replay invariant) and re-prompts when the
structure route is picked with nothing loaded."""

import os

import pytest
from rich.console import Console

import proprep.structure_prep.homology_searcher as hs
from proprep.structure_prep.homology_searcher import BLASTIntegrationModule
from proprep.utils.enhanced_menu import OptionStatus


class _Workspace(dict):
    def has(self, k):
        return k in self

    def set(self, k, v):
        self[k] = v


class _Processor:
    def __init__(self):
        self.console = Console(file=open(os.devnull, "w"), force_terminal=False)
        self.workspace = _Workspace()

    def _get_workspace(self):
        return self.workspace


def _module(workspace=None):
    m = BLASTIntegrationModule()
    m.processor = _Processor()
    if workspace is not None:
        m.processor.workspace = workspace
    m.initialize()
    return m


def _scripted(answers):
    """Return a prompt_with_context stand-in that pops answers in order and
    records the choice sets it was offered."""
    seen = []

    def fake(**kwargs):
        seen.append(kwargs.get("choices"))
        return answers.pop(0)

    fake.seen = seen
    return fake


# ---------------------------------------------------------------- menu ----


def test_run_blast_option_available_without_structure():
    m = _module()
    opts = {o.key: o for o in m.get_enhanced_menu_options(_Workspace())}
    assert opts["1"].status == OptionStatus.AVAILABLE
    assert opts["1"].dependency_text == ""
    # MODELLER build still needs a template structure.
    assert opts["5"].status == OptionStatus.BLOCKED


def test_module_never_gated_on_structure():
    m = _module()
    ws = _Workspace()
    assert m.can_process(ws) is True
    assert m.availability_note(ws) is None
    assert m.get_workspace_requirements() == []


# --------------------------------------------------------- choice set ----


def test_sequence_prompt_choices_constant_with_and_without_structure(monkeypatch):
    m = _module()

    fake = _scripted(["1", "MKV"])
    monkeypatch.setattr(hs, "prompt_with_context", fake)
    assert m._get_sequence_input(None, {}) == "MKV"
    without = fake.seen[0]

    fake = _scripted(["1", "MKV"])
    monkeypatch.setattr(hs, "prompt_with_context", fake)
    assert m._get_sequence_input("1ABC", {"A": "MKVL"}) == "MKV"
    with_structure = fake.seen[0]

    assert without == with_structure == ["1", "2", "3"]


def test_structure_route_reprompts_when_nothing_loaded(monkeypatch):
    m = _module()
    fake = _scripted(["3", "1", "MKVLA"])
    monkeypatch.setattr(hs, "prompt_with_context", fake)
    assert m._get_sequence_input(None, {}) == "MKVLA"
    # Offered the method prompt twice (3 rejected, then 1), then the sequence prompt.
    assert fake.seen == [["1", "2", "3"], ["1", "2", "3"], None]


def test_structure_route_uses_loaded_chain(monkeypatch):
    m = _module()
    fake = _scripted(["3", "B"])
    monkeypatch.setattr(hs, "prompt_with_context", fake)
    assert m._get_sequence_input("1ABC", {"A": "AAAA", "B": "BBBB"}) == "BBBB"
    assert fake.seen == [["1", "2", "3"], ["A", "B"]]


# ------------------------------------------------------- entry point ----


def test_run_blast_search_proceeds_without_structure(monkeypatch):
    m = _module()
    calls = {}

    def fake_interactive(structure_id, chain_sequences):
        calls["args"] = (structure_id, chain_sequences)
        return False

    monkeypatch.setattr(m, "_run_interactive_blast", fake_interactive)
    ws = _Workspace()
    m.run_blast_search(ws, interactive=True)
    assert calls["args"] == (None, {})


def test_automated_blast_without_structure_is_refused(monkeypatch):
    m = _module()
    monkeypatch.setattr(
        m, "_run_automated_blast", lambda *a: pytest.fail("should not run")
    )
    ws = _Workspace()
    assert m.run_blast_search(ws, interactive=False) is ws
