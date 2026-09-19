"""A rejected answer is named. A key pressed while ProPrep was printing (a ' next to Enter)
waits in the terminal's buffer and turns the next "4" into "'4"; Rich's "Please select one of
the available options" then reads as "option 4 is unavailable"."""

import io
import json
import sys
import types

from rich.console import Console

import proprep.utils.prompts as prompts
from proprep.utils.prompts import confirm_with_context, prompt_with_context

CHOICES = ["1", "2", "3", "4", "5", "6", "b", "s", "m", "x"]


def _ask(monkeypatch, typed, ask, processor=None):
    monkeypatch.setattr(sys, "stdin", io.StringIO(typed))
    console = Console(width=120, record=True, file=io.StringIO(), force_terminal=False)
    monkeypatch.setattr(prompts, "get_console", lambda: console)
    import rich.prompt
    monkeypatch.setattr(rich.prompt, "get_console", lambda: console)
    return ask(processor), console.export_text()


def test_a_stray_character_before_the_answer_is_shown(monkeypatch):
    answer, out = _ask(monkeypatch, "'4\n4\n", lambda p: prompt_with_context(
        p, "Enter your choice", choices=CHOICES, default="1", options_map={"4": "Generate prmtop/rst7"}))
    assert answer == "4"
    assert '''Received "'4", which is not one of the available options''' in out
    assert "Please select one of" not in out


def test_markup_characters_in_the_input_do_not_break_the_message(monkeypatch):
    answer, out = _ask(monkeypatch, "[/bold]4\n\x1b[A\n2\n", lambda p: prompt_with_context(
        p, "Enter your choice", choices=CHOICES))
    assert answer == "2"
    assert "Received '[/bold]4'" in out and r"Received '\x1b[A'" in out       # an arrow key, made visible


def test_a_valid_answer_in_either_case_or_with_spaces_is_not_reported(monkeypatch):
    answer, out = _ask(monkeypatch, " B \n", lambda p: prompt_with_context(p, "Enter your choice", choices=CHOICES))
    assert answer == "b" and "Received" not in out


def test_confirm_names_what_it_received(monkeypatch):
    answer, out = _ask(monkeypatch, "'y\ny\n", lambda p: confirm_with_context(p, "View the tLEaP input file?", default=True))
    assert answer is True
    assert '''Received "'y"; please enter y or n''' in out


def test_the_recording_is_unchanged_only_the_accepted_answer_is_stored(monkeypatch, tmp_path):
    from proprep.utils.session_recorder import SessionManager
    manager = SessionManager()
    log = str(tmp_path / "session.json")
    manager.start_recording(log, {})
    try:
        answer, out = _ask(monkeypatch, "'4\n4\n", lambda p: prompt_with_context(
            p, "Enter your choice", choices=CHOICES, module="Topology Generator Menu",
            description="Select option", options_map={"4": "Generate prmtop/rst7"}),
            processor=types.SimpleNamespace(session_manager=manager))
    finally:
        manager.stop()
    interactions = json.load(open(log))["interactions"]
    assert answer == "4" and "Received" in out
    assert [i["response"] for i in interactions] == ["4"]
    assert sorted(interactions[0]) == ["choices", "context", "index", "prompt", "response", "timestamp", "type"]
