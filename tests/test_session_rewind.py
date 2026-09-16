"""
Typing ``undo`` at any prompt rewinds the session in-process.

The manual recovery used to be: exit, relaunch, pick the log, find and edit
the wrong answer, replay. The rewind keyword does the same from inside the
program. These tests cover the four seams it touches: Rich's validation loop
(a prompt with ``choices`` must accept the word), the interceptor (the picker
must never be recorded), the signal (it must pass ``except Exception``), and
``main``'s rebuild after the unwind.
"""

import io
import json
import os
import types

import pytest
from rich.console import Console
from rich.prompt import Confirm, IntPrompt, Prompt

from proprep.utils import session_rewind as rw
from proprep.utils.session_recorder import (
    HybridInterceptor, InterceptedPrompt, SessionManager, SessionRecorder,
)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _quiet_console():
    return Console(file=io.StringIO(), force_terminal=False, width=100)


def _recorder(tmp_path, n=4):
    """A recording SessionRecorder holding n answers with rich context."""
    rec = SessionRecorder(str(tmp_path / "s.json"))
    rec.recording = True
    for i in range(n):
        rec.session_data["interactions"].append({
            "index": i, "type": "prompt", "prompt": f"Q{i}", "response": str(i % 2 + 1),
            "choices": ["1", "2"],
            "context": {"module": "M", "description": f"question {i}",
                        "options_map": {"1": "one", "2": "two"}, "option_label": "one"},
        })
    rec._interaction_count = n
    return rec


class _Scripted:
    """Stands in for the saved original Rich ``ask``.

    Rich raises RewindKeyword from inside its loop when the keyword is typed;
    the stub does the same so the interceptor path is exercised end to end.
    """

    def __init__(self, answers):
        self.answers = list(answers)
        self.asked = []

    def __call__(self, prompt, **kwargs):
        self.asked.append(prompt)
        answer = self.answers.pop(0)
        if rw.is_rewind_keyword(answer):
            raise rw.RewindKeyword()
        return answer


def _interceptor(rec, answers, cls=InterceptedPrompt):
    icp = cls(recorder=rec)
    icp._original_prompt_ask = _Scripted(answers)
    return icp


# ---------------------------------------------------------------------------
# Rich validation loop
# ---------------------------------------------------------------------------

def test_keyword_is_accepted_through_a_choices_prompt():
    """Rich re-asks on an invalid choice; the word must get past that."""
    with pytest.raises(rw.RewindKeyword):
        rw._call_with_class(Prompt.ask, rw.RewindPrompt, "Q", choices=["1", "2"],
                            console=_quiet_console(), stream=io.StringIO("undo\n"))


def test_keyword_is_accepted_by_confirm_and_int_prompts():
    with pytest.raises(rw.RewindKeyword):
        rw._call_with_class(Confirm.ask, rw.RewindConfirm, "Sure?",
                            console=_quiet_console(), stream=io.StringIO("UNDO\n"))
    with pytest.raises(rw.RewindKeyword):
        rw._call_with_class(IntPrompt.ask, rw.RewindIntPrompt, "N",
                            console=_quiet_console(), stream=io.StringIO(" undo \n"))


def test_ordinary_validation_is_unchanged():
    value = rw._call_with_class(Prompt.ask, rw.RewindPrompt, "Q", choices=["1", "2"],
                                console=_quiet_console(), stream=io.StringIO("9\n2\n"))
    assert value == "2"


def test_a_stubbed_original_is_called_directly():
    """Existing tests replace the originals with lambdas; keep that working."""
    assert rw._call_with_class(lambda p, **kw: "stub", rw.RewindPrompt, "Q") == "stub"


# ---------------------------------------------------------------------------
# signal
# ---------------------------------------------------------------------------

def test_signal_passes_through_except_exception():
    with pytest.raises(rw.RewindRequested):
        try:
            raise rw.RewindRequested(3, "reanswer", "s.json")
        except Exception:  # what most of ProPrep's handlers say
            pytest.fail("rewind signal was swallowed by except Exception")


def test_hybrid_args_for_both_modes():
    assert rw.RewindRequested(43, "change", "s.json", "-2").hybrid_args() == {
        "truncate_at": 43, "keep_following": True, "new_value": "-2"}
    assert rw.RewindRequested(43, "reanswer", "s.json").hybrid_args() == {
        "truncate_at": 42, "keep_following": False, "new_value": None}
    # Re-answering the very first question replays nothing.
    assert rw.RewindRequested(0, "reanswer", "s.json").hybrid_args()["truncate_at"] == -1


# ---------------------------------------------------------------------------
# picker through the interceptor
# ---------------------------------------------------------------------------

def test_change_mode_raises_with_validated_value_and_records_nothing(tmp_path):
    rec = _recorder(tmp_path)
    before = json.dumps(rec.session_data)
    # keyword, pick #2, change it, bad value rejected, then '2' accepted
    icp = _interceptor(rec, ["undo", "2", "c", "zzz", "2"])
    with pytest.raises(rw.RewindRequested) as info:
        icp._intercepted_prompt_ask("Live question", choices=["a", "b"])
    r = info.value
    assert (r.index, r.mode, r.new_value) == (2, "change", "2")
    assert r.session_file == str(tmp_path / "s.json")
    assert json.dumps(rec.session_data) == before, "picker prompts were recorded"
    assert icp._in_rich_prompt is False


def test_reanswer_mode_defaults_to_the_last_answer(tmp_path):
    rec = _recorder(tmp_path, n=5)
    icp = _interceptor(rec, ["undo", "4", "r"])
    with pytest.raises(rw.RewindRequested) as info:
        icp._intercepted_prompt_ask("Live question")
    assert (info.value.index, info.value.mode) == (4, "reanswer")


def test_confirm_answer_is_normalised_to_yes_no(tmp_path):
    rec = _recorder(tmp_path)
    rec.session_data["interactions"][1].update({"type": "confirm", "choices": None})
    icp = _interceptor(rec, ["undo", "1", "c", "maybe", "Y"])
    with pytest.raises(rw.RewindRequested) as info:
        icp._intercepted_prompt_ask("Live question")
    assert info.value.new_value == "yes"


def test_cancel_returns_to_the_question_and_records_the_real_answer(tmp_path):
    rec = _recorder(tmp_path)
    icp = _interceptor(rec, ["undo", "x", "b"])
    assert icp._intercepted_prompt_ask("Live question", choices=["a", "b"]) == "b"
    assert rec.session_data["interactions"][-1]["response"] == "b"
    assert len(rec.session_data["interactions"]) == 5


def test_more_pages_back_then_index_is_accepted(tmp_path):
    rec = _recorder(tmp_path, n=30)
    icp = _interceptor(rec, ["undo", "m", "3", "r"])
    with pytest.raises(rw.RewindRequested) as info:
        icp._intercepted_prompt_ask("Live question")
    assert info.value.index == 3


def test_out_of_range_index_is_re_asked(tmp_path):
    rec = _recorder(tmp_path)
    icp = _interceptor(rec, ["undo", "99", "abc", "1", "r"])
    with pytest.raises(rw.RewindRequested) as info:
        icp._intercepted_prompt_ask("Live question")
    assert info.value.index == 1


def test_without_a_recorder_the_keyword_explains_and_re_asks():
    icp = _interceptor(None, ["undo", "a"])
    assert icp._intercepted_prompt_ask("Live question", choices=["a"]) == "a"


def test_nothing_recorded_yet_re_asks(tmp_path):
    rec = _recorder(tmp_path, n=0)
    icp = _interceptor(rec, ["undo", "a"])
    assert icp._intercepted_prompt_ask("Live question") == "a"
    assert len(rec.session_data["interactions"]) == 1


def test_hybrid_interceptor_takes_the_same_path(tmp_path):
    rec = _recorder(tmp_path)
    icp = _interceptor(rec, ["undo", "0", "r"], cls=HybridInterceptor)
    with pytest.raises(rw.RewindRequested) as info:
        icp._intercepted_prompt_ask("Live question")
    assert info.value.hybrid_args()["truncate_at"] == -1


def test_builtin_input_path_honours_the_keyword(tmp_path):
    rec = _recorder(tmp_path)
    icp = InterceptedPrompt(recorder=rec)
    icp._original_prompt_ask = _Scripted(["2", "r"])   # picker answers
    icp._original_input = lambda prompt="": "undo"       # builtin input returns the text
    with pytest.raises(rw.RewindRequested) as info:
        icp._intercepted_input("raw> ")
    assert info.value.index == 2


# ---------------------------------------------------------------------------
# hybrid mode with the -1 sentinel
# ---------------------------------------------------------------------------

def test_hybrid_mode_truncate_minus_one_replays_nothing(tmp_path, monkeypatch):
    log = tmp_path / "s.json"
    log.write_text(json.dumps({"interactions": [
        {"index": 0, "type": "prompt", "prompt": "Q0", "response": "1", "context": {}},
        {"index": 1, "type": "prompt", "prompt": "Q1", "response": "2", "context": {}},
    ], "metadata": {}}))
    monkeypatch.chdir(tmp_path)
    mgr = SessionManager()
    try:
        mgr.start_hybrid_mode(str(log), truncate_at=-1, keep_following_interactions=False)
        assert mgr.replayer.session_data["interactions"] == []
        assert mgr.recorder.session_data["interactions"] == []
        assert list(tmp_path.glob("s.json.bak-*")), "log was not backed up"
    finally:
        mgr.stop()


# ---------------------------------------------------------------------------
# main: rebuild after the unwind
# ---------------------------------------------------------------------------

class _Workspace(dict):
    def set(self, k, v):
        self[k] = v


class _FakeProcessor:
    instances = []

    def __init__(self):
        self.workspace = _Workspace()
        self.cleaned = False
        self.hybrid_calls = []
        _FakeProcessor.instances.append(self)

    def cleanup(self):
        self.cleaned = True

    def start_session_hybrid(self, filename, **kw):
        self.hybrid_calls.append((filename, kw))

    def get_module_instance(self, name):
        raise AssertionError("no structure load expected without pdbid/pdbfile")

    def _get_workspace(self):
        return self.workspace


def test_rebuild_for_rewind_replaces_processor_and_replays(tmp_path, monkeypatch):
    from proprep import main as m
    _FakeProcessor.instances.clear()
    monkeypatch.setattr(m, "PDBProcessor", _FakeProcessor)
    monkeypatch.setattr(m, "integrate_session_manager", lambda p: None)
    monkeypatch.setattr(m, "reload_structure_from_session_metadata", lambda p, f: False)
    elsewhere = tmp_path / "elsewhere"; elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    project = str(tmp_path)
    args = types.SimpleNamespace(menu_grid=True, menu_list=False, analysis=False,
                                 pdbview=None, debug=True, pdbid=None, pdbfile=None)
    old = _FakeProcessor()
    rewind = rw.RewindRequested(7, "change", str(tmp_path / "s.json"), "x")

    new = m._rebuild_for_rewind(old, rewind, args, project, "workflow")

    assert old.cleaned and new is not old
    assert os.getcwd() == os.path.realpath(project)
    assert new.workspace["menu_mode"] == "workflow"
    assert new.workspace["menu_layout"] == "grid"
    assert new.workspace["project_directory"] == project
    assert new.workspace["debug"] is True
    assert new.hybrid_calls == [(str(tmp_path / "s.json"),
                                 {"truncate_at": 7, "keep_following": True, "new_value": "x"})]
