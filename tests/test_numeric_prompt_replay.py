"""
Numeric prompts must be replayed by question, not by position.

Replaying a session into a build that had gained two new integer prompts gave:

    Metal types already used (M positions to skip) (3): [REPLAY: ]
    Ligating-atom types already used (Y positions to skip) (11): [REPLAY: 48]

48 was the Gaussian memory in GB, recorded for an entirely different question.

Cause: IntPrompt and FloatPrompt are NOT subclasses of Prompt — they inherit
PromptBase.ask — so patching Prompt.ask never covered them. They fell through to
the builtin input() interception, where Rich has already printed the question
itself and passes nothing on, so every numeric answer was recorded as
type='input' with an EMPTY prompt string. Matching then had nothing to match on,
and any numeric question would take the next numeric answer in the file.

In the log from that run, all 15 'input' interactions had an empty prompt.
"""

import json

import pytest
from rich.prompt import Confirm, FloatPrompt, IntPrompt, Prompt

from proprep.utils.session_recorder import (
    HybridInterceptor, SessionRecorder, SessionReplayer,
)


def _replayer(tmp_path, interactions):
    path = tmp_path / "session.json"
    path.write_text(json.dumps({
        "interactions": [
            {"type": t, "prompt": p, "response": r, "context": {}}
            for t, p, r in interactions
        ]
    }))
    r = SessionReplayer(str(path))
    r.replaying = True
    return r


def _interceptor(tmp_path, interactions=(), recording=False):
    rec = SessionRecorder()
    rec.recording = recording
    h = HybridInterceptor(rec, _replayer(tmp_path, interactions) if interactions else None)
    # Never touch stdin: the live path is stubbed with a sentinel.
    h._original_int_ask = lambda prompt, **kw: -999
    h._original_float_ask = lambda prompt, **kw: -9.99
    return h, rec


MEMORY_Q = "Memory (GB)"
PROC_Q = "Number of processors"
NEW_Q = "Metal types already used (M positions to skip)"


def test_a_new_numeric_prompt_does_not_steal_another_answer(tmp_path):
    """The reported failure, reproduced."""
    h, _ = _interceptor(tmp_path, [
        ("prompt", MEMORY_Q, "48"),
        ("prompt", PROC_Q, "24"),
    ])

    got = h._intercepted_int_ask(NEW_Q)

    assert got == -999, f"a question absent from the log took answer {got!r}"


def test_recorded_numeric_answers_match_their_own_question(tmp_path):
    h, _ = _interceptor(tmp_path, [
        ("prompt", MEMORY_Q, "48"),
        ("prompt", PROC_Q, "24"),
    ])

    assert h._intercepted_int_ask(MEMORY_Q) == 48
    assert h._intercepted_int_ask(PROC_Q) == 24


def test_the_question_must_match_not_just_be_numeric(tmp_path):
    """Text is checked, so a numeric answer cannot serve a different question.

    Replay is strict about position as well: asking out of order diverges and
    falls through to live input rather than handing over the wrong number.
    Recording the question text is what makes that check possible at all —
    before, every numeric answer was logged with an empty prompt.
    """
    h, _ = _interceptor(tmp_path, [
        ("prompt", MEMORY_Q, "48"),
        ("prompt", PROC_Q, "24"),
    ])

    assert h._intercepted_int_ask(PROC_Q) == -999, "out of order: must not take 48"
    assert h._intercepted_int_ask(MEMORY_Q) == 48, "in order: resynchronises"
    assert h._intercepted_int_ask(PROC_Q) == 24


def test_a_replayed_int_is_an_int_not_a_string(tmp_path):
    h, _ = _interceptor(tmp_path, [("prompt", MEMORY_Q, "48")])

    got = h._intercepted_int_ask(MEMORY_Q)

    assert isinstance(got, int) and got == 48


def test_floats_replay_as_floats(tmp_path):
    h, _ = _interceptor(tmp_path, [("prompt", "pH value", "7.4")])

    got = h._intercepted_float_ask("pH value")

    assert isinstance(got, float) and got == pytest.approx(7.4)


def test_a_non_numeric_recorded_answer_is_not_forced(tmp_path):
    """A text answer under a matching question belongs to a different build."""
    h, _ = _interceptor(tmp_path, [("prompt", MEMORY_Q, "high")])

    got = h._intercepted_int_ask(MEMORY_Q)   # must not raise

    assert got == -999, "should ask rather than coerce a non-number"


def test_numeric_answers_are_recorded_with_their_question(tmp_path):
    """The root cause: these used to be logged with an empty prompt string."""
    h, rec = _interceptor(tmp_path, recording=True)

    h._intercepted_int_ask(MEMORY_Q)

    logged = rec.session_data["interactions"]
    assert len(logged) == 1
    assert logged[0]["prompt"] == MEMORY_Q, "the question text must be recorded"
    assert logged[0]["prompt"] != ""
    assert logged[0]["response"] == "-999"


def test_install_and_uninstall_cover_the_numeric_prompts(tmp_path):
    # NOT the stubbed helper: it overwrites _original_int_ask, which is exactly
    # what uninstall restores from.
    rec = SessionRecorder()
    rec.recording = False
    h = HybridInterceptor(rec, None)
    originals = (Prompt.ask, Confirm.ask, IntPrompt.ask, FloatPrompt.ask)

    h.install()
    try:
        # Bound methods are rebuilt on each attribute access, so compare the
        # underlying function and instance rather than identity.
        assert IntPrompt.ask.__func__ is HybridInterceptor._intercepted_int_ask
        assert IntPrompt.ask.__self__ is h
        assert FloatPrompt.ask.__func__ is HybridInterceptor._intercepted_float_ask
        assert FloatPrompt.ask.__self__ is h
    finally:
        h.uninstall()

    assert (Prompt.ask, Confirm.ask, IntPrompt.ask, FloatPrompt.ask) == originals


def test_intprompt_is_not_a_prompt_subclass():
    """Why patching Prompt.ask was never enough — pin the assumption."""
    assert not issubclass(IntPrompt, Prompt)
    assert not issubclass(FloatPrompt, Prompt)
