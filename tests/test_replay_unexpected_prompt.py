"""
An unexpected prompt must not destroy a replay.

Two defects made a leftover ``workflow_state.json`` block replay of a session
log:

1. ``WorkflowChecklist.run`` offered "Resume from saved state?" whenever the
   file existed. The recorded run began with no state file, so its log has no
   answer for that question.

2. ``SessionReplayer.get_next_response`` scanned FORWARD past non-matching
   interactions and, on a miss, left the position at the END of the log. So the
   one unanswerable prompt consumed every remaining recorded answer: replay
   reported itself exhausted and every later prompt fell through to live input,
   which stalls an unattended batch replay.

The forward scan itself is deliberate — it lets replay tolerate a recorded
prompt the current run does not ask — so the fix rewinds only when nothing
matched.
"""

import json

import pytest

from proprep.utils.session_recorder import SessionReplayer
from proprep.utils.workflow_checklist import WorkflowChecklist, WorkflowStep


def _write_log(path, interactions):
    path.write_text(json.dumps({
        "session_name": "test",
        "interactions": [
            {"type": t, "prompt": p, "response": r, "context": {}}
            for t, p, r in interactions
        ],
    }))
    return str(path)


def _replayer(tmp_path, interactions):
    r = SessionReplayer(_write_log(tmp_path / "session.json", interactions))
    r.replaying = True
    return r


# --------------------------------------------------------------------------- #
# SessionReplayer.get_next_response
# --------------------------------------------------------------------------- #

def test_unexpected_prompt_does_not_consume_the_log(tmp_path):
    """The core defect: one unmatched prompt used to eat every later answer."""
    r = _replayer(tmp_path, [
        ("prompt", "Pick a structure", "3"),
        ("confirm", "Proceed?", "yes"),
        ("prompt", "Total charge of small model", "-2"),
    ])

    # A prompt the recording never saw (the stale-state resume offer).
    assert r.get_next_response("confirm", "Resume from saved state?") is None

    # Every recorded answer is still available, in order.
    assert r.get_next_response("prompt", "Pick a structure") == "3"
    assert r.get_next_response("confirm", "Proceed?") == "yes"
    assert r.get_next_response("prompt", "Total charge of small model") == "-2"


def test_replay_is_not_reported_exhausted_after_a_miss(tmp_path):
    """has_more_interactions drives the callers' '[Replay complete]' latch."""
    r = _replayer(tmp_path, [
        ("prompt", "Pick a structure", "3"),
        ("confirm", "Proceed?", "yes"),
    ])

    r.get_next_response("confirm", "Resume from saved state?")

    assert r.has_more_interactions(), \
        "an unexpected prompt marked replay exhausted with answers still owed"


def test_a_recorded_prompt_this_run_skips_stalls_rather_than_leaping(tmp_path):
    """Matching is STRICT: only the next unconsumed interaction can answer.

    The forward scan this replaces could not tell "a recorded prompt the run
    does not ask" apart from "the same question asked at a different point in
    the workflow", and leapt for the latter — consuming the checklist decisions
    in between. Stalling is the honest outcome: the prompt falls through to live
    input, and the position is kept so replay can resynchronise.
    """
    r = _replayer(tmp_path, [
        ("confirm", "Overwrite existing files?", "yes"),   # not asked this run
        ("prompt", "Pick a structure", "3"),
    ])

    assert r.get_next_response("prompt", "Pick a structure") is None
    assert r.interaction_index == 0, "the position must be kept for resync"
    assert r.has_more_interactions()


def test_replay_resynchronises_once_the_recorded_question_returns(tmp_path):
    """Divergence is recoverable: answer live, then replay picks up again."""
    r = _replayer(tmp_path, [
        ("prompt", "Select action", "n"),
        ("prompt", "Select action", "11"),
    ])

    # A prompt the log never saw — no answer, no movement.
    assert r.get_next_response("confirm", "Add hydrogen(s) to FES?") is None
    assert r.interaction_index == 0

    # The recorded question comes round: replay resumes exactly where it was.
    assert r.get_next_response("prompt", "Select action") == "n"
    assert r.get_next_response("prompt", "Select action") == "11"


def test_repeated_misses_do_not_drift(tmp_path):
    """Several unexpected prompts in a row must each cost nothing."""
    r = _replayer(tmp_path, [("prompt", "Pick a structure", "3")])

    for _ in range(5):
        assert r.get_next_response("confirm", "Some new question?") is None

    assert r.get_next_response("prompt", "Pick a structure") == "3"


def test_type_mismatch_alone_is_still_a_miss(tmp_path):
    """Same text, wrong interaction type: a miss, and non-destructive."""
    r = _replayer(tmp_path, [("prompt", "Proceed?", "y")])

    assert r.get_next_response("confirm", "Proceed?") is None
    assert r.get_next_response("prompt", "Proceed?") == "y"


# --------------------------------------------------------------------------- #
# WorkflowChecklist resume offer
# --------------------------------------------------------------------------- #

class _Manager:
    def __init__(self, replaying):
        self._replaying = replaying

    def is_replaying(self):
        return self._replaying


class _Processor:
    def __init__(self, replaying):
        self.session_manager = _Manager(replaying)

    def _get_workspace(self):
        return None


def _checklist(tmp_path, replaying, monkeypatch):
    steps = [WorkflowStep(id="s1", name="Step 1", description="d",
                          handler="h", section="S")]
    cl = WorkflowChecklist(steps=steps, executor=object(),
                           processor=_Processor(replaying),
                           state_dir=tmp_path)
    # A state file left behind by an earlier run.
    (tmp_path / WorkflowChecklist.STATE_FILENAME).write_text(json.dumps({
        "workflow_id": "abcd1234",
        "workflow_name": "Workflow",
        "working_directory": str(tmp_path),
        "created_at": "2026-08-14T10:00:00",
        "updated_at": "2026-08-14T10:05:00",
        "step_statuses": {},
        "workspace_snapshot": {},
        "pdb_file": None,
    }))

    offered = []
    monkeypatch.setattr(WorkflowChecklist, "_offer_resume",
                        lambda self, f: offered.append(f) or True)
    monkeypatch.setattr(WorkflowChecklist, "_run_loop", lambda self: True)
    return cl, offered


def test_replay_skips_the_resume_offer(tmp_path, monkeypatch):
    cl, offered = _checklist(tmp_path, replaying=True, monkeypatch=monkeypatch)

    cl.run()

    assert offered == [], "the resume prompt was asked during replay"
    assert cl.state is not None, "a fresh state should have been initialised"


def test_interactive_run_still_offers_resume(tmp_path, monkeypatch):
    cl, offered = _checklist(tmp_path, replaying=False, monkeypatch=monkeypatch)

    cl.run()

    assert len(offered) == 1, "an interactive user lost the resume offer"


def test_no_session_manager_is_treated_as_live(tmp_path, monkeypatch):
    """A processor without a session manager must keep the resume offer."""
    cl, offered = _checklist(tmp_path, replaying=False, monkeypatch=monkeypatch)
    del cl.processor.session_manager

    cl.run()

    assert len(offered) == 1
