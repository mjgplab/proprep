"""
Regression test: optional checklist steps can actually be SKIPPED from the menu.

The step definitions carried optional=True and _skip_step() existed, but no
user-facing command was wired to it, so an optional step (e.g. Structure
Completeness) could not be skipped interactively. A "<num>s" / "s<num>" command
now maps to a skip action for that step.
"""

import proprep.utils.workflow_checklist as wc
from proprep.utils.workflow_checklist import WorkflowChecklist, WorkflowStep


def _checklist():
    steps = [
        WorkflowStep(id="prep-1", name="Filtering", description="",
                     handler="h1", section="Prep", dependencies=[], optional=True),
        WorkflowStep(id="prep-1b", name="Completeness", description="",
                     handler="h2", section="Prep", dependencies=[], optional=True),
        WorkflowStep(id="prep-2", name="Triage", description="",
                     handler="h3", section="Prep", dependencies=[]),
    ]
    return WorkflowChecklist(steps=steps, executor=object(), workflow_name="Test")


def _action_for(monkeypatch, typed):
    cl = _checklist()
    monkeypatch.setattr(wc, "prompt_with_context",
                        lambda *a, **k: typed)
    return cl._get_user_action()


def test_num_s_skips_that_step(monkeypatch):
    # Step 2 in display order is prep-1b (Completeness), which is optional.
    assert _action_for(monkeypatch, "2s") == "skip:prep-1b"
    assert _action_for(monkeypatch, "s2") == "skip:prep-1b"


def test_bare_s_is_still_save(monkeypatch):
    assert _action_for(monkeypatch, "s") == "save"


def test_num_still_runs_and_num_i_still_info(monkeypatch):
    assert _action_for(monkeypatch, "2") == "run:prep-1b"
    assert _action_for(monkeypatch, "2i") == "info:prep-1b"


def test_skip_only_marks_optional_steps():
    # _skip_step refuses a non-optional step (prep-2 / Triage).
    cl = _checklist()
    cl._initialize_state(pdb_file=None)
    assert cl._skip_step("prep-2") is False
    assert cl._skip_step("prep-1b") is True
    assert cl.state.get_step_status("prep-1b").status == "skipped"
