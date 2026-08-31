"""A checkpoint is a pause, not a completion.

A step that stops early to wait for a calculation run outside ProPrep has not
produced its artifacts. Recording it "completed" ticked it green, satisfied the
next step's dependency, and let the run walk past it -- so the failure surfaced
later, at the first step that consumed the artifact that was never made.
"""

from proprep.utils.workflow_checklist import WorkflowChecklist, WorkflowStep


class _Executor:
    def __init__(self, result):
        self.result = result

    def h_gate(self):
        return self.result

    def h_next(self):
        return {"summary": "done"}


def _checklist(tmp_path, result):
    steps = [
        WorkflowStep(id="s1", name="QM Gate", description="", handler="h_gate",
                     section="A", dependencies=[], checkpoint=True,
                     checkpoint_message="Run it externally, then resume."),
        WorkflowStep(id="s2", name="Consumer", description="", handler="h_next",
                     section="A", dependencies=["s1"]),
    ]
    cl = WorkflowChecklist(steps=steps, executor=_Executor(result),
                           workflow_name="Test", state_dir=tmp_path)
    cl._initialize_state()
    return cl


def test_a_checkpoint_is_recorded_in_progress(tmp_path):
    cl = _checklist(tmp_path, {"checkpoint": True})
    cl._run_step("s1")
    status = cl.state.get_step_status("s1")
    assert status.status == "in_progress"
    assert status.completed_at is None


def test_a_paused_step_is_what_next_offers(tmp_path):
    """Not the step after it."""
    cl = _checklist(tmp_path, {"checkpoint": True})
    cl._run_step("s1")
    assert cl._get_next_pending_step().id == "s1"


def test_a_paused_step_does_not_satisfy_a_dependency(tmp_path):
    cl = _checklist(tmp_path, {"checkpoint": True})
    cl._run_step("s1")
    assert cl._check_dependencies(cl._get_step_by_id("s2")) == ["s1"]


def test_the_pause_is_announced_rather_than_a_green_tick(tmp_path):
    cl = _checklist(tmp_path, {"checkpoint": True})
    assert cl._run_step("s1") is False
    assert "external calculation" in cl.state.get_step_status("s1").output_summary


def test_a_real_completion_still_completes(tmp_path):
    cl = _checklist(tmp_path, {"summary": "ESP generated"})
    assert cl._run_step("s1") is True
    assert cl.state.get_step_status("s1").status == "completed"
    assert cl._check_dependencies(cl._get_step_by_id("s2")) == []
    assert cl._get_next_pending_step().id == "s2"
