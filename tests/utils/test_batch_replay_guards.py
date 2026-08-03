"""Regression tests for the two batch-replay guards.

Batch replay used to have no way to fail. A prompt the recording could not
answer returned None from ``SessionReplayer.get_next_response``, which the
interceptor treats as a cue to ask the user; but batch redirects stdout into
``run.log`` and leaves stdin on the real terminal, so the question went where
nobody was looking and the run blocked forever behind a frozen progress bar.
Because the failed search also consumes the rest of the recording, every later
prompt in that run blocked too.

Two guards, tested here. Before the batch, ``validate_variable_coverage``
checks the runs against the template. During a run, ``strict_replay`` raises
``ReplayDivergenceError`` instead of prompting. Also covered is the accounting
that keeps a run the batch never reached distinct from one that failed, and
the retry list that lets the follow-up batch consume what is still owed.
"""

import json
import os
from datetime import datetime
from io import StringIO

import pytest
from rich.console import Console
from rich.prompt import Confirm, Prompt

from proprep.utils.batch_processor import (
    BatchProcessor,
    load_input_list,
    normalize_input_list,
)
from proprep.utils.session_recorder import (
    InterceptedPrompt,
    ReplayDivergenceError,
    SessionReplayer,
)
from proprep.utils.template_converter import validate_variable_coverage

TEMPLATE = {
    "template": True,
    "template_variables": {
        "input_protein": {"type": "string", "description": "PDB ID", "required": True},
        "unused_req": {"type": "string", "description": "never used", "required": True},
    },
    "interactions": [
        {"type": "prompt", "prompt": "Enter PDB ID", "response": "{{ input_protein }}"},
        {"type": "confirm", "prompt": "Proceed?", "response": "yes"},
        {"type": "input", "prompt": "Notes: ", "response": "none"},
        {"type": "prompt", "prompt": "Set cutoff", "response": "{{ undeclared_var }}"},
    ],
}

FULL_VARS = {"input_protein": "1ABC", "unused_req": "x", "undeclared_var": "5"}


# --------------------------------------------------------------------------
# Coverage check
# --------------------------------------------------------------------------

def test_satisfied_rows_produce_no_errors():
    errors, _ = validate_variable_coverage(TEMPLATE, [dict(FULL_VARS)])
    assert errors == []


def test_undeclared_placeholder_warns():
    """A placeholder absent from template_variables is invisible to the
    required-variable check, which iterates declarations rather than uses."""
    _, warnings = validate_variable_coverage(TEMPLATE, [dict(FULL_VARS)])
    assert any("undeclared_var" in w and "not declared" in w for w in warnings)


def test_missing_declared_required_var_is_an_error():
    """Previously raised from SessionReplayer construction mid-batch, before
    the restoring finally was in scope."""
    errors, _ = validate_variable_coverage(TEMPLATE, [{"input_protein": "1ABC"}])
    assert any("unused_req" in e and "never used" in e for e in errors)


def test_missing_undeclared_placeholder_is_an_error():
    errors, _ = validate_variable_coverage(TEMPLATE, [{"input_protein": "1ABC"}])
    assert any("undeclared_var" in e for e in errors)


def test_blank_cell_is_reported_with_its_run_number():
    rows = [dict(FULL_VARS), dict(FULL_VARS, input_protein="")]
    errors, _ = validate_variable_coverage(TEMPLATE, rows)
    assert any("input_protein" in e and "run(s) 2" in e for e in errors)


def test_unmatched_column_warns_but_does_not_error():
    rows = [dict(FULL_VARS, typo_col="z")]
    errors, warnings = validate_variable_coverage(TEMPLATE, rows)
    assert errors == []
    assert any("typo_col" in w for w in warnings)


def test_normalize_binds_plain_identifiers_and_drops_blanks():
    """Validation and execution must normalize identically, or they disagree
    about what a run supplies."""
    assert normalize_input_list(["1ABC", "  ", "2XYZ"]) == [
        {"input_protein": "1ABC"},
        {"input_protein": "2XYZ"},
    ]


# --------------------------------------------------------------------------
# Strict replay guard
# --------------------------------------------------------------------------

@pytest.fixture
def template_file(tmp_path):
    path = tmp_path / "template.json"
    path.write_text(json.dumps(TEMPLATE))
    return str(path)


@pytest.fixture
def interceptors(template_file):
    """Build installed interceptors, always uninstalling afterwards.

    install() monkey-patches Prompt.ask, Confirm.ask and builtins.input
    globally, so a leak here would corrupt every later test in the session.
    """
    installed = []

    def _make(strict=True, variables=None):
        replayer = SessionReplayer(
            replay_file=template_file, variables=variables or dict(FULL_VARS)
        )
        interceptor = InterceptedPrompt(replayer=replayer, strict_replay=strict)
        replayer.start_replay()
        interceptor.install()
        installed.append(interceptor)
        return replayer, interceptor

    yield _make

    for interceptor in installed:
        interceptor.uninstall()


def test_matching_prompt_still_replays(interceptors):
    interceptors()
    assert Prompt.ask("Enter PDB ID") == "1ABC"


def test_divergent_prompt_raises_instead_of_blocking(interceptors):
    interceptors()
    with pytest.raises(ReplayDivergenceError) as excinfo:
        Prompt.ask("A prompt that was never recorded")

    message = str(excinfo.value)
    assert "never recorded" in message           # the prompt reached
    assert "Enter PDB ID" in message             # what the recording expected


def test_confirm_path_raises(interceptors):
    interceptors()
    with pytest.raises(ReplayDivergenceError):
        Confirm.ask("An unrecorded question?")


def test_input_path_raises(interceptors):
    interceptors()
    with pytest.raises(ReplayDivergenceError):
        input("An unrecorded request: ")


def test_exhausted_recording_raises_with_a_distinct_message(interceptors):
    replayer, _ = interceptors()
    replayer.interaction_index = len(TEMPLATE["interactions"])

    with pytest.raises(ReplayDivergenceError, match="no interactions left"):
        Prompt.ask("Anything at all")


def test_non_strict_replay_still_falls_through(interceptors):
    """Interactive resume hands control to the user, which is correct when
    somebody is present to type."""
    _, interceptor = interceptors(strict=False)
    assert interceptor._replayed_response("prompt", "unrecorded") is None


def test_unresolved_placeholder_raises(interceptors):
    """A surviving {{ placeholder }} is not a usable answer, so replay reaches
    a prompt it cannot satisfy."""
    interceptors(variables={"input_protein": "1ABC", "unused_req": "x"})

    assert Prompt.ask("Enter PDB ID") == "1ABC"
    assert Confirm.ask("Proceed?") is True
    assert input("Notes: ") == "none"
    with pytest.raises(ReplayDivergenceError):
        Prompt.ask("Set cutoff")


def test_peek_pending_does_not_consume(interceptors):
    replayer, _ = interceptors()
    assert replayer.peek_pending()["prompt"] == "Enter PDB ID"
    assert replayer.peek_pending()["prompt"] == "Enter PDB ID"
    assert Prompt.ask("Enter PDB ID") == "1ABC"


# --------------------------------------------------------------------------
# Not-attempted accounting and the retry list
# --------------------------------------------------------------------------

PLANNED = [{"input_protein": f"P{i:03d}", "ph": "7.0"} for i in range(1, 501)]


def _stopped_batch(attempted=3, failed_at=3):
    """A 500-run batch that halted early: runs before failed_at succeeded."""
    processor = BatchProcessor(console=Console(file=StringIO(), width=100))
    processor.start_time = datetime(2026, 7, 22, 10, 0, 0)
    processor.end_time = datetime(2026, 7, 22, 10, 5, 30)
    for i, row in enumerate(PLANNED[:attempted], 1):
        ok = i != failed_at
        processor.results.append(
            {
                "protein": row["input_protein"],
                "variables": row,
                "success": ok,
                "output_dir": f"run_{i:02d}",
                "error": None if ok else "Replay diverged: no recorded answer",
                "timestamp": "2026-07-22T10:00:00",
            }
        )
    return processor


def _report_of(processor, base_dir):
    processor._save_batch_report(base_dir, "tpl.json", PLANNED)
    name = next(p for p in os.listdir(base_dir) if p.startswith("batch_report_"))
    with open(os.path.join(base_dir, name)) as f:
        return json.load(f), name


def test_summary_names_runs_that_were_never_started():
    """Reported as failures, 497 untouched runs would read as 497 broken
    structures rather than one batch that stopped after three."""
    processor = _stopped_batch()
    processor._display_batch_summary("tpl.json", PLANNED)
    output = processor.console.file.getvalue()

    assert "Planned: 500" in output
    assert "Not attempted: 497" in output
    assert "stopped after run 3 of 500" in output


def test_report_counts_reconcile_to_the_planned_total(tmp_path):
    info = _report_of(_stopped_batch(), str(tmp_path))[0]["batch_info"]

    assert info["total_proteins"] == 500
    assert info["attempted"] == 3
    assert (info["successful"], info["failed"], info["not_attempted"]) == (2, 1, 497)
    assert info["stopped_early"] is True
    assert (
        info["successful"] + info["failed"] + info["not_attempted"]
        == info["total_proteins"]
    )


def test_report_lists_the_untouched_runs(tmp_path):
    report = _report_of(_stopped_batch(), str(tmp_path))[0]
    assert len(report["not_attempted"]) == 497
    assert report["not_attempted"][0]["input_protein"] == "P004"


def test_retry_list_pairs_with_its_report(tmp_path):
    _, report_name = _report_of(_stopped_batch(), str(tmp_path))
    retry_name = next(
        p for p in os.listdir(str(tmp_path)) if p.startswith("batch_retry_")
    )
    assert retry_name[len("batch_retry_"):-4] == report_name[len("batch_report_"):-5]


def test_retry_list_round_trips_as_an_input_list(tmp_path):
    """The point of writing an ordinary input list rather than a bespoke
    manifest: it feeds back into --batch-list with no conversion step."""
    _report_of(_stopped_batch(), str(tmp_path))
    retry_file = os.path.join(
        str(tmp_path),
        next(p for p in os.listdir(str(tmp_path)) if p.startswith("batch_retry_")),
    )

    reloaded = normalize_input_list(load_input_list(retry_file))

    assert len(reloaded) == 498                                   # 1 failed + 497 unrun
    assert reloaded[0] == {"input_protein": "P003", "ph": "7.0"}  # failure first
    assert reloaded[1]["input_protein"] == "P004"                 # then original order
    assert reloaded[-1]["input_protein"] == "P500"


def test_no_retry_list_when_every_run_succeeded(tmp_path):
    processor = _stopped_batch(attempted=2, failed_at=0)
    processor._save_batch_report(str(tmp_path), "tpl.json", PLANNED[:2])
    assert not any(p.startswith("batch_retry_") for p in os.listdir(str(tmp_path)))
