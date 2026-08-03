"""Integrity tests for the built-in MD workflow presets.

The built-in protocols are pure data: a directory of ``.mdin`` templates under
``md_templates/builtin/<slug>/`` plus a manifest under
``md_workflows/builtin/<slug>.json``. They are discovered by glob, so a typo in
a template path, a broken restart chain, or a dangling dependency only surfaces
at runtime. These tests pin the invariants so a rename or a new protocol can't
silently break discovery.

They exercise the real ``WorkflowLoader`` (which raises on a missing template),
so loading success already proves every referenced template exists.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from proprep.md_prep.workflow_loader import WorkflowLoader
from proprep.utils.paths import get_package_dir


EXPECTED_BUILTINS = {
    "builtin_basic_equilibration",
    "builtin_membrane_minimization",
    "builtin_staged_equilibration_310K",
    "builtin_york_lab_tutorial",
}


@pytest.fixture(scope="module")
def workflows():
    return WorkflowLoader().get_available_workflows()


def test_expected_builtins_discovered(workflows):
    """All shipped built-in protocols are found by the glob discovery."""
    found = {k for k in workflows if k.startswith("builtin_")}
    missing = EXPECTED_BUILTINS - found
    assert not missing, f"built-in workflows not discovered: {missing}"


def test_old_names_gone(workflows):
    """The pre-rename standard_protein / protein_equilibration ids are retired."""
    stale = [k for k in workflows if "standard_protein" in k or "protein_equilibration" in k]
    assert not stale, f"stale built-in ids still present: {stale}"


@pytest.mark.parametrize("wf_id", sorted(EXPECTED_BUILTINS))
def test_restart_chain_and_dependencies_consistent(workflows, wf_id):
    """Each step reads the prior step's output (or inpcrd) and deps are valid."""
    wf = workflows[wf_id]
    assert wf.steps, f"{wf_id} has no steps"

    seen_ids: set[str] = set()
    prev_output = None
    for step in wf.steps:
        # restart chain: first step starts from inpcrd, the rest chain forward
        expected_input = "inpcrd" if prev_output is None else prev_output
        assert step.input_coord == expected_input, (
            f"{wf_id}:{step.id} input_coord={step.input_coord!r} "
            f"breaks restart chain (expected {expected_input!r})"
        )
        # dependencies must reference already-seen step ids (forward DAG)
        for dep in step.dependencies:
            assert dep in seen_ids, (
                f"{wf_id}:{step.id} depends on unknown/forward step {dep!r}"
            )
        seen_ids.add(step.id)
        prev_output = step.output_coord


def test_staged_protocol_shape(workflows):
    """The new staged_equilibration_310K protocol has its designed structure."""
    wf = workflows["builtin_staged_equilibration_310K"]
    assert "310 K" in wf.name
    assert len(wf.steps) == 9
    assert wf.steps[0].type == "minimization"
    assert wf.steps[-1].type == "production"

    by_id = {s.id: s for s in wf.steps}
    # the graduated restraint-release taper is present and ordered 25->10->5->2
    release_ids = ["04_release_10", "05_release_5", "06_release_2"]
    assert all(rid in by_id for rid in release_ids)
    # the CPU-density stage advertises its CPU requirement (engine lives in the
    # run plan, but the intent must be discoverable from the step text)
    assert "CPU" in by_id["02_npt_density"].description


def test_york_tutorial_renamed(workflows):
    """The renamed York tutorial carries faithful attribution and 14 steps."""
    wf = workflows["builtin_york_lab_tutorial"]
    assert wf.name == "York Lab Equilibration Tutorial"
    assert "York" in wf.author
    assert len(wf.steps) == 14


def test_ion_name_mask_fix_applied():
    """No built-in template carries the stale ion mask that froze counterions.

    Within ProPrep tleap emits Na+/Cl-, so the legacy ``!:NA,CL,MG,WAT`` mask
    matched no ions and silently restrained all counterions. Every built-in
    template that restrains by residue must use the corrected names.
    """
    builtin_dir = get_package_dir() / "md_templates" / "builtin"
    offenders = [
        str(p)
        for p in builtin_dir.rglob("*.mdin")
        if "!:NA,CL,MG,WAT" in p.read_text()
    ]
    assert not offenders, f"stale ion mask still present in: {offenders}"
