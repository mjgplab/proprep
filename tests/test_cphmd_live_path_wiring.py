#!/usr/bin/env python3
"""
Regression tests for constant-pH (CpHMD) propagation through the LIVE MD setup
path.

Background: the MD Manager has two setup implementations. The constant-pH offer
(`_check_and_offer_cpmd`) and the WorkflowConfig registration that carries
cpin_file/cpmd_settings lived ONLY in `_setup_workflows`, which is unreachable
from the menu. The live path (`_step1_workflow_centric_configuration`) builds
SimulationConfigs but never registered a WorkflowConfig, so every downstream
CpHMD guard resolved the workflow to None via `_get_workflow_for_step` and the
(otherwise complete) machinery silently no-opped: no CPIN staging, no
icnstph/solvph injection, no -cpin/-cpout/-cprestrt flags.

`_maybe_enable_cphmd_for_assignments` is the hook that fixes this. These tests
lock in:
  - it registers a WorkflowConfig (with cpmd_settings + cpin_file) when a CPIN
    in the workspace matches the selected structure's topology;
  - it does NOT offer/register when there is no production step;
  - it does NOT register when the CPIN belongs to a different topology;
  - once a workflow is registered, the downstream contract fires: production
    mdins get icnstph/solvph injected, and the engine command gains
    -cpin/-cpout/-cprestrt.

Run with: pytest tests/test_cphmd_live_path_wiring.py
"""

import io
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from rich.console import Console  # noqa: E402

from proprep.md_prep.molecular_dynamics_manager import (  # noqa: E402
    MolecularDynamicsManager,
    SimulationQueue,
    SimulationConfig,
    WorkflowConfig,
)


class FakeWorkspace:
    """Minimal dict-backed stand-in for the ProPrep workspace."""

    def __init__(self, data=None):
        self._data = dict(data or {})

    def get(self, key, default=None):
        return self._data.get(key, default)

    def set(self, key, value):
        self._data[key] = value

    def has(self, key):
        return key in self._data


class FakeAssignment:
    """Stand-in for the WorkflowCentricStep1Manager assignment object."""

    def __init__(self, workflow_id, prmtop, rst7, name):
        self.custom_workflow_id = workflow_id
        self.structure_pair = {"prmtop": prmtop, "rst7": rst7}
        self.structure_name = name


class FakeProcessor:
    """Minimal processor: console + workspace, as the manager's properties read."""

    def __init__(self, workspace):
        self.console = Console(file=io.StringIO())
        self._workspace = workspace

    def _get_workspace(self):
        return self._workspace


def _make_manager():
    """Build a MolecularDynamicsManager without running its heavy __init__.

    `console`, `workspace`, and `simulation_queue` are read-only properties
    backed by `processor` / `_simulation_queue`, so we wire those backing
    fields; all methods under test are the real, bound implementations.
    """
    m = object.__new__(MolecularDynamicsManager)
    ws = FakeWorkspace()
    m.processor = FakeProcessor(ws)
    m._simulation_queue = SimulationQueue(workspace=ws)
    return m


def _prod(workflow_id="wf1", name="prod", prmtop="sys.prmtop"):
    return SimulationConfig(
        name=name, template_id="t", mdin_path="", engine="pmemd",
        prmtop=prmtop, rst7="sys.rst7",
        workflow_id=workflow_id, workflow_step=1,
        simulation_type="production",
    )


def _min(workflow_id="wf1", name="min", prmtop="sys.prmtop"):
    return SimulationConfig(
        name=name, template_id="t", mdin_path="", engine="pmemd",
        prmtop=prmtop, rst7="sys.rst7",
        workflow_id=workflow_id, workflow_step=1,
        simulation_type="minimization",
    )


CPMD_INFO = {
    "cpin_file": "/tmp/sys.cpin",
    "cpin_config": {"prmtop_file": "sys.prmtop", "simulation_type": "implicit",
                    "igb": 2, "num_residues": 18},
    "cpmd_settings": {"icnstph": 1, "solvph": 7.0, "ntcnstph": 100,
                      "saltcon": 0.1, "igb": 2},
}


# ---------------------------------------------------------------------------
# _maybe_enable_cphmd_for_assignments — the fix itself
# ---------------------------------------------------------------------------

def test_registers_workflow_when_cpin_matches_structure():
    m = _make_manager()
    m.simulation_queue.add_simulation(_prod(prmtop="sys.prmtop"))
    m._check_and_offer_cpmd = lambda prmtop: CPMD_INFO

    assignments = [FakeAssignment("wf1", "sys.prmtop", "sys.rst7", "sys")]
    m._maybe_enable_cphmd_for_assignments(assignments)

    wf = m.simulation_queue._workflows.get("wf1")
    assert wf is not None, "WorkflowConfig should be registered for matching CPIN"
    assert wf.cpin_file == "/tmp/sys.cpin"
    assert wf.cpmd_settings["icnstph"] == 1
    assert wf.cpmd_settings["solvph"] == 7.0
    # Metadata-only registration: queue steps are tracked separately, so the
    # workflow must NOT carry its own copy (that would double-add on persist).
    assert wf.steps == []


def test_no_offer_without_production_step():
    """Minimization-only protocol: never even prompt for CpHMD."""
    m = _make_manager()
    m.simulation_queue.add_simulation(_min())
    called = {"offered": False}

    def _spy(prmtop):
        called["offered"] = True
        return CPMD_INFO

    m._check_and_offer_cpmd = _spy
    m._maybe_enable_cphmd_for_assignments(
        [FakeAssignment("wf1", "sys.prmtop", "sys.rst7", "sys")]
    )

    assert called["offered"] is False
    assert m.simulation_queue._workflows == {}


def test_no_registration_when_cpin_topology_mismatches():
    """A CPIN generated for another protein must not attach to this one."""
    m = _make_manager()
    m.simulation_queue.add_simulation(_prod(prmtop="mine.prmtop"))
    m._check_and_offer_cpmd = lambda prmtop: CPMD_INFO  # cpin is for sys.prmtop

    m._maybe_enable_cphmd_for_assignments(
        [FakeAssignment("wf1", "mine.prmtop", "mine.rst7", "mine")]
    )

    assert m.simulation_queue._workflows == {}


def test_matches_modified_prmtop_for_explicit_solvent():
    """Explicit-solvent CPIN: md pair is swapped to *_cpin.prmtop, which must
    still match via cpin_config['modified_prmtop']."""
    m = _make_manager()
    m.simulation_queue.add_simulation(_prod(prmtop="sys_cpin.prmtop"))
    info = {
        "cpin_file": "/tmp/sys.cpin",
        "cpin_config": {"prmtop_file": "sys.prmtop",
                        "modified_prmtop": "/abs/sys_cpin.prmtop",
                        "simulation_type": "explicit"},
        "cpmd_settings": {"icnstph": 2, "solvph": 7.0, "ntcnstph": 100,
                          "saltcon": 0.1, "ntrelax": 200},
    }
    m._check_and_offer_cpmd = lambda prmtop: info

    m._maybe_enable_cphmd_for_assignments(
        [FakeAssignment("wf1", "sys_cpin.prmtop", "sys.rst7", "sys")]
    )

    assert m.simulation_queue._workflows.get("wf1") is not None


def test_declining_offer_registers_nothing():
    m = _make_manager()
    m.simulation_queue.add_simulation(_prod())
    m._check_and_offer_cpmd = lambda prmtop: None  # user declined / no cpin

    m._maybe_enable_cphmd_for_assignments(
        [FakeAssignment("wf1", "sys.prmtop", "sys.rst7", "sys")]
    )

    assert m.simulation_queue._workflows == {}


# ---------------------------------------------------------------------------
# Downstream contract — proves the machinery fires once a workflow exists
# ---------------------------------------------------------------------------

def _register_cphmd_workflow(m, cpin_file, workflow_id="wf1"):
    m.simulation_queue._workflows[workflow_id] = WorkflowConfig(
        workflow_id=workflow_id, name="sys", description="Constant pH MD",
        system_prmtop="sys.prmtop", initial_rst7="sys.rst7", steps=[],
        cpin_file=cpin_file, cpin_config=CPMD_INFO["cpin_config"],
        cpmd_settings=CPMD_INFO["cpmd_settings"],
    )


def test_injection_adds_icnstph_solvph_for_production():
    m = _make_manager()
    _register_cphmd_workflow(m, "/tmp/sys.cpin")
    mdin = "title\n&cntrl\n  imin=0,\n  nstlim=1000,\n/\n"

    out = m._maybe_inject_cpmd_params(mdin, _prod(), silent=True)

    assert "icnstph=1" in out
    assert "solvph=7.0" in out


def test_injection_skips_non_production():
    m = _make_manager()
    _register_cphmd_workflow(m, "/tmp/sys.cpin")
    mdin = "title\n&cntrl\n  imin=1,\n/\n"

    out = m._maybe_inject_cpmd_params(mdin, _min(), silent=True)

    assert "icnstph" not in out
    assert out == mdin


def test_injection_noop_without_registered_workflow():
    """No WorkflowConfig (the pre-fix state) -> mdin untouched."""
    m = _make_manager()
    mdin = "title\n&cntrl\n  imin=0,\n/\n"

    out = m._maybe_inject_cpmd_params(mdin, _prod(), silent=True)

    assert out == mdin


def test_setup_cpmd_files_emits_cpin_flags(tmp_path):
    m = _make_manager()
    cpin = tmp_path / "sys.cpin"
    cpin.write_text("CPIN")
    _register_cphmd_workflow(m, str(cpin))

    sim_dir = tmp_path / "step1"
    sim_dir.mkdir()

    result = m._setup_cpmd_files(
        _prod(), sim_dir, tmp_path, "sys.prmtop", "prod", silent=True
    )

    assert result is not None
    flags = result["flags"]
    assert "-cpin" in flags
    assert "-cpout" in flags
    assert "-cprestrt" in flags
    # First production step copies the original CPIN into the step dir.
    assert (sim_dir / "sys.cpin").exists()


def test_setup_cpmd_files_noop_without_workflow():
    m = _make_manager()
    result = m._setup_cpmd_files(
        _prod(), Path("/tmp"), Path("/tmp"), "sys.prmtop", "prod", silent=True
    )
    assert result is None


# ---------------------------------------------------------------------------
# _check_and_offer_cpmd — the REAL method (the other tests mock it, so its own
# body was never exercised). Regression: it called self._get_workspace(), which
# MolecularDynamicsManager has no such method — the manager exposes a `workspace`
# property that delegates to processor._get_workspace(). The stray call raised
# AttributeError right after protocol assignment, and the setup wizard's
# try/except swallowed it as "Error in Step 1" and bounced back to file
# selection.
# ---------------------------------------------------------------------------

def test_check_and_offer_cpmd_empty_workspace_returns_none():
    """No cpin in the workspace: must return None, NOT raise AttributeError."""
    m = _make_manager()  # FakeWorkspace starts empty
    assert m._check_and_offer_cpmd("sys.prmtop") is None


def test_check_and_offer_cpmd_reads_workspace_via_property():
    """With a cpin recorded but the file missing on disk, the method still reads
    the workspace (through the `workspace` property) and returns None cleanly."""
    m = _make_manager()
    m.workspace.set("cpin_config", {"simulation_type": "implicit", "num_residues": 5})
    m.workspace.set("cpin_file", "/nonexistent/does_not_exist.cpin")
    assert m._check_and_offer_cpmd("sys.prmtop") is None
