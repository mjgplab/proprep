"""A simulation still running from an earlier session is restored, not dropped.

_restore_running_simulations() called a method that never existed
(_restore_workflow_pending_steps) for any running step that belonged to a
workflow; the AttributeError was caught, printed as a warning, and the
live process was left out of the kept list.
"""

from __future__ import annotations

import json
import os
import types
from datetime import datetime

from proprep.md_prep.molecular_dynamics_manager import MolecularDynamicsManager as MDM


class _Console:
    def __init__(self):
        self.lines = []

    def print(self, *a, **k):
        self.lines.append(" ".join(str(x) for x in a))


def test_running_workflow_step_is_restored(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    sims = tmp_path / "simulations"; sims.mkdir()
    sim_dir = sims / "02_System_Heating"; sim_dir.mkdir()
    (sims / ".running_processes.json").write_text(json.dumps({
        "structure_System Heating": {
            "pid": os.getpid(),  # alive for the duration of the test
            "sim_dir": str(sim_dir),
            "mdout_file": str(sim_dir / "heating.mdout"),
            "started_at": datetime.now().isoformat(),
        }
    }))

    m = object.__new__(MDM)
    console = _Console()
    m.processor = types.SimpleNamespace(console=console)
    m._load_sim_config_from_dir = lambda d: types.SimpleNamespace(
        name="System Heating", workflow_id="wf1", depends_on="Energy Minimization")
    m._restore_pending_workflow_steps = lambda: 0

    m._restore_running_simulations()

    assert "structure_System Heating" in m.running_simulations
    assert m.running_simulations["structure_System Heating"]["restored"] is True
    assert not any("Could not restore" in line for line in console.lines), console.lines
    kept = json.loads((sims / ".running_processes.json").read_text())
    assert "structure_System Heating" in kept
