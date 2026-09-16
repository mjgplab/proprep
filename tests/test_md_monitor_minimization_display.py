"""The monitor reports a minimization as cycles, and lists each step's output.

Seen during the 1.19.1 workshop check: "monitoring the minimization" showed
Step 2,000 / Time 4.000 ps / Temperature 75.89 K, which is the heating step.
"""

from __future__ import annotations

import os
import time
import types
from pathlib import Path

from rich.console import Console

from proprep.md_prep.amber_workflow_components import AMBERMonitor
from proprep.md_prep.molecular_dynamics_manager import MolecularDynamicsManager as MDM

MIN_OUT = """\
   imin    =       1, nmropt  =       0

   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER
      1      -1.5166E+05     1.0223E+03     2.2397E+05     HD21     5639
    500      -1.5400E+05     1.2000E+00     7.1612E+00     CA        12
"""
MD_OUT = """\
   imin    =       0, nmropt  =       1
 NSTEP =     2000   TIME(PS) =       4.000  TEMP(K) =    75.89  PRESS =     0.0
 Etot   =    -18893.1800  EKtot   =       900.0000  EPtot      =    -19793.1800
"""


def _render(path: Path) -> str:
    mon = AMBERMonitor(str(path)); mon.parse_amber_output()
    console = Console(width=80, record=True)
    m = object.__new__(MDM); m.processor = types.SimpleNamespace(console=console)
    m._show_simulation_status(mon)
    return console.export_text()


def test_minimization_shows_cycles_not_time(tmp_path):
    p = tmp_path / "min.mdout"; p.write_text(MIN_OUT)
    out = _render(p)
    assert "Cycle" in out and "500" in out
    assert "Time (ps)" not in out and "Temperature" not in out
    assert "Max gradient" in out and "7.1612" in out


def test_md_still_shows_time_and_temperature(tmp_path):
    p = tmp_path / "heat.mdout"; p.write_text(MD_OUT)
    out = _render(p)
    assert "Time (ps)" in out and "4.000" in out and "75.89" in out
    assert "Cycle" not in out


def test_recent_mdouts_listed_per_step_newest_first(tmp_path):
    batch = tmp_path / "simulations" / "batch_x"
    a = batch / "01_Energy_Minimization"; b = batch / "02_System_Heating"
    a.mkdir(parents=True); b.mkdir()
    (a / "min.mdout").write_text(MIN_OUT); time.sleep(0.02); (b / "heat.mdout").write_text(MD_OUT)
    os.utime(b / "heat.mdout", None)
    found = MDM._scan_recent_mdouts(tmp_path / "simulations", [])
    labels = [label for label, _ in found]
    assert labels == ["batch_x/02_System_Heating", "batch_x/01_Energy_Minimization"]
    # a tracked step is left to the tracked list
    tracked = [("x", {"mdout_file": b / "heat.mdout", "sim_dir": b}, "running")]
    assert [l for l, _ in MDM._scan_recent_mdouts(tmp_path / "simulations", tracked)] == ["batch_x/01_Energy_Minimization"]
