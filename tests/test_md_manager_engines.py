"""MD Manager engines: sander.MPI is offered and run, and the menu is constant.

AmberTools ships sander and sander.MPI only; pmemd, pmemd.MPI and pmemd.cuda
come with Amber. Up to 1.19.0 the MD Manager offered only pmemd.MPI as the
multi-core engine, so a workshop attendee on a constructor install had no
working parallel option and the default (pmemd) did not exist.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from proprep.md_prep.molecular_dynamics_manager import MolecularDynamicsManager as MDM


def _manager():
    return object.__new__(MDM)  # __init__ needs a processor; the methods under test do not


def test_engine_menu_is_constant_and_ordered():
    assert MDM.ENGINES == ["sander", "sander.MPI", "pmemd", "pmemd.MPI", "pmemd.cuda"]
    assert _manager()._get_available_engines() == MDM.ENGINES


def _cmd(engine, hw):
    m = _manager()
    f = Path("x")
    return m._build_amber_command_line(engine, f.with_name("s.mdin"), f.with_name("p.prmtop"),
                                       f.with_name("c.rst"), f.with_name("o.out"), f.with_name("r.rst"),
                                       f.with_name("t.nc"), f.with_name("ref.rst"), hw)


def test_sander_mpi_runs_under_mpirun():
    cmd = _cmd("sander.MPI", {"mpi_tasks": 4})
    assert cmd[:4] == ["mpirun", "-np", "4", "sander.MPI"]
    assert "-inf" not in cmd  # sander takes no -inf


def test_pmemd_mpi_still_runs_under_mpirun():
    cmd = _cmd("pmemd.MPI", {"mpi_tasks": 8})
    assert cmd[:4] == ["mpirun", "-np", "8", "pmemd.MPI"]
    assert "-inf" in cmd


def test_serial_engines_run_bare():
    assert _cmd("sander", {})[0] == "sander"
    assert _cmd("pmemd", {})[0] == "pmemd"


def test_recommended_cpu_engine_falls_back_to_sander_mpi(monkeypatch):
    installed = {"sander", "sander.MPI"}
    monkeypatch.setattr(shutil, "which", lambda name: f"/bin/{name}" if name in installed else None)
    m = _manager()
    assert m._cpu_parallel_engine() == "sander.MPI"
    assert m._default_engine() == "sander.MPI"
    assert m._engine_label("pmemd.MPI") == "pmemd.MPI (not installed)"
    assert m._engine_label("sander.MPI") == "sander.MPI"


def test_recommended_cpu_engine_prefers_pmemd_mpi_when_present(monkeypatch):
    monkeypatch.setattr(shutil, "which", lambda name: f"/bin/{name}")
    m = _manager()
    assert m._cpu_parallel_engine() == "pmemd.MPI"
    assert m._default_engine() == "pmemd.cuda"


def test_workshop_workflow_and_templates_exist():
    root = Path(__file__).resolve().parents[1] / "src/proprep"
    import json
    wf = json.load(open(root / "md_workflows/builtin/workshop_cpu.json"))
    for step in wf["steps"]:
        assert (root / "md_templates" / step["template"]).is_file(), step["template"]
    heat = (root / "md_templates/builtin/workshop_cpu/01_heating.mdin").read_text()
    # the &wt ramp must end inside nstlim, or the temperature never reaches 300 K
    assert "nstlim=10000" in heat and "ISTEP2=10000" in heat
