"""Every writer of simulation.mdin (batch, standalone, SLURM x3) must inject the
titration namelist for production steps, not only the live runner."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from test_cphmd_live_path_wiring import _make_manager, _prod, _min, _register_cphmd_workflow  # noqa: E402

MDIN = "Production\n&cntrl\n  ntt=3,\n  temp0=300.0,\n/\n"


def test_imported_mdin_gets_titration_namelist_for_production(tmp_path):
    m = _make_manager()
    _register_cphmd_workflow(m, "/tmp/sys.cpin", workflow_id="wf1")
    src = tmp_path / "prod.mdin"; src.write_text(MDIN)
    sim = _prod("wf1"); sim.mdin_path = str(src)
    dest = tmp_path / "step6" / "simulation.mdin"; dest.parent.mkdir()
    m._stage_mdin_from_path(sim, dest)
    text = dest.read_text()
    assert "icnstph=1," in text and "solvph=7.0," in text
    assert text.index("icnstph") < text.rindex("/")          # inside &cntrl


def test_non_production_and_unregistered_steps_are_untouched(tmp_path):
    m = _make_manager()
    _register_cphmd_workflow(m, "/tmp/sys.cpin", workflow_id="wf1")
    src = tmp_path / "min.mdin"; src.write_text(MDIN)
    sim = _min("wf1"); sim.mdin_path = str(src)
    dest = tmp_path / "simulation.mdin"
    m._stage_mdin_from_path(sim, dest)
    assert "icnstph" not in dest.read_text()
    sim2 = _prod("wf-none"); sim2.mdin_path = str(src)
    m._stage_mdin_from_path(sim2, dest)
    assert "icnstph" not in dest.read_text()


def test_every_mdin_writer_injects():
    """Static guard: each function that writes simulation.mdin from a template
    calls the injector (the raw-copy branches go through _stage_mdin_from_path)."""
    import inspect, re
    from proprep.md_prep.molecular_dynamics_manager import MolecularDynamicsManager
    src = inspect.getsource(MolecularDynamicsManager)
    for name in ("_prepare_workflow_files", "_prepare_standalone_simulation", "_run_amber_simulation",
                 "_generate_independent_slurm_scripts", "_build_workflow_sims_for_slurm",
                 "_generate_workflow_slurm_scripts"):
        body = inspect.getsource(getattr(MolecularDynamicsManager, name))
        assert "_maybe_inject_cpmd_params" in body, name
        assert 'shutil.copy2(sim_config.mdin_path, sim_dir / "simulation.mdin")' not in body, name
