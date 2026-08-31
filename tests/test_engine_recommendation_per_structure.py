"""Recommended engine assignment numbers NPT steps within each structure, so
every microstate's density equilibration is 'early NPT' (CPU), not just the
first structure's in the queue."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from test_cphmd_live_path_wiring import _make_manager, SimulationConfig  # noqa: E402

PARAMS = {'min': {'imin': 1}, 'heat': {'imin': 0, 'ntp': 0}, 'npt_density': {'imin': 0, 'ntp': 1},
          'npt_relax': {'imin': 0, 'ntp': 1}, 'nvt': {'imin': 0, 'ntp': 0}, 'prod': {'imin': 0, 'ntp': 0}}


def _sim(wf, step):
    return SimulationConfig(name=step, template_id="t", mdin_path="", engine="pmemd",
                            prmtop=f"{wf}.prmtop", rst7=f"{wf}.rst7", workflow_id=wf,
                            workflow_step=1, simulation_type="production" if step == 'prod' else step)


def test_first_npt_of_every_structure_is_cpu():
    m = _make_manager()
    m._get_simulation_mdin_params = lambda c: PARAMS[c.name]
    queue = [_sim(wf, step) for wf in ("ms001", "ms002", "ms003") for step in PARAMS]
    positions = m._npt_positions(queue)
    results = [(c.workflow_id, c.name, m._classify_for_engine(c, *positions[i])[0]) for i, c in enumerate(queue)]
    for wf in ("ms001", "ms002", "ms003"):
        by_step = {name: target for w, name, target in results if w == wf}
        assert by_step['min'] == 'cpu'
        assert by_step['npt_density'] == 'cpu', wf          # early NPT in EVERY structure
        assert by_step['npt_relax'] == 'gpu'
        assert by_step['heat'] == by_step['nvt'] == by_step['prod'] == 'gpu'


def test_single_npt_step_stays_on_cpu_and_non_npt_get_minus_one():
    m = _make_manager()
    m._get_simulation_mdin_params = lambda c: PARAMS[c.name]
    queue = [_sim("ms001", s) for s in ("min", "npt_density", "prod")]
    assert m._npt_positions(queue) == [(-1, 1), (0, 1), (-1, 1)]
    assert m._classify_for_engine(queue[1], 0, 1)[0] == 'cpu'
