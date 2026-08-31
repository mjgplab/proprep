"""Per-microstate tLEaP logs, multi-topology titration-file generation, and the
MD Manager's per-structure titration lookup."""
import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
from test_cphmd_live_path_wiring import FakeAssignment, _make_manager, _prod, CPMD_INFO  # noqa: E402

from proprep.tleap_prep.tleap_input_generator import TLeapInputGenerator
from proprep.md_prep.molecular_dynamics_manager import TITRATION_PH


# ---- (1) logs -----------------------------------------------------------------

def test_leap_log_for_reads_the_scripts_logfile_line(tmp_path):
    script = tmp_path / "microstate_003_tleap.in"
    script.write_text("logFile microstate_003_leap.log\nsource leaprc.protein.ff14SB\n")
    assert TLeapInputGenerator._leap_log_for(str(script)) == str(tmp_path / "microstate_003_leap.log")
    plain = tmp_path / "plain_tleap.in"
    plain.write_text("source leaprc.protein.ff14SB\n")
    assert TLeapInputGenerator._leap_log_for(str(plain)) == "leap.log"
    assert TLeapInputGenerator._leap_log_for(str(tmp_path / "missing.in")) == "leap.log"


# ---- (3) topology selection ---------------------------------------------------

@pytest.mark.parametrize("text,n,expected", [
    ("all", 5, [1, 2, 3, 4, 5]),
    ("3", 5, [3]),
    ("1,3,5", 5, [1, 3, 5]),
    ("1-4,7", 8, [1, 2, 3, 4, 7]),
    ("4-2", 5, [2, 3, 4]),
    ("", 5, None),
    ("0", 5, None),
    ("6", 5, None),
    ("a,b", 5, None),
])
def test_parse_topology_selection(text, n, expected):
    assert TLeapInputGenerator._parse_topology_selection(text, n) == expected


def test_titratable_signature_ignores_order_and_extra_keys():
    a = [{'resname': 'AS4', 'resnum': 12, 'pka': 4.0}, {'resname': 'PRN', 'resnum': 96, 'pka': 4.8}]
    b = [{'resname': 'PRN', 'resnum': 96, 'chain': 'A'}, {'resname': 'AS4', 'resnum': 12}]
    assert TLeapInputGenerator._titratable_signature(a) == TLeapInputGenerator._titratable_signature(b)
    assert TLeapInputGenerator._titratable_signature(a) != TLeapInputGenerator._titratable_signature(a[:1])


# ---- (4) MD manager: each structure gets its own titration files ---------------

def _workspace_with_two_topologies(tmp_path):
    cpins = {}
    for i in (1, 2):
        c = tmp_path / f"transformed_microstate_00{i}.cpin"; c.write_text("x"); cpins[i] = str(c)
    configs = [
        {'modes': [TITRATION_PH], 'files': {TITRATION_PH: cpins[i]},
         'prmtop_file': f"transformed_microstate_00{i}.prmtop",
         'modified_prmtop': f"transformed_microstate_00{i}_cpin.prmtop",
         'simulation_type': 'explicit', 'igb': None, 'intdiel': 1.0}
        for i in (1, 2)
    ]
    return configs, cpins


def test_each_structure_registers_its_own_cpin(tmp_path):
    m = _make_manager()
    configs, cpins = _workspace_with_two_topologies(tmp_path)
    m.workspace.set('titration_configs', configs)
    m.workspace.set('titration_config', configs[-1])          # legacy "last one" field
    m.simulation_queue.queue.append(_prod("wf1", prmtop="transformed_microstate_001_cpin.prmtop"))

    seen = {}
    def fake_offer(prmtop, resolved=None, applies_to=None):
        seen['resolved'] = resolved; seen['applies_to'] = applies_to
        return {'cpin_file': resolved[1][TITRATION_PH], 'titration_files': resolved[1],
                'cpin_config': resolved[0], 'cpmd_settings': CPMD_INFO['cpmd_settings']}
    m._check_and_offer_cpmd = fake_offer

    assignments = [
        FakeAssignment("wf1", "transformed_microstate_001_cpin.prmtop", "a.rst7", "ms001"),
        FakeAssignment("wf2", "transformed_microstate_002_cpin.prmtop", "b.rst7", "ms002"),
        FakeAssignment("wf3", "other_protein.prmtop", "c.rst7", "other"),
    ]
    m._maybe_enable_cphmd_for_assignments(assignments)

    wfs = m.simulation_queue._workflows
    assert wfs["wf1"].cpin_file == cpins[1]
    assert wfs["wf2"].cpin_file == cpins[2]           # not the last-generated file
    assert "wf3" not in wfs                            # no titration files for it
    assert [name for name, _ in seen['applies_to']] == ["ms001", "ms002"]   # prompt listed both


def test_legacy_single_config_path_is_unchanged():
    m = _make_manager()
    m.workspace.set('cpin_config', CPMD_INFO['cpin_config'])
    m.workspace.set('cpin_file', CPMD_INFO['cpin_file'])
    m.simulation_queue.queue.append(_prod("wf1"))
    m._check_and_offer_cpmd = lambda prmtop, **kw: CPMD_INFO
    m._maybe_enable_cphmd_for_assignments([FakeAssignment("wf1", "sys.prmtop", "sys.rst7", "sys")])
    assert m.simulation_queue._workflows["wf1"].cpin_file == "/tmp/sys.cpin"
