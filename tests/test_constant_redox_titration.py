#!/usr/bin/env python3
"""
Regression tests for constant-Redox (cein) and combined (cpein) titration.

Background: ProPrep could generate a `cpin` for constant-pH MD but had no path
to `ceinutil.py` or `cpeinutil.py`, so the constant_E heme sets it can now
build (residue name HEH) had no way to actually titrate in redox. The three
utilities differ in ways that fail *quietly* rather than loudly, so most of
what follows pins down those specific failure modes:

  - Flag triples. sander parses three parallel sets (mdfil.F90:290-341).
    Emitting the wrong triple for a generated file is a silent mismatch, not
    an error. Note that `-cpin` and `-cein` are COMPLEMENTARY and are passed
    together for a system that titrates in both (Amber tutorial 33); the only
    exclusion in mdfil.F90 (line 586) is one-directional and guards `-cpein`.
  - Reference-energy coverage. ceinutil accepts any of igb 1/2/5/7/8 and
    intdiel 1/2, but HEH has energies for far fewer than that; the unsupported
    combinations produce None-filled output rather than a diagnostic.
  - Atom counts. All three utilities filter termini by comparing a residue's
    atom count against ParmEd's definition and drop mismatches without
    comment. An 87-atom HEH that comes out at 86 simply stops titrating.

Several tests assert ProPrep's hardcoded tables against ParmEd itself, so an
AmberTools upgrade that changes the shipped reference energies fails here
rather than in someone's simulation.

Run with: pytest tests/test_constant_redox_titration.py
"""

import io
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from rich.console import Console  # noqa: E402

from proprep.utils.titration_modes import (  # noqa: E402
    PH,
    REDOX,
    PHREDOX,
    PH_RESIDUE_PKA,
    REDOX_RESIDUE_EO,
    PHREDOX_RESIDUE_NAMES,
    get_mode,
    modes_for_residue_names,
    partition_residues_by_mode,
    engine_flags_for_modes,
    mdin_keyword_sets,
    parmed_titratable_atom_count,
    parmed_residue_type,
)
from proprep.md_prep.molecular_dynamics_manager import (  # noqa: E402
    MolecularDynamicsManager,
    SimulationQueue,
    SimulationConfig,
    WorkflowConfig,
)

REPO_ROOT = Path(__file__).parent.parent
HEME_DIR = (REPO_ROOT / "src/proprep/forcefield_params/specialized_residues"
            / "heme/bis_his_c_type/oxidized/low_spin")

parmed = pytest.importorskip("parmed")


# ---------------------------------------------------------------------------
# Mode descriptors: the flag triples sander actually parses
# ---------------------------------------------------------------------------

def test_flag_triples_match_sander_parser():
    """Each mode carries the exact triple sander/pmemd accept.

    Hardcoded rather than derived: these strings are an external contract with
    the Fortran argument parser, so the test should fail if anyone "tidies"
    them.
    """
    assert get_mode(PH).flag_input == "-cpin"
    assert get_mode(PH).flag_output == "-cpout"
    assert get_mode(PH).flag_restart == "-cprestrt"

    assert get_mode(REDOX).flag_input == "-cein"
    assert get_mode(REDOX).flag_output == "-ceout"
    assert get_mode(REDOX).flag_restart == "-cerestrt"

    assert get_mode(PHREDOX).flag_input == "-cpein"
    assert get_mode(PHREDOX).flag_output == "-cpeout"
    assert get_mode(PHREDOX).flag_restart == "-cperestrt"


def test_modes_map_to_the_right_utility():
    assert get_mode(PH).utility == "cpinutil.py"
    assert get_mode(REDOX).utility == "ceinutil.py"
    assert get_mode(PHREDOX).utility == "cpeinutil.py"


def test_cpin_and_cein_are_passed_together():
    """A pH+redox run gets all six flags in one command.

    Amber tutorial 33 section 3 runs:
        -cpin mp8.cpin -cpout ... -cprestrt ... -cein mp8.cein -ceout ...
        -cerestrt ...
    sander never compares cpin_specified against cein_specified; the only
    exclusion (mdfil.F90:586) is one-directional and guards cpein alone.
    """
    flags = engine_flags_for_modes({PH: "sys.cpin", REDOX: "sys.cein"}, "prod")
    assert flags == [
        "-cpin", "sys.cpin", "-cpout", "prod.cpout",
        "-cprestrt", "prod.cprestrt",
        "-cein", "sys.cein", "-ceout", "prod.ceout",
        "-cerestrt", "prod.cerestrt",
    ]


def test_cpein_flags_do_not_overlap_the_others():
    """cpein IS exclusive with both others, so it shares no flag with them."""
    combined = set(get_mode(PHREDOX).engine_flags("sys.cpein", "step1"))
    ph = set(get_mode(PH).engine_flags("sys.cpin", "step1"))
    redox = set(get_mode(REDOX).engine_flags("sys.cein", "step1"))
    assert not (combined & {"-cpin", "-cein"})
    assert "-cpein" not in (ph | redox)


def test_combined_run_sets_both_keyword_families():
    """Tutorial 33 sets icnstph/solvph/ntcnstph AND icnste/solve/ntcnste in
    the same &cntrl."""
    families = mdin_keyword_sets([PH, REDOX])
    assert [f.mdin_flag_keyword for f in families] == ["icnstph", "icnste"]
    assert [f.mdin_setpoint_keyword for f in families] == ["solvph", "solve"]

    # A cpein drives both families from a single file.
    joint = mdin_keyword_sets([PHREDOX])
    assert [f.mdin_flag_keyword for f in joint] == ["icnstph", "icnste"]


def test_engine_flags_shape():
    flags = get_mode(REDOX).engine_flags("sys.cein", "step3")
    assert flags == [
        "-cein", "sys.cein",
        "-ceout", "step3.ceout",
        "-cerestrt", "step3.cerestrt",
    ]


def test_unknown_mode_key_is_rejected_but_none_defaults_to_ph():
    """None must mean constant pH: that is what every record written before
    the redox modes existed carries."""
    assert get_mode(None).key == PH
    with pytest.raises(ValueError):
        get_mode("nonsense")


def test_only_ceinutil_lacks_output_prmtop_support():
    """ceinutil.py imports changeRadii/change and never calls them; it has no
    -op flag. cpinutil and cpeinutil both do."""
    assert get_mode(PH).supports_output_prmtop is True
    assert get_mode(PHREDOX).supports_output_prmtop is True
    assert get_mode(REDOX).supports_output_prmtop is False


# ---------------------------------------------------------------------------
# Mode auto-detection
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("resnames,expected", [
    (["AS4", "GL4", "HIP"], [PH]),
    (["HEH"], [REDOX]),
    (["AS4", "HEH"], [PH, REDOX]),
    (["HEH", "PRN"], [PH, REDOX]),
    (["TYX"], [PHREDOX]),
    (["ALA", "WAT"], []),
    ([], []),
])
def test_mode_detection(resnames, expected):
    assert modes_for_residue_names(resnames) == expected


def test_detection_is_case_insensitive():
    assert modes_for_residue_names(["heh"]) == [REDOX]


def test_heme_plus_carboxylate_gets_two_complementary_files():
    """The constant_E + constant_pH heme set emits HEH *and* two PRN.

    That is exactly tutorial 33's microperoxidase case, which generates a cein
    for HEH and a cpin for PRN/GL4 and runs them together -- NOT a cpein.
    """
    assert modes_for_residue_names(["HEH", "PRN", "PRN", "AS4"]) == [PH, REDOX]


def test_proton_coupled_site_forces_the_joint_file():
    """TYX is typ='phredox': cpinutil drops it (not "ph") and ceinutil drops
    it (not "redox"), so only cpein can carry it -- and cpein cannot be
    combined with either of the others, so it absorbs everything."""
    assert modes_for_residue_names(["TYX", "AS4", "HEH"]) == [PHREDOX]


def test_partition_sends_each_residue_to_the_right_file():
    """Mixing families is a hard error at the utility, not a silent skip:
    ceinutil rejects a -resnums list containing an AS4 outright."""
    residues = [
        {"resname": "HEH", "resnum": 1},
        {"resname": "PRN", "resnum": 2},
        {"resname": "AS4", "resnum": 3},
    ]
    parts = partition_residues_by_mode(residues, [PH, REDOX])
    assert [r["resname"] for r in parts[PH]] == ["PRN", "AS4"]
    assert [r["resname"] for r in parts[REDOX]] == ["HEH"]


# ---------------------------------------------------------------------------
# Cross-checks against ParmEd (drift protection for AmberTools upgrades)
# ---------------------------------------------------------------------------

def test_residue_tables_match_parmed_types():
    """Each mode's residue table must agree with ParmEd's `typ`, because that
    is the field the utilities filter on."""
    for resname in PH_RESIDUE_PKA:
        assert parmed_residue_type(resname) == "ph", resname
    for resname in REDOX_RESIDUE_EO:
        assert parmed_residue_type(resname) == "redox", resname
    for resname in PHREDOX_RESIDUE_NAMES:
        assert parmed_residue_type(resname) == "phredox", resname


def test_heh_is_the_only_redox_residue():
    """If AmberTools ever ships a second redox residue, REDOX_RESIDUE_EO must
    learn about it or that residue silently never titrates."""
    from parmed.amber import titratable_residues as res
    redox = {n for n in res.titratable_residues
             if getattr(res, n).typ == "redox"}
    assert redox == set(REDOX_RESIDUE_EO)


def test_eo_values_match_parmed():
    from parmed.amber import titratable_residues as res
    for resname, eo in REDOX_RESIDUE_EO.items():
        assert getattr(res, resname).Eo == pytest.approx(eo)


def _igb_coverage(resname, explicit):
    """igb values for which EVERY state of `resname` has a reference energy.

    A state with `None` is what makes an unsupported combination produce a
    meaningless file, so a mode may only offer an igb where all states are set.
    """
    from parmed.amber import titratable_residues as res
    residue = getattr(res, resname)
    covered = set()
    for igb in (1, 2, 5, 7, 8):
        ok = True
        for state in residue.states:
            energy = state.refene
            if explicit:
                energy = energy.solvent_energies_obj if hasattr(
                    energy, 'solvent_energies_obj') else energy.solvent
            if getattr(energy, f'igb{igb}') is None:
                ok = False
                break
        if ok:
            covered.add(igb)
    return covered


def test_redox_igb_coverage_matches_parmed():
    """The igb sets ProPrep offers must match HEH's actual coverage.

    Verified against ParmEd 4.3.1: implicit {2,5,7,8}, explicit {2,5,7}.
    igb=1 has no reduced-state energy at all; igb=8 has one only for implicit.
    """
    mode = get_mode(REDOX)
    assert mode.allowed_igb('implicit') == _igb_coverage('HEH', explicit=False)
    assert mode.allowed_igb('explicit') == _igb_coverage('HEH', explicit=True)


def test_redox_rejects_intdiel_2():
    """HEH registers its dielc2 energies with no arguments, so every intdiel=2
    value is None -- the mode must not offer it."""
    from parmed.amber import titratable_residues as res
    reduced = res.HEH.states[1]
    assert all(getattr(reduced.refene.dielc2, f'igb{i}') is None
               for i in (1, 2, 5, 7, 8))
    assert 2.0 not in get_mode(REDOX).allowed_intdiel
    assert get_mode(REDOX).allowed_intdiel == {1.0}


def test_ph_mode_still_allows_everything_it_used_to():
    """Constant pH behaviour must not narrow as a side effect of this work."""
    mode = get_mode(PH)
    assert mode.allowed_igb('implicit') == {1, 2, 5, 7, 8}
    assert mode.allowed_igb('explicit') == {1, 2, 5, 7, 8}
    assert mode.allowed_intdiel == {1.0, 2.0}


# ---------------------------------------------------------------------------
# The silent-skip guard: HEH atom counts
# ---------------------------------------------------------------------------

def test_parmed_heh_atom_count_is_known():
    assert parmed_titratable_atom_count("HEH") is not None


@pytest.mark.parametrize("libname", [
    "Henriques_HEH_RESP_constE_FixedpH.lib",
    "Henriques_HEH_RESP_constE_constpH.lib",
])
def test_shipped_heh_libs_match_parmed_atom_count(libname):
    """ProPrep's HEH units must have exactly as many atoms as ParmEd's HEH.

    This is the whole point of naming the residue HEH. All three utilities
    filter termini by atom count and `continue` past any mismatch with no
    message and no non-zero exit, so a lib that is one atom off produces a
    heme that never titrates. Asserted against ParmEd rather than a literal so
    an AmberTools change is caught too.
    """
    lib = HEME_DIR / libname
    assert lib.exists(), f"missing library: {lib}"

    expected = parmed_titratable_atom_count("HEH")
    in_table = False
    count = 0
    for line in lib.read_text().splitlines():
        if line.startswith("!entry.HEH.unit.atoms table"):
            in_table = True
            continue
        if in_table:
            if line.startswith("!entry"):
                break
            if line.strip():
                count += 1

    assert count == expected, (
        f"{libname} declares {count} HEH atoms but ParmEd expects {expected}; "
        f"ceinutil would drop this residue silently"
    )


# ---------------------------------------------------------------------------
# Generator-side validation
# ---------------------------------------------------------------------------

def _make_generator():
    """Build a TleapInputGenerator without running its heavy __init__."""
    from proprep.tleap_prep.tleap_input_generator import TLeapInputGenerator

    class _Proc:
        def __init__(self):
            self.console = Console(file=io.StringIO())

    g = object.__new__(TLeapInputGenerator)
    g.processor = _Proc()
    return g


def test_validate_rejects_igb1_for_redox():
    g = _make_generator()
    ok, problems = g._validate_mode_constraints(
        get_mode(REDOX), 'implicit', igb=1, intdiel=1.0)
    assert not ok
    assert any("igb=1" in p for p in problems)


def test_validate_rejects_intdiel2_for_redox():
    g = _make_generator()
    ok, problems = g._validate_mode_constraints(
        get_mode(REDOX), 'implicit', igb=2, intdiel=2.0)
    assert not ok
    assert any("intdiel" in p for p in problems)


def test_validate_accepts_supported_redox_combination():
    g = _make_generator()
    ok, problems = g._validate_mode_constraints(
        get_mode(REDOX), 'implicit', igb=2, intdiel=1.0)
    assert ok, problems


def test_missing_op_does_not_block_explicit_redox():
    """ceinutil has no -op, but that never blocks a run.

    In a cpin+cein run cpinutil writes the radii-corrected prmtop and both
    files share it (tutorial 33 uses one mp8_es.new.prmtop for the whole
    simulation). In a redox-only run ProPrep writes it with ParmEd. Either
    way the missing flag is informational.
    """
    g = _make_generator()
    ok, problems = g._validate_mode_constraints(
        get_mode(REDOX), 'explicit', igb=None, intdiel=1.0)
    assert ok, problems

    ok, problems = g._validate_mode_constraints(
        get_mode(REDOX), 'explicit', igb=None, intdiel=1.0,
        companion_writes_prmtop=True)
    assert ok, problems


def test_validate_explicit_ph_is_unaffected():
    g = _make_generator()
    ok, problems = g._validate_mode_constraints(
        get_mode(PH), 'explicit', igb=None, intdiel=1.0)
    assert ok, problems


def test_heh_state_descriptions_are_redox_not_protonation():
    g = _make_generator()
    assert "oxidized" in g._get_cpin_state_description("HEH", 0)
    assert "reduced" in g._get_cpin_state_description("HEH", 1)


def test_titration_file_residue_count_parsing(tmp_path):
    g = _make_generator()
    cein = tmp_path / "sys.cein"
    cein.write_text(
        "&CNSTPHE_LIMITS ntres=1, maxh=2, natchrg=174, ntstates=2, /\n"
        "&CNSTE\n CHRGDAT= 0.666,\n TRESCNT=1,\n"
        " RESNAME='System: sys', 'Residue: HEH 1',\n/\n"
    )
    assert g._count_titration_file_residues(str(cein)) == 1

    empty = tmp_path / "empty.cein"
    empty.write_text("&CNSTE\n/\n")
    assert g._count_titration_file_residues(str(empty)) is None

    assert g._count_titration_file_residues(str(tmp_path / "nope")) is None


# ---------------------------------------------------------------------------
# MD-side plumbing
# ---------------------------------------------------------------------------

class FakeWorkspace:
    def __init__(self, data=None):
        self._data = dict(data or {})

    def get(self, key, default=None):
        return self._data.get(key, default)

    def set(self, key, value):
        self._data[key] = value

    def has(self, key):
        return key in self._data


class FakeProcessor:
    def __init__(self, workspace):
        self.console = Console(file=io.StringIO())
        self._workspace = workspace

    def _get_workspace(self):
        return self._workspace


def _make_manager(workspace_data=None):
    m = object.__new__(MolecularDynamicsManager)
    ws = FakeWorkspace(workspace_data)
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


def test_workflow_config_defaults_to_ph_when_mode_absent():
    """Workflows persisted before this change carry no titration_mode; they
    must keep behaving as constant pH."""
    wf = WorkflowConfig(
        workflow_id="wf1", name="n", description="d",
        system_prmtop="sys.prmtop", initial_rst7="sys.rst7", steps=[],
        cpin_file="/tmp/sys.cpin",
    )
    assert wf.titration_mode is None
    assert get_mode(wf.titration_mode).flag_input == "-cpin"


def test_titration_mode_round_trips_through_the_workspace():
    m = _make_manager()
    m.simulation_queue.add_simulation(_prod())
    m.simulation_queue._workflows["wf1"] = WorkflowConfig(
        workflow_id="wf1", name="n", description="Constant Redox potential MD",
        system_prmtop="sys.prmtop", initial_rst7="sys.rst7", steps=[],
        cpin_file="/tmp/sys.cein",
        cpin_config={"prmtop_file": "sys.prmtop"},
        cpmd_settings={"icnste": 1, "solve": -0.203, "ntcnste": 100},
        titration_mode=REDOX,
    )
    m.simulation_queue._sync_to_workspace()

    stored = m.workspace.get("md_workflows")["wf1"]
    assert stored["titration_mode"] == REDOX

    reloaded = SimulationQueue(workspace=m.workspace)
    _ = reloaded.queue  # triggers the load
    assert reloaded._workflows["wf1"].titration_mode == REDOX


def test_setup_titration_files_emits_redox_flags(tmp_path):
    """A workflow whose mode is redox must produce -cein/-ceout/-cerestrt."""
    cein = tmp_path / "sys.cein"
    cein.write_text("&CNSTE\n/\n")
    sim_dir = tmp_path / "prod"
    sim_dir.mkdir()

    m = _make_manager()
    config = _prod()
    workflow = WorkflowConfig(
        workflow_id="wf1", name="n", description="d",
        system_prmtop="sys.prmtop", initial_rst7="sys.rst7", steps=[config],
        cpin_file=str(cein), titration_mode=REDOX,
    )
    m._get_workflow_for_step = lambda c: workflow

    result = m._setup_cpmd_files(config, sim_dir, tmp_path, "sys.prmtop",
                                "prod", silent=True)
    assert result["flags"] == [
        "-cein", "sys.cein",
        "-ceout", "prod.ceout",
        "-cerestrt", "prod.cerestrt",
    ]
    assert (sim_dir / "sys.cein").exists()


def test_setup_titration_files_still_emits_cpin_flags_for_ph(tmp_path):
    cpin = tmp_path / "sys.cpin"
    cpin.write_text("&CNSTPH\n/\n")
    sim_dir = tmp_path / "prod"
    sim_dir.mkdir()

    m = _make_manager()
    config = _prod()
    workflow = WorkflowConfig(
        workflow_id="wf1", name="n", description="d",
        system_prmtop="sys.prmtop", initial_rst7="sys.rst7", steps=[config],
        cpin_file=str(cpin),  # no titration_mode -> constant pH
    )
    m._get_workflow_for_step = lambda c: workflow

    result = m._setup_cpmd_files(config, sim_dir, tmp_path, "sys.prmtop",
                                "prod", silent=True)
    assert result["flags"] == [
        "-cpin", "sys.cpin",
        "-cpout", "prod.cpout",
        "-cprestrt", "prod.cprestrt",
    ]


def test_setup_titration_files_emits_combined_flags(tmp_path):
    cpein = tmp_path / "sys.cpein"
    cpein.write_text("&CNSTPHE\n/\n")
    sim_dir = tmp_path / "prod"
    sim_dir.mkdir()

    m = _make_manager()
    config = _prod()
    workflow = WorkflowConfig(
        workflow_id="wf1", name="n", description="d",
        system_prmtop="sys.prmtop", initial_rst7="sys.rst7", steps=[config],
        cpin_file=str(cpein), titration_mode=PHREDOX,
    )
    m._get_workflow_for_step = lambda c: workflow

    flags = m._setup_cpmd_files(config, sim_dir, tmp_path, "sys.prmtop",
                                "prod", silent=True)["flags"]
    assert flags[0] == "-cpein"
    assert "-cpin" not in flags and "-cein" not in flags


def test_restart_chaining_uses_the_modes_extension(tmp_path):
    """A constant-redox chain must look for a .cerestrt, not a .cprestrt that
    sander will never write."""
    run_dir = tmp_path / "run"
    prev_dir = run_dir / "prod1"
    prev_dir.mkdir(parents=True)
    (prev_dir / "prod1.cerestrt").write_text("restart")
    cur_dir = run_dir / "prod2"
    cur_dir.mkdir()

    m = _make_manager()
    step1 = _prod(name="prod1")
    step2 = _prod(name="prod2")
    step2.workflow_step = 1
    workflow = WorkflowConfig(
        workflow_id="wf1", name="n", description="d",
        system_prmtop="sys.prmtop", initial_rst7="sys.rst7",
        steps=[step1, step2], cpin_file="/tmp/sys.cein", titration_mode=REDOX,
    )
    m._get_workflow_for_step = lambda c: workflow

    found = m._find_previous_cprestrt(step2, cur_dir, 'cerestrt')
    assert found is not None and found.name == "prod1.cerestrt"

    # The pH extension must not match the redox restart file.
    assert m._find_previous_cprestrt(step2, cur_dir, 'cprestrt') is None


def test_workspace_resolution_returns_every_generated_file():
    m = _make_manager({
        "titration_config": {
            "modes": [PH, REDOX],
            "files": {PH: "/tmp/sys.cpin", REDOX: "/tmp/sys.cein"},
        },
        "cpin_config": {"num_residues": 18},
        "cpin_file": "/tmp/sys.cpin",
    })
    config, files = m._resolve_titration_workspace_config()
    assert files == {PH: "/tmp/sys.cpin", REDOX: "/tmp/sys.cein"}


def test_workspace_resolution_falls_back_to_legacy_cpin_keys():
    """Workspaces written before the redox files have only cpin_config/
    cpin_file, and those are constant pH by definition."""
    m = _make_manager({
        "cpin_config": {"num_residues": 18},
        "cpin_file": "/tmp/sys.cpin",
    })
    config, files = m._resolve_titration_workspace_config()
    assert files == {PH: "/tmp/sys.cpin"}


def test_workspace_resolution_returns_none_when_nothing_generated():
    assert _make_manager()._resolve_titration_workspace_config() is None


# ---------------------------------------------------------------------------
# The generated run script (a separate emission site from _setup_cpmd_files)
# ---------------------------------------------------------------------------

def _write_run_script(tmp_path, mode_key, titration_filename, n_prod=2):
    """Drive _create_workflow_run_script for a workflow with `n_prod` steps."""
    m = _make_manager()
    steps = [_prod(name=f"prod{i}") for i in range(1, n_prod + 1)]
    for i, step in enumerate(steps):
        step.workflow_step = i
        step.engine = "pmemd"
        step.hardware_config = {}
    m.simulation_queue._workflows["wf1"] = WorkflowConfig(
        workflow_id="wf1", name="n", description="d",
        system_prmtop="sys.prmtop", initial_rst7="sys.rst7", steps=steps,
        cpin_file=str(tmp_path / titration_filename),
        titration_mode=mode_key,
    )
    m._get_workflow_for_step = lambda c: m.simulation_queue._workflows["wf1"]
    m._is_production_step = lambda c: True
    m._get_template_type = lambda t: "production"

    workflow_dir = tmp_path / "wf"
    workflow_dir.mkdir()
    m._create_workflow_run_script(
        workflow_dir, steps, [f"{s.name}/simulation.mdin" for s in steps],
        "sys.prmtop", "sys.rst7", extended_production_cycles=0)
    return (workflow_dir / "run_workflow.sh").read_text()


def test_run_script_emits_redox_flags(tmp_path):
    script = _write_run_script(tmp_path, REDOX, "sys.cein")
    assert "-cein" in script
    assert "-ceout" in script
    assert "-cerestrt" in script
    # A cein must never be handed to sander under the constant-pH flags.
    assert "-cpin" not in script
    assert "-cprestrt" not in script


def test_run_script_chains_redox_restarts(tmp_path):
    """The second production step must read the first step's .cerestrt --
    a .cprestrt would never exist for a constant-redox run.

    Steps are keyed step1/step2 in the generated script, independent of the
    SimulationConfig names.
    """
    script = _write_run_script(tmp_path, REDOX, "sys.cein", n_prod=2)
    # Step 1 reads the freshly staged cein; step 2 continues from step 1's
    # redox-state restart file.
    assert "sys.cein" in script
    assert "../step1/step1.cerestrt" in script
    assert "cprestrt" not in script


def test_run_script_emits_combined_flags(tmp_path):
    script = _write_run_script(tmp_path, PHREDOX, "sys.cpein")
    assert "-cpein" in script and "-cpeout" in script
    assert "-cperestrt" in script
    # -cpein is mutually exclusive with both others at the engine.
    assert " -cpin " not in script and " -cein " not in script


def test_run_script_unchanged_for_constant_ph(tmp_path):
    script = _write_run_script(tmp_path, None, "sys.cpin")
    assert "-cpin" in script and "-cpout" in script and "-cprestrt" in script
    assert "-cein" not in script and "-cpein" not in script


# ---------------------------------------------------------------------------
# The combined run: a cpin AND a cein together (Amber tutorial 33)
# ---------------------------------------------------------------------------

def _combined_workflow(tmp_path, config, steps=None):
    (tmp_path / "sys.cpin").write_text("&CNSTPH\n/\n")
    (tmp_path / "sys.cein").write_text("&CNSTE\n/\n")
    return WorkflowConfig(
        workflow_id="wf1", name="n", description="d",
        system_prmtop="sys.prmtop", initial_rst7="sys.rst7",
        steps=steps if steps is not None else [config],
        cpin_file=str(tmp_path / "sys.cpin"),
        titration_files={PH: str(tmp_path / "sys.cpin"),
                         REDOX: str(tmp_path / "sys.cein")},
    )


def test_setup_stages_and_flags_both_files(tmp_path):
    """The mp8 case: cpin for the propionates, cein for the heme, one run."""
    sim_dir = tmp_path / "prod"
    sim_dir.mkdir()
    m = _make_manager()
    config = _prod()
    workflow = _combined_workflow(tmp_path, config)
    m._get_workflow_for_step = lambda c: workflow

    result = m._setup_cpmd_files(config, sim_dir, tmp_path, "sys.prmtop",
                                "prod", silent=True)
    assert result["flags"] == [
        "-cpin", "sys.cpin", "-cpout", "prod.cpout",
        "-cprestrt", "prod.cprestrt",
        "-cein", "sys.cein", "-ceout", "prod.ceout",
        "-cerestrt", "prod.cerestrt",
    ]
    # Both files must actually reach the run directory.
    assert (sim_dir / "sys.cpin").exists()
    assert (sim_dir / "sys.cein").exists()


def test_each_file_chains_its_own_restart(tmp_path):
    """cpin chains from .cprestrt and cein from .cerestrt, independently."""
    run_dir = tmp_path / "run"
    prev_dir = run_dir / "prod1"
    prev_dir.mkdir(parents=True)
    (prev_dir / "prod1.cprestrt").write_text("ph restart")
    (prev_dir / "prod1.cerestrt").write_text("redox restart")
    cur_dir = run_dir / "prod2"
    cur_dir.mkdir()

    m = _make_manager()
    step1, step2 = _prod(name="prod1"), _prod(name="prod2")
    step2.workflow_step = 1
    workflow = _combined_workflow(tmp_path, step2, steps=[step1, step2])
    m._get_workflow_for_step = lambda c: workflow

    flags = m._setup_cpmd_files(step2, cur_dir, tmp_path, "sys.prmtop",
                                "prod2", silent=True)["flags"]
    assert flags[:2] == ["-cpin", "prod1.cprestrt"]
    assert flags[6:8] == ["-cein", "prod1.cerestrt"]


def test_run_script_emits_both_triples(tmp_path):
    m = _make_manager()
    steps = [_prod(name="prod1")]
    steps[0].workflow_step = 0
    steps[0].engine = "pmemd"
    steps[0].hardware_config = {}
    m.simulation_queue._workflows["wf1"] = _combined_workflow(
        tmp_path, steps[0], steps=steps)
    m._get_workflow_for_step = lambda c: m.simulation_queue._workflows["wf1"]
    m._is_production_step = lambda c: True
    m._get_template_type = lambda t: "production"

    workflow_dir = tmp_path / "wf"
    workflow_dir.mkdir()
    m._create_workflow_run_script(
        workflow_dir, steps, ["prod1/simulation.mdin"],
        "sys.prmtop", "sys.rst7", extended_production_cycles=0)
    script = (workflow_dir / "run_workflow.sh").read_text()

    for flag in ("-cpin", "-cpout", "-cprestrt",
                 "-cein", "-ceout", "-cerestrt"):
        assert flag in script, flag
    assert "-cpein" not in script


def test_combined_workflow_round_trips(tmp_path):
    m = _make_manager()
    m.simulation_queue.add_simulation(_prod())
    m.simulation_queue._workflows["wf1"] = _combined_workflow(
        tmp_path, _prod())
    m.simulation_queue._sync_to_workspace()

    stored = m.workspace.get("md_workflows")["wf1"]
    assert set(stored["titration_files"]) == {PH, REDOX}

    reloaded = SimulationQueue(workspace=m.workspace)
    _ = reloaded.queue
    assert set(reloaded._workflows["wf1"].active_titration_files()) == {PH, REDOX}


def test_legacy_single_file_workflow_still_resolves():
    """A workflow persisted before multi-file support has only cpin_file."""
    wf = WorkflowConfig(
        workflow_id="wf1", name="n", description="d",
        system_prmtop="sys.prmtop", initial_rst7="sys.rst7", steps=[],
        cpin_file="/tmp/sys.cpin",
    )
    assert wf.active_titration_files() == {PH: "/tmp/sys.cpin"}
