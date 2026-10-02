#!/usr/bin/env python3
"""
Molecular Docking module: registration, start-up without the chemistry
packages, and the whole menu driven with scripted answers.

The scripted run docks biotin back into streptavidin (PDB 1STP, unmodified
in tests/data/docking) the way a user would, in dashboard order: choose the
structure, take biotin as the ligand (which leaves the receptor) and make it
the carboxylate with two typed edits, keep the rest of the HETATM proposal,
protonate at pH 7, keep the torsions, box it, keep Vina's defaults, run, and look at the results. The top pose
must land within 2 A of the crystal.

Run with: pytest tests/test_docking_module.py
"""

import io
import os
import shutil
import subprocess
import sys

import pytest
from rich.console import Console

from proprep.application import menu_commands

SRC = os.path.join(os.path.dirname(__file__), "..", "src")
DATA = os.path.join(os.path.dirname(__file__), "data", "docking")


def test_the_module_is_listed_everywhere_a_module_must_be():
    with open(os.path.join(SRC, "proprep", "application", "pdbprocessor.py")) as handle:
        assert '"proprep.docking_prep.docking_module"' in handle.read()
    assert "Molecular Docking" in menu_commands.WorkflowMenuCommand.STAGE_TOOLS["simulate"]["tools"]
    with open(os.path.join(SRC, "proprep", "application", "menu_commands.py")) as handle:
        text = handle.read()
    assert text.count('"Molecular Docking",') == 2          # main-menu group and module order
    with open(os.path.join(SRC, "..", "proprep.spec")) as handle:
        assert "'proprep.docking_prep.docking_module'," in handle.read()


def test_proprep_starts_where_rdkit_and_meeko_are_missing():
    """The AmberTools system-Python build has no RDKit; importing the module must not fail
    (ProPrep would print an import warning at every start), and opening it says what is missing."""
    code = ("import sys\n"
            "for name in ('rdkit', 'gemmi', 'meeko', 'vina'):\n"
            "    sys.modules[name] = None\n"
            "from proprep.docking_prep.docking_module import MolecularDockingModule\n"
            "print(MolecularDockingModule.NAME)\n")
    run = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                         env={**os.environ, "PYTHONPATH": os.path.abspath(SRC)})
    assert run.returncode == 0, run.stderr
    assert run.stdout.strip() == "Molecular Docking"


# ------------------------------------------------------------ scripted run ---

class _Workspace(dict):
    def set(self, key, value):
        self[key] = value

    def has(self, key):
        return key in self


class _Processor:
    def __init__(self, folder):
        self.out = io.StringIO()
        self.console = Console(file=self.out, width=160)
        self.workspace = _Workspace(working_directory=str(folder))

    def _get_workspace(self):
        return self.workspace


def _scripted(monkeypatch, answers, asked):
    import proprep.utils.prompts as prompts

    def next_answer(processor, prompt, *args, **kwargs):
        asked.append(prompt)
        assert answers, f"no scripted answer left for: {prompt}"
        value = answers.pop(0)
        if value == "" and kwargs.get("default") is not None:
            return str(kwargs["default"])
        return value

    monkeypatch.setattr(prompts, "prompt_with_context", next_answer)
    monkeypatch.setattr(prompts, "prompt_float_with_retry",
                        lambda p, prompt, default, **k: float(next_answer(p, prompt, default=default) or default))
    monkeypatch.setattr(prompts, "prompt_int_with_retry",
                        lambda p, prompt, default, **k: int(next_answer(p, prompt, default=default) or default))
    monkeypatch.setattr(prompts, "confirm_with_context",
                        lambda p, prompt, default=None, **k: next_answer(p, prompt).lower().startswith("y"))


def test_whole_menu_redocks_biotin(monkeypatch, tmp_path):
    pytest.importorskip("meeko", reason="Meeko is required for docking")
    pytest.importorskip("vina", reason="AutoDock Vina is required for docking")
    if shutil.which("pdb2pqr") is None and not os.path.exists(os.path.join(os.path.dirname(sys.executable), "pdb2pqr")):
        pytest.skip("pdb2pqr is required for the protonation step")
    import proprep.utils.structure_selector as selector
    from proprep.docking_prep.docking_module import MolecularDockingModule

    pdb = tmp_path / "1STP.pdb"
    shutil.copy(os.path.join(DATA, "1STP.pdb.txt"), pdb)
    monkeypatch.setattr(selector.StructureSelector, "get_structure", lambda self, *a, **k: str(pdb))
    answers = [
        "1", "y", "all",                                   # receptor: this structure, every chain
        "2", "s", "1", "remove HO2", "charge O12 -1", "",   # biotin as ligand, made the carboxylate
                                                           # (one edit prompt: BTN's CCD entry is sound)
        "3", "",                                           # the ligand is already out of the receptor
        "4", "7.0", "",                                    # protonation at pH 7, no changes
        "7", "",                                           # rotatable bonds as proposed
        "8", "l", "8",                                     # box around the ligand, 8 A padding
        "9", "vina", "42", "8", "9", "3.0", "1.0", "0.375", "0",
        "r", "v", "0", "b",
    ]
    asked = []
    _scripted(monkeypatch, answers, asked)
    module = MolecularDockingModule()
    module.processor = _Processor(tmp_path)
    module.initialize()
    assert module.process(module.processor.workspace) is True
    assert answers == []

    workspace = module.processor.workspace
    state = workspace["docking_state"]
    assert state["receptor"]["keep_hetero"] == []
    assert state["ligand"]["ccd_edits"] == []
    assert state["ligand"]["edits"] == [{"op": "remove_atom", "atoms": ["HO2"]},
                                        {"op": "set_formal_charge", "atoms": ["O12"], "value": -1}]
    assert not [p for p in asked if p.startswith("Edit for BTN")]   # no CCD prompt for a sound entry
    results = workspace["docking_results"]
    assert results["rmsd"][0] < 2.0
    assert os.path.exists(workspace["docking_poses_sdf"])
    shown = module.processor.out.getvalue()
    assert "net charge -1; charged atoms: O12 -1" in shown
    assert "RMSD to crystal (A)" in shown

    # session replay matches prompt text exactly: no question may carry the current value
    assert not [p for p in asked if "now" in p.lower().split("(")[0] or "(now" in p.lower()]
    assert "Formal charge" not in " ".join(asked)          # no metals here, so no charge questions


def test_blocked_menu_entry_says_why(monkeypatch):
    """Other tools under a blocked ○ show '⚠ Needs ...'; the docking entry does too."""
    import proprep.utils.structure_selector as selector
    from proprep.docking_prep.docking_module import MolecularDockingModule
    module = MolecularDockingModule()
    monkeypatch.setattr(selector.StructureSelector, "get_structure", lambda self, *a, **k: None)
    assert module.availability_note({}) == "Needs a loaded structure"
    [option] = module.get_enhanced_menu_options({})
    assert option.status.name == "BLOCKED" and option.dependency_text == "Needs a loaded structure"
    monkeypatch.setattr(selector.StructureSelector, "get_structure", lambda self, *a, **k: "/x/prot.pdb")
    assert module.availability_note({}) is None


# ------------------------------------------------------------ viewer picks ---

class _PickModule:
    """The parts of the module the pick paths use."""
    def __init__(self, tmp_path, state):
        self.processor = _Processor(tmp_path)
        self.console = self.processor.console
        self.state = state
        self.saved = 0

    def save_state(self):
        self.saved += 1

    def show_ligand_in_viewer(self, build):
        pass


def _biotin(tmp_path):
    pytest.importorskip("rdkit")
    from proprep.docking_prep import ligand_sources as ls
    from proprep.docking_prep.ccd_chemistry import component_from_block
    with open(os.path.join(DATA, "ccd_BTN.cif.txt")) as handle:
        btn = component_from_block(handle.read(), "BTN")
    pdb = tmp_path / "1STP.pdb"
    shutil.copy(os.path.join(DATA, "1STP.pdb.txt"), pdb)
    return ls.from_structure_residue(str(pdb), "A", 300, "", "BTN", btn), str(pdb)


def _fake_pick(monkeypatch, results):
    from proprep.structure_prep.viewer_coordinator import viewer
    monkeypatch.setattr(viewer, "pick", lambda kind, prompt, structure_index=None, timeout=300.0: results.pop(0))


def test_a_picked_bond_or_atom_becomes_a_typed_edit(monkeypatch, tmp_path):
    """The viewer reports indices (an SDF has no atom names); they map back to the ligand's names."""
    from proprep.docking_prep.docking_menus_ligand import _picked_edit
    from proprep.docking_prep.docking_state import DockingState
    build, _ = _biotin(tmp_path)
    names = [a.GetProp("name") for a in build.mol.GetAtoms()]
    carbonyl = next((i, j) for i, j in ((b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in build.mol.GetBonds())
                    if build.mol.GetBondBetweenAtoms(i, j).GetBondTypeAsDouble() == 2.0)
    m = _PickModule(tmp_path, DockingState())
    _fake_pick(monkeypatch, [{"structure": 1, "atom1": {"index": carbonyl[0]}, "atom2": {"index": carbonyl[1]}},
                             {"structure": 1, "atom": {"index": carbonyl[1]}}])
    _scripted(monkeypatch, ["b", "1", "a", "c", "-1"], [])
    assert _picked_edit(m, build) == f"bond {names[carbonyl[0]]} {names[carbonyl[1]]} 1"
    assert _picked_edit(m, build) == f"charge {names[carbonyl[1]]} -1"


def test_no_viewer_means_no_pick_and_a_prompt_to_type(monkeypatch, tmp_path):
    from proprep.docking_prep.docking_menus_ligand import _picked_edit
    from proprep.docking_prep.docking_state import DockingState
    build, _ = _biotin(tmp_path)
    m = _PickModule(tmp_path, DockingState())
    _fake_pick(monkeypatch, [None])
    _scripted(monkeypatch, ["a"], [])
    assert _picked_edit(m, build) == ""
    assert "Type it instead" in m.processor.out.getvalue()


def test_a_picked_receptor_atom_adds_its_residue_to_the_flexible_ones(monkeypatch, tmp_path):
    from proprep.docking_prep.docking_menus_receptor import choose_flexible
    from proprep.docking_prep.docking_state import DockingState
    _, pdb = _biotin(tmp_path)
    state = DockingState()
    state.receptor.pdb_path, state.receptor.chains, state.receptor.flexible = pdb, ["A"], ["A:45"]
    m = _PickModule(tmp_path, state)
    _fake_pick(monkeypatch, [{"structure": 0, "atom": {"index": 1, "chain": "A", "resno": 88, "resname": "SER",
                                                       "inscode": ""}}])
    _scripted(monkeypatch, ["p", "n"], [])
    choose_flexible(m)
    assert state.receptor.flexible == ["A:45", "A:88"] and m.saved == 1


def test_a_picked_bond_switches_that_torsion(monkeypatch, tmp_path):
    pytest.importorskip("meeko", reason="Meeko is required for torsions")
    from proprep.docking_prep.docking_menus_ligand import choose_torsions
    from proprep.docking_prep.docking_state import DockingState
    from proprep.docking_prep.ligand_prep import torsion_bonds
    build, pdb = _biotin(tmp_path)
    bond = torsion_bonds(build.mol, rigid_macrocycles=False)[0]
    index = {a.GetProp("name"): a.GetIdx() for a in build.mol.GetAtoms()}
    state = DockingState()
    state.ligand.source = {"kind": "structure_residue"}
    m = _PickModule(tmp_path, state)
    m.build_ligand = lambda apply_edits_=True: build
    _fake_pick(monkeypatch, [{"structure": 1, "atom1": {"index": index[bond.atoms[0]]},
                              "atom2": {"index": index[bond.atoms[1]]}}])
    _scripted(monkeypatch, ["p", ""], [])
    choose_torsions(m)
    assert state.ligand.rigid_bonds == [sorted(bond.atoms)]


def test_option_letters_in_prompts_survive_rich_markup(monkeypatch):
    """Rich reads "[s]" as its strikethrough tag and printed "Ligand source:  residue of the
    structure,  SMILES,  file" in the running program; the fake prompts of the other tests never
    render, so this one renders through Rich itself."""
    import proprep.utils.prompts as prompts
    from proprep.docking_prep import docking_ui as ui
    rendered = []

    def render(processor, prompt, **kwargs):
        console = Console(file=io.StringIO(), width=200)
        console.print(prompt, end="")
        rendered.append(console.file.getvalue())
        return kwargs.get("default", "")
    monkeypatch.setattr(prompts, "prompt_with_context", render)
    ui.ask(None, "Ligand source: [s] residue of the structure, [m] SMILES, [f] file", default="s",
           description="test")
    assert rendered[-1] == "Ligand source: [s] residue of the structure, [m] SMILES, [f] file"
    table = ui.table("t", ["name"], [["mol_[s]"]])
    console = Console(file=io.StringIO(), width=80)
    console.print(table)
    assert "mol_[s]" in console.file.getvalue()


def test_the_settings_explain_the_scoring_functions_before_asking(monkeypatch, tmp_path):
    from proprep.docking_prep.docking_module import MolecularDockingModule
    from proprep.docking_prep import docking_menus_ligand as lig
    module = MolecularDockingModule()
    module.processor = _Processor(tmp_path)
    module.initialize()
    module.save_state = lambda: None
    asked = []
    _scripted(monkeypatch, ["vinardo", "42", "8", "9", "3.0", "1.0", "0.375", "0"], asked)
    lig.choose_settings(module)
    shown = module.processor.out.getvalue()
    for line in ("vina     AutoDock Vina's empirical function (Trott & Olson, 2010)",
                 "vinardo  A re-parameterisation of Vina's terms (Quiroga & Villarreal, 2016)",
                 "ad4      AutoDock 4's force-field function (Huey et al., 2007)",
                 "not comparable between functions"):
        assert line in shown
    assert shown.index("Scoring functions") < len(shown) and asked[0].startswith("Scoring function")
    assert module.state.settings["scoring"] == "vinardo"


def test_showing_a_pose_says_what_happened(monkeypatch, tmp_path):
    """After a run, pose 1 goes to an open viewer, or the terminal says how to see it; asked for
    from the results, a pose opens the viewer if none is open. Nothing is silent."""
    from unittest.mock import MagicMock
    import proprep.structure_prep.viewer_coordinator as vc
    from proprep.docking_prep.docking_module import MolecularDockingModule
    module = MolecularDockingModule()
    module.processor = _Processor(tmp_path)
    module.initialize()
    module.state.receptor.pdb_path = "receptor.pdb"
    module.state.last_run = {"files": {"pose_1": "pose_1.sdf", "pose_2": "pose_2.sdf"}}
    fake = MagicMock()
    monkeypatch.setattr(vc, "viewer", fake)
    monkeypatch.setattr(vc, "_is_web_shell_mode", lambda: False)

    fake.is_running.return_value = False
    module._show_pose(1)
    assert not fake.show_structures.called
    assert "The viewer is not open; 'v' shows any pose in it." in module.processor.out.getvalue()

    module._show_pose(2, asked=True)
    fake.show_structures.assert_called_with(["receptor.pdb", "pose_2.sdf"], force=True)
    assert "Viewer: pose 2 with the receptor and the search box (opened in your browser)." \
        in module.processor.out.getvalue()

    fake.is_running.return_value = True
    module._show_pose(1)
    fake.show_structures.assert_called_with(["receptor.pdb", "pose_1.sdf"], force=False)
    assert "Viewer: pose 1 with the receptor and the search box." in module.processor.out.getvalue()


def test_run_campaign_then_results_and_the_newer_results_are_offered_first(monkeypatch, tmp_path):
    """The dashboard reads r (one ligand), c (the same for every ligand of a library), then v; with
    both a run and a campaign, v names each, marks the newer, and offers it as the default."""
    from proprep.docking_prep.docking_module import MolecularDockingModule
    from proprep.docking_prep import docking_menus_batch
    module = MolecularDockingModule()
    module.processor = _Processor(tmp_path)
    module.initialize()
    module._dashboard()
    shown = module.processor.out.getvalue()
    assert shown.index("  r Run docking") < shown.index("  c Docking campaign: run docking for every ligand of") \
        < shown.index("  v View results: no run or campaign yet")

    module.state.last_run = {"folder": str(tmp_path / "docking" / "docking_20261002_083835"), "files": {}}
    module.state.last_campaign = {"folder": str(tmp_path / "docking" / "campaign_lib_20261002_091500")}
    assert module._newer(module.state.last_run["folder"], module.state.last_campaign["folder"]) == "c"
    assert "the run docking_20261002_083835 and the campaign campaign_lib_20261002_091500" in module._results_summary()
    seen = []
    monkeypatch.setattr(docking_menus_batch, "results", lambda m: seen.append("campaign"))
    asked = []
    _scripted(monkeypatch, [""], asked)
    module._view_results(docking_menus_batch)
    assert seen == ["campaign"]
    assert "The campaign: campaign_lib_20261002_091500  (newer)" in module.processor.out.getvalue()
