#!/usr/bin/env python3
"""Integration smoke tests for the pb_titrate ↔ ProPrep wiring.

No PBSA required — these verify the pieces that don't need an external
calculator:

  [1] Workflow surface: 8 steps declared, every handler resolves, the
      executor can be constructed and the entry point is callable.
  [2] sites.discover_sites maps cleanly between BPTI prmtop and a
      detected_redox_sites-like list (already covered in detail by
      envelope_retention/run_test.py; spot-check here).
  [3] Topology Generator menu now includes 'pb_titrate' between
      'generate_topology' and 'generate_cpin'; status logic gates on
      has_topology and titrate_recommendations correctly.
  [4] cpin step's option-3 ('Use PB Titrate recommendations') maps a
      titrate_recommendations dict to the same {resnum: state_id} shape
      that options 1 and 2 produce, with a sensible fallback for sites
      not in the recommendations.

Tests that *do* require PBSA (BPTI/HEWL/RNase pKa regression, bundle
smoke, ion rebalance) live in their own directories and are documented
in the project plan as the 'requires AmberTools' tier.
"""
from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import parmed as pmd

HERE = Path(__file__).resolve().parent


def test_workflow_surface():
    """Every WorkflowStep handler resolves on PBTitrateWorkflow."""
    from proprep.pb_titrate.workflow import (
        PB_TITRATE_STEPS, PBTitrateWorkflow, run_pb_titrate_workflow,
    )
    assert len(PB_TITRATE_STEPS) == 8, (
        f"Expected 8 steps, got {len(PB_TITRATE_STEPS)}")
    expected_ids = ["pbt-1", "pbt-2", "pbt-3", "pbt-4",
                    "pbt-5", "pbt-6", "pbt-7", "pbt-8"]
    actual_ids = [s.id for s in PB_TITRATE_STEPS]
    assert actual_ids == expected_ids, (
        f"Step IDs out of order: {actual_ids}")
    for step in PB_TITRATE_STEPS:
        assert hasattr(PBTitrateWorkflow, step.handler), (
            f"Handler {step.handler} missing on PBTitrateWorkflow")
    assert callable(run_pb_titrate_workflow)
    print(f"  ✓ 8 workflow steps, all handlers resolve, entry point callable")


def test_topology_generator_menu_has_pb_titrate():
    """The Topology Generator now lists 'pb_titrate' as a menu option."""
    from proprep.tleap_prep.tleap_input_generator import TLeapInputGenerator
    from proprep.tleap_prep.tleap_commands import RunPBTitrateCommand

    inst = TLeapInputGenerator()
    options = inst.get_menu_options()
    assert "pb_titrate" in options, (
        f"pb_titrate not in menu options: {list(options.keys())}")
    # And the new command class is importable.
    assert RunPBTitrateCommand is not None
    print(f"  ✓ menu options include pb_titrate ({len(options)} total)")


def test_enhanced_menu_status_logic():
    """get_enhanced_menu_options gates pb_titrate on topology, marks
    completed when titrate_recommendations is populated."""
    from proprep.tleap_prep.tleap_input_generator import TLeapInputGenerator
    from proprep.utils.enhanced_menu import OptionStatus

    inst = TLeapInputGenerator()

    # Case A: no topology → pb_titrate BLOCKED
    ws_a = MagicMock()
    ws_a.get = MagicMock(side_effect=lambda k, default=None: {
        "transformed_pdb_file": None,
        "generated_microstate_pdbs": None,
        "tleap_input_file": None,
        "generated_microstate_tleap_files": None,
        "output_dir": str(HERE),
        "cpin_file": None,
        "titrate_recommendations": None,
    }.get(k, default))

    opts = inst.get_enhanced_menu_options(ws_a)
    pbt = next(o for o in opts if o.description.startswith("Refine"))
    assert pbt.status == OptionStatus.BLOCKED, (
        f"Expected BLOCKED with no topology, got {pbt.status}")

    # Case B: topology present, recs populated → pb_titrate COMPLETED.
    # The method imports glob locally, so we patch glob.glob globally.
    import glob as _glob
    ws_b = MagicMock()
    ws_b.get = MagicMock(side_effect=lambda k, default=None: {
        "transformed_pdb_file": "x.pdb",
        "generated_microstate_pdbs": None,
        "tleap_input_file": "x.tleap",
        "generated_microstate_tleap_files": None,
        "output_dir": str(HERE),
        "cpin_file": None,
        "titrate_recommendations": {("AS4", 3): {"state_id": 1}},
    }.get(k, default))
    with patch.object(_glob, "glob", return_value=["fake.prmtop"]):
        opts = inst.get_enhanced_menu_options(ws_b)
    pbt = next(o for o in opts if o.description.startswith("Refine"))
    assert pbt.status == OptionStatus.COMPLETED, (
        f"Expected COMPLETED when recs present, got {pbt.status}")
    print(f"  ✓ pb_titrate menu status: BLOCKED without topology, "
          f"COMPLETED when recs present")


def test_cpin_option3_state_assignment():
    """Option 3 of _set_initial_protonation_states pulls state_id from
    titrate_recommendations and falls back to default (state 0) for
    selected residues that have no recommendation."""
    from proprep.tleap_prep.tleap_input_generator import TLeapInputGenerator

    inst = TLeapInputGenerator()
    inst.processor = MagicMock()
    inst.processor.console = MagicMock()
    inst.processor.workspace = MagicMock()

    # Mock workspace returns a recommendation dict covering 2 of 3 residues.
    titrate_recs = {
        ("AS4",  3): {"state_id": 1, "state_name": "PROT",
                       "prot_count": 1, "pka_corr": 4.0, "net_charge": 0.0},
        ("LYS", 15): {"state_id": 0, "state_name": "PROT",
                       "prot_count": 3, "pka_corr": 10.4, "net_charge": 1.0},
        # No entry for HIP-42 — should fall back to state 0.
    }
    inst.processor.workspace.get = MagicMock(
        side_effect=lambda k, default=None: (
            titrate_recs if k == "titrate_recommendations" else default))

    selected_residues = [
        {"resname": "AS4", "resnum": 3,  "chain": "A", "pka": 4.0},
        {"resname": "LYS", "resnum": 15, "chain": "A", "pka": 10.4},
        {"resname": "HIP", "resnum": 42, "chain": "A", "pka": 6.6},
    ]

    # Force the prompt to return "3" (option 3: titrate recs).
    with patch("proprep.tleap_prep.tleap_input_generator.prompt_with_context",
                return_value="3"):
        result = inst._set_initial_protonation_states(selected_residues, step=1)

    assert result is not None, "Option 3 returned None unexpectedly"
    # AS4-3 → state 1 (from recommendation)
    assert result[3] == 1, f"AS4-3 expected state 1, got {result.get(3)}"
    # LYS-15 → state 0 (from recommendation)
    assert result[15] == 0, f"LYS-15 expected state 0, got {result.get(15)}"
    # HIP-42 → state 0 (default fallback, no recommendation)
    assert result[42] == 0, f"HIP-42 expected state 0 (default), got {result.get(42)}"
    print(f"  ✓ option 3 maps recs correctly: "
          f"AS4-3→1 (rec), LYS-15→0 (rec), HIP-42→0 (fallback)")


def test_cpin_option3_only_offered_when_recs_exist():
    """When titrate_recommendations is empty, option 3 is not offered;
    options 1 and 2 are unchanged."""
    from proprep.tleap_prep.tleap_input_generator import TLeapInputGenerator

    inst = TLeapInputGenerator()
    inst.processor = MagicMock()
    inst.processor.console = MagicMock()
    inst.processor.workspace = MagicMock()
    # No recommendations
    inst.processor.workspace.get = MagicMock(return_value=None)

    captured_choices = []

    def fake_prompt(processor, prompt, **kwargs):
        if "Select option" in prompt:
            captured_choices.append(kwargs.get("choices"))
        return "1"  # cpinutil defaults

    with patch("proprep.tleap_prep.tleap_input_generator.prompt_with_context",
                side_effect=fake_prompt):
        result = inst._set_initial_protonation_states([], step=1)

    assert captured_choices, "No 'Select option' prompt issued"
    assert captured_choices[0] == ["1", "2"], (
        f"Expected only options 1, 2 when no recs; got {captured_choices[0]}")
    assert result is None  # option 1 returns None
    print(f"  ✓ option 3 hidden when no recommendations")


def test_sites_discover_with_redox_envelope_smoke():
    """sites.discover_sites with a synthetic RedoxSite-like object groups
    the titratable residue into a multi-residue envelope (full coverage
    in envelope_retention/run_test.py)."""
    from proprep.pb_titrate.sites import discover_sites
    PRMTOP = HERE.parent / "bpti_asp3" / "bpti.prmtop"
    RST7   = HERE.parent / "bpti_asp3" / "bpti.rst7"
    s = pmd.load_file(str(PRMTOP), str(RST7))

    rs = SimpleNamespace(
        site_id="fake",
        residue_groups={(str(getattr(s.residues[i], "chain", "") or ""),
                          int(s.residues[i].number) + 1,
                          ""): []
                         for i in [1, 2, 3, 4]},
    )
    sites_list = discover_sites(s, [rs])
    as4 = next(x for x in sites_list if x.resname == "AS4")
    assert as4.is_multi_residue
    assert as4.envelope_idxs == {1, 2, 3, 4}
    print(f"  ✓ AS4-3 wraps into envelope {sorted(as4.envelope_idxs)} "
          f"from synthetic redox site")


def main():
    print("=== pb_titrate integration smoke tests ===\n")
    print("[1] Workflow surface")
    test_workflow_surface()
    print("\n[2] Topology Generator menu has pb_titrate option")
    test_topology_generator_menu_has_pb_titrate()
    print("\n[3] Enhanced menu status logic")
    test_enhanced_menu_status_logic()
    print("\n[4] cpin option 3 — state assignment from titrate_recommendations")
    test_cpin_option3_state_assignment()
    print("\n[5] cpin option 3 — hidden when no recommendations exist")
    test_cpin_option3_only_offered_when_recs_exist()
    print("\n[6] sites.discover_sites with synthetic redox envelope")
    test_sites_discover_with_redox_envelope_smoke()
    print("\n=== ALL INTEGRATION TESTS PASSED ===")


if __name__ == "__main__":
    main()
