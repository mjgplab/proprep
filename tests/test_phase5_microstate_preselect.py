#!/usr/bin/env python3
"""
Phase 5 tests: microstate FF-set pre-selection helper.

_preferred_ff_set_for_combo reads workspace transformer_info and returns the set
a site's Stage-1 pH-treatment choice implies for an exact (transformer, redox,
spin) combo — gated on an actual ph_treatment so no-fork cofactors get no
pre-selection. (The dead, hardcoded-constph _generate_microstate_tleap_script is
removed; its absence is asserted too.)

Run with: pytest tests/test_phase5_microstate_preselect.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from proprep.tleap_prep.tleap_input_generator import TLeapInputGenerator  # noqa: E402


def _gen(transformer_info):
    gen = object.__new__(TLeapInputGenerator)
    gen.get_from_workspace = lambda k, d=None: transformer_info if k == "transformer_info" else d
    return gen


COMBO = ("heme_cys_axial_b_type", "oxidized", "high_spin")


def test_preferred_set_for_matching_fork_combo():
    gen = _gen([{
        "transformer_type": "heme_cys_axial_b_type",
        "redox_state": "oxidized", "spin_state": "high_spin",
        "ph_treatment": "fixed_pH", "forcefield_set": "Guberman_HTO_RESP_FixedpH",
    }])
    assert gen._preferred_ff_set_for_combo(*COMBO) == "Guberman_HTO_RESP_FixedpH"


def test_none_when_no_ph_treatment_choice():
    gen = _gen([{
        "transformer_type": "heme_cys_axial_b_type",
        "redox_state": "oxidized", "spin_state": "high_spin",
        "ph_treatment": None, "forcefield_set": "Guberman_HTO_RESP",
    }])
    assert gen._preferred_ff_set_for_combo(*COMBO) is None


def test_none_when_combo_does_not_match():
    gen = _gen([{
        "transformer_type": "heme_cys_axial_b_type",
        "redox_state": "reduced", "spin_state": "high_spin",
        "ph_treatment": "fixed_pH", "forcefield_set": "Guberman_HTR_RESP_FixedpH",
    }])
    assert gen._preferred_ff_set_for_combo(*COMBO) is None  # oxidized != reduced


def test_none_on_empty_transformer_info():
    assert _gen([])._preferred_ff_set_for_combo(*COMBO) is None


def test_dead_hardcoded_constph_function_removed():
    assert not hasattr(TLeapInputGenerator, "_generate_microstate_tleap_script")
