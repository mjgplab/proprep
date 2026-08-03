#!/usr/bin/env python3
"""
Regression tests for water-model -> divalent+ ion-set recommendation pairing.

The Topology Generator used to always recommend the 12-6-4 OPC ion set (the
catalog's static default) regardless of the water model the user picked. This
produced a contradictory menu: pick TIP3P water and the DIVALENT+ ION
PARAMETERS menu would still recommend an OPC-calibrated ion set. The pairing
now lives in the catalog (`recommended_ions_for_water`) and is consumed by
`_select_standard_forcefields_interactive` to re-mark the ions menu, mirroring
the protein-FF -> water-model logic in test_water_model_recommendation.py.

These tests lock in the rules:
  - OPC       -> 12-6-4 OPC (most accurate)   (prefers the 12-6-4 variant)
  - TIP3P     -> 12-6-4 TIP3P
  - SPC/E     -> 12-6-4 SPC/E
  - TIP4P-Ew  -> 12-6-4 TIP4P-Ew
  - OPC3      -> 12-6 OPC3                     (no 12-6-4 set exists for OPC3)
plus the fallback to "Default only" for water models with no dedicated Li/Merz
set (polarizable / 5-point / Force-Balance / SPC/Eb) and for no water selected,
and the invariant that every recommended ion name actually exists as an option.

Run with: pytest tests/test_ion_parameter_recommendation.py
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from proprep.forcefield_params.forcefield_catalog import (  # noqa: E402
    FORCEFIELD_OPTIONS,
    recommended_ions_for_water,
)


# Water names exactly as they appear in FORCEFIELD_OPTIONS['water'] -> the ion
# option name that should be recommended.
PREFER_1264_RULES = [
    ("OPC", "12-6-4 OPC (most accurate)"),
    ("TIP3P", "12-6-4 TIP3P"),
    ("SPC/E", "12-6-4 SPC/E"),
    ("TIP4P-Ew", "12-6-4 TIP4P-Ew"),
]

# Water models with a Li/Merz set but no 12-6-4 variant fall to the 12-6 set.
NON_1264_RULES = [
    ("OPC3", "12-6 OPC3"),
]

# Water models with no dedicated divalent+ set at all -> "Default only".
FALLBACK_WATERS = ["SPC/Eb", "TIP5P", "OPC3-pol", "TIP3P-FB", "TIP4P-FB"]


@pytest.mark.parametrize("water,expected", PREFER_1264_RULES + NON_1264_RULES)
def test_recommended_ions_for_water(water, expected):
    assert recommended_ions_for_water(water) == expected


@pytest.mark.parametrize("water", FALLBACK_WATERS)
def test_water_without_limerz_set_falls_back_to_default(water):
    assert recommended_ions_for_water(water) == "Default only"


def test_no_water_selected_falls_back_to_default():
    """Implicit solvent (water None / empty) -> Default only, never the static OPC."""
    assert recommended_ions_for_water(None) == "Default only"
    assert recommended_ions_for_water("") == "Default only"


def test_recommendation_never_mismatches_water():
    """The headline bug: a TIP3P (or any non-OPC) choice must never resolve to
    an OPC-tagged ion set. Verify the recommended set's for_water tag actually
    contains the chosen water (or is the tagless Default-only fallback)."""
    ion_by_name = {o["name"]: o for o in FORCEFIELD_OPTIONS["ions"]["options"]}
    for water, _ in PREFER_1264_RULES + NON_1264_RULES:
        rec = recommended_ions_for_water(water)
        opt = ion_by_name[rec]
        assert water in opt.get("for_water", []), (
            f"{water} resolved to {rec!r}, whose for_water is "
            f"{opt.get('for_water')!r} — mismatch"
        )


def test_every_recommended_ion_exists_in_menu():
    """A typo'd target would silently mark nothing in the UI — guard it."""
    ion_names = {o["name"] for o in FORCEFIELD_OPTIONS["ions"]["options"]}
    all_waters = [w["name"] for w in FORCEFIELD_OPTIONS["water"]["options"]]
    for water in all_waters + [None]:
        rec = recommended_ions_for_water(water)
        assert rec in ion_names, f"{rec!r} (for {water!r}) not an ion menu option"


def test_default_only_is_a_real_option():
    """The fallback target must exist so callers can mark it recommended."""
    ion_names = {o["name"] for o in FORCEFIELD_OPTIONS["ions"]["options"]}
    assert "Default only" in ion_names


def test_every_for_water_tag_names_a_real_water_model():
    """Each ion option's for_water entries must match actual water menu names,
    or the recommendation lookup can never fire for them."""
    water_names = {w["name"] for w in FORCEFIELD_OPTIONS["water"]["options"]}
    for opt in FORCEFIELD_OPTIONS["ions"]["options"]:
        for w in opt.get("for_water", []):
            assert w in water_names, (
                f"ion option {opt['name']!r} tags unknown water {w!r}"
            )
