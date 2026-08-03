#!/usr/bin/env python3
"""
Regression tests for protein-FF -> water-model recommendation pairing.

The Topology Generator used to always recommend OPC water regardless of the
protein force field the user picked. The pairing now lives in the catalog
(`recommended_water_for_protein`) and is consumed by
`_select_standard_forcefields_interactive` to re-mark the WATER MODEL menu.

These tests lock in the documented rules:
  - ff14SB           -> TIP3P
  - ff19SB           -> OPC
  - constant pH      -> TIP3P
  - constant redox   -> TIP3P
  - constant pH+redox -> TIP3P
plus the secondary pairings (ff15ipq->SPC/Eb, fb15->TIP3P-FB, legacy->TIP3P),
the unknown-FF fallback (None), and the invariant that every recommended water
name actually exists as an option in the water menu.

Run with: pytest tests/test_water_model_recommendation.py
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from proprep.forcefield_params.forcefield_catalog import (  # noqa: E402
    FORCEFIELD_OPTIONS,
    PROTEIN_WATER_RECOMMENDATION,
    recommended_water_for_protein,
)


# The headline rules the user specified, keyed by the protein option's leaprc
# exactly as it appears in FORCEFIELD_OPTIONS['protein'].
PRIMARY_RULES = [
    ("leaprc.protein.ff14SB", "TIP3P"),
    ("leaprc.protein.ff19SB", "OPC"),
    ("leaprc.constph", "TIP3P"),                       # constant pH
    ("leaprc.conste", "TIP3P"),                        # constant redox
    (["leaprc.constph", "leaprc.conste"], "TIP3P"),    # constant pH + redox
]

SECONDARY_RULES = [
    ("leaprc.protein.ff15ipq", "SPC/Eb"),
    ("leaprc.protein.fb15", "TIP3P-FB"),
    ("leaprc.protein.ff03.r1", "TIP3P"),
    ("oldff/leaprc.ff99SBildn", "TIP3P"),
]


@pytest.mark.parametrize("leaprc,expected", PRIMARY_RULES + SECONDARY_RULES)
def test_recommended_water_for_protein(leaprc, expected):
    assert recommended_water_for_protein(leaprc) == expected


def test_unknown_protein_ff_has_no_recommendation():
    """Unknown FF -> None so the caller keeps the catalog's default (OPC)."""
    assert recommended_water_for_protein("leaprc.protein.does_not_exist") is None


def test_none_input_is_safe():
    assert recommended_water_for_protein(None) is None
    assert recommended_water_for_protein([]) is None


def test_combined_list_first_known_wins():
    """A list resolves on the first constituent with a known pairing."""
    # Unknown leaprc ahead of conste should be skipped, not block resolution.
    assert recommended_water_for_protein(
        ["leaprc.protein.unknown", "leaprc.conste"]
    ) == "TIP3P"


def test_every_recommended_water_exists_in_menu():
    """A typo'd target would silently mark nothing in the UI — guard it."""
    water_names = {o["name"] for o in FORCEFIELD_OPTIONS["water"]["options"]}
    for target in set(PROTEIN_WATER_RECOMMENDATION.values()):
        assert target in water_names, f"{target!r} not a water menu option"


def test_every_catalog_protein_ff_is_mapped():
    """Every protein FF the user can pick should have an explicit pairing, so
    the WATER MODEL menu is never left on the static OPC default by omission."""
    for opt in FORCEFIELD_OPTIONS["protein"]["options"]:
        leaprc = opt["leaprc"]
        assert recommended_water_for_protein(leaprc) is not None, (
            f"protein option {opt['name']!r} ({leaprc!r}) has no water pairing"
        )
