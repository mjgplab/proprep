#!/usr/bin/env python3
"""Tests for tLEaP/PDB-safe residue-name validation (proprep.utils.tleap_utils).

Background: tLEaP matches ``loadpdb`` residues to templates by the OFF/lib
entry name (not the residue name inside a mol2), lexes a digit-leading bare
token like ``9E2`` as scientific notation, and PDB ``resName`` is only three
columns wide. So a residue name that must survive into a reusable library and
match a structure residue has to be 1-3 chars, letter-leading, alphanumeric.

These pin the rules behind rejecting names like ``9E2`` at the small-molecule
parameterizer's naming prompt.

Run with: pytest tests/utils/test_tleap_resname_validation.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from proprep.utils.tleap_utils import (
    is_tleap_safe_resname,
    suggest_tleap_safe_resname,
    tleap_safe_unit_var,
)


class TestIsTleapSafeResname:
    def test_accepts_plain_letter_leading(self):
        assert is_tleap_safe_resname("LIG")
        assert is_tleap_safe_resname("lig")
        assert is_tleap_safe_resname("A1")
        assert is_tleap_safe_resname("E92")   # letter-leading digits are fine
        assert is_tleap_safe_resname("Cl1")

    def test_rejects_digit_leading(self):
        # The 9E2 / 0G6 / 1N7 CCD-code family that started this whole saga.
        assert not is_tleap_safe_resname("9E2")
        assert not is_tleap_safe_resname("0G6")
        assert not is_tleap_safe_resname("1N7")

    def test_rejects_too_long(self):
        # PDB resName is 3 columns; a 4-char name truncates on read.
        assert not is_tleap_safe_resname("ABCD")
        assert not is_tleap_safe_resname("m9e2")

    def test_rejects_empty_and_nonalnum(self):
        assert not is_tleap_safe_resname("")
        assert not is_tleap_safe_resname(None)  # type: ignore[arg-type]
        assert not is_tleap_safe_resname("A-1")
        assert not is_tleap_safe_resname("A B")


class TestSuggestTleapSafeResname:
    def test_digit_leading_gets_letter_prefix_and_truncates(self):
        assert suggest_tleap_safe_resname("9E2") == "X9E"
        assert suggest_tleap_safe_resname("0G6") == "X0G"

    def test_strips_nonalnum_and_uppercases(self):
        assert suggest_tleap_safe_resname("lig-1") == "LIG"

    def test_empty_falls_back_to_lig(self):
        assert suggest_tleap_safe_resname("") == "LIG"
        assert suggest_tleap_safe_resname("---") == "LIG"

    def test_suggestions_are_themselves_safe(self):
        for name in ["9E2", "0G6", "1N7", "lig-1", "", "ABCD", "m9e2"]:
            assert is_tleap_safe_resname(suggest_tleap_safe_resname(name))


class TestUnitVarVsResnameContract:
    """tleap_safe_unit_var is a throwaway-handle helper; it does NOT produce a
    name safe to persist as a lib entry (it can be 4+ chars / mangled)."""

    def test_mangled_unit_var_is_not_a_safe_resname(self):
        # 9E2 -> m9E2 is fine as a script variable but 4 chars: unusable as a
        # PDB/lib residue name. This is exactly the trap the naming prompt avoids.
        assert tleap_safe_unit_var("9E2") == "m9E2"
        assert not is_tleap_safe_resname(tleap_safe_unit_var("9E2"))

    def test_already_safe_name_passes_through_unit_var(self):
        assert tleap_safe_unit_var("LIG") == "LIG"
