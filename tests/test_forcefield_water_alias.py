#!/usr/bin/env python3
"""Regression tests for bare LEaP unit-alias handling (e.g. ``HOH = TP3``).

Crystallographic waters are named ``HOH`` in the PDB, but no AMBER ``.lib``
defines a ``HOH`` unit. tleap makes ``HOH`` usable at runtime via a variable
assignment in ``leaprc.water.*`` -- ``HOH = TP3`` binds the alias name to the
water-model UNIT. The Force Field Explorer must replay that alias to type or
browse waters. These tests cover the two halves of that:

1. ``parse_leaprc`` captures the bare ``NAME = UNIT`` form (unit test, hermetic).
2. ``ForceFieldData`` registers the alias only when the target unit loaded, and
   the aliased name then both types and browses (integration test, needs
   AMBERHOME).
"""

import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from proprep.forcefield_prep.ff_parsers import parse_leaprc


# ---------------------------------------------------------------------------
# 1. Parser unit tests (hermetic -- no AMBER install required)
# ---------------------------------------------------------------------------

def _parse_text(tmp_path, text):
    p = tmp_path / "leaprc.test"
    p.write_text(text)
    return parse_leaprc(p)


def test_parse_leaprc_captures_bare_unit_alias(tmp_path):
    contents = _parse_text(tmp_path, "loadOff solvents.lib\nHOH = TP3\nWAT = TP3\n")
    assert contents.unit_aliases == {"HOH": "TP3", "WAT": "TP3"}


def test_parse_leaprc_unit_alias_ignores_comments_and_calls(tmp_path):
    text = (
        "# HOH = SHOULD_NOT_MATCH\n"      # commented out
        "HIS = HIE\n"                     # real alias
        "x = createUnit foo\n"            # function-call RHS (has a space) -> skip
        "WAT = TP3   # inline comment\n"  # real alias with trailing comment
    )
    contents = _parse_text(tmp_path, text)
    assert contents.unit_aliases == {"HIS": "HIE", "WAT": "TP3"}


def test_parse_leaprc_addpdbresmap_not_treated_as_unit_alias(tmp_path):
    # The multi-line addPdbResMap block must not leak into unit_aliases.
    text = 'addPdbResMap {\n  { 0 "HID" "NHID" }\n}\nHOH = TP3\n'
    contents = _parse_text(tmp_path, text)
    assert contents.unit_aliases == {"HOH": "TP3"}
    assert contents.pdb_res_map == {"HID": "NHID"}


# ---------------------------------------------------------------------------
# 2. End-to-end against a real AMBER install
# ---------------------------------------------------------------------------

requires_amber = pytest.mark.skipif(
    not os.environ.get("AMBERHOME"),
    reason="AMBERHOME not set; skipping integration test against real leaprc data",
)


@requires_amber
def test_hoh_types_and_browses_with_water_model():
    from proprep.forcefield_prep.forcefield_data import ForceFieldData

    ff = ForceFieldData()
    ff.load_leaprc("leaprc.protein.ff14SB")
    ff.load_leaprc("leaprc.water.tip3p")

    # HOH aliases onto the TIP3P unit and its oxygen types as OW.
    assert ff.residue_aliases.get("HOH") == "TP3"
    atom_def = ff.get_atom_definition("HOH", "O")
    assert atom_def is not None
    assert atom_def.atom_type == "OW"
    assert atom_def.charge == pytest.approx(-0.834)

    # add_residue_alias copies the atoms, so HOH is browsable too.
    assert "HOH" in ff.get_available_residues()
    assert ff.has_residue("HOH")


@requires_amber
def test_hoh_untyped_without_a_water_model():
    """Gating: a protein-only load has no water unit, so HOH must not resolve."""
    from proprep.forcefield_prep.forcefield_data import ForceFieldData

    ff = ForceFieldData()
    ff.load_leaprc("leaprc.protein.ff14SB")

    assert ff.residue_aliases.get("HOH") is None
    assert ff.get_atom_definition("HOH", "O") is None


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
