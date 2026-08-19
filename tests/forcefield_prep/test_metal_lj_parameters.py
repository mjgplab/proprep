"""
Metal Lennard-Jones terms: MCPB.py's table, and the right water model.

Two defects behind "there is no LJ parameter for Mo", reported from a 4UHX run
whose frcmod carried ``r=1.5, eps=0.01, VERIFY MANUALLY`` for molybdenum.

1. ProPrep looked only at the IOD sets, which cover ions with a measured
   ion-oxygen distance. Mo(VI) is not a free aqueous ion, so it has no IOD
   entry -- but MCPB.py's own IonLJParaDict carries it, taken from UFF:

       Mo6: (1.506, 0.056, 'Adopted from atom type Mo3+6/Mo6+6 from UFF
                            (Rappe et al. JACS, 114, 10024)')

   The placeholder radius was close; epsilon was low by 5.6x.

2. The database was constructed with ``water_model='tip3p'`` hardcoded at every
   call site, while that run had selected ``leaprc.water.opc``. LJ parameters
   are fitted per water model, so every metal in a non-TIP3P system -- not just
   the exotic ones -- got terms from the wrong set.

The leaprc mapping is on the exact suffix. The ad-hoc test it replaces,
``if "opc" in wm.lower()``, sent leaprc.water.opc3 to the OPC set.
"""

import logging

import pytest

from proprep.forcefield_prep.mcpb.metal_ion_database import (
    MetalIonDatabase, SUPPORTED_WATER_MODELS, water_model_from_leaprc,
)


pymsmt = pytest.importorskip(
    "pymsmt.mol.element", reason="AmberTools/pymsmt not installed")


# --------------------------------------------------------------------------- #
# leaprc -> LJ table
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("leaprc,expected", [
    ("leaprc.water.tip3p", "tip3p"),
    ("leaprc.water.opc", "opc"),
    ("leaprc.water.opc3", "opc3"),      # NOT "opc"
    ("leaprc.water.tip4pew", "tip4pew"),
    ("leaprc.water.spce", "spce"),
])
def test_supported_models_map_to_themselves(leaprc, expected):
    assert water_model_from_leaprc(leaprc) == expected


def test_opc3_is_not_swallowed_by_opc():
    """The substring bug: "opc" is inside "opc3"."""
    assert water_model_from_leaprc("leaprc.water.opc3") != "opc"


@pytest.mark.parametrize("leaprc,expected", [
    ("leaprc.water.opc3pol", "opc3"),   # polarizable OPC3 -> OPC3
    ("leaprc.water.spceb", "spce"),
    ("leaprc.water.tip5p", "tip3p"),    # no counterpart at all
    ("", "tip3p"),
    (None, "tip3p"),
])
def test_models_without_their_own_set_borrow_one(leaprc, expected):
    assert water_model_from_leaprc(leaprc) == expected


def test_borrowing_is_announced(caplog):
    """Silently substituting a water model would be the same class of bug."""
    with caplog.at_level(logging.WARNING):
        water_model_from_leaprc("leaprc.water.tip5p")

    assert "tip5p" in caplog.text and "tip3p" in caplog.text


def test_every_mapping_target_is_a_table_the_database_accepts():
    for leaprc in ("leaprc.water.opc3pol", "leaprc.water.spceb",
                   "leaprc.water.tip5p", "leaprc.water.opc"):
        assert water_model_from_leaprc(leaprc) in SUPPORTED_WATER_MODELS


# --------------------------------------------------------------------------- #
# Mo, and the rest of MCPB.py's table
# --------------------------------------------------------------------------- #

def test_molybdenum_resolves_instead_of_falling_back():
    """The reported symptom."""
    radius, epsilon, source = MetalIonDatabase(water_model='tip3p')._get_vdw_params('MO', 6)

    assert (radius, epsilon) == pytest.approx((1.506, 0.056))
    assert "VERIFY MANUALLY" not in source
    assert "UFF" in source


def test_the_source_carries_its_citation():
    _r, _e, source = MetalIonDatabase()._get_vdw_params('MO', 6)

    assert "Rappe" in source


@pytest.mark.parametrize("element", ["MO", "Mo", "mo"])
def test_the_element_symbol_is_matched_case_insensitively(element):
    """Elements arrive uppercase from PDB files and mixed-case elsewhere."""
    radius, _e, _s = MetalIonDatabase()._get_vdw_params(element, 6)

    assert radius == pytest.approx(1.506)


def test_an_element_in_no_table_still_falls_back():
    """The placeholder must survive for genuinely unparameterized metals."""
    radius, epsilon, source = MetalIonDatabase()._get_vdw_params('MO', 3)

    assert (radius, epsilon) == (1.5, 0.01)
    assert "VERIFY MANUALLY" in source


def test_an_empty_element_does_not_raise():
    radius, _e, source = MetalIonDatabase()._get_vdw_params('', 2)

    assert radius == 1.5 and "VERIFY MANUALLY" in source


# --------------------------------------------------------------------------- #
# the water model actually changes the numbers
# --------------------------------------------------------------------------- #

def test_iron_differs_between_tip3p_and_opc():
    """If these agreed, the hardcoded tip3p would have been harmless."""
    tip3p = MetalIonDatabase(water_model='tip3p')._get_vdw_params('FE', 3)
    opc = MetalIonDatabase(water_model='opc')._get_vdw_params('FE', 3)

    assert tip3p[0] != opc[0]
    assert tip3p[1] != opc[1]


def test_the_iod_sets_are_unchanged_by_the_new_tier():
    """MCPB.py's table is consulted after IOD, so IOD entries must still win."""
    db = MetalIonDatabase(water_model='tip3p')

    radius, epsilon, source = db._get_vdw_params('ZN', 2)

    assert (radius, epsilon) == pytest.approx((1.395, 0.014917))
    assert "IOD" in source


def test_uff_entries_do_not_vary_with_water_model():
    """Mo6 is UFF-derived, not fitted against water; it should be identical."""
    a = MetalIonDatabase(water_model='tip3p')._get_vdw_params('MO', 6)
    b = MetalIonDatabase(water_model='opc')._get_vdw_params('MO', 6)

    assert a == b


def test_an_unsupported_model_is_downgraded_at_construction():
    db = MetalIonDatabase(water_model='tip5p')

    assert db.water_model == 'tip3p'
