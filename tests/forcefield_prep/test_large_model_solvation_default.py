"""
The large model's solvation is asked for, not inherited.

The prompt for implicit solvation on the large model sat entirely under
`large_charge < 0`. Two things followed:

* an anionic large model silently adopted the small model's SCRF, announcing it
  only after the fact ("applying same solvation for consistency") — the choice
  was never offered;
* a neutral or cationic large model dropped the small model's SCRF entirely,
  with no message at all, leaving the two calculations at different levels of
  theory for no stated reason.

Solvation is now asked in every case. This pins what the answer defaults to.
"""

import pytest

from proprep.forcefield_prep.metal_site_parameterizer import MetalSiteWorkflowManager

_default = MetalSiteWorkflowManager._default_large_model_solvation

WATER = "SCRF=(Solvent=Water)"


@pytest.mark.parametrize("charge", [-3, -2, -1, 0, 1, 2])
def test_matches_the_small_model_whenever_it_was_solvated(charge):
    """Diverging from the small model should take a deliberate 'no'."""
    assert _default(WATER, charge) is True


@pytest.mark.parametrize("charge", [-1, -2, -3])
def test_unsolvated_small_model_still_recommends_it_for_an_anion(charge):
    assert _default("", charge) is True


@pytest.mark.parametrize("charge", [0, 1, 2])
def test_unsolvated_small_model_and_no_anion_defaults_to_no(charge):
    assert _default("", charge) is False


def test_the_moco_case():
    """4UHX MoCo: charge -3, and the small model was run in water."""
    assert _default(WATER, -3) is True


def test_a_cationic_site_still_follows_the_small_model():
    """The case the old anion-only gate dropped silently."""
    assert _default(WATER, 2) is True
    assert _default("", 2) is False
