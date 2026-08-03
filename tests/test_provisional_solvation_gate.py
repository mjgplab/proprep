"""Regression: the Topology Generator must clear a PROVISIONAL (preprocessing-
origin) implicit solvation so the user is prompted for a real solvation choice
on the production build — but must NOT do so while preprocessing runs its own
metal-free tleap (that build is intentionally box-less), nor for membrane
systems, nor when the user genuinely chose implicit.

This guards both the original bug (box-less production topology) and its first
regression (being asked about solvent during the MCPB preprocessing tleap).
"""

from proprep.tleap_prep.tleap_input_generator import TLeapInputGenerator


def _gen(preprocessing_active=False):
    g = TLeapInputGenerator.__new__(TLeapInputGenerator)
    flag = {"_preprocessing_tleap_active": preprocessing_active}
    g.get_from_workspace = lambda k, d=None: flag.get(k, d)
    return g


PROVISIONAL = {"solvent_model": "implicit", "provisional": True}
USER_IMPLICIT = {"solvent_model": "implicit"}
USER_EXPLICIT = {"solvent_model": "explicit", "use_octahedron": True}


def test_production_run_clears_provisional():
    # User-initiated production build inherited the placeholder -> prompt.
    assert _gen()._provisional_solvation_needs_prompt(PROVISIONAL, False) is True


def test_preprocessing_run_keeps_provisional():
    # Preprocessing's own metal-free tleap must stay box-less and silent.
    assert _gen(preprocessing_active=True)._provisional_solvation_needs_prompt(PROVISIONAL, False) is False


def test_membrane_never_prompts():
    assert _gen()._provisional_solvation_needs_prompt(PROVISIONAL, True) is False


def test_genuine_user_choices_are_respected():
    assert _gen()._provisional_solvation_needs_prompt(USER_IMPLICIT, False) is False
    assert _gen()._provisional_solvation_needs_prompt(USER_EXPLICIT, False) is False


def test_missing_or_nondict_solvation():
    assert _gen()._provisional_solvation_needs_prompt(None, False) is False
    assert _gen()._provisional_solvation_needs_prompt({}, False) is False
