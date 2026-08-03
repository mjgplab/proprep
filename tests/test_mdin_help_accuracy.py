"""The MD wizard's parameter help must not drift back from the Amber manual.

An audit of this database against the manual (chapters 23 sander / 24 pmemd; see
docs/mdin_help_audit.md) found ten parameters carrying outright errors. None was
a typo -- each was fluent, plausible prose that happened to be wrong, which is
exactly the kind of thing that survives review. These tests pin the corrected
facts so a future edit has to argue with the manual rather than with nothing.

They deliberately assert the *fact*, not the wording: defaults, option sets and
namelist assignments, plus the absence of the specific invented claims.
"""

import re
from pathlib import Path

import pytest

from proprep.md_prep.amber_parameter_database import (
    Parameter,
    build_parameter_database,
)

DB = build_parameter_database()


# ── Schema ──────────────────────────────────────────────────────────────

def test_advisory_fields_are_gone():
    """why_default / when_to_change asked questions the manual does not answer.

    They could only be filled with unsourceable advice, so the fields were
    removed rather than emptied -- an empty field is an invitation to refill it.
    """
    fields = set(Parameter.__dataclass_fields__)
    assert "why_default" not in fields
    assert "when_to_change" not in fields
    assert "manual_notes" in fields


def test_no_advisory_fields_in_the_source():
    src = Path(__file__).parent.parent / "src/proprep/md_prep/amber_parameter_database.py"
    text = src.read_text()
    body = text.split('"""', 2)[-1]          # skip the module docstring, which explains them
    assert "why_default=" not in body
    assert "when_to_change=" not in body


# ── Defaults the manual states explicitly ───────────────────────────────

@pytest.mark.parametrize("name,expected", [
    ("ew_type", 0),                  # "Standard use is to have EW_TYPE = 0"
    ("nsnb", 25),                    # "Default is 25."
    ("mdout_flush_interval", 300),   # "a default minimum interval of 300 seconds"
    ("baroscalingdir", 0),           # "= 0 box size scales randomly ... (default)"
    ("ntp", 0),
    ("ntc", 1),
    ("ntf", 1),
    ("imin", 0),
    ("dt", 0.001),
    ("temp0", 300.0),
    ("comp", 44.6),
    ("vlimit", 20.0),
    ("order", 4),
    ("tol", 1e-5),
    ("dsum_tol", 1e-5),
    ("drms", 1e-4),
])
def test_default_matches_manual(name, expected):
    assert DB[name].default == expected


def test_ntwr_default_is_nstlim_not_zero():
    """The manual says "Default = nstlim." Showing 0 implies restarts are off."""
    assert DB["ntwr"].default == "nstlim"


# ── Namelist assignment ─────────────────────────────────────────────────

@pytest.mark.parametrize("name,namelist", [
    ("ew_type", "ewald"),            # introduced under "The &ewald namelist has..."
    ("fft_grids_per_ang", "ewald"),  # manual: "In &ewald."
    ("dsum_tol", "ewald"),
    ("order", "ewald"),
    ("vdw_cutoff", "cntrl"),         # manual: "In &cntrl, these variables..."
    ("es_cutoff", "cntrl"),
    ("mdout_flush_interval", "cntrl"),
])
def test_namelist_matches_manual(name, namelist):
    assert DB[name].namelist == namelist


# ── Option sets ─────────────────────────────────────────────────────────

def test_ew_type_options_are_not_inverted():
    """The database had 0='Regular Ewald (slow)' and 1='PME (fast)' -- backwards."""
    opts = DB["ew_type"].options
    assert "mesh" in opts[0].lower() or "pme" in opts[0].lower()
    assert "exact" in opts[1].lower()


@pytest.mark.parametrize("name,expected", [
    ("ntf", {1, 2, 3, 4, 5, 6, 7, 8}),
    ("ivcap", {0, 1, 2, 5}),
    ("baroscalingdir", {0, 1, 2, 3}),
    ("ntt", {0, 1, 2, 3, 9, 10, 11}),
    ("ntmin", {0, 1, 2, 3, 4, 5}),
    ("imin", {0, 1, 5, 6, 7}),
    ("csurften", {0, 1, 2, 3}),
])
def test_option_set_matches_manual(name, expected):
    assert set(DB[name].options) == expected


def test_ithermostat_has_no_undocumented_option():
    """The manual: "Two types of thermostats are currently available" -- 1 and 2.

    An option 0 ("No thermostat in middle scheme") had been invented; it was the
    only fabricated option value in the audit.
    """
    assert set(DB["ithermostat"].options) == {1, 2}


def test_every_parameter_can_express_its_own_default():
    """ivcap and baroscalingdir both omitted option 0 -- which was the default.

    A few defaults are genuinely conditional (lj1264 follows the prmtop; ntwr
    defaults to nstlim), so no option value can be "the default". Those are
    exempt from the membership check but must say so in manual_notes -- the
    exemption is not a licence to leave the default unexplained.
    """
    offenders, conditional = [], []
    for p in DB.values():
        if not p.options:
            continue
        if p.default in p.options:
            continue
        (conditional if isinstance(p.default, str) else offenders).append(p.name)

    assert offenders == [], f"default not selectable: {offenders}"
    for name in conditional:
        assert DB[name].manual_notes, (
            f"{name} has a conditional default and must explain it in manual_notes"
        )


# ── Specific claims the manual contradicts ──────────────────────────────

def test_middle_scheme_is_leapfrog_not_velocity_verlet():
    """Manual: "=1 'middle' scheme based on the leapfrog algorithm"."""
    text = (DB["ischeme"].what_it_does + " " + " ".join(DB["ischeme"].options.values())).lower()
    assert "leapfrog" in text
    assert "verlet" not in text


def test_ithermostat_points_at_therm_par_not_gamma_ln():
    """The middle scheme's coupling comes from therm_par."""
    text = " ".join(DB["ithermostat"].options.values()).lower()
    assert "therm_par" in text
    assert "gamma_ln" not in text


def test_ithermostat_option_2_is_andersen():
    """It had been labelled "Velocity rescaling (Bussi-like)"; the manual says Andersen."""
    assert "andersen" in DB["ithermostat"].options[2].lower()


def test_idistr_is_a_frequency_not_a_distribution_type():
    """Manual: "the frequency at which the thermostat velocity distribution
    functions are accumulated" -- not a choice of distribution."""
    p = DB["idistr"]
    assert "frequency" in (p.brief + p.what_it_does).lower()
    assert p.options is None, "0=Uniform / 1=Gaussian were invented"


def test_mdout_flush_interval_is_in_seconds():
    """The help said "in steps"; the manual says integer seconds."""
    assert "second" in DB["mdout_flush_interval"].what_it_does.lower()


def test_baroscalingdir_is_not_described_as_a_membrane_setting():
    p = DB["baroscalingdir"]
    text = (p.brief + " " + p.what_it_does).lower()
    assert "membrane" not in text
    assert "barostat" in text


# ── Invented specifics ──────────────────────────────────────────────────

@pytest.mark.parametrize("name,fragment", [
    ("therm_par", "2-50"),      # invented range; the manual gives no fixed range
    ("ioutfm", "50%"),          # invented figure; the manual says only "smaller"
    ("ntc", "10 fs"),           # invented vibration period
    ("jfastw", "SETTLE"),       # the manual never names SETTLE
])
def test_invented_specific_is_absent(name, fragment):
    p = DB[name]
    haystack = " ".join(filter(None, [p.brief, p.what_it_does, p.manual_notes]))
    assert fragment.lower() not in haystack.lower()


# ── manual_notes hygiene ────────────────────────────────────────────────

def test_manual_notes_exist_where_the_manual_gives_guidance():
    """These are parameters whose manual entry carries explicit guidance."""
    for name in ["dt", "vlimit", "taup", "gamma_ln", "tol", "comp", "fswitch",
                 "ntb", "ntr", "ncyc", "dielc", "vdw_cutoff", "es_cutoff"]:
        assert DB[name].manual_notes, f"{name} should carry manual guidance"


def test_pressure_keywords_the_manual_documents_are_present():
    """Both were missing from the database entirely.

    baro_stochastic matters because the manual says it "improves on the Berendsen
    barostat and produces the correct isothermal-isobaric ensemble" -- i.e. it is
    the fix for the known flaw in the default barostat. ninterface is required
    alongside csurften.
    """
    assert DB["baro_stochastic"].default == 0
    assert set(DB["baro_stochastic"].options) == {0, 1}
    assert DB["ninterface"].default == 2


def test_keywords_the_wizard_referenced_but_did_not_offer():
    """The help pointed at parameters a user could not reach.

    ntr's guidance says to set netfrc=0; nsnb's says it is only consulted when
    nbflag=0; fswitch's says it is incompatible with the 12-6-4 models. All three
    referents were absent from the database.
    """
    for name in ["netfrc", "nbflag", "skinnb", "lj1264", "plj1264"]:
        assert name in DB, f"{name} is documented and referenced but not offered"


def test_ewald_namelist_assignment_for_the_added_keywords():
    """netfrc, nbflag and skinnb are defined under "The &ewald namelist has the
    following variables" (section 23.7.2), not in &cntrl."""
    for name in ["netfrc", "nbflag", "skinnb"]:
        assert DB[name].namelist == "ewald"
    for name in ["lj1264", "plj1264"]:
        assert DB[name].namelist == "cntrl"


def test_skinnb_default_is_two_angstrom():
    assert DB["skinnb"].default == 2.0


def test_netfrc_records_the_restraint_manager_seam():
    """Positional restraints are owned by the MD Manager's restraint manager, not
    this wizard. Anyone reading netfrc's note should know that."""
    assert "restraint manager" in DB["netfrc"].manual_notes.lower()


def test_ivcap_option_2_inactivates_the_cap():
    """It had been labelled "Orthorhombic virtual box", which is not what it does."""
    text = DB["ivcap"].options[2].lower()
    assert "inactivat" in text
    assert "orthorhombic" not in text


def test_surften_default_is_verified_against_the_manual():
    """The sander/pmemd excerpt did not cover it; the full manual states it."""
    assert DB["surften"].default == 0.005
    assert DB["surften"].manual_notes


def test_manual_notes_are_absent_where_the_manual_is_silent():
    """A note that could not be sourced is worse than no note.

    These parameters get no explicit guidance in the supplied chapters, so the
    field must stay empty rather than be filled to look complete.
    """
    for name in ["restraint_wt", "cutcap", "fcap", "nfft2", "nfft3"]:
        assert not DB[name].manual_notes, f"{name} has no manual guidance to cite"


# ── Feature completeness ────────────────────────────────────────────────

def test_water_cap_can_actually_be_positioned():
    """ivcap=1 needs cutcap AND a centre. The centre was not offered at all."""
    for name in ["xcap", "ycap", "zcap"]:
        assert name in DB
    assert DB["ivcap"] and DB["cutcap"] and DB["fcap"]


def test_electric_field_is_not_stuck_static():
    """efx/efy/efz were offered without efn, effreq or efphase, so the field
    could only ever be a static, absolutely-scaled one."""
    for name in ["efn", "effreq", "efphase"]:
        assert name in DB


def test_flush_interval_pair_is_complete_and_distinct():
    """mdout defaults to 300 s, mdinfo to 60 s. Only mdout was offered."""
    assert DB["mdout_flush_interval"].default == 300
    assert DB["mdinfo_flush_interval"].default == 60


def test_ewald_tolerance_pair_is_complete():
    """dsum_tol governs the direct sum, rsum_tol the reciprocal sum."""
    assert DB["dsum_tol"].default == 1e-5
    assert DB["rsum_tol"].default == 5e-5
    assert DB["rsum_tol"].namelist == "ewald"


# ── Second help surface: inline mdin comments ───────────────────────────
# amber_annotated_templates writes these comments INTO the generated mdin, so a
# wrong one persists in the user's own input file rather than just on screen.

def _inline_comment(param):
    from proprep.md_prep.amber_annotated_templates import AmberAnnotatedTemplate
    t = AmberAnnotatedTemplate.__new__(AmberAnnotatedTemplate)
    return AmberAnnotatedTemplate._get_param_comment(t, param)


def test_ntmin_inline_comment_is_not_inverted():
    """It read "1=steepest descent only, 2=steepest+conjugate" -- exactly backwards.

    Manual: "=1 For NCYC cycles the steepest descent method is used then conjugate
    gradient is switched on (default). =2 Only the steepest descent method is used."
    """
    c = _inline_comment("ntmin").lower()
    assert "1=steepest descent then conjugate" in c
    assert "2=steepest descent only" in c


@pytest.mark.parametrize("param,fragment", [
    ("imin", "0=md"),
    ("ntx", "1=coordinates only"),
    ("irest", "0=new simulation"),
    ("ntxo", "2=netcdf"),
    ("ioutfm", "1=netcdf"),
    ("ntt", "3=langevin"),
    ("ntb", "2=constant pressure"),
    ("barostat", "2=monte carlo"),
    ("ntc", "2=bonds with h"),
    ("nmropt", "1=&wt blocks enabled"),
])
def test_inline_comments_that_check_out_against_the_manual(param, fragment):
    assert fragment in _inline_comment(param).lower()
