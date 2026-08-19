"""
Importing external parameters: declared atom types, and the right site.

Two defects found while importing an externally-obtained FAD set.

1. The wizard never asked about atom types. The rationale is stated in
   library_promotion's docstring -- "only metal sites can introduce new atom
   types; small molecules (GAFF) and modified AAs never do" -- and it holds for
   parameters ProPrep GENERATES, where antechamber assigns existing GAFF types.
   It does not hold for an IMPORT: those files came from a collaborator or a
   paper and may declare anything. A declared type with no addAtomTypes entry
   fails at tleap with nothing in the wizard having mentioned it.

2. After depositing, the transformer offer listed every detected site.
   Importing FAD offered the Fe2S2 and MoCo clusters, neither of which contains
   a FAD residue. A transformer is built against a site, so only a site holding
   the imported residue is a sensible target.
"""

import pytest

from proprep.forcefield_prep.library_promotion import (
    _known_atom_types, parse_frcmod_mass_types,
)
from proprep.forcefield_prep.forcefield_parameterizer import (
    _rank_sites_for_library,
)


# --------------------------------------------------------------------------- #
# reading a frcmod's MASS section
# --------------------------------------------------------------------------- #

def _frcmod(tmp_path, body):
    path = tmp_path / "x.frcmod"
    path.write_text(body)
    return path


def test_declared_types_are_read(tmp_path):
    path = _frcmod(tmp_path, "title\n\nMASS\nM1 55.85\nY7 32.06\n\nBOND\nM1-Y7 70.0 2.3\n")

    assert parse_frcmod_mass_types(path) == [("M1", 55.85), ("Y7", 32.06)]


def test_an_empty_mass_section_declares_nothing(tmp_path):
    """The common case for GAFF output -- and for the reported FAD set."""
    path = _frcmod(tmp_path, "FAD parameters\n\nMASS\n\nBOND\ncd-ne 381.80 1.414\n")

    assert parse_frcmod_mass_types(path) == []


def test_the_section_ends_at_the_blank_line(tmp_path):
    """BOND entries below must not be mistaken for types."""
    path = _frcmod(tmp_path, "t\n\nMASS\nM1 55.85\n\nBOND\ncd-ne 381.8 1.414\n")

    assert [t for t, _m in parse_frcmod_mass_types(path)] == ["M1"]


def test_a_missing_file_is_not_an_error(tmp_path):
    assert parse_frcmod_mass_types(tmp_path / "nope.frcmod") == []


def test_a_frcmod_with_no_mass_header_declares_nothing(tmp_path):
    assert parse_frcmod_mass_types(_frcmod(tmp_path, "t\n\nBOND\na-b 1.0 1.0\n")) == []


# --------------------------------------------------------------------------- #
# which types already exist
# --------------------------------------------------------------------------- #

known = pytest.mark.skipif(not _known_atom_types(),
                           reason="Amber parameter files not locatable")


@known
@pytest.mark.parametrize("atom_type", ["c3", "cd", "ne", "os", "p5", "HO", "N", "CT"])
def test_standard_types_are_recognised(atom_type):
    assert atom_type in _known_atom_types()


@known
@pytest.mark.parametrize("atom_type", ["2C", "3C"])
def test_types_declared_in_a_frcmod_are_recognised(atom_type):
    """
    ff19SB's 2C/3C live in frcmod.ff19SB, not parm19.dat. Scanning only the
    .dat files reported them as new and would have asked the user to declare
    types the force field already has.
    """
    assert atom_type in _known_atom_types()


@known
@pytest.mark.parametrize("atom_type", ["M1", "Y7", "YA", "YB"])
def test_mcpb_style_types_are_not_recognised(atom_type):
    """These genuinely need an addAtomTypes entry."""
    assert atom_type not in _known_atom_types()


@known
def test_only_frcmods_a_leaprc_sources_are_consulted():
    """
    Scanning every shipped frcmod is over-broad. frcmod.tumuc, a niche
    nucleic-acid set, declares YA as a hydrogen (mass 1.008) -- the same two
    letters ProPrep gives a metal-bound oxygen. Letting an unsourced file vouch
    for a type would suppress a prompt the user needs.
    """
    assert "YA" not in _known_atom_types()


# --------------------------------------------------------------------------- #
# identifying the site an imported library belongs to
# --------------------------------------------------------------------------- #

def _site(site_id, residues):
    """residues: {resname: [atom_name, ...]}"""
    from types import SimpleNamespace
    atoms = [SimpleNamespace(resname=r, atom_name=a)
             for r, names in residues.items() for a in names]
    return SimpleNamespace(site_id=site_id, atoms=atoms)


FAD_ATOMS = ["PA", "O1A", "O2A", "O5B", "C5B", "C4B", "N1", "C2", "O2", "N5"]

SITES = [
    _site("site_1", {"CYM": ["N", "CA", "CB", "SG"], "FES": ["FE1", "S1"]}),
    _site("site_2", {"MOS": ["MO", "S", "O1"]}),
    _site("site_3", {"FAD": FAD_ATOMS}),
]
LIB = {a.upper() for a in FAD_ATOMS}


def test_a_name_match_ranks_first():
    ranked = _rank_sites_for_library(SITES, "FAD", LIB)

    assert ranked[0][0].site_id == "site_3"
    assert ranked[0][2] == 1.0


def test_a_renamed_library_still_finds_its_residue():
    """
    THE case a transformer exists for. Parameters for an oxidized FAD named
    FAO, against a structure that says FAD: matching on the residue NAME finds
    nothing, because renaming FAD -> FAO is the transformer's whole purpose.
    Atom names survive the rename.
    """
    ranked = _rank_sites_for_library(SITES, "FAO", LIB)

    assert ranked[0][0].site_id == "site_3"
    assert ranked[0][1] == "FAD"
    assert ranked[0][2] == pytest.approx(1.0)


def test_unrelated_sites_score_low():
    ranked = _rank_sites_for_library(SITES, "FAO", LIB)
    scores = {site.site_id: score for site, _r, score in ranked}

    assert scores.get("site_1", 0.0) < 0.8
    assert scores.get("site_2", 0.0) < 0.8


def test_a_partial_overlap_is_reported_proportionally():
    """A residue sharing half its names is a weak match, not a binary miss."""
    half = _site("s", {"XYZ": ["PA", "O1A", "ZZ1", "ZZ2"]})

    ranked = _rank_sites_for_library([half], "FAO", LIB)

    assert ranked[0][2] == pytest.approx(0.5)


def test_no_library_atoms_falls_back_to_the_name():
    """An unreadable lib must not make every site look like a match."""
    ranked = _rank_sites_for_library(SITES, "FAD", set())

    assert [s.site_id for s, _r, _sc in ranked] == ["site_3"]


def test_nothing_matches_when_the_residue_is_absent():
    metal_only = SITES[:2]

    ranked = _rank_sites_for_library(metal_only, "FAO", LIB)

    assert all(score < 0.8 for _s, _r, score in ranked)


def test_hydrogens_are_excluded_from_the_comparison(tmp_path):
    """
    A crystal structure usually has none, so counting them would depress every
    comparison. The reported FAD library declared 84 atoms, 53 non-hydrogen,
    and the structure's residue had exactly those 53.
    """
    from proprep.forcefield_prep.library_promotion import library_atom_names

    lib = tmp_path / "x.lib"
    lib.write_text(
        '!entry.X.unit.atoms table  str name  str type\n'
        ' "PA" "p5" 0 1 0 1 15 0.0\n'
        ' "H1" "ho" 0 1 0 2 1 0.0\n'
        '!entry.X.unit.atomspertinfo table\n')

    assert library_atom_names(lib) == {"PA"}


# --------------------------------------------------------------------------- #
# what the wizard tells the user to do next
# --------------------------------------------------------------------------- #

def _offer(monkeypatch, sites, result):
    """Run the post-import offer with prompts declined, capturing the output."""
    from rich.console import Console
    from types import SimpleNamespace
    from proprep.forcefield_prep import forcefield_parameterizer as fp

    monkeypatch.setattr(fp, "confirm_with_context", lambda *a, **k: False)

    workspace = SimpleNamespace(
        get=lambda key, default=None: sites if key == "detected_redox_sites" else default)
    parameterizer = fp.ForcefieldParameterizer.__new__(fp.ForcefieldParameterizer)
    parameterizer.console = Console(record=True, width=110)
    parameterizer.processor = SimpleNamespace(_get_workspace=lambda: workspace)
    parameterizer._offer_transformer_for_import(result)
    return parameterizer.console.export_text()


REAL_RESULT = {
    "library_path": "/u/.proprep/forcefield_params/specialized_residues/small_molecules/FAD",
    "state_dir": "/u/.proprep/forcefield_params/specialized_residues/"
                 "small_molecules/FAD/single_state/default",
}


def test_the_link_target_is_the_entry_not_the_state_directory(monkeypatch):
    """
    FORCEFIELD_PATH names the entry; the state is carried separately. Existing
    transformers show this: FORCEFIELD_PATH = 'metal_sites/4hux_fe2s2'.
    """
    from proprep.forcefield_prep.forcefield_parameterizer import (
        ForcefieldParameterizer,
    )

    seed = ForcefieldParameterizer._import_forcefield_seed(REAL_RESULT)

    assert seed["path"] == "small_molecules/FAD"
    assert seed["redox_state"] == "single_state"
    assert seed["spin_state"] == "default"


def test_the_message_describes_binding_not_only_renaming(monkeypatch):
    """
    It used to say a library is inert "until something renames its residues",
    which makes an already-correctly-named cofactor look like it needs nothing.
    Binding is the part that is always required.
    """
    text = _offer(monkeypatch, [], REAL_RESULT)

    assert "binds" in text
    assert "pass-through" in text


def test_it_names_the_pass_through_route_when_no_site_matches(monkeypatch):
    """The reported situation: sites exist, none of them is the imported residue."""
    text = _offer(monkeypatch,
                  [_site("site_1", {"CYM": ["N", "CA"], "FES": ["FE1"]})],
                  REAL_RESULT)

    assert "No detected site resembles these parameters" in text
    assert "Redox Site Detector" in text
    assert "save" in text and "pass-through" in text


def test_it_still_names_the_library_to_link(monkeypatch):
    text = _offer(monkeypatch, [], REAL_RESULT)

    assert "small_molecules/FAD" in text
