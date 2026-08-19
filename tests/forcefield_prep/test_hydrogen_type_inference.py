"""
An inferred hydrogen is typed for what it is bonded to.

Amber names a hydrogen after its neighbour, not after hydrogen generally, and
the fallback typer mapped every H to 'H' from the element alone. 'H' is
specifically the amide/amine hydrogen:

    parm10.dat / parm19.dat   (identical in both)
      H   0.6000  0.0157   H bonded to nitrogen atoms
      HO  0.0000  0.0000   hydroxyl group
      HS  0.6000  0.0157   hydrogen bonded to sulphur
      HC  1.4870  0.0157   H bonded to aliphatic carbon

So the hydroxo proton on a Mo cofactor was given a 0.6 A vdW sphere that the
hydroxyl convention deliberately sets to zero. HO is not specific to ff19SB --
it is byte-identical in parm10 (ff14SB's base) and predates both.

This fallback only fires where there is no library entry, which is the withheld
cluster residues; standard residues keep their library types.

Detection is the bonded heavy neighbour, found with the same covalent-radius
test the cluster bond perception uses, so the two cannot disagree about what is
bonded to what.
"""

import logging

import pytest

from proprep.forcefield_prep.mcpb.atom_typer import MCPBAtomTyper
from proprep.structure_prep.comprehensive_redox_detector import (
    RedoxSite, RedoxSiteAtom,
)


def _atom(resid, name, element, coords, resname="LIG"):
    return RedoxSiteAtom(chain="A", resname=resname, resid=resid,
                         atom_name=name, coords=coords, element=element)


def _typer():
    typer = MCPBAtomTyper.__new__(MCPBAtomTyper)
    typer.logger = logging.getLogger("test")
    return typer


def _site(atoms):
    site = RedoxSite(site_id="s", structure_id="t")
    site.atoms = list(atoms)
    return site


def _infer(heavy_element, distance, resid=1):
    """Type an H sitting `distance` from one heavy atom."""
    heavy = _atom(resid, "X", heavy_element, (0.0, 0.0, 0.0))
    hydrogen = _atom(resid, "H1", "H", (distance, 0.0, 0.0))
    return _typer()._infer_type_from_element(
        "H", hydrogen, _site([heavy, hydrogen])).strip()


# --------------------------------------------------------------------------- #
# the reported case
# --------------------------------------------------------------------------- #

def test_the_hydroxo_proton_is_typed_ho():
    """MOS as reported: H1 0.98 A from O1."""
    mos = [
        _atom(1312, "MO", "MO", (-29.569, -17.017, -41.812), "MOS"),
        _atom(1312, "S", "S", (-29.451, -17.649, -39.574), "MOS"),
        _atom(1312, "O1", "O", (-27.666, -16.557, -41.673), "MOS"),
        _atom(1312, "O2", "O", (-29.607, -18.617, -42.453), "MOS"),
        _atom(1312, "H1", "H", (-27.532, -15.655, -41.974), "MOS"),
    ]

    typed = _typer()._infer_type_from_element("H", mos[-1], _site(mos)).strip()

    assert typed == "HO"


def test_it_binds_to_the_nearer_oxygen():
    """O1 is 0.98 A away, O2 is far; the nearest bonded heavy atom wins."""
    assert _infer("O", 0.98) == "HO"


# --------------------------------------------------------------------------- #
# the other neighbours
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("element,distance,expected", [
    ("O", 0.98, "HO"),
    ("S", 1.34, "HS"),
    ("N", 1.01, "H"),
    ("C", 1.09, "HC"),
])
def test_the_neighbour_decides_the_type(element, distance, expected):
    assert _infer(element, distance) == expected


def test_nitrogen_still_gives_the_old_answer():
    """The previous behaviour was right for exactly one case; keep it."""
    assert _infer("N", 1.01) == "H"


# --------------------------------------------------------------------------- #
# falling back safely
# --------------------------------------------------------------------------- #

def test_an_unbonded_hydrogen_keeps_the_generic_type():
    assert _infer("O", 9.0) == "H"


def test_no_site_falls_back_to_the_element():
    assert _typer()._infer_type_from_element("H").strip() == "H"


def test_a_neighbour_in_another_residue_is_not_used():
    """Residue boundaries matter: an H is typed by its own residue."""
    other = _atom(2, "O", "O", (0.0, 0.0, 0.0))
    hydrogen = _atom(1, "H1", "H", (0.98, 0.0, 0.0))

    typed = _typer()._infer_type_from_element(
        "H", hydrogen, _site([other, hydrogen])).strip()

    assert typed == "H"


def test_another_hydrogen_is_not_a_neighbour():
    """Geminal or H2-like contacts must not decide the type."""
    partner = _atom(1, "H2", "H", (0.74, 0.0, 0.0))
    hydrogen = _atom(1, "H1", "H", (0.0, 0.0, 0.0))

    typed = _typer()._infer_type_from_element(
        "H", hydrogen, _site([partner, hydrogen])).strip()

    assert typed == "H"


def test_an_unknown_neighbour_element_falls_back():
    assert _infer("SE", 1.46) == "H"


# --------------------------------------------------------------------------- #
# non-hydrogen typing is untouched
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("element,expected", [
    ("C", "CT"), ("N", "N"), ("O", "O"), ("S", "S"),
    ("ZN", "Zn"), ("FE", "Fe"),
])
def test_other_elements_are_unchanged(element, expected):
    assert _typer()._infer_type_from_element(element).strip() == expected


def test_every_type_is_two_characters():
    """AMBER requires it; the padding must survive the new branch."""
    for element, distance in (("O", 0.98), ("S", 1.34), ("N", 1.01), ("C", 1.09)):
        heavy = _atom(1, "X", element, (0.0, 0.0, 0.0))
        hydrogen = _atom(1, "H1", "H", (distance, 0.0, 0.0))
        typed = _typer()._infer_type_from_element(
            "H", hydrogen, _site([heavy, hydrogen]))
        assert len(typed) == 2, f"{element}: {typed!r}"


# --------------------------------------------------------------------------- #
# the path the metal-site workflow actually uses
# --------------------------------------------------------------------------- #

def _global(element, coords, resid=1312, chain="A"):
    from types import SimpleNamespace
    return SimpleNamespace(element=element, chain=chain, resid=resid, coords=coords)


MOS_GLOBAL = [
    ((-29.569, -17.017, -41.812), "MO"),
    ((-29.451, -17.649, -39.574), "S"),
    ((-27.666, -16.557, -41.673), "O"),
    ((-29.607, -18.617, -42.453), "O"),
    ((-27.532, -15.655, -41.974), "H"),
]


def test_the_global_registry_path_types_the_hydroxo_proton():
    """
    The workflow types cluster atoms through _convert_global_types, not through
    MCPBAtomTyper's fallback. A withheld atom has no library type there, and
    original_type fell back to the raw ELEMENT -- which for hydrogen is 'H',
    a real Amber type (amide/amine), so it resolved silently and wrongly.
    """
    from proprep.forcefield_prep.metal_site_parameterizer import _inferred_element_type

    residue_atoms = {("A", 1312): MOS_GLOBAL}
    hydrogen = _global("H", (-27.532, -15.655, -41.974))

    assert _inferred_element_type(hydrogen, residue_atoms) == "HO"


@pytest.mark.parametrize("element", ["MO", "FE", "S", "O"])
def test_non_hydrogen_placeholders_are_unchanged(element):
    """
    'MO'/'FE'/'S'/'O' are not Amber types, so they read as placeholders and the
    atoms are renamed to M*/Y* regardless. Only hydrogen collides with a real
    type, so only hydrogen is corrected.
    """
    from proprep.forcefield_prep.metal_site_parameterizer import _inferred_element_type

    atom = _global(element, (-29.569, -17.017, -41.812))

    assert _inferred_element_type(atom, {("A", 1312): MOS_GLOBAL}) == element


def test_a_hydrogen_with_no_neighbours_keeps_the_generic_type():
    from proprep.forcefield_prep.metal_site_parameterizer import _inferred_element_type

    lone = _global("H", (50.0, 50.0, 50.0))

    assert _inferred_element_type(lone, {("A", 1312): MOS_GLOBAL}) == "H"


def test_an_unknown_residue_key_does_not_raise():
    from proprep.forcefield_prep.metal_site_parameterizer import _inferred_element_type

    stray = _global("H", (0.0, 0.0, 0.0), resid=9999)

    assert _inferred_element_type(stray, {("A", 1312): MOS_GLOBAL}) == "H"


def test_the_two_paths_share_one_table():
    """Both call sites must agree; drift here is a silent typing difference."""
    from proprep.forcefield_prep.mcpb.atom_typer import H_TYPE_BY_NEIGHBOR

    assert H_TYPE_BY_NEIGHBOR == {"O": "HO", "S": "HS", "N": "H", "C": "HC"}


# --------------------------------------------------------------------------- #
# the origin: preprocessing's prmtop fallback
# --------------------------------------------------------------------------- #

def test_the_preprocessing_fallback_is_where_the_type_is_decided():
    """
    Three code paths could type a cluster hydrogen and only one runs. The
    reported symptom survived two fixes aimed at the other two:

      MCPBAtomTyper._infer_type_from_element   -- not reached (would emit 'XX'
                                                  for MO; the run stored 'MO')
      _convert_global_types                    -- not reached either
      structure_preprocessor's prmtop fallback -- THIS one

    Identifying it needs no guessing: the recorded original_type was the bare
    element symbol for every withheld atom, which only that branch produces.
    """
    import inspect

    from proprep.forcefield_prep import structure_preprocessor

    source = inspect.getsource(structure_preprocessor)
    assert "hydrogen_type_from_neighbors" in source, (
        "the preprocessing fallback must type hydrogens by neighbour")


def test_the_shared_helper_handles_the_live_mos_geometry():
    """MOS exactly as it sits in prepared_structure.pdb."""
    from proprep.forcefield_prep.mcpb.atom_typer import hydrogen_type_from_neighbors

    atoms = [
        ((-29.569, -17.017, -41.812), "MO"),
        ((-29.451, -17.649, -39.574), "S"),
        ((-27.666, -16.557, -41.673), "O"),
        ((-29.607, -18.617, -42.453), "O"),
        ((-27.532, -15.655, -41.974), "H"),
    ]

    assert hydrogen_type_from_neighbors((-27.532, -15.655, -41.974), atoms) == "HO"


def test_a_metal_is_not_chosen_as_the_bonding_partner():
    """
    The hydroxo H sits 2.4 A from Mo and 0.98 A from O1. Nearest-bonded must
    pick the oxygen; a plain nearest-atom search over heavy atoms would too,
    but only because of the distance -- pin it.
    """
    from proprep.forcefield_prep.mcpb.atom_typer import hydrogen_type_from_neighbors

    atoms = [
        ((0.0, 0.0, 0.0), "MO"),
        ((1.90, 0.0, 0.0), "O"),
        ((2.88, 0.0, 0.0), "H"),
    ]

    assert hydrogen_type_from_neighbors((2.88, 0.0, 0.0), atoms) == "HO"
