"""
Coordinate-keyed charge lookups must survive a float32/float64 round trip.

Reported as "the small model's suggested QM charge is 0". Every atom in the
metal site printed 0.0000 in the charge breakdown, including CYM side chains
that plainly had charges in atom_type_assignments:

    CYM 114 CB  charge=-0.2413
    CYM 114 SG  charge=-0.8843999989024437

The lookup key is the atom's coordinate tuple. BioPython hands back float32;
assignments that have been through JSON -- a resumed session, a saved
atom_type_assignments.json -- come back as float64. They are not equal and do
not hash alike:

    tuple(atom.get_coord())  ==  (-41.416, -16.455, -50.088)   float32
    key from JSON            ==  (-41.416, -16.455, -50.088)   float64
    equal?  False

numpy prints the shortest repr that round-trips in float32, so the two display
identically while comparing unequal. Nothing looked wrong in the output either:
a row of 0.0000 is exactly what a withheld cluster atom legitimately shows.

Consequence: the writer summed 0.0 over all 44 site atoms and proposed 0 for
[Fe2S2(SCys)4]2-, whose charge is -2.
"""

import pytest

from proprep.forcefield_prep.pdb_writer import PDBWriter


np = pytest.importorskip("numpy")


def _writer():
    """A PDBWriter with only the state the lookup touches (no PDB parse)."""
    w = PDBWriter.__new__(PDBWriter)
    w._assignment_index = None
    w._assignment_index_src = None
    return w


# The reported site, as float64 keys (post-JSON) with real ff14SB charges.
SITE = {
    (-41.416, -16.455, -50.088): {"atom_name": "N", "charge": -0.4157},
    (-40.945, -17.763, -50.578): {"atom_name": "CA", "charge": -0.0351},
    (-41.874, -18.430, -51.585): {"atom_name": "CB", "charge": -0.2413},
    (-42.093, -17.443, -53.093): {"atom_name": "SG", "charge": -0.8844},
}


def _as_float32(coords):
    """What BioPython's atom.get_coord() yields."""
    return tuple(np.array(coords, dtype=np.float32))


# --------------------------------------------------------------------------- #
# the defect
# --------------------------------------------------------------------------- #

def test_float32_and_float64_keys_are_not_equal():
    """Pins the premise -- if this ever became true the fallback is dead code."""
    exact = (-41.416, -16.455, -50.088)
    from_biopython = _as_float32(exact)

    assert repr(from_biopython[0]) == repr(exact[0])      # they LOOK identical
    assert from_biopython != exact                        # but are not equal
    assert exact not in {from_biopython: 1}


def test_a_float32_coordinate_finds_a_float64_assignment():
    writer = _writer()

    for coords, expected in SITE.items():
        charge = writer._get_charge_from_coords(_as_float32(coords), SITE)

        assert charge == pytest.approx(expected["charge"])


def test_the_reported_total_is_recovered():
    """0.0000 for every site atom was the symptom; the sum is what matters."""
    writer = _writer()

    total = sum(writer._get_charge_from_coords(_as_float32(c), SITE) for c in SITE)

    assert total == pytest.approx(sum(v["charge"] for v in SITE.values()))
    assert total != 0.0


def test_the_float64_key_direction_also_works():
    """Assignments keyed by float32, looked up with float64."""
    assignments = {_as_float32(c): v for c, v in SITE.items()}
    writer = _writer()

    charge = writer._get_charge_from_coords((-41.874, -18.430, -51.585), assignments)

    assert charge == pytest.approx(-0.2413)


# --------------------------------------------------------------------------- #
# the exact path must be unchanged
# --------------------------------------------------------------------------- #

def test_an_exact_key_still_matches():
    writer = _writer()

    assert writer._get_charge_from_coords(
        (-41.416, -16.455, -50.088), SITE) == pytest.approx(-0.4157)


def test_an_absent_coordinate_returns_zero():
    writer = _writer()

    assert writer._get_charge_from_coords((999.0, 999.0, 999.0), SITE) == 0.0


def test_a_withheld_atom_with_charge_none_is_still_zero():
    """
    A bridging sulfide before its formal charge is collected. This must stay
    0.0 rather than crash the sum -- and must not be confused with a miss.
    """
    assignments = {(1.0, 2.0, 3.0): {"atom_name": "S1", "charge": None}}
    writer = _writer()

    assert writer._get_charge_from_coords(_as_float32((1.0, 2.0, 3.0)),
                                          assignments) == 0.0


def test_a_dataclass_style_assignment_works_through_the_fallback():
    class Assignment:
        charge = -0.5

    assignments = {(1.0, 2.0, 3.0): Assignment()}
    writer = _writer()

    assert writer._get_charge_from_coords(
        _as_float32((1.0, 2.0, 3.0)), assignments) == pytest.approx(-0.5)


# --------------------------------------------------------------------------- #
# the index
# --------------------------------------------------------------------------- #

def test_the_index_is_rebuilt_for_a_different_dict():
    """The writer is reused for the small model and then the large one."""
    writer = _writer()
    writer._get_charge_from_coords(_as_float32((-41.416, -16.455, -50.088)), SITE)

    other = {(1.0, 2.0, 3.0): {"charge": 7.5}}
    assert writer._get_charge_from_coords(
        _as_float32((1.0, 2.0, 3.0)), other) == pytest.approx(7.5)
    # and the old dict still resolves when passed again
    assert writer._get_charge_from_coords(
        _as_float32((-41.874, -18.430, -51.585)), SITE) == pytest.approx(-0.2413)


def test_non_coordinate_keys_do_not_break_the_index():
    """Assignment dicts are not guaranteed to be coordinate-keyed throughout."""
    assignments = {**SITE, "not-a-coord": {"charge": 1.0}, 7: {"charge": 2.0}}
    writer = _writer()

    assert writer._get_charge_from_coords(
        _as_float32((-41.416, -16.455, -50.088)), assignments) == pytest.approx(-0.4157)


@pytest.mark.parametrize("bad", [None, "xyz", (1.0, 2.0), (1.0, 2.0, "z")])
def test_unusable_coordinates_return_zero_rather_than_raise(bad):
    writer = _writer()

    assert writer._get_charge_from_coords(bad, SITE) == 0.0


def test_coord_key_rounds_to_pdb_precision():
    """3 decimals is what a PDB file carries, so this is lossless for them."""
    assert PDBWriter._coord_key((-41.4164999, -16.455, -50.088)) == (
        -41.416, -16.455, -50.088)
    assert PDBWriter._coord_key(_as_float32((-41.416, -16.455, -50.088))) == (
        -41.416, -16.455, -50.088)
