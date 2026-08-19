"""
Checklist state must round-trip numbers as numbers.

Resuming a metal-site run failed while reinserting the metals:

    Error when writing atom ((' ', 1311, ' '), ('FE1', ' ')):
    must be real number, not dict

``_serialize_value`` tested ``isinstance(value, (str, int, float, bool))`` and
fell through to a ``str()`` fallback for anything else. Coordinates read from
BioPython are numpy scalars: ``np.float32`` is NOT a Python float and has no
``__dict__``, so each coordinate was stored as

    {"__type__": "str", "value": "-46.078"}

and came back as a dict containing a string. PDBIO then refused to write it.

numpy registers its scalar types with the ``numbers`` ABCs, so they are
recognised there. ``MetalInfo.from_dict`` also unwraps and coerces, so a state
file written before this still resumes.
"""

import numpy as np
import pytest

from proprep.forcefield_prep.structure_preprocessor import MetalInfo
from proprep.utils.workflow_checklist import _deserialize_value, _serialize_value


# --------------------------------------------------------------------------- #
# serialization
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("value,expected", [
    (np.float32(-46.078), -46.078),
    (np.float64(1.5), 1.5),
    (np.int32(7), 7),
    (np.int64(-3), -3),
])
def test_numpy_scalars_serialize_as_numbers(value, expected):
    got = _serialize_value(value)

    assert not isinstance(got, dict), f"fell through to the str fallback: {got!r}"
    assert got == pytest.approx(expected)


def test_python_numbers_are_unchanged():
    assert _serialize_value(1.5) == 1.5
    assert _serialize_value(3) == 3
    assert _serialize_value("x") == "x"


def test_bool_stays_a_bool_not_an_int():
    """bool is numbers.Integral, so order of checks matters."""
    got = _serialize_value(True)

    assert got is True and isinstance(got, bool)


def test_numpy_coordinates_survive_a_round_trip():
    coords = tuple(np.array([-46.078, -17.593, -47.380], dtype=np.float32))

    restored = _deserialize_value(_serialize_value(coords))

    assert all(isinstance(c, float) for c in restored)
    assert restored == pytest.approx((-46.078, -17.593, -47.380), abs=1e-3)


# --------------------------------------------------------------------------- #
# reading state files written before the fix
# --------------------------------------------------------------------------- #

LEGACY_ATOM = {
    "atom_name": "FE1",
    "element": "Fe",
    "coords": [
        {"__type__": "str", "value": "-46.078"},
        {"__type__": "str", "value": "-17.593"},
        {"__type__": "str", "value": "-47.38"},
    ],
    "original_chain": "A",
    "original_resid": 1310,
    "original_resname": "FES",
    "is_isolated": True,
    "cluster_id": "A:1310:FES",
}


def test_legacy_string_wrapped_coordinates_are_recovered():
    """The exact shape in the reported state file."""
    metal = MetalInfo.from_dict(LEGACY_ATOM)

    assert metal.atom_name == "FE1"
    assert all(isinstance(c, float) for c in metal.coords)
    assert metal.coords == pytest.approx((-46.078, -17.593, -47.38))


def test_current_shape_still_reads():
    current = {**LEGACY_ATOM, "coords": [-46.078, -17.593, -47.38]}

    assert MetalInfo.from_dict(current).coords == pytest.approx(
        (-46.078, -17.593, -47.38))


def test_a_tuple_wrapped_coordinate_list_is_recovered():
    wrapped = {**LEGACY_ATOM,
               "coords": {"__type__": "tuple", "value": [-1.0, -2.0, -3.0]}}

    assert MetalInfo.from_dict(wrapped).coords == pytest.approx((-1.0, -2.0, -3.0))


def test_an_unparseable_coordinate_does_not_abort_the_resume():
    broken = {**LEGACY_ATOM, "coords": ["x", -17.593, -47.38]}

    coords = MetalInfo.from_dict(broken).coords

    assert coords[0] == 0.0
    assert coords[1] == pytest.approx(-17.593)


def test_coordinates_are_writable_by_pdbio():
    """The failure was a TypeError inside PDBIO's format string."""
    metal = MetalInfo.from_dict(LEGACY_ATOM)

    line = "%8.3f%8.3f%8.3f" % metal.coords   # must not raise

    assert "-46.078" in line
