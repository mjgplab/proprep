"""
A RedoxSite must rebuild from checklist state, not only from exported JSON.

Resuming a metal-site run failed at structure recombination with

    tLEaP error: 'site_id'

and the next step then died on the prepared structure that step never produced.
``dict_to_redox_site`` is the documented normalizer for "a resumed session can
hand object-expecting code dict-form sites", but it only understood the flat
dict the exporters write. Checklist state stores objects as
``{"__type__": ..., "value": {...}}`` and spells three fields differently,
because it dumps the dataclass field names:

    exported                      checklist state
    coordinates                   coords
    atom1_coordinates             atom1_coords
    atom1                         atom1_residue_info

Both shapes are now accepted.
"""

import pytest

from proprep.structure_prep.comprehensive_redox_detector import (
    RedoxSite, dict_to_redox_site,
)
from proprep.utils.workspace import unwrap_serialized


EXPORTED = {
    "site_id": "site_1",
    "structure_id": "4UHX",
    "site_type": "metal",
    "centers": [{
        "chain": "A", "resname": "FES", "resid": 1310,
        "atom_name": "FE1", "coordinates": [1.0, 2.0, 3.0],
        "center_type": "organometallic_cofactor", "element": "FE",
    }],
    "atoms": [{
        "chain": "A", "resname": "FES", "resid": 1310,
        "atom_name": "FE1", "coordinates": [1.0, 2.0, 3.0], "element": "FE",
    }],
    "bonds": [{
        "atom1_coordinates": [1.0, 2.0, 3.0], "atom2_coordinates": [4.0, 5.0, 6.0],
        "bond_type": "coordinate", "chemical_type": "coordinate", "distance": 2.3,
        "atom1": {"resid": 1310}, "atom2": {"resid": 114},
    }],
}

# The same site as checklist state writes it: wrapped, dataclass field names.
STATE = {
    "__type__": "RedoxSite",
    "value": {
        "site_id": "site_1",
        "structure_id": "4UHX",
        "site_type": "metal",
        "centers": [{
            "chain": "A", "resname": "FES", "resid": 1310,
            "atom_name": "FE1",
            "coords": {"__type__": "tuple", "value": [1.0, 2.0, 3.0]},
            "center_type": {"__type__": "CenterType",
                            "value": {"_value_": "organometallic_cofactor"}},
            "element": "FE",
        }],
        "atoms": [{
            "chain": "A", "resname": "FES", "resid": 1310,
            "atom_name": "FE1",
            "coords": {"__type__": "tuple", "value": [1.0, 2.0, 3.0]},
            "element": "FE",
        }],
        "bonds": [{
            "atom1_coords": {"__type__": "tuple", "value": [1.0, 2.0, 3.0]},
            "atom2_coords": {"__type__": "tuple", "value": [4.0, 5.0, 6.0]},
            "bond_type": "coordinate", "chemical_type": "coordinate",
            "distance": 2.3,
            "atom1_residue_info": {"resid": 1310},
            "atom2_residue_info": {"resid": 114},
        }],
    },
}


def test_rebuilds_from_checklist_state():
    """The reported failure: KeyError 'site_id' on the wrapper."""
    site = dict_to_redox_site(STATE)

    assert isinstance(site, RedoxSite)
    assert site.site_id == "site_1"
    assert len(site.centers) == 1
    assert len(site.atoms) == 1
    assert len(site.bonds) == 1


def test_the_two_shapes_agree():
    """Same site, two serializations, same object."""
    a = dict_to_redox_site(EXPORTED)
    b = dict_to_redox_site(STATE)

    assert a.site_id == b.site_id
    assert a.centers[0].coords == b.centers[0].coords == (1.0, 2.0, 3.0)
    assert a.centers[0].center_type == b.centers[0].center_type
    assert a.atoms[0].coords == b.atoms[0].coords
    assert a.bonds[0].atom1_coords == b.bonds[0].atom1_coords
    assert a.bonds[0].atom2_coords == b.bonds[0].atom2_coords
    assert a.bonds[0].atom1_residue_info == b.bonds[0].atom1_residue_info


def test_still_reads_the_exported_shape():
    """The shape it always handled must keep working."""
    site = dict_to_redox_site(EXPORTED)

    assert site.centers[0].coords == (1.0, 2.0, 3.0)
    assert site.bonds[0].distance == pytest.approx(2.3)


def test_is_idempotent_on_objects():
    """It doubles as a normalizer, so objects pass through untouched."""
    site = dict_to_redox_site(STATE)

    assert dict_to_redox_site(site) is site


def test_a_site_saved_before_bonds_exist_rebuilds_empty():
    """State snapshotted at step 8 has no bonds yet; that is not an error."""
    state = {"__type__": "RedoxSite",
             "value": {**STATE["value"], "bonds": []}}

    site = dict_to_redox_site(state)

    assert site.bonds == []
    assert len(site.atoms) == 1


# --------------------------------------------------------------------------- #
# the shared unwrapper
# --------------------------------------------------------------------------- #

def test_unwrap_handles_the_nesting_checklist_state_produces():
    assert unwrap_serialized({"__type__": "tuple", "value": [1, 2, 3]}) == (1, 2, 3)
    assert unwrap_serialized(
        {"__type__": "CenterType", "value": {"_value_": "metal_ion"}}) == "metal_ion"
    assert unwrap_serialized({"__type__": "Path", "value": "/tmp/x"}) == "/tmp/x"
    assert unwrap_serialized(
        {"a": {"__type__": "tuple", "value": [1]}}) == {"a": (1,)}
    assert unwrap_serialized(
        [{"__type__": "tuple", "value": [1]}]) == [(1,)]


def test_unwrap_leaves_ordinary_values_alone():
    plain = {"site_id": "s", "n": 3, "xs": [1, 2], "nested": {"k": "v"}}

    assert unwrap_serialized(plain) == plain
