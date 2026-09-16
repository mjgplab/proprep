"""Site templates resolve bonds and boundary atoms by identity, not position.

Row numbers in the bond table and atom indices in the site list are both
site-specific: same-type residues are ordered by a near-tie distance sort,
and transformed HCO hemes carry grafted atoms in per-site order. Replaying
either as a position produced 6-12 Å "bonds" and the wrong boundary on
the maxtilt structure (64 hemes, 48 with a different HCO atom order)."""

import os

from rich.console import Console

from proprep.structure_prep.comprehensive_redox_detector import (
    BoundaryDefinition,
    DetectionConfig,
    RedoxSite,
    RedoxSiteAtom,
    SiteRefinementInterface,
    SiteTemplate,
    _boundary_atom_identities,
    _resolve_boundary_atoms,
)


def _atom(chain, resname, resid, name, xyz, element="C"):
    return RedoxSiteAtom(chain, resname, resid, name, tuple(float(v) for v in xyz), element)


def _refiner():
    return SiteRefinementInterface(
        DetectionConfig(), console=Console(file=open(os.devnull, "w")), template_mode=False
    )


def _template(bond_selections, boundary_atoms=None, boundary_indices=None):
    return SiteTemplate(
        site_type="heme_bis_his_c_type",
        residue_types=["CYO"], residue_counts=[2],
        atom_filtering_metals=[], atom_filtering_nonmetals=[],
        boundary_definition=BoundaryDefinition.MIN_DISTANCE_CUSTOM,
        residue_selection_choice="1-2",
        bond_pairs=list(bond_selections),
        bond_selections=bond_selections,
        continue_searching=False,
        custom_boundary_atom_indices=boundary_indices,
        custom_boundary_atoms=boundary_atoms,
    )


def _heme_with_two_cys(cys_order):
    """HCO with thioether carbons CBB2 at x=0 and CBC1 at x=6; two CYO whose
    CA sit 1.5 Å from one each. ``cys_order`` lists which CYO is added first."""
    site = RedoxSite("s", "t")
    site.atoms += [
        _atom("A", "HCO", 300, "FE", (3, 3, 0), "FE"),
        _atom("A", "HCO", 300, "CBB2", (0, 0, 0)),
        _atom("A", "HCO", 300, "CBC1", (6, 0, 0)),
    ]
    cys = {
        "near_CBB2": _atom("A", "CYO", 25, "CA", (0, 1.5, 0)),
        "near_CBC1": _atom("A", "CYO", 28, "CA", (6, 1.5, 0)),
    }
    for key in cys_order:
        site.atoms.append(cys[key])
    return site


def _bond_map(site):
    return {
        (b.atom1_residue_info["resid"], b.atom1_residue_info["atom_name"],
         b.atom2_residue_info["atom_name"]): round(b.distance, 2)
        for b in site.bonds
    }


# ------------------------------------------------------------ bonds ----


def test_same_type_residues_matched_by_geometry_not_row_order():
    # Captured on a site where row 2 was the CBB2 Cys and row 3 the CBC1 Cys.
    bonds = {
        "2-1": [{"source_atom": "CA", "target_atom": "CBB2",
                 "source_resname": "CYO", "target_resname": "HCO"}],
        "3-1": [{"source_atom": "CA", "target_atom": "CBC1",
                 "source_resname": "CYO", "target_resname": "HCO"}],
    }
    for order in (["near_CBB2", "near_CBC1"], ["near_CBC1", "near_CBB2"]):
        site = _heme_with_two_cys(order)
        _refiner()._apply_template_bonds(site, _template(bonds))
        assert _bond_map(site) == {(25, "CA", "CBB2"): 1.5, (28, "CA", "CBC1"): 1.5}, order


def test_legacy_bond_without_resnames_falls_back_to_row_order():
    bonds = {"2-1": [{"source_atom": "CA", "target_atom": "CBB2"}]}
    site = _heme_with_two_cys(["near_CBC1", "near_CBB2"])  # row 2 is the far Cys
    _refiner()._apply_template_bonds(site, _template(bonds))
    assert _bond_map(site) == {(28, "CA", "CBB2"): 6.18}


def test_identical_bonds_land_on_distinct_residues():
    """Two 'HIS NE2 - FE' bonds must not both pick the closer His."""
    site = RedoxSite("s", "t")
    site.atoms += [
        _atom("A", "HEM", 300, "FE", (0, 0, 0), "FE"),
        _atom("A", "HIS", 18, "NE2", (0, 0, 2.0), "N"),
        _atom("A", "HIS", 87, "NE2", (0, 0, -2.1), "N"),
    ]
    bonds = {
        "2-1": [{"source_atom": "NE2", "target_atom": "FE",
                 "source_resname": "HIS", "target_resname": "HEM"}],
        "3-1": [{"source_atom": "NE2", "target_atom": "FE",
                 "source_resname": "HIS", "target_resname": "HEM"}],
    }
    _refiner()._apply_template_bonds(site, _template(bonds))
    assert sorted(b.atom1_residue_info["resid"] for b in site.bonds) == [18, 87]


def test_missing_atom_is_skipped_not_misbonded():
    site = _heme_with_two_cys(["near_CBB2", "near_CBC1"])
    bonds = {"2-1": [{"source_atom": "SG", "target_atom": "CBB2",
                      "source_resname": "CYO", "target_resname": "HCO"}]}
    _refiner()._apply_template_bonds(site, _template(bonds))
    assert site.bonds == []


# --------------------------------------------------------- boundary ----


def _heme_atoms(order):
    coords = {"FE": (0, 0, 0), "CBB2": (1, 0, 0), "CB2": (2, 0, 0), "NB": (3, 0, 0), "C3D": (4, 0, 0)}
    return [_atom("A", "HCO", 300, n, coords[n]) for n in order]


def test_boundary_atoms_resolved_by_name_across_atom_orders():
    template_site = RedoxSite("s1", "t")
    template_site.atoms = _heme_atoms(["FE", "CBB2", "CB2", "NB", "C3D"])
    indices = [1, 2]  # CBB2, CB2 as the user picked them
    identities = _boundary_atom_identities(template_site, indices)
    assert identities == [("HCO", 0, "CBB2"), ("HCO", 0, "CB2")]

    other = RedoxSite("s2", "t")
    other.atoms = _heme_atoms(["FE", "CBB2", "NB", "C3D", "CB2"])  # index 2 is now NB
    coords, warnings = _resolve_boundary_atoms(other, identities, indices)
    assert coords == [(1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]
    assert warnings == []


def test_boundary_falls_back_to_index_when_name_absent():
    other = RedoxSite("s2", "t")
    other.atoms = _heme_atoms(["FE", "CBB2", "NB"])
    coords, warnings = _resolve_boundary_atoms(other, [("HCO", 0, "CB2")], [2])
    assert coords == [(3.0, 0.0, 0.0)]
    assert any("not found by name" in w for w in warnings)


def test_boundary_index_only_template_still_works():
    other = RedoxSite("s2", "t")
    other.atoms = _heme_atoms(["FE", "CBB2", "NB"])
    coords, warnings = _resolve_boundary_atoms(other, None, [1, 9])
    assert coords == [(1.0, 0.0, 0.0)]
    assert any("out of range" in w for w in warnings)


def test_ordinal_separates_two_residues_of_one_name():
    site = RedoxSite("s", "t")
    site.atoms = [_atom("A", "HEM", 1, "FE", (0, 0, 0), "FE"), _atom("A", "HEM", 2, "FE", (5, 0, 0), "FE")]
    assert _boundary_atom_identities(site, [1]) == [("HEM", 1, "FE")]
    coords, _ = _resolve_boundary_atoms(site, [("HEM", 1, "FE")], None)
    assert coords == [(5.0, 0.0, 0.0)]
