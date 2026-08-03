"""Heme transformer compatibility uses exact ligand criteria.

A heme's axial ligand set IS its site type: a bis-His c-type heme has two
His and two thioether Cys, a His/Met c-type has one His and one Met, a
b-type has no thioether Cys at all. A site missing one of those is a
different cofactor needing different parameters, not a near miss.

These tests pin that down, because the transformers previously accepted a
site meeting 80% of their requirements, which let each heme transformer
match its siblings' sites.
"""

import pytest

from proprep.structure_prep.comprehensive_redox_detector import (
    CenterType,
    RedoxCenter,
    RedoxSite,
    RedoxSiteAtom,
    RedoxSiteBond,
)
# Importing the package registers every transformer.
import proprep.redoxsite_prep.transformation.transformers  # noqa: F401
from proprep.redoxsite_prep.transformation.redox_transformer_framework import (
    RedoxSiteTransformerBase,
    redox_transformer_registry,
)

HEME_TRANSFORMERS = [
    "heme_bis_his_c_type",
    "heme_his_met_axial_c_type",
    "heme_bis_his_b_type",
    "heme_cys_axial_b_type",
]


def build_heme_site(heme_resname, ligands):
    """Build a heme RedoxSite.

    ``ligands`` is a list of ``(resname, chain, resid)``. Coordinates are
    synthetic but unique, which is all the validators need — they resolve
    atoms through the site's coordinate map, never through geometry.
    """
    site = RedoxSite(site_id="site_1", structure_id="test")
    counter = [0]

    def next_coord():
        counter[0] += 1
        return (float(counter[0]), 0.0, 0.0)

    coords = {}
    fe = next_coord()
    coords[(heme_resname, "FE")] = fe
    site.add_center(RedoxCenter(
        chain="A", resname=heme_resname, resid=501, coords=fe,
        center_type=CenterType.ORGANOMETALLIC_COFACTOR,
        atom_name="FE", element="FE",
    ))

    # Propionate carve-out atoms plus the two thioether-forming carbons.
    for name in ("C2A", "CAA", "C3D", "CAD", "CAB", "CAC"):
        c = next_coord()
        coords[(heme_resname, name)] = c
        site.add_atom(RedoxSiteAtom(
            chain="A", resname=heme_resname, resid=501,
            atom_name=name, coords=c, element="C",
        ))

    for resname, chain, resid in ligands:
        for name in ("CA", "CB", "NE2", "SG", "SD"):
            c = next_coord()
            coords[(resname, name, chain, resid)] = c
            site.add_atom(RedoxSiteAtom(
                chain=chain, resname=resname, resid=resid,
                atom_name=name, coords=c, element="C",
            ))

    def add_bond(c1, c2, chemical_type):
        site.bonds.append(RedoxSiteBond(
            atom1_coords=c1, atom2_coords=c2, bond_type="interresidue",
            chemical_type=chemical_type, distance=0.0,
            atom1_element="C", atom2_element="C",
        ))

    add_bond(coords[(heme_resname, "C2A")], coords[(heme_resname, "CAA")], "covalent")
    add_bond(coords[(heme_resname, "C3D")], coords[(heme_resname, "CAD")], "covalent")
    for resname, chain, resid in ligands:
        add_bond(coords[(resname, "CA", chain, resid)],
                 coords[(resname, "CB", chain, resid)], "covalent")
        axial_atom = {"HIS": "NE2", "MET": "SD", "CYS": "SG"}.get(resname)
        if axial_atom:
            add_bond(coords[(heme_resname, "FE")],
                     coords[(resname, axial_atom, chain, resid)], "coordinate")
    return site


BIS_HIS_C = [("HIS", "A", 332), ("HIS", "B", 16), ("CYS", "A", 328), ("CYS", "A", 331)]
HIS_MET_C = [("HIS", "A", 332), ("MET", "A", 60), ("CYS", "A", 328), ("CYS", "A", 331)]
BIS_HIS_B = [("HIS", "A", 60), ("HIS", "A", 90)]
CYS_AXIAL_B = [("CYS", "A", 400)]


def matching_transformers(site):
    """Names of the heme transformers that accept this site."""
    return {
        name for name in HEME_TRANSFORMERS
        if redox_transformer_registry.get_transformer(name)
        .evaluate_redox_site(site).is_valid
    }


@pytest.mark.parametrize("heme,ligands,expected", [
    ("HEC", BIS_HIS_C, "heme_bis_his_c_type"),
    ("HEM", BIS_HIS_C, "heme_bis_his_c_type"),
    ("HEC", HIS_MET_C, "heme_his_met_axial_c_type"),
    ("HEM", BIS_HIS_B, "heme_bis_his_b_type"),
    ("HEM", CYS_AXIAL_B, "heme_cys_axial_b_type"),
])
def test_exactly_one_heme_transformer_matches(heme, ligands, expected):
    """Each ligand set selects its own transformer and no other."""
    assert matching_transformers(build_heme_site(heme, ligands)) == {expected}


def test_bis_his_site_rejects_his_met_transformer():
    """The original symptom: no Met present, so no His/Met match."""
    site = build_heme_site("HEC", BIS_HIS_C)
    ev = (redox_transformer_registry
          .get_transformer("heme_his_met_axial_c_type")
          .evaluate_redox_site(site))
    assert not ev.is_valid
    failed = {d.description for d in ev.details if not d.passed}
    assert any("MET" in d for d in failed)


def test_c_type_heme_misnamed_hem_still_rejects_b_type():
    """A c-type heme deposited as HEM must not match the b-type transformer.

    HEM-vs-HEC naming was the only thing separating them; a mis-annotated
    c-type passed every b-type check. The b-type requirements now state
    that a b-type heme has no thioether Cys.
    """
    site = build_heme_site("HEM", BIS_HIS_C)
    ev = (redox_transformer_registry
          .get_transformer("heme_bis_his_b_type")
          .evaluate_redox_site(site))
    assert not ev.is_valid
    failed = {d.description for d in ev.details if not d.passed}
    assert any("CYS" in d for d in failed)


def test_confidence_still_reports_partial_credit():
    """is_valid is exact, but confidence keeps the ratio for ranking."""
    site = build_heme_site("HEC", BIS_HIS_C)
    ev = (redox_transformer_registry
          .get_transformer("heme_his_met_axial_c_type")
          .evaluate_redox_site(site))
    assert not ev.is_valid
    assert 0.0 < ev.confidence < 1.0


def test_matching_site_reports_full_confidence():
    site = build_heme_site("HEC", BIS_HIS_C)
    ev = (redox_transformer_registry
          .get_transformer("heme_bis_his_c_type")
          .evaluate_redox_site(site))
    assert ev.is_valid
    assert ev.confidence == pytest.approx(1.0)
    # The transformer declares a HEM-variant and a HEC-variant bond group
    # under require_one_group, so on a HEC site the HEM group is expected to
    # come up short; one satisfied group is the requirement.
    group_totals = [d for d in ev.details if "(total)" in d.description]
    assert group_totals and any(d.passed for d in group_totals)


# ── Bond credit is clamped ───────────────────────────────────────────────
# Callers add bond counts into the same met/total tally they use for
# pass/fail requirement checks, so surplus bonds must not buy credit that
# offsets a failed composition check.

class _BondOnlyTransformer(RedoxSiteTransformerBase):
    """Bare transformer used to exercise the bond validator directly."""


def _bonds_section(min_count, pairs, require_one_group=False, extra_group=None):
    groups = [{
        "description": "group A",
        "min_count": min_count,
        "bond_types": {"coordinate": pairs},
    }]
    if extra_group:
        groups.append(extra_group)
    return {"required_bond_groups": groups, "require_one_group": require_one_group}


def test_surplus_bonds_do_not_exceed_the_requirement():
    """Two matching bonds against a requirement of one credits one, not two."""
    site = build_heme_site("HEM", [("CYS", "A", 400), ("CYS", "A", 410)])
    found, required, _ = _BondOnlyTransformer._validate_bond_groups(
        site, _bonds_section(1, [(("HEM", "FE"), ("CYS", "SG"))]),
    )
    assert (found, required) == (1, 1)


def test_shortfall_is_reported_honestly():
    site = build_heme_site("HEM", [("CYS", "A", 400)])
    found, required, _ = _BondOnlyTransformer._validate_bond_groups(
        site, _bonds_section(3, [(("HEM", "FE"), ("CYS", "SG"))]),
    )
    assert (found, required) == (1, 3)


def test_require_one_group_prefers_a_group_that_passed():
    """A failing group with more bonds must not displace a passing one.

    Selecting the group with the largest raw count could return a group
    that fell short of a larger min_count, sinking an evaluation that an
    already-satisfied group had met.
    """
    site = build_heme_site("HEM", [("CYS", "A", 400), ("CYS", "A", 410)])
    big_failing_group = {
        "description": "group B (unreachable min_count)",
        "min_count": 99,
        "bond_types": {"covalent": [(("CYS", "CA"), ("CYS", "CB"))]},
    }
    found, required, _ = _BondOnlyTransformer._validate_bond_groups(
        site,
        _bonds_section(
            2, [(("HEM", "FE"), ("CYS", "SG"))],
            require_one_group=True, extra_group=big_failing_group,
        ),
    )
    assert (found, required) == (2, 2)


def test_no_group_satisfied_returns_zero_credit():
    site = build_heme_site("HEM", [("CYS", "A", 400)])
    found, required, _ = _BondOnlyTransformer._validate_bond_groups(
        site, _bonds_section(5, [(("HEM", "FE"), ("CYS", "SG"))],
                             require_one_group=True),
    )
    assert found == 0
    assert required == 5
