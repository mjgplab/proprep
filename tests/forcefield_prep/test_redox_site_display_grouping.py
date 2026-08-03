"""
Regression tests for _group_residues_into_display_units: what ends up in ONE
parameterization unit.

Two behaviors are pinned here.

1. A multi-member redox site collapses to ONE unit, not one per center.
   7K0W has a binuclear Mn center bridged by the baloxavir inhibitor E4Z (one
   RedoxSite: MN 202 + MN 203 + E4Z), plus two independent SO4 ions. The residue
   list used to flatten that to three rows reading as three independent jobs.

2. A covalent AA-ligand adduct arrives COMPLETE.
   Extraction seeds a unit only from a site's redox CENTERS and skips centers
   that are standard amino acids, so a covalent inhibitor (KRAS 6OIM: CYS A:12
   SG bonded to MOV A:303 C25) reached classification as the ligand alone and
   read as a lone small molecule. Grouping pulls the covalently bonded partner
   back in, so the unit classifies as a from-structure modified amino acid.
   A metal site is NOT expanded that way: MCPB re-detects its own coordination
   sphere, and a c-type heme's Cys-heme thioether links are covalent too.

These tests are hermetic: every residue carries the CCD/element data the
classifier needs, so no network lookup happens (the ccd_parser stub raises).
"""

import pytest
from Bio.PDB.StructureBuilder import StructureBuilder

from proprep.forcefield_prep.forcefield_parameterizer import ForcefieldParameterizer
from proprep.forcefield_prep.forcefield_worker import NonStandardResidue
from proprep.structure_prep.comprehensive_redox_detector import (
    RedoxSiteAtom,
    RedoxSiteBond,
)


class _FakeSite:
    """Stands in for a RedoxSite: grouping reads site_id, atoms and bonds."""

    def __init__(self, site_id, atoms=None, bonds=None):
        self.site_id = site_id
        self.atoms = atoms or []
        self.bonds = bonds or []


class _ExplodingCCD:
    def get_residue_data(self, name):  # pragma: no cover - must never be hit
        raise AssertionError(f"unexpected CCD network lookup for {name}")


def _res(name, chain, resid, site=None, *, ccd_data=None, atom_count=0,
         elements=None):
    r = NonStandardResidue(name=name, chain_id=chain, resid=resid,
                           category="unknown")
    if site is not None:
        r.source_redox_site = site
        r.redox_site_id = site.site_id
    if ccd_data is not None:
        r.ccd_data = ccd_data
    if atom_count:
        r.atom_count = atom_count
    if elements:
        r.redox_site_elements = set(elements)
    return r


def _bare_parameterizer(residues, structure=None):
    # Bypass the heavy __init__ (CCD parser, FF tables); pin only what the
    # grouping + classification path reads.
    p = object.__new__(ForcefieldParameterizer)
    p.non_standard_residues = residues
    p.user_residue_classifications = {}
    p.classification_settings = {
        "min_small_molecule_atoms": 2,
        "max_small_molecule_atoms": 100,
    }
    p.ccd_parser = _ExplodingCCD()
    p._get_structure_object = lambda: structure
    return p


def _atom(chain, resname, resid, name, element, coords):
    return RedoxSiteAtom(chain=chain, resname=resname, resid=resid,
                         atom_name=name, coords=coords, element=element)


def _bond(a1, a2, chemical_type):
    return RedoxSiteBond(
        atom1_coords=a1.coords, atom2_coords=a2.coords,
        bond_type="interresidue", chemical_type=chemical_type, distance=1.8,
        atom1_element=a1.element, atom2_element=a2.element,
        atom1_residue_info={"chain": a1.chain, "resname": a1.resname,
                            "resid": a1.resid, "atom_name": a1.atom_name},
        atom2_residue_info={"chain": a2.chain, "resname": a2.resname,
                            "resid": a2.resid, "atom_name": a2.atom_name},
    )


def _structure(residues):
    """Build a BioPython structure: residues = [(resname, resid, [atom names])]."""
    sb = StructureBuilder()
    sb.init_structure("t")
    sb.init_model(0)
    sb.init_chain("A")
    sb.init_seg(" ")
    serial = 1
    for resname, resid, atom_names in residues:
        sb.init_residue(resname, " ", resid, " ")
        for i, name in enumerate(atom_names):
            sb.init_atom(name, (float(i), float(resid), 0.0), 0.0, 1.0, " ",
                         name, serial, name[0])
            serial += 1
    return sb.get_structure()


# ── 1. site collapses to one unit ──────────────────────────────────────────

def test_binuclear_site_collapses_so4_stay_loose():
    site = _FakeSite("Mn-binuclear-A202")
    so4_ccd = {"type": "NON-POLYMER"}
    residues = [
        _res("E4Z", "A", 201, site, atom_count=30),      # bridging inhibitor
        _res("MN", "A", 202, site),
        _res("MN", "A", 203, site),
        _res("SO4", "A", 204, ccd_data=so4_ccd, atom_count=5),   # independent
        _res("SO4", "A", 205, ccd_data=so4_ccd, atom_count=5),   # independent
    ]
    units = _bare_parameterizer(residues)._group_residues_into_display_units()

    # Three site members -> ONE metal-site unit; two loose SO4 -> two units.
    assert len(units) == 3

    site_units = [u for u in units if u["site_id"] == "Mn-binuclear-A202"]
    assert len(site_units) == 1
    su = site_units[0]
    assert su["category"] == "metal_site"
    assert {m.resid for m in su["metals"]} == {202, 203}
    assert [l.name for l in su["ligands"]] == ["E4Z"]

    loose = [u["members"][0].name for u in units if u["site_id"] is None]
    assert loose == ["SO4", "SO4"]
    assert all(u["category"] == "small_molecule"
               for u in units if u["site_id"] is None)


def test_label_reads_as_one_binuclear_site():
    site = _FakeSite("Mn-binuclear-A202")
    residues = [
        _res("MN", "A", 202, site),
        _res("MN", "A", 203, site),
        _res("E4Z", "A", 201, site, atom_count=30),
    ]
    p = _bare_parameterizer(residues)
    unit = p._group_residues_into_display_units()[0]
    name, category, status = p._format_unit(unit)

    assert "Metal site" in name
    assert "MN(A:202)" in name and "MN(A:203)" in name
    assert "E4Z" in name           # ligand shown in the row
    assert category == "Metal Site"
    assert "binuclear" in status


def test_single_member_site_stays_one_residue():
    # The overwhelmingly common case: one lone metal site -> unchanged behavior.
    site = _FakeSite("Zn-mono-A300")
    residues = [_res("ZN", "A", 300, site)]
    units = _bare_parameterizer(residues)._group_residues_into_display_units()
    assert len(units) == 1
    assert [m.name for m in units[0]["members"]] == ["ZN"]
    assert units[0]["category"] == "metal_site"


def test_no_redox_site_all_loose():
    residues = [
        _res("SO4", "A", 204, ccd_data={"type": "NON-POLYMER"}, atom_count=5),
        _res("HEM", "A", 1, elements=["FE"]),
    ]
    units = _bare_parameterizer(residues)._group_residues_into_display_units()
    assert len(units) == 2
    assert all(len(u["members"]) == 1 and u["site_id"] is None for u in units)
    assert [u["category"] for u in units] == ["small_molecule", "metal_site"]


# ── 2. covalent adduct arrives complete ────────────────────────────────────

def _kras_like_site():
    """CYS A:12 SG covalently bonded to ligand MOV A:303 C25 (6OIM)."""
    sg = _atom("A", "CYS", 12, "SG", "S", (-6.344, -3.260, 0.409))
    cys_atoms = [
        _atom("A", "CYS", 12, "N", "N", (0.0, 12.0, 0.0)),
        _atom("A", "CYS", 12, "CA", "C", (1.0, 12.0, 0.0)),
        _atom("A", "CYS", 12, "C", "C", (2.0, 12.0, 0.0)),
        _atom("A", "CYS", 12, "O", "O", (3.0, 12.0, 0.0)),
        sg,
    ]
    c25 = _atom("A", "MOV", 303, "C25", "C", (-5.364, -2.168, -0.643))
    mov_atoms = [c25, _atom("A", "MOV", 303, "C21", "C", (-2.61, -7.064, -0.371))]
    site = _FakeSite("site_1", atoms=cys_atoms + mov_atoms,
                     bonds=[_bond(sg, c25, "covalent")])
    return site


def test_covalent_adduct_pulls_in_the_bonded_amino_acid():
    site = _kras_like_site()
    # Extraction yields the ligand ONLY: CYS is a standard AA and not a center.
    # MOV's own CCD entry says non-polymer, so on its own it reads as a small
    # molecule -- that was the reported bug, and is what this test pins.
    residues = [_res("MOV", "A", 303, site, ccd_data={"type": "NON-POLYMER"},
                     atom_count=41)]
    structure = _structure([
        ("ALA", 11, ["N", "CA", "C", "O", "CB"]),
        ("CYS", 12, ["N", "CA", "C", "O", "CB", "SG"]),
        ("GLY", 13, ["N", "CA", "C", "O"]),
        ("MOV", 303, ["C25", "C21"]),
    ])

    units = _bare_parameterizer(residues, structure)._group_residues_into_display_units()

    assert len(units) == 1
    unit = units[0]
    assert unit["category"] == "modified_amino_acid"
    assert unit["procedure"] == "from_structure"
    assert {(m.name, m.resid) for m in unit["members"]} == {("MOV", 303), ("CYS", 12)}

    # The pulled-in Cys is a unit member only; it is not a non-standard residue
    # and must not leak into the analyzed list.
    assert [r.name for r in residues] == ["MOV"]
    cys = next(m for m in unit["members"] if m.name == "CYS")
    assert cys.is_covalent_partner is True


def test_metal_site_is_not_expanded_by_covalent_bonds():
    # A c-type heme is covalently anchored by CYS-SG thioether links. Those must
    # not pull the Cys into the unit: MCPB re-detects the site itself.
    fe = _atom("A", "HEM", 1, "FE", "FE", (0.0, 0.0, 0.0))
    cab = _atom("A", "HEM", 1, "CAB", "C", (2.0, 0.0, 0.0))
    sg14 = _atom("A", "CYS", 14, "SG", "S", (3.5, 0.0, 0.0))
    ne2 = _atom("A", "HIS", 18, "NE2", "N", (0.0, 0.0, 2.1))
    site = _FakeSite("cheme", atoms=[fe, cab, sg14, ne2],
                     bonds=[_bond(sg14, cab, "covalent"),
                            _bond(fe, ne2, "coordinate")])
    residues = [_res("HEM", "A", 1, site, elements=["FE"])]

    units = _bare_parameterizer(residues)._group_residues_into_display_units()

    assert len(units) == 1
    assert units[0]["category"] == "metal_site"
    assert [m.name for m in units[0]["members"]] == ["HEM"]
    assert units[0]["ligands"] == []


def test_nonstandard_partner_is_absorbed_not_duplicated():
    # If the covalent partner is itself NON-standard, the structure scan picks
    # it up as its own residue while the site also bonds it to the center. It
    # must end up in the conjugate unit ONCE, not parameterized twice.
    c1 = _atom("A", "LIG", 500, "C1", "C", (0.0, 0.0, 0.0))
    n1 = _atom("A", "NAG", 600, "N1", "N", (1.5, 0.0, 0.0))
    site = _FakeSite("s", atoms=[c1, n1], bonds=[_bond(c1, n1, "covalent")])

    lig = _res("LIG", "A", 500, site, ccd_data={"type": "NON-POLYMER"},
               atom_count=20)
    nag = _res("NAG", "A", 600, ccd_data={"type": "NON-POLYMER"}, atom_count=14)
    units = _bare_parameterizer([lig, nag])._group_residues_into_display_units()

    assert len(units) == 1
    assert {m.name for m in units[0]["members"]} == {"LIG", "NAG"}
    # The scanned object itself is reused, so its CCD data is not thrown away.
    assert any(m is nag for m in units[0]["members"])


def test_coordinate_bond_alone_does_not_expand_a_metal_free_unit():
    # Only covalent bonds promote a partner into the unit.
    lig = _atom("A", "LIG", 300, "O1", "O", (0.0, 0.0, 0.0))
    od1 = _atom("A", "ASP", 50, "OD1", "O", (2.5, 0.0, 0.0))
    site = _FakeSite("s", atoms=[lig, od1],
                     bonds=[_bond(lig, od1, "coordinate")])
    residues = [_res("LIG", "A", 300, site, ccd_data={"type": "NON-POLYMER"},
                     atom_count=10)]

    units = _bare_parameterizer(residues)._group_residues_into_display_units()

    assert len(units) == 1
    assert [m.name for m in units[0]["members"]] == ["LIG"]
    assert units[0]["category"] == "small_molecule"


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-v"]))
