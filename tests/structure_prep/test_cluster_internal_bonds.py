"""
A withheld cluster's own bonds have to reach the Seminario step.

The bonded-parameter step draws bonds from two sources, and a pure inorganic
cluster is in neither:

- LINK records, which are ``chemical_type == 'coordinate'`` bonds on the site,
  i.e. metal-to-ligand ACROSS residues
- the prmtop, which never saw the residue -- it is withheld whole

So Fe-S inside FES, and Mo-S/Mo-O inside MOS plus the O-H of the hydroxo, never
got a force constant. tleap builds the residue from the deposited library's
connectivity table, finds those bonds and reports a missing parameter for each
(``M3-YB``, ``YA-H``).

Perceived from geometry rather than asked for, and tagged ``cluster_internal``
rather than ``coordinate``. Three consumers key on that exact string and must
not see these:

- the tleap generator emits ``bond`` commands for coordinate bonds; an
  intra-residue bond is redundant with the library unit
- the auto-transformer's WL fingerprint walks coordinate bonds; an
  intra-residue bond is a self-edge there
- the LINK writer, which is widened deliberately

Selecting on ``redox_site.centers`` does NOT work: for a cluster the center
describes the RESIDUE, with element None and centroid coords matching no atom.
"""

import pytest

from proprep.structure_prep.comprehensive_redox_detector import (
    CLUSTER_INTERNAL, RedoxSite, RedoxSiteAtom, perceive_cluster_internal_bonds,
)


def _atom(resname, resid, name, element, coords):
    return RedoxSiteAtom(chain="A", resname=resname, resid=resid,
                         atom_name=name, coords=coords, element=element)


# The reported clusters, at their real geometries.
FES_ATOMS = [
    _atom("FES", 1311, "FE1", "FE", (-46.078, -17.593, -47.380)),
    _atom("FES", 1311, "FE2", "FE", (-43.100, -18.216, -47.517)),
    _atom("FES", 1311, "S1", "S", (-44.425, -16.836, -48.664)),
    _atom("FES", 1311, "S2", "S", (-44.818, -19.213, -46.534)),
]

MOS_ATOMS = [
    _atom("MOS", 1312, "MO", "MO", (-29.569, -17.017, -41.812)),
    _atom("MOS", 1312, "S", "S", (-29.451, -17.649, -39.574)),
    _atom("MOS", 1312, "O1", "O", (-27.666, -16.557, -41.673)),
    _atom("MOS", 1312, "O2", "O", (-29.607, -18.617, -42.453)),
    _atom("MOS", 1312, "H1", "H", (-27.532, -15.655, -41.974)),
]


def _site(atoms, bonds=None):
    site = RedoxSite(site_id="site_1", structure_id="4UHX")
    site.atoms = list(atoms)
    site.bonds = list(bonds or [])
    site.centers = []
    site.coord_to_pdb = {
        a.coords: {"chain": a.chain, "resid": a.resid, "resname": a.resname,
                   "atom_name": a.atom_name, "element": a.element}
        for a in atoms
    }
    return site


def _pairs(site):
    return {
        tuple(sorted((b.atom1_residue_info["atom_name"],
                      b.atom2_residue_info["atom_name"])))
        for b in site.bonds if b.chemical_type == CLUSTER_INTERNAL
    }


# --------------------------------------------------------------------------- #
# the reported clusters
# --------------------------------------------------------------------------- #

def test_the_molybdenum_cofactor_gets_its_four_bonds():
    site = _site(MOS_ATOMS)

    assert perceive_cluster_internal_bonds(site) == 4
    assert _pairs(site) == {("MO", "S"), ("MO", "O1"), ("MO", "O2"), ("H1", "O1")}


def test_the_hydroxo_o_h_is_included():
    """Not a metal bond, but it is inside the withheld residue and unparameterized."""
    site = _site(MOS_ATOMS)
    perceive_cluster_internal_bonds(site)

    assert ("H1", "O1") in _pairs(site)


def test_the_iron_sulfur_cluster_gets_its_four_fe_s_bonds():
    site = _site(FES_ATOMS)

    assert perceive_cluster_internal_bonds(site) == 4
    assert _pairs(site) == {
        ("FE1", "S1"), ("FE1", "S2"), ("FE2", "S1"), ("FE2", "S2")}


def test_no_metal_metal_bond_is_added():
    """Fe...Fe clears any radius cutoff; MCPB bridges them through the sulfides."""
    site = _site(FES_ATOMS)
    perceive_cluster_internal_bonds(site)

    assert ("FE1", "FE2") not in _pairs(site)


# --------------------------------------------------------------------------- #
# the tag, which three consumers filter on
# --------------------------------------------------------------------------- #

def test_the_bonds_are_not_tagged_coordinate():
    site = _site(MOS_ATOMS)
    perceive_cluster_internal_bonds(site)

    assert all(b.chemical_type == CLUSTER_INTERNAL for b in site.bonds)
    assert not any(b.chemical_type == "coordinate" for b in site.bonds)


def test_the_bonds_are_marked_intraresidue():
    site = _site(MOS_ATOMS)
    perceive_cluster_internal_bonds(site)

    assert all(b.bond_type == "intraresidue" for b in site.bonds)


def test_residue_info_is_populated_for_the_link_writer():
    """The LINK writer reads atom_name off both endpoints."""
    site = _site(MOS_ATOMS)
    perceive_cluster_internal_bonds(site)

    for bond in site.bonds:
        assert bond.atom1_residue_info.get("atom_name")
        assert bond.atom2_residue_info.get("atom_name")


# --------------------------------------------------------------------------- #
# what must be left alone
# --------------------------------------------------------------------------- #

def test_an_organometallic_cofactor_is_skipped():
    """
    A heme's internal bonds come from its own library. Perceiving forty of them
    would hand Seminario bonds it has no business parameterizing. Carbon in the
    residue is what distinguishes it from a pure inorganic cluster.
    """
    hem = [
        _atom("HEM", 1, "FE", "FE", (0.0, 0.0, 0.0)),
        _atom("HEM", 1, "NA", "N", (2.0, 0.0, 0.0)),
        _atom("HEM", 1, "C1A", "C", (3.2, 0.6, 0.0)),
    ]
    site = _site(hem)

    assert perceive_cluster_internal_bonds(site) == 0


def test_a_residue_with_no_metal_is_skipped():
    cys = [
        _atom("CYM", 114, "CB", "C", (0.0, 0.0, 0.0)),
        _atom("CYM", 114, "SG", "S", (1.8, 0.0, 0.0)),
    ]

    assert perceive_cluster_internal_bonds(_site(cys)) == 0


def test_a_lone_metal_ion_has_nothing_to_add():
    site = _site([_atom("ZN", 261, "ZN", "ZN", (0.0, 0.0, 0.0))])

    assert perceive_cluster_internal_bonds(site) == 0


def test_an_existing_bond_is_not_duplicated():
    site = _site(MOS_ATOMS)
    perceive_cluster_internal_bonds(site)
    before = len(site.bonds)

    assert perceive_cluster_internal_bonds(site) == 0
    assert len(site.bonds) == before


def test_a_coordination_bond_already_present_is_preserved():
    """The Mo-S(MTE) links must survive untouched."""
    from proprep.structure_prep.comprehensive_redox_detector import RedoxSiteBond

    existing = RedoxSiteBond(
        atom1_coords=MOS_ATOMS[0].coords, atom2_coords=(-31.0, -18.0, -40.0),
        bond_type="interresidue", chemical_type="coordinate", distance=2.4,
        atom1_element="MO", atom2_element="S")
    site = _site(MOS_ATOMS, bonds=[existing])

    perceive_cluster_internal_bonds(site)

    assert existing in site.bonds
    assert sum(1 for b in site.bonds if b.chemical_type == "coordinate") == 1


def test_a_site_with_no_atoms_is_not_an_error():
    assert perceive_cluster_internal_bonds(_site([])) == 0


# --------------------------------------------------------------------------- #
# the LINK writer must carry them
# --------------------------------------------------------------------------- #

def test_the_link_writer_accepts_the_new_type(tmp_path):
    from proprep.forcefield_prep.mcpb.fingerprint_generator import FingerprintGenerator

    site = _site(MOS_ATOMS)
    perceive_cluster_internal_bonds(site)
    coord_to_id = {a.coords: i + 1 for i, a in enumerate(MOS_ATOMS)}

    out = tmp_path / "standard.fingerprint"
    with open(out, "w") as fh:
        FingerprintGenerator()._write_link_records(fh, site, coord_to_id)

    text = out.read_text()
    assert text.count("LINK") == 4
    assert "MO" in text and "O1" in text
