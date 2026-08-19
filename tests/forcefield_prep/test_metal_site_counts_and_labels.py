"""
Metal-site counts and labels must describe sites and atoms, not residues.

Four defects in one screen of 4UHX output, all the same species — a count taken
over residues where the meaningful unit is a site or an atom:

  "MCPB (mononuclear)" for an Fe2S2 cluster
      nuclearity counted metal-bearing RESIDUES. FES is one residue holding two
      Fe, so a binuclear cluster read as mononuclear.

  "Found 4 metal sites" for three
      every member of a metal-site unit is stamped category="metal_site", the
      coordinating ligand included, so MoCo's MOS + MTE counted twice.

  "Multi-Site Metal Parameterization: 4 sites" after selecting three
      len(combined_residues) rather than the number of selected units.

  "Parameterize with option 2, view status with option 3"
      the menu is 1 analyze / 2 import / 3 parameterize / 4 status / 5 help.
      The suggestions kept the numbering from before "import" was added.
"""

import pytest

from proprep.forcefield_prep.forcefield_parameterizer import ForcefieldParameterizer


class _Atom:
    def __init__(self, chain, resid, resname, element):
        self.chain = chain
        self.resid = resid
        self.resname = resname
        self.element = element


class _Site:
    """Stand-in RedoxSite: the parameterizer only reads .atoms and .site_id."""

    def __init__(self, site_id, atoms):
        self.site_id = site_id
        self.atoms = atoms


class _Res:
    def __init__(self, name, chain_id, resid, site=None, category=None):
        self.name = name
        self.chain_id = chain_id
        self.resid = resid
        self.category = category
        self.source_redox_site = site
        self.redox_site_id = site.site_id if site else None
        self.biopython_residue = None
        self.ccd_data = None


def _param():
    return ForcefieldParameterizer.__new__(ForcefieldParameterizer)


def _fes_unit():
    """4UHX site_1: one FES residue holding FE1 + FE2 and two bridging S."""
    site = _Site("site_1", [
        _Atom("A", 3001, "FES", "FE"),
        _Atom("A", 3001, "FES", "FE"),
        _Atom("A", 3001, "FES", "S"),
        _Atom("A", 3001, "FES", "S"),
    ])
    fes = _Res("FES", "A", 3001, site, "metal_site")
    return {"members": [fes], "metals": [fes], "ligands": [],
            "category": "metal_site", "site_id": "site_1"}


def _mos_unit():
    """4UHX site_2: one Mo, plus an MTE ligand residue in the same site."""
    site = _Site("site_2", [
        _Atom("A", 3004, "MOS", "MO"),
        _Atom("A", 3004, "MOS", "S"),
        _Atom("A", 3003, "MTE", "C"),
    ])
    mos = _Res("MOS", "A", 3004, site, "metal_site")
    mte = _Res("MTE", "A", 3003, site, "metal_site")
    return {"members": [mos, mte], "metals": [mos], "ligands": [mte],
            "category": "metal_site", "site_id": "site_2"}


# --------------------------------------------------------------------------- #
# nuclearity
# --------------------------------------------------------------------------- #

def test_fe2s2_is_binuclear_not_mononuclear():
    name, category, status = _param()._format_unit(_fes_unit())

    assert "binuclear" in status, status
    assert "mononuclear" not in status
    assert category == "Metal Site"
    assert "FES(A:3001)" in name


def test_single_metal_site_is_still_mononuclear():
    name, _category, status = _param()._format_unit(_mos_unit())

    assert "mononuclear" in status, status
    assert "+ligand" in status
    assert "MTE" in name


def test_fe4s4_is_tetranuclear():
    site = _Site("site_1", [_Atom("A", 500, "SF4", "FE") for _ in range(4)]
                 + [_Atom("A", 500, "SF4", "S") for _ in range(4)])
    sf4 = _Res("SF4", "A", 500, site, "metal_site")
    unit = {"members": [sf4], "metals": [sf4], "ligands": [],
            "category": "metal_site", "site_id": "site_1"}

    _name, _category, status = _param()._format_unit(unit)

    assert "tetranuclear" in status, status


def test_metal_atom_count_reads_atoms_not_residues():
    p = _param()
    fes = _fes_unit()["metals"][0]
    mos = _mos_unit()["metals"][0]

    assert p._count_metal_atoms(fes) == 2   # FE1 + FE2, one residue
    assert p._count_metal_atoms(mos) == 1   # the S is not a metal


# --------------------------------------------------------------------------- #
# site counting
# --------------------------------------------------------------------------- #

def test_ligand_member_does_not_count_as_its_own_site():
    """The 4UHX case: 3 sites, 4 residues stamped metal_site."""
    p = _param()
    residues = _fes_unit()["members"] + _mos_unit()["members"] + [
        _Res("FES", "A", 3002, _Site("site_3", [_Atom("A", 3002, "FES", "FE"),
                                                _Atom("A", 3002, "FES", "FE")]),
             "metal_site"),
    ]
    assert len(residues) == 4, "fixture should reproduce the 4-residue layout"

    n_sites = len({p._unit_key(r) for r in residues
                   if r.category == "metal_site"})

    assert n_sites == 3


def test_standalone_residues_are_their_own_units():
    p = _param()
    a = _Res("FAD", "A", 3006, None, "small_molecule")
    b = _Res("MLI", "A", 3010, None, "small_molecule")
    c = _Res("MLI", "A", 3011, None, "small_molecule")

    assert len({p._unit_key(r) for r in (a, b, c)}) == 3, \
        "two MLI copies at different resids are two units"


# --------------------------------------------------------------------------- #
# menu suggestion numbering
# --------------------------------------------------------------------------- #

class _WS:
    def __init__(self, d):
        self.d = d

    def get(self, k, default=None):
        return self.d.get(k, default)


def _suggestion(residues, pending=None):
    p = _param()
    return p.get_menu_suggestion(_WS({
        "non_standard_residues": residues,
        "pending_parameterizations": pending or {},
    }))


def test_suggestion_points_at_the_actual_menu_options():
    residues = _fes_unit()["members"] + _mos_unit()["members"]
    text = _suggestion(residues)

    # 3 parameterize, 4 status, 5 help — matching get_menu_options.
    assert "option 3" in text and "option 4" in text and "option 5" in text
    assert "option 2" not in text, f"stale pre-import numbering: {text}"


def test_suggestion_counts_sites_not_residues():
    residues = _fes_unit()["members"] + _mos_unit()["members"]
    assert len(residues) == 3          # FES, MOS, MTE

    text = _suggestion(residues)

    assert "2 metal sites" in text, text


def test_pending_suggestion_points_at_parameterize():
    text = _suggestion([], pending={"FES": {}})

    assert "option 3" in text and "option 4" in text
    assert "option 2" not in text, f"stale numbering: {text}"


def test_singular_site_is_not_pluralised():
    text = _suggestion(_fes_unit()["members"])

    assert "1 metal site." in text or "1 metal site," in text, text
