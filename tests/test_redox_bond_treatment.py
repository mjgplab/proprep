"""Regression tests for the RedoxSiteBond.treatment field.

`treatment` distinguishes how a metal-ligand contact is REALIZED in the force
field ("bonded" = MCPB bonded model; "restrained" = kept as a nonbonded residue
held by an MD distance restraint), independently of what the contact IS
chemically (`chemical_type`). The field must:

  * default to "bonded" so all existing behavior and already-serialized sites
    are unchanged (the change is additive/non-breaking), and
  * round-trip through export -> import and through a legacy dict that predates
    the field (missing key -> "bonded").
"""

import io
import json
import os
import tempfile

from proprep.structure_prep.comprehensive_redox_detector import (
    RedoxSite,
    RedoxSiteAtom,
    RedoxSiteBond,
    dict_to_redox_site,
    _export_to_json,
)


def _site_with_two_waters():
    """A metal (MN 185) coordinated by two waters: 183 restrained, 184 bonded."""
    site = RedoxSite("site_1", "test")
    mn = (0.0, 0.0, 0.0)
    o183 = (2.2, 0.0, 0.0)
    o184 = (-2.2, 0.0, 0.0)
    site.add_atom(RedoxSiteAtom("A", "MN", 185, "MN", mn, "Mn"))
    site.add_atom(RedoxSiteAtom("A", "WAT", 183, "O", o183, "O"))
    site.add_atom(RedoxSiteAtom("A", "WAT", 184, "O", o184, "O"))
    site.add_bond_with_classification(mn, o183, 2.2, treatment="restrained")
    site.add_bond_with_classification(mn, o184, 2.2, treatment="bonded")
    return site, mn, o183, o184


def _bond_dict(treatment=None):
    d = {
        "atom1": {"chain": "A", "resname": "MN", "resid": 1,
                  "atom_name": "MN", "element": "Mn"},
        "atom2": {"chain": "A", "resname": "HOH", "resid": 2,
                  "atom_name": "O", "element": "O"},
        "bond_type": "interresidue",
        "chemical_type": "coordinate",
        "distance": 2.2,
        "atom1_element": "Mn",
        "atom2_element": "O",
        "atom1_coordinates": [0, 0, 0],
        "atom2_coordinates": [1, 0, 0],
    }
    if treatment is not None:
        d["treatment"] = treatment
    return d


def test_default_is_bonded():
    b = RedoxSiteBond((0, 0, 0), (1, 0, 0), "interresidue", "coordinate",
                      2.1, "Mn", "O")
    assert b.treatment == "bonded"


def test_explicit_restrained_is_honored():
    b = RedoxSiteBond((0, 0, 0), (1, 0, 0), "interresidue", "coordinate",
                      2.1, "Mn", "O", treatment="restrained")
    assert b.treatment == "restrained"


def test_legacy_dict_without_treatment_defaults_to_bonded():
    """A site JSON written before this field existed must still load."""
    site = dict_to_redox_site(
        {"site_id": "site_1", "centers": [], "atoms": [],
         "bonds": [_bond_dict(treatment=None)]}
    )
    assert site.bonds[0].treatment == "bonded"


def test_treatment_survives_export_import_round_trip():
    site = dict_to_redox_site(
        {"site_id": "site_1", "centers": [], "atoms": [],
         "bonds": [_bond_dict(treatment="restrained")]}
    )
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as d:
        try:
            os.chdir(d)
            _export_to_json([site], "dummy.pdb")  # writes dummy_redox_sites.json
            data = json.load(open("dummy_redox_sites.json"))
        finally:
            os.chdir(cwd)
    assert data["sites"][0]["bonds"][0]["treatment"] == "restrained"
    reimported = dict_to_redox_site(data["sites"][0])
    assert reimported.bonds[0].treatment == "restrained"


# ---------------------------------------------------------------------------
# §3 gate integration: a restrained metal-water contact must be excluded from
# the bonded model (LINK record, tLEaP bond, RESP metal-site set) while the
# bonded one is kept.
# ---------------------------------------------------------------------------

def test_collect_restrained_ligands_identifies_only_restrained_water():
    from proprep.forcefield_prep.metal_site_parameterizer import _collect_restrained_ligands
    site, mn, o183, o184 = _site_with_two_waters()
    coords, resids = _collect_restrained_ligands(site)
    assert resids == {183}                         # only the restrained water
    assert tuple(round(x, 3) for x in o183) in coords
    assert tuple(round(x, 3) for x in o184) not in coords  # bonded water excluded


def test_link_records_skip_restrained_bond():
    from proprep.forcefield_prep.mcpb.fingerprint_generator import FingerprintGenerator
    site, mn, o183, o184 = _site_with_two_waters()
    coord_to_id = {mn: 1, o183: 2, o184: 3}
    buf = io.StringIO()
    FingerprintGenerator()._write_link_records(buf, site, coord_to_id)
    out = buf.getvalue()
    assert "1-MN 3-O" in out          # bonded water 184 -> LINK present
    assert "1-MN 2-O" not in out      # restrained water 183 -> no LINK


def test_tleap_bond_commands_skip_restrained_bond():
    from rich.console import Console
    from proprep.forcefield_prep.structure_preprocessor import StructurePreprocessor
    site, mn, o183, o184 = _site_with_two_waters()
    pre = StructurePreprocessor(processor=None, console=Console())
    cmds = pre._convert_site_bonds_to_tleap_commands([site])
    assert any("mol.184.O" in c for c in cmds)      # bonded water -> bond command
    assert not any("mol.183.O" in c for c in cmds)  # restrained water -> skipped


def test_mcpb4_naming_excludes_restrained_water():
    """mcpb-4 must not rename a restrained water: it has no mol2/library unit,
    so a phantom WT* name in the PDB would have no tLEaP unit to build from.
    This mirrors the exclusion filter in _checklist_mcpb_4_integration.
    """
    from proprep.forcefield_prep.metal_site_parameterizer import _collect_restrained_ligands
    from proprep.forcefield_prep.mcpb.integration_utils import generate_unique_residue_names

    site, mn, o183, o184 = _site_with_two_waters()

    # Fingerprint-style residue keys (int resid, resname) for every site residue,
    # including BOTH waters — as the fingerprint actually contains them.
    combined_residue_keys = [(185, "MN"), (183, "WAT"), (184, "WAT")]

    restrained_resids = set()
    _, _resids = _collect_restrained_ligands(site)
    for _r in _resids:
        restrained_resids.add(int(_r))

    filtered = [k for k in combined_residue_keys if k[0] not in restrained_resids]
    name_map = generate_unique_residue_names(filtered)

    assert (183, "WAT") not in name_map     # restrained water: not renamed
    assert (184, "WAT") in name_map         # bonded water: still gets an MCPB name
    assert (185, "MN") in name_map


def test_global_types_path_excludes_restrained_water_from_renaming():
    """The multi-site global-types path must also skip Y-renaming for a
    restrained water, matching coords tolerantly (rounding-safe)."""
    from rich.console import Console
    from proprep.forcefield_prep.metal_site_parameterizer import MetalSiteWorkflowManager
    from proprep.forcefield_prep.mcpb.global_atom_registry import GlobalTypeAssignment

    site, mn, o183, o184 = _site_with_two_waters()
    mgr = MetalSiteWorkflowManager(console=Console())
    mgr.provided_redox_site = site

    def gta(coords, renamed):
        return GlobalTypeAssignment(
            site_index=0, coords=coords, chain="A", resname="WAT", resid=0,
            atom_name="O", element="O", original_type="OW",
            global_renamed_type=renamed, is_metal=False, is_ligand=True,
            site_id="site_1",
        )

    # global registry flagged BOTH water O's as ligands (is_ligand=True)
    global_types = {o183: gta(o183, "YA"), o184: gta(o184, "YB")}
    ta = mgr._convert_global_types(global_types)

    assert ta[o183]["is_metal_ligand"] is False   # restrained -> NOT renamed
    assert ta[o184]["is_metal_ligand"] is True     # bonded -> still a ligand
