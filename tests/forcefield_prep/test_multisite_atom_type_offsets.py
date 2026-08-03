"""
Regression test: M*/Y* atom-type names must be globally unique across multiple
metal sites in one protein.

Every metal site's frcmod/mol2 is loaded into the SAME tLEaP session, so if two
sites each restart their renaming at M1/Y1 the second site's types silently
overwrite the first's. The fix threads a running offset out of
_apply_systematic_renaming (stashed as type_offset_metal_end /
type_offset_ligand_end) and back into the next site's metal_start/ligand_start.
"""

from proprep.forcefield_prep.metal_site_parameterizer import MetalSiteWorkflowManager


def _assignment(is_center=False, is_metal_ligand=False):
    return {
        "is_center": is_center,
        "is_metal_ligand": is_metal_ligand,
        "renamed_type": None,
        "renamed": False,
    }


def _make_site_assignments(n_metals, n_ligands):
    """Build a coords->assignment dict with n_metals centers + n_ligands ligands."""
    assignments = {}
    i = 0
    for _ in range(n_metals):
        assignments[(float(i), 0.0, 0.0)] = _assignment(is_center=True)
        i += 1
    for _ in range(n_ligands):
        assignments[(float(i), 0.0, 0.0)] = _assignment(is_metal_ligand=True)
        i += 1
    return assignments


def test_offsets_advance_across_sites():
    wf = MetalSiteWorkflowManager(console=None)

    # --- Site 1: 1 metal, 3 ligands ---
    site1 = _make_site_assignments(n_metals=1, n_ligands=3)
    wf._apply_systematic_renaming(redox_site=None, type_assignments=site1,
                                  metal_start=0, ligand_start=0)
    site1_metals = sorted(a["renamed_type"] for a in site1.values() if a["is_center"])
    site1_ligands = sorted(a["renamed_type"] for a in site1.values() if a["is_metal_ligand"])
    assert site1_metals == ["M1"]
    assert site1_ligands == ["Y1", "Y2", "Y3"]
    assert wf.type_offset_metal_end == 1
    assert wf.type_offset_ligand_end == 3

    # --- Site 2: continue from site 1's ending offsets ---
    metal_off = wf.type_offset_metal_end
    ligand_off = wf.type_offset_ligand_end
    site2 = _make_site_assignments(n_metals=1, n_ligands=2)
    wf._apply_systematic_renaming(redox_site=None, type_assignments=site2,
                                  metal_start=metal_off, ligand_start=ligand_off)
    site2_metals = sorted(a["renamed_type"] for a in site2.values() if a["is_center"])
    site2_ligands = sorted(a["renamed_type"] for a in site2.values() if a["is_metal_ligand"])
    assert site2_metals == ["M2"]
    assert site2_ligands == ["Y4", "Y5"]

    # --- The whole point: no name collides across the two sites ---
    all_names = site1_metals + site1_ligands + site2_metals + site2_ligands
    assert len(all_names) == len(set(all_names)), f"collision: {all_names}"


def test_default_offsets_restart_at_one():
    """With no offsets passed, numbering starts at M1/Y1 (single-site behavior)."""
    wf = MetalSiteWorkflowManager(console=None)
    site = _make_site_assignments(n_metals=1, n_ligands=2)
    wf._apply_systematic_renaming(redox_site=None, type_assignments=site)
    metals = sorted(a["renamed_type"] for a in site.values() if a["is_center"])
    ligands = sorted(a["renamed_type"] for a in site.values() if a["is_metal_ligand"])
    assert metals == ["M1"]
    assert ligands == ["Y1", "Y2"]
