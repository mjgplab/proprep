"""
The large model's MK ReadRadii block must be populated for cluster sites too.

large_resp.gjf was written with a bare ``Pop=MK`` and no radii block, so the
Merz-Kollman ESP sampled the metals with Gaussian's defaults instead of the
radii the metal database supplies.

The radii were collected by walking ``redox_site.centers`` and looking each
center's ``coords`` up in ``type_assignments``. That holds only for a lone
``metal_ion``, where the center IS the metal atom. For an organometallic
cofactor or a pure cluster the center describes the RESIDUE: in the reported
structure its ``element`` was None and its ``coords`` empty, so nothing matched
and the dict came back empty — silently, since an empty dict just selects the
``Pop=MK`` spelling.

They now come from the metal atoms in ``type_assignments``, which already carry
``is_center``, the element, and the charge the user entered.
"""

import io

import pytest
from rich.console import Console

from proprep.forcefield_prep.metal_site_parameterizer import MetalSiteWorkflowManager


def _mgr():
    m = MetalSiteWorkflowManager.__new__(MetalSiteWorkflowManager)
    m.console = Console(file=io.StringIO(), width=200)
    m.logger = None
    return m


def _atom(element, charge, is_center):
    return {"element": element, "charge": charge, "is_center": is_center}


def test_fe2s2_cluster_resolves_its_iron_radius():
    """A pure cluster: the reported case."""
    m = _mgr()
    assignments = {
        (0.0, 0.0, 0.0): _atom("FE", 3.0, True),
        (1.0, 0.0, 0.0): _atom("FE", 3.0, True),
        (2.0, 0.0, 0.0): _atom("S", -2.0, False),
        (3.0, 0.0, 0.0): _atom("S", -2.0, False),
    }

    radii = m._collect_metal_radii(assignments)

    assert set(radii) == {"Fe"}
    # MCPB.py's own MK table (Smith et al. JCTC 2023, 19, 2064), not the IOD
    # force-field value — see test_esp_radius_is_not_the_forcefield_radius.
    assert radii["Fe"] == pytest.approx(1.383, abs=1e-3)


def test_moco_resolves_its_molybdenum_radius():
    m = _mgr()
    assignments = {
        (0.0, 0.0, 0.0): _atom("MO", 6.0, True),
        (1.0, 0.0, 0.0): _atom("S", -2.0, False),
        (2.0, 0.0, 0.0): _atom("O", -2.0, False),
    }

    radii = m._collect_metal_radii(assignments)

    assert set(radii) == {"Mo"}
    assert radii["Mo"] > 0


def test_element_case_is_normalised():
    """PDB elements are upper case; the database keys are title case."""
    m = _mgr()

    assert set(m._collect_metal_radii({(0.0, 0.0, 0.0): _atom("ZN", 2.0, True)})) == {"Zn"}


def test_non_centre_atoms_are_ignored():
    m = _mgr()
    assignments = {
        (0.0, 0.0, 0.0): _atom("S", -2.0, False),
        (1.0, 0.0, 0.0): _atom("O", -2.0, False),
    }

    assert m._collect_metal_radii(assignments) == {}


def test_a_metal_with_no_charge_yet_is_skipped_not_crashed():
    """charge is None until the user is asked."""
    m = _mgr()

    assert m._collect_metal_radii({(0.0, 0.0, 0.0): _atom("FE", None, True)}) == {}


def test_empty_result_says_so():
    """An empty dict silently selects 'Pop=MK'; that must be visible."""
    m = _mgr()

    m._collect_metal_radii({(0.0, 0.0, 0.0): _atom("S", -2.0, False)})

    assert "No metal radii resolved" in m.console.file.getvalue()


def test_same_element_at_two_charges_is_reported():
    """ReadRadii is per element, so the second oxidation state cannot differ."""
    m = _mgr()
    assignments = {
        (0.0, 0.0, 0.0): _atom("FE", 3.0, True),
        (1.0, 0.0, 0.0): _atom("FE", 2.0, True),
    }

    radii = m._collect_metal_radii(assignments)

    assert set(radii) == {"Fe"}, "one radius per element"
    out = m.console.file.getvalue()
    assert "ReadRadii is per element" in out
    assert "+3" in out and "+2" in out


def test_two_different_metals_both_appear():
    m = _mgr()
    assignments = {
        (0.0, 0.0, 0.0): _atom("FE", 3.0, True),
        (1.0, 0.0, 0.0): _atom("ZN", 2.0, True),
    }

    assert set(m._collect_metal_radii(assignments)) == {"Fe", "Zn"}


def test_object_style_assignments_are_supported():
    """type_assignments values are dicts on one route and objects on another."""
    class _A:
        def __init__(self, element, charge, is_center):
            self.element, self.charge, self.is_center = element, charge, is_center

    m = _mgr()

    radii = m._collect_metal_radii({(0.0, 0.0, 0.0): _A("FE", 3.0, True)})

    assert set(radii) == {"Fe"}


def test_esp_radius_is_not_the_forcefield_radius():
    """Two different quantities from two different sources.

    The MK ReadRadii value decides how close to the nucleus ESP grid points may
    fall; the IOD value is a Lennard-Jones parameter for the frcmod. They must
    not be conflated — pointing the ESP lookup at MCPB.py's table must leave the
    force-field parameters alone.
    """
    from proprep.forcefield_prep.mcpb.metal_ion_database import MetalIonDatabase

    db = MetalIonDatabase(water_model="tip3p")

    esp = db.get_vdw_radius("Fe", 3)
    ff_radius, _eps, ff_source = db._get_vdw_params("Fe", 3)

    assert esp == pytest.approx(1.383, abs=1e-3), "ESP radius: MCPB.py table"
    assert ff_radius == pytest.approx(1.386, abs=1e-3), "force field: IOD, unchanged"
    assert "IOD" in ff_source


def test_mo6_resolves_from_the_mcpb_table():
    """The reported gap: no IOD table runs past tetravalent, so Mo(VI) fell to
    an unsourced generic 1.5."""
    from proprep.forcefield_prep.mcpb.metal_ion_database import MetalIonDatabase

    radius = MetalIonDatabase(water_model="tip3p").get_vdw_radius("Mo", 6)

    assert radius == pytest.approx(1.337, abs=1e-3)
    assert radius != 1.5, "should no longer be the generic fallback"


def test_readradii_is_selected_only_when_radii_exist():
    """The spelling in the route line is driven by the dict being non-empty."""
    m = _mgr()
    empty = m._collect_metal_radii({(0.0, 0.0, 0.0): _atom("S", -2.0, False)})
    filled = m._collect_metal_radii({(0.0, 0.0, 0.0): _atom("FE", 3.0, True)})

    assert ("Pop(MK,ReadRadii)" if filled else "Pop=MK") == "Pop(MK,ReadRadii)"
    assert ("Pop(MK,ReadRadii)" if empty else "Pop=MK") == "Pop=MK"
