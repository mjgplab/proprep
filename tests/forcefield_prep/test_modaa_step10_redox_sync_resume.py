"""
Regression: modified-AA Route B step 10 (FF integration) syncs RedoxSite objects
to the renamed PDB even after a session RESUME.

detected_redox_sites round-trips through workflow_state.json. The checklist
serializer wraps complex objects as {"__type__": "RedoxSite", "value": {...}} and
the deserializer deliberately leaves them as dicts for the caller to reconstruct
(see _deserialize_value). So on a resumed run the workspace hands step 10 dict-
form sites, and sync_redox_sites_from_pdb — which reads site.atoms — used to fail
with "'dict' object has no attribute 'atoms'", skipping the sync and leaving
tLEaP bond commands pointing at pre-rename atom names.

Step 10 now normalizes via _ensure_redox_site_objects before syncing. These tests
pin the full serialize -> deserialize -> normalize -> sync chain on a synthetic
covalent adduct (a Cys SG covalently bonded to a ligand C, renamed to lib names).
"""

import pytest

from proprep.utils.workflow_checklist import _serialize_value, _deserialize_value
from proprep.forcefield_prep.structure_preprocessor import _ensure_redox_site_objects
from proprep.structure_prep.comprehensive_redox_detector import (
    RedoxSite, RedoxSiteAtom, RedoxSiteBond, sync_redox_sites_from_pdb,
)

# (element, resname, resid, atom_name, x, y, z) — coords are the permanent key.
_SG = ("S", "CYS", 12, "SG", 1.000, 2.000, 3.000)
_C25 = ("C", "MOV", 303, "C25", 1.500, 2.500, 3.500)


def _build_site():
    site = RedoxSite("site_1", "t")
    for el, rn, ri, an, x, y, z in (_SG, _C25):
        site.add_atom(RedoxSiteAtom(chain="A", resname=rn, resid=ri,
                                    atom_name=an, coords=(x, y, z), element=el))
    site.bonds.append(RedoxSiteBond(
        atom1_coords=_SG[4:], atom2_coords=_C25[4:],
        bond_type="interresidue", chemical_type="covalent", distance=1.8,
        atom1_element="S", atom2_element="C",
        atom1_residue_info={"chain": "A", "resname": "CYS", "resid": 12,
                            "atom_name": "SG"},
        atom2_residue_info={"chain": "A", "resname": "MOV", "resid": 303,
                            "atom_name": "C25"},
    ))
    return site


def _pdb_atom(serial, name, resname, resid, x, y, z, element):
    return (f"ATOM  {serial:>5} {name:<4} {resname:>3} A{resid:>4}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}\n")


def _renamed_pdb(tmp_path):
    """Same coords as the site, but post-rename residue/atom names (CS1 / MV1)."""
    p = tmp_path / "renamed.pdb"
    p.write_text(
        _pdb_atom(1, "S1", "CS1", 12, *_SG[4:], "S")
        + _pdb_atom(2, "C6", "MV1", 303, *_C25[4:], "C")
    )
    return str(p)


def _resumed_sites():
    """A site as it re-emerges from workflow_state.json on resume: wrapper dicts."""
    restored = _deserialize_value(_serialize_value([_build_site()]))
    # Precondition: this is the broken form — a dict, not a RedoxSite.
    assert isinstance(restored[0], dict)
    assert not hasattr(restored[0], "atoms")
    return restored


def test_raw_resumed_sites_would_break_sync(tmp_path):
    # Documents the bug: syncing the un-normalized resume form raises exactly the
    # error the user saw. This is what the caller must prevent.
    with pytest.raises(AttributeError):
        sync_redox_sites_from_pdb(_renamed_pdb(tmp_path), _resumed_sites())


def test_normalize_then_sync_rekeys_after_resume(tmp_path):
    sites = _ensure_redox_site_objects(_resumed_sites())
    assert isinstance(sites[0], RedoxSite)
    assert len(sites[0].atoms) == 2

    changed = sync_redox_sites_from_pdb(_renamed_pdb(tmp_path), sites)
    assert changed is True

    by_coords = {a.coords: a for a in sites[0].atoms}
    sg = by_coords[_SG[4:]]
    c25 = by_coords[_C25[4:]]
    # Atoms re-keyed to the renamed residue/atom names (bond commands stay valid).
    assert (sg.resname, sg.atom_name) == ("CS1", "S1")
    assert (c25.resname, c25.atom_name) == ("MV1", "C6")


def test_object_form_passes_through_unchanged(tmp_path):
    # A fresh (non-resumed) run hands objects straight through; normalization is
    # a no-op and the sync still re-keys.
    sites = _ensure_redox_site_objects([_build_site()])
    assert isinstance(sites[0], RedoxSite)
    assert sync_redox_sites_from_pdb(_renamed_pdb(tmp_path), sites) is True
    assert {a.resname for a in sites[0].atoms} == {"CS1", "MV1"}


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-v"]))
