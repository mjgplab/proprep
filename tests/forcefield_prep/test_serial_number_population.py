"""
Regression test for _populate_serial_numbers coordinate matching.

RedoxSite atom coords can come back from a persisted workflow_state as
full-precision *strings* ('30.838055'). The population matched atoms to the PDB
via a dict keyed by floats rounded to 3 dp, comparing the atom's RAW coords —
so a string coord never matched, no serial got stamped, and the MCPB fingerprint
silently fell back to synthetic 10000+ atom IDs, breaking bond/angle extraction
(0 bonds -> Seminario has nothing to do).

Fix: normalize both sides with round(float(v), 3).
"""

from types import SimpleNamespace
import logging

from proprep.forcefield_prep.metal_site_parameterizer import MetalSiteWorkflowManager
from proprep.structure_prep.comprehensive_redox_detector import RedoxSiteAtom


class _Console:
    def print(self, *a, **k):
        pass


def _wf():
    wf = object.__new__(MetalSiteWorkflowManager)
    wf.console = _Console()
    wf.logger = logging.getLogger("test_serial_pop")
    return wf


def _pdb(tmp_path):
    # Two atoms at known coords / serials.
    p = tmp_path / "prepared.pdb"
    p.write_text(
        "ATOM    678  N   HID A  45      30.838  -4.626  14.299  1.00  0.00           N\n"
        "ATOM   3014 MN   MN  A 185      27.166  -4.943  19.905  1.00  0.00          MN\n"
        "END\n"
    )
    return str(p)


def _atom(coords):
    return RedoxSiteAtom(
        chain="A", resname="HID", resid=45, atom_name="N",
        coords=coords, element="N",
    )


def test_string_coords_still_populate_serials(tmp_path):
    pdb = _pdb(tmp_path)
    # Full-precision STRING coords, as restored from workflow_state.json.
    site = SimpleNamespace(atoms=[
        _atom(("30.838055", "-4.6255007", "14.299155")),
        _atom(("27.166", "-4.943", "19.905")),
    ])
    wf = _wf()
    wf._populate_serial_numbers(site, pdb)

    assert site.atoms[0].properties.get("serial_number") == 678
    assert site.atoms[1].properties.get("serial_number") == 3014


def test_float_coords_still_work(tmp_path):
    pdb = _pdb(tmp_path)
    # Full-precision FLOAT coords (fresh detection, pre-serialization).
    site = SimpleNamespace(atoms=[
        _atom((30.838055, -4.6255007, 14.299155)),
        _atom((27.166, -4.943, 19.905)),
    ])
    wf = _wf()
    wf._populate_serial_numbers(site, pdb)

    assert site.atoms[0].properties.get("serial_number") == 678
    assert site.atoms[1].properties.get("serial_number") == 3014


def test_unmatched_atom_left_without_serial(tmp_path):
    pdb = _pdb(tmp_path)
    site = SimpleNamespace(atoms=[_atom(("999.0", "999.0", "999.0"))])
    wf = _wf()
    wf._populate_serial_numbers(site, pdb)
    assert "serial_number" not in site.atoms[0].properties
