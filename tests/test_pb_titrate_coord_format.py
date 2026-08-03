"""Regression: pb_titrate must write inpcrd coordinates that fit Amber's
fixed-column ``(6F12.7)`` format, so pbsa's formatted read never SIGSEGVs.

Background: large multi-chain assemblies can place a sliced PB cluster's
atoms past ±1000 Å. ParmEd's ``%12.7f`` writer then spills a 13th column
(e.g. ``-1234.5678901``), collapsing the inter-field gap; pbsa reads
misaligned tokens and gfortran crashes in ``formatted_transfer_scalar_read``
*before* any solve — deterministic, so every retry repeats it. ``run_pbsa``
recenters the cluster (a translation, hence PB-energy-invariant) so coords
stay small and in-range, EXCEPT on the phiout path where the grid frame must
match the original-frame coordinates coupling.py interpolates against.
"""
import numpy as np
import parmed as pmd
import pytest

from proprep.pb_titrate.pb_backend import _recenter_for_format

# A clean inpcrd coordinate line is exactly 6 fields x 12 columns.
MAX_CLEAN_WIDTH = 72


def _build_structure(n_atoms: int) -> pmd.Structure:
    s = pmd.Structure()
    res = pmd.Residue("AS4", number=0)
    for i in range(n_atoms):
        a = pmd.Atom(name=f"A{i}", atomic_number=6)
        s.add_atom(a, res.name, res.number)
    return s


def _max_coord_line_width(structure: pmd.Structure, tmp_path) -> int:
    out = tmp_path / "test.inpcrd"
    structure.save(str(out), overwrite=True, format="rst7")
    # line 0: title, line 1: natom; coordinate rows follow
    rows = out.read_text().splitlines()[2:]
    return max(len(r) for r in rows)


def test_negative_far_coords_overflow_without_recenter(tmp_path):
    """Sanity-check the failure mode the fix targets actually exists."""
    s = _build_structure(6)
    # ~ -1500 Å on every axis: 13 columns each, fields collide
    s.coordinates = np.full((6, 3), -1500.123456, dtype=np.float64)
    assert _max_coord_line_width(s, tmp_path) > MAX_CLEAN_WIDTH


def test_recenter_brings_far_coords_into_format(tmp_path):
    s = _build_structure(6)
    s.coordinates = np.random.RandomState(0).rand(6, 3) * 40.0 - 1520.0
    _recenter_for_format(s)
    assert _max_coord_line_width(s, tmp_path) <= MAX_CLEAN_WIDTH
    # all-positive, small-magnitude after the shift
    assert s.coordinates.min() >= 0.0
    assert s.coordinates.max() < 1000.0


def test_recenter_preserves_internal_geometry(tmp_path):
    """A pure translation: all pairwise distances are unchanged, so PB
    energies (and the cluster-minus-model ΔΔG) are invariant."""
    s = _build_structure(8)
    xyz = np.random.RandomState(1).rand(8, 3) * 30.0 - 1510.0
    s.coordinates = xyz.copy()
    before = xyz - xyz[0]
    _recenter_for_format(s)
    after = s.coordinates - s.coordinates[0]
    np.testing.assert_allclose(before, after, atol=1e-6)


def test_recenter_noop_for_inrange_coords():
    """In-range coordinates are still translated to the +margin corner, but
    geometry is preserved and values stay tiny — never makes things worse."""
    s = _build_structure(4)
    xyz = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0],
                    [-2.0, 0.0, 1.0], [3.0, -1.0, 2.0]])
    s.coordinates = xyz.copy()
    _recenter_for_format(s)
    assert s.coordinates.min() >= 0.0
    assert s.coordinates.max() < 1000.0


def test_recenter_handles_empty_coords():
    s = pmd.Structure()
    # no atoms -> coordinates is None/empty; must not raise
    _recenter_for_format(s)


# ---------------------------------------------------------------------------
# Box strip: a focused FD-PB calc is non-periodic, and a sliced cluster cut
# from a solvated box keeps IFBOX=1 with many disconnected molecules, which
# triggers pbsa's 1-slot ATOMS_PER_MOLECULE overrun (SIGSEGV in rdparm2's
# formatted read). run_pbsa sets structure.box = None so pbsa skips that read.
# These tests verify the assumption the fix relies on: box=None clears IFBOX.
# ---------------------------------------------------------------------------
import os
import proprep.pb_titrate as _pbt

_BUNDLED_PRMTOP = os.path.join(
    os.path.dirname(_pbt.__file__), "tests", "bpti_multi", "coupling_work",
    "self", "AS40003", "cluster_state0.prmtop")


def _ifbox(parm) -> int:
    return parm.parm_data["POINTERS"][27]


@pytest.mark.skipif(not os.path.exists(_BUNDLED_PRMTOP),
                    reason="bundled test prmtop not installed")
def test_box_none_clears_ifbox(tmp_path):
    parm = pmd.load_file(_BUNDLED_PRMTOP)
    # give it a periodic box -> IFBOX becomes 1 on save (mimics a solvated slice)
    parm.box = [80.0, 80.0, 80.0, 90.0, 90.0, 90.0]
    boxed = tmp_path / "boxed.prmtop"
    parm.save(str(boxed), overwrite=True)
    assert _ifbox(pmd.load_file(str(boxed))) == 1

    # the fix: dropping the box clears IFBOX, so pbsa skips SOLVENT_POINTERS /
    # ATOMS_PER_MOLECULE (the read that overruns for nspm >= 2).
    parm.box = None
    nobox = tmp_path / "nobox.prmtop"
    parm.save(str(nobox), overwrite=True)
    reloaded = pmd.load_file(str(nobox))
    assert _ifbox(reloaded) == 0
    txt = open(str(nobox)).read()
    assert "ATOMS_PER_MOLECULE" not in txt
    assert "SOLVENT_POINTERS" not in txt
