"""
Regression tests for the ONIOM single-bond backbone cut + HX cap parameters.

Covers the redesign of the full-residue boundary (single-bond C-cut/N-cut,
flanking Cα left in LOW, one HX link per cut) and the GAFF-supplied formamide
cap parameters (carbonyl C–H, absent from ff14SB).
"""
import io
import logging
import re
from collections import namedtuple

import parmed as pmd
import pytest

import proprep.oniom_prep.oniom_qmmm_preparator as M
import proprep.oniom_prep.oniom_writer as W
from proprep.oniom_prep.data_structures import ONIOMLayer


# ── helpers ──────────────────────────────────────────────────────────────
def _preparator_cls():
    src = open(M.__file__).read()
    name = next(c for c in re.findall(r"^class (\w+)", src, re.M)
                if hasattr(getattr(M, c, None), "_detect_boundary_atoms"))
    return getattr(M, name)


def _build_peptide(n):
    """An n-residue poly-ALA with backbone + Cβ and peptide bonds."""
    s = pmd.Structure()
    for ri in range(n):
        for an in ("N", "H", "CA", "HA", "CB", "C", "O"):
            s.add_atom(pmd.Atom(name=an, type=an), "ALA", ri + 1)
    byres = [{a.name: a for a in r.atoms} for r in s.residues]
    for r in byres:
        for x, y in (("N", "H"), ("N", "CA"), ("CA", "HA"),
                     ("CA", "CB"), ("CA", "C"), ("C", "O")):
            s.bonds.append(pmd.Bond(r[x], r[y]))
    for i in range(n - 1):
        s.bonds.append(pmd.Bond(byres[i]["C"], byres[i + 1]["N"]))
    return s


def _run_boundary(n, high_residues):
    s = _build_peptide(n)
    cls = _preparator_cls()
    m = cls.__new__(cls)
    m._parm = s
    m.boundary_atom_indices = []
    m.high_residue_indices = list(high_residues)
    m.high_atom_indices = set()
    m.layer_membership = {a.idx: ONIOMLayer.LOW for a in s.atoms}
    for ri in high_residues:
        for a in s.residues[ri].atoms:
            m.layer_membership[a.idx] = ONIOMLayer.HIGH

    class _Console:
        def print(self, *a, **k):
            pass

    class _P:
        console = _Console()

    m.processor = _P()
    m._detect_boundary_atoms()
    return s, m


def _layer(m, s, res_idx, atom_name):
    a = next(a for a in s.residues[res_idx].atoms if a.name == atom_name)
    return m.layer_membership[a.idx]


# ── boundary cut tests ───────────────────────────────────────────────────
def test_single_bond_cut_leaves_flanking_ca_low():
    """Middle residue in QM: each flanking Cα stays LOW and is the boundary."""
    s, m = _run_boundary(3, [1])

    # Flanking Cα are the boundary atoms, both in LOW.
    boundary_names = {(s.atoms[i].residue.idx, s.atoms[i].name)
                      for i in m.boundary_atom_indices}
    assert boundary_names == {(0, "CA"), (2, "CA")}
    assert _layer(m, s, 0, "CA") == ONIOMLayer.LOW
    assert _layer(m, s, 2, "CA") == ONIOMLayer.LOW


def test_c_cut_promotes_only_carbonyl():
    """N-terminal flank (prev): C-cut promotes {C,O} only — formamide cap."""
    s, m = _run_boundary(3, [1])
    assert _layer(m, s, 0, "C") == ONIOMLayer.HIGH
    assert _layer(m, s, 0, "O") == ONIOMLayer.HIGH
    # Cα/Cβ/N of the flank stay LOW (single bond cut, not acetamide).
    for nm in ("N", "CA", "CB"):
        assert _layer(m, s, 0, nm) == ONIOMLayer.LOW


def test_n_cut_promotes_only_amide():
    """C-terminal flank (next): N-cut promotes {N,H} only — primary-amide cap."""
    s, m = _run_boundary(3, [1])
    assert _layer(m, s, 2, "N") == ONIOMLayer.HIGH
    assert _layer(m, s, 2, "H") == ONIOMLayer.HIGH
    for nm in ("C", "CA", "CB"):
        assert _layer(m, s, 2, nm) == ONIOMLayer.LOW


def test_one_link_per_boundary():
    """Exactly one boundary atom per cut (no 2-caps-on-one-Cα)."""
    s, m = _run_boundary(5, [2])
    assert len(m.boundary_atom_indices) == 2


def test_lone_gap_residue_promoted_whole():
    """A LOW residue flanked by HIGH on both sides is promoted whole."""
    s, m = _run_boundary(5, [1, 3])
    # res2 (sandwiched) becomes entirely HIGH.
    for a in s.residues[2].atoms:
        assert m.layer_membership[a.idx] == ONIOMLayer.HIGH
    # Boundaries move outward to res0.CA and res4.CA.
    boundary_names = {(s.atoms[i].residue.idx, s.atoms[i].name)
                      for i in m.boundary_atom_indices}
    assert boundary_names == {(0, "CA"), (4, "CA")}


# ── HX cap parameter tests ───────────────────────────────────────────────
_Bond = namedtuple("Bond", "force_constant eq_length")
_Angle = namedtuple("Angle", "force_constant eq_angle")
_VDW = namedtuple("VDW", "radius well_depth")


class _FFReader:
    """ff14SB-like: has N–H but no carbonyl C–H (and no C-based angles)."""
    def get_bond_parameter(self, t1, t2):
        return _Bond(434.0, 1.0100) if {t1, t2} == {"N", "H"} else None

    def get_angle_parameter(self, a, b, c):
        return _Angle(50.0, 118.0) if b == "N" else None

    def get_nonbonded_parameter(self, t):
        return _VDW(0.6, 0.0157) if t in ("H", "H1", "HC", "HA") else None


def _writer_with_caps(parent_types):
    _Atom = namedtuple("Atom", "type")
    _Link = namedtuple("Link", "qm_parent_idx mm_parent_idx")

    class _Parm:
        atoms = [_Atom(t) for t in parent_types]

    class _Setup:
        link_atoms = [_Link(i, 90 + i) for i in range(len(parent_types))]
        parm = _Parm()

    cls = next(getattr(W, n) for n in dir(W)
               if isinstance(getattr(W, n), type)
               and hasattr(getattr(W, n), "_write_link_atom_parameters"))
    w = cls.__new__(cls)
    w.oniom_setup = _Setup()
    w.ff_reader = _FFReader()
    w.type_remap = {}
    w.logger = logging.getLogger("test")
    return w


def test_formamide_cap_bond_supplied_from_gaff():
    """Carbonyl C parent → C–HX bond from GAFF c-h4 (ff14SB lacks it)."""
    w = _writer_with_caps(["C"])
    buf = io.StringIO()
    w._write_link_atom_parameters(buf, {"HX", "C", "O", "N"}, include_comments=False)
    out = buf.getvalue()
    assert "HrmStr1  C     HX" in out
    assert "310.7" in out and "1.1121" in out          # GAFF c-h4


def test_formamide_cap_angles_supplied_from_gaff():
    """Carbonyl C parent → O–C–HX and N–C–HX angles from GAFF."""
    w = _writer_with_caps(["C"])
    buf = io.StringIO()
    w._write_link_atom_parameters(buf, {"HX", "C", "O", "N"}, include_comments=False)
    out = buf.getvalue()
    assert "54.2000" in out and "120.7000" in out       # GAFF h4-c-o
    assert "50.7000" in out and "113.4400" in out       # GAFF h4-c2-n


def test_amide_cap_uses_ff14sb_n_h():
    """Amide N parent → N–HX bond copied from the FF's N–H (not GAFF)."""
    w = _writer_with_caps(["N"])
    buf = io.StringIO()
    w._write_link_atom_parameters(buf, {"HX", "N", "C"}, include_comments=False)
    out = buf.getvalue()
    assert "HrmStr1  N     HX" in out
    assert "434.0000" in out and "1.0100" in out
