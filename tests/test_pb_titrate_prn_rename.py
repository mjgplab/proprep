"""Regression: pb_titrate's modern-FF PDB writer must rename the heme
propionate PRN to its static fixed-pH counterpart (PRP/PRD), not leave it as
the titratable constant-pH/redox PRN.

Background: PB-Titrate titrates each heme propionate as PRN (defined only in
conste.lib). The modern-FF rebuild path (pdb_rename.build_rename_map) renames
every titratable residue to a standard fixed-charge name so it can be rebuilt
under ff14SB/ff19SB *without* leaprc.constph/conste. PRN was missing from the
rename table, so it was silently skipped (``if rn not in _RENAME: continue``)
and written out still named PRN — which then has no template under the modern
FF (conste.lib is not loaded) and, even if it built, would stay titratable,
defeating the purpose of baking states in. The fix maps
PRN → PRP (protonated, COOH) / PRD (deprotonated, COO⁻), which ProPrep's heme
``*_FixedpH`` libraries bundle, and surfaces a build note that the rebuild must
use the matching FixedpH heme set.

These tests mock get_residue so they need no AmberTools cpinutil.py.
"""
import parmed as pmd
import pytest

import proprep.pb_titrate.pdb_rename as pdb_rename
from proprep.pb_titrate.pdb_rename import build_rename_map, _RENAME, _FALLBACK_CLASS
from proprep.pb_titrate.residues import ProtonationState, TitratableResidue


# A minimal PRN with the two states cpinutil reports: STATE 0 = deprotonated
# (PRD, COO⁻, prot_count 0), STATE 1 = protonated (PRP, COOH, prot_count 1).
def _fake_prn() -> TitratableResidue:
    s0 = ProtonationState(name="STATE 0",
                          charges={"CG": 0.7, "O1": -0.85, "O2": -0.85},
                          prot_count=0, pka_corr=0.0)
    s1 = ProtonationState(name="STATE 1",
                          charges={"CG": 0.6, "O1": -0.55, "O2": -0.65, "HO": 0.45},
                          prot_count=1, pka_corr=4.5)
    return TitratableResidue(name="PRN", pka_model=4.5, states=[s0, s1])


@pytest.fixture(autouse=True)
def _mock_get_residue(monkeypatch):
    monkeypatch.setattr(pdb_rename, "get_residue",
                        lambda name: _fake_prn() if name == "PRN"
                        else pytest.fail(f"unexpected residue {name}"))


def _struct_with_prns(numbers):
    """A bare structure with one PRN residue per (0-based) number given."""
    s = pmd.Structure()
    for n in numbers:
        s.add_atom(pmd.Atom(name="CG", atomic_number=6), "PRN", n)
    return s


def test_prn_is_registered_in_rename_table():
    """The core regression: PRN must be in the rename table (was missing)."""
    assert "PRN" in _RENAME
    assert _RENAME["PRN"] == {"prot": "PRP", "deprot": "PRD"}
    # standard-pH fallback for a propionate is the deprotonated carboxylate
    assert _FALLBACK_CLASS["PRN"] == "deprot"


def test_recommended_protonated_prn_becomes_prp():
    s = _struct_with_prns([200])
    state_map = {("PRN", 201): "STATE 1"}   # resnum is 1-based (number+1)
    rename, report = build_rename_map(s, state_map)
    assert rename[("PRN", 201)] == "PRP"
    r = report[0]
    assert r["new_resname"] == "PRP" and r["source"] == "recommended"
    assert r["warning"] is None


def test_recommended_deprotonated_prn_becomes_prd():
    s = _struct_with_prns([200])
    state_map = {("PRN", 201): "STATE 0"}
    rename, report = build_rename_map(s, state_map)
    assert rename[("PRN", 201)] == "PRD"
    assert report[0]["source"] == "recommended"


def test_unrecommended_prn_falls_back_to_prd():
    """No PB recommendation → standard-pH default (deprotonated COO⁻ = PRD)."""
    s = _struct_with_prns([200])
    rename, report = build_rename_map(s, {})
    assert rename[("PRN", 201)] == "PRD"
    assert report[0]["source"] == "fallback"
    assert report[0]["warning"]  # fallback emits a warning


def test_each_propionate_renamed_independently():
    """A heme has two PRN residues (A-ring, D-ring); each gets its own state."""
    s = _struct_with_prns([200, 201])
    state_map = {("PRN", 201): "STATE 1",   # A-ring protonated
                 ("PRN", 202): "STATE 0"}   # D-ring deprotonated
    rename, _ = build_rename_map(s, state_map)
    assert rename[("PRN", 201)] == "PRP"
    assert rename[("PRN", 202)] == "PRD"


def test_prn_rename_emits_fixedph_build_note():
    """PRP/PRD live only in the heme FixedpH lib — the report must flag that."""
    s = _struct_with_prns([200])
    _, report = build_rename_map(s, {("PRN", 201): "STATE 0"})
    note = report[0].get("note")
    assert note and "FixedpH" in note
    # it's a note, not a warning (expected behaviour, not a problem)
    assert report[0]["warning"] is None


def test_metal_fixed_prn_is_deprotonated():
    """A propionate held as a metal coordinator fixes at PRD (carboxylate)."""
    s = _struct_with_prns([200])
    metal_fixed = [{"resname": "PRN", "resnum": 201, "fix_mode": "deprotonated"}]
    rename, report = build_rename_map(s, {}, metal_fixed)
    assert rename[("PRN", 201)] == "PRD"
    assert report[0]["source"] == "metal-fixed"
