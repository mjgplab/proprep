"""Regression: PB-Titrate's pre-PB strip must keep STRUCTURAL ions and remove
BULK neutralizing ions — even when they share an element and residue name.

Amber writes calcium as residue ``CA`` for both a structural site and a bulk
counter-ion, so the old name-based strip (``:Ca2+`` / ``:CA``) could not tell
them apart: it either kept both (bulk Ca²⁺ leaked into the protonation-renamed
PDB) or, if the name were added, stripped the structural ones too. The fix
classifies each monatomic ion by RedoxSite membership + a coordination
backstop. These tests build a tiny structure (structural Ca²⁺ coordinated to an
Asp carboxylate, bulk Ca²⁺ off in solvent) and check the split.
"""
import types

import numpy as np
import parmed as pmd

import proprep.pb_titrate.workflow as wf
from proprep.pb_titrate.workflow import PBTitrateWorkflow
from proprep.pb_titrate.sites import _parmed_reskey

_classify = PBTitrateWorkflow._classify_ions_for_strip
_confirm = PBTitrateWorkflow._confirm_ion_classification


class _StubWS:
    """Carries only what _classify_ions_for_strip touches."""
    WATER_RESNAMES = PBTitrateWorkflow.WATER_RESNAMES


def _build_structure():
    """Asp carboxylate + structural Ca²⁺ (2.4 Å from OD1) + bulk Ca²⁺ (50 Å) +
    a water. Atom add order must match the coordinate rows below."""
    s = pmd.Structure()
    s.add_atom(pmd.Atom(name="N", atomic_number=7), "ASP", 1, chain="A")
    s.add_atom(pmd.Atom(name="CA", atomic_number=6), "ASP", 1, chain="A")   # alpha C
    s.add_atom(pmd.Atom(name="CG", atomic_number=6), "ASP", 1, chain="A")
    s.add_atom(pmd.Atom(name="OD1", atomic_number=8), "ASP", 1, chain="A")
    s.add_atom(pmd.Atom(name="CA", atomic_number=20), "CA", 100, chain="A")  # structural Ca²⁺
    s.add_atom(pmd.Atom(name="CA", atomic_number=20), "CA", 200, chain="A")  # bulk Ca²⁺
    s.add_atom(pmd.Atom(name="O", atomic_number=8), "WAT", 300, chain="A")
    # OD1 is the sole nearby coordinator (2.4 Å); the rest of Asp is held off so
    # the cutoff test can cross the 2.4 Å contact cleanly.
    s.coordinates = np.array([
        [5.0, 5.0, 0.0],     # N
        [4.0, 5.0, 0.0],     # alpha C
        [1.0, 5.0, 0.0],     # CG   (~5.2 Å from the structural Ca)
        [0.0, 0.0, 0.0],     # OD1  (carboxylate O)
        [2.4, 0.0, 0.0],     # structural Ca²⁺ — 2.4 Å from OD1
        [50.0, 0.0, 0.0],    # bulk Ca²⁺ — far, no coordinator nearby
        [10.0, 10.0, 10.0],  # water O
    ], dtype=float)
    return s


def _by_resnum(classified):
    return {e["resnum"]: e for e in classified}


def test_structural_ca_kept_via_coordination_bulk_stripped():
    s = _build_structure()
    out = _classify(_StubWS(), s, redox_keys=set(), cutoff=3.0)
    d = _by_resnum(out)
    # only the two single-atom Ca residues are ion candidates (not ASP, not WAT)
    assert set(d) == {101, 201}
    assert d[101]["decision"] == "keep"
    assert "coordinated" in d[101]["reason"]
    assert d[201]["decision"] == "strip"
    assert "bulk" in d[201]["reason"]


def test_redox_site_membership_overrides_to_keep():
    """A bulk-looking ion that the user declared as a RedoxSite is kept."""
    s = _build_structure()
    bulk = next(r for r in s.residues if r.name == "CA" and r.number == 200)
    redox_keys = {_parmed_reskey(bulk)}
    out = _classify(_StubWS(), s, redox_keys=redox_keys, cutoff=3.0)
    d = _by_resnum(out)
    assert d[201]["decision"] == "keep"
    assert d[201]["reason"] == "RedoxSite"


def test_tight_cutoff_drops_coordination_backstop():
    """With a cutoff below the 2.4 Å contact, the structural Ca²⁺ is no longer
    'coordinated' — and with no RedoxSite it falls through to bulk."""
    s = _build_structure()
    out = _classify(_StubWS(), s, redox_keys=set(), cutoff=2.0)
    d = _by_resnum(out)
    assert d[101]["decision"] == "strip"
    assert d[201]["decision"] == "strip"


def test_no_ions_returns_empty():
    s = pmd.Structure()
    s.add_atom(pmd.Atom(name="N", atomic_number=7), "ALA", 1, chain="A")
    s.add_atom(pmd.Atom(name="CA", atomic_number=6), "ALA", 1, chain="A")
    s.coordinates = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
    assert _classify(_StubWS(), s, redox_keys=set(), cutoff=3.0) == []


# --- interactive override -------------------------------------------------

class _Console:
    def print(self, *a, **k):
        pass


def _stub_self():
    return types.SimpleNamespace(console=_Console(), processor=None)


def test_confirm_toggle_reclassifies_by_resnum(monkeypatch):
    classified = [
        {"resname": "CA", "resnum": 100, "chain": "A", "atom_idx": 4,
         "decision": "keep", "reason": "coordinated"},
        {"resname": "CA", "resnum": 200, "chain": "A", "atom_idx": 5,
         "decision": "strip", "reason": "bulk (not coordinated)"},
    ]
    monkeypatch.setattr(wf, "prompt_with_context", lambda **k: "200")
    _confirm(_stub_self(), classified)
    d = _by_resnum(classified)
    assert d[200]["decision"] == "keep"   # toggled strip → keep
    assert d[100]["decision"] == "keep"   # untouched


def test_confirm_empty_input_accepts_as_is(monkeypatch):
    classified = [
        {"resname": "CA", "resnum": 200, "chain": "A", "atom_idx": 5,
         "decision": "strip", "reason": "bulk (not coordinated)"},
    ]
    monkeypatch.setattr(wf, "prompt_with_context", lambda **k: "")
    _confirm(_stub_self(), classified)
    assert classified[0]["decision"] == "strip"


def test_confirm_ignores_unknown_resnum(monkeypatch):
    classified = [
        {"resname": "CA", "resnum": 100, "chain": "A", "atom_idx": 4,
         "decision": "keep", "reason": "coordinated"},
    ]
    monkeypatch.setattr(wf, "prompt_with_context", lambda **k: "999 foo")
    _confirm(_stub_self(), classified)
    assert classified[0]["decision"] == "keep"  # unchanged
