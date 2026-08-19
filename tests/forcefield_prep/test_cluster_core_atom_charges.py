"""
A withheld cluster's core atoms get a formal charge, so the QM total is complete.

A pure inorganic cluster is withheld from the force-field pass as a whole
residue, so nothing in it arrives with a charge. _collect_metal_charges asked
only for the metals (``is_center``) and left the rest at None, which
PDBWriter._get_charge_from_coords counts as 0.0. The suggested QM charge was
therefore short by the core atoms' formal charge.

The numbers below are from a real 4UHX run (Fe2S2 + Mo cofactor):

    site_1 FES1309   FE1 +3, FE2 +3, S1 None, S2 None   -> suggested +2
                     with four CYM at -1: 6 - 4 + 0 = +2
                     true [Fe2S2(SCys)4]2- is -2, short by the two S2-

    site_2 MOS1310   MO -2 (a compensating value: Mo is never -2), S/O1/O2 None

The Mo row is the second-order damage: the metal's formal charge is also the
van der Waals radius key, so compensating there asks the database for a radius
that does not exist.
"""

import io

import pytest
from rich.console import Console

from proprep.forcefield_prep import metal_site_parameterizer as msp_mod
from proprep.forcefield_prep.metal_site_parameterizer import MetalSiteWorkflowManager


def _manager(answers):
    """A bare manager whose prompts come from `answers` (a list of strings)."""
    mgr = MetalSiteWorkflowManager.__new__(MetalSiteWorkflowManager)
    mgr.console = Console(file=io.StringIO(), width=200)
    mgr.processor = None
    mgr.logger = None
    mgr._asked = []
    return mgr


def _install_prompt(monkeypatch, mgr, answers):
    queue = list(answers)

    def fake_prompt(processor, prompt, **kwargs):
        mgr._asked.append(prompt)
        if not queue:
            raise AssertionError(f"unexpected extra prompt: {prompt}")
        return queue.pop(0)

    monkeypatch.setattr(msp_mod, "prompt_with_context", fake_prompt)
    return queue


def _fes_assignments():
    """The Fe2S2 site as preprocessing hands it over: no charges anywhere."""
    return {
        (1.0, 0.0, 0.0): {"resname": "FES", "resid": 1309, "chain": "A",
                          "atom_name": "FE1", "element": "FE",
                          "charge": None, "is_center": True},
        (2.0, 0.0, 0.0): {"resname": "FES", "resid": 1309, "chain": "A",
                          "atom_name": "FE2", "element": "FE",
                          "charge": None, "is_center": True},
        (1.5, 1.0, 0.0): {"resname": "FES", "resid": 1309, "chain": "A",
                          "atom_name": "S1", "element": "S",
                          "charge": None, "is_center": False},
        (1.5, -1.0, 0.0): {"resname": "FES", "resid": 1309, "chain": "A",
                           "atom_name": "S2", "element": "S",
                           "charge": None, "is_center": False},
        # Four cysteinates, already charged from the prmtop.
        (5.0, 0.0, 0.0): {"resname": "CYM", "resid": 112, "chain": "A",
                          "atom_name": "SG", "element": "S",
                          "charge": -1.0, "is_center": False},
        (6.0, 0.0, 0.0): {"resname": "CYM", "resid": 115, "chain": "A",
                          "atom_name": "SG", "element": "S",
                          "charge": -1.0, "is_center": False},
        (7.0, 0.0, 0.0): {"resname": "CYM", "resid": 147, "chain": "A",
                          "atom_name": "SG", "element": "S",
                          "charge": -1.0, "is_center": False},
        (8.0, 0.0, 0.0): {"resname": "CYM", "resid": 149, "chain": "A",
                          "atom_name": "SG", "element": "S",
                          "charge": -1.0, "is_center": False},
    }


def _total(assignments):
    """What PDBWriter sums: None counts as 0.0."""
    return sum((a["charge"] or 0.0) for a in assignments.values())


def test_bridging_sulfides_are_asked_for_and_complete_the_total(monkeypatch):
    mgr = _manager(None)
    # Two metals (charge, spin each), then the two bridging sulfides.
    _install_prompt(monkeypatch, mgr, ["+3", "5", "+3", "5", "-2", "-2"])

    assignments = _fes_assignments()
    result = mgr._collect_metal_charges(assignments, redox_site=None, interactive=True)

    assert result[(1.5, 1.0, 0.0)]["charge"] == -2.0
    assert result[(1.5, -1.0, 0.0)]["charge"] == -2.0

    # 2 Fe(+3) + 2 S(-2) + 4 CYM(-1) = -2, the real [Fe2S2(SCys)4]2-.
    assert _total(result) == pytest.approx(-2.0)


def test_pre_fix_total_was_short_by_the_sulfides(monkeypatch):
    """Pin the arithmetic the bug produced, so the regression is unambiguous."""
    assignments = _fes_assignments()
    # Metals only, as the old code did.
    assignments[(1.0, 0.0, 0.0)]["charge"] = 3.0
    assignments[(2.0, 0.0, 0.0)]["charge"] = 3.0
    assert _total(assignments) == pytest.approx(2.0)   # what the run reported
    # The gap is exactly the two sulfides.
    assert _total(assignments) - (-2.0) == pytest.approx(4.0)


def test_no_charge_is_left_as_none(monkeypatch):
    mgr = _manager(None)
    _install_prompt(monkeypatch, mgr, ["+3", "5", "+3", "5", "-2", "-2"])

    result = mgr._collect_metal_charges(_fes_assignments(), redox_site=None,
                                        interactive=True)

    unset = [(a["resname"], a["atom_name"])
             for a in result.values() if a["charge"] is None]
    assert unset == [], f"atoms left with no charge: {unset}"


def test_core_atoms_are_asked_by_name(monkeypatch):
    """The prompt has to identify WHICH core atom, like the metal prompt does."""
    mgr = _manager(None)
    _install_prompt(monkeypatch, mgr, ["+3", "5", "+3", "5", "-2", "-2"])

    mgr._collect_metal_charges(_fes_assignments(), redox_site=None, interactive=True)

    core_prompts = [p for p in mgr._asked if "Formal charge for" in p]
    assert len(core_prompts) == 2
    assert any("S1" in p for p in core_prompts)
    assert any("S2" in p for p in core_prompts)
    assert all("FES" in p and "A:1309" in p for p in core_prompts)


def test_site_with_no_cluster_core_asks_nothing_extra(monkeypatch):
    """A lone Zn ion has no core atoms — the flow must be unchanged."""
    mgr = _manager(None)
    _install_prompt(monkeypatch, mgr, ["+2", "0"])

    assignments = {
        (1.0, 0.0, 0.0): {"resname": "ZN", "resid": 300, "chain": "A",
                          "atom_name": "ZN", "element": "ZN",
                          "charge": None, "is_center": True},
        (5.0, 0.0, 0.0): {"resname": "CYM", "resid": 10, "chain": "A",
                          "atom_name": "SG", "element": "S",
                          "charge": -1.0, "is_center": False},
    }
    result = mgr._collect_metal_charges(assignments, redox_site=None, interactive=True)

    assert not [p for p in mgr._asked if "Formal charge for" in p]
    assert _total(result) == pytest.approx(1.0)   # Zn(+2) + one CYM(-1)


def test_non_interactive_zeroes_core_atoms_and_says_so(monkeypatch):
    """Automation keeps working, but the shortfall is stated, not hidden."""
    mgr = _manager(None)
    _install_prompt(monkeypatch, mgr, [])

    result = mgr._collect_metal_charges(_fes_assignments(), redox_site=None,
                                        interactive=False)

    assert result[(1.5, 1.0, 0.0)]["charge"] == 0.0
    output = mgr.console.file.getvalue()
    assert "Non-interactive" in output


def test_nothing_to_collect_returns_untouched(monkeypatch):
    mgr = _manager(None)
    _install_prompt(monkeypatch, mgr, [])

    assignments = {
        (5.0, 0.0, 0.0): {"resname": "CYM", "resid": 10, "chain": "A",
                          "atom_name": "SG", "element": "S",
                          "charge": -1.0, "is_center": False},
    }
    result = mgr._collect_metal_charges(assignments, redox_site=None, interactive=True)

    assert mgr._asked == []
    assert _total(result) == pytest.approx(-1.0)
