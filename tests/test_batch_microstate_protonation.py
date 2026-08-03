#!/usr/bin/env python3
"""
Batch microstate generation: propionate/pH treatment is bound per redox state.

Regression for the bug where the batch path (option 3, "Generate all redox
microstates") never asked about propionate protonation and silently fell through
to the cofactor's is_default forcefield set (PRD/fixed_pH for c-type heme), while
the single path (option 2) prompts for it.

The fix adds Step 2b (`_configure_state_protonation`) which binds the choice to
each (site, redox/spin state) so a site's reduced and oxidized forms may carry
different protonation, and `_generate_single_microstate` merges the per-state
choice into that microstate's transformer parameters.

Run with: pytest tests/test_batch_microstate_protonation.py
"""

import io
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from rich.console import Console  # noqa: E402

from proprep.redoxsite_prep.transformation.redox_transformation_manager import (  # noqa: E402
    RedoxTransformationManager,
)


class _Site:
    def __init__(self, site_id):
        self.site_id = site_id


def _manager():
    """A manager with __init__ bypassed and just enough state wired up."""
    m = object.__new__(RedoxTransformationManager)
    m.console = Console(file=io.StringIO(), force_terminal=False)
    m.transformation_parameters = {}
    m.site_state_protonation = {}
    return m


def _script_choices(m):
    """Stand in for the interactive configurator: script a redox-DEPENDENT
    protonation so reduced != oxidized, proving states are independent.

    reduced/low_spin  -> fixed_pH, both rings protonated (PRP)
    oxidized/low_spin -> constant_pH (titratable PRN, no per-ring choice)
    """
    def fake_config(pname, pdef, site_ids, site_numbers, tname):
        for sid in site_ids:
            cur = m.transformation_parameters[sid]
            redox = cur["redox_state"]
            if pname == "ph_treatment":
                cur[pname] = "fixed_pH" if redox == "reduced" else "constant_pH"
            elif pname.startswith("protonation_"):
                # per-ring protomer only applies under fixed_pH (real gating)
                if cur.get("ph_treatment") == "fixed_pH":
                    cur[pname] = "protonated"
    m._configure_choice_parameter_interactive = fake_config


def _ms(**site_states):
    """Build a microstate dict: _ms(site_1=('reduced','low_spin'), ...)."""
    return {
        sid: {"redox_state": r, "spin_state": s}
        for sid, (r, s) in site_states.items()
    }


def test_protonation_bound_per_redox_state():
    m = _manager()
    m.redox_sites = [_Site("site_1"), _Site("site_2")]
    m.site_transformer_assignments = {
        "site_1": "heme_bis_his_c_type",
        "site_2": "heme_bis_his_c_type",
    }
    _script_choices(m)

    # Both sites use both reduced and oxidized across the selected microstates.
    microstates = [
        _ms(site_1=("reduced", "low_spin"), site_2=("reduced", "low_spin")),
        _ms(site_1=("oxidized", "low_spin"), site_2=("oxidized", "low_spin")),
    ]
    assert m._configure_state_protonation(microstates) is True

    for sid in ("site_1", "site_2"):
        reduced = m.site_state_protonation[sid]["reduced_low_spin"]
        oxidized = m.site_state_protonation[sid]["oxidized_low_spin"]

        # reduced -> fixed_pH with both rings chosen
        assert reduced["ph_treatment"] == "fixed_pH"
        assert reduced["protonation_propionate_a"] == "protonated"
        assert reduced["protonation_propionate_d"] == "protonated"

        # oxidized -> constant_pH, and the per-ring keys are NOT present
        # (gated off outside fixed_pH), proving the two states are independent.
        assert oxidized["ph_treatment"] == "constant_pH"
        assert "protonation_propionate_a" not in oxidized
        assert "protonation_propionate_d" not in oxidized

    # scratch store is cleared so it cannot leak into generation
    assert m.transformation_parameters == {}


def test_no_fork_transformer_is_silent_and_unconfigured():
    """A transformer with no ph_treatment fork must not be prompted and must
    leave site_state_protonation empty (falls back to its single treatment)."""
    m = _manager()
    m.redox_sites = [_Site("site_1")]
    m.site_transformer_assignments = {"site_1": "no_transformation"}

    called = {"n": 0}

    def fake_config(*a, **k):
        called["n"] += 1
    m._configure_choice_parameter_interactive = fake_config

    microstates = [_ms(site_1=("unchanged", "unchanged"))]
    assert m._configure_state_protonation(microstates) is True
    assert m.site_state_protonation == {}
    assert called["n"] == 0


def test_only_used_states_are_configured():
    """Scoping: a state a site never uses in the selected microstates must not
    be configured, and a per-state prompt must only include sites that use it."""
    m = _manager()
    m.redox_sites = [_Site("site_1"), _Site("site_2")]
    m.site_transformer_assignments = {
        "site_1": "heme_bis_his_c_type",
        "site_2": "heme_bis_his_c_type",
    }
    _script_choices(m)

    # site_1 varies (reduced + oxidized); site_2 is always reduced.
    microstates = [
        _ms(site_1=("reduced", "low_spin"), site_2=("reduced", "low_spin")),
        _ms(site_1=("oxidized", "low_spin"), site_2=("reduced", "low_spin")),
    ]
    assert m._configure_state_protonation(microstates) is True

    # site_1 configured for both states; site_2 only for reduced.
    assert set(m.site_state_protonation["site_1"]) == {"reduced_low_spin",
                                                       "oxidized_low_spin"}
    assert set(m.site_state_protonation["site_2"]) == {"reduced_low_spin"}


def test_merge_selects_correct_state_key():
    """The param-merge logic used by _generate_single_microstate must pick the
    protonation for the site's chosen state, and fall back to {} when absent."""
    m = _manager()
    m.site_state_protonation = {
        "site_1": {
            "reduced_low_spin": {"ph_treatment": "fixed_pH",
                                 "protonation_propionate_a": "protonated"},
            "oxidized_low_spin": {"ph_treatment": "constant_pH"},
        }
    }

    def merged(site_id, states):
        params = {"redox_state": states["redox_state"],
                  "spin_state": states["spin_state"]}
        key = f"{states['redox_state']}_{states['spin_state']}"
        params.update(m.site_state_protonation.get(site_id, {}).get(key, {}))
        return params

    red = merged("site_1", {"redox_state": "reduced", "spin_state": "low_spin"})
    assert red["ph_treatment"] == "fixed_pH"
    assert red["protonation_propionate_a"] == "protonated"

    ox = merged("site_1", {"redox_state": "oxidized", "spin_state": "low_spin"})
    assert ox["ph_treatment"] == "constant_pH"
    assert "protonation_propionate_a" not in ox

    # A site with no configured protonation merges nothing (transformer default).
    none = merged("site_2", {"redox_state": "reduced", "spin_state": "low_spin"})
    assert none == {"redox_state": "reduced", "spin_state": "low_spin"}


class _Center:
    def __init__(self):
        self.resname, self.chain, self.resid = "HEC", "A", 801


class _SiteC(_Site):
    def __init__(self, site_id):
        super().__init__(site_id)
        self.centers = [_Center()]


def test_batch_transformer_info_carries_ph_treatment_and_ff_set():
    """The workspace transformer_info the batch flow stores must let the Topology
    Generator's combo readers resolve the chosen treatment + implied set, so the
    FF-set picker pre-selects/filters instead of offering every set (the bug)."""
    m = _manager()
    m.redox_sites = [_SiteC("site_1")]
    m.site_transformer_assignments = {"site_1": "heme_bis_his_c_type"}
    m.site_state_protonation = {
        "site_1": {
            "reduced_low_spin": {"ph_treatment": "fixed_pH",
                                 "protonation_propionate_a": "protonated",
                                 "protonation_propionate_d": "deprotonated"},
            "oxidized_low_spin": {"ph_treatment": "fixed_pH",
                                  "protonation_propionate_a": "deprotonated",
                                  "protonation_propionate_d": "deprotonated"},
        }
    }
    microstates = [
        _ms(site_1=("reduced", "low_spin")),
        _ms(site_1=("oxidized", "low_spin")),
    ]
    info = m._build_batch_transformer_info(microstates)

    by_combo = {(e["redox_state"], e["spin_state"]): e for e in info}
    assert set(by_combo) == {("reduced", "low_spin"), ("oxidized", "low_spin")}

    red = by_combo[("reduced", "low_spin")]
    assert red["ph_treatment"] == "fixed_pH"
    assert red["forcefield_set"] == "Henriques_HCR_RESP_FixedpH"

    ox = by_combo[("oxidized", "low_spin")]
    assert ox["ph_treatment"] == "fixed_pH"
    assert ox["forcefield_set"] == "Henriques_HCO_RESP_FixedpH"


def test_batch_transformer_info_constant_ph_resolves_constph_set():
    """constant_pH choice must resolve the conste (constant-pH) set, not fixed."""
    m = _manager()
    m.redox_sites = [_SiteC("site_1")]
    m.site_transformer_assignments = {"site_1": "heme_bis_his_c_type"}
    m.site_state_protonation = {
        "site_1": {"reduced_low_spin": {"ph_treatment": "constant_pH"}}
    }
    info = m._build_batch_transformer_info([_ms(site_1=("reduced", "low_spin"))])
    assert info[0]["ph_treatment"] == "constant_pH"
    assert info[0]["forcefield_set"] == "Henriques_HCR_RESP_constpH"


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-v"]))
