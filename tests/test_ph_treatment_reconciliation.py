"""Regression: when a structure's heme propionates have been baked from the
titratable PRN into static PRP/PRD (the PB-Titrate modern-FF rebuild), the
Topology Generator must switch the heme site from its recorded constant-pH set
to the fixed-pH companion — otherwise it loads the constant-pH lib (which
defines PRN, not PRP/PRD) under leaprc.constph/conste and tleap fails.

Two layers are covered:
  1. forcefield_params.find_companion_set — resolve the opposite-treatment twin
     of a forcefield set by shared base name (real bundle metadata).
  2. TLeapInputGenerator._reconcile_ph_treatment_with_structure — flip stale
     constant_pH heme entries to fixed_pH when PRP/PRD are present (logic tested
     against a lightweight stub so it needs no real workspace/PDB).
"""
import pytest

from proprep.forcefield_params import find_companion_set
from proprep.forcefield_params.loader import _set_base_name
from proprep.tleap_prep.tleap_input_generator import TLeapInputGenerator


# --- Layer 1: companion-set resolution (real bundle) -----------------------

_CP = "heme/bis_his_c_type"


def test_base_name_strips_ph_treatment_suffix():
    assert _set_base_name("Henriques_HCO_RESP_constpH") == "Henriques_HCO_RESP"
    assert _set_base_name("Henriques_HCO_RESP_FixedpH") == "Henriques_HCO_RESP"
    # no suffix → unchanged
    assert _set_base_name("Foo_Bar") == "Foo_Bar"


def test_constph_resolves_to_fixedph_companion():
    assert find_companion_set(_CP, "oxidized", "low_spin",
                              "Henriques_HCO_RESP_constpH",
                              "fixed_pH") == "Henriques_HCO_RESP_FixedpH"


def test_fixedph_resolves_back_to_constph_companion():
    assert find_companion_set(_CP, "oxidized", "low_spin",
                              "Henriques_HCO_RESP_FixedpH",
                              "constant_pH") == "Henriques_HCO_RESP_constpH"


def test_unknown_set_falls_back_to_default_of_target_treatment():
    # No base-name match → the default fixed-pH set (FixedpH is is_default here).
    assert find_companion_set(_CP, "oxidized", "low_spin",
                              "Nonexistent_Set", "fixed_pH") \
        == "Henriques_HCO_RESP_FixedpH"


def test_none_set_falls_back_to_default():
    assert find_companion_set(_CP, "oxidized", "low_spin", None, "fixed_pH") \
        == "Henriques_HCO_RESP_FixedpH"


def test_missing_cofactor_returns_none():
    assert find_companion_set("heme/does_not_exist", "oxidized", "low_spin",
                              "X", "fixed_pH") is None


# --- Layer 2: generator reconciliation logic (stubbed self) ----------------

class _StubConsole:
    def print(self, *a, **k):
        pass


class _StubGen:
    """Minimal stand-in exposing only what the reconciliation method touches."""
    def __init__(self, has_fixed_ph):
        self._has_fixed_ph = has_fixed_ph
        self.persisted = None

        class _P:
            console = _StubConsole()
        self.processor = _P()

    def get_from_workspace(self, key, default=None):
        return default

    def update_workspace(self, key, value):
        self.persisted = value

    def _build_uses_fixed_ph_propionate(self):
        return self._has_fixed_ph


def _heme_site(ph_treatment="constant_pH",
               set_name="Henriques_HCO_RESP_constpH"):
    return {
        "site_id": "site_1",
        "residue_name": "HEC",
        "has_transformer": True,
        "cofactor_path": _CP,
        "redox_state": "oxidized",
        "spin_state": "low_spin",
        "ph_treatment": ph_treatment,
        "forcefield_set": set_name,
    }


# Bind the real method onto the stub so we test the actual logic.
_reconcile = TLeapInputGenerator._reconcile_ph_treatment_with_structure


def test_reconcile_flips_constph_to_fixedph_when_prp_prd_present():
    info = [_heme_site()]
    out = _reconcile(_StubGen(has_fixed_ph=True), info)
    site = out[0]
    assert site["ph_treatment"] == "fixed_pH"
    assert site["forcefield_set"] == "Henriques_HCO_RESP_FixedpH"


def test_reconcile_persists_to_workspace_on_change():
    gen = _StubGen(has_fixed_ph=True)
    info = [_heme_site()]
    _reconcile(gen, info)
    assert gen.persisted is info  # re-persisted the corrected list


def test_reconcile_is_noop_without_fixed_ph_propionates():
    """The ordinary constant-pH rebuild must be left completely untouched."""
    gen = _StubGen(has_fixed_ph=False)
    info = [_heme_site()]
    out = _reconcile(gen, info)
    assert out[0]["ph_treatment"] == "constant_pH"
    assert out[0]["forcefield_set"] == "Henriques_HCO_RESP_constpH"
    assert gen.persisted is None  # nothing changed → nothing persisted


def test_reconcile_leaves_already_fixedph_sites_alone():
    gen = _StubGen(has_fixed_ph=True)
    info = [_heme_site(ph_treatment="fixed_pH",
                       set_name="Henriques_HCO_RESP_FixedpH")]
    _reconcile(gen, info)
    assert info[0]["forcefield_set"] == "Henriques_HCO_RESP_FixedpH"
    assert gen.persisted is None  # no constant_pH entry to flip


def test_reconcile_skips_non_transformer_entries():
    gen = _StubGen(has_fixed_ph=True)
    info = [{"site_id": "s", "has_transformer": False,
             "ph_treatment": "constant_pH"}]
    _reconcile(gen, info)
    assert gen.persisted is None
