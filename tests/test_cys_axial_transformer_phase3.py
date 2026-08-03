#!/usr/bin/env python3
"""
Phase 3 acceptance tests for the cys-axial b-type heme transformer.

Two halves:
1. Backward-compat against REAL on-disk metadata (no mocks): pre-Phase-6 the
   cofactor has no protonation_model, so the transformer must behave exactly as
   before — propionates → PRN, HEM → HTO/HTR, CYS → CTO/CTR, and no ph_treatment
   parameter.
2. Fixed-pH path with metadata injected (monkeypatched), proving the four
   residue codes are now resolved from metadata: propionates → PRP/PRD per ring,
   with the structural steps (selectors, atom renames) unchanged.

Run with: pytest tests/test_cys_axial_transformer_phase3.py
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import proprep.forcefield_params as ffp  # noqa: E402
from proprep.forcefield_params import loader as ff_loader  # noqa: E402
from proprep.redoxsite_prep.transformation.transformers.cys_axial_b_type_heme import (  # noqa: E402
    CysAxialBTypeHemeTransformer as T,
)

COMPONENTS = {
    "center_chain": "A", "center_id": 100,
    "prop_a_id": 101, "prop_d_id": 102,
    "proximal_cys_id": 50, "proximal_cys_chain": "A",
}


def _seq_by_id(parameters):
    seq = T.get_transformation_sequence(COMPONENTS, parameters)
    return {t["id"]: t for t in seq}


# ---- 1. Live behavior against REAL (Phase-6-patched) metadata, no mocks -----

def test_live_param_defs_include_fork():
    defs = T.get_parameter_definitions()
    assert defs["ph_treatment"]["options"] == ["constant_pH", "fixed_pH"]
    assert "protonation_propionate_a" in defs
    assert "protonation_propionate_d" in defs


@pytest.mark.parametrize("redox,heme,cys", [
    ("oxidized", "HTO", "CTO"),
    ("reduced", "HTR", "CTR"),
])
def test_constant_ph_sequence_uses_prn(redox, heme, cys):
    steps = _seq_by_id({"redox_state": redox, "spin_state": "high_spin",
                        "ph_treatment": "constant_pH"})
    assert steps["extract_propionate_a"]["action"]["change_residue_name"] == "PRN"
    assert steps["extract_propionate_d"]["action"]["change_residue_name"] == "PRN"
    assert steps["rename_heme_redox_specific"]["action"]["change_residue_name"] == heme
    assert steps["rename_proximal_cys"]["action"]["change_residue_name"] == cys
    # structural details unchanged
    assert steps["extract_propionate_a"]["action"]["rename_atoms"] == {
        "CAA": "CA", "CBA": "CB", "CGA": "CG", "O1A": "O1", "O2A": "O2"}
    assert steps["extract_propionate_a"]["action"]["change_residue_id"] == 101


@pytest.mark.parametrize("redox,heme,cys", [
    ("oxidized", "HTO", "CTO"),
    ("reduced", "HTR", "CTR"),
])
def test_fixed_ph_sequence_real_metadata(redox, heme, cys):
    steps = _seq_by_id({"redox_state": redox, "spin_state": "high_spin",
                        "ph_treatment": "fixed_pH",
                        "protonation_propionate_a": "protonated",
                        "protonation_propionate_d": "deprotonated"})
    assert steps["extract_propionate_a"]["action"]["change_residue_name"] == "PRP"
    assert steps["extract_propionate_d"]["action"]["change_residue_name"] == "PRD"
    assert steps["rename_heme_redox_specific"]["action"]["change_residue_name"] == heme
    assert steps["rename_proximal_cys"]["action"]["change_residue_name"] == cys


def test_default_no_treatment_resolves_to_constant_prn():
    # No ph_treatment given -> representative default set (constant_pH) -> PRN
    steps = _seq_by_id({"redox_state": "oxidized", "spin_state": "high_spin"})
    assert steps["extract_propionate_a"]["action"]["change_residue_name"] == "PRN"


def test_validate_requires_ph_treatment():
    ok, _ = T.validate_parameters({"redox_state": "oxidized", "spin_state": "high_spin",
                                   "ph_treatment": "constant_pH"})
    assert ok
    bad, _ = T.validate_parameters({"redox_state": "oxidized", "spin_state": "high_spin"})
    assert not bad  # ph_treatment is now a forked param -> required


# ---- 2. Fixed-pH path (injected metadata) ----------------------------------

def _metadata():
    fixed_sites = [
        {"role": "propionate_a", "label": "A-ring propionate",
         "variants": {"deprotonated": "PRD", "protonated": "PRP"}, "default": "deprotonated"},
        {"role": "propionate_d", "label": "D-ring propionate",
         "variants": {"deprotonated": "PRD", "protonated": "PRP"}, "default": "deprotonated"},
    ]
    constant_sites = [
        {"role": "propionate_a", "label": "A-ring propionate", "titratable_residue": "PRN"},
        {"role": "propionate_d", "label": "D-ring propionate", "titratable_residue": "PRN"},
    ]
    return {
        "redox_states": {"oxidized": {"spin_states": {"high_spin": {
            "residue_name": "HTO",
            "ligand_residue_names": {"proximal_cys": "CTO"},
            "forcefield_sets": {
                "Guberman_HTO_RESP": {
                    "files": {"frcmod": "x.frcmod", "lib": "x.lib"},
                    "ph_treatment": "constant_pH",
                    "protonation_model": {"sites": constant_sites}, "is_default": True},
                "Guberman_HTO_RESP_FixedpH": {
                    "files": {"frcmod": "y.frcmod", "lib": "y.lib"},
                    "ph_treatment": "fixed_pH",
                    "protonation_model": {"sites": fixed_sites}, "is_default": False},
            },
        }}}}
    }


@pytest.fixture
def fixed_ph(monkeypatch):
    meta = _metadata()

    def fake_discover(path, redox, spin):
        if (redox, spin) != ("oxidized", "high_spin"):
            return []
        return [{"name": n, "ph_treatment": i.get("ph_treatment"),
                 "protonation_model": i.get("protonation_model"),
                 "is_default": i.get("is_default", False)}
                for n, i in meta["redox_states"][redox]["spin_states"][spin]["forcefield_sets"].items()]

    monkeypatch.setattr(ff_loader, "load_forcefield_metadata", lambda p: meta)
    monkeypatch.setattr(ffp, "load_forcefield_metadata", lambda p: meta)
    monkeypatch.setattr(ffp, "discover_forcefield_files", fake_discover)


def test_param_defs_gain_fork_and_site_params(fixed_ph):
    defs = T.get_parameter_definitions()
    assert defs["ph_treatment"]["options"] == ["constant_pH", "fixed_pH"]
    assert "protonation_propionate_a" in defs and "protonation_propionate_d" in defs


def test_fixed_ph_sequence_resolves_protomers(fixed_ph):
    steps = _seq_by_id({
        "redox_state": "oxidized", "spin_state": "high_spin",
        "ph_treatment": "fixed_pH",
        "protonation_propionate_a": "protonated",
        "protonation_propionate_d": "deprotonated"})
    assert steps["extract_propionate_a"]["action"]["change_residue_name"] == "PRP"
    assert steps["extract_propionate_d"]["action"]["change_residue_name"] == "PRD"
    assert steps["rename_heme_redox_specific"]["action"]["change_residue_name"] == "HTO"
    assert steps["rename_proximal_cys"]["action"]["change_residue_name"] == "CTO"
    # structural details still unchanged under fixed-pH
    assert steps["extract_propionate_a"]["action"]["rename_atoms"] == {
        "CAA": "CA", "CBA": "CB", "CGA": "CG", "O1A": "O1", "O2A": "O2"}


def test_constant_ph_sequence_uses_prn(fixed_ph):
    steps = _seq_by_id({
        "redox_state": "oxidized", "spin_state": "high_spin",
        "ph_treatment": "constant_pH"})
    assert steps["extract_propionate_a"]["action"]["change_residue_name"] == "PRN"
    assert steps["extract_propionate_d"]["action"]["change_residue_name"] == "PRN"


def test_validate_fixed_ph(fixed_ph):
    base = {"redox_state": "oxidized", "spin_state": "high_spin", "ph_treatment": "fixed_pH"}
    ok, _ = T.validate_parameters({**base, "protonation_propionate_a": "protonated",
                                   "protonation_propionate_d": "deprotonated"})
    assert ok
    bad, msg = T.validate_parameters({**base, "protonation_propionate_a": "bogus",
                                      "protonation_propionate_d": "deprotonated"})
    assert not bad
    # constant_pH needs no protonation params (they're gated off)
    ok2, _ = T.validate_parameters({"redox_state": "oxidized", "spin_state": "high_spin",
                                    "ph_treatment": "constant_pH"})
    assert ok2
