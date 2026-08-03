#!/usr/bin/env python3
"""
Phase 2 acceptance tests for the generic protonation / pH-treatment hooks on
RedoxSiteTransformerBase.

Exercised in isolation: forcefield_params is patched (no disk metadata needed,
no Phase-6 patch needed). A dummy transformer subclass drives the base hooks:

- protonation_parameter_definitions(): ph_treatment fork + per-site params
  (empty static options so they're gated), and {} when only one treatment.
- available_ph_treatments() and get_valid_options() gating (protonation params
  visible only under fixed_pH).
- select_forcefield_set_name() per treatment.
- resolve_output_residue_names(): role -> residue name via the resolver.

Run with: pytest tests/test_transformer_protonation_hooks.py
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import proprep.forcefield_params as ffp  # noqa: E402
from proprep.forcefield_params import loader as ff_loader  # noqa: E402
from proprep.redoxsite_prep.transformation.redox_transformer_framework import (  # noqa: E402
    RedoxSiteTransformerBase,
)


class DummyHeme(RedoxSiteTransformerBase):
    TRANSFORMER_NAME = "dummy_heme"
    FORCEFIELD_PATH = "fake/heme"


# ---- fixtures: crafted metadata + matching discover output -----------------

def _metadata(with_fixed=True):
    constant_sites = [
        {"role": "propionate_a", "label": "A-ring propionate", "titratable_residue": "PRN"},
        {"role": "propionate_d", "label": "D-ring propionate", "titratable_residue": "PRN"},
    ]
    fixed_sites = [
        {"role": "propionate_a", "label": "A-ring propionate",
         "variants": {"deprotonated": "PRD", "protonated": "PRP"}, "default": "deprotonated"},
        {"role": "propionate_d", "label": "D-ring propionate",
         "variants": {"deprotonated": "PRD", "protonated": "PRP"}, "default": "deprotonated"},
    ]
    sets = {
        "Guberman_HTO_RESP": {
            "files": {"frcmod": "x.frcmod", "lib": "x.lib"},
            "ph_treatment": "constant_pH",
            "protonation_model": {"sites": constant_sites},
            "is_default": True,
        }
    }
    if with_fixed:
        sets["Guberman_HTO_RESP_FixedpH"] = {
            "files": {"frcmod": "y.frcmod", "lib": "y.lib"},
            "ph_treatment": "fixed_pH",
            "protonation_model": {"sites": fixed_sites},
            "is_default": False,
        }
    return {
        "redox_states": {"oxidized": {"spin_states": {"high_spin": {
            "residue_name": "HTO",
            "ligand_residue_names": {"proximal_cys": "CTO"},
            "forcefield_sets": sets,
        }}}}
    }


def _install(monkeypatch, with_fixed=True):
    meta = _metadata(with_fixed)

    def fake_discover(path, redox, spin):
        if (redox, spin) != ("oxidized", "high_spin"):
            return []
        out = []
        for name, info in meta["redox_states"][redox]["spin_states"][spin]["forcefield_sets"].items():
            out.append({
                "name": name,
                "ph_treatment": info.get("ph_treatment"),
                "protonation_model": info.get("protonation_model"),
                "is_default": info.get("is_default", False),
            })
        return out

    # resolver + get_protonation_model call loader's module-global; the
    # framework's _all_redox_spin_pairs imports load_forcefield_metadata from
    # the package — patch both bindings.
    monkeypatch.setattr(ff_loader, "load_forcefield_metadata", lambda p: meta)
    monkeypatch.setattr(ffp, "load_forcefield_metadata", lambda p: meta)
    # set enumeration reads files on disk -> patch to crafted output
    monkeypatch.setattr(ffp, "discover_forcefield_files", fake_discover)


STATE = {"redox_state": "oxidized", "spin_state": "high_spin"}


# ---- parameter definitions -------------------------------------------------

def test_param_defs_include_fork_and_gated_site_params(monkeypatch):
    _install(monkeypatch)
    defs = DummyHeme.protonation_parameter_definitions()
    assert defs["ph_treatment"]["options"] == ["constant_pH", "fixed_pH"]
    assert defs["ph_treatment"]["default"] == "constant_pH"
    for role in ("protonation_propionate_a", "protonation_propionate_d"):
        assert defs[role]["type"] == "choice"
        assert defs[role]["options"] == []          # gated (empty static options)
        assert defs[role]["default"] == "deprotonated"


def test_no_fork_when_single_treatment(monkeypatch):
    _install(monkeypatch, with_fixed=False)
    assert DummyHeme.protonation_parameter_definitions() == {}


# ---- option gating ---------------------------------------------------------

def test_available_ph_treatments(monkeypatch):
    _install(monkeypatch)
    assert DummyHeme.available_ph_treatments("oxidized", "high_spin") == ["constant_pH", "fixed_pH"]


def test_ph_treatment_valid_options(monkeypatch):
    _install(monkeypatch)
    assert DummyHeme.get_valid_options("ph_treatment", STATE) == ["constant_pH", "fixed_pH"]


def test_protonation_gated_off_under_constant(monkeypatch):
    _install(monkeypatch)
    params = {**STATE, "ph_treatment": "constant_pH"}
    assert DummyHeme.get_valid_options("protonation_propionate_a", params) == []


def test_protonation_options_under_fixed(monkeypatch):
    _install(monkeypatch)
    params = {**STATE, "ph_treatment": "fixed_pH"}
    assert DummyHeme.get_valid_options("protonation_propionate_a", params) == ["deprotonated", "protonated"]


# ---- set selection + resolution -------------------------------------------

def test_select_set_name_by_treatment(monkeypatch):
    _install(monkeypatch)
    assert DummyHeme.select_forcefield_set_name(
        {**STATE, "ph_treatment": "fixed_pH"}) == "Guberman_HTO_RESP_FixedpH"
    assert DummyHeme.select_forcefield_set_name(
        {**STATE, "ph_treatment": "constant_pH"}) == "Guberman_HTO_RESP"


def test_resolve_output_residue_names_fixed(monkeypatch):
    _install(monkeypatch)
    params = {**STATE, "ph_treatment": "fixed_pH",
              "protonation_propionate_a": "protonated",
              "protonation_propionate_d": "deprotonated"}
    assert DummyHeme.resolve_output_residue_names(params) == {
        "center": "HTO", "proximal_cys": "CTO",
        "propionate_a": "PRP", "propionate_d": "PRD"}


def test_resolve_output_residue_names_constant(monkeypatch):
    _install(monkeypatch)
    params = {**STATE, "ph_treatment": "constant_pH"}
    assert DummyHeme.resolve_output_residue_names(params) == {
        "center": "HTO", "proximal_cys": "CTO",
        "propionate_a": "PRN", "propionate_d": "PRN"}
