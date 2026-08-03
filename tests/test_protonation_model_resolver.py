#!/usr/bin/env python3
"""
Phase 1 acceptance tests for the constant-pH / fixed-pH cofactor FF infrastructure.

Covers the forcefield_params data layer in isolation (metadata is injected via
monkeypatch, so these don't depend on any on-disk metadata.json patch):

- resolve_residue_names: role -> residue name across center / ligand /
  protonation-site roles, for fixed-pH (variant pick + default) and constant-pH
  (titratable) sets, plus the error paths.
- get_protonation_model: returns the set's block or None.
- get_prerequisite_leaprc_groups(set_name=...): per-set override with fallback
  to the cofactor-global block.
- _parse_leaprc_groups: groups + legacy flat list.

Run with: pytest tests/test_protonation_model_resolver.py
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from proprep.forcefield_params import loader  # noqa: E402
from proprep.forcefield_params.loader import (  # noqa: E402
    resolve_residue_names,
    get_protonation_model,
    get_prerequisite_leaprc_groups,
    _parse_leaprc_groups,
    InvalidForcefieldError,
)

COFACTOR = "heme/cys_axial_b_type"


def _metadata():
    """Crafted metadata mirroring the planned cys-axial schema:
    a constant-pH default set (PRN) and a fixed-pH sibling set (PRD/PRP),
    each with its own prerequisites; plus a cofactor-global fallback.
    """
    fixed_sites = [
        {"role": "propionate_a", "label": "A-ring propionate",
         "variants": {"deprotonated": "PRD", "protonated": "PRP"},
         "default": "deprotonated"},
        {"role": "propionate_d", "label": "D-ring propionate",
         "variants": {"deprotonated": "PRD", "protonated": "PRP"},
         "default": "deprotonated"},
    ]
    constant_sites = [
        {"role": "propionate_a", "label": "A-ring propionate", "titratable_residue": "PRN"},
        {"role": "propionate_d", "label": "D-ring propionate", "titratable_residue": "PRN"},
    ]
    return {
        "cofactor_type": "cys_axial_b_type",
        "prerequisites": {"leaprc_groups": [{"satisfied_by": ["leaprc.constph"]}]},
        "redox_states": {
            "oxidized": {
                "spin_states": {
                    "high_spin": {
                        "residue_name": "HTO",
                        "ligand_residue_names": {"proximal_cys": "CTO"},
                        "forcefield_sets": {
                            "Guberman_HTO_RESP": {
                                "files": {"frcmod": "x.frcmod", "lib": "x.lib"},
                                "ph_treatment": "constant_pH",
                                "protonation_model": {"sites": constant_sites},
                                "prerequisites": {
                                    "leaprc_groups": [{"satisfied_by": ["leaprc.constph"]}]
                                },
                            },
                            "Guberman_HTO_RESP_FixedpH": {
                                "files": {"frcmod": "y.frcmod", "lib": "y.lib"},
                                "ph_treatment": "fixed_pH",
                                "protonation_model": {"sites": fixed_sites},
                                "prerequisites": {
                                    "leaprc_groups": [{"satisfied_by": [
                                        "leaprc.protein.ff14SB", "leaprc.protein.ff19SB"]}]
                                },
                            },
                        },
                    }
                }
            }
        },
    }


@pytest.fixture(autouse=True)
def _patch_metadata(monkeypatch):
    monkeypatch.setattr(loader, "load_forcefield_metadata", lambda cofactor_path: _metadata())


# ---- resolve_residue_names -------------------------------------------------

def test_fixed_ph_mixed_protonation():
    out = resolve_residue_names(
        COFACTOR, "oxidized", "high_spin", "Guberman_HTO_RESP_FixedpH",
        {"propionate_a": "protonated", "propionate_d": "deprotonated"},
    )
    assert out == {"center": "HTO", "proximal_cys": "CTO",
                   "propionate_a": "PRP", "propionate_d": "PRD"}


def test_fixed_ph_defaults_when_no_choice():
    out = resolve_residue_names(
        COFACTOR, "oxidized", "high_spin", "Guberman_HTO_RESP_FixedpH")
    assert out["propionate_a"] == "PRD" and out["propionate_d"] == "PRD"


def test_constant_ph_uses_titratable_residue():
    out = resolve_residue_names(
        COFACTOR, "oxidized", "high_spin", "Guberman_HTO_RESP")
    assert out == {"center": "HTO", "proximal_cys": "CTO",
                   "propionate_a": "PRN", "propionate_d": "PRN"}


def test_invalid_protonation_choice_raises():
    with pytest.raises(InvalidForcefieldError):
        resolve_residue_names(
            COFACTOR, "oxidized", "high_spin", "Guberman_HTO_RESP_FixedpH",
            {"propionate_a": "bogus"})


def test_unknown_set_raises():
    with pytest.raises(InvalidForcefieldError):
        resolve_residue_names(COFACTOR, "oxidized", "high_spin", "NoSuchSet")


# ---- get_protonation_model -------------------------------------------------

def test_get_protonation_model_present_and_absent():
    model = get_protonation_model(
        COFACTOR, "oxidized", "high_spin", "Guberman_HTO_RESP_FixedpH")
    assert [s["role"] for s in model["sites"]] == ["propionate_a", "propionate_d"]
    assert get_protonation_model(
        COFACTOR, "oxidized", "high_spin", "NoSuchSet") is None


# ---- per-set prerequisites -------------------------------------------------

def test_per_set_prereqs_fixed_ph_avoids_constph():
    groups = get_prerequisite_leaprc_groups(COFACTOR, "Guberman_HTO_RESP_FixedpH")
    assert groups == [["leaprc.protein.ff14SB", "leaprc.protein.ff19SB"]]
    assert not any("leaprc.constph" in g for g in groups)


def test_per_set_prereqs_constant_ph_requires_constph():
    groups = get_prerequisite_leaprc_groups(COFACTOR, "Guberman_HTO_RESP")
    assert groups == [["leaprc.constph"]]


def test_global_prereqs_when_no_set_name():
    assert get_prerequisite_leaprc_groups(COFACTOR) == [["leaprc.constph"]]


def test_unknown_set_falls_back_to_global_prereqs():
    # A set name with no own prerequisites block -> cofactor-global fallback.
    assert get_prerequisite_leaprc_groups(COFACTOR, "NoSuchSet") == [["leaprc.constph"]]


# ---- _parse_leaprc_groups --------------------------------------------------

def test_parse_groups_and_legacy_flat_list():
    assert _parse_leaprc_groups(
        {"leaprc_groups": [{"satisfied_by": ["a", "b"]}, ["c"]]}) == [["a", "b"], ["c"]]
    # legacy flat list -> each its own required (size-1) group
    assert _parse_leaprc_groups({"leaprcs": ["a", "b"]}) == [["a"], ["b"]]
    assert _parse_leaprc_groups({}) == []


# ---- backward-compat: cofactor with no protonation_model -------------------

def test_no_protonation_model_resolves_center_and_ligands_only(monkeypatch):
    meta = {
        "redox_states": {"oxidized": {"spin_states": {"low_spin": {
            "residue_name": "HCO",
            "ligand_residue_names": {"axial_his_stub": "HIO"},
            "forcefield_sets": {"Henriques_HCO_RESP": {
                "files": {"frcmod": "a.frcmod", "lib": "a.lib"}}}}}}},
        "prerequisites": {"leaprc_groups": [{"satisfied_by": ["leaprc.conste"]}]},
    }
    monkeypatch.setattr(loader, "load_forcefield_metadata", lambda cofactor_path: meta)
    out = resolve_residue_names("heme/bis_his_c_type", "oxidized", "low_spin",
                                "Henriques_HCO_RESP")
    assert out == {"center": "HCO", "axial_his_stub": "HIO"}
    # no set_name still returns global; unknown set falls back to global too
    assert get_prerequisite_leaprc_groups("heme/bis_his_c_type",
                                          "Henriques_HCO_RESP") == [["leaprc.conste"]]
