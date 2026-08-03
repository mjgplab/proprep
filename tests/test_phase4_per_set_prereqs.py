#!/usr/bin/env python3
"""
Phase 4 acceptance tests: the Topology Generator's per-set leaprc prerequisites.

The functional core of Phase 4 is that _collect_cofactor_prereq_groups reads each
site's Stage-1-implied forcefield set (from transformer_info) and applies that
SET's prerequisites — so a fixed-pH heme requires ff14SB/ff19SB and does NOT drag
in leaprc.constph, while a constant-pH heme still requires constph.

Metadata is injected (monkeypatch) so the fixed-pH set + per-set prereqs exist
without the Phase-6 on-disk patch. The generator instance is built with
object.__new__ (these methods only touch get_from_workspace + the resolver).

Run with: pytest tests/test_phase4_per_set_prereqs.py
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from proprep.forcefield_params import loader as ff_loader  # noqa: E402
from proprep.tleap_prep.tleap_input_generator import TLeapInputGenerator  # noqa: E402

COFACTOR = "heme/cys_axial_b_type"


def _metadata():
    return {
        "prerequisites": {"leaprc_groups": [{"satisfied_by": ["leaprc.constph"]}]},
        "redox_states": {"oxidized": {"spin_states": {"high_spin": {
            "residue_name": "HTO",
            "ligand_residue_names": {"proximal_cys": "CTO"},
            "forcefield_sets": {
                "Guberman_HTO_RESP": {
                    "files": {"frcmod": "x.frcmod", "lib": "x.lib"},
                    "ph_treatment": "constant_pH",
                    "prerequisites": {"leaprc_groups": [{"satisfied_by": ["leaprc.constph"]}]},
                },
                "Guberman_HTO_RESP_FixedpH": {
                    "files": {"frcmod": "y.frcmod", "lib": "y.lib"},
                    "ph_treatment": "fixed_pH",
                    "prerequisites": {"leaprc_groups": [{"satisfied_by": [
                        "leaprc.protein.ff14SB", "leaprc.protein.ff19SB"]}]},
                },
            },
        }}}},
    }


@pytest.fixture(autouse=True)
def _patch_metadata(monkeypatch):
    monkeypatch.setattr(ff_loader, "load_forcefield_metadata", lambda p: _metadata())


def _gen(transformer_info):
    gen = object.__new__(TLeapInputGenerator)
    gen.get_from_workspace = lambda key, default=None: (
        transformer_info if key == "transformer_info" else default)
    return gen


def _site(set_name, ph_treatment, resname="HTO"):
    return {
        "has_transformer": True,
        "transformer_type": "heme_cys_axial_b_type",
        "cofactor_path": COFACTOR,
        "residue_name": resname,
        "redox_state": "oxidized",
        "spin_state": "high_spin",
        "forcefield_set": set_name,
        "ph_treatment": ph_treatment,
        "parameters": {"redox_state": "oxidized", "spin_state": "high_spin",
                       "ph_treatment": ph_treatment},
    }


# ---- _resolve_site_forcefield_set -----------------------------------------

def test_resolve_uses_recorded_set():
    gen = _gen([])
    assert gen._resolve_site_forcefield_set(
        _site("Guberman_HTO_RESP_FixedpH", "fixed_pH")) == "Guberman_HTO_RESP_FixedpH"


def test_resolve_none_on_unknown_transformer():
    gen = _gen([])
    # no recorded set + unknown transformer type -> None (caller falls back to global)
    assert gen._resolve_site_forcefield_set(
        {"transformer_type": "does_not_exist", "parameters": {}}) is None


# ---- _collect_cofactor_prereq_groups (the functional core) -----------------

def test_fixed_ph_site_requires_ff_not_constph():
    gen = _gen([_site("Guberman_HTO_RESP_FixedpH", "fixed_pH")])
    cof = gen._collect_cofactor_prereq_groups()
    assert len(cof) == 1
    groups = cof[0]["groups"]
    assert groups == [["leaprc.protein.ff14SB", "leaprc.protein.ff19SB"]]
    assert not any("leaprc.constph" in g for g in groups)


def test_constant_ph_site_requires_constph():
    gen = _gen([_site("Guberman_HTO_RESP", "constant_pH")])
    cof = gen._collect_cofactor_prereq_groups()
    assert cof[0]["groups"] == [["leaprc.constph"]]


def test_no_recorded_set_falls_back_to_global():
    # An older transformer_info with no forcefield_set + an unknown transformer
    # type -> set unresolved -> cofactor-global prereqs (constph here).
    site = _site(None, None)
    site["transformer_type"] = "does_not_exist"
    gen = _gen([site])
    cof = gen._collect_cofactor_prereq_groups()
    assert cof[0]["groups"] == [["leaprc.constph"]]


def test_site_without_cofactor_path_skipped():
    gen = _gen([{"has_transformer": True, "transformer_type": "x"}])  # no cofactor_path
    assert gen._collect_cofactor_prereq_groups() == []
