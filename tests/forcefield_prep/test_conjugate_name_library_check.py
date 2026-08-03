"""
Regression tests: naming a covalent AA<->ligand adduct consults the forcefield
library so two different conjugates cannot both land on the same tLEaP unit name.

Two behaviors are pinned:

  1. loader.get_registered_residue_names enumerates every unit name the bundled +
     user library already defines (str residue_name, {source: code} maps, and
     ligand_residue_names), and never raises — an unreadable library yields {}.

  2. ForcefieldParameterizer._name_conjugate_residue advances the DEFAULT past a
     library-occupied name (a second Cys adduct does not default to CS1 again),
     and, when the user TYPES a name a library entry already owns, warns and
     requires explicit confirmation before reusing it.
"""

from types import SimpleNamespace

import pytest

from proprep.forcefield_params import loader
import proprep.forcefield_prep.forcefield_parameterizer as FP
from proprep.forcefield_prep.forcefield_parameterizer import ForcefieldParameterizer


# ── loader.get_registered_residue_names ─────────────────────────────────────

def test_registered_names_collects_str_map_and_ligands(monkeypatch):
    fake = {
        "heme/bis_his_c_type": {"redox_states": {
            "fe3": {"spin_states": {"low": {
                # single-residue: residue_name is a bare string
                "residue_name": "HTO",
                "ligand_residue_names": {"prox": "HIP", "dist": "MET"},
            }}}}},
        "modified_aa/CS1": {"redox_states": {
            "single_state": {"spin_states": {"default": {
                # covalent adduct: residue_name is a {source: code} map
                "residue_name": {"CYS": "CS1", "MOV": "MO1"},
            }}}}},
    }
    monkeypatch.setattr(loader, "get_available_cofactor_types", lambda: fake)

    names = loader.get_registered_residue_names()

    assert names["HTO"] == "heme/bis_his_c_type"
    assert names["HIP"] == "heme/bis_his_c_type"
    assert names["MET"] == "heme/bis_his_c_type"
    assert names["CS1"] == "modified_aa/CS1"
    assert names["MO1"] == "modified_aa/CS1"


def test_registered_names_never_raises(monkeypatch):
    def boom():
        raise RuntimeError("library on fire")
    monkeypatch.setattr(loader, "get_available_cofactor_types", boom)

    assert loader.get_registered_residue_names() == {}


def test_registered_names_uppercases_and_skips_blanks(monkeypatch):
    fake = {"x/y": {"redox_states": {"s": {"spin_states": {"d": {
        "residue_name": "abc", "ligand_residue_names": {"r": "  "},
    }}}}}}
    monkeypatch.setattr(loader, "get_available_cofactor_types", lambda: fake)

    names = loader.get_registered_residue_names()
    assert "ABC" in names          # uppercased
    assert all(k.strip() for k in names)  # blank ligand name dropped


# ── _name_conjugate_residue ─────────────────────────────────────────────────

def _namer():
    p = object.__new__(ForcefieldParameterizer)
    p.console = FP._console if hasattr(FP, "_console") else _Console()
    p.processor = None
    return p


class _Console:
    def print(self, *a, **k):
        pass


def _residues():
    return [SimpleNamespace(name="CYS"), SimpleNamespace(name="MOV")]


def test_default_advances_past_library_name(monkeypatch):
    # CS1 already owned by an earlier adduct -> the default must become CS2.
    monkeypatch.setattr(loader, "get_registered_residue_names",
                        lambda: {"CS1": "modified_aa/CS1"})
    seen = {}

    def fake_prompt(processor, prompt, default=None, **k):
        seen["default"] = default
        return default  # user accepts the offered default

    monkeypatch.setattr(FP, "prompt_with_context", fake_prompt)
    p = _namer()
    p.console = _Console()

    name = p._name_conjugate_residue("CYS", 12, _residues())

    assert seen["default"] == "CS2"
    assert name == "CS2"


def test_default_is_cs1_when_library_empty(monkeypatch):
    monkeypatch.setattr(loader, "get_registered_residue_names", lambda: {})
    seen = {}
    monkeypatch.setattr(FP, "prompt_with_context",
                        lambda processor, prompt, default=None, **k:
                        seen.setdefault("default", default) or default)
    p = _namer()
    p.console = _Console()

    name = p._name_conjugate_residue("CYS", 12, _residues())
    assert seen["default"] == "CS1"
    assert name == "CS1"


def test_typed_library_collision_declined_reprompts(monkeypatch):
    monkeypatch.setattr(loader, "get_registered_residue_names",
                        lambda: {"ZZZ": "small_molecules/ZZZ"})
    replies = iter(["ZZZ", "QQQ"])  # first collides, second is free
    monkeypatch.setattr(FP, "prompt_with_context",
                        lambda *a, **k: next(replies))
    confirm_calls = []
    monkeypatch.setattr(FP, "confirm_with_context",
                        lambda *a, **k: confirm_calls.append(a) or False)
    p = _namer()
    p.console = _Console()

    name = p._name_conjugate_residue("CYS", 12, _residues())

    assert name == "QQQ"                 # advanced past the declined name
    assert len(confirm_calls) == 1       # confirmation asked once, for ZZZ


def test_typed_library_collision_confirmed_reuses(monkeypatch):
    monkeypatch.setattr(loader, "get_registered_residue_names",
                        lambda: {"ZZZ": "small_molecules/ZZZ"})
    monkeypatch.setattr(FP, "prompt_with_context", lambda *a, **k: "ZZZ")
    monkeypatch.setattr(FP, "confirm_with_context", lambda *a, **k: True)
    p = _namer()
    p.console = _Console()

    name = p._name_conjugate_residue("CYS", 12, _residues())
    assert name == "ZZZ"                 # explicit reuse honored


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-v"]))
