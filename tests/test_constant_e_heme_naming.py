#!/usr/bin/env python3
"""
Acceptance tests for the constant_E redox-treatment fork on the bis-His c-type
heme.

AMBER's constant-redox machinery is name-keyed: ceinutil.py selects titratable
residues by matching RESIDUE_LABEL in the prmtop against parmed's
titratable_residues, where the c-type heme is HEH. A topology naming it HCO or
HCR is skipped in silence, so the heme never titrates. These tests pin the
naming end to end.

Run against the REAL shipped metadata (not mocks), because the thing under test
is what the library actually emits. Framework-level hook behaviour is checked
too, including that cofactors without a redox fork are unaffected.

Run with: pytest tests/test_constant_e_heme_naming.py
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from proprep.forcefield_params import (  # noqa: E402
    discover_forcefield_files,
    resolve_residue_names,
)
from proprep.forcefield_params.loader import (  # noqa: E402
    find_companion_set,
    get_registered_residue_names,
)
from proprep.redoxsite_prep.transformation.transformers.bis_his_c_type_heme import (  # noqa: E402
    BisHisCTypeHemeTransformer as BisHis,
)

COFACTOR = "heme/bis_his_c_type"

CONST_E_SETS = {
    "Henriques_HEH_RESP_constE_constpH": "constant_pH",
    "Henriques_HEH_RESP_constE_FixedpH": "fixed_pH",
}


def _sets(redox="oxidized"):
    return {s["name"]: s for s in discover_forcefield_files(COFACTOR, redox, "low_spin")}


# ---- the 2x2 matrix exists and is labelled on both axes --------------------

def test_oxidized_branch_ships_all_four_combinations():
    sets = _sets("oxidized")
    got = {(s["redox_treatment"], s["ph_treatment"]) for s in sets.values()}
    assert got == {
        ("fixed_E", "constant_pH"),
        ("fixed_E", "fixed_pH"),
        ("constant_E", "constant_pH"),
        ("constant_E", "fixed_pH"),
    }


def test_reduced_branch_is_fixed_E_only():
    """constant_E carries no committed oxidation state, so its sets live under a
    single branch. A reduced duplicate would be the same library reachable by a
    second path, which would emit a second loadoff of the same HEH unit."""
    treatments = {s["redox_treatment"] for s in _sets("reduced").values()}
    assert treatments == {"fixed_E"}


# ---- residue naming: the actual bug this fixes -----------------------------

@pytest.mark.parametrize("set_name,ph_treatment", sorted(CONST_E_SETS.items()))
def test_constant_E_sets_name_the_centre_HEH(set_name, ph_treatment):
    resolved = resolve_residue_names(COFACTOR, "oxidized", "low_spin", set_name)
    assert resolved["center"] == "HEH", (
        f"{set_name} must emit HEH so ceinutil.py recognises the heme"
    )


@pytest.mark.parametrize("redox,expected", [("oxidized", "HCO"), ("reduced", "HCR")])
def test_fixed_E_sets_keep_the_state_specific_code(redox, expected):
    for name, info in _sets(redox).items():
        if info["redox_treatment"] != "fixed_E":
            continue
        resolved = resolve_residue_names(COFACTOR, redox, "low_spin", name)
        assert resolved["center"] == expected


def test_per_set_override_does_not_leak_to_siblings():
    """residue_name sits above forcefield_sets, so the constant_E override must
    be per-set — otherwise it would rename the fixed_E siblings too."""
    centres = {
        name: resolve_residue_names(COFACTOR, "oxidized", "low_spin", name)["center"]
        for name in _sets("oxidized")
    }
    assert centres["Henriques_HCO_RESP_constpH"] == "HCO"
    assert centres["Henriques_HCO_RESP_FixedpH"] == "HCO"
    assert centres["Henriques_HEH_RESP_constE_constpH"] == "HEH"
    assert centres["Henriques_HEH_RESP_constE_FixedpH"] == "HEH"


def test_stubs_and_propionates_are_unchanged_by_the_redox_fork():
    """Only the centre is renamed: the His/Cys stubs are treatment-independent
    and the propionates stay on the pH axis."""
    for name, ph in CONST_E_SETS.items():
        r = resolve_residue_names(COFACTOR, "oxidized", "low_spin", name)
        assert r["axial_his_stub"] == "HIO"
        assert r["thioether_cys_stub"] == "CYO"
        assert r["propionate_a"] == ("PRN" if ph == "constant_pH" else "PRD")


def test_HEH_is_claimed_in_the_residue_name_registry():
    """A newly parameterized residue must not collide with HEH."""
    assert get_registered_residue_names().get("HEH") == COFACTOR


# ---- the shipped library really declares HEH -------------------------------

def _lib_units(lib_path):
    units = []
    with open(lib_path) as fh:
        for line in fh:
            if line.startswith("!entry.") and ".unit.name single str" in line:
                units.append(line.split(".")[1])
    return units


def _lib_atoms(lib_path, unit):
    """(name, type, charge) for each atom of `unit`."""
    header = (
        f"!entry.{unit}.unit.atoms table  str name  str type  int typex  "
        f"int resx  int flags  int seq  int elmnt  dbl chg"
    )
    out, inside = [], False
    with open(lib_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == header:
                inside = True
                continue
            if inside and line.startswith("!"):
                break
            if inside:
                f = line.replace('"', "").split()
                out.append((f[0], f[1], float(f[-1])))
    return out


@pytest.mark.parametrize("set_name", sorted(CONST_E_SETS))
def test_constant_E_library_declares_a_HEH_unit(set_name):
    lib = _sets("oxidized")[set_name]["lib"]
    units = _lib_units(lib)
    assert "HEH" in units
    assert "HCO" not in units and "HCR" not in units


@pytest.mark.parametrize("set_name", sorted(CONST_E_SETS))
def test_constant_E_HEH_carries_the_oxidized_charges(set_name):
    """AMBER ships one HEH unit holding the ferric charges; the state is applied
    from the cein file at step 0. The constant_E libs must match that
    convention, and must be atom-for-atom identical to the fixed_E oxidized
    unit they were derived from."""
    heh = _lib_atoms(_sets("oxidized")[set_name]["lib"], "HEH")
    hco = _lib_atoms(_sets("oxidized")["Henriques_HCO_RESP_constpH"]["lib"], "HCO")
    assert len(heh) == 87
    assert heh == hco


def test_constant_E_and_fixed_E_share_one_frcmod_content():
    """Bonded terms are redox-independent in the Henriques/Crespo model (AMBER
    likewise ships a single frcmod.conste), so the fork must not have introduced
    a divergent parameter file."""
    sets = _sets("oxidized")
    bodies = {
        name: Path(info["frcmod"]).read_text().split("\n", 1)[1]
        for name, info in sets.items()
    }
    assert len(set(bodies.values())) == 1


# ---- the redox_state prompt is suppressed under constant_E -----------------

def test_redox_treatment_is_offered_before_redox_state():
    """The manager configures choice parameters in dict insertion order, and
    redox_treatment gates redox_state, so it must come first."""
    order = list(BisHis.get_parameter_definitions())
    assert order.index("redox_treatment") < order.index("redox_state")


def test_redox_state_is_gated_off_under_constant_E():
    assert BisHis.is_parameter_gated_off(
        "redox_state", {"redox_treatment": "constant_E", "spin_state": "low_spin"}
    ) is True


def test_redox_state_is_not_gated_under_fixed_E():
    assert BisHis.is_parameter_gated_off(
        "redox_state", {"redox_treatment": "fixed_E", "spin_state": "low_spin"}
    ) is False


def test_gated_redox_state_is_stamped_with_the_constant_E_branch():
    """A gated parameter is skipped AND stamped: redox_state is structurally
    required because it keys the metadata tree."""
    value = BisHis.gated_parameter_value(
        "redox_state", {"redox_treatment": "constant_E"}
    )
    assert value == "oxidized"
    assert BisHis.constant_E_redox_state() == "oxidized"


def test_stale_redox_state_is_redirected_under_constant_E():
    """A replayed session or a microstate combo may still carry 'reduced'.
    Without the redirect it would silently resolve back to a fixed_E set and
    re-emit HCR."""
    params = {
        "redox_treatment": "constant_E",
        "redox_state": "reduced",
        "spin_state": "low_spin",
        "ph_treatment": "constant_pH",
    }
    assert BisHis.effective_redox_state(params) == "oxidized"
    assert BisHis.select_forcefield_set_name(params) == "Henriques_HEH_RESP_constE_constpH"
    assert BisHis.get_parameter_mappings(params)["heme_name"] == "HEH"


# ---- end to end through the transformation sequence ------------------------

def _components():
    return {
        "center_chain": "A", "center_id": 100,
        "b_ring_cys_chain": "A", "b_ring_cys_id": 50,
        "c_ring_cys_chain": "A", "c_ring_cys_id": 53,
        "proximal_ligand_chain": "A", "proximal_ligand_id": 54,
        "distal_ligand_chain": "A", "distal_ligand_id": 80,
        "prop_a_chain": "A", "prop_a_id": 101,
        "prop_d_chain": "A", "prop_d_id": 102,
    }


def _final_heme_name(parameters):
    seq = BisHis.get_transformation_sequence(_components(), parameters)
    step = next(s for s in seq if s["id"] == "apply_redox_specific_heme_name")
    return step["action"]["change_residue_name"]


@pytest.mark.parametrize("ph_treatment", ["constant_pH", "fixed_pH"])
def test_sequence_emits_HEH_under_constant_E(ph_treatment):
    name = _final_heme_name({
        "redox_treatment": "constant_E", "redox_state": "oxidized",
        "spin_state": "low_spin", "ph_treatment": ph_treatment,
    })
    assert name == "HEH"


@pytest.mark.parametrize("redox,expected", [("oxidized", "HCO"), ("reduced", "HCR")])
def test_sequence_still_emits_HCO_HCR_under_fixed_E(redox, expected):
    name = _final_heme_name({
        "redox_treatment": "fixed_E", "redox_state": redox,
        "spin_state": "low_spin", "ph_treatment": "constant_pH",
    })
    assert name == expected


def test_validate_parameters_accepts_a_constant_E_configuration():
    ok, msg = BisHis.validate_parameters({
        "redox_treatment": "constant_E", "redox_state": "oxidized",
        "spin_state": "low_spin", "ph_treatment": "constant_pH",
    })
    assert ok, msg


# ---- companion pairing stays inside its own redox treatment ----------------

@pytest.mark.parametrize("start,target,expected", [
    ("Henriques_HEH_RESP_constE_constpH", "fixed_pH", "Henriques_HEH_RESP_constE_FixedpH"),
    ("Henriques_HEH_RESP_constE_FixedpH", "constant_pH", "Henriques_HEH_RESP_constE_constpH"),
    ("Henriques_HCO_RESP_constpH", "fixed_pH", "Henriques_HCO_RESP_FixedpH"),
])
def test_companion_set_pairs_within_the_redox_treatment(start, target, expected):
    """The PB-Titrate rebuild swaps a set's pH treatment. _set_base_name strips
    only the pH suffix, and _constE sits before it, so the bases still match and
    the swap can't jump across the redox axis."""
    assert find_companion_set(COFACTOR, "oxidized", "low_spin", start, target) == expected


# ---- cofactors without a redox fork are untouched --------------------------

@pytest.mark.parametrize("cofactor,redox,spin", [
    ("heme/cys_axial_b_type", "oxidized", "high_spin"),
    ("heme/his_met_axial_c_type", "oxidized", "low_spin"),
    ("heme/bis_his_b_type", "oxidized", "low_spin"),
])
def test_other_heme_leaves_declare_no_constant_E_set(cofactor, redox, spin):
    """Only the Henriques bis-His c-heme is AMBER's HEH reference compound. The
    Guberman sets are different parameterizations with different atom sets, so
    naming them HEH would hand ceinutil the wrong reference energies."""
    treatments = {
        s.get("redox_treatment") for s in discover_forcefield_files(cofactor, redox, spin)
    }
    assert "constant_E" not in treatments


def test_no_redox_fork_means_no_redox_treatment_parameter():
    """redox_treatment_parameter_definitions() returns {} below two treatments,
    so a single-treatment cofactor gains no new prompt."""
    from proprep.redoxsite_prep.transformation.transformers.cys_axial_b_type_heme import (
        CysAxialBTypeHemeTransformer as CysAxial,
    )
    assert CysAxial.redox_treatment_parameter_definitions() == {}
    assert "redox_treatment" not in CysAxial.get_parameter_definitions()


# ---- the manager actually skips the prompt (not just narrows it) -----------

def _manager():
    """A manager with __init__ bypassed and just enough state wired up."""
    import io
    from rich.console import Console
    from proprep.redoxsite_prep.transformation.redox_transformation_manager import (
        RedoxTransformationManager,
    )
    m = object.__new__(RedoxTransformationManager)
    m.console = Console(file=io.StringIO(), force_terminal=False)
    m.transformation_parameters = {}
    m.processor = None  # prompt_with_context is patched out; only the arg is read
    return m


def _run_redox_state_prompt(redox_treatment, monkeypatch):
    """Drive the real configurator for redox_state with one site already set to
    `redox_treatment`, and report whether the user was prompted."""
    import proprep.redoxsite_prep.transformation.redox_transformation_manager as mod

    prompted = []

    def _no_prompt(*a, **kw):
        prompted.append(kw.get("prompt", a[1] if len(a) > 1 else ""))
        return ""

    # The configurator imports the registry from the framework module inside the
    # function body, so patch it there rather than on the manager module.
    import proprep.redoxsite_prep.transformation.redox_transformer_framework as fw

    monkeypatch.setattr(mod, "prompt_with_context", _no_prompt)
    monkeypatch.setattr(
        fw.redox_transformer_registry, "get_transformer", lambda _n: BisHis
    )

    m = _manager()
    m.transformation_parameters = {
        "site_1": {"redox_treatment": redox_treatment, "spin_state": "low_spin"}
    }
    pdef = BisHis.get_parameter_definitions()["redox_state"]
    m._configure_choice_parameter_interactive(
        "redox_state", pdef, ["site_1"], [1], "bis_his_c_type_heme"
    )
    return prompted, m.transformation_parameters["site_1"]


def test_manager_skips_and_stamps_redox_state_under_constant_E(monkeypatch):
    prompted, params = _run_redox_state_prompt("constant_E", monkeypatch)
    assert prompted == [], "redox_state must not be asked under constant_E"
    assert params["redox_state"] == "oxidized", (
        "a gated parameter must still be stamped: redox_state keys the metadata tree"
    )


def test_manager_still_prompts_for_redox_state_under_fixed_E(monkeypatch):
    prompted, _ = _run_redox_state_prompt("fixed_E", monkeypatch)
    assert prompted, "redox_state must still be asked under fixed_E"
