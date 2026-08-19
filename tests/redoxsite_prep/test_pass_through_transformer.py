"""
A transformer with no edits, carrying only a force-field link.

Renaming and parameter-binding are separate jobs. MCPB output needs both (CYS
-> CM1, plus the library); an externally parameterized cofactor already named
correctly in the PDB needs only the second. There was no way to express the
second: the transformer creator refused to save without operations, and the
built-in no_transformation pass-through declares FORCEFIELD_PATH = None and is
skipped by the Topology Generator's force-field selection.

So "needs no renaming" and "has parameters" were mutually exclusive, and a
cofactor like an imported FAD had nowhere to live -- the deposited library was
unreachable because nothing bound it to the site.

Roles are seeded from the site when the builder is constructed, so a spec with
no operations still records WHICH site it applies to. That is what separates it
from no_transformation, which matches every site and therefore could never
safely carry parameters -- it would apply one residue's library to all of them.
"""

from types import SimpleNamespace

import pytest

from proprep.redoxsite_prep.transformation.spec_transformer import (
    build_spec_transformer_class,
)


def _spec(roles, operations=(), forcefield=None, **extra):
    spec = {
        "schema_version": 1,
        "name": extra.pop("name", "fad_cofactor"),
        "description": "FAD cofactor",
        "source": "table_transformer_creator",
        "roles": list(roles),
        "operations": list(operations),
    }
    if forcefield:
        spec["forcefield"] = forcefield
    spec.update(extra)
    return spec


def _role(resname, label=None):
    return {"label": label or resname.lower(), "resname": resname,
            "fingerprint": None, "discriminators": []}


def _site(site_id, residues):
    atoms = [SimpleNamespace(chain="A", resid=100 + i, resname=r,
                             atom_name="X", coords=(0.0, 0.0, 0.0))
             for i, r in enumerate(residues)]
    return SimpleNamespace(site_id=site_id, atoms=atoms)


PASS_THROUGH = _spec([_role("FAD")], forcefield={"path": "small_molecules/FAD"},
                     pass_through=True)


# --------------------------------------------------------------------------- #
# it binds parameters without renaming
# --------------------------------------------------------------------------- #

def test_it_carries_the_library_path():
    """FORCEFIELD_PATH is what the Topology Generator resolves parameters through."""
    cls = build_spec_transformer_class(PASS_THROUGH)

    assert cls.FORCEFIELD_PATH == "small_molecules/FAD"


def test_it_renames_nothing():
    cls = build_spec_transformer_class(PASS_THROUGH)

    assert cls.SPEC["operations"] == []


def test_it_matches_the_site_it_was_built_from():
    cls = build_spec_transformer_class(PASS_THROUGH)

    evaluation = cls.evaluate_redox_site(_site("fad_site", ["FAD"]))

    assert evaluation.is_valid
    assert evaluation.confidence == 1.0


@pytest.mark.parametrize("residues", [["CYM", "FES"], ["MTE", "MOS"], ["HEM"]])
def test_it_does_not_claim_unrelated_sites(residues):
    """
    The critical difference from no_transformation, which is is_valid=True for
    everything. A pass-through that matched every site would apply one
    residue's library to all of them.
    """
    cls = build_spec_transformer_class(PASS_THROUGH)

    assert not cls.evaluate_redox_site(_site("other", residues)).is_valid


def test_the_built_in_pass_through_still_carries_no_parameters():
    """It matches everything, so it must not gain a library path."""
    from proprep.redoxsite_prep.transformation.transformers.no_transformation import (
        NoTransformationTransformer,
    )

    assert NoTransformationTransformer.FORCEFIELD_PATH is None
    assert NoTransformationTransformer.evaluate_redox_site(
        _site("anything", ["FAD"])).is_valid


# --------------------------------------------------------------------------- #
# roles are what make it specific
# --------------------------------------------------------------------------- #

def test_roles_survive_into_the_match_criteria():
    cls = build_spec_transformer_class(PASS_THROUGH)

    assert cls.RENAMES == [{"resname": "FAD", "target": "fad", "signature": None}]


def test_a_spec_with_no_roles_matches_nothing():
    """
    Which is why the creator refuses to save one: total == 0 makes is_valid
    False, so it would be an inert file.
    """
    cls = build_spec_transformer_class(
        _spec([], forcefield={"path": "small_molecules/FAD"}))

    assert not cls.evaluate_redox_site(_site("fad_site", ["FAD"])).is_valid


def test_a_multi_residue_site_requires_all_of_them():
    cls = build_spec_transformer_class(
        _spec([_role("MTE"), _role("MOS")],
              forcefield={"path": "metal_sites/moco"}))

    assert cls.evaluate_redox_site(_site("s", ["MTE", "MOS"])).is_valid
    assert not cls.evaluate_redox_site(_site("s", ["MTE"])).is_valid


# --------------------------------------------------------------------------- #
# a renaming transformer is unaffected
# --------------------------------------------------------------------------- #

def test_a_renaming_spec_still_works():
    """The MCPB case: renames AND a library, which already worked."""
    spec = _spec(
        [_role("CYS", "cys_1")],
        operations=[{"op": "rename_residue", "selector": {"role": "cys_1"},
                     "action": {"change_residue_name": "CM1"}}],
        forcefield={"path": "metal_sites/4hux_fe2s2"}, name="fe2s2")
    cls = build_spec_transformer_class(spec)

    assert cls.FORCEFIELD_PATH == "metal_sites/4hux_fe2s2"
    assert cls.SPEC["operations"]
    assert cls.evaluate_redox_site(_site("s", ["CYS"])).is_valid


# --------------------------------------------------------------------------- #
# saving one
# --------------------------------------------------------------------------- #

def _creator(monkeypatch, tmp_path, *, prompts, confirms):
    """A TableTransformerCreator with its prompts scripted and output redirected."""
    from rich.console import Console
    from proprep.redoxsite_prep.transformation import table_transformer_creator as ttc

    pq, cq = list(prompts), list(confirms)
    monkeypatch.setattr(ttc, "prompt_with_context",
                        lambda *a, default=None, **k: (
                            pq.pop(0) if pq and pq[0] is not None
                            else (pq.pop(0), default)[1]))
    monkeypatch.setattr(ttc, "confirm_with_context",
                        lambda *a, default=None, **k: cq.pop(0))
    monkeypatch.setattr(ttc, "DEFAULT_USER_TRANSFORMER_DIR", tmp_path)

    creator = ttc.TableTransformerCreator.__new__(ttc.TableTransformerCreator)
    creator.console = Console(quiet=True)
    creator.processor = None
    creator.forcefield_default = {}
    return creator, ttc


class _Builder:
    """Minimal RecipeBuilder stand-in: operations plus seeded roles."""

    def __init__(self, operations, role_meta):
        self.operations = list(operations)
        self.role_meta = dict(role_meta)
        self.parameters = []
        self.reference_state = {}

    def to_recipe(self, name, description="", forcefield=None):
        recipe = {"name": name, "description": description,
                  "roles": [{"label": k, **v} for k, v in self.role_meta.items()],
                  "operations": self.operations}
        if forcefield:
            recipe["forcefield"] = forcefield
        return recipe


ROLE_META = {"fad": {"resname": "FAD", "fingerprint": None, "discriminators": []}}


def test_a_pass_through_with_a_library_link_is_saved(monkeypatch, tmp_path):
    """The reported case: FAD already named correctly, parameters imported."""
    import json

    creator, _ttc = _creator(
        monkeypatch, tmp_path,
        # template name, forcefield path, redox_state, spin_state
        prompts=["FAD cofactor", "small_molecules/FAD", "", ""],
        # save-as-pass-through?, link a library?
        confirms=[True, True])

    out = creator._save(_Builder([], ROLE_META))

    assert out is not None
    recipe = json.loads((tmp_path / "fad_cofactor.json").read_text())
    assert recipe["operations"] == []
    assert recipe["forcefield"]["path"] == "small_molecules/FAD"
    assert recipe["pass_through"] is True
    assert [r["resname"] for r in recipe["roles"]] == ["FAD"]


def test_a_pass_through_without_a_link_is_refused(monkeypatch, tmp_path):
    """
    It would rename nothing and bind nothing -- which no_transformation already
    does, and more broadly. Writing an inert file helps no one.
    """
    creator, _ttc = _creator(
        monkeypatch, tmp_path,
        prompts=["FAD cofactor"],
        confirms=[True, False])          # save as pass-through, but no link

    assert creator._save(_Builder([], ROLE_META)) is None
    assert list(tmp_path.glob("*.json")) == []


def test_declining_the_pass_through_saves_nothing(monkeypatch, tmp_path):
    creator, _ttc = _creator(monkeypatch, tmp_path, prompts=[], confirms=[False])

    assert creator._save(_Builder([], ROLE_META)) is None


def test_no_operations_and_no_roles_is_still_refused(monkeypatch, tmp_path):
    """Nothing to match on: the spec could never apply to anything."""
    creator, _ttc = _creator(monkeypatch, tmp_path, prompts=[], confirms=[])

    assert creator._save(_Builder([], {})) is None


# --------------------------------------------------------------------------- #
# choosing the library, rather than typing its path blind
# --------------------------------------------------------------------------- #

def _picker(monkeypatch, answers, *, entries=None, metadata=None):
    from rich.console import Console
    from proprep.redoxsite_prep.transformation import table_transformer_creator as ttc

    queue = list(answers)
    monkeypatch.setattr(
        ttc, "prompt_with_context",
        lambda *a, default=None, **k: (queue.pop(0) if queue else (default or "")))
    if entries is not None:
        monkeypatch.setattr(
            "proprep.forcefield_params.loader.get_available_cofactor_types",
            lambda: {e: {} for e in entries})
    if metadata is not None:
        monkeypatch.setattr(
            "proprep.forcefield_params.loader.load_forcefield_metadata",
            lambda path: metadata)

    creator = ttc.TableTransformerCreator.__new__(ttc.TableTransformerCreator)
    creator.console = Console(quiet=True)
    creator.processor = None
    creator.forcefield_default = {}
    return creator


ENTRIES = ["heme/bis_his_c_type", "metal_sites/4hux_fe2s2", "small_molecules/FAD"]
FAD_META = {"redox_states": {"single_state": {"spin_states": {"default": {}}}}}
HEME_META = {"redox_states": {
    "oxidized": {"spin_states": {"low_spin": {}, "high_spin": {}}},
    "reduced": {"spin_states": {"low_spin": {}}}}}


def test_a_library_can_be_chosen_by_number(monkeypatch):
    creator = _picker(monkeypatch, ["3"], entries=ENTRIES)

    assert creator._pick_forcefield_path() == "small_molecules/FAD"


def test_a_path_can_still_be_typed(monkeypatch):
    """Entries outside the library, or a fresh deposit the listing missed."""
    creator = _picker(monkeypatch, ["some/other/path"], entries=ENTRIES)

    assert creator._pick_forcefield_path() == "some/other/path"


def test_the_just_deposited_entry_is_the_default(monkeypatch):
    creator = _picker(monkeypatch, [], entries=ENTRIES)

    assert creator._pick_forcefield_path(
        default="small_molecules/FAD") == "small_molecules/FAD"


def test_typing_is_the_fallback_when_nothing_is_listed(monkeypatch):
    creator = _picker(monkeypatch, ["heme/bis_his_c_type"], entries=[])

    assert creator._pick_forcefield_path() == "heme/bis_his_c_type"


def test_a_single_state_entry_needs_no_question(monkeypatch):
    """
    The reported confusion: small_molecules/FAD lives at single_state/default,
    so the redox slot is "single_state" -- not the "default" a reasonable
    person guesses. With one state there is nothing to ask.
    """
    creator = _picker(monkeypatch, [], metadata=FAD_META)

    assert creator._pick_forcefield_state("small_molecules/FAD") == (
        "single_state", "default")


def test_a_multi_state_entry_lists_its_leaves(monkeypatch):
    creator = _picker(monkeypatch, ["2"], metadata=HEME_META)

    assert creator._pick_forcefield_state("heme/bis_his_c_type") == (
        "oxidized", "low_spin")


def test_the_seeded_state_is_the_default(monkeypatch):
    creator = _picker(monkeypatch, [], metadata=HEME_META)

    assert creator._pick_forcefield_state(
        "heme/bis_his_c_type", "reduced", "low_spin") == ("reduced", "low_spin")


def test_unreadable_metadata_falls_back_to_typing(monkeypatch):
    creator = _picker(monkeypatch, ["oxidized", "high_spin"], metadata={})

    assert creator._pick_forcefield_state("x/y") == ("oxidized", "high_spin")
