"""
Importing parameters offers to build the transformer they need.

A deposited library is inert on its own: something has to recognize the site in
a structure and rename its residues to the names the library uses. That is the
transformer's job, so an import with no transformer leaves the parameters
unreachable.

The catch is that the transformer editor works on a DETECTED redox site, while
the import wizard is deliberately usable with no structure loaded. So the offer
is made when sites are present and replaced by instructions when they are not,
rather than launching into a guaranteed "No redox sites detected".
"""

import io

import pytest
from rich.console import Console

from proprep.forcefield_prep import forcefield_parameterizer as fp_mod
from proprep.forcefield_prep.forcefield_parameterizer import ForcefieldParameterizer


class _WS:
    def __init__(self, d=None):
        self.d = dict(d or {})

    def get(self, k, default=None):
        return self.d.get(k, default)

    def set(self, k, v):
        self.d[k] = v


class _RedoxModule:
    def __init__(self, boom=False):
        self.calls = []
        self.boom = boom

    def create_custom_transformer(self, forcefield_default=None):
        if self.boom:
            raise RuntimeError("editor unavailable")
        self.calls.append(forcefield_default)
        return True


class _Processor:
    def __init__(self, workspace, module=None):
        self._ws = workspace
        self._module = module

    def _get_workspace(self):
        return self._ws

    def get_module_instance(self, name):
        assert name == "Redox Site Preparer", name
        return self._module


IMPORT_RESULT = {
    "library_path": "/home/u/.proprep/forcefield_params/specialized_residues/"
                    "metal_sites/fe2s2_cys4",
    "metadata_path": "/home/u/.proprep/.../metadata.json",
    "state_dir": "/home/u/.proprep/forcefield_params/specialized_residues/"
                 "metal_sites/fe2s2_cys4/oxidized/high_spin",
    "set_name": "user_fe2s2_cys4",
    "copied_files": [],
}


def _param(workspace, module=None):
    p = ForcefieldParameterizer.__new__(ForcefieldParameterizer)
    p.console = Console(file=io.StringIO(), width=200)
    p.processor = _Processor(workspace, module)
    return p


# --------------------------------------------------------------------------- #
# the force-field link derived from the deposit
# --------------------------------------------------------------------------- #

def test_seed_is_relative_to_specialized_residues():
    seed = ForcefieldParameterizer._import_forcefield_seed(IMPORT_RESULT)

    assert seed["path"] == "metal_sites/fe2s2_cys4"
    assert seed["redox_state"] == "oxidized"
    assert seed["spin_state"] == "high_spin"


def test_seed_tolerates_a_library_path_it_cannot_place():
    seed = ForcefieldParameterizer._import_forcefield_seed(
        {"library_path": "/somewhere/else/fe2s2"})

    assert "path" not in seed


# --------------------------------------------------------------------------- #
# the offer
# --------------------------------------------------------------------------- #

def test_transformer_creator_launches_with_the_imported_library(monkeypatch):
    module = _RedoxModule()
    p = _param(_WS({"detected_redox_sites": [object()]}), module)
    monkeypatch.setattr(fp_mod, "confirm_with_context",
                        lambda processor, prompt, **kw: True)

    p._offer_transformer_for_import(IMPORT_RESULT)

    assert len(module.calls) == 1
    assert module.calls[0] == {"path": "metal_sites/fe2s2_cys4",
                               "redox_state": "oxidized",
                               "spin_state": "high_spin"}


def test_declining_leaves_instructions(monkeypatch):
    module = _RedoxModule()
    p = _param(_WS({"detected_redox_sites": [object()]}), module)
    monkeypatch.setattr(fp_mod, "confirm_with_context",
                        lambda processor, prompt, **kw: False)

    p._offer_transformer_for_import(IMPORT_RESULT)

    assert module.calls == []
    out = p.console.file.getvalue()
    assert "Redox Site Preparer" in out
    assert "metal_sites/fe2s2_cys4" in out


def test_no_detected_sites_explains_instead_of_failing(monkeypatch):
    """The import wizard runs without a structure; the editor needs a site."""
    module = _RedoxModule()
    p = _param(_WS({}), module)
    called = []
    monkeypatch.setattr(fp_mod, "confirm_with_context",
                        lambda processor, prompt, **kw: called.append(prompt) or True)

    p._offer_transformer_for_import(IMPORT_RESULT)

    assert module.calls == [], "must not launch an editor with no site to edit"
    assert called == [], "must not even ask when it cannot be done"
    out = p.console.file.getvalue()
    assert "Redox Site Detector" in out
    assert "metal_sites/fe2s2_cys4" in out


def test_editor_failure_does_not_undo_the_import(monkeypatch):
    module = _RedoxModule(boom=True)
    p = _param(_WS({"detected_redox_sites": [object()]}), module)
    monkeypatch.setattr(fp_mod, "confirm_with_context",
                        lambda processor, prompt, **kw: True)

    p._offer_transformer_for_import(IMPORT_RESULT)   # must not raise

    out = p.console.file.getvalue()
    assert "parameters are imported" in out


def test_missing_module_is_reported_not_raised(monkeypatch):
    p = _param(_WS({"detected_redox_sites": [object()]}), module=None)
    monkeypatch.setattr(fp_mod, "confirm_with_context",
                        lambda processor, prompt, **kw: True)

    p._offer_transformer_for_import(IMPORT_RESULT)

    assert "Redox Site Preparer" in p.console.file.getvalue()


# --------------------------------------------------------------------------- #
# wiring from the menu option
# --------------------------------------------------------------------------- #

def test_cancelled_import_makes_no_offer(monkeypatch):
    p = _param(_WS({"detected_redox_sites": [object()]}), _RedoxModule())
    offered = []
    p._offer_transformer_for_import = lambda result: offered.append(result)
    monkeypatch.setattr(
        "proprep.forcefield_prep.library_promotion.run_import_wizard",
        lambda console, processor: None)

    assert p._import_into_user_library() is True
    assert offered == [], "a cancelled import has nothing to build a transformer for"


def test_successful_import_makes_the_offer(monkeypatch):
    p = _param(_WS({"detected_redox_sites": [object()]}), _RedoxModule())
    offered = []
    p._offer_transformer_for_import = lambda result: offered.append(result)
    monkeypatch.setattr(
        "proprep.forcefield_prep.library_promotion.run_import_wizard",
        lambda console, processor: IMPORT_RESULT)

    assert p._import_into_user_library() is True
    assert offered == [IMPORT_RESULT]


def test_import_error_still_returns_to_the_menu(monkeypatch):
    p = _param(_WS({}), _RedoxModule())
    offered = []
    p._offer_transformer_for_import = lambda result: offered.append(result)

    def boom(console, processor):
        raise OSError("disk gone")

    monkeypatch.setattr(
        "proprep.forcefield_prep.library_promotion.run_import_wizard", boom)

    assert p._import_into_user_library() is True
    assert offered == []
    assert "Could not import parameters" in p.console.file.getvalue()


# --------------------------------------------------------------------------- #
# the seed reaches the editor
# --------------------------------------------------------------------------- #

def test_preparer_forwards_the_seed_to_the_editor(monkeypatch):
    from proprep.redoxsite_prep.redoxsite_integration import RedoxSitePreparationModule

    mod = RedoxSitePreparationModule.__new__(RedoxSitePreparationModule)
    site = type("S", (), {"site_id": "site_1", "atoms": [object()]})()
    mod.processor = _Processor(_WS({"detected_redox_sites": [site]}))
    mod.processor.console = Console(file=io.StringIO(), width=200)

    seen = {}

    class _Creator:
        def __init__(self, processor, redox_site, forcefield_default=None):
            seen["forcefield_default"] = forcefield_default

        def create(self):
            return "/transformers/x.json"

    monkeypatch.setattr(
        "proprep.redoxsite_prep.transformation.table_transformer_creator."
        "TableTransformerCreator", _Creator)

    seed = {"path": "metal_sites/fe2s2_cys4", "redox_state": "oxidized"}
    assert mod.create_custom_transformer(forcefield_default=seed) is True
    assert seen["forcefield_default"] == seed


def test_editor_stores_the_seed():
    from proprep.redoxsite_prep.transformation.table_transformer_creator import (
        TableTransformerCreator,
    )

    processor = _Processor(_WS({}))
    processor.console = Console(file=io.StringIO())
    seed = {"path": "metal_sites/fe2s2_cys4"}

    creator = TableTransformerCreator(processor, redox_site=None,
                                      forcefield_default=seed)
    assert creator.forcefield_default == seed

    # Absent seed must stay falsy, so the link prompt keeps defaulting to no.
    assert TableTransformerCreator(processor, redox_site=None).forcefield_default == {}
