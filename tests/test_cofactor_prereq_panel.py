"""
The cofactor-prerequisite panel must actually run.

It was rewritten to report declared prerequisites instead of narrating
chemistry, and the rewrite called ``self._collect_cofactor_prereqs()`` -- a
method that does not exist. The real one is
``_collect_cofactor_prereq_groups``. The whole suite passed, because it
exercised the inference function and nothing ever invoked the panel:

    Error executing Show Generate tLEaP input for single state menu:
    'TLeapInputGenerator' object has no attribute '_collect_cofactor_prereqs'

An AttributeError on a method name is invisible until the line runs, so a
panel with no test has no coverage of its own wiring.
"""

from types import SimpleNamespace

import pytest
from rich.console import Console

from proprep.tleap_prep.tleap_input_generator import TLeapInputGenerator


def _generator(cofactors, width=100):
    generator = TLeapInputGenerator.__new__(TLeapInputGenerator)
    generator.processor = SimpleNamespace(console=Console(record=True, width=width))
    generator._collect_cofactor_prereq_groups = lambda: cofactors
    return generator


COFACTORS = [
    {"residue_name": "FES", "location": "A.1310",
     "groups": [["leaprc.protein.ff19SB", "leaprc.protein.ff14SB"]]},
    {"residue_name": "MTE", "location": "A.1311",
     "groups": [["leaprc.gaff2"]]},
]


# --------------------------------------------------------------------------- #
# it runs at all
# --------------------------------------------------------------------------- #

def test_the_panel_runs():
    """The regression: it raised AttributeError before producing any output."""
    generator = _generator(COFACTORS)

    generator._show_cofactor_ff_prerequisites_panel()   # must not raise

    assert generator.processor.console.export_text()


def test_it_calls_the_method_that_exists():
    assert hasattr(TLeapInputGenerator, "_collect_cofactor_prereq_groups")
    assert not hasattr(TLeapInputGenerator, "_collect_cofactor_prereqs")


# --------------------------------------------------------------------------- #
# what it reports
# --------------------------------------------------------------------------- #

def test_each_cofactor_and_its_requirement_appear():
    generator = _generator(COFACTORS)
    generator._show_cofactor_ff_prerequisites_panel()
    text = generator.processor.console.export_text()

    for fragment in ("FES", "A.1310", "MTE", "leaprc.gaff2"):
        assert fragment in text


def test_alternatives_within_a_group_are_shown_as_alternatives():
    generator = _generator(COFACTORS)
    generator._show_cofactor_ff_prerequisites_panel()

    assert " or " in generator.processor.console.export_text()


def test_a_cofactor_with_two_groups_gets_a_row_each():
    """MoCo needs a protein FF AND gaff2; both must be visible."""
    generator = _generator([{
        "residue_name": "MOS", "location": "A.1312",
        "groups": [["leaprc.protein.ff19SB"], ["leaprc.gaff2"]],
    }])
    generator._show_cofactor_ff_prerequisites_panel()
    text = generator.processor.console.export_text()

    assert "leaprc.protein.ff19SB" in text
    assert "leaprc.gaff2" in text


def test_nothing_is_printed_when_no_cofactor_declares_anything():
    """Silence is right: there is nothing to satisfy."""
    generator = _generator([])

    generator._show_cofactor_ff_prerequisites_panel()

    assert generator.processor.console.export_text().strip() == ""


def test_a_missing_location_does_not_break_the_row():
    generator = _generator([{"residue_name": "LIG", "location": "",
                             "groups": [["leaprc.gaff2"]]}])

    generator._show_cofactor_ff_prerequisites_panel()

    assert "LIG" in generator.processor.console.export_text()


# --------------------------------------------------------------------------- #
# the prose it replaced must stay gone
# --------------------------------------------------------------------------- #

def test_no_invented_chemistry_remains():
    """
    The old panel explained a ribitol tail and which bond would fail -- neither
    knowable about an arbitrary parameter set, and both wrong for an imported
    one.
    """
    generator = _generator(COFACTORS)
    generator._show_cofactor_ff_prerequisites_panel()
    text = generator.processor.console.export_text().lower()

    for claim in ("ribitol", "sugar chain", "ch2/ch/oh"):
        assert claim not in text


# --------------------------------------------------------------------------- #
# no other self.<method> call is a typo
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("module_name,class_name", [
    ("proprep.tleap_prep.tleap_input_generator", "TLeapInputGenerator"),
    ("proprep.forcefield_prep.structure_preprocessor", "StructurePreprocessor"),
    ("proprep.forcefield_prep.metal_site_parameterizer", "MetalSiteWorkflowManager"),
    ("proprep.redoxsite_prep.transformation.table_transformer_creator",
     "TableTransformerCreator"),
    ("proprep.redoxsite_prep.transformation.redox_transformation_manager",
     "RedoxTransformationManager"),
    ("proprep.forcefield_prep.mcpb.mol2_writer", "Mol2Writer"),
    ("proprep.forcefield_prep.mcpb.atom_typer", "MCPBAtomTyper"),
    # Added once its three unresolved calls were fixed: _get_site_attribute
    # (defined on PDBProcessor, not here) and the two resume methods that were
    # dispatched to but never written.
    ("proprep.forcefield_prep.forcefield_parameterizer", "ForcefieldParameterizer"),
])
def test_every_self_method_call_resolves(module_name, class_name):
    """
    Catches the mistake statically instead of waiting for the line to run.

    A misspelled self.<method> is a clean AttributeError at call time and
    invisible before it, which is how a panel shipped calling a method that
    never existed. Attribute ACCESS is not checked -- only calls -- because
    instance attributes are legitimately set outside __init__ here.
    """
    import ast
    import importlib
    import inspect

    cls = getattr(importlib.import_module(module_name), class_name)
    tree = ast.parse(inspect.getsource(cls))

    called = {
        node.func.attr
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and isinstance(node.func.value, ast.Name)
        and node.func.value.id == "self"
    }

    missing = sorted(name for name in called if not hasattr(cls, name))
    assert not missing, f"{class_name} calls undefined method(s): {missing}"
