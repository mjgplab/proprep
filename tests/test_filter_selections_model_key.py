"""filter_selections must stay {chain_id: {component_type: [resid, ...]}}.

Up to 1.19.0 PDB Filter also wrote the chosen model index into that dict
as an int under "selected_model" for multi-model structures. The
Protonation State Analyzer iterated the values and crashed with
"'int' object has no attribute 'items'" on every NMR ensemble where one
model was picked (1UAO, 18 models). The index now lives under its own
workspace key; chain_selections() shields consumers from the old shape.
"""

from __future__ import annotations

import re
from pathlib import Path

from proprep.utils.workspace import chain_selections

LEGACY = {"A": {"protein": [1, 2, 3]}, "B": {"protein": [4]}, "selected_model": 0}


def test_chain_selections_drops_legacy_model_key():
    assert chain_selections(LEGACY) == {"A": {"protein": [1, 2, 3]}, "B": {"protein": [4]}}


def test_chain_selections_is_identity_on_clean_dict():
    clean = {"A": {"protein": [1]}}
    assert chain_selections(clean) == clean


def test_chain_selections_handles_missing():
    assert chain_selections(None) == {}
    assert chain_selections({}) == {}


def test_consumers_can_iterate_legacy_value():
    """The exact loop shape from the analyzer, run on the legacy value."""
    for chain_id, component_selections in chain_selections(LEGACY).items():
        for comp_type, residues in component_selections.items():
            assert isinstance(residues, list), (chain_id, comp_type)


def test_pdb_filter_no_longer_writes_the_model_into_filter_selections():
    """The producer stores the index under selected_model_idx, not inside
    the chain dict. Source-level check: the assignment must be gone and the
    workspace key present."""
    src = Path(__file__).resolve().parents[1] / "src/proprep/structure_prep/pdb_filter.py"
    text = src.read_text()
    assert not re.search(r'filter_selections\[\s*["\']selected_model["\']\s*\]\s*=', text)
    assert 'workspace.set("selected_model_idx", selected_model_idx)' in text
