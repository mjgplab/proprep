"""The alt-loc picker prints each alternate's viewer colour next to its occupancy,
and the viewer reps are labelled with the same colour name and occupancy, so the
red/blue on screen can be matched to the numbers in the prompt."""
from unittest.mock import MagicMock

import pytest

from proprep.structure_prep import structure_completeness as sc
from proprep.structure_prep import viewer_coordinator as vc

Owner = next(cls for cls in vars(sc).values()
             if isinstance(cls, type) and hasattr(cls, '_focus_altloc_viewer'))


@pytest.fixture
def fixer():
    return Owner.__new__(Owner)


def test_altloc_colors_are_named_and_stable(fixer):
    assert fixer._altloc_color('A', 0) == ('#e74c3c', 'red')
    assert fixer._altloc_color('b', 1) == ('#3498db', 'blue')
    # letters beyond the palette cycle the fallback set, by position
    assert fixer._altloc_color('G', 6) == fixer._ALTLOC_FALLBACK_PALETTE[6 % 4]


def test_viewer_reps_carry_color_name_and_occupancy(fixer, monkeypatch):
    fake = MagicMock()
    monkeypatch.setattr(vc, 'viewer', fake)
    state = {'prev_labels': ['altloc_scaffold', 'altloc_A'], 'env_distance': None, 'neighbor_search': None}

    fixer._focus_altloc_viewer(state, 'A', 'SER', 37, ['A', 'B'], occupancies={'A': '0.60', 'B': '0.40'})

    labels = {c.kwargs['label']: c.kwargs for c in fake.highlight.call_args_list}
    assert labels['altloc_A']['display_label'] == 'Alt A (red), occ 0.60'
    assert labels['altloc_B']['display_label'] == 'Alt B (blue), occ 0.40'
    assert labels['altloc_A']['color'] == '#e74c3c'
    assert labels['altloc_A']['selection'] if 'selection' in labels['altloc_A'] else True
    # stale reps from the previous residue are cleared first
    assert [c.args[0] for c in fake.unhighlight.call_args_list] == ['altloc_scaffold', 'altloc_A']
