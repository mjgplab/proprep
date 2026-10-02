"""The alt-loc picker prints each alternate's viewer colour next to its occupancy,
and the viewer reps are labelled with the same colour name and occupancy, so the
red/blue on screen can be matched to the numbers in the prompt."""
from unittest.mock import MagicMock

from proprep.structure_prep import altloc_picker
from proprep.structure_prep import viewer_coordinator as vc


def test_altloc_colors_are_named_and_stable():
    assert altloc_picker.altloc_color('A', 0) == ('#e74c3c', 'red')
    assert altloc_picker.altloc_color('b', 1) == ('#3498db', 'blue')
    # letters beyond the palette cycle the fallback set, by position
    assert altloc_picker.altloc_color('G', 6) == altloc_picker.FALLBACK_PALETTE[6 % 4]


def test_viewer_reps_carry_color_name_and_occupancy(monkeypatch):
    fake = MagicMock()
    monkeypatch.setattr(vc, 'viewer', fake)
    state = {'prev_labels': ['altloc_scaffold', 'altloc_A'], 'env_distance': None, 'points': None, 'owners': []}

    altloc_picker.focus_viewer(state, ('A', 37, ''), ['A', 'B'], occupancies={'A': '0.60', 'B': '0.40'})

    labels = {c.kwargs['label']: c for c in fake.highlight.call_args_list}
    assert labels['altloc_A'].kwargs['display_label'] == 'Alt A (red), occ 0.60'
    assert labels['altloc_B'].kwargs['display_label'] == 'Alt B (blue), occ 0.40'
    assert labels['altloc_A'].kwargs['color'] == '#e74c3c'
    assert labels['altloc_A'].args[0] == ':A and 37 and %A'
    # stale reps from the previous residue are cleared first
    assert [c.args[0] for c in fake.unhighlight.call_args_list] == ['altloc_scaffold', 'altloc_A']


def test_insertion_codes_reach_the_viewer_selection_and_the_shell(monkeypatch):
    fake = MagicMock()
    monkeypatch.setattr(vc, 'viewer', fake)
    owners = [('A', 40, 'A'), ('A', 40, 'A'), ('A', 40, ''), ('B', 7, ''), ('A', 90, '')]
    points = [[0, 0, 0], [1, 0, 0], [3, 0, 0], [0, 4, 0], [20, 0, 0]]
    state = {'prev_labels': [], 'env_distance': 5.0, 'points': __import__('numpy').array(points, float),
             'owners': owners}

    altloc_picker.focus_viewer(state, ('A', 40, 'A'), ['A', 'B'])

    labels = {c.kwargs['label']: c for c in fake.highlight.call_args_list}
    assert labels['altloc_B'].args[0] == ':A and 40^A and %B'
    # residue 40 (no insertion code) and B:7 are within 5 A; A:90 is not; 40A itself is the target
    assert labels['altloc_environment'].args[0] == '(:A and (40)) or (:B and (7))'
