"""Regression: the per-metric water analysis tables must read analysis parameters
from the WaterAnalyzer, not from PDBFilterWorker (which has no ``parameters``).

Metric 4 crashed with ``'PDBFilterWorker' object has no attribute 'parameters'``
because ``_display_sasa_table`` was the one table that was not handed the analyzer.
"""
from unittest.mock import MagicMock

import pytest
from rich.console import Console

from proprep.structure_prep.pdb_filter_worker import PDBFilterWorker, WaterAnalyzer


@pytest.fixture
def worker():
    w = PDBFilterWorker.__new__(PDBFilterWorker)  # display methods need no __init__ state
    w.processor = MagicMock()
    return w


@pytest.fixture
def analyzer():
    a = WaterAnalyzer.__new__(WaterAnalyzer)
    a.parameters = {'sasa_probe_radius': 1.4, 'burial_atom_types': 'protein,hetero'}
    return a


def _analysis(resnum, sasa, category, access='bulk'):
    return {
        'residue_number': resnum, 'residue_name': 'HOH', 'chain_id': 'A',
        'burial_sasa': sasa, 'burial_sasa_isolated': 98.5, 'burial_access': access,
        'burial_closest_distance': 2.8, 'burial_closest_atom': 'SER72 OG',
        'burial_category': category,
    }


def test_burial_table_uses_analyzer_parameters(worker, analyzer):
    console = Console(file=open('/dev/null', 'w'), width=120)
    analyses = [_analysis(301, 0.0, 'Enclosed', 'enclosed'), _analysis(302, 35.0, 'Exposed')]

    # Must not raise; before the fix this was an AttributeError on self.parameters
    worker._display_water_analysis_table(analyses, ['4'], console, analyzer)


def test_burial_table_reports_the_analyzer_probe(worker, analyzer):
    analyzer.parameters['sasa_probe_radius'] = 1.6
    console = Console(record=True, file=open('/dev/null', 'w'), width=120)

    worker._display_sasa_table([_analysis(301, 12.0, 'Exposed')], console, analyzer)

    text = console.export_text()
    assert '1.6 Å probe' in text
    assert 'SER72 OG 2.80 Å' in text


def test_burial_table_shows_covered_percent_and_the_rules(worker, analyzer):
    console = Console(record=True, file=open('/dev/null', 'w'), width=120)
    row = _analysis(301, 0.6, 'Exposed')          # 0.6 of 98.5 Å² -> 99% covered

    worker._display_sasa_table([row], console, analyzer)

    text = console.export_text()
    assert '99%' in text
    assert 'Clash' in text and '2.2 Å' in text          # the clash cutoff is stated
    assert '2.80 Å' in text                              # r_water + probe touching distance (1.4 + 1.4)
    assert '3.30 Å' in text                              # contact distance + one grid cell: the Access reach
    assert 'SASA = 0.0' in text                          # the buried rule
    assert 'two different distances' in text
    assert 'first that matches wins' in text


def test_profile_and_directional_tables_state_their_label_conventions(worker, monkeypatch):
    import proprep.utils.prompts as prompts
    monkeypatch.setattr(prompts, 'confirm_with_context', lambda *a, **k: False)   # no per-water chart
    monkeypatch.setattr(prompts, 'prompt_with_context', lambda *a, **k: 'done')
    console = Console(record=True, file=open('/dev/null', 'w'), width=140)
    profile_row = {'residue_number': 301, 'residue_name': 'HOH',
                   'burial_profile': {'final_count': 30.0, 'saturation_radius': 5.0, 'saturation_count': 28.0,
                                      'steepest_start': 3.0, 'steepest_end': 3.5, 'radii': [2.0, 2.5],
                                      'burial_counts': [1.0, 2.0]}}
    worker._display_burial_profile_table([profile_row], console)
    text = console.export_text()
    assert '<10%/step' in text and 'less than 10%' in text

    console = Console(record=True, file=open('/dev/null', 'w'), width=160)
    directional_row = {'residue_number': 302, 'residue_name': 'HOH',
                       'directional_burial': {'total_weight': 12.0, 'primary_direction': 'N',
                                              'pocket_opening': 'S',
                                              'pattern_type': 'Highly directional (range > 1.5 × mean)'}}
    worker._display_directional_analysis_table([directional_row], console)
    text = console.export_text()
    assert 'range > 1.5 × mean' in text          # the label carries its rule
    assert '0.5×' in text and '1.5×' in text     # and the legend restates the convention
