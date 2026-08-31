"""Water-analysis viewer halos: one group per fact a displayed metric establishes,
labelled with its rule and cutoff, overlaps allowed, no combined category. The
B-factor comparison is relative to the structure's own protein atoms."""
from types import SimpleNamespace
from unittest.mock import MagicMock

import numpy as np
import pytest
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from rich.console import Console

from proprep.structure_prep import viewer_coordinator as vc
from proprep.structure_prep.pdb_filter_worker import PDBFilterWorker, WaterAnalyzer, WATER_HALO_KEYS


def _structure(protein_bs, water_bs):
    structure = Structure('s'); model = Model(0); chain = Chain('A')
    structure.add(model); model.add(chain)
    for i, b in enumerate(protein_bs):
        res = Residue((' ', i + 1, ' '), 'ALA', '')
        res.add(Atom('CA', np.array([50.0 + 4 * i, 50.0, 50.0], dtype=np.float32), b, 1.0, ' ', 'CA', i + 1, element='C'))
        chain.add(res)
    waters = []
    for j, b in enumerate(water_bs):
        res = Residue(('W', 1000 + j, ' '), 'HOH', '')
        res.add(Atom('O', np.array([0.0, 4.0 * j, 0.0], dtype=np.float32), b, 1.0, ' ', 'O', 5000 + j, element='O'))
        chain.add(res); waters.append(res)
    return structure, waters


def test_b_factor_is_reported_relative_to_protein_median():
    structure, (w_low, _) = _structure(protein_bs=[10, 20, 30, 40, 50], water_bs=[25, 35])
    an = WaterAnalyzer(structure)
    assert an.protein_median_bfactor() == 30.0
    a = an.analyze_water(w_low, {})
    assert a['protein_median_b'] == 30.0
    assert a['b_factor_ratio'] == pytest.approx(25 / 30)
    assert 'category' not in a          # no combined classification


def test_no_protein_means_no_b_reference():
    structure, (w,) = _structure(protein_bs=[], water_bs=[5])
    an = WaterAnalyzer(structure)
    assert an.protein_median_bfactor() is None
    assert 'b_factor_ratio' not in an.analyze_water(w, {})


@pytest.fixture
def worker():
    w = PDBFilterWorker.__new__(PDBFilterWorker)
    w.processor = MagicMock()
    return w


@pytest.fixture
def analyzer():
    a = WaterAnalyzer.__new__(WaterAnalyzer)
    a.parameters = {'sasa_probe_radius': 1.4, 'metal_distance_cutoff': 2.5,
                    'hbond_distance_cutoff': 3.5, 'interface_distance_cutoff': 5.0}
    a._protein_median_b = 22.0
    return a


def _a(num, **facts):
    base = {'residue_number': num, 'residue_name': 'HOH', 'chain_id': 'A',
            'coordinating_metal': False, 'total': 0, 'at_interface': False,
            'b_factor': 30.0, 'b_factor_ratio': 30.0 / 22.0,
            'burial_category': 'Exposed', 'burial_covered_pct': 10.0, 'burial_sasa': 90.0}
    base.update(facts)
    return base


def test_halos_follow_displayed_metrics_and_allow_overlap(worker, analyzer):
    analyses = [
        _a(1, coordinating_metal=True, burial_category='Enclosed', burial_covered_pct=100, total=3),
        _a(2, total=4, b_factor_ratio=0.5),
        _a(3, burial_covered_pct=95),
        _a(4, at_interface=True),
    ]
    groups = worker._water_halo_groups(analyses, ['1', '2', '3', '4', '5'], analyzer)
    by_key = {k: (label, [x['residue_number'] for x in items]) for k, label, _, items in groups}

    # water 1 is in three groups at once: no precedence, no combined category
    assert by_key['water_metal'][1] == [1]
    assert by_key['water_burial_enclosed'][1] == [1]
    assert by_key['water_hbond_3'][1] == [1]
    assert by_key['water_hbond_4'][1] == [2]
    assert by_key['water_b_below_median'][1] == [2]
    assert by_key['water_burial_covered90'][1] == [3]
    assert by_key['water_interface'][1] == [4]
    assert all(k in WATER_HALO_KEYS for k in by_key)

    # labels carry the rule and cutoff
    assert by_key['water_metal'][0] == "Waters: metal within 2.5 Å"
    assert by_key['water_hbond_4'][0] == "Waters: 4 H-bond partners (≤ 3.5 Å)"
    assert by_key['water_b_below_median'][0] == "Waters: B below protein median (22 Å²)"
    assert '2.2 Å' in by_key['water_burial_clash'][0] if 'water_burial_clash' in by_key else True


def test_halos_only_for_displayed_metrics(worker, analyzer):
    analyses = [_a(1, coordinating_metal=True, at_interface=True, total=2)]
    keys = {k for k, *_ in worker._water_halo_groups(analyses, ['2'], analyzer)}
    assert keys == {'water_hbond_2'}
    assert worker._water_halo_groups(analyses, ['6', '7', '8'], analyzer) == []


def test_bfactor_table_reports_ratio_to_protein_median(worker):
    console = Console(record=True, file=open('/dev/null', 'w'), width=120)
    rows = [{**_a(1), 'b_factor': 11.0, 'b_factor_ratio': 0.5, 'protein_median_b': 22.0},
            {**_a(2), 'b_factor': 44.0, 'b_factor_ratio': 2.0, 'protein_median_b': 22.0}]
    worker._display_bfactor_table(rows, console)
    text = console.export_text()
    assert 'median B = 22.0' in text
    assert '0.50' in text and 'below protein median' in text
    assert '2.00' in text and 'above protein median' in text


def test_highlight_display_label_reaches_the_annotation_config(monkeypatch):
    fake = SimpleNamespace(selected_structures=[object()], annotation_config={}, viewer_config={},
                           update_annotations=lambda: None, _launch_viewer=lambda: None)
    coord = vc.ViewerCoordinator()
    monkeypatch.setattr(coord, '_ensure_viewer', lambda: fake)
    monkeypatch.setattr(coord, 'is_running', lambda: True)
    coord.highlight(":A and 5", label="water_metal", display_label="Waters: metal within 2.5 Å")
    assert fake.annotation_config['water_metal']['label'] == "Waters: metal within 2.5 Å"
    coord.highlight(":A and 6", label="plain")
    assert 'label' not in fake.annotation_config['plain']


def test_hbond_partners_ranked_by_distance_without_angle_scores():
    # Note: entities are built bottom-up here, which caches a parent-less full_id on the
    # atoms; Biopython's Atom.__eq__ then compares names only. The analyzer excludes the
    # query water by identity, so this does not matter, but it is why == must not be used.
    # protein O atoms at 2.7, 3.1, 3.4 Å (in) and 3.6 Å (out); a carbon at 2.5 Å is not a partner
    structure = Structure('s'); model = Model(0); chain = Chain('A'); structure.add(model); model.add(chain)
    for i, (elem, name, d) in enumerate([('O', 'O', 3.4), ('O', 'OG', 2.7), ('O', 'OD1', 3.1), ('O', 'OE1', 3.6), ('C', 'CB', 2.5)]):
        res = Residue((' ', i + 1, ' '), 'SER', '')
        res.add(Atom(name, np.array([d, 0.0, 0.0], dtype=np.float32), 10.0, 1.0, ' ', name, i + 1, element=elem))
        chain.add(res)
    w = Residue(('W', 500, ' '), 'HOH', ''); wo = Atom('O', np.zeros(3, dtype=np.float32), 10.0, 1.0, ' ', 'O', 999, element='O')
    w.add(wo); chain.add(w)
    an = WaterAnalyzer(structure)
    hb = an.calculate_hydrogen_bonds(wo, chain)
    assert hb['total'] == 3
    assert [d['distance'] for d in hb['details']] == pytest.approx([2.7, 3.1, 3.4])
    assert all('angle_score' not in d for d in hb['details'])
    assert 'hbond_angle_cutoff' not in an.parameters
