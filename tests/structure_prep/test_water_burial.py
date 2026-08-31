"""Water burial from geometry: accessible area of the oxygen against the protein,
bulk connectivity by flood fill, and the derived neighbourhood cutoff.

Structures are built synthetically so the expected answers follow from the
geometry with no reference data.
"""
import math

import numpy as np
import pytest
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from proprep.structure_prep.pdb_filter_worker import WaterAnalyzer

PROBE = 1.4
R_W = WaterAnalyzer.WATER_OXYGEN_RADIUS
ISOLATED = 4 * math.pi * (R_W + PROBE) ** 2


def _fibonacci_sphere(n, radius):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = math.pi * (1 + 5 ** 0.5) * i
    return radius * np.stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)], axis=1)


def _build(protein_coords, water_coords, resname='ALA', element='C'):
    """One chain of single-atom 'protein' residues plus HOH residues."""
    structure = Structure('s'); model = Model(0); chain = Chain('A')
    structure.add(model); model.add(chain)
    for i, c in enumerate(protein_coords):
        res = Residue((' ', i + 1, ' '), resname, '')
        res.add(Atom('CA', np.asarray(c, dtype=np.float32), 20.0, 1.0, ' ', 'CA', i + 1, element=element))
        chain.add(res)
    waters = []
    for j, c in enumerate(water_coords):
        res = Residue(('W', 1000 + j, ' '), 'HOH', '')
        res.add(Atom('O', np.asarray(c, dtype=np.float32), 20.0, 1.0, ' ', 'O', 5000 + j, element='O'))
        chain.add(res)
        waters.append(res)
    return structure, waters


def test_isolated_water_has_full_sphere_and_is_bulk():
    structure, (w,) = _build([[50.0, 50.0, 50.0]], [[0.0, 0.0, 0.0]])
    a = WaterAnalyzer(structure).calculate_burial_analysis(w)
    assert a['burial_sasa'] == pytest.approx(ISOLATED, rel=1e-3)
    assert a['burial_covered_pct'] == pytest.approx(0.0, abs=0.1)
    assert a['burial_access'] == 'bulk'
    assert a['burial_category'] == 'Exposed'


def test_water_in_closed_shell_is_enclosed_with_zero_area():
    # 400 carbons on a 5 Å sphere: neighbour spacing ~0.9 Å, far below what a 1.4 Å probe needs
    structure, (w,) = _build(_fibonacci_sphere(400, 5.0), [[0.0, 0.0, 0.0]])
    a = WaterAnalyzer(structure).calculate_burial_analysis(w)
    assert a['burial_sasa'] == 0.0
    assert a['burial_covered_pct'] == 100.0
    assert a['burial_access'] == 'enclosed'
    assert a['burial_category'] == 'Enclosed'
    assert a['burial_closest_distance'] == pytest.approx(5.0, abs=1e-3)


def test_shell_with_a_probe_sized_opening_is_bulk():
    # Same shell with the cap z > 2.5 removed: opening radius sqrt(25 - 6.25) = 4.3 Å,
    # so a probe centre (needs 1.7 + 1.4 = 3.1 Å clearance from the rim atoms) passes.
    shell = _fibonacci_sphere(400, 5.0)
    structure, (w,) = _build(shell[shell[:, 2] <= 2.5], [[0.0, 0.0, 0.0]])
    a = WaterAnalyzer(structure).calculate_burial_analysis(w)
    assert a['burial_access'] == 'bulk'
    assert a['burial_sasa'] > 0.0
    assert a['burial_category'] == 'Exposed'


def test_water_on_a_flat_slab_is_partly_exposed():
    xs, ys = np.meshgrid(np.arange(-12, 12.1, 1.5), np.arange(-12, 12.1, 1.5))
    slab = np.stack([xs.ravel(), ys.ravel(), np.zeros(xs.size)], axis=1)
    structure, (w,) = _build(slab, [[0.0, 0.0, 3.3]])
    a = WaterAnalyzer(structure).calculate_burial_analysis(w)
    assert 0.0 < a['burial_sasa'] < ISOLATED
    assert a['burial_access'] == 'bulk'


def test_clash_flag_uses_wwpdb_close_contact_distance_and_ignores_metals():
    structure, (w,) = _build([[2.0, 0.0, 0.0]], [[0.0, 0.0, 0.0]])
    assert WaterAnalyzer(structure).calculate_burial_analysis(w)['burial_category'] == 'Clash'
    structure, (w,) = _build([[2.0, 0.0, 0.0]], [[0.0, 0.0, 0.0]], resname='ZN', element='ZN')
    assert WaterAnalyzer(structure).calculate_burial_analysis(w)['burial_category'] != 'Clash'


def test_neighbourhood_cutoff_is_derived_from_the_radii():
    structure, _ = _build(_fibonacci_sphere(50, 8.0), [[0.0, 0.0, 0.0]], element='S')
    ctx = WaterAnalyzer(structure)._burial_context()
    # r_water + r_max + 2 * probe: beyond this no atom's accessible sphere meets the water's
    assert ctx['sasa_cutoff'] == pytest.approx(R_W + WaterAnalyzer.SASA_RADII['S'] + 2 * PROBE)


def test_neighbourhood_sasa_matches_the_full_calculation():
    freesasa = pytest.importorskip('freesasa')
    rng = np.random.default_rng(7)
    coords = rng.uniform(-15, 15, size=(300, 3))
    structure, (w,) = _build(coords, [[0.0, 0.0, 0.0]])
    an = WaterAnalyzer(structure)
    local = an.calculate_water_sasa(w['O'].coord)
    radii = [WaterAnalyzer.SASA_RADII['C']] * len(coords) + [R_W]
    full = freesasa.calcCoord(np.vstack([coords, [[0, 0, 0]]]).ravel().tolist(), radii,
                              freesasa.Parameters({'probe-radius': PROBE, 'algorithm': freesasa.LeeRichards,
                                                   'n-slices': WaterAnalyzer.SASA_SLICES}))
    assert local == pytest.approx(full.atomArea(len(coords)), abs=1e-6)


def test_other_waters_do_not_occlude_by_default():
    shell = _fibonacci_sphere(400, 5.0)
    structure, waters = _build([[50.0, 50.0, 50.0]], np.vstack([[[0.0, 0.0, 0.0]], shell]))
    an = WaterAnalyzer(structure)
    assert an.calculate_burial_analysis(waters[0])['burial_sasa'] == pytest.approx(ISOLATED, rel=1e-3)
    an.set_parameters(burial_atom_types='protein,hetero,water')
    assert an.calculate_burial_analysis(waters[0])['burial_sasa'] == 0.0
