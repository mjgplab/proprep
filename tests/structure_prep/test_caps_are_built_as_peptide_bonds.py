"""ACE and NME caps are peptide bonds onto the residue, and are built as such.

6R2Q: after the structure fixer, tLEaP reported the carbonyl O of an ACE 0.41 A
from the next residue's CB, and another cap's CH3 0.75 A from a backbone C. The
caps were placed at a fixed offset along the x axis of the file's frame from
the residue's first atom (for NME, from its LAST atom, usually a side-chain
atom), whatever way the residue faced, on the belief that "tLEaP rebuilds these
coordinates regardless". It does not: it builds only missing atoms. Measured on
the real structure, the caps had C-N-CA angles of 166, 54 and 50 degrees (121.7
is a peptide's) and omega torsions of -25, -157 and 4 (180 is trans).

The residue below is real: LYS 2 of chain B of the repaired 6R2Q, backbone and
CB, with the neighbouring atoms that the old cap ran into.
"""

import io

import numpy as np
import pytest
from rich.console import Console

from proprep.structure_prep.structure_completeness import CappingHandler

RESIDUE = """\
ATOM      1  N   LYS B   2      19.188  63.329  26.931  1.00  0.00           N
ATOM      2  CA  LYS B   2      20.503  63.959  26.903  1.00  0.00           C
ATOM      3  C   LYS B   2      21.594  62.905  26.748  1.00  0.00           C
ATOM      4  O   LYS B   2      21.455  61.951  25.981  1.00  0.00           O
ATOM      5  CB  LYS B   2      20.574  64.989  25.772  1.00  0.00           C
""".splitlines(keepends=True)


def _xyz(line):
    return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])


def _atoms(lines):
    return {l[12:16].strip(): _xyz(l) for l in lines}


def _angle(a, b, c):
    u, v = a - b, c - b
    return float(np.degrees(np.arccos(u @ v / np.linalg.norm(u) / np.linalg.norm(v))))


def _torsion(a, b, c, d):
    b1, b2, b3 = b - a, c - b, d - c
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    return float(np.degrees(np.arctan2(np.cross(n1, n2) @ (b2 / np.linalg.norm(b2)), n1 @ n2)))


@pytest.fixture
def handler():
    return CappingHandler(Console(file=io.StringIO(), width=200))


def test_ace_is_a_trans_peptide_bond_onto_the_residues_nitrogen(handler):
    residue = _atoms(RESIDUE)
    everything = np.array([_xyz(l) for l in RESIDUE])
    cap = _atoms(handler._create_ace_cap("B", 1, 100, RESIDUE[0], residue_lines=RESIDUE, all_coords=everything))
    assert set(cap) == {"CH3", "C", "O"}
    assert np.linalg.norm(cap["C"] - residue["N"]) == pytest.approx(1.329, abs=0.002)
    assert np.linalg.norm(cap["O"] - cap["C"]) == pytest.approx(1.231, abs=0.002)
    assert np.linalg.norm(cap["CH3"] - cap["C"]) == pytest.approx(1.525, abs=0.002)
    assert _angle(cap["C"], residue["N"], residue["CA"]) == pytest.approx(121.7, abs=0.2)
    assert abs(_torsion(cap["CH3"], cap["C"], residue["N"], residue["CA"])) == pytest.approx(180.0, abs=0.5)   # trans
    assert abs(_torsion(cap["O"], cap["C"], residue["N"], residue["CA"])) == pytest.approx(0.0, abs=0.5)
    # planar carbonyl: the three angles at C add up to 360
    at_c = _angle(residue["N"], cap["C"], cap["O"]) + _angle(residue["N"], cap["C"], cap["CH3"]) + _angle(cap["O"], cap["C"], cap["CH3"])
    assert at_c == pytest.approx(360.0, abs=0.5)


def test_ace_turns_away_from_whatever_is_in_its_way(handler):
    residue = _atoms(RESIDUE)
    base = np.array([_xyz(l) for l in RESIDUE])
    free = _atoms(handler._create_ace_cap("B", 1, 1, RESIDUE[0], residue_lines=RESIDUE, all_coords=base))
    # put an atom exactly where that cap's oxygen went: the cap must go elsewhere
    crowded = np.vstack([base, free["O"], free["CH3"]])
    moved = _atoms(handler._create_ace_cap("B", 1, 1, RESIDUE[0], residue_lines=RESIDUE, all_coords=crowded))
    assert np.linalg.norm(moved["O"] - free["O"]) > 1.0
    assert min(np.linalg.norm(crowded[5:] - v, axis=1).min() for v in moved.values()) > 1.5
    # and it is still the same peptide bond
    assert _angle(moved["C"], residue["N"], residue["CA"]) == pytest.approx(121.7, abs=0.2)


def test_ace_does_not_depend_on_which_way_the_file_faces(handler):
    """The old cap was laid along the file's x axis: turn the structure and it changed shape."""
    turn = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    def turned(line):
        x, y, z = _xyz(line) @ turn
        return f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
    rotated = [turned(l) for l in RESIDUE]
    a = _atoms(handler._create_ace_cap("B", 1, 1, RESIDUE[0], residue_lines=RESIDUE, all_coords=np.array([_xyz(l) for l in RESIDUE])))
    b = _atoms(handler._create_ace_cap("B", 1, 1, rotated[0], residue_lines=rotated, all_coords=np.array([_xyz(l) for l in rotated])))
    for name in a:
        assert np.allclose(a[name] @ turn, b[name], atol=0.01)


def test_nme_is_bonded_to_the_carbonyl_carbon_not_to_the_residues_last_atom(handler):
    residue = _atoms(RESIDUE)
    cap = _atoms(handler._create_nme_cap("B", 3, 200, RESIDUE[-1], residue_lines=RESIDUE))     # last atom is CB
    assert set(cap) == {"N", "H"}                  # tLEaP builds the methyl under its own name
    assert np.linalg.norm(cap["N"] - residue["C"]) == pytest.approx(1.329, abs=0.002)
    assert np.linalg.norm(cap["N"] - residue["CB"]) > 2.0
    assert _angle(residue["CA"], residue["C"], cap["N"]) == pytest.approx(116.2, abs=0.2)
    assert abs(_torsion(residue["O"], residue["C"], cap["N"], cap["H"])) == pytest.approx(180.0, abs=0.5)      # H trans to O
    assert np.linalg.norm(cap["H"] - cap["N"]) == pytest.approx(1.010, abs=0.002)


def test_a_residue_without_a_backbone_falls_back_and_says_so():
    out = io.StringIO()
    handler = CappingHandler(Console(file=out, width=200))
    only_n = RESIDUE[:1]
    cap = handler._create_ace_cap("B", 1, 1, only_n[0], residue_lines=only_n, all_coords=np.zeros((0, 3)))
    assert len(cap) == 3 and "no N, CA and C to build on" in out.getvalue()


def test_a_backbone_on_one_line_falls_back_instead_of_writing_nan(handler):
    """N, CA and C collinear define no plane: the placement would divide by zero."""
    straight = [f"ATOM  {i:5d}  {name:<3s} ALA A   2    {x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00\n"
                for i, (name, x) in enumerate((("N", 0.0), ("CA", 1.5), ("C", 3.0), ("O", 4.2)), 1)]
    for lines in (handler._create_ace_cap("A", 1, 1, straight[0], residue_lines=straight, all_coords=np.zeros((0, 3))),
                  handler._create_nme_cap("A", 3, 1, straight[-1], residue_lines=straight)):
        for line in lines:
            assert np.isfinite([float(line[30:38]), float(line[38:46]), float(line[46:54])]).all()
