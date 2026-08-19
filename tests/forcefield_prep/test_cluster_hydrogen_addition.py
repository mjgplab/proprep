"""
An inorganic cluster can be given a hydrogen, e.g. the Mo-OH of a Mo cofactor.

Hydrogen addition covers protein (category A) and organic residues (category
B); category F, the pure inorganic clusters, was never offered to it, and
`reduce` has no chemistry for a Mo-S-O core in any case. So a molybdenum
cofactor whose resting state is Mo(=O)(=S)(OH) reached the QM model as a bare
oxo, with the wrong electron count and the wrong charge, and there was no way to
fix it.

Editing the generated .gjf by hand is not the alternative: the Gaussian input
and the model PDB are matched by index, so an atom added to only one of them
shifts Seminario's indices and leaves the deposited residue template without the
hydrogen its charges were fitted with. The hydrogen has to enter the structure
before the models are built, which is what this step does.
"""

import numpy as np
import pytest

from proprep.forcefield_prep.hydrogen_editor import HydrogenEditor
from proprep.forcefield_prep.structure_preprocessor import StructurePreprocessor


# --------------------------------------------------------------------------- #
# placement geometry
# --------------------------------------------------------------------------- #

def _angle(a, b, c):
    """Angle a-b-c in degrees."""
    v1, v2 = np.array(a) - np.array(b), np.array(c) - np.array(b)
    cos = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))


def test_single_neighbour_hydrogen_is_bent_not_linear():
    """A hydroxo H opposite its only bond would be a linear Mo-O-H."""
    o = np.array([0.0, 0.0, 0.0])
    mo = np.array([1.9, 0.0, 0.0])
    bonds = [(mo - o) / np.linalg.norm(mo - o)]

    d = HydrogenEditor._place_hydrogen(o, bonds, np.array([mo]), 0.96, 109.5)
    h = o + d * 0.96

    assert _angle(mo, o, h) == pytest.approx(109.5, abs=0.5)
    assert np.linalg.norm(h - o) == pytest.approx(0.96, abs=1e-6)


@pytest.mark.parametrize("bend", [104.5, 109.5, 120.0])
def test_the_opening_angle_is_honoured(bend):
    o = np.array([0.0, 0.0, 0.0])
    mo = np.array([0.0, 2.1, 0.0])
    bonds = [(mo - o) / np.linalg.norm(mo - o)]

    d = HydrogenEditor._place_hydrogen(o, bonds, np.array([mo]), 0.96, bend)

    assert _angle(mo, o, o + d * 0.96) == pytest.approx(bend, abs=0.5)


def test_rotation_avoids_the_other_ligands():
    """The one free degree of freedom is spent on clearance, not arbitrarily.

    The bonded neighbour is excluded from the clearance set on purpose: it is
    equidistant from every azimuth, so including it makes the minimum identical
    all the way round and the scan silently degenerates to the first angle.
    """
    o = np.array([0.0, 0.0, 0.0])
    mo = np.array([1.9, 0.0, 0.0])
    bonds = [(mo - o) / np.linalg.norm(mo - o)]
    # Two further Mo ligands crowd the +y side. Mo itself is NOT in this set.
    crowd = np.array([[1.9, 1.6, 0.0], [2.4, 2.2, 0.0]])

    d = HydrogenEditor._place_hydrogen(o, bonds, crowd, 0.96, 109.5)
    h = o + d * 0.96

    assert h[1] < 0.0, f"H should swing away from the +y ligands, got {h.round(3)}"
    # And it beats the naive first-azimuth choice it would otherwise take.
    naive = o + HydrogenEditor._place_hydrogen(
        o, bonds, np.zeros((0, 3)), 0.96, 109.5) * 0.96
    assert (float(np.min(np.linalg.norm(crowd - h, axis=1)))
            >= float(np.min(np.linalg.norm(crowd - naive, axis=1))))


def test_two_neighbours_keep_completing_the_coordination():
    """The determined case must be unchanged."""
    b = [np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])]

    d = HydrogenEditor._place_hydrogen(np.zeros(3), b, np.zeros((0, 3)), 1.01, 109.5)

    assert d == pytest.approx([-0.7071, -0.7071, 0.0], abs=1e-3)


def test_cancelling_neighbours_are_reported_not_guessed():
    b = [np.array([1.0, 0.0, 0.0]), np.array([-1.0, 0.0, 0.0])]

    assert HydrogenEditor._place_hydrogen(
        np.zeros(3), b, np.zeros((0, 3)), 1.01, 109.5) is None


# --------------------------------------------------------------------------- #
# merging the hydrogen into the structure
# --------------------------------------------------------------------------- #

from proprep.utils.pdb_format import atom_name_field  # noqa: E402


def _pdb_line(serial, name, element, resname, chain, resid, xyz, record="HETATM"):
    """Build a column-correct record (hand-counting these is how names drift)."""
    return (
        f"{record:<6.6s}{serial:5d} {atom_name_field(name, element)} "
        f"{resname:>3.3s} {chain}{resid:>4d}    "
        f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
        f"  1.00  0.00          {element:>2s}"
    )


MOCO_PDB = [
    _pdb_line(1, "N", "N", "ALA", "A", 10, (10.0, 10.0, 10.0), "ATOM"),
    _pdb_line(2, "CA", "C", "ALA", "A", 10, (11.0, 10.0, 10.0), "ATOM"),
    _pdb_line(3, "MO", "MO", "MOS", "A", 3004, (0.0, 0.0, 0.0)),
    _pdb_line(4, "S", "S", "MOS", "A", 3004, (1.9, 0.0, 0.0)),
    _pdb_line(5, "O1", "O", "MOS", "A", 3004, (-1.7, 0.0, 0.0)),
    _pdb_line(6, "O2", "O", "MOS", "A", 3004, (0.0, -1.7, 0.0)),
    _pdb_line(7, "C1", "C", "FAD", "A", 3006, (20.0, 20.0, 20.0)),
    "END",
]


def _preprocessor(tmp_path):
    sp = StructurePreprocessor.__new__(StructurePreprocessor)
    pdb = tmp_path / "structure.pdb"
    pdb.write_text("\n".join(MOCO_PDB) + "\n")
    sp._pdb_file = str(pdb)
    return sp, pdb


def _residue_pdb_with_h(tmp_path):
    """What HydrogenEditor leaves behind: the residue plus a new H."""
    p = tmp_path / "cluster_MOS_A_3004.pdb"
    p.write_text("\n".join([
        _pdb_line(1, "MO", "MO", "MOS", "A", 3004, (0.0, 0.0, 0.0)),
        _pdb_line(2, "S", "S", "MOS", "A", 3004, (1.9, 0.0, 0.0)),
        _pdb_line(3, "O1", "O", "MOS", "A", 3004, (-1.7, 0.0, 0.0)),
        _pdb_line(4, "O2", "O", "MOS", "A", 3004, (0.0, -1.7, 0.0)),
        _pdb_line(5, "H1", "H", "MOS", "A", 3004, (-2.02, 0.905, 0.0)),
        "END",
    ]) + "\n")
    return p


def test_hydrogen_lands_inside_the_cluster_residue(tmp_path):
    sp, pdb = _preprocessor(tmp_path)
    edited = _residue_pdb_with_h(tmp_path)

    added = sp._merge_cluster_hydrogens(str(edited), "A", 3004, "MOS")

    assert added == 1
    lines = [l for l in pdb.read_text().splitlines() if l.startswith(("ATOM", "HETATM"))]
    mos = [l for l in lines if l[17:20].strip() == "MOS"]
    assert len(mos) == 5, "the H should join the cluster residue"
    assert mos[-1][12:16].strip() == "H1"
    # Contiguous: the residue's atoms are not split by the following FAD.
    idx = [i for i, l in enumerate(lines) if l[17:20].strip() == "MOS"]
    assert idx == list(range(idx[0], idx[0] + 5))


def test_serials_are_renumbered_without_gaps(tmp_path):
    sp, pdb = _preprocessor(tmp_path)

    sp._merge_cluster_hydrogens(str(_residue_pdb_with_h(tmp_path)), "A", 3004, "MOS")

    lines = [l for l in pdb.read_text().splitlines() if l.startswith(("ATOM", "HETATM"))]
    serials = [int(l[6:11]) for l in lines]
    assert serials == list(range(1, len(lines) + 1))


def test_the_merged_hydrogen_reads_back_as_hydrogen(tmp_path):
    """Columns matter: a name in the wrong column changes the element."""
    sp, pdb = _preprocessor(tmp_path)
    sp._merge_cluster_hydrogens(str(_residue_pdb_with_h(tmp_path)), "A", 3004, "MOS")

    h = [l for l in pdb.read_text().splitlines()
         if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() == "MOS"][-1]

    assert h[12:16] == " H1 ", repr(h[12:16])   # one-letter element, indented
    assert h[76:78] == " H"
    assert h[21] == "A"
    assert h[22:26].strip() == "3004"


def test_coordinates_survive_the_merge(tmp_path):
    sp, pdb = _preprocessor(tmp_path)
    sp._merge_cluster_hydrogens(str(_residue_pdb_with_h(tmp_path)), "A", 3004, "MOS")

    h = [l for l in pdb.read_text().splitlines()
         if l.startswith(("ATOM", "HETATM")) and l[12:16].strip() == "H1"][0]

    assert float(h[30:38]) == pytest.approx(-2.020)
    assert float(h[38:46]) == pytest.approx(0.905)


def test_no_hydrogens_is_a_no_op(tmp_path):
    sp, pdb = _preprocessor(tmp_path)
    before = pdb.read_text()
    plain = tmp_path / "plain.pdb"
    plain.write_text("HETATM    1 MO    MOS A3004       0.000   0.000   0.000  1.00  0.00          MO\n")

    assert sp._merge_cluster_hydrogens(str(plain), "A", 3004, "MOS") == 0
    assert pdb.read_text() == before


# --------------------------------------------------------------------------- #
# the hydrogen must not be typed as a metal ligand
# --------------------------------------------------------------------------- #

class _Atom:
    def __init__(self, name, element, coords, chain="A", resid=3004, resname="MOS"):
        self.atom_name, self.element, self.coords = name, element, coords
        self.chain, self.resid, self.resname = chain, resid, resname


class _Site:
    def __init__(self, atoms):
        self.atoms = atoms


def test_a_hydroxo_hydrogen_is_not_given_a_metal_ligand_type():
    """Y* types are for atoms coordinating the metal; this H bonds an oxygen."""
    from proprep.forcefield_prep.metal_site_parameterizer import MetalSiteWorkflowManager

    mo = _Atom("MO", "MO", (0.0, 0.0, 0.0))
    atoms = [
        mo,
        _Atom("S", "S", (1.9, 0.0, 0.0)),
        _Atom("O1", "O", (-1.7, 0.0, 0.0)),
        _Atom("O2", "O", (0.0, -1.7, 0.0)),
        _Atom("H1", "H", (-2.02, 0.905, 0.0)),
    ]
    mgr = MetalSiteWorkflowManager.__new__(MetalSiteWorkflowManager)

    ligands = mgr._cluster_internal_ligand_coords(_Site(atoms), {mo.coords})

    assert (-1.7, 0.0, 0.0) in ligands, "the oxo/hydroxo O still coordinates Mo"
    assert (1.9, 0.0, 0.0) in ligands, "the sulfido still coordinates Mo"
    assert (-2.02, 0.905, 0.0) not in ligands, "the H does not coordinate the metal"
    assert mo.coords not in ligands, "metals stay M*"
