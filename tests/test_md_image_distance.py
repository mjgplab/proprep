"""How close a solute comes to its own periodic images, for any Amber box.

Asked for after a membrane build in which the protein looked off-centre in its
patch. Where the solute sits in the box means nothing under periodic
boundaries; its distance to its own copies does, and it changes along a run.

The numbers are checked against two things that share no code with ProPrep's
k-d tree search: cpptraj's ``minimage`` (through pytraj), and a brute-force
search over two shells of neighbouring cells. Checked live when written, too:
every frame of a real truncated-octahedral run (chignolin, 300 frames) agreed
with cpptraj to 0.00000 A, and the monitor read a trajectory that sander was
still writing (4, 9, 13, ... frames on successive calls).
"""

import io
import types

import numpy as np
import pytest
from rich.console import Console

from proprep.md_prep import image_distance as imd

pt = pytest.importorskip("pytraj")

from proprep.md_prep import molecular_dynamics_manager as mdm_module  # noqa: E402
from proprep.md_prep.molecular_dynamics_manager import MolecularDynamicsManager as MDM  # noqa: E402
from proprep.md_prep.trajectory_analyzer import TrajectoryAnalyzer  # noqa: E402
from proprep.utils import prompts as prompts_module  # noqa: E402

TOP, NC = pt.datafiles.tz2_ortho_parm7, pt.datafiles.tz2_ortho_nc
OCT = imd.TRUNCATED_OCTAHEDRON_ANGLE
BOXES = {
    "rectangular": [31.0, 36.0, 40.0, 90.0, 90.0, 90.0],
    "truncated octahedron": [30.0, 30.0, 30.0, OCT, OCT, OCT],
    "triclinic": [30.0, 34.0, 38.0, 80.0, 100.0, 70.0],
}


def _brute_force(xyz, box, reach=2):
    """Every atom against every atom of every copy within ``reach`` cells. No tree, no half-space."""
    vectors, best = imd.lattice_vectors(box), np.inf
    for i in range(-reach, reach + 1):
        for j in range(-reach, reach + 1):
            for k in range(-reach, reach + 1):
                if i == j == k == 0:
                    continue
                copy = xyz + i * vectors[0] + j * vectors[1] + k * vectors[2]
                best = min(best, np.sqrt(((xyz[:, None, :] - copy[None, :, :]) ** 2).sum(-1)).min())
    return best


# ── The box ─────────────────────────────────────────────────────────────

@pytest.mark.parametrize("shape, box", BOXES.items())
def test_box_shape_is_named_from_its_angles_and_edges(shape, box):
    assert imd.box_shape(box) == shape


def test_a_trajectory_files_two_decimal_angle_is_still_a_truncated_octahedron():
    assert imd.box_shape([40.13, 40.13, 40.13, 109.47, 109.47, 109.47]) == "truncated octahedron"
    # equal angles but unequal edges is not one
    assert imd.box_shape([40.0, 41.0, 40.0, 109.47, 109.47, 109.47]) == "triclinic"


def test_lattice_vectors_have_the_box_lengths_and_angles():
    for box in BOXES.values():
        a, b, c = imd.lattice_vectors(box)
        assert np.allclose([np.linalg.norm(v) for v in (a, b, c)], box[:3])
        angle = lambda u, v: np.degrees(np.arccos(u @ v / np.linalg.norm(u) / np.linalg.norm(v)))
        assert np.allclose([angle(b, c), angle(a, c), angle(a, b)], box[3:], atol=1e-6)


@pytest.mark.parametrize("box", [None, [0.0, 0.0, 0.0, 90.0, 90.0, 90.0], [30.0, 30.0]])
def test_no_box_means_no_images(box):
    with pytest.raises(imd.NotPeriodic):
        imd.minimum_image_distance(np.zeros((2, 3)), box)


# ── The distance ────────────────────────────────────────────────────────

def test_a_single_atom_is_one_cell_edge_from_its_nearest_copy():
    """Known without any calculation: the shortest lattice vector."""
    atom = np.array([[3.0, -7.0, 11.0]])
    assert imd.minimum_image_distance(atom, BOXES["rectangular"])[0] == pytest.approx(31.0)
    assert imd.minimum_image_distance(atom, BOXES["truncated octahedron"])[0] == pytest.approx(30.0)


@pytest.mark.parametrize("shape", BOXES)
def test_agrees_with_a_brute_force_search_over_two_shells_of_cells(shape):
    rng = np.random.default_rng(7)
    xyz = rng.normal(scale=6.0, size=(60, 3))
    distance, first, second, shift = imd.minimum_image_distance(xyz, BOXES[shape])
    assert distance == pytest.approx(_brute_force(xyz, BOXES[shape]), abs=1e-9)
    # the pair and the cell it names really are that far apart
    copy = xyz[second] + np.asarray(shift, float) @ imd.lattice_vectors(BOXES[shape])
    assert np.linalg.norm(xyz[first] - copy) == pytest.approx(distance, abs=1e-9)


def test_where_the_solute_sits_in_the_box_changes_nothing():
    rng = np.random.default_rng(3)
    xyz = rng.normal(scale=5.0, size=(40, 3))
    for box in BOXES.values():
        here = imd.minimum_image_distance(xyz, box)[0]
        moved = imd.minimum_image_distance(xyz + np.array([11.0, -14.5, 9.25]), box)[0]
        assert here == pytest.approx(moved, abs=1e-9)


# ── Against cpptraj, on a real solute ───────────────────────────────────

def _analyzer(box=None):
    traj = pt.load(NC, top=TOP)
    if box is not None:
        traj.unitcells[:] = box
    return TrajectoryAnalyzer(TOP, traj_object=traj)


@pytest.mark.parametrize("shape, box", [("rectangular", None), ("truncated octahedron", BOXES["truncated octahedron"]),
                                        ("triclinic", BOXES["triclinic"])])
def test_every_frame_agrees_with_cpptraj_minimage(shape, box):
    analyzer = _analyzer(box)
    mask = analyzer.solute_mask()
    before = analyzer.traj.xyz.copy()
    result = analyzer.calculate_image_distance(mask)
    reference = np.asarray(pt.compute(f"minimage MI {mask} {mask}", analyzer.traj)["MI"])
    assert result["box_shape"] == shape and result["frames"] == list(range(10))
    assert np.allclose(result["distances"], reference, atol=1e-3)
    assert np.array_equal(before, analyzer.traj.xyz), "the analysis moved the shared trajectory"


def test_the_solute_is_found_from_the_topologys_molecules():
    """trpzip2: one 220-atom molecule, then water. No residue names involved."""
    assert _analyzer().solute_mask() == "@1-220"


def test_the_closest_pair_is_named_and_the_stride_is_honoured():
    analyzer = _analyzer()
    result = analyzer.calculate_image_distance("@1-220", stride=3)
    assert result["frames"] == [0, 3, 6, 9] and len(result["distances"]) == 4
    closest = result["closest"]
    assert closest["distance"] == min(result["distances"]) and closest["frame"] in result["frames"]
    for label in (closest["atom"], closest["image_atom"]):
        name, number, atom = label.split()
        assert 1 <= int(number) <= 13 and name.isalpha() and atom


def test_image_distance_does_not_depend_on_what_ran_before():
    """rmsd and pca superpose in place without rotating the box; this must not see that."""
    analyzer = _analyzer()
    fresh = analyzer.calculate_image_distance("@1-220")["distances"]
    analyzer.calculate_rmsd("@CA", 0)
    analyzer.calculate_pca("@CA", 2)
    assert np.allclose(fresh, analyzer.calculate_image_distance("@1-220", label="again")["distances"])


def test_bad_input_is_refused():
    analyzer = _analyzer()
    with pytest.raises(ValueError, match="selects no atom"):
        analyzer.calculate_image_distance(":NOPE")
    with pytest.raises(ValueError, match="stride"):
        analyzer.calculate_image_distance("@1-220", stride=0)


# ── The cutoff the run used ─────────────────────────────────────────────

def test_cutoff_is_read_from_the_mdout_in_preference_to_the_input(tmp_path):
    (tmp_path / "simulation.mdin").write_text(" &cntrl\n  cut=9.0,     ! Nonbonded cutoff\n /\n")
    assert imd.read_cutoff(tmp_path) == (9.0, "simulation.mdin")
    (tmp_path / "prod.mdout").write_text("     dielc   =   1.00000, cut     =  10.00000, intdiel =   1.00000\n")
    assert imd.read_cutoff(tmp_path) == (10.0, "prod.mdout")


def test_a_commented_out_cutoff_is_not_read(tmp_path):
    (tmp_path / "simulation.mdin").write_text(" &cntrl\n  ! cut=12.0 was tried first\n  ntb=1,\n /\n")
    assert imd.read_cutoff(tmp_path) is None
    assert imd.read_cutoff(tmp_path / "missing") is None


# ── In the MD Manager ───────────────────────────────────────────────────

def _drive(monkeypatch, call, answers):
    def prompt(processor, prompt=None, *args, **kwargs):
        for fragment, answer in answers:
            if fragment in str(prompt).lower():
                return answer
        return kwargs["default"]

    monkeypatch.setattr(mdm_module, "prompt_with_context", prompt)
    monkeypatch.setattr(prompts_module, "prompt_with_context", prompt)
    console = Console(width=140, record=True, file=io.StringIO(), force_terminal=False)
    manager = object.__new__(MDM)
    manager.processor = types.SimpleNamespace(console=console)
    call(manager)
    return console.export_text()


def test_a_solute_inside_the_cutoff_of_its_image_is_said_plainly(monkeypatch):
    """A 30 A truncated octahedron is too small for trpzip2: about 11 A to its copy."""
    analyzer = _analyzer(BOXES["truncated octahedron"])
    out = _drive(monkeypatch, lambda m: m._analyze_image_distance(analyzer),
                 [("select the solute", "2"), ("cutoff", "12")])
    assert "truncated octahedron" in out and "10 of 10" in out
    assert "interacts directly with itself" in out


def test_a_solute_clear_of_its_image_is_said_plainly_too(monkeypatch, tmp_path):
    analyzer = _analyzer()
    (tmp_path / "prod.mdout").write_text("     dielc   =   1.00000, cut     =  10.00000, intdiel =   1.00000\n")
    out = _drive(monkeypatch, lambda m: m._analyze_image_distance(analyzer, sim_dir=tmp_path), [])
    assert "read from prod.mdout" in out and "rectangular" in out
    assert "never comes within the 10 Å cutoff" in out and "0 of 10" in out


def test_without_a_cutoff_beside_the_run_ambers_default_is_named_as_such(monkeypatch, tmp_path):
    out = _drive(monkeypatch, lambda m: m._analyze_image_distance(_analyzer(), sim_dir=tmp_path), [])
    assert "Amber's own default" in out and "8 Å cutoff" in out


def test_the_monitor_says_why_when_there_is_nothing_to_read(monkeypatch, tmp_path):
    out = _drive(monkeypatch, lambda m: m._monitor_image_distance(tmp_path), [])
    assert "No topology" in out
    (tmp_path / "system.prmtop").write_text("x")
    out = _drive(monkeypatch, lambda m: m._monitor_image_distance(tmp_path), [])
    assert "No trajectory (.nc)" in out and "minimization writes no trajectory" in out


def test_the_monitor_reads_a_runs_own_files(monkeypatch, tmp_path):
    import shutil
    shutil.copy(TOP, tmp_path / "system.prmtop")
    shutil.copy(NC, tmp_path / "prod.nc")
    (tmp_path / "prod.mdout").write_text("     dielc   =   1.00000, cut     =  10.00000, intdiel =   1.00000\n")
    out = _drive(monkeypatch, lambda m: m._monitor_image_distance(tmp_path), [("select the solute", "2")])
    assert "10 frames written so far" in out and "read from prod.mdout" in out
    assert "Closest approach" in out and "Error" not in out
