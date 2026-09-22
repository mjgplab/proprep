"""PPM3 orients the protein; ProPrep keeps the protein.

Seen on 6R2Q. MEMEMBED (alpha-helical or beta-barrel mode) left the
transmembrane barrel complex under the membrane. PPM3 places it as OPM does,
but packmol-memgen's own ``--ppm`` packs PPM3's output file, which holds only
residues PPM3 knows and no hydrogens: 24,057 atoms became 10,769, without the
20 hemes (HCO, PRD), HIO, CYO, HIE, HIP, CYX, LYN, GLH, waters or ions. It
also fails outright when ``-p`` has a directory part ("./protein_with_h.pdb").

ProPrep therefore takes only a rigid transform from PPM3 and passes the whole
structure as ``--preoriented``. Checked live on 6R2Q with real PPM3 and
packmol-memgen: 24,057 atoms in PROT0.pdb, 20 Fe, protein volume in both
leaflets (30,191 and 24,985 A^3), OPM's orientation to within 1 A per chain.
"""

import io
import os
import types

import numpy as np
import pytest
from rich.console import Console

from proprep.membrane_prep import membrane_builder as mb
from proprep.membrane_prep import orientation, packmol_runner
from proprep.membrane_prep.orientation import (
    OrientationError,
    apply_transform,
    rigid_transform,
    transform_from_oriented_copy,
)

# A rotation about an oblique axis plus a shift: nothing a bug could hit by luck.
ANGLE = np.radians(218.4)
AXIS = np.array([1.0, 2.0, -0.5]) / np.linalg.norm([1.0, 2.0, -0.5])
K = np.array([[0, -AXIS[2], AXIS[1]], [AXIS[2], 0, -AXIS[0]], [-AXIS[1], AXIS[0], 0]])
ROTATION = np.eye(3) + np.sin(ANGLE) * K + (1 - np.cos(ANGLE)) * K @ K
SHIFT = np.array([12.5, -40.25, 7.0])

ATOMS = [  # name, residue name, residue number, element, x, y, z
    ("N", "CYX", 1, "N", 66.740, 9.069, 4.563), ("CA", "CYX", 1, "C", 68.175, 9.270, 4.711),
    ("HA", "CYX", 1, "H", 68.601, 9.301, 3.708), ("SG", "CYX", 1, "S", 70.176, 10.988, 5.768),
    ("N", "HIO", 2, "N", 71.755, 9.256, 0.938), ("CA", "HIO", 2, "C", 70.793, 10.165, 1.533),
    ("NE2", "HIO", 2, "N", 71.073, 10.260, 3.037), ("FE", "HCO", 3, "FE", 23.639, 107.015, 3.850),
    ("NA", "HCO", 3, "N", 25.100, 106.000, 4.900),
]


def _line(i, name, resname, resid, element, xyz):
    return (f"ATOM  {i:5d} {name:<4s} {resname:>3s}  {resid:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {element:>2s}\n")


def _write(path, keep=lambda atom: True, move=lambda xyz: xyz):
    rows = [a for a in ATOMS if keep(a)]
    path.write_text("REMARK made for a test\n" + "".join(
        _line(i, n, rn, ri, el, move(np.array([x, y, z]))) for i, (n, rn, ri, el, x, y, z) in enumerate(rows, 1))
        + "END\n")
    return path


def test_rigid_transform_recovers_a_known_move_and_is_never_a_mirror_image():
    points = np.array([[a[4], a[5], a[6]] for a in ATOMS])
    rotation, translation, rmsd = rigid_transform(points, points @ ROTATION + SHIFT)
    assert np.allclose(rotation, ROTATION, atol=1e-9) and np.allclose(translation, SHIFT, atol=1e-7)
    assert rmsd < 1e-9 and np.isclose(np.linalg.det(rotation), 1.0)
    # Even when the target IS a mirror image, the answer stays a proper rotation.
    mirrored = points * np.array([1.0, 1.0, -1.0])
    assert np.isclose(np.linalg.det(rigid_transform(points, mirrored)[0]), 1.0)


def test_whole_structure_is_moved_from_the_atoms_the_tool_kept(tmp_path):
    """The oriented copy has lost its hydrogens, its HIO and its heme, as PPM3's does."""
    original = _write(tmp_path / "protein.pdb")
    oriented_copy = _write(tmp_path / "ppm3out.pdb",
                           keep=lambda a: a[1] == "CYX" and a[3] != "H",
                           move=lambda xyz: xyz @ ROTATION + SHIFT)
    rotation, translation, rmsd, shared = transform_from_oriented_copy(original, oriented_copy)
    assert shared == 3 and rmsd < 0.01
    moved = apply_transform(original, tmp_path / "out.pdb", rotation, translation)

    before, after = original.read_text().splitlines(), (tmp_path / "out.pdb").read_text().splitlines()
    assert moved == len(ATOMS) and len(before) == len(after)
    for old, new in zip(before, after):
        assert old[:30] == new[:30] and old[54:] == new[54:]          # only coordinates change
    fe = next(l for l in after if l[12:16].strip() == "FE")
    expected = np.array([23.639, 107.015, 3.850]) @ ROTATION + SHIFT
    # Three atoms within 3 A of each other, rounded to PDB's three decimals, carry
    # the fit out to an iron 100 A away: a long lever arm, hence 0.05 A here. A real
    # run fits thousands of atoms across the whole protein (6R2Q: 12,374).
    assert np.allclose([float(fe[30:38]), float(fe[38:46]), float(fe[46:54])], expected, atol=0.05)


def test_numbering_and_chain_names_do_not_matter_only_that_both_files_agree(tmp_path):
    """The structure is matched with the orientation tool's copy of ITSELF, so whatever a
    repair did to its chain names and residue numbers is the same in both."""
    def renumbered(atom_line_args):
        i, n, rn, ri, el, xyz = atom_line_args
        return (f"ATOM  {i:5d} {n:<4s} {rn:>3s} Q{ri + 7000:4d}    "
                f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {el:>2s}\n")
    def write(path, keep, move):
        rows = [a for a in ATOMS if keep(a)]
        path.write_text("".join(renumbered((i, n, rn, ri, el, move(np.array([x, y, z]))))
                                for i, (n, rn, ri, el, x, y, z) in enumerate(rows, 1)))
        return path
    original = write(tmp_path / "repaired.pdb", lambda a: True, lambda xyz: xyz)
    copy = write(tmp_path / "tool_output.pdb", lambda a: a[1] == "CYX" and a[3] != "H",
                 lambda xyz: xyz @ ROTATION + SHIFT)
    rotation, translation, rmsd, shared = transform_from_oriented_copy(original, copy)
    assert shared == 3 and rmsd < 0.01 and np.allclose(rotation, ROTATION, atol=1e-3)


def test_a_copy_that_is_not_rigid_is_refused(tmp_path):
    original = _write(tmp_path / "protein.pdb")
    stretched = _write(tmp_path / "bad.pdb", move=lambda xyz: xyz * np.array([1.0, 1.0, 1.3]))
    with pytest.raises(OrientationError, match="not a rigid copy"):
        transform_from_oriented_copy(original, stretched)


def test_too_few_shared_atoms_is_refused(tmp_path):
    original = _write(tmp_path / "protein.pdb")
    sparse = _write(tmp_path / "sparse.pdb", keep=lambda a: a[0] == "SG")
    with pytest.raises(OrientationError, match="at least 3"):
        transform_from_oriented_copy(original, sparse)


# ── The build ───────────────────────────────────────────────────────────

def _build(monkeypatch, tmp_path, method):
    pdb = _write(tmp_path / "protein_with_h.pdb")
    module = mb.MembraneBuilderModule()
    module.processor = types.SimpleNamespace(console=Console(file=io.StringIO(), width=300))
    module.config.protein_pdb = str(pdb)
    module.config.skip_protonation = True
    module.config.orientation_method = method
    calls, oriented = [], []

    def fake_orient(protein_pdb, oriented_pdb, scratch_dir, n_ter="in"):
        oriented.append((protein_pdb, oriented_pdb))
        return {"atoms_moved": 9, "atoms_matched": 3, "fit_rmsd": 0.0005,
                "ppm3_result": "emin= -74.9 thickn= 22.9", "ppm3_log": "x"}

    def run_packmol(args, work_dir, console, **kwargs):
        calls.append(list(args))
        return packmol_runner.PackmolResult(success=False, error_message="stop here")

    monkeypatch.setattr(mb, "confirm_with_context", lambda *a, **k: True)
    monkeypatch.setattr(module, "_show_review", lambda workspace: None)
    # tLEaP's pass returns the path the builder itself makes: "<work_dir>/protein_with_h.pdb".
    monkeypatch.setattr(module, "_run_pre_tleap_hydrogen_pass", lambda ws, wd: os.path.join(wd, pdb.name))
    monkeypatch.setattr(orientation, "orient_with_ppm3", fake_orient)
    monkeypatch.setattr(packmol_runner, "find_packmol_memgen", lambda: "/x/packmol-memgen")
    monkeypatch.setattr(packmol_runner, "run_packmol_memgen", run_packmol)
    module._run_build(workspace=None)
    return calls[0], oriented, module


def test_ppm3_is_run_by_proprep_and_passed_as_preoriented(monkeypatch, tmp_path):
    args, oriented, module = _build(monkeypatch, tmp_path, "ppm3")
    assert len(oriented) == 1 and oriented[0][1].endswith("protein_oriented.pdb")
    assert "--preoriented" in args and "--ppm" not in args
    assert args[args.index("-p") + 1] == "protein_oriented.pdb"
    # The user's choice is untouched: the next build orients again.
    assert module.config.preoriented is False and module.config.orientation_method == "ppm3"


def test_packmol_memgen_is_given_a_bare_file_name(monkeypatch, tmp_path):
    """A directory part in -p breaks packmol-memgen's scratch-folder names."""
    args, oriented, _ = _build(monkeypatch, tmp_path, "memembed")
    assert oriented == []
    assert args[args.index("-p") + 1] == "protein_with_h.pdb"
    assert "--preoriented" not in args
