"""
Orient a protein in the membrane frame without losing any of it.

packmol-memgen can run PPM3 itself (``--ppm``), but it then packs PPM3's
*output file*, and PPM3 writes only the residues in its own library, without
hydrogens. On a ProPrep-prepared 6R2Q that dropped all 20 hemes (HCO and their
PRD propionates), every heme-ligating HIO and CYO, the HIE, HIP, CYX, LYN and
GLH residues, the waters and the calcium ions: 24,057 atoms in, 10,769 out.
Every ProPrep structure carries such residue names, so that route is lossy
here as a rule, not as an exception.

Here PPM3 is asked for one thing only, where the membrane is. Its output gives
a rigid transform (the atoms it did keep, matched to the same atoms of the
input), that transform is applied to the complete input structure, and
packmol-memgen is given the result as ``--preoriented``. Nothing is dropped
because nothing but coordinates is taken from PPM3.

PPM is the method behind the OPM database. For 6R2Q this reproduces OPM's
deposited orientation to within 1 A on every chain; MEMEMBED, in either its
alpha-helical or its beta-barrel mode, leaves that complex under the membrane.
"""

import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

# The fit is between two copies of the same atoms, so it is exact up to the
# three decimals of a PDB file. Anything worse means the wrong atoms were
# paired, and the transform must not be used.
MAX_FIT_RMSD = 0.01  # A

AtomKey = Tuple[str, int, str, str]      # chain, residue number, insertion code, atom name


class OrientationError(RuntimeError):
    """The protein could not be oriented; the message says why."""


def _atoms(pdb_path) -> Dict[AtomKey, np.ndarray]:
    """Heavy atoms of a PDB file by identity; an atom named twice keeps its first record."""
    found: Dict[AtomKey, np.ndarray] = {}
    with open(pdb_path, errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")) or line[17:20] == "DUM":
                continue
            name = line[12:16].strip()
            if name.startswith("H") or line[76:78].strip() == "H":
                continue
            try:
                key = (line[21], int(line[22:26]), line[26], name)
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue
            found.setdefault(key, xyz)
    return found


def rigid_transform(source: np.ndarray, target: np.ndarray):
    """Rotation and translation taking ``source`` onto ``target`` (Kabsch), and the RMSD left.

    Returns ``(rotation, translation, rmsd)`` for row vectors:
    ``fitted = points @ rotation + translation``.
    """
    source_centre, target_centre = source.mean(axis=0), target.mean(axis=0)
    u, _, vt = np.linalg.svd((source - source_centre).T @ (target - target_centre))
    # A proper rotation: never a mirror image of the protein.
    rotation = u @ np.diag([1.0, 1.0, np.sign(np.linalg.det(u @ vt))]) @ vt
    translation = target_centre - source_centre @ rotation
    fitted = source @ rotation + translation
    rmsd = float(np.sqrt(((fitted - target) ** 2).sum(axis=1).mean()))
    return rotation, translation, rmsd


def transform_from_oriented_copy(original_pdb, oriented_pdb):
    """The rigid transform that carries ``original_pdb`` into the frame of ``oriented_pdb``.

    ``oriented_pdb`` is a copy of the same structure, moved, and possibly
    missing atoms. Raises OrientationError when too few atoms are shared or
    the copy is not a rigid image of the original.
    """
    original, oriented = _atoms(original_pdb), _atoms(oriented_pdb)
    shared = sorted(set(original) & set(oriented))
    if len(shared) < 3:
        raise OrientationError(
            f"Only {len(shared)} atoms of {Path(original_pdb).name} are found again in "
            f"{Path(oriented_pdb).name}; at least 3 are needed to place it.")
    rotation, translation, rmsd = rigid_transform(
        np.array([original[k] for k in shared]), np.array([oriented[k] for k in shared]))
    if rmsd > MAX_FIT_RMSD:
        raise OrientationError(
            f"{Path(oriented_pdb).name} is not a rigid copy of {Path(original_pdb).name}: "
            f"{len(shared)} shared atoms fit with an RMSD of {rmsd:.3f} A "
            f"(limit {MAX_FIT_RMSD} A). The orientation was not applied.")
    return rotation, translation, rmsd, len(shared)


def apply_transform(pdb_in, pdb_out, rotation: np.ndarray, translation: np.ndarray) -> int:
    """Write ``pdb_in`` with every atom moved; all other text is copied as it is."""
    moved = 0
    with open(pdb_in, errors="replace") as source, open(pdb_out, "w") as out:
        for line in source:
            if line.startswith(("ATOM", "HETATM")):
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                x, y, z = xyz @ rotation + translation
                line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
                moved += 1
            out.write(line)
    return moved


def find_ppm3() -> Optional[str]:
    """PPM3's executable (``immers``), as AmberTools installs it."""
    found = shutil.which("immers")
    if found:
        return found
    amberhome = os.environ.get("AMBERHOME")
    if amberhome and os.path.isfile(os.path.join(amberhome, "bin", "immers")):
        return os.path.join(amberhome, "bin", "immers")
    return None


def _ppm3_residue_library() -> Optional[str]:
    try:
        import packmol_memgen
    except ImportError:
        return None
    path = os.path.join(os.path.dirname(packmol_memgen.__file__), "lib", "ppm3", "res.lib")
    return path if os.path.isfile(path) else None


def run_ppm3(protein_pdb, scratch_dir, n_ter: str = "in") -> Tuple[Path, Path]:
    """Run PPM3 on a copy of the protein; return its output PDB and its log.

    The input deck is the one packmol-memgen writes for a single planar
    bilayer. PPM3 reads fixed file names from its working directory.
    """
    executable, library = find_ppm3(), _ppm3_residue_library()
    if executable is None:
        raise OrientationError("PPM3 (immers) was not found on PATH or in $AMBERHOME/bin.")
    if library is None:
        raise OrientationError("PPM3's residue library (packmol_memgen/lib/ppm3/res.lib) was not found.")

    scratch = Path(scratch_dir)
    scratch.mkdir(parents=True, exist_ok=True)
    shutil.copy(protein_pdb, scratch / "ppm3in.pdb")
    shutil.copy(library, scratch / "res.lib")
    with open(scratch / "ppm3in.pdb", errors="replace") as handle:
        chains = sorted({line[21] for line in handle if line.startswith("ATOM")})
    (scratch / "ppm.inp").write_text(
        "2\nno\nppm3in.pdb\n1\n    \nplanar\n" + f"{n_ter}\n" + ",".join(chains) + "\n")

    log = scratch / "ppm3.log"
    with open(scratch / "ppm.inp") as deck, open(log, "w") as out:
        subprocess.run([executable], stdin=deck, stdout=out, stderr=subprocess.STDOUT, cwd=scratch)
    output = scratch / "ppm3inout.pdb"
    if not output.exists() or output.stat().st_size == 0:
        raise OrientationError(f"PPM3 wrote no oriented structure. Its log is {log}.")
    return output, log


def ppm3_summary(log_path) -> Optional[str]:
    """PPM3's own result line (transfer energy, hydrophobic thickness, tilt)."""
    try:
        lines = Path(log_path).read_text(errors="replace").splitlines()
    except OSError:
        return None
    for line in reversed(lines):
        if "emin=" in line and "thickn=" in line:
            return " ".join(line.split())
    return None


def orient_with_ppm3(protein_pdb, oriented_pdb, scratch_dir, n_ter: str = "in") -> Dict[str, object]:
    """Write ``oriented_pdb``: all of ``protein_pdb``, in the membrane frame PPM3 finds.

    Returns what was done, for the caller to report.
    """
    ppm3_output, log = run_ppm3(protein_pdb, scratch_dir, n_ter=n_ter)
    rotation, translation, rmsd, shared = transform_from_oriented_copy(protein_pdb, ppm3_output)
    moved = apply_transform(protein_pdb, oriented_pdb, rotation, translation)
    return {"atoms_moved": moved, "atoms_matched": shared, "fit_rmsd": rmsd,
            "ppm3_result": ppm3_summary(log), "ppm3_log": str(log)}
