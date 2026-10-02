"""
The search box, the docking settings, and one docking run with AutoDock Vina.

Three scoring functions: Vina and Vinardo (maps computed by Vina itself), and
AutoDock4 (maps computed by autogrid4, then read by Vina). Vina ignores partial
charges; AD4's electrostatic and desolvation terms use them, and AD4 was
calibrated with Gasteiger charges on receptor and ligand, which is what the
Meeko preparation writes.

Every setting that changes the result is a field of ``DockingSettings`` with no
default here; the caller shows Vina's defaults and records the values used. The
seed must be given and must not be 0: Vina reads 0 as "pick a random seed",
and a random seed makes the run impossible to replay.

Vina's own console output (for example its warning when the search space is
very large) is captured and returned in the run's notes instead of being
restated here.
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import time
from contextlib import contextmanager
from dataclasses import asdict, dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

SCORING_FUNCTIONS = ("vina", "vinardo", "ad4")

_ENERGY_COLUMNS = {
    "vina": ("total", "inter", "intra", "torsions", "intra best pose"),
    "vinardo": ("total", "inter", "intra", "torsions", "intra best pose"),
    "ad4": ("total", "inter", "intra", "torsions", "-intra"),
}


@dataclass(frozen=True)
class SearchBox:
    center: Tuple[float, float, float]
    size: Tuple[float, float, float]           # edge lengths, Angstrom
    source: str                                # how it was defined, for the record

    def to_dict(self) -> dict:
        return {"center": list(self.center), "size": list(self.size), "source": self.source}


def box_around(coordinates: np.ndarray, padding: float, source: str) -> SearchBox:
    """The axis-aligned box enclosing ``coordinates`` plus ``padding`` Angstrom on every side."""
    coordinates = np.asarray(coordinates, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3 or len(coordinates) == 0:
        raise ValueError("A box needs at least one xyz coordinate.")
    if padding < 0:
        raise ValueError("Padding cannot be negative.")
    low, high = coordinates.min(axis=0), coordinates.max(axis=0)
    center = (low + high) / 2.0
    size = (high - low) + 2.0 * padding
    return SearchBox(tuple(round(float(x), 3) for x in center), tuple(round(float(x), 3) for x in size),
                     f"{source}, padding {padding:g} A")


@dataclass(frozen=True)
class DockingSettings:
    scoring: str
    seed: int
    exhaustiveness: int
    n_poses: int               # poses generated and returned
    energy_range: float        # kcal/mol above the best pose
    min_rmsd: float            # Angstrom between returned poses
    max_evals: int             # 0: Vina's own heuristic
    spacing: float             # map grid spacing, Angstrom
    cpu: int                   # 0: all available cores

    def validate(self) -> None:
        if self.scoring not in SCORING_FUNCTIONS:
            raise ValueError(f"Scoring function must be one of {', '.join(SCORING_FUNCTIONS)}.")
        if self.seed == 0:
            raise ValueError("Seed 0 tells Vina to pick a random seed; give a nonzero seed so the run can be replayed.")
        for name in ("exhaustiveness", "n_poses"):
            if getattr(self, name) < 1:
                raise ValueError(f"{name} must be at least 1.")
        if self.energy_range <= 0 or self.min_rmsd < 0 or self.spacing <= 0 or self.max_evals < 0 or self.cpu < 0:
            raise ValueError("energy_range and spacing must be positive; min_rmsd, max_evals and cpu non-negative.")

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class DockingRun:
    poses_pdbqt: str
    energies: List[Dict[str, float]]           # one dict per pose, columns named for the scoring function
    settings: DockingSettings
    box: SearchBox
    seconds: float
    notes: List[str] = field(default_factory=list)


def _flush_c_stdio() -> None:
    """Flush the C library's stdio buffers, where Vina's C++ output waits (fflush(NULL))."""
    import ctypes
    try:
        ctypes.CDLL(None).fflush(None)
    except (OSError, AttributeError):
        pass


@contextmanager
def _capture_output(path: str):
    """Send the process's stdout and stderr file descriptors (Vina's C++ output) to ``path``."""
    sys.stdout.flush()
    sys.stderr.flush()
    saved = (os.dup(1), os.dup(2))
    with open(path, "w") as handle:
        os.dup2(handle.fileno(), 1)
        os.dup2(handle.fileno(), 2)
        try:
            yield
        finally:
            _flush_c_stdio()
            sys.stdout.flush()
            sys.stderr.flush()
            os.dup2(saved[0], 1)
            os.dup2(saved[1], 2)
            os.close(saved[0])
            os.close(saved[1])


def _pdbqt_types(pdbqt: str) -> List[str]:
    from .ligand_prep import pdbqt_atom_type
    types: List[str] = []
    for line in pdbqt.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            atom_type = pdbqt_atom_type(line)
            if atom_type and atom_type not in types:
                types.append(atom_type)
    return types


def _autogrid(workdir: str, receptor_file: str, receptor_pdbqt: str, ligand_types: Sequence[str],
              box: SearchBox, spacing: float) -> Tuple[str, List[str]]:
    """AD4 maps for ``ligand_types``; returns the map prefix and notes (warnings AutoGrid printed)."""
    from meeko.gridbox import get_gpf_string
    from .dependencies import find_executable
    executable = find_executable("autogrid4")
    if executable is None:
        raise RuntimeError("AD4 scoring needs autogrid4 (conda-forge package 'autogrid').")
    receptor_types = _pdbqt_types(receptor_pdbqt)
    # dielectric: any negative value selects AutoGrid's distance-dependent (Mehler-Solmajer)
    # function and its magnitude is ignored (identical maps for -42 and -0.1465)
    gpf, npts = get_gpf_string(list(box.center), list(box.size), os.path.basename(receptor_file),
                               receptor_types, list(ligand_types), map_prefix="receptor", spacing=spacing)
    with open(os.path.join(workdir, "receptor.gpf"), "w") as handle:
        handle.write(gpf)
    run = subprocess.run([executable, "-p", "receptor.gpf", "-l", "receptor.glg"], cwd=workdir,
                         capture_output=True, text=True, timeout=3600)
    log = open(os.path.join(workdir, "receptor.glg")).read() if os.path.exists(os.path.join(workdir, "receptor.glg")) else ""
    if run.returncode != 0 or "Successful Completion" not in log:
        raise RuntimeError(f"autogrid4 failed (exit {run.returncode}): {(run.stderr or log)[-500:]}")
    warnings = [line.strip() for line in (log + run.stdout + run.stderr).splitlines() if "WARNING" in line.upper()]
    notes = [f"AutoGrid maps: {npts[0]}x{npts[1]}x{npts[2]} points at {spacing} A, "
             f"types {' '.join(ligand_types)}; distance-dependent dielectric."]
    if warnings:
        notes.append(f"autogrid4 printed {len(warnings)} warnings, e.g. {warnings[0]}")
    return os.path.join(workdir, "receptor"), notes


def build_ad4_maps(workdir: str, receptor_rigid_pdbqt: str, ligand_types: Sequence[str], box: SearchBox,
                   spacing: float) -> Tuple[str, List[str]]:
    """AutoGrid maps for every type in ``ligand_types``, once, for reuse by many ligands.

    Returns the map prefix and notes. A campaign builds them for the union of
    its ligands' types (and any flexible side chains') and passes the prefix
    to ``run_docking(ad4_maps=...)``.
    """
    os.makedirs(workdir, exist_ok=True)
    rigid_file = os.path.join(workdir, "receptor_rigid.pdbqt")
    with open(rigid_file, "w") as handle:
        handle.write(receptor_rigid_pdbqt)
    return _autogrid(workdir, rigid_file, receptor_rigid_pdbqt, list(ligand_types), box, spacing)


def run_docking(receptor_rigid_pdbqt: str, receptor_flex_pdbqt: str, ligand_pdbqt: str,
                box: SearchBox, settings: DockingSettings, workdir: Optional[str] = None,
                ad4_maps: Optional[str] = None) -> DockingRun:
    """Dock one prepared ligand into one prepared receptor.

    ``ad4_maps``: prefix of AutoGrid maps built beforehand (``build_ad4_maps``)
    to reuse instead of running autogrid4 for this ligand; AD4 scoring only.
    """
    from vina import Vina
    settings.validate()
    workdir = workdir or tempfile.mkdtemp(prefix="proprep_docking_")
    os.makedirs(workdir, exist_ok=True)
    rigid_file = os.path.join(workdir, "receptor_rigid.pdbqt")
    with open(rigid_file, "w") as handle:
        handle.write(receptor_rigid_pdbqt)
    flex_file = None
    if receptor_flex_pdbqt:
        flex_file = os.path.join(workdir, "receptor_flex.pdbqt")
        with open(flex_file, "w") as handle:
            handle.write(receptor_flex_pdbqt)
    with open(os.path.join(workdir, "ligand.pdbqt"), "w") as handle:
        handle.write(ligand_pdbqt)

    notes: List[str] = []
    log_path = os.path.join(workdir, "vina_output.txt")
    started = time.time()
    with _capture_output(log_path):
        v = Vina(sf_name=settings.scoring, cpu=settings.cpu, seed=settings.seed, verbosity=1)
        if settings.scoring == "ad4":
            ligand_types = _pdbqt_types(ligand_pdbqt)
            for extra in _pdbqt_types(receptor_flex_pdbqt or ""):
                if extra not in ligand_types:
                    ligand_types.append(extra)       # flexible side chains are scored on the maps too
            if ad4_maps:
                absent = [t for t in ligand_types if not os.path.exists(f"{ad4_maps}.{t}.map")]
                if absent:
                    raise ValueError(f"The AutoGrid maps at {ad4_maps} have no map for type(s) {', '.join(absent)}.")
                prefix = ad4_maps
            else:
                prefix, grid_notes = _autogrid(workdir, rigid_file, receptor_rigid_pdbqt, ligand_types, box,
                                               settings.spacing)
                notes.extend(grid_notes)
            if flex_file:
                v.set_receptor(flex_pdbqt_filename=flex_file)
            v.set_ligand_from_string(ligand_pdbqt)
            v.load_maps(prefix)
        else:
            v.set_receptor(rigid_pdbqt_filename=rigid_file, flex_pdbqt_filename=flex_file)
            v.set_ligand_from_string(ligand_pdbqt)
            v.compute_vina_maps(center=list(box.center), box_size=list(box.size), spacing=settings.spacing)
        v.dock(exhaustiveness=settings.exhaustiveness, n_poses=settings.n_poses,
               min_rmsd=settings.min_rmsd, max_evals=settings.max_evals)
        poses = v.poses(n_poses=settings.n_poses, energy_range=settings.energy_range)
        energies = v.energies(n_poses=settings.n_poses, energy_range=settings.energy_range)
    seconds = time.time() - started
    with open(log_path) as handle:
        vina_output = handle.read()
    for line in vina_output.splitlines():
        if "WARNING" in line.upper():
            notes.append(f"Vina: {line.strip()}")
    columns = _ENERGY_COLUMNS[settings.scoring]
    table = [{name: round(float(value), 3) for name, value in zip(columns, row)} for row in np.atleast_2d(energies)]
    notes.append(f"{len(table)} poses within {settings.energy_range} kcal/mol of the best, "
                 f"{settings.scoring} scoring, seed {settings.seed}, exhaustiveness {settings.exhaustiveness}.")
    return DockingRun(poses_pdbqt=poses, energies=table, settings=settings, box=box, seconds=round(seconds, 1),
                      notes=notes)
