"""
How close a solute comes to its own periodic images.

Under periodic boundary conditions the simulated system is one cell of an
infinite lattice, so the solute has copies of itself all around it. Where the
solute sits in the box means nothing (the whole system can be slid without
changing the physics). What does mean something is the shortest distance from
any atom of the solute to any atom of any copy: once that falls below the
nonbonded cutoff, the solute interacts directly with itself. It changes along a
run, because the solute tumbles, unfolds or extends, and because the box itself
shrinks during constant-pressure equilibration.

Any Amber box is handled: a rectangular one, a truncated octahedron (which
Amber stores as a triclinic cell with all three angles 109.47 degrees) and a
general triclinic one. The copies are generated from the three lattice vectors,
so nothing here is specific to a shape.

The search uses a k-d tree (scikit-learn, already a ProPrep dependency), so its
cost grows close to linearly with the number of atoms. cpptraj's ``minimage``
gives the same numbers (checked on an orthorhombic and a truncated-octahedral
trajectory, to 0.001 A) but compares every atom with every atom of every
image: 0.9 s a frame for 6,000 atoms, against 0.2 s here, and the gap widens
with the square of the size.
"""

import re
from itertools import product
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import numpy as np

# The angle of Amber's truncated octahedron is acos(-1/3) = 109.4712 degrees;
# trajectories record it to two decimals. Used only to NAME the shape.
TRUNCATED_OCTAHEDRON_ANGLE = float(np.degrees(np.arccos(-1.0 / 3.0)))
ANGLE_TOLERANCE = 0.02        # degrees
LENGTH_TOLERANCE = 1.0e-3     # relative, for "all three edges equal"

# Half of the 26 neighbouring cells: the distance to the copy at +v is the
# distance to the copy at -v, so the other half would repeat the work.
_HALF_SHIFTS = [s for s in product((-1, 0, 1), repeat=3) if s > (0, 0, 0)]


class NotPeriodic(ValueError):
    """The frame has no periodic box, so it has no images."""


def lattice_vectors(box: Sequence[float]) -> np.ndarray:
    """The three cell vectors (rows) from ``a, b, c, alpha, beta, gamma``, Amber's convention:
    a along x, b in the xy plane."""
    if box is None or len(box) < 6:
        raise NotPeriodic("no box information")
    a, b, c, alpha, beta, gamma = (float(v) for v in box[:6])
    if min(a, b, c) <= 0.0:
        raise NotPeriodic("the box has no size")
    alpha, beta, gamma = np.radians([alpha, beta, gamma])
    cx = c * np.cos(beta)
    cy = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    cz = np.sqrt(max(c * c - cx * cx - cy * cy, 0.0))
    return np.array([[a, 0.0, 0.0],
                     [b * np.cos(gamma), b * np.sin(gamma), 0.0],
                     [cx, cy, cz]])


def box_shape(box: Sequence[float]) -> str:
    """"rectangular", "truncated octahedron" or "triclinic", from the box angles and edges."""
    a, b, c, alpha, beta, gamma = (float(v) for v in box[:6])
    angles = np.array([alpha, beta, gamma])
    if np.all(np.abs(angles - 90.0) < ANGLE_TOLERANCE):
        return "rectangular"
    equal_edges = max(a, b, c) - min(a, b, c) < LENGTH_TOLERANCE * max(a, b, c)
    if equal_edges and np.all(np.abs(angles - TRUNCATED_OCTAHEDRON_ANGLE) < ANGLE_TOLERANCE):
        return "truncated octahedron"
    return "triclinic"


def minimum_image_distance(xyz: np.ndarray, box: Sequence[float]) -> Tuple[float, int, int, Tuple[int, int, int]]:
    """
    Shortest distance from any of the atoms to any atom of any periodic copy of them.

    Returns:
        (distance, index of the atom, index of the atom of the copy, the copy's cell
        as multiples of the three lattice vectors)
    """
    from sklearn.neighbors import KDTree

    xyz = np.asarray(xyz, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or len(xyz) == 0:
        raise ValueError("expected an (n_atoms, 3) array with at least one atom")
    vectors = lattice_vectors(box)
    tree = KDTree(xyz)
    best = (np.inf, -1, -1, (0, 0, 0))
    for shift in _HALF_SHIFTS:
        distances, nearest = tree.query(xyz + np.asarray(shift, dtype=float) @ vectors, k=1)
        in_copy = int(distances[:, 0].argmin())
        if distances[in_copy, 0] < best[0]:
            best = (float(distances[in_copy, 0]), int(nearest[in_copy, 0]), in_copy, shift)
    return best


def read_cutoff(sim_dir) -> Optional[Tuple[float, str]]:
    """
    The nonbonded cutoff a simulation used: ``(value in A, file it was read from)``, or None.

    The mdout is preferred, because it echoes what the engine actually ran with
    ("cut     =  10.00000"); the input file is the fallback.
    """
    pattern = re.compile(r"\bcut\s*=\s*([0-9]*\.?[0-9]+)", re.IGNORECASE)
    folder = Path(sim_dir)
    if not folder.is_dir():
        return None
    for glob in ("*.mdout", "*.out", "*.mdin", "*.in"):
        for path in sorted(folder.glob(glob)):
            try:
                with open(path, errors="replace") as handle:
                    for line in handle:
                        found = pattern.search(line.split("!")[0])
                        if found:
                            return float(found.group(1)), path.name
            except OSError:
                continue
    return None


def describe_frames_below(distances: Sequence[float], frames: Sequence[int], cutoff: float) -> List[int]:
    """Frame numbers (as given) whose image distance is below the cutoff."""
    return [int(f) for d, f in zip(distances, frames) if d < cutoff]
