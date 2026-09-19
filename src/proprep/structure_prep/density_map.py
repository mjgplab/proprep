"""
Electron Density Maps - a structure against the density it was built into

For an X-ray entry whose structure factors were deposited, PDBe serves two maps
in CCP4 format: the 2mFo-DFc map and the mFo-DFc difference map. They are
CALCULATED by PDBe from the deposited structure factors and the deposited model,
so the 2mFo-DFc map is biased towards that model; the difference map is where
disagreement between model and data shows.

Two things stand between those files and a correct picture, and both are dealt
with here rather than in the browser, where they cannot be tested:

1. A PDBe map covers exactly one unit cell starting at the origin, and a
   deposited model usually lies partly or wholly in neighbouring cells (none of
   the cofactor atoms of 1LTZ or 1J8U are inside the map's box). Density is
   periodic, so the map is re-cut around the model by taking grid indices
   modulo the cell sampling: an exact copy of grid values, no interpolation.

2. A map belongs to the coordinate frame of the deposited entry. A structure
   that ProPrep superposed onto another, or a symmetry copy in a biological
   assembly, no longer lies in it. `frame_check` compares positions, not atom
   names, because ProPrep renames residues and atoms along the way.

Contour levels in sigma refer to the mean and standard deviation of the WHOLE
unit cell, the crystallographic convention. A region cut around a protein holds
more density than the cell on average, so its own statistics would put "1 sigma"
somewhere else; the whole-cell values travel with the re-cut map.
"""

import logging
import os
import struct
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import requests

logger = logging.getLogger(__name__)

# kind -> (what it is, PDBe URL)
MAP_KINDS = {
    "2fofc": ("2mFo-DFc", "https://www.ebi.ac.uk/pdbe/entry-files/{pdb_id}.ccp4"),
    "fofc": ("mFo-DFc difference", "https://www.ebi.ac.uk/pdbe/entry-files/{pdb_id}_diff.ccp4"),
}

# A structure is shown against a map only if its heavy atoms sit, on average, above this
# level of the 2mFo-DFc map. 1 sigma is the level at which that map is conventionally
# contoured to judge whether a model is in density. Measured: models in their own map
# average +2.7 to +3.7 sigma, the same atoms moved elsewhere 0.0 (the map's mean). This
# is the test that does not depend on file headers: a superposed structure that kept its
# CRYST1 and HEADER records agrees with itself perfectly and still fails it.
IN_DENSITY_SIGMA = 1.0

# CCP4 data modes ProPrep reads
_MODE_DTYPES = {0: "i1", 1: "i2", 2: "f4", 6: "u2"}


class DensityMapError(Exception):
    """A map could not be fetched, read or used; the message says why."""


@dataclass
class DensityMap:
    """A CCP4 map, indexed grid[x, y, z] along the cell axes."""
    grid: np.ndarray
    cell: Tuple[float, float, float, float, float, float]
    sampling: Tuple[int, int, int]      # grid intervals along a, b, c of the whole cell
    start: Tuple[int, int, int]         # grid index of grid[0, 0, 0]
    mean: float                         # of the whole-cell map this came from
    sigma: float                        # likewise

    @property
    def covers_whole_cell(self) -> bool:
        return tuple(self.grid.shape) == tuple(self.sampling)

    @property
    def spacing(self) -> float:
        """Largest grid spacing along the cell axes, in Angstroms."""
        return max(self.cell[i] / self.sampling[i] for i in range(3))


# ---------------------------------------------------------------------------
# Fetching
# ---------------------------------------------------------------------------

def fetch_density_maps(pdb_id: str, out_dir: str, timeout: int = 120) -> Optional[Dict[str, str]]:
    """
    Download an entry's two maps from PDBe, next to the structure. Files already
    there are reused.

    Returns:
        {"2fofc": path, "fofc": path}, or None if PDBe has no maps for the entry
        (NMR and cryo-EM entries, and X-ray entries deposited without structure factors)

    Raises:
        DensityMapError: PDBe could not be reached
    """
    paths = {}
    for kind, (_, url) in MAP_KINDS.items():
        path = os.path.join(out_dir, f"{pdb_id.lower()}_{kind}.ccp4")
        if not (os.path.exists(path) and os.path.getsize(path) > 1024):
            try:
                response = requests.get(url.format(pdb_id=pdb_id.lower()), timeout=timeout)
            except requests.exceptions.RequestException as e:
                raise DensityMapError(f"PDBe could not be reached: {e}") from e
            if response.status_code == 404:
                return None
            if response.status_code != 200:
                raise DensityMapError(f"PDBe answered HTTP {response.status_code} for the {kind} map of {pdb_id}")
            with open(path, "wb") as f:
                f.write(response.content)
        paths[kind] = path
    return paths


# ---------------------------------------------------------------------------
# CCP4 format
# ---------------------------------------------------------------------------

def read_ccp4(path: str) -> DensityMap:
    """Read a CCP4/MRC map into grid[x, y, z]."""
    with open(path, "rb") as f:
        raw = f.read()
    if len(raw) < 1024 or raw[208:212] != b"MAP ":
        raise DensityMapError(f"{os.path.basename(path)} is not a CCP4 map")

    # Byte order from the machine stamp; 0x44 is little-endian, 0x11 big-endian
    order = ">" if raw[212] == 0x11 else "<"
    nc, nr, ns, mode, c0, r0, s0, nx, ny, nz = struct.unpack(order + "10i", raw[0:40])
    cell = struct.unpack(order + "6f", raw[40:64])
    mapc, mapr, maps = struct.unpack(order + "3i", raw[64:76])
    nsymbt = struct.unpack(order + "i", raw[92:96])[0]

    if mode not in _MODE_DTYPES:
        raise DensityMapError(f"CCP4 data mode {mode} is not supported")
    if sorted((mapc, mapr, maps)) != [1, 2, 3]:
        raise DensityMapError(f"CCP4 axis order {mapc},{mapr},{maps} is not a permutation of 1,2,3")

    data = np.frombuffer(raw, dtype=order + _MODE_DTYPES[mode], offset=1024 + nsymbt, count=nc * nr * ns)
    # In the file the column index runs fastest, then row, then section; mapc/mapr/maps
    # say which cell axis (1=a, 2=b, 3=c) each of them follows
    by_file_axis = data.reshape(ns, nr, nc).transpose(2, 1, 0)      # [column, row, section]
    file_axis_of = {mapc: 0, mapr: 1, maps: 2}
    grid = by_file_axis.transpose(file_axis_of[1], file_axis_of[2], file_axis_of[3]).astype(np.float32)
    file_start = (c0, r0, s0)
    start = tuple(file_start[file_axis_of[axis]] for axis in (1, 2, 3))

    return DensityMap(grid=grid, cell=tuple(cell), sampling=(nx, ny, nz), start=start,
                      mean=float(grid.mean()), sigma=float(grid.std()))


def write_ccp4(density: DensityMap, path: str) -> None:
    """Write grid[x, y, z] as a CCP4 map (mode 2, axis order a, b, c, little-endian)."""
    nc, nr, ns = density.grid.shape
    header = bytearray(1024)
    struct.pack_into("<10i", header, 0, nc, nr, ns, 2, *density.start, *density.sampling)
    struct.pack_into("<6f", header, 40, *density.cell)
    struct.pack_into("<3i", header, 64, 1, 2, 3)
    struct.pack_into("<3f", header, 76, float(density.grid.min()), float(density.grid.max()), float(density.grid.mean()))
    struct.pack_into("<2i", header, 88, 1, 0)       # P1: a box cut from the cell has no symmetry of its own
    header[208:212] = b"MAP "
    header[212:216] = bytes([0x44, 0x41, 0x00, 0x00])
    struct.pack_into("<f", header, 216, float(density.grid.std()))
    label = b"ProPrep: PDBe map re-cut around the model"
    struct.pack_into("<i", header, 220, 1)
    header[224:224 + len(label)] = label
    with open(path, "wb") as f:
        f.write(bytes(header))
        f.write(np.ascontiguousarray(density.grid.transpose(2, 1, 0), dtype="<f4").tobytes())


# ---------------------------------------------------------------------------
# Coordinates
# ---------------------------------------------------------------------------

def fractionalization_matrix(cell) -> np.ndarray:
    """The matrix taking orthogonal Angstrom coordinates to fractional ones (PDB convention)."""
    a, b, c = cell[:3]
    alpha, beta, gamma = np.radians(cell[3:6])
    cos_a, cos_b, cos_g, sin_g = np.cos(alpha), np.cos(beta), np.cos(gamma), np.sin(gamma)
    volume = np.sqrt(1 - cos_a ** 2 - cos_b ** 2 - cos_g ** 2 + 2 * cos_a * cos_b * cos_g)
    orthogonalization = np.array([
        [a, b * cos_g, c * cos_b],
        [0.0, b * sin_g, c * (cos_a - cos_b * cos_g) / sin_g],
        [0.0, 0.0, c * volume / sin_g],
    ])
    return np.linalg.inv(orthogonalization)


def read_pdb_frame(pdb_file: str) -> Tuple[Optional[tuple], Optional[np.ndarray], Optional[np.ndarray], Optional[str]]:
    """
    Read what ties a PDB file to its crystal.

    Returns:
        (cell from CRYST1, SCALE matrix, SCALE translation, idcode from HEADER); each None when absent
    """
    cell, idcode, rows = None, None, {}
    with open(pdb_file) as f:
        for line in f:
            record = line[:6]
            if record == "HEADER" and len(line) >= 66 and line[62:66].strip():
                idcode = line[62:66].strip().upper()
            elif record == "CRYST1":
                try:
                    cell = tuple(float(line[a:b]) for a, b in ((6, 15), (15, 24), (24, 33), (33, 40), (40, 47), (47, 54)))
                except ValueError:
                    cell = None
            elif record in ("SCALE1", "SCALE2", "SCALE3"):
                try:
                    rows[int(record[5])] = [float(line[10:20]), float(line[20:30]), float(line[30:40]), float(line[45:55])]
                except ValueError:
                    pass
            elif record in ("ATOM  ", "HETATM"):
                break
    if len(rows) == 3:
        scale = np.array([rows[i][:3] for i in (1, 2, 3)])
        shift = np.array([rows[i][3] for i in (1, 2, 3)])
        # A file written without crystal information carries the identity here
        if not np.allclose(scale, np.eye(3)):
            return cell, scale, shift, idcode
    return cell, None, None, idcode


def heavy_atom_coordinates(pdb_file: str) -> np.ndarray:
    """Coordinates of the non-hydrogen atoms of the first model."""
    coordinates = []
    with open(pdb_file) as f:
        for line in f:
            record = line[:6]
            if record == "ENDMDL":
                break
            if record not in ("ATOM  ", "HETATM"):
                continue
            element = line[76:78].strip().upper() if len(line) >= 78 else ""
            if not element:
                element = line[12:14].strip().lstrip("0123456789")[:1].upper()
            if element in ("H", "D"):
                continue
            try:
                coordinates.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                continue
    return np.array(coordinates, dtype=np.float64).reshape(-1, 3)


def nearest_distances(points: np.ndarray, reference: np.ndarray, chunk: int = 512) -> np.ndarray:
    """For each point, the distance to the nearest reference point."""
    nearest = np.empty(len(points))
    for begin in range(0, len(points), chunk):
        block = points[begin:begin + chunk]
        squared = ((block[:, None, :] - reference[None, :, :]) ** 2).sum(axis=2)
        nearest[begin:begin + chunk] = np.sqrt(squared.min(axis=1))
    return nearest


@dataclass
class FrameCheck:
    """Whether a structure still lies where its deposited entry does."""
    atoms: int                  # heavy atoms of the structure
    coincident_fraction: float  # of them, the fraction within `tolerance` of a deposited atom
    median_distance: float      # median distance to the nearest deposited atom, Angstroms
    tolerance: float            # half the map's grid spacing

    @property
    def same_frame(self) -> bool:
        return self.median_distance <= self.tolerance


def frame_check(structure_file: str, deposited_file: str, grid_spacing: float) -> FrameCheck:
    """
    Compare a structure with the deposited entry by position.

    A structure that was only filtered, protonated or renamed keeps its heavy atoms
    exactly where they were deposited. One that was superposed onto another
    structure, or minimized far from the crystal, does not. The test is the median
    distance from each heavy atom to the nearest deposited heavy atom, against half
    the grid spacing of the map: a displacement below that cannot be seen in the map.
    """
    structure = heavy_atom_coordinates(structure_file)
    deposited = heavy_atom_coordinates(deposited_file)
    if not len(structure) or not len(deposited):
        raise DensityMapError("no atoms to compare")
    distances = nearest_distances(structure, deposited)
    tolerance = grid_spacing / 2.0
    return FrameCheck(atoms=len(structure), coincident_fraction=float((distances <= tolerance).mean()),
                      median_distance=float(np.median(distances)), tolerance=tolerance)


# ---------------------------------------------------------------------------
# Re-cutting
# ---------------------------------------------------------------------------

def recut_around(density: DensityMap, coordinates: np.ndarray, margin: float,
                 scale: Optional[np.ndarray] = None, shift: Optional[np.ndarray] = None) -> DensityMap:
    """
    Cut the periodic density around a set of atoms.

    Args:
        density: a map covering one whole unit cell
        coordinates: atom positions in the deposited frame, Angstroms
        margin: how far beyond the outermost atoms the box extends, Angstroms
        scale, shift: the entry's SCALE records; the cell's own matrix when absent

    Returns:
        A map over the box, its grid values copied from the periodic cell
    """
    if not density.covers_whole_cell:
        raise DensityMapError("the map does not cover a whole unit cell, so it cannot be extended periodically")
    if scale is None:
        scale, shift = fractionalization_matrix(density.cell), np.zeros(3)

    fractional = coordinates @ scale.T + shift
    sampling = np.array(density.sampling)
    # The margin in grid units along each axis, from the spacing between grid planes
    plane_spacing = 1.0 / (np.linalg.norm(scale, axis=1) * sampling)
    pad = np.ceil(margin / plane_spacing).astype(int)
    low = np.floor(fractional.min(axis=0) * sampling).astype(int) - pad
    high = np.ceil(fractional.max(axis=0) * sampling).astype(int) + pad

    index = [(np.arange(low[i], high[i] + 1) - density.start[i]) % sampling[i] for i in range(3)]
    region = density.grid[np.ix_(index[0], index[1], index[2])]
    return DensityMap(grid=np.ascontiguousarray(region), cell=density.cell, sampling=density.sampling,
                      start=(int(low[0]), int(low[1]), int(low[2])), mean=density.mean, sigma=density.sigma)


def density_at(density: DensityMap, coordinates: np.ndarray,
               scale: Optional[np.ndarray] = None, shift: Optional[np.ndarray] = None) -> np.ndarray:
    """Map value at the grid point nearest each atom, in sigma above the whole-cell mean."""
    if scale is None:
        scale, shift = fractionalization_matrix(density.cell), np.zeros(3)
    index = np.rint((coordinates @ scale.T + shift) * np.array(density.sampling)).astype(int) - np.array(density.start)
    if density.covers_whole_cell:
        index %= np.array(density.sampling)
    elif (index < 0).any() or (index >= np.array(density.grid.shape)).any():
        raise DensityMapError("an atom lies outside the map")
    values = density.grid[index[:, 0], index[:, 1], index[:, 2]]
    return (values - density.mean) / density.sigma


def mean_density_at_atoms(map_path: str, structure_file: str, frame_file: str) -> float:
    """
    Mean 2mFo-DFc density at a structure's heavy atoms, in whole-cell sigma.

    Args:
        map_path: the whole-cell map as fetched
        structure_file: the structure whose atoms are sampled
        frame_file: the PDB file whose SCALE records define the crystal frame
    """
    _, scale, shift, _ = read_pdb_frame(frame_file)
    coordinates = heavy_atom_coordinates(structure_file)
    if not len(coordinates):
        raise DensityMapError("no atoms to sample the map at")
    return float(density_at(read_ccp4(map_path), coordinates, scale, shift).mean())


def check_standard_frame(density: DensityMap, coordinates: np.ndarray,
                         scale: Optional[np.ndarray], shift: Optional[np.ndarray]) -> None:
    """
    Refuse an entry whose coordinates are not in the standard crystal orientation.

    ProPrep places the re-cut box with the entry's SCALE records. NGL places the map
    from the unit cell alone (a along x, b in the xy plane), which is the orientation
    nearly every entry uses. Where an entry's SCALE says otherwise, the two disagree
    and the density would be drawn off the atoms; half the grid spacing is again the
    displacement below which the map cannot tell.

    Raises:
        DensityMapError: the two frames differ by more than half the grid spacing
    """
    if scale is None:
        return
    standard = fractionalization_matrix(density.cell)
    difference = (coordinates @ scale.T + shift) - (coordinates @ standard.T)
    displacement = float(np.linalg.norm(difference @ np.linalg.inv(standard).T, axis=1).max())
    if displacement > density.spacing / 2.0:
        raise DensityMapError(
            f"this entry's coordinates are not in the standard crystal orientation: its SCALE records place atoms "
            f"up to {displacement:.2f} A from where the unit cell alone does, and the viewer places maps by the unit "
            f"cell. The density would be drawn off the atoms, so it is not shown")


def prepare_maps_for_structure(deposited_file: str, map_paths: Dict[str, str], out_dir: str,
                               margin: float) -> List[Dict]:
    """
    Re-cut an entry's maps around its deposited model.

    Returns:
        One dict per map: kind, label, path of the re-cut file, whole-cell mean and sigma
    """
    _, scale, shift, _ = read_pdb_frame(deposited_file)
    coordinates = heavy_atom_coordinates(deposited_file)
    prepared = []
    for kind, path in map_paths.items():
        density = read_ccp4(path)
        check_standard_frame(density, coordinates, scale, shift)
        region = recut_around(density, coordinates, margin, scale, shift)
        recut_path = os.path.join(out_dir, os.path.basename(path).replace(".ccp4", "_around_model.ccp4"))
        write_ccp4(region, recut_path)
        prepared.append({"kind": kind, "label": MAP_KINDS[kind][0], "path": recut_path,
                         "mean": density.mean, "sigma": density.sigma, "spacing": density.spacing})
    return prepared
