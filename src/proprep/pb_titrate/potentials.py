"""Read pbsa potential map output (`phiform=1`, ASCII format).

The pbsa ASCII potential file (typically `pbsa.phi`) is a flat dump
of the linearized PB potential on the focused fine grid:

    # comment block describing the format
    h gox goy goz                    # grid spacing + origin
    xm ym zm                         # grid dimensions
    phi[1] phi[2] ... phi[xm·ym·zm]  # values, mapped via i + xm*(j-1 + ym*(k-1))

Coordinates of grid point (i, j, k), with i,j,k in 1..{xm,ym,zm}:
    x = gox + h*i
    y = goy + h*j
    z = goz + h*k

Units: kcal/(mol·e) — i.e., potential per unit charge, so that
    Σ q_atom · φ(r_atom)
gives energy in kcal/mol when q is in elementary-charge units.

Used by the O(N) coupling path in `coupling.py` to evaluate
W_ij = Σ over j's atoms of Δq_j · Δφ_i(r_atom).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Tuple

import numpy as np


@dataclass
class PhiGrid:
    """Linearized-PB potential map on a uniform Cartesian grid."""
    h:       float                # grid spacing (Å)
    origin:  Tuple[float, float, float]  # (gox, goy, goz)
    dims:    Tuple[int, int, int]        # (xm, ym, zm)
    phi:     np.ndarray            # shape (xm, ym, zm), units kcal/(mol·e)

    @property
    def xmin(self) -> float:
        return self.origin[0] + self.h        # i=1 → gox + h
    @property
    def ymin(self) -> float:
        return self.origin[1] + self.h
    @property
    def zmin(self) -> float:
        return self.origin[2] + self.h
    @property
    def xmax(self) -> float:
        return self.origin[0] + self.h * self.dims[0]
    @property
    def ymax(self) -> float:
        return self.origin[1] + self.h * self.dims[1]
    @property
    def zmax(self) -> float:
        return self.origin[2] + self.h * self.dims[2]


def read_pbsa_phi(path: Path) -> PhiGrid:
    """Parse a pbsa ASCII phi file (phiform=1)."""
    text = Path(path).read_text()
    lines = [L for L in text.splitlines() if L.strip() and not L.lstrip().startswith("#")]
    # First non-comment non-blank line: h, gox, goy, goz
    parts = lines[0].split()
    if len(parts) < 4:
        raise ValueError(f"Expected 4 numbers (h, gox, goy, goz) on first data "
                          f"line of {path}; got: {parts}")
    h = float(parts[0])
    origin = (float(parts[1]), float(parts[2]), float(parts[3]))
    # Second data line: xm, ym, zm
    parts = lines[1].split()
    if len(parts) < 3:
        raise ValueError(f"Expected 3 dims on second data line of {path}; "
                          f"got: {parts}")
    xm, ym, zm = int(parts[0]), int(parts[1]), int(parts[2])
    # Remaining: flat phi values
    flat = []
    for L in lines[2:]:
        flat.extend(float(x) for x in L.split())
    expected = xm * ym * zm
    if len(flat) != expected:
        raise ValueError(f"Expected {expected} phi values in {path}, "
                          f"got {len(flat)}")
    arr = np.array(flat, dtype=np.float64)
    # Map index = i + xm*(j-1 + ym*(k-1)) with i,j,k starting at 1
    # In numpy: arr[idx-1] for 1-based linear index
    # Reshape: phi[i-1, j-1, k-1] = arr[(i-1) + xm*((j-1) + ym*(k-1))]
    # This is Fortran-order (column-major) reshape with shape (xm, ym, zm)
    phi = arr.reshape((zm, ym, xm)).transpose(2, 1, 0)
    return PhiGrid(h=h, origin=origin, dims=(xm, ym, zm), phi=phi)


def interpolate(grid: PhiGrid, xyz: np.ndarray) -> np.ndarray:
    """Trilinear interpolation of phi at (x, y, z) coordinates.

    `xyz` shape (N, 3) returns shape (N,). Points outside the grid are
    clamped to the boundary (PB potential extrapolates to 0 at infinity
    for finite systems; clamping is safer than extrapolating far outside).
    """
    h = grid.h
    gox, goy, goz = grid.origin
    xm, ym, zm = grid.dims

    # Continuous indices: i_c = (x - gox) / h - 1  (since x = gox + h*i)
    # Wait: with x = gox + h*i and i ∈ [1, xm], x ∈ [gox+h, gox+h*xm].
    # i = (x - gox)/h, so 0-based index = i - 1 = (x - gox)/h - 1
    # Cleaner: use 0-based index directly.
    pts = np.atleast_2d(xyz).astype(np.float64)
    ic = (pts[:, 0] - gox) / h - 1.0
    jc = (pts[:, 1] - goy) / h - 1.0
    kc = (pts[:, 2] - goz) / h - 1.0

    # Clamp to valid range [0, dim-1]
    ic = np.clip(ic, 0.0, xm - 1)
    jc = np.clip(jc, 0.0, ym - 1)
    kc = np.clip(kc, 0.0, zm - 1)

    i0 = np.floor(ic).astype(np.int64); i1 = np.minimum(i0 + 1, xm - 1)
    j0 = np.floor(jc).astype(np.int64); j1 = np.minimum(j0 + 1, ym - 1)
    k0 = np.floor(kc).astype(np.int64); k1 = np.minimum(k0 + 1, zm - 1)
    fi = ic - i0; fj = jc - j0; fk = kc - k0

    p = grid.phi
    c000 = p[i0, j0, k0]; c001 = p[i0, j0, k1]
    c010 = p[i0, j1, k0]; c011 = p[i0, j1, k1]
    c100 = p[i1, j0, k0]; c101 = p[i1, j0, k1]
    c110 = p[i1, j1, k0]; c111 = p[i1, j1, k1]

    c00 = c000 * (1 - fi) + c100 * fi
    c01 = c001 * (1 - fi) + c101 * fi
    c10 = c010 * (1 - fi) + c110 * fi
    c11 = c011 * (1 - fi) + c111 * fi
    c0  = c00 * (1 - fj) + c10 * fj
    c1  = c01 * (1 - fj) + c11 * fj
    out = c0 * (1 - fk) + c1 * fk
    return out


def coupling_from_phi(phi_grid: PhiGrid,
                       atom_xyz: np.ndarray,
                       delta_charges: np.ndarray) -> float:
    """W_ij from a precomputed Δφ_i grid evaluated at site j's atoms.

        W_ij = Σ over atoms a of j: Δq_j(a) · Δφ_i(r_a)

    Returns scalar in kcal/mol (matches φ units of kcal/(mol·e) × q in e).
    """
    phi_at_atoms = interpolate(phi_grid, atom_xyz)
    return float(np.dot(delta_charges, phi_at_atoms))
