"""
Hessian Matrix Parser
Extracts Hessian matrix and optimized coordinates from QM output files.

© 2024 ProPrep Developer. All rights reserved.

This module parses quantum mechanical output files to extract:
- Cartesian Hessian matrix (second derivatives of energy)
- Optimized atomic coordinates
- Force constants in internal coordinates

Supports:
- Gaussian formatted checkpoint files (.fchk)
- GAMESS-US output files (.log)
- ORCA Hessian files (.hess)

Based on MCPB's get_matrix_from_fchk and get_crds_from_fchk functions.

PROPRIETARY SOFTWARE: This file contains proprietary code belonging to the ProPrep project.
Unauthorized copying, distribution, or modification is strictly prohibited.
"""

import numpy as np
from pathlib import Path
from typing import Tuple, Optional, List
import linecache
import re


class HessianParseError(Exception):
    """Exception raised when Hessian parsing fails."""
    pass


class HessianParser:
    """
    Parse Hessian matrix from QM output files.

    The Hessian matrix contains second derivatives of energy with respect
    to Cartesian coordinates and is essential for computing force constants
    via the Seminario method.
    """

    @staticmethod
    def parse_gaussian_fchk(fchk_file: Path) -> Tuple[np.ndarray, np.ndarray]:
        """
        Extract Hessian matrix and optimized coordinates from Gaussian fchk file.

        Gaussian stores the lower triangle of the symmetric Hessian matrix in
        the formatted checkpoint file. This function reconstructs the full matrix.

        Args:
            fchk_file: Path to Gaussian formatted checkpoint file

        Returns:
            Tuple of (hessian_matrix, coordinates):
                - hessian_matrix: (3N, 3N) symmetric matrix in Hartree/Bohr²
                - coordinates: (3N,) flat array in Bohr

        Raises:
            HessianParseError: If file cannot be parsed or is invalid
        """
        if not fchk_file.exists():
            raise HessianParseError(f"fchk file not found: {fchk_file}")

        # Get number of atoms
        n_atoms = HessianParser._get_n_atoms_from_fchk(fchk_file)

        # Extract optimized coordinates
        coords = HessianParser._extract_coordinates_from_fchk(fchk_file, n_atoms)

        # Extract Hessian matrix
        hessian = HessianParser._extract_hessian_from_fchk(fchk_file, n_atoms)

        return hessian, coords

    @staticmethod
    def parse_gamess_log(log_file: Path) -> Tuple[np.ndarray, np.ndarray]:
        """
        Extract Hessian matrix and coordinates from GAMESS-US log file.

        Args:
            log_file: Path to GAMESS output log file

        Returns:
            Tuple of (hessian_matrix, coordinates)

        Raises:
            HessianParseError: If file cannot be parsed
        """
        if not log_file.exists():
            raise HessianParseError(f"GAMESS log file not found: {log_file}")

        # Extract coordinates
        coords = HessianParser._extract_coordinates_from_gamess(log_file)
        n_atoms = len(coords) // 3

        # Extract Hessian
        hessian = HessianParser._extract_hessian_from_gamess(log_file, n_atoms)

        return hessian, np.array(coords)

    @staticmethod
    def _get_n_atoms_from_fchk(fchk_file: Path) -> int:
        """Get number of atoms from fchk file."""
        with open(fchk_file, 'r') as f:
            for line in f:
                if 'Number of atoms' in line:
                    try:
                        return int(line.split()[-1])
                    except (ValueError, IndexError):
                        raise HessianParseError("Could not parse number of atoms")

        raise HessianParseError("'Number of atoms' not found in fchk file")

    @staticmethod
    def _extract_coordinates_from_fchk(fchk_file: Path, n_atoms: int) -> np.ndarray:
        """
        Extract optimized Cartesian coordinates from fchk file.

        Coordinates are stored in Bohr (atomic units).

        Args:
            fchk_file: Path to fchk file
            n_atoms: Number of atoms

        Returns:
            (3N,) array of coordinates in Bohr
        """
        coords = []
        in_coords_section = False

        with open(fchk_file, 'r') as f:
            for line in f:
                if 'Current cartesian coordinates' in line:
                    in_coords_section = True
                    continue

                if in_coords_section:
                    # Check if line contains numeric data (may start with spaces)
                    stripped = line.strip()
                    if stripped and (stripped[0].isdigit() or stripped[0] == '-'):
                        try:
                            values = line.split()
                            coords.extend([float(v) for v in values])
                        except ValueError:
                            pass  # Skip lines that can't be parsed as floats

                    # Stop when we have all coordinates
                    if len(coords) >= 3 * n_atoms:
                        break

                    # Stop at next section header (lines that don't start with space/digit)
                    if stripped and not (stripped[0].isdigit() or stripped[0] == '-' or line[0] == ' '):
                        if len(coords) > 0:  # Only break if we've found some coords
                            break

        if len(coords) < 3 * n_atoms:
            raise HessianParseError(
                f"Expected {3*n_atoms} coordinates, found {len(coords)}"
            )

        return np.array(coords[:3*n_atoms])

    @staticmethod
    def _extract_hessian_from_fchk(fchk_file: Path, n_atoms: int) -> np.ndarray:
        """
        Extract Hessian matrix from fchk file.

        Gaussian stores the lower triangle of the symmetric Hessian matrix.
        This function reconstructs the full symmetric matrix.

        The Hessian is stored in the order:
        H[0,0], H[1,0], H[1,1], H[2,0], H[2,1], H[2,2], ...

        Args:
            fchk_file: Path to fchk file
            n_atoms: Number of atoms

        Returns:
            (3N, 3N) symmetric Hessian matrix in Hartree/Bohr²
        """
        msize = 3 * n_atoms
        n_elements = msize * (msize + 1) // 2

        # Find Cartesian Force Constants section
        begin_line = None
        n_elements_in_file = None

        with open(fchk_file, 'r') as f:
            for i, line in enumerate(f, 1):
                if 'Cartesian Force Constants' in line:
                    begin_line = i + 1
                    line_parts = line.split()
                    try:
                        n_elements_in_file = int(line_parts[-1])
                    except (ValueError, IndexError):
                        raise HessianParseError("Could not parse number of force constants")
                    break

        if begin_line is None:
            raise HessianParseError(
                "No 'Cartesian Force Constants' section found in fchk file. "
                "Ensure the QM calculation included frequency analysis."
            )

        if n_elements != n_elements_in_file:
            raise HessianParseError(
                f"Matrix size mismatch: expected {n_elements} elements "
                f"for {n_atoms} atoms, found {n_elements_in_file}"
            )

        # Read force constants (5 values per line)
        fc_values = []
        n_lines = (n_elements + 4) // 5  # Ceiling division

        for line_num in range(begin_line, begin_line + n_lines):
            line = linecache.getline(str(fchk_file), line_num)
            if not line:
                raise HessianParseError(f"Unexpected end of file at line {line_num}")

            values = line.split()
            fc_values.extend([float(v) for v in values])

        linecache.clearcache()

        if len(fc_values) < n_elements:
            raise HessianParseError(
                f"Not enough force constants: expected {n_elements}, found {len(fc_values)}"
            )

        # Build symmetric matrix from lower triangle
        hessian = np.zeros((msize, msize))
        idx = 0

        for j in range(msize):
            for k in range(j + 1):
                hessian[j, k] = fc_values[idx]
                hessian[k, j] = fc_values[idx]  # Symmetry
                idx += 1

        return hessian

    @staticmethod
    def _extract_coordinates_from_gamess(log_file: Path) -> List[float]:
        """
        Extract optimized coordinates from GAMESS log file.

        Args:
            log_file: Path to GAMESS log file

        Returns:
            Flat list of coordinates in Bohr
        """
        coords = []
        in_coords_section = False

        with open(log_file, 'r') as f:
            for line in f:
                # Look for final coordinates section
                if 'COORDINATES OF ALL ATOMS ARE' in line:
                    in_coords_section = True
                    next(f)  # Skip header line
                    next(f)  # Skip blank line
                    continue

                if in_coords_section:
                    if line.strip() == '':
                        break

                    parts = line.split()
                    if len(parts) >= 5:
                        try:
                            # GAMESS format: ATOM CHARGE X Y Z
                            x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                            coords.extend([x, y, z])
                        except (ValueError, IndexError):
                            continue

        if not coords:
            raise HessianParseError("Could not find coordinates in GAMESS log file")

        return coords

    @staticmethod
    def _extract_hessian_from_gamess(log_file: Path, n_atoms: int) -> np.ndarray:
        """
        Extract Hessian matrix from GAMESS log file.

        Args:
            log_file: Path to GAMESS log file
            n_atoms: Number of atoms

        Returns:
            (3N, 3N) Hessian matrix
        """
        msize = 3 * n_atoms
        hessian = np.zeros((msize, msize))

        in_hessian_section = False
        current_col_block = 0

        with open(log_file, 'r') as f:
            for line in f:
                if 'HESSIAN MATRIX' in line:
                    in_hessian_section = True
                    continue

                if in_hessian_section:
                    # Parse column headers
                    if line.strip().startswith('1') and len(line.split()) >= 5:
                        col_indices = [int(x) - 1 for x in line.split()]
                        current_col_block = col_indices
                        continue

                    # Parse matrix row
                    parts = line.split()
                    if len(parts) > 1 and parts[0].isdigit():
                        row_idx = int(parts[0]) - 1
                        values = [float(x) for x in parts[1:]]

                        for col_offset, value in enumerate(values):
                            if col_offset < len(current_col_block):
                                col_idx = current_col_block[col_offset]
                                hessian[row_idx, col_idx] = value

        if np.all(hessian == 0):
            raise HessianParseError("Could not extract Hessian from GAMESS log file")

        # Make symmetric
        hessian = (hessian + hessian.T) / 2

        return hessian

    @staticmethod
    def validate_hessian(hessian: np.ndarray, coords: np.ndarray) -> Tuple[bool, str]:
        """
        Validate extracted Hessian matrix.

        Checks:
        - Matrix is square
        - Matrix is symmetric
        - Size matches coordinates
        - No NaN or infinite values

        Args:
            hessian: Hessian matrix to validate
            coords: Coordinate array

        Returns:
            (is_valid, message) tuple
        """
        n_coords = len(coords)
        expected_size = n_coords

        # Check square
        if hessian.shape[0] != hessian.shape[1]:
            return False, f"Hessian is not square: {hessian.shape}"

        # Check size matches coordinates
        if hessian.shape[0] != expected_size:
            return False, (
                f"Hessian size {hessian.shape[0]} does not match "
                f"coordinate size {expected_size}"
            )

        # Check symmetry
        if not np.allclose(hessian, hessian.T, rtol=1e-5):
            max_asymmetry = np.max(np.abs(hessian - hessian.T))
            return False, f"Hessian is not symmetric (max diff: {max_asymmetry:.2e})"

        # Check for invalid values
        if np.any(np.isnan(hessian)):
            return False, "Hessian contains NaN values"

        if np.any(np.isinf(hessian)):
            return False, "Hessian contains infinite values"

        # Check for all zeros
        if np.all(hessian == 0):
            return False, "Hessian is all zeros"

        return True, "Hessian is valid"

    @staticmethod
    def get_hessian_statistics(hessian: np.ndarray) -> dict:
        """
        Get statistical information about the Hessian matrix.

        Useful for debugging and quality assessment.

        Args:
            hessian: Hessian matrix

        Returns:
            Dict with statistics
        """
        eigenvalues = np.linalg.eigvalsh(hessian)

        return {
            "shape": hessian.shape,
            "min_value": float(np.min(hessian)),
            "max_value": float(np.max(hessian)),
            "mean_value": float(np.mean(hessian)),
            "min_eigenvalue": float(np.min(eigenvalues)),
            "max_eigenvalue": float(np.max(eigenvalues)),
            "n_negative_eigenvalues": int(np.sum(eigenvalues < -1e-6)),
            "condition_number": float(np.linalg.cond(hessian))
        }
