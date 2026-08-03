"""
Seminario Method Implementation
Computes bond and angle force constants from QM Hessian matrix.

© 2024 ProPrep Developer. All rights reserved.

This module implements the Seminario method for extracting force constants
from quantum mechanical Hessian matrices. The method projects sub-blocks of
the Hessian onto internal coordinate directions (bonds and angles) to obtain
force constants suitable for molecular mechanics force fields.

Reference:
    J. M. Seminario, "Calculation of intramolecular force fields from
    second-derivative tensors", International Journal of Quantum Chemistry,
    1996, 60(7), 1271-1277.

PROPRIETARY SOFTWARE: This file contains proprietary code belonging to the ProPrep project.
Unauthorized copying, distribution, or modification is strictly prohibited.
"""

import numpy as np
from typing import Tuple, List, Dict, Optional
from dataclasses import dataclass
from rich.console import Console

# Unit conversion constants (from atomic units to AMBER units)
B_TO_A = 0.529177249           # Bohr to Angstrom
H_TO_KCAL_MOL = 627.5095       # Hartree to kcal/mol
HB2_TO_KCAL_MOL_A2 = H_TO_KCAL_MOL / (B_TO_A ** 2)  # Hartree/Bohr² to kcal/(mol·Å²) = 2240.87


@dataclass
class BondParameter:
    """
    Bond force constant and equilibrium length.

    Attributes:
        atom1_idx: 0-indexed first atom
        atom2_idx: 0-indexed second atom
        atom1_type: Atom type name (e.g., "M1", "L1")
        atom2_type: Atom type name
        force_constant: Force constant in kcal/(mol·Å²)
        eq_length: Equilibrium bond length in Angstrom
        std_dev: Standard deviation if averaged over matrix orderings
    """
    atom1_idx: int
    atom2_idx: int
    atom1_type: str
    atom2_type: str
    force_constant: float
    eq_length: float
    std_dev: float = 0.0


@dataclass
class AngleParameter:
    """
    Angle force constant and equilibrium angle.

    Attributes:
        atom1_idx: 0-indexed first atom
        atom2_idx: 0-indexed second atom (central)
        atom3_idx: 0-indexed third atom
        atom1_type: Atom type name
        atom2_type: Atom type name (central)
        atom3_type: Atom type name
        force_constant: Force constant in kcal/(mol·rad²)
        eq_angle: Equilibrium angle in degrees
        std_dev: Standard deviation if averaged
    """
    atom1_idx: int
    atom2_idx: int
    atom3_idx: int
    atom1_type: str
    atom2_type: str
    atom3_type: str
    force_constant: float
    eq_angle: float
    std_dev: float = 0.0


class SeminarioMethod:
    """
    Seminario method for extracting force constants from Hessian matrix.

    Projects sub-blocks of the Hessian matrix onto internal coordinate
    directions to obtain force constants. The method is particularly
    well-suited for metal-ligand interactions where quantum mechanical
    calculations provide accurate second derivatives.

    The Hessian matrix contains second derivatives of energy with respect
    to Cartesian coordinates:
        H[i,j] = ∂²E / ∂x_i ∂x_j

    The Seminario method extracts force constants by:
    1. Selecting appropriate sub-blocks of H
    2. Diagonalizing to find eigenvalues (vibrational modes)
    3. Projecting onto the internal coordinate direction
    """

    def __init__(self, hessian: np.ndarray, coords: np.ndarray,
                 scale_factor: float = 1.0, console: Optional[Console] = None):
        """
        Initialize Seminario method.

        Args:
            hessian: (3N, 3N) Hessian matrix in Hartree/Bohr²
            coords: (3N,) coordinates in Bohr
            scale_factor: Frequency scaling factor (default 1.0)
                         Applied as (scale_factor)² to force constants
            console: Rich console for output (optional)
        """
        self.hessian = hessian
        self.coords = coords
        self.scale_factor = scale_factor
        self.n_atoms = len(coords) // 3
        self.console = console or Console()

        # Validate inputs
        if hessian.shape[0] != hessian.shape[1]:
            raise ValueError(f"Hessian must be square, got {hessian.shape}")

        if hessian.shape[0] != 3 * self.n_atoms:
            raise ValueError(
                f"Hessian size {hessian.shape[0]} does not match "
                f"coordinate size {len(coords)}"
            )

    def compute_bond_force_constant(
        self,
        atom1_idx: int,
        atom2_idx: int,
        atom1_type: str,
        atom2_type: str,
        crystal_coords: Optional[Tuple[float, ...]] = None,
        average_matrices: bool = True
    ) -> BondParameter:
        """
        Compute bond force constant using Seminario method.

        The method extracts a 3×3 sub-block of the Hessian between atoms i and j,
        diagonalizes it, and projects the eigenvalues onto the bond direction.

        Args:
            atom1_idx: 0-indexed atom 1
            atom2_idx: 0-indexed atom 2
            atom1_type: Atom type name (e.g., "M1", "L1")
            atom2_type: Atom type name
            crystal_coords: Optional (x1,y1,z1,x2,y2,z2) in Angstrom
                          If provided, used for equilibrium length
            average_matrices: Average over H[i,j] and H[j,i] (recommended)

        Returns:
            BondParameter with force constant and equilibrium length
        """
        # Extract coordinates (in Bohr)
        c1 = self.coords[3*atom1_idx:3*atom1_idx+3]
        c2 = self.coords[3*atom2_idx:3*atom2_idx+3]

        # Bond vector and distance
        vec12 = c2 - c1
        dist_bohr = np.linalg.norm(vec12)

        if dist_bohr < 1e-6:
            raise ValueError(f"Atoms {atom1_idx} and {atom2_idx} have identical coordinates")

        unit_vec = vec12 / dist_bohr

        # Equilibrium distance (use crystal structure if provided)
        if crystal_coords:
            c1_cryst = np.array(crystal_coords[:3])
            c2_cryst = np.array(crystal_coords[3:])
            eq_dist = np.linalg.norm(c2_cryst - c1_cryst)  # Already in Angstrom
        else:
            eq_dist = dist_bohr * B_TO_A

        # Extract 3×3 sub-Hessian: H[3i:3i+3, 3j:3j+3]
        # Negative sign for force constant definition
        sub_hessian = -self.hessian[3*atom1_idx:3*atom1_idx+3,
                                     3*atom2_idx:3*atom2_idx+3]

        # Eigendecomposition
        eigvals, eigvecs = np.linalg.eig(sub_hessian)

        # Project onto bond direction
        # Sum of eigenvalue × |projection of eigenvector onto bond|
        # Note: For asymmetric sub-Hessians, eigenvalues/vectors may be complex.
        # We use full complex values and take .real at the end (matching MCPB).
        fc1 = sum(eigvals[i] * abs(np.dot(eigvecs[:,i], unit_vec))
                  for i in range(3))

        if average_matrices:
            # Average over transposed matrix H[j,i]
            sub_hessian_T = -self.hessian[3*atom2_idx:3*atom2_idx+3,
                                           3*atom1_idx:3*atom1_idx+3]
            eigvals_T, eigvecs_T = np.linalg.eig(sub_hessian_T)

            fc2 = sum(eigvals_T[i] * abs(np.dot(eigvecs_T[:,i], unit_vec))
                     for i in range(3))

            # Take real part after averaging (complex parts should cancel for conjugate pairs)
            fc_mean = np.real((fc1 + fc2) / 2)
            fc_std = np.std([np.real(fc1), np.real(fc2)])
        else:
            fc_mean = np.real(fc1)
            fc_std = 0.0

        # Convert to kcal/(mol·Å²) and apply scaling
        # Factor 0.5 because AMBER uses k(r-r0)² not (1/2)k(r-r0)²
        fc_final = fc_mean * HB2_TO_KCAL_MOL_A2 * 0.5 * self.scale_factor**2
        fc_std_final = fc_std * HB2_TO_KCAL_MOL_A2 * 0.5 * self.scale_factor**2

        # Ensure non-negative force constant
        if fc_final < 0:
            self.console.print(
                f"[yellow]⚠️  Warning: Negative force constant for {atom1_type}-{atom2_type} "
                f"({fc_final:.1f}), setting to 0.0[/yellow]"
            )
            fc_final = 0.0

        return BondParameter(
            atom1_idx=atom1_idx,
            atom2_idx=atom2_idx,
            atom1_type=atom1_type,
            atom2_type=atom2_type,
            force_constant=fc_final,
            eq_length=eq_dist,
            std_dev=fc_std_final
        )

    def compute_angle_force_constant(
        self,
        atom1_idx: int,
        atom2_idx: int,
        atom3_idx: int,
        atom1_type: str,
        atom2_type: str,
        atom3_type: str,
        crystal_coords: Optional[np.ndarray] = None,
        average_matrices: bool = True
    ) -> AngleParameter:
        """
        Compute angle force constant using Seminario method.

        The method uses two 3×3 sub-Hessians (for bonds A-B and C-B) and
        projects them onto directions perpendicular to the bonds within
        the plane of the angle.

        Args:
            atom1_idx: 0-indexed atom A
            atom2_idx: 0-indexed atom B (central)
            atom3_idx: 0-indexed atom C
            atom1_type: Atom type name A
            atom2_type: Atom type name B (central)
            atom3_type: Atom type name C
            crystal_coords: Optional (9,) array [x1,y1,z1,x2,y2,z2,x3,y3,z3]
                          in Angstrom for equilibrium angle
            average_matrices: Average over 4 matrix orderings (recommended)

        Returns:
            AngleParameter with force constant and equilibrium angle
        """
        # Extract coordinates (in Bohr)
        c1 = self.coords[3*atom1_idx:3*atom1_idx+3]
        c2 = self.coords[3*atom2_idx:3*atom2_idx+3]
        c3 = self.coords[3*atom3_idx:3*atom3_idx+3]

        # Vectors and distances
        vec12 = c2 - c1
        vec32 = c2 - c3
        dist12 = np.linalg.norm(vec12)
        dist32 = np.linalg.norm(vec32)

        if dist12 < 1e-6 or dist32 < 1e-6:
            raise ValueError(f"Degenerate angle: atoms {atom1_idx}-{atom2_idx}-{atom3_idx}")

        # Unit vectors
        u12 = vec12 / dist12
        u32 = vec32 / dist32

        # Angle value
        cos_angle = np.dot(u12, u32)
        cos_angle = np.clip(cos_angle, -1.0, 1.0)  # Handle numerical errors
        angle_rad = np.arccos(cos_angle)
        angle_deg = np.degrees(angle_rad)

        # Use crystal structure angle if provided
        if crystal_coords is not None:
            c1_cryst = crystal_coords[:3]
            c2_cryst = crystal_coords[3:6]
            c3_cryst = crystal_coords[6:]
            v12_cryst = c2_cryst - c1_cryst
            v32_cryst = c2_cryst - c3_cryst
            cos_cryst = np.dot(v12_cryst, v32_cryst) / (
                np.linalg.norm(v12_cryst) * np.linalg.norm(v32_cryst)
            )
            cos_cryst = np.clip(cos_cryst, -1.0, 1.0)
            angle_deg = np.degrees(np.arccos(cos_cryst))

        # Normal to plane ABC
        normal = np.cross(u32, u12)
        normal_norm = np.linalg.norm(normal)

        if normal_norm < 1e-6:
            raise ValueError(f"Linear angle detected: {atom1_idx}-{atom2_idx}-{atom3_idx}")

        normal = normal / normal_norm

        # Perpendicular vectors in plane
        u_PA = np.cross(normal, u12)
        u_PC = np.cross(u32, normal)

        # Compute force constant (average over 4 matrix orderings)
        fc_values = []

        # Ordering 1: H[A,B] and H[C,B]
        fc_values.append(self._compute_angle_fc_single(
            atom1_idx, atom2_idx, atom3_idx,
            u_PA, u_PC, dist12, dist32
        ))

        if average_matrices:
            # Ordering 2: H[B,A] and H[C,B]
            fc_values.append(self._compute_angle_fc_single(
                atom1_idx, atom2_idx, atom3_idx,
                u_PA, u_PC, dist12, dist32,
                transpose_12=True
            ))

            # Ordering 3: H[A,B] and H[B,C]
            fc_values.append(self._compute_angle_fc_single(
                atom1_idx, atom2_idx, atom3_idx,
                u_PA, u_PC, dist12, dist32,
                transpose_32=True
            ))

            # Ordering 4: H[B,A] and H[B,C]
            fc_values.append(self._compute_angle_fc_single(
                atom1_idx, atom2_idx, atom3_idx,
                u_PA, u_PC, dist12, dist32,
                transpose_12=True, transpose_32=True
            ))

        fc_mean = np.mean(fc_values)
        fc_std = np.std(fc_values) if len(fc_values) > 1 else 0.0

        # Convert to kcal/(mol·rad²) and scale
        # Factor 0.5 because AMBER uses k(θ-θ0)² not (1/2)k(θ-θ0)²
        fc_final = fc_mean * H_TO_KCAL_MOL * 0.5 * self.scale_factor**2
        fc_std_final = fc_std * H_TO_KCAL_MOL * 0.5 * self.scale_factor**2

        # Ensure non-negative force constant
        if fc_final < 0:
            self.console.print(
                f"[yellow]⚠️  Warning: Negative force constant for "
                f"{atom1_type}-{atom2_type}-{atom3_type} ({fc_final:.2f}), setting to 0.0[/yellow]"
            )
            fc_final = 0.0

        return AngleParameter(
            atom1_idx=atom1_idx,
            atom2_idx=atom2_idx,
            atom3_idx=atom3_idx,
            atom1_type=atom1_type,
            atom2_type=atom2_type,
            atom3_type=atom3_type,
            force_constant=fc_final,
            eq_angle=angle_deg,
            std_dev=fc_std_final
        )

    def _compute_angle_fc_single(
        self,
        atom1_idx: int,
        atom2_idx: int,
        atom3_idx: int,
        u_PA: np.ndarray,
        u_PC: np.ndarray,
        dist12: float,
        dist32: float,
        transpose_12: bool = False,
        transpose_32: bool = False
    ) -> float:
        """
        Compute angle force constant for one matrix ordering.

        Uses the compliance formula:
            k_angle = 1 / (1/k_AB + 1/k_CB)

        where k_AB and k_CB are computed from projections of Hessian eigenvalues.

        Args:
            atom1_idx: First atom index
            atom2_idx: Central atom index
            atom3_idx: Third atom index
            u_PA: Unit vector perpendicular to bond A-B in plane
            u_PC: Unit vector perpendicular to bond C-B in plane
            dist12: Distance A-B in Bohr
            dist32: Distance C-B in Bohr
            transpose_12: Use H[B,A] instead of H[A,B]
            transpose_32: Use H[B,C] instead of H[C,B]

        Returns:
            Force constant in Hartree/radian²
        """
        # Extract sub-Hessians with optional transposition
        if transpose_12:
            H_12 = -self.hessian[3*atom2_idx:3*atom2_idx+3,
                                 3*atom1_idx:3*atom1_idx+3]
        else:
            H_12 = -self.hessian[3*atom1_idx:3*atom1_idx+3,
                                 3*atom2_idx:3*atom2_idx+3]

        if transpose_32:
            H_32 = -self.hessian[3*atom2_idx:3*atom2_idx+3,
                                 3*atom3_idx:3*atom3_idx+3]
        else:
            H_32 = -self.hessian[3*atom3_idx:3*atom3_idx+3,
                                 3*atom2_idx:3*atom2_idx+3]

        # Eigendecompositions
        eigvals_12, eigvecs_12 = np.linalg.eig(H_12)
        eigvals_32, eigvecs_32 = np.linalg.eig(H_32)

        # Project onto perpendicular directions
        # Note: For asymmetric sub-Hessians, eigenvalues/vectors may be complex.
        # We use full complex values and take .real at the end (matching MCPB).
        contrib_12 = sum(eigvals_12[i] * abs(np.dot(u_PA, eigvecs_12[:,i]))
                        for i in range(3))
        contrib_32 = sum(eigvals_32[i] * abs(np.dot(u_PC, eigvecs_32[:,i]))
                        for i in range(3))

        # Take real part before compliance calculation
        contrib_12 = np.real(contrib_12)
        contrib_32 = np.real(contrib_32)

        # Avoid division by zero
        if abs(contrib_12) < 1e-10 or abs(contrib_32) < 1e-10:
            return 0.0

        # Compliance formula
        compliance_12 = 1.0 / (contrib_12 * dist12**2)
        compliance_32 = 1.0 / (contrib_32 * dist32**2)

        fc = 1.0 / (compliance_12 + compliance_32)

        return fc

    def compute_all_parameters(
        self,
        bonds: List[Tuple[int, int, str, str]],
        angles: List[Tuple[int, int, int, str, str, str]],
        crystal_structure: Optional[Dict] = None
    ) -> Tuple[List[BondParameter], List[AngleParameter]]:
        """
        Compute force constants for all bonds and angles.

        Args:
            bonds: List of (atom1_idx, atom2_idx, atom1_type, atom2_type)
            angles: List of (atom1_idx, atom2_idx, atom3_idx,
                            atom1_type, atom2_type, atom3_type)
            crystal_structure: Optional dict with crystal coordinates for
                             equilibrium values

        Returns:
            (bond_parameters, angle_parameters) tuple
        """
        bond_params = []
        angle_params = []

        # Compute bond parameters
        for atom1_idx, atom2_idx, atom1_type, atom2_type in bonds:
            crystal_coords = None
            if crystal_structure and (atom1_idx, atom2_idx) in crystal_structure.get('bonds', {}):
                crystal_coords = crystal_structure['bonds'][(atom1_idx, atom2_idx)]

            param = self.compute_bond_force_constant(
                atom1_idx, atom2_idx, atom1_type, atom2_type,
                crystal_coords=crystal_coords
            )
            bond_params.append(param)

        # Compute angle parameters
        for atom1_idx, atom2_idx, atom3_idx, atom1_type, atom2_type, atom3_type in angles:
            crystal_coords = None
            if crystal_structure and (atom1_idx, atom2_idx, atom3_idx) in crystal_structure.get('angles', {}):
                crystal_coords = crystal_structure['angles'][(atom1_idx, atom2_idx, atom3_idx)]

            param = self.compute_angle_force_constant(
                atom1_idx, atom2_idx, atom3_idx,
                atom1_type, atom2_type, atom3_type,
                crystal_coords=crystal_coords
            )
            angle_params.append(param)

        return bond_params, angle_params
