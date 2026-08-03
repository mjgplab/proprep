"""
Structure Orientation Module

Aligns protein structures along Cartesian axes using various methods:
- Principal axes alignment (PCA/inertia tensor)
- Two residues/atoms to axis alignment
- Plane fitting with normal to Z
- Center at origin
- Custom rotation (Euler angles)

Useful for: membrane protein setup, MD with electric fields, steered MD,
visualization consistency, QM/MM calculations.
"""

import logging
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Union
from dataclasses import dataclass
from datetime import datetime

import numpy as np

from proprep.utils.module_registry import ProcessingModule, register_module
from proprep.utils.prompts import (
    prompt_with_context,
    confirm_with_context,
    int_prompt_with_context,
    float_prompt_with_context
)

logger = logging.getLogger(__name__)


@dataclass
class TransformationRecord:
    """Record of applied transformation for reproducibility."""
    method: str
    rotation_matrix: np.ndarray
    translation_vector: np.ndarray
    parameters: Dict[str, Any]
    timestamp: str

    def to_dict(self) -> Dict:
        """Convert to JSON-serializable dictionary."""
        return {
            'method': self.method,
            'rotation_matrix': self.rotation_matrix.tolist(),
            'translation_vector': self.translation_vector.tolist(),
            'parameters': self.parameters,
            'timestamp': self.timestamp
        }

    @classmethod
    def from_dict(cls, data: Dict) -> 'TransformationRecord':
        """Create from dictionary."""
        return cls(
            method=data['method'],
            rotation_matrix=np.array(data['rotation_matrix']),
            translation_vector=np.array(data['translation_vector']),
            parameters=data['parameters'],
            timestamp=data['timestamp']
        )


def calculate_center_of_mass(coords: np.ndarray, masses: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Calculate center of mass (or geometric center if no masses provided).

    Args:
        coords: Nx3 array of coordinates
        masses: Optional N array of atomic masses

    Returns:
        3-element array with center coordinates
    """
    if masses is None:
        return np.mean(coords, axis=0)
    else:
        total_mass = np.sum(masses)
        return np.sum(coords * masses[:, np.newaxis], axis=0) / total_mass


def calculate_inertia_tensor(coords: np.ndarray, masses: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Calculate the inertia tensor for a set of coordinates.

    Args:
        coords: Nx3 array of coordinates (should be centered at origin)
        masses: Optional N array of atomic masses (uses unit masses if None)

    Returns:
        3x3 inertia tensor
    """
    if masses is None:
        masses = np.ones(len(coords))

    I = np.zeros((3, 3))

    for i, (coord, m) in enumerate(zip(coords, masses)):
        x, y, z = coord
        I[0, 0] += m * (y**2 + z**2)
        I[1, 1] += m * (x**2 + z**2)
        I[2, 2] += m * (x**2 + y**2)
        I[0, 1] -= m * x * y
        I[0, 2] -= m * x * z
        I[1, 2] -= m * y * z

    # Symmetric tensor
    I[1, 0] = I[0, 1]
    I[2, 0] = I[0, 2]
    I[2, 1] = I[1, 2]

    return I


def get_principal_axes(coords: np.ndarray, masses: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate principal axes of inertia.

    Args:
        coords: Nx3 array of coordinates
        masses: Optional atomic masses

    Returns:
        Tuple of (eigenvalues, eigenvectors) where eigenvectors are columns
        Sorted by eigenvalue (smallest to largest, so shortest axis first)
    """
    # Center coordinates
    com = calculate_center_of_mass(coords, masses)
    centered = coords - com

    # Calculate inertia tensor
    I = calculate_inertia_tensor(centered, masses)

    # Diagonalize
    eigenvalues, eigenvectors = np.linalg.eigh(I)

    # Sort by eigenvalue (ascending - smallest moment = longest axis)
    idx = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    return eigenvalues, eigenvectors


def build_rotation_matrix_from_axes(
    eigenvectors: np.ndarray,
    axis_mapping: Dict[str, str]
) -> np.ndarray:
    """
    Build rotation matrix to align principal axes to Cartesian axes.

    Args:
        eigenvectors: 3x3 matrix where columns are principal axes
                     (column 0 = shortest axis, column 2 = longest axis)
        axis_mapping: Dict mapping principal axis to Cartesian axis
                     e.g., {'longest': 'X', 'middle': 'Y', 'shortest': 'Z'}

    Returns:
        3x3 rotation matrix
    """
    # Map axis names to indices
    principal_idx = {'shortest': 0, 'middle': 1, 'longest': 2}
    cartesian_vectors = {
        'X': np.array([1, 0, 0]),
        'Y': np.array([0, 1, 0]),
        'Z': np.array([0, 0, 1])
    }

    # Build target coordinate system
    target = np.zeros((3, 3))
    for principal_name, cartesian_name in axis_mapping.items():
        p_idx = principal_idx[principal_name]
        target[:, p_idx] = cartesian_vectors[cartesian_name]

    # Rotation matrix: R * eigenvectors = target
    # So R = target @ eigenvectors.T
    R = target @ eigenvectors.T

    # Ensure proper rotation (det = +1)
    if np.linalg.det(R) < 0:
        # Flip one axis to get proper rotation
        target[:, 0] = -target[:, 0]
        R = target @ eigenvectors.T

    return R


def rotation_matrix_from_vectors(vec1: np.ndarray, vec2: np.ndarray) -> np.ndarray:
    """
    Create rotation matrix that rotates vec1 to align with vec2.

    Args:
        vec1: Source vector (will be rotated)
        vec2: Target vector (to align with)

    Returns:
        3x3 rotation matrix
    """
    # Normalize vectors
    a = vec1 / np.linalg.norm(vec1)
    b = vec2 / np.linalg.norm(vec2)

    # Check if vectors are parallel or anti-parallel
    dot = np.dot(a, b)
    if np.abs(dot - 1.0) < 1e-10:
        return np.eye(3)
    if np.abs(dot + 1.0) < 1e-10:
        # 180 degree rotation - find perpendicular axis
        perp = np.array([1, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1, 0])
        perp = perp - np.dot(perp, a) * a
        perp = perp / np.linalg.norm(perp)
        # Rotation matrix for 180 degrees around perp
        return 2 * np.outer(perp, perp) - np.eye(3)

    # Rodrigues' rotation formula
    v = np.cross(a, b)
    s = np.linalg.norm(v)
    c = dot

    vx = np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

    R = np.eye(3) + vx + vx @ vx * (1 - c) / (s * s)
    return R


def euler_to_rotation_matrix(angles: Tuple[float, float, float], convention: str = 'ZYX') -> np.ndarray:
    """
    Convert Euler angles to rotation matrix.

    Args:
        angles: (alpha, beta, gamma) in degrees
        convention: Rotation order (default 'ZYX' - common for aerospace)

    Returns:
        3x3 rotation matrix
    """
    alpha, beta, gamma = np.radians(angles)

    # Individual rotation matrices
    Rx = lambda a: np.array([
        [1, 0, 0],
        [0, np.cos(a), -np.sin(a)],
        [0, np.sin(a), np.cos(a)]
    ])

    Ry = lambda a: np.array([
        [np.cos(a), 0, np.sin(a)],
        [0, 1, 0],
        [-np.sin(a), 0, np.cos(a)]
    ])

    Rz = lambda a: np.array([
        [np.cos(a), -np.sin(a), 0],
        [np.sin(a), np.cos(a), 0],
        [0, 0, 1]
    ])

    if convention == 'ZYX':
        return Rz(alpha) @ Ry(beta) @ Rx(gamma)
    elif convention == 'XYZ':
        return Rx(alpha) @ Ry(beta) @ Rz(gamma)
    elif convention == 'ZXZ':
        return Rz(alpha) @ Rx(beta) @ Rz(gamma)
    else:
        raise ValueError(f"Unknown Euler convention: {convention}")


def fit_plane_to_points(coords: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit a plane to a set of 3D points using SVD.

    Args:
        coords: Nx3 array of coordinates

    Returns:
        Tuple of (centroid, normal_vector)
    """
    centroid = np.mean(coords, axis=0)
    centered = coords - centroid

    # SVD to find the normal (smallest singular value direction)
    _, _, Vt = np.linalg.svd(centered)
    normal = Vt[-1]  # Last row = direction of smallest variance = normal

    # Ensure normal points "up" (positive Z component)
    if normal[2] < 0:
        normal = -normal

    return centroid, normal


def apply_transformation(
    pdb_path: str,
    output_path: str,
    rotation: np.ndarray,
    translation: np.ndarray
) -> Dict[str, int]:
    """
    Apply rotation and translation to a PDB file.

    Args:
        pdb_path: Input PDB file
        output_path: Output PDB file
        rotation: 3x3 rotation matrix
        translation: 3-element translation vector

    Returns:
        Statistics dict with atom counts
    """
    output_lines = []
    atom_count = 0

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                # Parse coordinates
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                # Apply transformation: R @ (coords - translation) or R @ coords + translation
                # For centering + rotation: first translate, then rotate
                coords = np.array([x, y, z])
                new_coords = rotation @ (coords + translation)

                # Rebuild line with new coordinates
                new_line = (
                    line[:30] +
                    f"{new_coords[0]:8.3f}" +
                    f"{new_coords[1]:8.3f}" +
                    f"{new_coords[2]:8.3f}" +
                    line[54:]
                )
                output_lines.append(new_line)
                atom_count += 1
            else:
                output_lines.append(line)

    # Add remark about transformation
    remark_lines = [
        "REMARK   1 STRUCTURE ORIENTED BY PROPREP\n",
        f"REMARK   1 ROTATION MATRIX:\n",
        f"REMARK   1   {rotation[0,0]:10.6f} {rotation[0,1]:10.6f} {rotation[0,2]:10.6f}\n",
        f"REMARK   1   {rotation[1,0]:10.6f} {rotation[1,1]:10.6f} {rotation[1,2]:10.6f}\n",
        f"REMARK   1   {rotation[2,0]:10.6f} {rotation[2,1]:10.6f} {rotation[2,2]:10.6f}\n",
        f"REMARK   1 TRANSLATION: {translation[0]:10.4f} {translation[1]:10.4f} {translation[2]:10.4f}\n",
    ]

    # Insert remarks after HEADER/TITLE if present, otherwise at beginning
    insert_idx = 0
    for i, line in enumerate(output_lines):
        if line.startswith(('HEADER', 'TITLE', 'COMPND', 'SOURCE', 'REMARK')):
            insert_idx = i + 1
        elif line.startswith(('ATOM', 'HETATM', 'MODEL')):
            break

    output_lines = output_lines[:insert_idx] + remark_lines + output_lines[insert_idx:]

    with open(output_path, 'w') as f:
        f.writelines(output_lines)

    return {'atoms_transformed': atom_count}


def extract_coords_from_pdb(pdb_path: str, selection: str = 'all') -> Tuple[np.ndarray, List[Dict]]:
    """
    Extract coordinates from PDB file.

    Args:
        pdb_path: Path to PDB file
        selection: 'all', 'backbone', 'CA', or specific residue selection

    Returns:
        Tuple of (Nx3 coordinates array, list of atom info dicts)
    """
    coords = []
    atom_info = []

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21]
                res_num = int(line[22:26])

                # Apply selection filter
                if selection == 'all':
                    pass
                elif selection == 'backbone':
                    if atom_name not in ('N', 'CA', 'C', 'O'):
                        continue
                elif selection == 'CA':
                    if atom_name != 'CA':
                        continue
                else:
                    # Could add more complex selection parsing here
                    pass

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                coords.append([x, y, z])
                atom_info.append({
                    'atom_name': atom_name,
                    'res_name': res_name,
                    'chain_id': chain_id,
                    'res_num': res_num
                })

    return np.array(coords), atom_info


def get_residue_center(pdb_path: str, chain: str, res_num: int) -> Optional[np.ndarray]:
    """Get the center of mass of a specific residue."""
    coords = []

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                if len(line) > 26:
                    line_chain = line[21]
                    try:
                        line_res = int(line[22:26])
                    except ValueError:
                        continue

                    if line_chain == chain and line_res == res_num:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append([x, y, z])

    if coords:
        return np.mean(coords, axis=0)
    return None


def get_atom_coords(pdb_path: str, chain: str, res_num: int, atom_name: str) -> Optional[np.ndarray]:
    """Get coordinates of a specific atom."""
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                if len(line) > 54:
                    line_atom = line[12:16].strip()
                    line_chain = line[21]
                    try:
                        line_res = int(line[22:26])
                    except ValueError:
                        continue

                    if (line_chain == chain and line_res == res_num and
                        line_atom == atom_name):
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        return np.array([x, y, z])
    return None


@register_module
class StructureOrientationModule(ProcessingModule):
    """Module for orienting structures along Cartesian axes."""

    NAME = "Structure Orientator"
    DESCRIPTION = "Align structure along Cartesian axes"
    VERSION = "1.0.0"
    CATEGORY = "preparation"
    PRIORITY = 50  # After structure loading/fixing

    @property
    def console(self):
        """Get console from processor or create new one"""
        if self.processor and hasattr(self.processor, 'console'):
            return self.processor.console
        else:
            from rich.console import Console
            return Console()

    def get_menu_options(self) -> Dict[str, str]:
        """Get module menu options."""
        return {
            "orient": "Orient structure along Cartesian axes",
        }

    def get_enhanced_menu_options(self, workspace):
        """Get menu options with enhanced status information."""
        from proprep.utils.enhanced_menu import MenuOption, OptionStatus
        from proprep.utils.structure_selector import StructureSelector

        options = []

        # Check if a structure is loaded
        selector = StructureSelector(workspace, self.console)
        has_structure = selector.get_structure(silent=True) is not None

        # Check if structure already oriented
        has_oriented = workspace.get("oriented_pdb_file") is not None

        if has_oriented:
            status = OptionStatus.COMPLETED
            desc = "Orient structure along Cartesian axes"
        elif has_structure:
            status = OptionStatus.AVAILABLE
            desc = "Orient structure along Cartesian axes"
        else:
            status = OptionStatus.BLOCKED
            desc = "Orient structure along Cartesian axes"

        options.append(MenuOption(
            key="1",
            description=desc,
            status=status,
            dependency_text="[Load a structure first] ○" if not has_structure else ""
        ))

        return options

    def handle_menu_option(self, option: str) -> bool:
        """Handle a menu option selection."""
        if option == "orient":
            return self.process(self.processor.workspace)
        return False

    def get_workspace_outputs(self) -> list:
        """Get workspace outputs"""
        return ["oriented_pdb_file", "orientation_record"]

    def availability_note(self, workspace):
        """Menu note when unavailable (○). Mirrors can_process exactly."""
        return None if self.can_process(workspace) else "Needs a loaded structure"

    def can_process(self, workspace) -> bool:
        """Check if a structure is available."""
        from proprep.utils.structure_selector import StructureSelector
        selector = StructureSelector(workspace, self.console)
        return selector.get_structure(silent=True) is not None

    def process(self, workspace) -> bool:
        """Orient structure along Cartesian axes."""
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(workspace, self.console)

        # Get input structure
        result = selector.get_structure(silent=False, return_key=True)
        if not result:
            self.console.print("[red]No structure available in workspace.[/red]")
            return False

        pdb_path, source_key = result
        self.console.print(f"\n[cyan]Input structure:[/cyan] {Path(pdb_path).name}")
        self.console.print(f"[grey50]Source: {source_key}[/grey50]")

        # Establish the input view in the viewer and drop XYZ axis
        # markers as a fixed reference frame the user can compare
        # against after the orientation runs (NGL auto-centres on
        # every structure load, so without these markers a rotated
        # structure looks identical to the original).
        self._show_orientation_view(pdb_path)

        # Show alignment options
        self.console.print("\n[bold]Alignment options:[/bold]")
        self.console.print("  1. Align principal axes to Cartesian axes (PCA)")
        self.console.print("  2. Align vector between two residues to an axis")
        self.console.print("  3. Align vector between two atoms to an axis")
        self.console.print("  4. Align plane (from residues) normal to Z")
        self.console.print("  5. Center structure at origin")
        self.console.print("  6. Apply custom rotation (Euler angles)")
        self.console.print("  7. Cancel")

        choice = prompt_with_context(
            self.processor,
            "\nSelect alignment method",
            choices=["1", "2", "3", "4", "5", "6", "7"],
            default="1",
            module=self.NAME,
            description="Select alignment method",
            options_map={
                "1": "Principal axes (PCA)",
                "2": "Two residues to axis",
                "3": "Two atoms to axis",
                "4": "Plane normal to Z",
                "5": "Center at origin",
                "6": "Custom Euler angles",
                "7": "Cancel"
            }
        )

        if choice == "7":
            self.console.print("[yellow]Cancelled.[/yellow]")
            return False

        # Dispatch to appropriate method
        method_map = {
            "1": self._align_principal_axes,
            "2": self._align_residues_to_axis,
            "3": self._align_atoms_to_axis,
            "4": self._align_plane_to_z,
            "5": self._center_at_origin,
            "6": self._apply_custom_rotation,
        }

        return method_map[choice](pdb_path, workspace)

    def _align_principal_axes(self, pdb_path: str, workspace: Dict) -> bool:
        """Align principal axes of inertia to Cartesian axes."""
        self.console.print("\n[bold]Principal Axes Alignment[/bold]")

        # Extract coordinates (use CA atoms for efficiency)
        coords, _ = extract_coords_from_pdb(pdb_path, selection='CA')

        if len(coords) < 3:
            self.console.print("[red]Need at least 3 CA atoms for principal axis alignment.[/red]")
            return False

        # Calculate principal axes
        eigenvalues, eigenvectors = get_principal_axes(coords)

        self.console.print("\n[cyan]Principal moments of inertia:[/cyan]")
        self.console.print(f"  Longest axis  (smallest moment): λ = {eigenvalues[0]:.1f}")
        self.console.print(f"  Middle axis:                      λ = {eigenvalues[1]:.1f}")
        self.console.print(f"  Shortest axis (largest moment):   λ = {eigenvalues[2]:.1f}")

        # Get user's axis mapping
        self.console.print("\n[bold]Map principal axes to Cartesian axes:[/bold]")

        longest_target = prompt_with_context(
            self.processor,
            "Align LONGEST axis to",
            choices=["X", "Y", "Z"],
            default="Z",
            module=self.NAME,
            description="Target for longest axis"
        ).upper()

        remaining = [a for a in ["X", "Y", "Z"] if a != longest_target]
        middle_target = prompt_with_context(
            self.processor,
            "Align MIDDLE axis to",
            choices=remaining,
            default=remaining[0],
            module=self.NAME,
            description="Target for middle axis"
        ).upper()

        shortest_target = [a for a in ["X", "Y", "Z"] if a not in [longest_target, middle_target]][0]
        self.console.print(f"[grey50]Shortest axis will align to {shortest_target}[/grey50]")

        axis_mapping = {
            'longest': longest_target,
            'middle': middle_target,
            'shortest': shortest_target
        }

        # Ask about centering
        center = confirm_with_context(
            self.processor,
            "Center structure at origin?",
            default=True,
            module=self.NAME,
            description="Center at origin"
        )

        # Build transformation
        com = calculate_center_of_mass(coords)
        translation = -com if center else np.zeros(3)

        # Note: eigenvectors are sorted by eigenvalue ascending
        # So column 0 = smallest eigenvalue = longest axis
        # We need to reverse the mapping
        rotation = build_rotation_matrix_from_axes(eigenvectors, {
            'shortest': axis_mapping['longest'],   # smallest eigenvalue = longest axis
            'middle': axis_mapping['middle'],
            'longest': axis_mapping['shortest']    # largest eigenvalue = shortest axis
        })

        # Apply and save
        return self._apply_and_save(
            pdb_path, workspace, rotation, translation,
            method="principal_axes",
            parameters={'axis_mapping': axis_mapping, 'centered': center}
        )

    def _align_residues_to_axis(self, pdb_path: str, workspace: Dict) -> bool:
        """Align vector between two residue centers to an axis."""
        self.console.print("\n[bold]Align Two Residues to Axis[/bold]")

        # Get first residue
        self.console.print("\n[cyan]First residue (will be at negative end of axis):[/cyan]")
        chain1 = prompt_with_context(
            self.processor, "Chain ID", default="A",
            module=self.NAME, description="First residue chain"
        ).strip()
        res1 = int_prompt_with_context(
            self.processor, "Residue number",
            module=self.NAME, description="First residue number"
        )

        # Get second residue
        self.console.print("\n[cyan]Second residue (will be at positive end of axis):[/cyan]")
        chain2 = prompt_with_context(
            self.processor, "Chain ID", default=chain1,
            module=self.NAME, description="Second residue chain"
        ).strip()
        res2 = int_prompt_with_context(
            self.processor, "Residue number",
            module=self.NAME, description="Second residue number"
        )

        # Get residue centers
        center1 = get_residue_center(pdb_path, chain1, res1)
        center2 = get_residue_center(pdb_path, chain2, res2)

        if center1 is None:
            self.console.print(f"[red]Residue {chain1}:{res1} not found.[/red]")
            return False
        if center2 is None:
            self.console.print(f"[red]Residue {chain2}:{res2} not found.[/red]")
            return False

        # Calculate vector
        vec = center2 - center1
        distance = np.linalg.norm(vec)
        self.console.print(f"\n[grey50]Vector length: {distance:.2f} Å[/grey50]")

        # Visual confirmation: residue 1 (negative end) red, residue 2
        # (positive end) green. Lets the user verify their pick before
        # committing to an axis.
        self._highlight_endpoint_pair(
            f":{chain1} and {res1}",
            f":{chain2} and {res2}",
        )

        # Get target axis
        target_axis = prompt_with_context(
            self.processor,
            "Align this vector to axis",
            choices=["X", "Y", "Z"],
            default="Z",
            module=self.NAME,
            description="Target axis"
        ).upper()

        # Ask about centering
        center = confirm_with_context(
            self.processor,
            "Center structure at origin?",
            default=True,
            module=self.NAME,
            description="Center at origin"
        )

        # Build transformation
        target_vec = {'X': np.array([1,0,0]), 'Y': np.array([0,1,0]), 'Z': np.array([0,0,1])}[target_axis]
        rotation = rotation_matrix_from_vectors(vec, target_vec)

        if center:
            coords, _ = extract_coords_from_pdb(pdb_path, selection='all')
            com = calculate_center_of_mass(coords)
            translation = -com
        else:
            translation = np.zeros(3)

        return self._apply_and_save(
            pdb_path, workspace, rotation, translation,
            method="residues_to_axis",
            parameters={
                'residue1': f"{chain1}:{res1}",
                'residue2': f"{chain2}:{res2}",
                'target_axis': target_axis,
                'centered': center
            }
        )

    def _align_atoms_to_axis(self, pdb_path: str, workspace: Dict) -> bool:
        """Align vector between two specific atoms to an axis."""
        self.console.print("\n[bold]Align Two Atoms to Axis[/bold]")

        # Get first atom
        self.console.print("\n[cyan]First atom (will be at negative end of axis):[/cyan]")
        chain1 = prompt_with_context(
            self.processor, "Chain ID", default="A",
            module=self.NAME, description="First atom chain"
        ).strip()
        res1 = int_prompt_with_context(
            self.processor, "Residue number",
            module=self.NAME, description="First atom residue"
        )
        atom1 = prompt_with_context(
            self.processor, "Atom name (e.g., CA, FE, ZN)", default="CA",
            module=self.NAME, description="First atom name"
        ).strip().upper()

        # Get second atom
        self.console.print("\n[cyan]Second atom (will be at positive end of axis):[/cyan]")
        chain2 = prompt_with_context(
            self.processor, "Chain ID", default=chain1,
            module=self.NAME, description="Second atom chain"
        ).strip()
        res2 = int_prompt_with_context(
            self.processor, "Residue number",
            module=self.NAME, description="Second atom residue"
        )
        atom2 = prompt_with_context(
            self.processor, "Atom name", default=atom1,
            module=self.NAME, description="Second atom name"
        ).strip().upper()

        # Get atom coordinates
        coord1 = get_atom_coords(pdb_path, chain1, res1, atom1)
        coord2 = get_atom_coords(pdb_path, chain2, res2, atom2)

        if coord1 is None:
            self.console.print(f"[red]Atom {chain1}:{res1}:{atom1} not found.[/red]")
            return False
        if coord2 is None:
            self.console.print(f"[red]Atom {chain2}:{res2}:{atom2} not found.[/red]")
            return False

        # Calculate vector
        vec = coord2 - coord1
        distance = np.linalg.norm(vec)
        self.console.print(f"\n[grey50]Atom-atom distance: {distance:.2f} Å[/grey50]")

        # Visual confirmation: atom 1 (negative end) red, atom 2
        # (positive end) green. Selection includes the atom-name
        # filter so the marker pinpoints the specific atom rather
        # than the whole residue.
        self._highlight_endpoint_pair(
            f":{chain1} and {res1} and .{atom1}",
            f":{chain2} and {res2} and .{atom2}",
        )

        # Get target axis
        target_axis = prompt_with_context(
            self.processor,
            "Align this vector to axis",
            choices=["X", "Y", "Z"],
            default="Z",
            module=self.NAME,
            description="Target axis"
        ).upper()

        # Ask about centering
        center = confirm_with_context(
            self.processor,
            "Center structure at origin?",
            default=True,
            module=self.NAME,
            description="Center at origin"
        )

        # Build transformation
        target_vec = {'X': np.array([1,0,0]), 'Y': np.array([0,1,0]), 'Z': np.array([0,0,1])}[target_axis]
        rotation = rotation_matrix_from_vectors(vec, target_vec)

        if center:
            coords, _ = extract_coords_from_pdb(pdb_path, selection='all')
            com = calculate_center_of_mass(coords)
            translation = -com
        else:
            translation = np.zeros(3)

        return self._apply_and_save(
            pdb_path, workspace, rotation, translation,
            method="atoms_to_axis",
            parameters={
                'atom1': f"{chain1}:{res1}:{atom1}",
                'atom2': f"{chain2}:{res2}:{atom2}",
                'target_axis': target_axis,
                'centered': center
            }
        )

    def _align_plane_to_z(self, pdb_path: str, workspace: Dict) -> bool:
        """Fit plane to selected residues and align normal to Z."""
        self.console.print("\n[bold]Align Plane Normal to Z[/bold]")
        self.console.print("[grey50]Select residues that define a plane (e.g., membrane-facing residues)[/grey50]")

        # Get residue selection
        self.console.print("\nEnter residues (format: chain:resnum, comma-separated)")
        self.console.print("[grey50]Example: A:10, A:20, A:30, A:40[/grey50]")

        residue_str = prompt_with_context(
            self.processor,
            "Residues",
            module=self.NAME,
            description="Residues defining plane"
        )

        # Parse residues and get coordinates
        plane_coords = []
        accepted_specs: List[str] = []
        for res_spec in residue_str.split(','):
            res_spec = res_spec.strip()
            if ':' in res_spec:
                chain, res_num = res_spec.split(':')
                try:
                    center = get_residue_center(pdb_path, chain.strip(), int(res_num.strip()))
                    if center is not None:
                        plane_coords.append(center)
                        accepted_specs.append(f"{chain.strip()}:{res_num.strip()}")
                except ValueError:
                    continue

        if len(plane_coords) < 3:
            self.console.print("[red]Need at least 3 residues to define a plane.[/red]")
            return False

        plane_coords = np.array(plane_coords)

        # Fit plane
        centroid, normal = fit_plane_to_points(plane_coords)
        self.console.print(f"\n[grey50]Plane fitted to {len(plane_coords)} residues[/grey50]")
        self.console.print(f"[grey50]Normal vector: ({normal[0]:.3f}, {normal[1]:.3f}, {normal[2]:.3f})[/grey50]")

        # Highlight the plane-defining residues so the user can verify
        # they've picked the right set before centering. Palette per
        # residue so each is independently identifiable in the rep
        # manager.
        self._highlight_plane_residues(accepted_specs)

        # Ask about centering
        center = confirm_with_context(
            self.processor,
            "Center structure at origin?",
            default=True,
            module=self.NAME,
            description="Center at origin"
        )

        # Build rotation to align normal with Z
        z_axis = np.array([0, 0, 1])
        rotation = rotation_matrix_from_vectors(normal, z_axis)

        if center:
            coords, _ = extract_coords_from_pdb(pdb_path, selection='all')
            com = calculate_center_of_mass(coords)
            translation = -com
        else:
            translation = np.zeros(3)

        return self._apply_and_save(
            pdb_path, workspace, rotation, translation,
            method="plane_to_z",
            parameters={
                'residues': residue_str,
                'normal': normal.tolist(),
                'centered': center
            }
        )

    def _center_at_origin(self, pdb_path: str, workspace: Dict) -> bool:
        """Center structure at origin without rotation."""
        self.console.print("\n[bold]Center at Origin[/bold]")

        # Get coordinates
        coords, _ = extract_coords_from_pdb(pdb_path, selection='all')

        # Calculate center
        use_mass = confirm_with_context(
            self.processor,
            "Use center of mass (vs geometric center)?",
            default=True,
            module=self.NAME,
            description="Use center of mass"
        )

        if use_mass:
            # Approximate masses by atom type would be better, but geometric is fine for now
            com = np.mean(coords, axis=0)  # Geometric for simplicity
            self.console.print("[grey50]Using geometric center (mass-weighted requires atom masses)[/grey50]")
        else:
            com = np.mean(coords, axis=0)

        self.console.print(f"\n[grey50]Current center: ({com[0]:.2f}, {com[1]:.2f}, {com[2]:.2f})[/grey50]")

        rotation = np.eye(3)
        translation = -com

        return self._apply_and_save(
            pdb_path, workspace, rotation, translation,
            method="center_origin",
            parameters={'original_center': com.tolist()}
        )

    def _apply_custom_rotation(self, pdb_path: str, workspace: Dict) -> bool:
        """Apply custom Euler angle rotation."""
        self.console.print("\n[bold]Custom Euler Angle Rotation[/bold]")

        # Get Euler angles
        self.console.print("\nEnter Euler angles in degrees (ZYX convention):")

        alpha = float_prompt_with_context(
            self.processor, "Rotation around Z (alpha)", default=0.0,
            module=self.NAME, description="Z rotation angle"
        )
        beta = float_prompt_with_context(
            self.processor, "Rotation around Y (beta)", default=0.0,
            module=self.NAME, description="Y rotation angle"
        )
        gamma = float_prompt_with_context(
            self.processor, "Rotation around X (gamma)", default=0.0,
            module=self.NAME, description="X rotation angle"
        )

        # Ask about centering
        center = confirm_with_context(
            self.processor,
            "Center structure at origin before rotation?",
            default=True,
            module=self.NAME,
            description="Center before rotation"
        )

        rotation = euler_to_rotation_matrix((alpha, beta, gamma), 'ZYX')

        if center:
            coords, _ = extract_coords_from_pdb(pdb_path, selection='all')
            com = calculate_center_of_mass(coords)
            translation = -com
        else:
            translation = np.zeros(3)

        return self._apply_and_save(
            pdb_path, workspace, rotation, translation,
            method="custom_euler",
            parameters={
                'alpha': alpha,
                'beta': beta,
                'gamma': gamma,
                'convention': 'ZYX',
                'centered': center
            }
        )

    def _apply_and_save(
        self,
        pdb_path: str,
        workspace: Dict,
        rotation: np.ndarray,
        translation: np.ndarray,
        method: str,
        parameters: Dict
    ) -> bool:
        """Apply transformation and save to workspace."""

        # Generate output filename
        input_name = Path(pdb_path).stem
        output_name = f"{input_name}_oriented.pdb"
        output_path = Path(pdb_path).parent / output_name

        # Check for existing file
        if output_path.exists():
            if not confirm_with_context(
                self.processor,
                f"{output_name} already exists. Overwrite?",
                default=True,
                module=self.NAME,
                description="Overwrite existing file"
            ):
                # Generate unique name
                i = 1
                while output_path.exists():
                    output_name = f"{input_name}_oriented_{i}.pdb"
                    output_path = Path(pdb_path).parent / output_name
                    i += 1

        # Apply transformation
        self.console.print(f"\n[cyan]Applying transformation...[/cyan]")

        try:
            stats = apply_transformation(
                pdb_path=pdb_path,
                output_path=str(output_path),
                rotation=rotation,
                translation=translation
            )

            self.console.print(f"[green]✓ Transformed {stats['atoms_transformed']:,} atoms[/green]")
            self.console.print(f"[green]✓ Saved: {output_path}[/green]")

            # Save transformation record
            record = TransformationRecord(
                method=method,
                rotation_matrix=rotation,
                translation_vector=translation,
                parameters=parameters,
                timestamp=datetime.now().isoformat()
            )

            # Save record to JSON file
            record_path = output_path.with_suffix('.orientation.json')
            with open(record_path, 'w') as f:
                json.dump(record.to_dict(), f, indent=2)
            self.console.print(f"[grey50]Transformation saved: {record_path.name}[/grey50]")

            # Save to workspace
            workspace.set("oriented_pdb_file", str(output_path))
            workspace.set("orientation_record", record.to_dict())
            self.console.print("[green]✓ Registered in workspace as 'oriented_pdb_file'[/green]")

            # Show the oriented structure with the same XYZ axis
            # markers from Hook 1, so the user can compare orientations
            # against the same fixed reference frame. Per-method
            # endpoint/plane labels are dropped first because the
            # show_structure swap re-anchors the scene around the
            # rotated coordinates and those highlights would land on
            # the same residue numbers but in different positions.
            self._clear_per_method_orientation_overlays()
            # Auto-fired post-orientation refresh. Default force=False on the
            # helper keeps it silent in CLI when no viewer is open; pushes
            # live to one that is. Web shell auto-launches via env flag.
            self._show_orientation_view(str(output_path))

            return True

        except Exception as e:
            self.console.print(f"[red]Error applying transformation: {e}[/red]")
            logger.exception("Error in structure orientation")
            return False

    # ── Viewer hooks ─────────────────────────────────────────────────────

    # Standard X/Y/Z colour convention (red/green/blue) at 15 Å from
    # origin, plus a white "you are here" marker at the origin itself
    # so the user can see where the protein lands when centering is
    # enabled. Stable labels so re-firing replaces.
    _AXIS_MARKERS = (
        ((0.0, 0.0, 0.0), 1.2, "#ffffff", "orient_axis_origin"),
        ((15.0, 0.0, 0.0), 2.0, "#e31a1c", "orient_axis_X"),
        ((0.0, 15.0, 0.0), 2.0, "#33a02c", "orient_axis_Y"),
        ((0.0, 0.0, 15.0), 2.0, "#1f78b4", "orient_axis_Z"),
    )

    def _show_orientation_view(self, pdb_path: str, *, force: bool = False) -> None:
        """Show the structure plus XYZ axis sphere markers.

        Used at module entry (baseline view) and after a successful
        transformation (oriented view), so the user can see the
        orientation change against a fixed reference frame. NGL's
        auto-centre on structure load otherwise hides any rotation.
        """
        if not pdb_path:
            return
        try:
            from proprep.structure_prep.viewer_coordinator import (
                viewer as _viewer,
            )
            _viewer.show_structure(pdb_path, force=force)
            for coords, radius, color, label in self._AXIS_MARKERS:
                _viewer.show_sphere(
                    coords, radius=radius, label=label,
                    color=color, opacity=0.9, force=force,
                )
        except Exception as exc:
            logger.debug("orientation viewer hook silenced: %s", exc)

    def _highlight_endpoint_pair(
        self,
        neg_selection: str,
        pos_selection: str,
    ) -> None:
        """Highlight a vector's two endpoints — negative end red,
        positive end green — matching the prompt copy that says the
        first pick is the "negative end" and the second is the
        "positive end".
        """
        try:
            from proprep.structure_prep.viewer_coordinator import (
                viewer as _viewer,
            )
            _viewer.unhighlight("orient_endpoint_neg")
            _viewer.unhighlight("orient_endpoint_pos")
            _viewer.highlight(
                neg_selection, style="ball+stick",
                color="#e31a1c", label="orient_endpoint_neg",
            )
            _viewer.highlight(
                pos_selection, style="ball+stick",
                color="#33a02c", label="orient_endpoint_pos",
            )
        except Exception as exc:
            logger.debug("endpoint viewer hook silenced: %s", exc)

    def _highlight_plane_residues(self, residue_specs: List[str]) -> None:
        """Highlight a list of plane-defining residues, palette-coloured
        by index so each residue is independently identifiable.

        ``residue_specs`` is a list of ``"chain:resnum"`` strings
        already validated by the caller.
        """
        try:
            from proprep.structure_prep.viewer_coordinator import (
                viewer as _viewer,
            )
            applied: List[str] = []
            for idx, spec in enumerate(residue_specs, 1):
                if ':' not in spec:
                    continue
                chain, res_num = spec.split(':', 1)
                chain = chain.strip()
                res_num = res_num.strip()
                if not chain or not res_num:
                    continue
                label = f"orient_plane_{idx}_{chain}_{res_num}"
                _viewer.highlight(
                    f":{chain} and {res_num}",
                    style="ball+stick",
                    color=f"palette:{idx}",
                    label=label,
                )
                applied.append(label)
            self._plane_viewer_labels = applied
        except Exception as exc:
            logger.debug("plane viewer hook silenced: %s", exc)
            self._plane_viewer_labels = []

    def _clear_per_method_orientation_overlays(self) -> None:
        """Drop endpoint and plane labels left behind by Hooks 3/4/5
        before swapping to the oriented structure. Axis spheres are
        kept across the swap (re-asserted by _show_orientation_view).
        """
        try:
            from proprep.structure_prep.viewer_coordinator import (
                viewer as _viewer,
            )
            _viewer.unhighlight("orient_endpoint_neg")
            _viewer.unhighlight("orient_endpoint_pos")
            for lbl in getattr(self, "_plane_viewer_labels", None) or []:
                _viewer.unhighlight(lbl)
            self._plane_viewer_labels = []
        except Exception as exc:
            logger.debug("orientation overlay cleanup silenced: %s", exc)
