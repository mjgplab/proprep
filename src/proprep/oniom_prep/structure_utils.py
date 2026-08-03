"""
Structure Utilities for ONIOM Preparation

Helper functions for PDB file selection, distance calculations, and other utilities.

© 2024 ProPrep Developer. All rights reserved.
"""

import os
import math
import logging
from typing import Optional, Tuple, List
from rich.console import Console


logger = logging.getLogger(__name__)


def get_priority_pdb_file(workspace, console: Optional[Console] = None, processor=None) -> Optional[str]:
    """
    Select the appropriate PDB file with priority ranking.

    Priority order (highest to lowest):
    1. protonation-updated (protonation_pdb_file)
    2. transformed (transformed_pdb_file)
    3. repaired (repaired_pdb_file)
    4. filtered (filtered_pdb_file)
    5. Interactive selection from available structures

    Args:
        workspace: ProPrep workspace object
        console: Optional Rich console for output
        processor: Optional processor for session recording

    Returns:
        Path to selected PDB file, or None if no valid file found
    """
    if console is None:
        console = Console()

    # Priority order matches other modules
    priority_keys = [
        ("protonation_pdb_file", "protonation-updated"),
        ("structure_with_prot_resnames", "protonation-updated"),  # Legacy key
        ("transformed_pdb_file", "transformed"),
        ("repaired_pdb_file", "repaired"),
        ("filtered_pdb_file", "filtered"),
    ]

    for workspace_key, description in priority_keys:
        pdb_file = workspace.get(workspace_key)
        if pdb_file and os.path.exists(pdb_file):
            console.print(f"[green]Using {description} PDB file: {pdb_file}[/green]")
            logger.info(f"Selected {description} PDB file: {pdb_file}")
            return pdb_file

    # No processed structure found - ask user to select interactively
    from proprep.utils.structure_selector import get_interactive_pdb_file

    pdb_file, file_key = get_interactive_pdb_file(
        workspace,
        console=console,
        processor=processor,
        return_key=True
    )

    if pdb_file:
        # Get human-readable source name
        source_mapping = {
            "rcsb_pdb_file": "RCSB PDB",
            "local_pdb_file": "local",
            "alphafold_pdb_file": "AlphaFold",
            "alphafold_homolog_pdb_file": "AlphaFold homolog",
            "aligned_target_pdb_file": "aligned target",
            "aligned_ref_pdb_file": "aligned reference"
        }
        source = source_mapping.get(file_key, "selected")
        console.print(f"[green]Using {source} PDB file: {pdb_file}[/green]")
        logger.info(f"Selected {source} PDB file: {pdb_file}")
        return pdb_file

    # No valid structure found
    console.print("[red]No valid PDB file found in workspace[/red]")
    console.print("[yellow]Please load a PDB file first[/yellow]")
    logger.error("No valid PDB file found in workspace")
    return None


def calculate_distance(coords1: Tuple[float, float, float],
                      coords2: Tuple[float, float, float]) -> float:
    """
    Calculate Euclidean distance between two points.

    Args:
        coords1: First point (x, y, z)
        coords2: Second point (x, y, z)

    Returns:
        Distance in Angstroms
    """
    dx = coords1[0] - coords2[0]
    dy = coords1[1] - coords2[1]
    dz = coords1[2] - coords2[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def calculate_geometric_center(coords_list: List[Tuple[float, float, float]]) -> Tuple[float, float, float]:
    """
    Calculate geometric center of a list of coordinates.

    Args:
        coords_list: List of (x, y, z) tuples

    Returns:
        (x, y, z) tuple of geometric center

    Raises:
        ValueError: If coords_list is empty
    """
    if not coords_list:
        raise ValueError("coords_list cannot be empty")

    n = len(coords_list)
    center_x = sum(c[0] for c in coords_list) / n
    center_y = sum(c[1] for c in coords_list) / n
    center_z = sum(c[2] for c in coords_list) / n

    return (center_x, center_y, center_z)


def normalize_vector(vector: Tuple[float, float, float]) -> Tuple[float, float, float]:
    """
    Normalize a 3D vector to unit length.

    Args:
        vector: (x, y, z) vector

    Returns:
        Normalized (x, y, z) unit vector

    Raises:
        ValueError: If vector has zero length
    """
    length = math.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)

    if length < 1e-10:
        raise ValueError("Cannot normalize zero-length vector")

    return (vector[0] / length, vector[1] / length, vector[2] / length)
