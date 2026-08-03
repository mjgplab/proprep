"""
ONIOM Data Structures

Shared data structures for ONIOM QM/MM setup.
All structures use parmed integer atom indices as primary keys.
Coordinates are stored as metadata for Gaussian output.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple
import math

from proprep.structure_prep.comprehensive_redox_detector import RedoxSite


# ===== ENUMS =====

class ONIOMLayer(Enum):
    """ONIOM layer designations."""
    HIGH = "H"      # QM region (highest level of theory)
    MEDIUM = "M"    # Optional intermediate layer (2-layer ONIOM doesn't use this)
    LOW = "L"       # MM region (molecular mechanics)


class FreezeFlag(Enum):
    """ONIOM freeze/active flags for geometry optimization."""
    ACTIVE = 0      # Atom can move during optimization
    FROZEN = -1     # Atom is frozen in place


# ===== DATA CLASSES =====

@dataclass
class LayerAssignment:
    """Layer assignment for one atom in ONIOM setup."""

    # Primary identifier (parmed atom index, 0-based)
    atom_idx: int

    # Coordinates for Gaussian output
    coords: Tuple[float, float, float]

    # Layer assignment
    layer: ONIOMLayer              # HIGH, MEDIUM, or LOW
    freeze: FreezeFlag             # ACTIVE or FROZEN

    # Residue metadata (from parmed)
    residue_idx: int               # parmed residue index (0-based)
    residue_name: str              # e.g. "ALA", "HEM"
    residue_number: int            # parmed residue.number (sequential after tLEaP)

    # Atom metadata
    atom_name: str
    element: str

    # Assignment reason (for user reporting)
    assignment_reason: str         # e.g., "Selected residue", "Alpha carbon bridge"


@dataclass
class LinkAtom:
    """Hydrogen link atom at QM/MM boundary."""

    # Parent atom indices (parmed, 0-based)
    qm_parent_idx: int             # QM atom (closer to QM region)
    mm_parent_idx: int             # MM atom (in MM region)

    # Link atom position and geometry (for Gaussian output)
    coords: Tuple[float, float, float]              # Calculated position
    qm_parent_coords: Tuple[float, float, float]    # QM parent position
    mm_parent_coords: Tuple[float, float, float]    # MM parent position
    bond_vector: Tuple[float, float, float]         # Unit vector from QM to MM parent
    bond_length: float                              # Original bond length (Angstroms)

    # Optional fields with defaults
    element: str = "H"                              # Always hydrogen for link atoms
    scale_factor: float = 0.723                     # Default C-H scale factor
    boundary_type: str = ""                         # e.g., "Backbone bond", "Sidechain bond"


@dataclass
class ConnectivityEntry:
    """One atom's connectivity information for Gaussian connectivity table."""

    atom_index: int                        # 1-based index in ONIOM atom list
    connected_indices: List[int] = field(default_factory=list)  # 1-based indices of bonded atoms
    bond_orders: List[float] = field(default_factory=list)      # Bond order for each connection

    # Optional: store element for validation
    element: str = ""


@dataclass
class ONIOMSetup:
    """Complete ONIOM calculation configuration."""

    # Source data
    redox_site: RedoxSite                  # Original RedoxSite
    parm: Any = None                       # parmed.Structure for downstream access

    # Layer assignments (keyed by parmed atom index, 0-based)
    layer_assignments: Dict[int, LayerAssignment] = field(default_factory=dict)

    # Link atoms
    link_atoms: List[LinkAtom] = field(default_factory=list)

    # Connectivity table (OPTIONAL — Gaussian can auto-generate)
    connectivity: List[ConnectivityEntry] = field(default_factory=list)
    use_explicit_connectivity: bool = False

    # Atom types and charges (keyed by parmed atom index, 0-based)
    atom_types: Dict[int, str] = field(default_factory=dict)
    charges: Dict[int, float] = field(default_factory=dict)

    # QM settings
    qm_functional: str = "B3LYP"
    qm_basis_set: str = "6-31G*"
    qm_charge: int = 0
    qm_multiplicity: int = 1

    # MM settings
    mm_forcefield: str = "AMBER"

    # Medium layer settings (for 3-layer ONIOM)
    medium_method: str = "HF/3-21G"
    medium_charge: int = 0
    medium_multiplicity: int = 1

    # Job settings
    job_type: str = "Opt"
    n_processors: int = 4
    memory_gb: int = 8
    additional_keywords: str = ""

    # Metadata
    n_layers: int = 2
    setup_timestamp: str = ""
    validation_passed: bool = False
    validation_warnings: List[str] = field(default_factory=list)


@dataclass
class LayerStatistics:
    """Summary statistics for user reporting."""

    # Atom counts by layer
    n_high: int = 0
    n_medium: int = 0
    n_low: int = 0
    n_link: int = 0

    # Freeze counts
    n_frozen: int = 0
    n_active: int = 0

    # Residue counts by layer (parmed residue indices)
    residues_high: List[int] = field(default_factory=list)
    residues_medium: List[int] = field(default_factory=list)
    residues_low: List[int] = field(default_factory=list)

    # Charge information
    total_charge_high: float = 0.0
    total_charge_medium: float = 0.0
    total_charge_low: float = 0.0


# ===== UTILITY FUNCTIONS =====

def calculate_layer_statistics(oniom_setup: ONIOMSetup) -> LayerStatistics:
    """
    Calculate statistics from ONIOMSetup for user reporting.

    Args:
        oniom_setup: Complete ONIOM setup

    Returns:
        LayerStatistics with counts and charge info
    """
    stats = LayerStatistics()

    residues_seen_high: set = set()
    residues_seen_medium: set = set()
    residues_seen_low: set = set()

    for atom_idx, assignment in oniom_setup.layer_assignments.items():
        if assignment.layer == ONIOMLayer.HIGH:
            stats.n_high += 1
            stats.total_charge_high += oniom_setup.charges.get(atom_idx, 0.0)
            if assignment.residue_idx not in residues_seen_high:
                residues_seen_high.add(assignment.residue_idx)
                stats.residues_high.append(assignment.residue_idx)

        elif assignment.layer == ONIOMLayer.MEDIUM:
            stats.n_medium += 1
            stats.total_charge_medium += oniom_setup.charges.get(atom_idx, 0.0)
            if assignment.residue_idx not in residues_seen_medium:
                residues_seen_medium.add(assignment.residue_idx)
                stats.residues_medium.append(assignment.residue_idx)

        elif assignment.layer == ONIOMLayer.LOW:
            stats.n_low += 1
            stats.total_charge_low += oniom_setup.charges.get(atom_idx, 0.0)
            if assignment.residue_idx not in residues_seen_low:
                residues_seen_low.add(assignment.residue_idx)
                stats.residues_low.append(assignment.residue_idx)

        if assignment.freeze == FreezeFlag.FROZEN:
            stats.n_frozen += 1
        else:
            stats.n_active += 1

    stats.n_link = len(oniom_setup.link_atoms)

    return stats
