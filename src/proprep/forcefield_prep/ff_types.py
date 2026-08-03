"""
Unified Force Field Type Definitions

This module provides the SINGLE source of truth for all force field-related
dataclasses used throughout the metal site parameterization workflow.

These types replace the previously inconsistent mix of:
- Plain dicts with string keys
- AtomTypeAssignment from mcpb/atom_typer.py
- GlobalTypeAssignment from mcpb/global_atom_registry.py
- Various parameter dicts

All code should import types from this module.

Author: ProPrep Development Team
"""

from dataclasses import dataclass, field
from typing import Optional, Tuple, List, Dict, Any
from enum import Enum


# =============================================================================
# Enumerations
# =============================================================================

class TerminalType(Enum):
    """Classification of residue position in protein chain."""
    N_TERMINAL = "N_TERMINAL"
    C_TERMINAL = "C_TERMINAL"
    INTERNAL = "INTERNAL"


class AtomSource(Enum):
    """Source of atom type assignment."""
    FORCE_FIELD = "force_field"        # From lib file (ff14SB, conste, etc.)
    METAL_DATABASE = "metal_database"  # From MetalIonDatabase
    ANTECHAMBER = "antechamber"        # From antechamber/GAFF2
    MANUAL = "manual"                  # User-entered
    GLOBAL_REGISTRY = "global_registry"  # From GlobalAtomTypeRegistry (legacy)


# =============================================================================
# Atom Type Assignment
# =============================================================================

@dataclass
class AtomTypeAssignment:
    """
    Unified atom type assignment - THE ONLY representation used everywhere.

    This dataclass captures everything needed to describe an atom's type
    assignment during metal site parameterization:

    - Identity: coordinates, serial number, residue info
    - Type assignment: original FF type and MCPB-renamed type (M1, Y1, etc.)
    - Classification: is it a metal center? coordinating ligand?
    - Parameters: charge, mass, spin state
    - Selection: can be excluded from QM model
    - Hydrogen: tracking for added capping hydrogens

    Usage:
        # Create assignment
        assignment = AtomTypeAssignment(
            coords=(1.234, 5.678, 9.012),
            chain='A',
            resname='HIS',
            resid=93,
            atom_name='NE2',
            element='N',
            original_type='NB',
            renamed_type='Y1',
            is_metal_ligand=True,
            charge=-0.5683,
            source=AtomSource.FORCE_FIELD
        )

        # Access attributes
        print(assignment.renamed_type)  # 'Y1'
        print(assignment.is_metal_ligand)  # True
    """
    # === Identity ===
    coords: Tuple[float, float, float]
    serial_number: Optional[int] = None  # PDB serial number (for fingerprint matching)

    # === Residue Info ===
    chain: str = ""
    resname: str = ""          # Original PDB residue name (e.g., "HIS", "HEM")
    ff_resname: str = ""       # Force field residue name (e.g., "HID", "HEH")
    resid: int = 0
    atom_name: str = ""
    element: str = ""

    # === Type Assignment ===
    original_type: str = ""    # From force field (e.g., "NB", "FE", "CA")
    renamed_type: str = ""     # After MCPB renaming (e.g., "Y1", "M1", or unchanged)

    # === Classification ===
    is_center: bool = False           # Is metal center?
    is_metal_ligand: bool = False     # Directly coordinates to metal?
    terminal_type: TerminalType = TerminalType.INTERNAL

    # === Parameters ===
    charge: Optional[float] = None
    mass: Optional[float] = None
    spin: Optional[int] = None        # For metals: unpaired electrons

    # === VDW Parameters (for metals) ===
    vdw_radius: Optional[float] = None
    vdw_epsilon: Optional[float] = None
    vdw_source: str = ""

    # === Source Tracking ===
    source: AtomSource = AtomSource.FORCE_FIELD
    source_detail: str = ""    # e.g., "conste.lib", "gaff2", "manual_entry"
    library_source: str = ""   # Legacy compatibility

    # === Atom Selection (for QM model) ===
    include_in_qm: bool = True
    exclude_reason: str = ""   # Why excluded (e.g., "propionate - separate residue")

    # === Hydrogen Addition ===
    needs_hydrogen: bool = False
    added_hydrogens: List[Tuple[float, float, float]] = field(default_factory=list)

    # === Legacy Compatibility ===
    renamed: bool = False      # True if renamed_type differs from original_type

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert to dictionary for JSON serialization.

        Returns:
            Dict suitable for json.dump()
        """
        return {
            'coords': list(self.coords),
            'serial_number': self.serial_number,
            'chain': self.chain,
            'resname': self.resname,
            'ff_resname': self.ff_resname,
            'resid': self.resid,
            'atom_name': self.atom_name,
            'element': self.element,
            'original_type': self.original_type,
            'renamed_type': self.renamed_type,
            'is_center': self.is_center,
            'is_metal_ligand': self.is_metal_ligand,
            'terminal_type': self.terminal_type.value,
            'charge': self.charge,
            'mass': self.mass,
            'spin': self.spin,
            'vdw_radius': self.vdw_radius,
            'vdw_epsilon': self.vdw_epsilon,
            'vdw_source': self.vdw_source,
            'source': self.source.value,
            'source_detail': self.source_detail,
            'library_source': self.library_source,
            'include_in_qm': self.include_in_qm,
            'exclude_reason': self.exclude_reason,
            'needs_hydrogen': self.needs_hydrogen,
            'added_hydrogens': [list(h) for h in self.added_hydrogens],
            'renamed': self.renamed
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'AtomTypeAssignment':
        """
        Create from dictionary (JSON deserialization).

        Args:
            data: Dict from json.load()

        Returns:
            AtomTypeAssignment instance
        """
        # Handle coords
        coords = tuple(data.get('coords', (0.0, 0.0, 0.0)))

        # Handle enums
        terminal_type = TerminalType(data.get('terminal_type', 'INTERNAL'))
        try:
            source = AtomSource(data.get('source', 'force_field'))
        except ValueError:
            source = AtomSource.FORCE_FIELD

        # Handle added_hydrogens
        added_h = [tuple(h) for h in data.get('added_hydrogens', [])]

        return cls(
            coords=coords,
            serial_number=data.get('serial_number'),
            chain=data.get('chain', ''),
            resname=data.get('resname', ''),
            ff_resname=data.get('ff_resname', ''),
            resid=data.get('resid', 0),
            atom_name=data.get('atom_name', ''),
            element=data.get('element', ''),
            original_type=data.get('original_type', ''),
            renamed_type=data.get('renamed_type', ''),
            is_center=data.get('is_center', False),
            is_metal_ligand=data.get('is_metal_ligand', False),
            terminal_type=terminal_type,
            charge=data.get('charge'),
            mass=data.get('mass'),
            spin=data.get('spin'),
            vdw_radius=data.get('vdw_radius'),
            vdw_epsilon=data.get('vdw_epsilon'),
            vdw_source=data.get('vdw_source', ''),
            source=source,
            source_detail=data.get('source_detail', ''),
            library_source=data.get('library_source', ''),
            include_in_qm=data.get('include_in_qm', True),
            exclude_reason=data.get('exclude_reason', ''),
            needs_hydrogen=data.get('needs_hydrogen', False),
            added_hydrogens=added_h,
            renamed=data.get('renamed', False)
        )


# =============================================================================
# Force Field Parameters
# =============================================================================

@dataclass
class MassParameter:
    """Atomic mass parameter from parmXX.dat."""
    atom_type: str
    mass: float
    polarizability: float = 0.0
    comment: str = ""
    source: str = ""  # e.g., "parm10.dat"


@dataclass
class BondParameter:
    """
    Bond stretching parameter.

    The key property provides a canonical lookup key that is
    order-independent (sorted alphabetically).
    """
    type1: str
    type2: str
    force_constant: float  # kcal/mol/Å²
    eq_length: float       # Å
    source: str = ""       # e.g., "parm10.dat", "frcmod.conste"

    @property
    def key(self) -> str:
        """Canonical key for lookup (alphabetically sorted)."""
        types = sorted([self.type1, self.type2])
        return f"{types[0]}-{types[1]}"


@dataclass
class AngleParameter:
    """
    Angle bending parameter.

    The key property provides a canonical lookup key that handles
    the symmetry of angles (A-B-C == C-B-A).
    """
    type1: str
    type2: str  # Central atom
    type3: str
    force_constant: float  # kcal/mol/rad²
    eq_angle: float        # degrees
    source: str = ""

    @property
    def key(self) -> str:
        """Canonical key (reversible: A-B-C == C-B-A)."""
        if self.type1 <= self.type3:
            return f"{self.type1}-{self.type2}-{self.type3}"
        return f"{self.type3}-{self.type2}-{self.type1}"


@dataclass
class DihedralParameter:
    """
    Dihedral (torsion) parameter.

    Dihedrals can have multiple terms with different periodicities.
    """
    type1: str
    type2: str
    type3: str
    type4: str
    divider: int           # IDIVF - divide barrier by this
    barrier: float         # PK - barrier height in kcal/mol
    phase: float           # PHASE - phase angle in degrees
    periodicity: float     # PN - periodicity (negative = more terms follow)
    source: str = ""

    @property
    def key(self) -> str:
        """Canonical key (reversible: A-B-C-D == D-C-B-A)."""
        forward = f"{self.type1}-{self.type2}-{self.type3}-{self.type4}"
        reverse = f"{self.type4}-{self.type3}-{self.type2}-{self.type1}"
        return min(forward, reverse)


@dataclass
class ImproperParameter:
    """
    Improper dihedral parameter (for planarity).

    In AMBER, improper dihedrals maintain planarity of sp2 centers.
    """
    type1: str  # Central atom
    type2: str
    type3: str
    type4: str
    barrier: float    # PK in kcal/mol
    phase: float      # PHASE in degrees
    periodicity: int  # PN (usually 2)
    source: str = ""


@dataclass
class NonbondedParameter:
    """
    Van der Waals (Lennard-Jones) parameter.

    AMBER uses the form: E = eps * [(R*/r)^12 - 2*(R*/r)^6]
    where R* = radius (Rmin/2) and eps = well_depth.
    """
    atom_type: str
    radius: float       # Å (Rmin/2)
    well_depth: float   # kcal/mol (epsilon)
    comment: str = ""
    source: str = ""


# =============================================================================
# Convenience Type Aliases
# =============================================================================

# Type alias for atom type assignment dictionaries
# Key is coordinate tuple, value is AtomTypeAssignment
TypeAssignmentDict = Dict[Tuple[float, float, float], AtomTypeAssignment]


# =============================================================================
# Conversion Utilities
# =============================================================================

def convert_legacy_dict_to_assignment(
    coords: Tuple[float, float, float],
    legacy_dict: Dict[str, Any]
) -> AtomTypeAssignment:
    """
    Convert a legacy plain-dict type assignment to AtomTypeAssignment.

    This handles the old format used in _assign_atom_types() and
    _apply_systematic_renaming():

        {
            'chain': 'A',
            'resname': 'HIS',
            'renamed_type': 'Y1',
            'is_center': False,
            'is_metal_ligand': True,
            ...
        }

    Args:
        coords: Atom coordinates
        legacy_dict: Old-style dict assignment

    Returns:
        AtomTypeAssignment instance
    """
    # Handle terminal_type which might be string or enum
    terminal_type = legacy_dict.get('terminal_type', 'middle')
    if isinstance(terminal_type, str):
        terminal_map = {
            'N_TERMINAL': TerminalType.N_TERMINAL,
            'C_TERMINAL': TerminalType.C_TERMINAL,
            'INTERNAL': TerminalType.INTERNAL,
            'middle': TerminalType.INTERNAL,
            'n_terminal': TerminalType.N_TERMINAL,
            'c_terminal': TerminalType.C_TERMINAL,
        }
        terminal_type = terminal_map.get(terminal_type, TerminalType.INTERNAL)

    # Handle source
    source_str = legacy_dict.get('source', 'force_field')
    try:
        source = AtomSource(source_str)
    except ValueError:
        if 'metal' in source_str.lower():
            source = AtomSource.METAL_DATABASE
        elif 'antechamber' in source_str.lower():
            source = AtomSource.ANTECHAMBER
        elif 'manual' in source_str.lower():
            source = AtomSource.MANUAL
        else:
            source = AtomSource.FORCE_FIELD

    return AtomTypeAssignment(
        coords=coords,
        chain=legacy_dict.get('chain', ''),
        resname=legacy_dict.get('resname', ''),
        ff_resname=legacy_dict.get('ff_resname', ''),
        resid=legacy_dict.get('resid', 0),
        atom_name=legacy_dict.get('atom_name', ''),
        element=legacy_dict.get('element', ''),
        original_type=legacy_dict.get('original_type', ''),
        renamed_type=legacy_dict.get('renamed_type', ''),
        is_center=legacy_dict.get('is_center', False),
        is_metal_ligand=legacy_dict.get('is_metal_ligand', False),
        terminal_type=terminal_type,
        charge=legacy_dict.get('charge'),
        mass=legacy_dict.get('mass'),
        spin=legacy_dict.get('spin'),
        vdw_radius=legacy_dict.get('vdw_params', {}).get('radius') if isinstance(legacy_dict.get('vdw_params'), dict) else None,
        vdw_epsilon=legacy_dict.get('vdw_params', {}).get('epsilon') if isinstance(legacy_dict.get('vdw_params'), dict) else None,
        vdw_source=legacy_dict.get('vdw_params', {}).get('source', '') if isinstance(legacy_dict.get('vdw_params'), dict) else '',
        source=source,
        source_detail=legacy_dict.get('source', ''),
        library_source=legacy_dict.get('library_source', ''),
        renamed=legacy_dict.get('renamed', False)
    )


def convert_assignment_to_legacy_dict(assignment: AtomTypeAssignment) -> Dict[str, Any]:
    """
    Convert AtomTypeAssignment to legacy plain-dict format.

    For backwards compatibility with code that expects the old format.

    Args:
        assignment: AtomTypeAssignment instance

    Returns:
        Dict in legacy format
    """
    result = {
        'chain': assignment.chain,
        'resname': assignment.resname,
        'ff_resname': assignment.ff_resname,
        'resid': assignment.resid,
        'atom_name': assignment.atom_name,
        'element': assignment.element,
        'original_type': assignment.original_type,
        'renamed_type': assignment.renamed_type,
        'is_center': assignment.is_center,
        'is_metal_ligand': assignment.is_metal_ligand,
        'terminal_type': assignment.terminal_type.value,
        'charge': assignment.charge,
        'renamed': assignment.renamed,
        'source': assignment.source_detail or assignment.source.value,
        'library_source': assignment.library_source,
    }

    # Add VDW params if present
    if assignment.vdw_radius is not None:
        result['vdw_params'] = {
            'radius': assignment.vdw_radius,
            'epsilon': assignment.vdw_epsilon,
            'source': assignment.vdw_source
        }

    # Add metal-specific params
    if assignment.spin is not None:
        result['spin'] = assignment.spin
    if assignment.mass is not None:
        result['mass'] = assignment.mass

    return result
