"""
Data models for Transformer Creator v3.

Role-based transformer specification using relationship phrases
and demonstrated PDB edits.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple


# ===== RELATIONSHIP PHRASES =====

class PhraseCategory(str, Enum):
    CONNECTIVITY = "connectivity"
    SPATIAL = "spatial"
    SEQUENCE = "sequence"


class SpatialComparison(str, Enum):
    CLOSEST = "closest"
    FARTHEST = "farthest"
    WITHIN = "within"


class SequenceComparison(str, Enum):
    LOWER = "lower"
    HIGHER = "higher"


class ChainRelation(str, Enum):
    SAME = "same"
    DIFFERENT = "different"


@dataclass
class ConnectivityPhrase:
    """[role] bonds [through [atom]] to [[atom] on] [role] [with [bond_type]]"""
    category: str = "connectivity"
    subject_role: str = ""
    target_role: str = ""
    source_atom: Optional[str] = None   # atom on subject role
    target_atom: Optional[str] = None   # atom on target role
    bond_type: Optional[str] = None     # coordinate, covalent, thioether, disulfide
    raw_text: str = ""


@dataclass
class SpatialPhrase:
    """[[atom] of] [role] is [closest to / farthest from / within N Å of] [[atom] on] [role]"""
    category: str = "spatial"
    subject_role: str = ""
    target_role: str = ""
    source_atom: Optional[str] = None   # atom on subject role
    target_atom: Optional[str] = None   # atom on target role
    comparison: str = ""                # "closest", "farthest", "within"
    distance: Optional[float] = None    # for "within N Å"
    raw_text: str = ""


@dataclass
class SequencePhrase:
    """
    Three sub-forms:
    - [role] is +N/-N from [role]
    - [role] has a lower/higher resid than [role]
    - [role] is on the same/different chain as [role]
    """
    category: str = "sequence"
    subject_role: str = ""
    raw_text: str = ""
    # Arithmetic form
    offset: Optional[int] = None
    offset_target_role: Optional[str] = None
    # Comparison form
    resid_comparison: Optional[str] = None  # "lower" or "higher"
    comparison_target_role: Optional[str] = None
    # Chain form
    chain_relation: Optional[str] = None    # "same" or "different"
    chain_target_role: Optional[str] = None


# Union type for all phrases
RelationshipPhrase = ConnectivityPhrase | SpatialPhrase | SequencePhrase


# ===== ROLE DEFINITIONS =====

@dataclass
class RoleDefinition:
    """A labeled residue in the site with its relationship constraints."""
    label: str                          # user's label (snake_case)
    resname: str                        # primary resname in current structure
    alt_resnames: List[str] = field(default_factory=list)  # alternatives for other structures
    chain: str = ""                     # in current structure
    resid: int = 0                      # in current structure
    is_center: bool = False
    atom_names: List[str] = field(default_factory=list)  # atoms in current structure
    phrases: List[RelationshipPhrase] = field(default_factory=list)


# ===== EDITING COMMANDS =====

class CommandType(str, Enum):
    RENAME_RESIDUE = "rename_residue"
    RENAME_ATOMS = "rename_atoms"
    MOVE_ATOMS = "move_atoms"
    MOVE_ATOMS_NEW = "move_atoms_new"
    HETATM = "hetatm"
    ATOM = "atom"


@dataclass
class EditCommand:
    """A single editing command in a pass."""
    command_type: CommandType
    raw_text: str
    source_role: str
    # For rename_residue
    new_resname: Optional[str] = None
    # For rename_atoms
    old_atom_names: Optional[List[str]] = None
    new_atom_names: Optional[List[str]] = None
    # For move_atoms / move_atoms_new
    atom_names: Optional[List[str]] = None
    target_role: Optional[str] = None
    # For move_atoms_new only
    new_residue_name: Optional[str] = None
    new_role_label: Optional[str] = None


@dataclass
class EditingPass:
    """A coherent batch of editing commands."""
    description: str
    commands: List[EditCommand] = field(default_factory=list)


# ===== PARAMETERS =====

@dataclass
class ParameterOption:
    """A single parameter definition."""
    name: str
    description: str
    param_type: str             # "choice" or "fixed"
    options: List[str] = field(default_factory=list)  # for choice type
    default: Optional[str] = None                      # for choice type
    fixed_value: Optional[str] = None                  # for fixed type
    note: Optional[str] = None                         # for fixed type


@dataclass
class ParameterSpec:
    """Parameter definitions and name mappings."""
    parameters: List[ParameterOption] = field(default_factory=list)
    # Maps parameter combo string → {abstract_name: concrete_resname}
    # e.g., {"reduced_low_spin": {"heme_name": "HCR", "proximal_name": "PHR"}}
    name_mappings: Dict[str, Dict[str, str]] = field(default_factory=dict)
    # Which role labels have state-dependent names
    state_dependent_roles: List[str] = field(default_factory=list)


# ===== TOP-LEVEL SPEC =====

@dataclass
class TransformerSpecV3:
    """Complete specification for a transformer, v3 format."""
    # Metadata
    name: str = ""
    description: str = ""
    forcefield_path: Optional[str] = None

    # Roles and relationships (steps 2-3)
    roles: List[RoleDefinition] = field(default_factory=list)

    # Editing passes (step 4)
    editing_passes: List[EditingPass] = field(default_factory=list)

    # Parameters (step 5)
    parameters: Optional[ParameterSpec] = None

    # Derived from RedoxSite (step 1, automatic)
    site_requirements: Optional[Dict[str, Any]] = None

    # Derived from new residues created during editing
    required_residue_count: int = 1
    residue_space_plan: Dict[str, int] = field(default_factory=dict)

    def get_role(self, label: str) -> Optional[RoleDefinition]:
        """Look up a role by label."""
        for role in self.roles:
            if role.label == label:
                return role
        return None

    def get_all_role_labels(self) -> List[str]:
        """Get all defined role labels."""
        return [r.label for r in self.roles]
