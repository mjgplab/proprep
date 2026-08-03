#!/usr/bin/env python3
"""
Data Models for Transformer Creator v2

Clean, declarative data structures using Pydantic to represent
all components of a RedoxSite transformer specification.
"""

from typing import List, Dict, Any, Optional, Union, Literal
from pydantic import BaseModel, Field, field_validator, model_validator
from enum import Enum


# ===== ENUMS =====

class CenterTypeEnum(str, Enum):
    """Redox center types"""
    METAL_ION = "metal_ion"
    ORGANOMETALLIC_COFACTOR = "organometallic_cofactor"  # Cofactors with embedded metal (heme, Fe-S)
    ORGANIC_COFACTOR = "organic_cofactor"                 # Purely organic cofactors
    REDOX_AMINO_ACID = "redox_amino_acid"


class StageType(str, Enum):
    """Component matching stage types"""
    INTRINSIC_FILTER = "intrinsic_filter"
    POSITIONAL = "positional"
    ARITHMETIC = "arithmetic"
    STRUCTURAL_INFO = "structural_info"
    EXCLUDE_MATCHED = "exclude_matched"
    CUSTOM = "custom"


class ParameterType(str, Enum):
    """Parameter types"""
    CHOICE = "choice"
    FIXED = "fixed"


# ===== SITE REQUIREMENTS =====

class CenterRequirements(BaseModel):
    """Requirements for redox centers"""
    required_count: int = Field(ge=1, description="Number of centers required")
    center_types: List[CenterTypeEnum] = Field(description="Allowed center types")
    elements: Optional[List[str]] = Field(default=None, description="Required elements (e.g., ['FE', 'CU'])")


class ResidueRequirement(BaseModel):
    """Requirements for a specific residue type"""
    min_count: int = Field(ge=0, description="Minimum count")
    max_count: Optional[int] = Field(default=None, description="Maximum count (None = unlimited)")


class AtomRequirements(BaseModel):
    """Requirements for atoms/residues in site"""
    required_residues: Dict[str, ResidueRequirement] = Field(
        default_factory=dict,
        description="Residue name -> count requirements"
    )
    alternative_groups: Optional[List[List[str]]] = Field(
        default=None,
        description="Groups of residues where only one from each group is needed"
    )


class BondGroup(BaseModel):
    """A group of required bonds"""
    description: str = Field(description="Description of this bond group")
    min_count: int = Field(ge=1, description="Minimum bonds that must match")
    bond_types: Dict[str, List[tuple]] = Field(
        description="Map of chemical_type -> list of (atom1_spec, atom2_spec) tuples"
    )


class BondRequirements(BaseModel):
    """Requirements for bonds in site"""
    required_bond_groups: List[BondGroup] = Field(
        default_factory=list,
        description="Groups of bonds (at least one group must be satisfied)"
    )
    require_one_group: bool = Field(
        default=False,
        description="If True, at least one bond group must meet min_count"
    )


class SiteRequirements(BaseModel):
    """Complete site matching requirements"""
    centers: CenterRequirements
    atoms: Optional[AtomRequirements] = None
    bonds: Optional[BondRequirements] = None


# ===== COMPONENT MATCHING =====

class IntrinsicFilterStage(BaseModel):
    """Filter by intrinsic properties"""
    type: Literal["intrinsic_filter"] = "intrinsic_filter"
    resname: Optional[str] = None
    resnames: Optional[List[str]] = None  # Alternative: list of acceptable names
    coord_atom: Optional[str] = None
    coord_atoms: Optional[List[str]] = None  # For multi-dentate ligands
    min_distance: Optional[float] = None
    max_distance: Optional[float] = None


class PositionalStage(BaseModel):
    """Sort and select by position"""
    type: Literal["positional"] = "positional"
    sort_by: str = Field(description="'resid', 'distance', or 'chain'")
    position: int = Field(description="Position in sorted list (0=first, -1=last)")
    reverse: bool = Field(default=False, description="Reverse sort order")
    reference_component: Optional[str] = Field(
        default=None,
        description="Component to measure distance from (if sort_by='distance')"
    )
    reference_atom: Optional[str] = Field(
        default=None,
        description="Atom in reference component (if applicable)"
    )


class ArithmeticStage(BaseModel):
    """Calculate from another component"""
    type: Literal["arithmetic"] = "arithmetic"
    base_component: str = Field(description="Component to calculate from")
    operation: str = Field(description="'add' or 'subtract'")
    value: int = Field(description="Value to add/subtract")


class StructuralInfoStage(BaseModel):
    """Use structural annotations"""
    type: Literal["structural_info"] = "structural_info"
    info_type: str = Field(description="Type of structural info to use")
    # Can be extended with specific fields as needed


class ExcludeMatchedStage(BaseModel):
    """Exclude already-matched components"""
    type: Literal["exclude_matched"] = "exclude_matched"


class CustomStage(BaseModel):
    """Custom Python logic"""
    type: Literal["custom"] = "custom"
    code: str = Field(description="Python code to execute")
    description: str = Field(description="Human-readable description")


# Union of all stage types
MatchingStage = Union[
    IntrinsicFilterStage,
    PositionalStage,
    ArithmeticStage,
    StructuralInfoStage,
    ExcludeMatchedStage,
    CustomStage
]


class ComponentMatchingRule(BaseModel):
    """Rule for matching a single component"""
    description: str = Field(description="Human-readable description")
    stages: List[MatchingStage] = Field(
        default_factory=list,
        description="Sequential stages to identify component"
    )


class ComponentMatching(BaseModel):
    """Complete component matching specification"""
    rules: Dict[str, ComponentMatchingRule] = Field(
        default_factory=dict,
        description="Component name -> matching rule"
    )


# ===== TRANSFORMATIONS =====

class TransformationSelector(BaseModel):
    """Selector for transformation target"""
    chain_id: Optional[str] = None  # Can be a component field like "center_chain"
    residue_id: Optional[Union[int, str]] = None  # Can be component field like "center_id"
    residue_name: Optional[str] = None
    atom_names: Optional[List[str]] = None
    insertion_code: Optional[str] = None


class TransformationAction(BaseModel):
    """Actions to apply in transformation"""
    change_residue_name: Optional[str] = None
    change_residue_id: Optional[Union[int, str]] = None
    change_chain_id: Optional[str] = None
    change_insertion_code: Optional[str] = None
    rename_atoms: Optional[Dict[str, str]] = None  # old_name -> new_name
    convert_to_hetatm: Optional[bool] = None


class Transformation(BaseModel):
    """Single transformation step"""
    id: str = Field(description="Unique transformation ID")
    description: str = Field(description="Human-readable description")
    selector: TransformationSelector
    action: TransformationAction


class TransformationSequence(BaseModel):
    """Sequence of transformations"""
    transformations: List[Transformation] = Field(default_factory=list)


# ===== PARAMETERS =====

class ParameterDefinition(BaseModel):
    """Definition of a user-configurable parameter"""
    description: str
    type: ParameterType
    options: Optional[List[str]] = None  # For CHOICE type
    default: Optional[str] = None  # For CHOICE type
    value: Optional[str] = None  # For FIXED type
    note: Optional[str] = None


class Parameters(BaseModel):
    """Parameter definitions and mappings"""
    definitions: Dict[str, ParameterDefinition] = Field(
        default_factory=dict,
        description="Parameter name -> definition"
    )
    mappings: Dict[str, Dict[str, str]] = Field(
        default_factory=dict,
        description="Parameter combination -> residue name mappings"
    )


# ===== COMPLETE TRANSFORMER SPECIFICATION =====

class TransformerSpec(BaseModel):
    """Complete specification for a RedoxSite transformer"""

    # Basic metadata
    name: str = Field(description="Transformer name (snake_case)")
    class_name: Optional[str] = Field(default=None, description="Python class name (CamelCase)")
    description: str = Field(description="Human-readable description")
    supported_site_types: Optional[List[str]] = Field(default=None, description="Site types this transformer handles")

    # Core specifications
    site_requirements: SiteRequirements
    component_matching: ComponentMatching
    transformation_sequence: TransformationSequence

    # Optional
    parameters: Optional[Parameters] = None
    forcefield_path: Optional[str] = None

    # Metadata
    template_source: Optional[str] = Field(
        default=None,
        description="Name of transformer used as template (if any)"
    )
    required_residue_count: int = Field(
        default=1,
        description="Number of residue IDs needed (for ID allocation)"
    )

    @model_validator(mode='after')
    def auto_generate_fields(self):
        """Auto-generate class_name and supported_site_types from name if not provided"""
        # Generate class_name if not provided
        if not self.class_name and self.name:
            # Convert snake_case to CamelCase
            class_name = ''.join(word.capitalize() for word in self.name.split('_'))
            self.class_name = f"{class_name}Transformer"

        # Generate supported_site_types if not provided
        if not self.supported_site_types and self.name:
            self.supported_site_types = [self.name]

        return self


# ===== HELPER FUNCTIONS =====

def create_empty_spec(name: str, description: str) -> TransformerSpec:
    """Create an empty transformer spec with minimal required fields"""
    return TransformerSpec(
        name=name,
        class_name='',  # Will be auto-generated
        description=description,
        supported_site_types=[],  # Will be auto-generated
        site_requirements=SiteRequirements(
            centers=CenterRequirements(
                required_count=1,
                center_types=[CenterTypeEnum.METAL_ION]
            )
        ),
        component_matching=ComponentMatching(),
        transformation_sequence=TransformationSequence()
    )
