"""
ProPrep Forcefield Parameters Module

This module provides centralized management of forcefield parameters for redox centers,
supporting multiple oxidation states and spin states for accurate molecular dynamics
simulations.

Directory Structure:
    specialized_residues/
    ├── cofactor_family/
    │   └── specific_type/
    │       ├── metadata.json
    │       └── redox_state/
    │           └── spin_state/
    │               ├── files.frcmod
    │               └── files.lib

Example:
    specialized_residues/
    ├── heme/
    │   └── bis_his_c_type/
    │       ├── metadata.json
    │       ├── fe2/
    │       │   ├── high_spin/
    │       │   │   ├── HEH_hs.frcmod
    │       │   │   └── HEH_hs.lib
    │       │   └── low_spin/
    │       │       ├── HEL_ls.frcmod
    │       │       └── HEL_ls.lib
    │       └── fe3/
    │           ├── high_spin/
    │           └── low_spin/
"""

from .loader import (
    load_forcefield_metadata,
    discover_forcefield_files,
    validate_forcefield_structure,
    get_available_cofactor_types,
    get_atom_types_for_state,
    get_residue_name_for_state,
    get_protonation_model,
    resolve_residue_names,
    get_prerequisite_leaprcs,
    get_prerequisite_leaprc_groups,
    find_companion_set,
    get_forcefield_base_path,
    get_user_forcefield_base_path,
    residue_has_template_connectivity,
    reset_template_connectivity_cache,
    ForcefieldNotFoundError,
    InvalidForcefieldError
)
from .user_library import (
    promote_state,
    PromotionRequest,
    UserLibraryError,
    LibraryCollisionError,
)

__all__ = [
    'load_forcefield_metadata',
    'discover_forcefield_files',
    'validate_forcefield_structure',
    'get_available_cofactor_types',
    'get_atom_types_for_state',
    'get_residue_name_for_state',
    'get_protonation_model',
    'resolve_residue_names',
    'get_prerequisite_leaprcs',
    'get_prerequisite_leaprc_groups',
    'find_companion_set',
    'get_forcefield_base_path',
    'get_user_forcefield_base_path',
    'residue_has_template_connectivity',
    'reset_template_connectivity_cache',
    'ForcefieldNotFoundError',
    'InvalidForcefieldError',
    'promote_state',
    'PromotionRequest',
    'UserLibraryError',
    'LibraryCollisionError',
]