"""
ORCA QM/MM Preparation Module

Provides tools for setting up ORCA QM/MM calculations from RedoxSite objects
and AMBER topology files. Integrates with ProPrep's workflow for automated
QM/MM input generation.
"""

import logging
logging.getLogger(__name__).setLevel(logging.WARNING)

from .orca_data_structures import (
    ORCAQMMMSetup,
    compress_indices,
    map_redox_atoms_to_prmtop,
)

# Import the preparator module (registered via the unified QM/MM Preparator)
from . import orca_qmmm_preparator

__all__ = [
    "ORCAQMMMSetup",
    "compress_indices",
    "map_redox_atoms_to_prmtop",
]

__version__ = "1.0.0"
