"""
ONIOM QM/MM Preparation Module

Provides tools for setting up Gaussian ONIOM calculations from RedoxSite objects.
Integrates with ProPrep's workflow for automated QM/MM input generation.

© 2024 ProPrep Developer. All rights reserved.
"""

import logging
logging.getLogger(__name__).setLevel(logging.WARNING)

from .data_structures import (
    ONIOMLayer,
    FreezeFlag,
    LayerAssignment,
    LinkAtom,
    ConnectivityEntry,
    ONIOMSetup,
    LayerStatistics,
)

# Import the preparator module (registered via the unified QM/MM Preparator)
from . import oniom_qmmm_preparator

__all__ = [
    "ONIOMLayer",
    "FreezeFlag",
    "LayerAssignment",
    "LinkAtom",
    "ConnectivityEntry",
    "ONIOMSetup",
    "LayerStatistics",
]

__version__ = "0.1.0"
