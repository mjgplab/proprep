"""
MPSA Structure Preparation Module

This module contains tools for preparing protein structures for molecular dynamics simulations.
It includes utilities for loading, filtering, analyzing and modifying PDB structures.
"""

from .altloc_selector import AltLocSelector, AltLocSelectorModule
from .chem_comp_dict_fetcher import CCDParser, display_residue_info, get_residue_info
from .disulfide_bond_detector_worker import DisulfideBondDetector
from .disulfide_bond_detector import DisulfideBondModule
from .homology_searcher import BLASTIntegrationModule
from .pdb_filter import PDBFilterModule
from .pdb_filter_worker import ComponentClassifier, PDBFilterTool

# Import all the modules to make them available when importing the package
from .pdb_loader import PDBLoaderModule, PDBMetadataExtractor, download_pdb
from .protonation_state_analyzer import ProtonationStateAnalyzer, ProtonationStateModule
from .structure_alignment import StructureAlignmentModule
from .interactive_structure_viewer import InteractiveStructureViewer
from .biological_assembly import BiologicalAssemblyGenerator
from .structure_orientation import StructureOrientationModule

# Define what gets imported with "from proprep.structure_prep import *"
__all__ = [
    # Classes
    "PDBLoaderModule",
    "PDBMetadataExtractor",
    "PDBFilterModule",
    "PDBFilterTool",
    "ComponentClassifier",
    "AltLocSelector",
    "AltLocSelectorModule",
    "DisulfideBondDetector",
    "DisulfideBondModule",
    "ProtonationStateAnalyzer",
    "ProtonationStateModule",
    "BLASTIntegrationModule",
    "StructureAlignmentModule",
    "InteractiveStructureViewer",
    "BiologicalAssemblyGenerator",
    "StructureOrientationModule",
    "CCDParser",
    "download_pdb",
    "get_residue_info",
    "display_residue_info",
]

# Version info
__version__ = "0.1.0"
