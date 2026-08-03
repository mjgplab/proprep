"""
Forcefield Parameterization Module

This module provides functionality for handling forcefield parameterization
of non-standard residues for molecular dynamics simulations, including
modified amino acids, small molecules, and metal sites.
"""

# Existing imports - keep exactly as they were
from proprep.forcefield_prep.forcefield_parameterizer import (
    ForcefieldParameterizer, NonStandardResidue)
from proprep.forcefield_prep.modified_amino_acid_parameterizer import (
    check_and_resume_pending_workflows, check_workflow_can_resume,
    parameterize_modified_residues, resume_paused_workflow, run_workflow)
from proprep.forcefield_prep.small_molecule_parameterizer import SmallMoleculeParameterizer

# Metal site parameterizer (step methods called by structure_preprocessor)
try:
    from proprep.forcefield_prep import metal_site_parameterizer
    _metal_site_available = True
except ImportError:
    metal_site_parameterizer = None
    _metal_site_available = False

# MCPB interface - for atom typing and fingerprint generation
try:
    from proprep.forcefield_prep import mcpb
    _mcpb_available = True
except ImportError:
    mcpb = None
    _mcpb_available = False

# ForceFieldLoader - for loading R* values for solvation estimation
from proprep.forcefield_prep.forcefield_loader import ForceFieldLoader

# Updated __all__ list - preserve existing exports and add new ones
__all__ = [
    # Existing exports (unchanged)
    'ForcefieldParameterizer',
    'NonStandardResidue',
    'SmallMoleculeParameterizer',
    'run_workflow',
    'parameterize_modified_residues',
    'resume_paused_workflow',
    'check_workflow_can_resume',
    'check_and_resume_pending_workflows',
    # ForceFieldLoader for solvation estimation
    'ForceFieldLoader',
]

# Add metal site exports if available
if _metal_site_available:
    __all__.append('metal_site_parameterizer')

if _mcpb_available:
    __all__.append('mcpb')

# Utility function to check metal site availability
def is_metal_site_available():
    """
    Check if metal site parameterization is available.
    
    Returns:
        bool: True if metal site parameterization can be used
    """
    return _metal_site_available

def get_available_parameterizers():
    """
    Get list of available parameterization methods.
    
    Returns:
        List[str]: Names of available parameterizers
    """
    available = ["Modified Amino Acids", "Small Molecules"]
    
    if _metal_site_available:
        available.append("Metal Sites")
    
    return available
