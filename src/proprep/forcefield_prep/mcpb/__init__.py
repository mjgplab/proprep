"""
MCPB Interface Module

Provides atom typing and fingerprint generation for MCPB (Metal Center Parameter Builder)
integration with ProPrep RedoxSite objects.

Main components:
- GlobalAtomTypeRegistry: Coordinates atom types across multiple metal sites
- AtomTyper: Assigns AMBER atom types to RedoxSite atoms
- FingerprintGenerator: Generates MCPB-compatible fingerprint files
- ForceFieldLibrary: Manages force field atom type libraries
- TerminalClassifier: Detects protein terminal residues
- QMInterface: Manages QM calculation workflow
- HessianParser: Extracts Hessian matrices from QM outputs
- SeminarioMethod: Computes force constants from Hessian
- FrcmodBuilder: Generates AMBER frcmod files
"""

from .global_atom_registry import GlobalAtomTypeRegistry, GlobalTypeAssignment
from .atom_typer import MCPBAtomTyper, AtomTypeRenamingMode, AtomTypeAssignment
from .fingerprint_generator import FingerprintGenerator
from .terminal_classifier import TerminalClassifier, TerminalType, TerminalDetectionResult
from .ff_library import ForceFieldLibrary
from .qm_interface import QMInterface, QMSoftware, QMCalculationMode
from .hessian_parser import HessianParser, HessianParseError
from .seminario import SeminarioMethod, BondParameter, AngleParameter
from .frcmod_builder import FrcmodBuilder, MassEntry, NonbondedEntry

__all__ = [
    # Multi-site coordination (NEW)
    'GlobalAtomTypeRegistry',
    'GlobalTypeAssignment',
    # Step 1 - Atom typing and fingerprints
    'MCPBAtomTyper',
    'AtomTypeRenamingMode',
    'AtomTypeAssignment',
    'FingerprintGenerator',
    'TerminalClassifier',
    'TerminalType',
    'TerminalDetectionResult',
    'ForceFieldLibrary',
    # Step 2 - QM and bonded parameters
    'QMInterface',
    'QMSoftware',
    'QMCalculationMode',
    'HessianParser',
    'HessianParseError',
    'SeminarioMethod',
    'BondParameter',
    'AngleParameter',
    'FrcmodBuilder',
    'MassEntry',
    'NonbondedEntry',
]
