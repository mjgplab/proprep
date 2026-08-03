"""
Force Field Loader for Solvation Estimation

Loads R* (van der Waals radius) values from AMBER force field files for accurate
solvation box estimation before running tLEaP.

Uses the existing ff_parsers module to read leaprc, parm*.dat, frcmod, and lib files.
"""

import logging
import os
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from .ff_parsers import (
    parse_leaprc,
    parse_parm_dat,
    parse_frcmod,
    parse_lib_file,
    resolve_leaprc_path,
    find_amber_file,
)
from .ff_types import NonbondedParameter

logger = logging.getLogger(__name__)


class ForceFieldLoader:
    """
    Load force field parameters for accurate solvation estimation.

    Uses existing ff_parsers to read R* values from parm*.dat and frcmod files,
    and atom type mappings from lib files.

    Supports both standard AMBER forcefields (loaded via leaprc) and custom
    frcmod files from redox site parameterization.
    """

    # Default R* values from tLEaP source (atom.h)
    ATOM_DEFAULT_RADIUS = 1.5  # For unknown heavy atom types
    HYDROGEN_DEFAULT_RADIUS = 1.0  # For unknown hydrogen types

    def __init__(self, amberhome: Optional[str] = None):
        """
        Initialize the force field loader.

        Args:
            amberhome: Path to AMBER installation. If None, uses $AMBERHOME.
        """
        if amberhome is None:
            amberhome = os.environ.get('AMBERHOME', '')

        self.amberhome = Path(amberhome) if amberhome else None

        if self.amberhome and not self.amberhome.exists():
            logger.warning(f"AMBERHOME directory does not exist: {self.amberhome}")
            self.amberhome = None

        # Parameter storage
        self.nonbonded_params: Dict[str, NonbondedParameter] = {}
        self.atom_type_map: Dict[Tuple[str, str], str] = {}  # (resname, atomname) -> atom_type

        # Track loaded files to avoid duplicates
        self._loaded_files: Set[str] = set()

        # Track warnings for missing types
        self._missing_types: Set[str] = set()
        self._unmapped_atoms: Set[Tuple[str, str]] = set()

    def load_from_leaprc(self, leaprc_name: str) -> bool:
        """
        Load parameters from a leaprc file, following source directives recursively.

        Args:
            leaprc_name: Name like 'leaprc.protein.ff19SB' or full path

        Returns:
            True if loaded successfully, False otherwise
        """
        if not self.amberhome:
            logger.warning("AMBERHOME not set, cannot load leaprc files")
            return False

        leaprc_path = resolve_leaprc_path(leaprc_name, self.amberhome)
        if not leaprc_path:
            logger.warning(f"Could not find leaprc file: {leaprc_name}")
            return False

        self._load_leaprc_recursive(leaprc_path)
        return True

    def _load_leaprc_recursive(self, leaprc_path: Path) -> None:
        """Recursively load leaprc and all sourced files."""
        path_str = str(leaprc_path)
        if path_str in self._loaded_files:
            return
        self._loaded_files.add(path_str)

        logger.debug(f"Loading leaprc: {leaprc_path}")
        contents = parse_leaprc(leaprc_path)

        # Follow source directives first (base parameters)
        for source_file in contents.source_files:
            source_path = resolve_leaprc_path(source_file, self.amberhome)
            if source_path:
                self._load_leaprc_recursive(source_path)

        # Load parm.dat files (base parameters)
        for parm_file in contents.parm_files:
            self._load_parm_file(parm_file)

        # Load frcmod files (these override parm.dat)
        for frcmod_file in contents.frcmod_files:
            self._load_frcmod_file(frcmod_file)

        # Load lib files for atom type mappings
        for lib_file in contents.lib_files:
            self._load_lib_file(lib_file)

    def _load_parm_file(self, filename: str) -> None:
        """Load parameters from a parm*.dat file."""
        if not self.amberhome:
            return

        path = find_amber_file(filename, self.amberhome, 'parm')
        if not path:
            logger.debug(f"Could not find parm file: {filename}")
            return

        path_str = str(path)
        if path_str in self._loaded_files:
            return
        self._loaded_files.add(path_str)

        logger.debug(f"Loading parm file: {path}")
        contents = parse_parm_dat(path)

        # Base parameters - don't overwrite existing (frcmod takes precedence)
        for atom_type, param in contents.nonbonded_parameters.items():
            if atom_type not in self.nonbonded_params:
                self.nonbonded_params[atom_type] = param

    def _load_frcmod_file(self, filename: str) -> None:
        """Load parameters from a frcmod file (overrides base params)."""
        if not self.amberhome:
            return

        path = find_amber_file(filename, self.amberhome, 'frcmod')
        if not path:
            logger.debug(f"Could not find frcmod file: {filename}")
            return

        path_str = str(path)
        if path_str in self._loaded_files:
            return
        self._loaded_files.add(path_str)

        logger.debug(f"Loading frcmod file: {path}")
        contents = parse_frcmod(path)

        # frcmod overrides existing parameters
        self.nonbonded_params.update(contents.nonbonded_parameters)

    def _load_lib_file(self, filename: str) -> None:
        """Load atom type mappings from a lib file."""
        if not self.amberhome:
            return

        path = find_amber_file(filename, self.amberhome, 'lib')
        if not path:
            logger.debug(f"Could not find lib file: {filename}")
            return

        path_str = str(path)
        if path_str in self._loaded_files:
            return
        self._loaded_files.add(path_str)

        logger.debug(f"Loading lib file: {path}")
        atoms = parse_lib_file(path)

        for atom in atoms:
            key = (atom.residue_name, atom.atom_name)
            self.atom_type_map[key] = atom.atom_type

    def load_custom_frcmod(self, frcmod_path: str) -> bool:
        """
        Load custom frcmod file (e.g., from redox site parameterization).

        These override any existing parameters (custom takes precedence).

        Args:
            frcmod_path: Full path to custom frcmod file

        Returns:
            True if loaded successfully, False otherwise
        """
        path = Path(frcmod_path)
        if not path.exists():
            logger.warning(f"Custom frcmod file not found: {frcmod_path}")
            return False

        path_str = str(path)
        if path_str in self._loaded_files:
            return True  # Already loaded
        self._loaded_files.add(path_str)

        logger.debug(f"Loading custom frcmod: {path}")
        contents = parse_frcmod(path)
        self.nonbonded_params.update(contents.nonbonded_parameters)
        return True

    def load_custom_lib(self, lib_path: str) -> bool:
        """
        Load custom lib file for atom type mappings.

        Args:
            lib_path: Full path to custom lib file

        Returns:
            True if loaded successfully, False otherwise
        """
        path = Path(lib_path)
        if not path.exists():
            logger.warning(f"Custom lib file not found: {lib_path}")
            return False

        path_str = str(path)
        if path_str in self._loaded_files:
            return True
        self._loaded_files.add(path_str)

        logger.debug(f"Loading custom lib: {path}")
        atoms = parse_lib_file(path)

        for atom in atoms:
            key = (atom.residue_name, atom.atom_name)
            self.atom_type_map[key] = atom.atom_type

        return True

    def load_from_selected_forcefields(self, selected_forcefields: dict) -> None:
        """
        Load parameters from selected standard forcefields.

        Args:
            selected_forcefields: Dict from forcefield selection UX, e.g.:
                {
                    'protein': {'leaprc': 'leaprc.protein.ff19SB', ...},
                    'water': {'leaprc': 'leaprc.water.opc', ...},
                    'modified_aa': [{'leaprc': '...', ...}, ...],  # multi-select is list
                    ...
                }
        """
        for category, ff_info in selected_forcefields.items():
            if ff_info is None:
                continue

            # Handle multi-select categories (list of dicts)
            if isinstance(ff_info, list):
                for item in ff_info:
                    if isinstance(item, dict):
                        self._load_single_forcefield(item)
            elif isinstance(ff_info, dict):
                self._load_single_forcefield(ff_info)

    def _load_single_forcefield(self, ff_info: dict) -> None:
        """Load a single forcefield from its info dict."""
        # Load leaprc if specified
        if 'leaprc' in ff_info:
            self.load_from_leaprc(ff_info['leaprc'])

        # Load any additional frcmod files
        for frcmod in ff_info.get('frcmods', []):
            if self.amberhome:
                path = find_amber_file(frcmod, self.amberhome, 'frcmod')
                if path:
                    self._load_frcmod_file(frcmod)

    def load_from_custom_forcefield_dict(self, custom_forcefields: dict) -> None:
        """
        Load parameters from custom forcefield selection (redox sites, etc.).

        Args:
            custom_forcefields: Dict keyed by (transformer_type, redox_state, spin_state),
                values contain {'frcmod': '/path/to/file.frcmod', 'lib': '/path/to/file.lib'}
        """
        for key, ff_info in custom_forcefields.items():
            if not isinstance(ff_info, dict):
                continue

            if 'frcmod' in ff_info:
                self.load_custom_frcmod(ff_info['frcmod'])

            if 'lib' in ff_info:
                self.load_custom_lib(ff_info['lib'])

    def get_r_star(self, atom_type: str) -> float:
        """
        Get R* (van der Waals radius) for an atom type.

        Args:
            atom_type: AMBER atom type (e.g., 'CT', 'N', 'H1')

        Returns:
            R* value in Angstroms
        """
        param = self.nonbonded_params.get(atom_type)
        if param:
            return param.radius

        # Track missing type for warnings
        if atom_type not in self._missing_types:
            self._missing_types.add(atom_type)
            logger.debug(f"No R* value for atom type '{atom_type}', using default")

        # Use element-based default
        if atom_type.startswith('H'):
            return self.HYDROGEN_DEFAULT_RADIUS
        return self.ATOM_DEFAULT_RADIUS

    # Atoms that only exist in terminal residues
    N_TERMINAL_ATOMS = {'H1', 'H2', 'H3'}  # N-terminal amino hydrogens
    C_TERMINAL_ATOMS = {'OXT'}  # C-terminal carboxyl oxygen

    # Protonation state variants (the only mappings we need)
    PROTONATION_VARIANTS = {
        'HIS': ['HID', 'HIE', 'HIP'],
        'CYS': ['CYX', 'CYM'],
        'ASP': ['ASH', 'AS4'],
        'GLU': ['GLH', 'GL4'],
        'LYS': ['LYN'],
        'TYR': ['TYM'],
    }

    # Water name variants
    WATER_VARIANTS = {
        'HOH': ['WAT', 'TIP3', 'TP3', 'TP4', 'OPC', 'SPC', 'SPCE'],
        'WAT': ['HOH', 'TIP3'],
    }

    def get_atom_type(self, resname: str, atom_name: str) -> Optional[str]:
        """
        Get AMBER atom type for a residue/atom combination.

        Lookup order:
        1. Direct (resname, atom_name) lookup
        2. Protonation state variants (HIS->HID/HIE/HIP, etc.)
        3. N-terminal variant (NXXX) for H1/H2/H3
        4. C-terminal variant (CXXX) for OXT
        5. Water variants

        Args:
            resname: Residue name from PDB
            atom_name: Atom name from PDB

        Returns:
            AMBER atom type, or None if not found
        """
        # 1. Try direct lookup
        key = (resname, atom_name)
        if key in self.atom_type_map:
            return self.atom_type_map[key]

        # 2. Try protonation state variants
        if resname in self.PROTONATION_VARIANTS:
            for variant in self.PROTONATION_VARIANTS[resname]:
                key = (variant, atom_name)
                if key in self.atom_type_map:
                    return self.atom_type_map[key]

        # 3. Try N-terminal variant for terminal-specific atoms
        if atom_name in self.N_TERMINAL_ATOMS:
            n_term_resname = f'N{resname}'
            key = (n_term_resname, atom_name)
            if key in self.atom_type_map:
                return self.atom_type_map[key]

        # 4. Try C-terminal variant for OXT
        if atom_name in self.C_TERMINAL_ATOMS:
            c_term_resname = f'C{resname}'
            key = (c_term_resname, atom_name)
            if key in self.atom_type_map:
                return self.atom_type_map[key]

        # 5. Try water variants
        if resname in self.WATER_VARIANTS:
            for variant in self.WATER_VARIANTS[resname]:
                key = (variant, atom_name)
                if key in self.atom_type_map:
                    return self.atom_type_map[key]

        # Track unmapped atom for warnings
        self._unmapped_atoms.add((resname, atom_name))
        return None

    def get_r_star_for_atom(self, resname: str, atom_name: str, element: str = '') -> float:
        """
        Get R* value for an atom, with fallback to element-based defaults.

        Args:
            resname: Residue name from PDB
            atom_name: Atom name from PDB
            element: Element symbol (used for fallback)

        Returns:
            R* value in Angstroms
        """
        atom_type = self.get_atom_type(resname, atom_name)

        if atom_type:
            return self.get_r_star(atom_type)

        # Fallback based on element
        if element == 'H' or atom_name.startswith('H'):
            return self.HYDROGEN_DEFAULT_RADIUS
        return self.ATOM_DEFAULT_RADIUS

    def get_warnings(self) -> List[str]:
        """
        Get warnings about missing atom types and unmapped atoms.

        Returns:
            List of warning messages
        """
        warnings = []

        if self._missing_types:
            warnings.append(
                f"Missing R* values for {len(self._missing_types)} atom type(s): "
                f"{', '.join(sorted(self._missing_types)[:10])}"
                f"{'...' if len(self._missing_types) > 10 else ''}"
            )

        if self._unmapped_atoms:
            sample = list(self._unmapped_atoms)[:5]
            sample_str = ', '.join(f"{r}:{a}" for r, a in sample)
            warnings.append(
                f"Could not determine atom type for {len(self._unmapped_atoms)} atom(s): "
                f"{sample_str}{'...' if len(self._unmapped_atoms) > 5 else ''}"
            )

        return warnings

    def get_statistics(self) -> dict:
        """Get statistics about loaded parameters."""
        return {
            'n_nonbonded_params': len(self.nonbonded_params),
            'n_atom_type_mappings': len(self.atom_type_map),
            'n_loaded_files': len(self._loaded_files),
            'n_missing_types': len(self._missing_types),
            'n_unmapped_atoms': len(self._unmapped_atoms),
        }

    def __repr__(self):
        stats = self.get_statistics()
        return (
            f"ForceFieldLoader("
            f"nonbonded={stats['n_nonbonded_params']}, "
            f"mappings={stats['n_atom_type_mappings']}, "
            f"files={stats['n_loaded_files']})"
        )
