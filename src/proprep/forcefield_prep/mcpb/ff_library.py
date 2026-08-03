"""
Force Field Library Manager

Manages AMBER force field atom type libraries.
Loads from AMBER mol2 files and custom mol2 files.
"""

import os
import logging
from typing import Dict, Optional, Tuple


class ForceFieldLibrary:
    """
    Manages force field atom type libraries.

    Loads atom types and charges from:
    1. AMBER force field mol2 files (ff19SB, ff14SB, etc.)
    2. Custom mol2 files (for non-standard residues)

    Library format: "ResName-AtomName" -> {atom_type, charge, source}
    """

    # Force field to mol2 file mapping
    FF_MOL2_MAP = {
        'ff19SB': 'parm19.mol2',
        'ff14SB': 'parm12.mol2',
        'ff14SB.redq': 'parm12_redq.mol2',
        'ff14SBonlysc': 'parm12.mol2',
        'ff10': 'parm10.mol2',
        'ff03': 'parm03.mol2',
        'ff03.r1': 'parm03_r1.mol2',
        'ff99SB': 'parm94.mol2',
        'ff99': 'parm94.mol2',
        'ff94': 'parm94.mol2',
        'ff15ipq': 'parm15ipq_10.0.mol2',
        'ff15ipq-vac': 'parm15ipq-vac_10.0.mol2',
        'fb15': 'parm_fb15.mol2',
    }

    def __init__(self, force_field: str, amberhome: Optional[str] = None, logger=None):
        """
        Initialize force field library.

        Args:
            force_field: Force field name (e.g., 'ff19SB')
            amberhome: Path to AMBER installation (auto-detect if None)
            logger: Optional logger instance
        """
        self.force_field = force_field
        self.logger = logger or logging.getLogger(__name__)

        # Detect AMBERHOME
        if amberhome:
            self.amberhome = amberhome
        else:
            self.amberhome = self._detect_amberhome()

        # Library dictionary: "ResName-AtomName" -> {atom_type, charge, source}
        self.library: Dict[str, Dict] = {}

        # Load force field library
        self._load_force_field_library()

    def _detect_amberhome(self) -> str:
        """Auto-detect AMBER installation from environment"""
        amberhome = os.getenv('AMBERHOME')
        if not amberhome:
            raise ValueError(
                "AMBERHOME environment variable not set. "
                "Please set AMBERHOME or provide amberhome parameter."
            )
        if not os.path.exists(amberhome):
            raise ValueError(f"AMBERHOME directory does not exist: {amberhome}")

        self.logger.info(f"Using AMBERHOME: {amberhome}")
        return amberhome

    def _load_force_field_library(self):
        """Load force field mol2 file into library dictionary"""
        mol2_path = self._get_ff_mol2_path()

        if not os.path.exists(mol2_path):
            raise FileNotFoundError(
                f"Force field mol2 file not found: {mol2_path}\n"
                f"Force field: {self.force_field}\n"
                f"Expected path: {mol2_path}"
            )

        self.logger.info(f"Loading force field library: {mol2_path}")
        self._parse_mol2_to_library(mol2_path, source=self.force_field)
        self.logger.info(f"Loaded {len(self.library)} atom type definitions")

    def _get_ff_mol2_path(self) -> str:
        """
        Map force field name to mol2 file path.

        Returns:
            Full path to mol2 file
        """
        mol2_file = self.FF_MOL2_MAP.get(self.force_field)
        if not mol2_file:
            raise ValueError(
                f"Unknown force field: {self.force_field}\n"
                f"Supported force fields: {list(self.FF_MOL2_MAP.keys())}"
            )

        return os.path.join(self.amberhome, 'dat', 'pymsmt', mol2_file)

    def _parse_mol2_to_library(self, mol2_path: str, source: str):
        """
        Parse mol2 file and extract atom type information.

        Mol2 format:
            @<TRIPOS>ATOM
            AtomID AtomName X Y Z AtomType ResID ResName Charge

        Args:
            mol2_path: Path to mol2 file
            source: Source identifier (e.g., 'ff19SB', 'custom:ZN.mol2')
        """
        try:
            # Use pymsmt's mol2io if available
            from pymsmt.mol.mol2io import get_atominfo

            mol, atids, resids = get_atominfo(mol2_path)

            for resid in resids:
                resname = mol.residues[resid].resname
                for atomid in mol.residues[resid].resconter:
                    atom = mol.atoms[atomid]
                    key = f"{resname}-{atom.atname}"

                    # Store atom type info
                    self.library[key] = {
                        'atom_type': atom.atomtype.strip(),
                        'charge': atom.charge,
                        'source': source,
                        'resname': resname,
                        'atom_name': atom.atname
                    }

        except ImportError:
            # Fallback: manual parsing if pymsmt not available
            self.logger.warning("pymsmt not available, using manual mol2 parsing")
            self._parse_mol2_manual(mol2_path, source)

    def _parse_mol2_manual(self, mol2_path: str, source: str):
        """
        Manual mol2 parsing (fallback if pymsmt not available).

        Args:
            mol2_path: Path to mol2 file
            source: Source identifier
        """
        in_atom_section = False
        current_resname = None

        with open(mol2_path, 'r') as f:
            for line in f:
                line = line.strip()

                if line == '@<TRIPOS>ATOM':
                    in_atom_section = True
                    continue
                elif line.startswith('@<TRIPOS>'):
                    in_atom_section = False
                    continue

                if in_atom_section and line:
                    parts = line.split()
                    if len(parts) >= 9:
                        # AtomID AtomName X Y Z AtomType ResID ResName Charge
                        atom_name = parts[1]
                        atom_type = parts[5]
                        resname = parts[7]
                        try:
                            charge = float(parts[8])
                        except (ValueError, IndexError):
                            charge = 0.0

                        key = f"{resname}-{atom_name}"

                        self.library[key] = {
                            'atom_type': atom_type.strip(),
                            'charge': charge,
                            'source': source,
                            'resname': resname,
                            'atom_name': atom_name
                        }

    def add_custom_library(self, mol2_path: str):
        """
        Add custom mol2 file (overwrites existing entries).

        Args:
            mol2_path: Path to custom mol2 file

        Raises:
            FileNotFoundError: If mol2 file doesn't exist
        """
        if not os.path.exists(mol2_path):
            raise FileNotFoundError(f"Custom mol2 file not found: {mol2_path}")

        source = f"custom:{os.path.basename(mol2_path)}"
        self.logger.info(f"Loading custom library: {mol2_path}")

        original_count = len(self.library)
        self._parse_mol2_to_library(mol2_path, source)
        new_count = len(self.library)

        self.logger.info(
            f"Added {new_count - original_count} new entries "
            f"from custom library (total: {new_count})"
        )

    def get_atom_type(self, lookup_key: str) -> Optional[Dict]:
        """
        Get atom type info for a lookup key.

        Args:
            lookup_key: Format: "ResName-AtomName", "NResName-AtomName", or "CResName-AtomName"
                       Examples: "CYS-SG", "NCYS-N", "CCYS-C"

        Returns:
            Dictionary with atom_type, charge, source, resname, atom_name
            or None if not found
        """
        return self.library.get(lookup_key)

    def get_charge(self, lookup_key: str) -> Optional[float]:
        """
        Get atom charge for a lookup key.

        Args:
            lookup_key: Format: "ResName-AtomName"

        Returns:
            Charge value or None if not found
        """
        info = self.get_atom_type(lookup_key)
        return info['charge'] if info else None

    def has_residue(self, resname: str) -> bool:
        """
        Check if residue type exists in library.

        Args:
            resname: Residue name (e.g., 'CYS', 'HIS')

        Returns:
            True if residue has at least one atom defined
        """
        for key in self.library.keys():
            if key.startswith(f"{resname}-"):
                return True
        return False

    def get_residue_atoms(self, resname: str) -> Dict[str, Dict]:
        """
        Get all atoms for a specific residue type.

        Args:
            resname: Residue name (e.g., 'CYS', 'NCYS')

        Returns:
            Dictionary mapping atom names to atom type info
        """
        residue_atoms = {}
        prefix = f"{resname}-"

        for key, info in self.library.items():
            if key.startswith(prefix):
                atom_name = key[len(prefix):]
                residue_atoms[atom_name] = info

        return residue_atoms

    def get_library_statistics(self) -> Dict:
        """
        Get statistics about the loaded library.

        Returns:
            Dictionary with library statistics
        """
        # Count residues
        residue_names = set()
        for key in self.library.keys():
            resname = key.split('-')[0]
            residue_names.add(resname)

        # Count sources
        sources = set()
        for info in self.library.values():
            sources.add(info['source'])

        return {
            'total_entries': len(self.library),
            'unique_residues': len(residue_names),
            'sources': list(sources),
            'force_field': self.force_field
        }

    def __repr__(self):
        stats = self.get_library_statistics()
        return (
            f"ForceFieldLibrary("
            f"ff={self.force_field}, "
            f"entries={stats['total_entries']}, "
            f"residues={stats['unique_residues']})"
        )
