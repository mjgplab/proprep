"""
Force Field Parameter Reader

Reads bond, angle, and dihedral parameters from AMBER force field files.
Supports both .dat (main parameter files) and .frcmod (modification files).

© 2024 ProPrep Developer. All rights reserved.
"""

import os
import logging
from typing import Dict, Optional, Tuple, List
from dataclasses import dataclass


@dataclass
class BondParam:
    """Bond parameter from force field library."""
    atom_type1: str
    atom_type2: str
    force_constant: float  # kcal/mol/Å²
    eq_length: float       # Å
    source: str = ""


@dataclass
class AngleParam:
    """Angle parameter from force field library."""
    atom_type1: str
    atom_type2: str
    atom_type3: str
    force_constant: float  # kcal/mol/rad²
    eq_angle: float        # degrees
    source: str = ""


@dataclass
class MassParam:
    """Mass parameter from force field library."""
    atom_type: str
    mass: float            # amu
    comment: str = ""


@dataclass
class NonbondedParam:
    """Nonbonded (VDW) parameter from force field library."""
    atom_type: str
    radius: float          # Rmin/2 in Angstrom
    well_depth: float      # epsilon in kcal/mol
    comment: str = ""


@dataclass
class DihedralParam:
    """Dihedral parameter from force field library."""
    atom_type1: str
    atom_type2: str
    atom_type3: str
    atom_type4: str
    idivf: int
    pk: float              # barrier height (kcal/mol)
    phase: float           # phase angle (degrees)
    pn: float              # periodicity
    source: str = ""


class FFParameterReader:
    """
    Read force field parameters from AMBER .dat and .frcmod files.

    Parses BOND, ANGL, and DIHE sections from AMBER parameter files.
    Parameters are stored in lookup dictionaries for fast access.
    """

    # Force field to parameter file mapping
    FF_PARAM_MAP = {
        'ff19SB': 'parm19.dat',
        'ff14SB': 'parm10.dat',  # ff14SB uses parm10 base
        'ff14SBonlysc': 'parm10.dat',
        'ff10': 'parm10.dat',
        'ff03': 'parm03.dat',
        'ff99SB': 'parm99.dat',
        'ff99': 'parm99.dat',
        'ff94': 'parm94.dat',
        'gaff2': 'gaff2.dat',  # General AMBER Force Field 2
        'gaff': 'gaff.dat',    # General AMBER Force Field 1
    }

    # Additional frcmod files for specific force fields
    FF_FRCMOD_MAP = {
        'ff19SB': ['frcmod.ff19SB'],
        'ff14SB': ['frcmod.ff14SB'],
        'ff14SBonlysc': ['frcmod.ff14SBonlysc'],
    }

    def __init__(self, force_field: str, amberhome: Optional[str] = None, logger=None):
        """
        Initialize parameter reader.

        Args:
            force_field: Force field name (e.g., 'ff19SB', 'ff14SB')
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

        # Parameter storage
        self.mass_params: Dict[str, MassParam] = {}
        self.nonbonded_params: Dict[str, NonbondedParam] = {}
        self.bond_params: Dict[str, BondParam] = {}
        self.angle_params: Dict[str, AngleParam] = {}
        self.dihedral_params: Dict[str, List[DihedralParam]] = {}

        # Load parameters
        self._load_parameters()

    def _detect_amberhome(self) -> str:
        """Auto-detect AMBER installation from environment."""
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

    def _load_parameters(self):
        """Load all parameters from main .dat file, frcmod files, and GAFF2."""
        # Load main parameter file for primary forcefield
        param_file = self._get_param_file_path()
        if os.path.exists(param_file):
            self.logger.info(f"Loading parameters from: {param_file}")
            self._parse_parameter_file(param_file, source=self.force_field)
        else:
            self.logger.warning(f"Parameter file not found: {param_file}")

        # Load frcmod files for primary forcefield (these override main parameters)
        frcmod_files = self.FF_FRCMOD_MAP.get(self.force_field, [])
        for frcmod_name in frcmod_files:
            frcmod_path = os.path.join(
                self.amberhome, 'dat', 'leap', 'parm', frcmod_name
            )
            if os.path.exists(frcmod_path):
                self.logger.info(f"Loading frcmod: {frcmod_name}")
                self._parse_parameter_file(frcmod_path, source=frcmod_name)
            else:
                self.logger.warning(f"frcmod file not found: {frcmod_path}")

        # Always load GAFF2 parameters for ligand/small molecule support
        # GAFF2 parameters are loaded AFTER primary forcefield, so primary FF takes precedence
        # for any overlapping atom types (protein forcefields override GAFF for standard residues)
        if self.force_field != 'gaff2' and self.force_field != 'gaff':
            gaff2_path = os.path.join(self.amberhome, 'dat', 'leap', 'parm', 'gaff2.dat')
            if os.path.exists(gaff2_path):
                self.logger.info("Loading GAFF2 parameters for ligand support")
                # Parse with lower priority - existing parameters won't be overwritten
                self._parse_parameter_file(gaff2_path, source='gaff2', overwrite=False)
            else:
                self.logger.debug(f"GAFF2 file not found: {gaff2_path} (optional)")

        self.logger.info(
            f"Loaded {len(self.mass_params)} mass parameters, "
            f"{len(self.nonbonded_params)} nonbonded parameters, "
            f"{len(self.bond_params)} bond parameters, "
            f"{len(self.angle_params)} angle parameters, "
            f"{len(self.dihedral_params)} dihedral parameters"
        )

    def _get_param_file_path(self) -> str:
        """Get path to main parameter .dat file."""
        param_file = self.FF_PARAM_MAP.get(self.force_field)
        if not param_file:
            raise ValueError(
                f"Unknown force field: {self.force_field}\n"
                f"Supported force fields: {list(self.FF_PARAM_MAP.keys())}"
            )

        return os.path.join(self.amberhome, 'dat', 'leap', 'parm', param_file)

    def _parse_parameter_file(self, file_path: str, source: str, overwrite: bool = True):
        """
        Parse AMBER parameter file (.dat or .frcmod).

        Two formats:
        1. .dat files (parm10.dat, etc.): Implicit sections without headers
           - Mass parameters first
           - Blank line + atom type list
           - Bond parameters (TYPE1-TYPE2 format)
           - Angle parameters (TYPE1-TYPE2-TYPE3 format)
           - Dihedral parameters (TYPE1-TYPE2-TYPE3-TYPE4 format)
           - Nonbonded parameters (indented lines)

        2. .frcmod files: Explicit section headers (BOND, ANGL, DIHE, etc.)

        Args:
            file_path: Path to parameter file
            source: Source identifier (e.g., 'ff19SB', 'frcmod.ff14SB')
            overwrite: If False, don't overwrite existing parameters (for GAFF2 loading)
        """
        is_frcmod = file_path.endswith('.frcmod')

        # Some .dat files (e.g., amberdyes.dat) use explicit section headers
        # like frcmod files. Detect this by checking for a MASS header.
        if not is_frcmod:
            with open(file_path, 'r') as f:
                for line in f:
                    stripped = line.strip()
                    if not stripped or stripped.startswith('#'):
                        continue
                    # Skip title line (first non-empty, non-comment line)
                    # Check second non-empty line for section header
                    break
                for line in f:
                    stripped = line.strip()
                    if stripped in ('MASS', 'BOND', 'ANGL', 'ANGLE', 'NONB', 'NONBON', 'DIHE', 'DIHEDRAL'):
                        is_frcmod = True
                    break

        if is_frcmod:
            self._parse_frcmod_file(file_path, source, overwrite)
        else:
            self._parse_dat_file(file_path, source, overwrite)

    def _parse_frcmod_file(self, file_path: str, source: str, overwrite: bool = True):
        """Parse .frcmod file with explicit section headers."""
        section = None

        with open(file_path, 'r') as f:
            for line in f:
                line_stripped = line.strip()

                # Skip empty lines and comments
                if not line_stripped or line_stripped.startswith('#'):
                    continue

                # Check for section headers
                if line_stripped == 'MASS':
                    section = 'MASS'
                    continue
                elif line_stripped in ['NONB', 'NONBON']:
                    section = 'NONB'
                    continue
                elif line_stripped == 'BOND':
                    section = 'BOND'
                    continue
                elif line_stripped in ['ANGLE', 'ANGL']:
                    section = 'ANGL'
                    continue
                elif line_stripped in ['DIHEDRAL', 'DIHE']:
                    section = 'DIHE'
                    continue
                elif line_stripped in ['IMPROPER', 'IMPR']:
                    section = None  # Skip impropers
                    continue
                elif line_stripped.startswith('LJEDIT') or line_stripped.startswith('NBON'):
                    section = None
                    continue

                # Parse parameter lines based on current section
                if section == 'MASS':
                    self._parse_mass_line(line, source, overwrite)
                elif section == 'NONB':
                    self._parse_nonb_line(line, source, overwrite)
                elif section == 'BOND':
                    self._parse_bond_line(line, source, overwrite)
                elif section == 'ANGL':
                    self._parse_angle_line(line, source, overwrite)
                elif section == 'DIHE':
                    self._parse_dihedral_line(line, source, overwrite)

    def _parse_dat_file(self, file_path: str, source: str, overwrite: bool = True):
        """
        Parse .dat file with implicit sections.

        AMBER .dat files use fixed-width format for bonded parameters:
        - Bonds: FORMAT(A2,1X,A2,2F10.2) - chars 0-1, 3-4 are atom types
        - Angles: FORMAT(A2,1X,A2,1X,A2,2F10.2) - chars 0-1, 3-4, 6-7 are atom types
        - Dihedrals: FORMAT(A2,1X,A2,1X,A2,1X,A2,I4,3F15.2) - chars 0-1, 3-4, 6-7, 9-10

        Format heuristics:
        - Lines with dash at position 2: bonded parameters (bond/angle/dihedral)
        - Lines starting with spaces: Nonbonded parameters
        - Other lines: Mass parameters or VDW equivalences
        """
        # Collect equivalence lines before "MOD4 RE" for VDW parameters
        equivalences = []

        with open(file_path, 'r') as f:
            lines = f.readlines()

        # Find MOD4 RE marker (marks VDW equivalence section)
        mod4_idx = None
        for i, line in enumerate(lines):
            if 'MOD4' in line and 'RE' in line:
                mod4_idx = i
                break

        # Parse equivalence lines (3-4 lines before MOD4 RE)
        if mod4_idx is not None:
            for i in range(max(0, mod4_idx - 4), mod4_idx):
                line = lines[i]
                stripped = line.strip()
                # Equivalence lines: not indented, multiple atom types separated by spaces
                if stripped and not line.startswith(' ') and 'MOD4' not in stripped:
                    atom_types = stripped.split()
                    # Must have at least 2 types and all should be 1-2 chars
                    if len(atom_types) >= 2 and all(len(t) <= 2 for t in atom_types):
                        equivalences.append(atom_types)

        # Parse all parameter lines
        for line in lines:
            stripped = line.strip()

            # Skip empty lines and comments
            if not stripped or len(stripped) < 3:
                continue

            # Skip MOD4 RE and END lines
            if 'MOD4' in stripped or 'END' in stripped:
                continue

            # Nonbonded parameters (lines starting with spaces, no dash at pos 2)
            if line.startswith('  ') and (len(line) < 3 or line[2] != '-'):
                self._parse_nonb_line(line, source, overwrite)
                continue

            # Bonded parameters: check for dash at positions 2, 5, 8
            # Position 2: always present for bonds/angles/dihedrals
            # Position 5: present for angles/dihedrals
            # Position 8: present for dihedrals only
            if len(line) >= 3 and line[2] == '-':
                # Count dashes to determine type
                if len(line) >= 9 and line[8] == '-':
                    # Has dash at position 8 → dihedral (4 atom types)
                    self._parse_dihedral_line(line, source, overwrite)
                elif len(line) >= 6 and line[5] == '-':
                    # Has dash at position 5 but not 8 → angle (3 atom types)
                    self._parse_angle_line(line, source, overwrite)
                else:
                    # Has dash at position 2 only → bond (2 atom types)
                    self._parse_bond_line(line, source, overwrite)
                continue

            # Mass parameters (TYPE MASS [POLARIZABILITY] [COMMENT])
            # Check if it looks like a mass line: first field is 1-2 chars, second is float
            parts = stripped.split()
            if len(parts) >= 2 and len(parts[0]) <= 3:
                try:
                    float(parts[1])  # Second field should be mass
                    self._parse_mass_line(line, source, overwrite)
                except ValueError:
                    pass  # Not a mass line or equivalence line

        # Apply VDW equivalences: copy parameters from reference type to equivalent types
        for equiv_list in equivalences:
            if not equiv_list:
                continue

            ref_type = equiv_list[0]  # First type is the reference
            equiv_types = equiv_list[1:]  # Rest are equivalent types

            # Get VDW parameters for reference type (keys are already stripped)
            if ref_type in self.nonbonded_params:
                ref_param = self.nonbonded_params[ref_type]

                # Copy to all equivalent types
                for equiv_type in equiv_types:
                    # Create new parameter for equivalent type
                    equiv_param = NonbondedParam(
                        atom_type=equiv_type,
                        radius=ref_param.radius,
                        well_depth=ref_param.well_depth,
                        comment=f"equivalent to {ref_type}"
                    )
                    self.nonbonded_params[equiv_type] = equiv_param

                    self.logger.debug(
                        f"VDW equivalence: {equiv_type} -> {ref_type} "
                        f"(R={ref_param.radius}, eps={ref_param.well_depth})"
                    )

    def _parse_mass_line(self, line: str, source: str, overwrite: bool = True):
        """
        Parse mass parameter line.

        Format: ATOMTYPE  MASS  COMMENT
        Example: C  12.01  0.616  !  sp2 C carbonyl group
        """
        parts = line.split()
        if len(parts) < 2:
            return

        try:
            atom_type = parts[0].strip()

            # Skip if already exists and overwrite=False
            if not overwrite and atom_type in self.mass_params:
                return

            mass = float(parts[1])

            # Extract comment (everything after the mass, skip polarizability)
            # Format can be: TYPE MASS POLARIZABILITY COMMENT
            # or: TYPE MASS COMMENT
            comment_parts = []
            for i in range(2, len(parts)):
                # Skip polarizability value (typically a decimal number)
                try:
                    float(parts[i])
                    continue
                except ValueError:
                    # This and everything after is comment
                    comment_parts = parts[i:]
                    break

            comment = ' '.join(comment_parts).replace('!', '').strip()

            mass_param = MassParam(atom_type, mass, comment)
            self.mass_params[atom_type] = mass_param

        except (ValueError, IndexError) as e:
            self.logger.debug(f"Failed to parse mass line: {line.strip()}")

    def _parse_nonb_line(self, line: str, source: str, overwrite: bool = True):
        """
        Parse nonbonded parameter line.

        Format: ATOMTYPE  RADIUS  EPSILON  COMMENT
        Example:   C   1.9080  0.1094
        """
        parts = line.split()
        if len(parts) < 3:
            return

        try:
            atom_type = parts[0].strip()

            # Skip if already exists and overwrite=False
            if not overwrite and atom_type in self.nonbonded_params:
                return

            radius = float(parts[1])      # Rmin/2
            well_depth = float(parts[2])  # epsilon

            # Extract comment (everything after epsilon)
            comment = ' '.join(parts[3:]).strip() if len(parts) > 3 else ""

            nonb_param = NonbondedParam(atom_type, radius, well_depth, comment)
            self.nonbonded_params[atom_type] = nonb_param

        except (ValueError, IndexError) as e:
            self.logger.debug(f"Failed to parse nonb line: {line.strip()}")

    def _parse_bond_line(self, line: str, source: str, overwrite: bool = True):
        """
        Parse bond parameter line.

        Fixed-width format: FORMAT(A2,1X,A2,2F10.2)
        - Chars 0-1: First atom type
        - Char 2: Space (dash in visualization)
        - Chars 3-4: Second atom type
        - Chars 5-14: Force constant
        - Chars 15-24: Equilibrium length

        Example: "C -N   490.0    1.335"
        """
        try:
            # Parse using fixed-width columns
            if len(line) < 15:
                return

            type1 = line[0:2].strip()
            type2 = line[3:5].strip()

            # Parse numeric values - extract from the rest of the line
            values_str = line[5:].strip()
            values = values_str.split()

            if len(values) < 2:
                return

            force_constant = float(values[0])
            eq_length = float(values[1])

            if not type1 or not type2:
                return

            # Skip if already exists and overwrite=False
            key1 = f"{type1}-{type2}"
            key2 = f"{type2}-{type1}"
            if not overwrite and key1 in self.bond_params:
                return

            # Store both orderings
            bond_param = BondParam(type1, type2, force_constant, eq_length, source)
            self.bond_params[key1] = bond_param
            self.bond_params[key2] = bond_param

        except (ValueError, IndexError) as e:
            self.logger.debug(f"Failed to parse bond line: {line.strip()}")

    def _parse_angle_line(self, line: str, source: str, overwrite: bool = True):
        """
        Parse angle parameter line.

        Fixed-width format: FORMAT(A2,1X,A2,1X,A2,2F10.2)
        - Chars 0-1: First atom type
        - Char 2: Space (dash in visualization)
        - Chars 3-4: Second atom type
        - Char 5: Space (dash in visualization)
        - Chars 6-7: Third atom type
        - Chars 8-17: Force constant
        - Chars 18-27: Equilibrium angle

        Example: "C -C -O     80.0      120.00"
        """
        try:
            # Parse using fixed-width columns
            if len(line) < 18:
                return

            type1 = line[0:2].strip()
            type2 = line[3:5].strip()
            type3 = line[6:8].strip()

            # Parse numeric values - extract from the rest of the line
            values_str = line[8:].strip()
            values = values_str.split()

            if len(values) < 2:
                return

            force_constant = float(values[0])
            eq_angle = float(values[1])

            if not type1 or not type2 or not type3:
                return

            # Skip if already exists and overwrite=False
            key1 = f"{type1}-{type2}-{type3}"
            key2 = f"{type3}-{type2}-{type1}"
            if not overwrite and key1 in self.angle_params:
                return

            # Store both orderings (reverse order is equivalent)
            angle_param = AngleParam(type1, type2, type3, force_constant, eq_angle, source)
            self.angle_params[key1] = angle_param
            self.angle_params[key2] = angle_param

        except (ValueError, IndexError) as e:
            self.logger.debug(f"Failed to parse angle line: {line.strip()}")

    def _parse_dihedral_line(self, line: str, source: str, overwrite: bool = True):
        """
        Parse dihedral parameter line.

        Fixed-width format: FORMAT(A2,1X,A2,1X,A2,1X,A2,I4,3F15.2)
        - Chars 0-1: First atom type
        - Char 2: Space (dash in visualization)
        - Chars 3-4: Second atom type
        - Char 5: Space (dash in visualization)
        - Chars 6-7: Third atom type
        - Char 8: Space (dash in visualization)
        - Chars 9-10: Fourth atom type
        - Chars 11-14: IDIVF (I4)
        - Chars 15-29: PK (F15.2)
        - Chars 30-44: PHASE (F15.2)
        - Chars 45-59: PN (F15.2)

        Example: "C -N -CX-C    1    0.00          0.0            -4."
        """
        try:
            # Parse using fixed-width columns
            if len(line) < 15:
                return

            type1 = line[0:2].strip()
            type2 = line[3:5].strip()
            type3 = line[6:8].strip()
            type4 = line[9:11].strip()

            # Parse numeric values - extract from the rest of the line
            values_str = line[11:].strip()
            values = values_str.split()

            if len(values) < 4:
                return

            idivf = int(values[0])
            pk = float(values[1])
            phase = float(values[2])
            pn = float(values[3])

            if not type1 or not type2 or not type3 or not type4:
                return

            # Skip if already exists and overwrite=False
            key1 = f"{type1}-{type2}-{type3}-{type4}"
            key2 = f"{type4}-{type3}-{type2}-{type1}"
            if not overwrite and key1 in self.dihedral_params:
                return

            # Store dihedral (note: dihedrals can have multiple terms)
            dihedral_param = DihedralParam(
                type1, type2, type3, type4, idivf, pk, phase, pn, source
            )

            # Dihedrals can have multiple terms, so store as list
            if key1 not in self.dihedral_params:
                self.dihedral_params[key1] = []
            if key2 not in self.dihedral_params:
                self.dihedral_params[key2] = []

            self.dihedral_params[key1].append(dihedral_param)
            self.dihedral_params[key2].append(dihedral_param)

        except (ValueError, IndexError) as e:
            self.logger.debug(f"Failed to parse dihedral line: {line.strip()}")

    def get_mass_parameter(self, atom_type: str) -> Optional[MassParam]:
        """
        Get mass parameter for atom type.

        Args:
            atom_type: Atom type

        Returns:
            MassParam or None if not found
        """
        return self.mass_params.get(atom_type)

    def get_nonbonded_parameter(self, atom_type: str) -> Optional[NonbondedParam]:
        """
        Get nonbonded (VDW) parameter for atom type.

        Args:
            atom_type: Atom type

        Returns:
            NonbondedParam or None if not found
        """
        return self.nonbonded_params.get(atom_type)

    def get_bond_parameter(self, type1: str, type2: str) -> Optional[BondParam]:
        """
        Get bond parameter for atom type pair.

        Args:
            type1: First atom type
            type2: Second atom type

        Returns:
            BondParam or None if not found
        """
        # Strip whitespace from atom types (fingerprint pads to 2 chars with spaces)
        type1 = type1.strip()
        type2 = type2.strip()
        key = f"{type1}-{type2}"
        return self.bond_params.get(key)

    def get_angle_parameter(self, type1: str, type2: str, type3: str) -> Optional[AngleParam]:
        """
        Get angle parameter for atom type triple.

        Args:
            type1: First atom type
            type2: Second atom type (central)
            type3: Third atom type

        Returns:
            AngleParam or None if not found
        """
        # Strip whitespace from atom types (fingerprint pads to 2 chars with spaces)
        type1 = type1.strip()
        type2 = type2.strip()
        type3 = type3.strip()
        key = f"{type1}-{type2}-{type3}"
        return self.angle_params.get(key)

    def get_dihedral_parameters(self, type1: str, type2: str, type3: str,
                               type4: str) -> List[DihedralParam]:
        """
        Get dihedral parameters for atom type quadruple.

        Note: Dihedrals can have multiple terms (Fourier series).

        Args:
            type1: First atom type
            type2: Second atom type
            type3: Third atom type
            type4: Fourth atom type

        Returns:
            List of DihedralParam (can be empty if not found)
        """
        key = f"{type1}-{type2}-{type3}-{type4}"
        return self.dihedral_params.get(key, [])

    def get_statistics(self) -> Dict:
        """Get statistics about loaded parameters."""
        # Count unique dihedral atom type combinations
        unique_dihedrals = len(set(
            f"{d.atom_type1}-{d.atom_type2}-{d.atom_type3}-{d.atom_type4}"
            for dihedrals in self.dihedral_params.values()
            for d in dihedrals
        )) // 2  # Divide by 2 because we store both orderings

        return {
            'force_field': self.force_field,
            'n_masses': len(self.mass_params),
            'n_nonbonded': len(self.nonbonded_params),
            'n_bonds': len(self.bond_params) // 2,  # Divide by 2 (both orderings stored)
            'n_angles': len(self.angle_params) // 2,
            'n_dihedrals': unique_dihedrals,
        }

    def __repr__(self):
        stats = self.get_statistics()
        return (
            f"FFParameterReader("
            f"ff={stats['force_field']}, "
            f"bonds={stats['n_bonds']}, "
            f"angles={stats['n_angles']}, "
            f"dihedrals={stats['n_dihedrals']})"
        )
