"""
Unified Force Field Parsers

This module provides the SINGLE set of parsing functions for all AMBER force field
file formats. All code should use these parsers instead of implementing their own.

Consolidates parsing logic from:
- ForcefieldDataCollector (leaprc, lib/off files)
- ForceFieldLibrary (mol2 files)
- FFParameterReader (parmXX.dat, frcmod files)

File Formats Supported:
- leaprc: Force field loading scripts (source, loadOff, loadAmberParams directives)
- lib/off: AMBER library files (atom definitions with types and charges)
- mol2: Tripos mol2 files (atom coordinates, types, charges)
- parmXX.dat: Base AMBER parameter files (mass, bond, angle, dihedral, VDW)
- frcmod: Force field modification files (parameter overrides)

Author: ProPrep Development Team
"""

import re
import os
from pathlib import Path
from typing import List, Dict, Optional, Tuple, NamedTuple
from dataclasses import dataclass, field

from .ff_types import (
    MassParameter,
    BondParameter,
    AngleParameter,
    DihedralParameter,
    ImproperParameter,
    NonbondedParameter,
)


# =============================================================================
# Return Types for Parsers
# =============================================================================

@dataclass
class AtomDefinition:
    """
    Atom definition from lib or mol2 file.

    Represents a single atom's type and charge assignment within a residue.
    """
    residue_name: str
    atom_name: str
    atom_type: str
    charge: float
    source_file: str = ""

    @property
    def key(self) -> Tuple[str, str]:
        """Lookup key: (residue_name, atom_name)."""
        return (self.residue_name, self.atom_name)


@dataclass
class LeaprcContents:
    """Contents extracted from a leaprc file."""
    # Files to load
    lib_files: List[str] = field(default_factory=list)      # loadOff, loadAmberPrep
    frcmod_files: List[str] = field(default_factory=list)   # loadAmberParams (frcmod.*)
    parm_files: List[str] = field(default_factory=list)     # loadAmberParams (parm*.dat)
    mol2_files: List[str] = field(default_factory=list)     # loadMol2
    source_files: List[str] = field(default_factory=list)   # source (other leaprc files)

    # Residue mappings from addPdbResMap
    pdb_res_map: Dict[str, str] = field(default_factory=dict)  # {PDB_name: FF_name}

    # Bare LEaP unit aliases, e.g. `HOH = TP3` / `WAT = TP3` / `HIS = HIE`.
    # tleap resolves these by variable assignment (the alias name is bound to
    # the same UNIT object), not via addPdbResMap. Kept separate so the loader
    # can gate registration on the target unit actually being loaded.
    unit_aliases: Dict[str, str] = field(default_factory=dict)  # {alias: unit}

    # Source file
    source_path: Optional[Path] = None


@dataclass
class ParmContents:
    """Contents extracted from a parmXX.dat file."""
    mass_parameters: Dict[str, MassParameter] = field(default_factory=dict)
    bond_parameters: Dict[str, BondParameter] = field(default_factory=dict)
    angle_parameters: Dict[str, AngleParameter] = field(default_factory=dict)
    dihedral_parameters: Dict[str, List[DihedralParameter]] = field(default_factory=dict)
    improper_parameters: List[ImproperParameter] = field(default_factory=list)
    nonbonded_parameters: Dict[str, NonbondedParameter] = field(default_factory=dict)
    # VDW equivalences: atom types that share VDW parameters
    vdw_equivalences: List[List[str]] = field(default_factory=list)
    source_file: str = ""


@dataclass
class FrcmodContents:
    """Contents extracted from a frcmod file."""
    mass_parameters: Dict[str, MassParameter] = field(default_factory=dict)
    bond_parameters: Dict[str, BondParameter] = field(default_factory=dict)
    angle_parameters: Dict[str, AngleParameter] = field(default_factory=dict)
    dihedral_parameters: Dict[str, List[DihedralParameter]] = field(default_factory=dict)
    improper_parameters: List[ImproperParameter] = field(default_factory=list)
    nonbonded_parameters: Dict[str, NonbondedParameter] = field(default_factory=dict)
    source_file: str = ""


# =============================================================================
# Leaprc Parser
# =============================================================================

def parse_leaprc(path: Path) -> LeaprcContents:
    """
    Parse a leaprc file to extract loading directives.

    Extracts:
    - source directives (chained leaprc files)
    - loadOff/loadAmberPrep directives (lib files)
    - loadAmberParams directives (frcmod and parm.dat files)
    - loadMol2 directives
    - addPdbResMap directives (residue name mappings)
    - bare unit-alias assignments (e.g. ``HOH = TP3``)

    Args:
        path: Path to leaprc file

    Returns:
        LeaprcContents with all extracted information
    """
    contents = LeaprcContents(source_path=path)

    try:
        with open(path, 'r') as f:
            text = f.read()
    except (FileNotFoundError, IOError):
        return contents

    # Patterns for different load commands
    patterns = {
        'loadOff': r'loadOff\s+([^\s#]+)',
        'loadoff': r'loadoff\s+([^\s#]+)',
        'loadAmberParams': r'loadAmberParams\s+([^\s#]+)',
        'loadamberparams': r'loadamberparams\s+([^\s#]+)',
        'source': r'source\s+([^\s#]+)',
        'loadAmberPrep': r'loadAmberPrep\s+([^\s#]+)',
        'loadMol2': r'loadMol2\s+([^\s#]+)',
        'loadPdb': r'loadPdb\s+([^\s#]+)',
    }

    for cmd, pattern in patterns.items():
        matches = re.finditer(pattern, text, re.IGNORECASE | re.MULTILINE)
        for match in matches:
            filename = match.group(1).strip()

            if cmd.lower() == 'source':
                contents.source_files.append(filename)
            elif cmd.lower() in ('loadoff', 'loadamberprep'):
                contents.lib_files.append(filename)
            elif cmd.lower() == 'loadamberparams':
                # Distinguish between frcmod and parm.dat files
                if 'frcmod' in filename.lower():
                    contents.frcmod_files.append(filename)
                elif filename.endswith('.dat'):
                    contents.parm_files.append(filename)
                else:
                    # Default: treat as frcmod
                    contents.frcmod_files.append(filename)
            elif cmd.lower() == 'loadmol2':
                contents.mol2_files.append(filename)

    # Parse addPdbResMap directives
    # Format: addPdbResMap { { 0 "PDB_NAME" "FF_NAME" } }
    pdb_res_pattern = r'addPdbResMap\s*\{\s*\{\s*\d+\s+"([^"]+)"\s+"([^"]+)"\s*\}\s*\}'
    for match in re.finditer(pdb_res_pattern, text, re.IGNORECASE):
        pdb_name = match.group(1)
        ff_name = match.group(2)
        contents.pdb_res_map[pdb_name] = ff_name

    # Parse bare unit-alias assignments, e.g. `HOH = TP3`, `WAT = TP3`,
    # `HIS = HIE`. This is the LEaP variable-assignment form (the alias name is
    # bound to the same UNIT), used by every leaprc.water.* to map the
    # crystallographic name HOH onto the selected water model. It is NOT
    # addPdbResMap, so it needs its own pattern. The full-line anchors plus an
    # optional trailing comment keep this from matching function-call RHSs or
    # commented-out lines; across the shipped AMBER leaprc corpus every match is
    # a genuine residue alias (water models, HIS defaults, spin-label MTN).
    unit_alias_pattern = r'^[ \t]*([A-Za-z0-9+\-]+)[ \t]*=[ \t]*([A-Za-z0-9+\-]+)[ \t]*(?:#.*)?$'
    for match in re.finditer(unit_alias_pattern, text, re.MULTILINE):
        alias_name = match.group(1)
        unit_name = match.group(2)
        contents.unit_aliases[alias_name] = unit_name

    return contents


# =============================================================================
# Lib/OFF File Parser
# =============================================================================

def parse_lib_file(path: Path, source_name: str = "") -> List[AtomDefinition]:
    """
    Parse an AMBER lib (OFF format) file.

    Extracts atom definitions from unit.atoms tables.

    Format:
        !entry.RESIDUE.unit.atoms table  str name  str type  int typex ...
        "CA" "CT" 0 1 131072 1 6 0.033700

    Args:
        path: Path to lib file
        source_name: Source identifier for tracking

    Returns:
        List of AtomDefinition objects
    """
    atoms = []
    source = source_name or path.name

    try:
        with open(path, 'r') as f:
            content = f.read()
    except (FileNotFoundError, IOError):
        return atoms

    # Pattern to find unit.atoms tables
    # Captures: residue name and the atoms section
    atom_table_pattern = r'!entry\.([^.]+)\.unit\.atoms\s+table[^\n]*\n((?:^[^!].*\n)*)'

    for match in re.finditer(atom_table_pattern, content, re.MULTILINE):
        residue_name = match.group(1)
        atoms_section = match.group(2)

        # Parse each atom line
        # Format: "atom_name" "atom_type" typex resx flags seq elmnt charge
        # The charge is the 8th field (0-indexed: field 7)
        atom_line_pattern = r'"([^"]+)"\s+"([^"]+)"\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+([-\d.]+)'

        for atom_match in re.finditer(atom_line_pattern, atoms_section):
            atom_name = atom_match.group(1)
            atom_type = atom_match.group(2)
            try:
                charge = float(atom_match.group(3))
            except ValueError:
                charge = 0.0

            atoms.append(AtomDefinition(
                residue_name=residue_name,
                atom_name=atom_name,
                atom_type=atom_type,
                charge=charge,
                source_file=source
            ))

    return atoms


# =============================================================================
# Mol2 File Parser
# =============================================================================

def parse_mol2_file(path: Path, source_name: str = "") -> List[AtomDefinition]:
    """
    Parse a Tripos mol2 file.

    Extracts atom definitions from @<TRIPOS>ATOM sections.

    Format:
        @<TRIPOS>ATOM
        atom_id atom_name x y z atom_type res_id res_name charge
        1 N -0.123 0.456 0.789 N 1 ALA -0.4157

    Args:
        path: Path to mol2 file
        source_name: Source identifier for tracking

    Returns:
        List of AtomDefinition objects
    """
    atoms = []
    source = source_name or path.name

    try:
        # Try using pymsmt if available (more robust)
        from pymsmt.mol.mol2io import get_atominfo

        mol, atids, resids = get_atominfo(str(path))

        for resid in resids:
            resname = mol.residues[resid].resname
            for atomid in mol.residues[resid].resconter:
                atom = mol.atoms[atomid]
                atoms.append(AtomDefinition(
                    residue_name=resname,
                    atom_name=atom.atname,
                    atom_type=atom.atomtype.strip(),
                    charge=atom.charge,
                    source_file=source
                ))

        return atoms

    except ImportError:
        # Fallback: manual parsing
        pass

    # Manual mol2 parsing
    try:
        with open(path, 'r') as f:
            lines = f.readlines()
    except (FileNotFoundError, IOError):
        return atoms

    in_atom_section = False

    for line in lines:
        line_stripped = line.strip()

        if line_stripped == '@<TRIPOS>ATOM':
            in_atom_section = True
            continue
        elif line_stripped.startswith('@<TRIPOS>'):
            in_atom_section = False
            continue

        if in_atom_section and line_stripped:
            parts = line_stripped.split()
            if len(parts) >= 9:
                # atom_id atom_name x y z atom_type res_id res_name charge
                atom_name = parts[1]
                atom_type = parts[5]
                resname = parts[7]
                try:
                    charge = float(parts[8])
                except (ValueError, IndexError):
                    charge = 0.0

                atoms.append(AtomDefinition(
                    residue_name=resname,
                    atom_name=atom_name,
                    atom_type=atom_type.strip(),
                    charge=charge,
                    source_file=source
                ))

    return atoms


# =============================================================================
# Parm.dat File Parser
# =============================================================================

def parse_parm_dat(path: Path, source_name: str = "") -> ParmContents:
    """
    Parse an AMBER parmXX.dat base parameter file.

    These files use fixed-width format with implicit sections:
    - Mass parameters (at start)
    - VDW equivalence lines (before MOD4 RE marker)
    - Bond parameters (TYPE1-TYPE2 at positions 0-1, 3-4)
    - Angle parameters (TYPE1-TYPE2-TYPE3 at positions 0-1, 3-4, 6-7)
    - Dihedral parameters (4 atom types)
    - Nonbonded parameters (indented lines)

    Args:
        path: Path to parm.dat file
        source_name: Source identifier for tracking

    Returns:
        ParmContents with all parsed parameters
    """
    contents = ParmContents(source_file=source_name or path.name)

    try:
        with open(path, 'r') as f:
            lines = f.readlines()
    except (FileNotFoundError, IOError):
        return contents

    source = contents.source_file

    # Find MOD4 RE marker (marks VDW equivalence section)
    mod4_idx = None
    for i, line in enumerate(lines):
        if 'MOD4' in line and 'RE' in line:
            mod4_idx = i
            break

    # Parse VDW equivalence lines (3-4 lines before MOD4 RE)
    if mod4_idx is not None:
        for i in range(max(0, mod4_idx - 4), mod4_idx):
            line = lines[i]
            stripped = line.strip()
            # Equivalence lines: not indented, multiple atom types separated by spaces
            if stripped and not line.startswith(' ') and 'MOD4' not in stripped:
                atom_types = stripped.split()
                # Must have at least 2 types and all should be 1-2 chars
                if len(atom_types) >= 2 and all(len(t) <= 2 for t in atom_types):
                    contents.vdw_equivalences.append(atom_types)

    # Parse all parameter lines
    for line in lines:
        stripped = line.strip()

        # Skip empty lines and short lines
        if not stripped or len(stripped) < 3:
            continue

        # Skip MOD4 RE and END markers
        if 'MOD4' in stripped or stripped == 'END':
            continue

        # Nonbonded parameters (lines starting with spaces, no dash at pos 2)
        if line.startswith('  ') and (len(line) < 3 or line[2] != '-'):
            param = _parse_nonbonded_line(line, source)
            if param:
                contents.nonbonded_parameters[param.atom_type] = param
            continue

        # Bonded parameters: check for dash at positions 2, 5, 8
        if len(line) >= 3 and line[2] == '-':
            if len(line) >= 9 and line[8] == '-':
                # Dihedral (4 atom types)
                param = _parse_dihedral_line(line, source)
                if param:
                    key = param.key
                    if key not in contents.dihedral_parameters:
                        contents.dihedral_parameters[key] = []
                    contents.dihedral_parameters[key].append(param)
            elif len(line) >= 6 and line[5] == '-':
                # Angle (3 atom types)
                param = _parse_angle_line(line, source)
                if param:
                    contents.angle_parameters[param.key] = param
            else:
                # Bond (2 atom types)
                param = _parse_bond_line(line, source)
                if param:
                    contents.bond_parameters[param.key] = param
            continue

        # Mass parameters (TYPE MASS [POLARIZABILITY] [COMMENT])
        parts = stripped.split()
        if len(parts) >= 2 and len(parts[0]) <= 3:
            try:
                float(parts[1])  # Second field should be mass
                param = _parse_mass_line(line, source)
                if param:
                    contents.mass_parameters[param.atom_type] = param
            except ValueError:
                pass

    # Apply VDW equivalences: copy parameters from reference type to equivalent types
    for equiv_list in contents.vdw_equivalences:
        if not equiv_list:
            continue

        ref_type = equiv_list[0]  # First type is the reference
        equiv_types = equiv_list[1:]  # Rest are equivalent types

        if ref_type in contents.nonbonded_parameters:
            ref_param = contents.nonbonded_parameters[ref_type]

            for equiv_type in equiv_types:
                equiv_param = NonbondedParameter(
                    atom_type=equiv_type,
                    radius=ref_param.radius,
                    well_depth=ref_param.well_depth,
                    comment=f"equivalent to {ref_type}",
                    source=source
                )
                contents.nonbonded_parameters[equiv_type] = equiv_param

    return contents


# =============================================================================
# Frcmod File Parser
# =============================================================================

def parse_frcmod(path: Path, source_name: str = "") -> FrcmodContents:
    """
    Parse an AMBER frcmod (force field modification) file.

    These files use explicit section headers:
    - MASS
    - BOND
    - ANGLE/ANGL
    - DIHEDRAL/DIHE
    - IMPROPER/IMPR
    - NONBON/NONB

    Args:
        path: Path to frcmod file
        source_name: Source identifier for tracking

    Returns:
        FrcmodContents with all parsed parameters
    """
    contents = FrcmodContents(source_file=source_name or path.name)

    try:
        with open(path, 'r') as f:
            lines = f.readlines()
    except (FileNotFoundError, IOError):
        return contents

    source = contents.source_file
    section = None

    for line in lines:
        line_stripped = line.strip()

        # Skip empty lines and comments
        if not line_stripped or line_stripped.startswith('#'):
            continue

        # Check for section headers
        if line_stripped == 'MASS':
            section = 'MASS'
            continue
        elif line_stripped in ('NONB', 'NONBON'):
            section = 'NONB'
            continue
        elif line_stripped == 'BOND':
            section = 'BOND'
            continue
        elif line_stripped in ('ANGLE', 'ANGL'):
            section = 'ANGL'
            continue
        elif line_stripped in ('DIHEDRAL', 'DIHE'):
            section = 'DIHE'
            continue
        elif line_stripped in ('IMPROPER', 'IMPR'):
            section = 'IMPR'
            continue
        elif line_stripped.startswith('LJEDIT') or line_stripped.startswith('NBON'):
            section = None
            continue

        # Parse based on current section
        if section == 'MASS':
            param = _parse_mass_line(line, source)
            if param:
                contents.mass_parameters[param.atom_type] = param
        elif section == 'NONB':
            param = _parse_nonbonded_line(line, source)
            if param:
                contents.nonbonded_parameters[param.atom_type] = param
        elif section == 'BOND':
            param = _parse_bond_line(line, source)
            if param:
                contents.bond_parameters[param.key] = param
        elif section == 'ANGL':
            param = _parse_angle_line(line, source)
            if param:
                contents.angle_parameters[param.key] = param
        elif section == 'DIHE':
            param = _parse_dihedral_line(line, source)
            if param:
                key = param.key
                if key not in contents.dihedral_parameters:
                    contents.dihedral_parameters[key] = []
                contents.dihedral_parameters[key].append(param)
        elif section == 'IMPR':
            param = _parse_improper_line(line, source)
            if param:
                contents.improper_parameters.append(param)

    return contents


# =============================================================================
# Line Parsing Helpers
# =============================================================================

def _parse_mass_line(line: str, source: str) -> Optional[MassParameter]:
    """
    Parse a mass parameter line.

    Format: ATOMTYPE  MASS  [POLARIZABILITY]  [COMMENT]
    Example: C  12.01  0.616  !  sp2 C carbonyl group
    """
    parts = line.split()
    if len(parts) < 2:
        return None

    try:
        atom_type = parts[0].strip()
        mass = float(parts[1])

        # Extract polarizability if present
        polarizability = 0.0
        comment_start = 2
        if len(parts) >= 3:
            try:
                polarizability = float(parts[2])
                comment_start = 3
            except ValueError:
                pass

        # Everything after mass/polarizability is comment
        comment_parts = []
        for i in range(comment_start, len(parts)):
            if parts[i] == '!' or parts[i].startswith('!'):
                comment_parts = parts[i:]
                break
            comment_parts.append(parts[i])

        comment = ' '.join(comment_parts).replace('!', '').strip()

        return MassParameter(
            atom_type=atom_type,
            mass=mass,
            polarizability=polarizability,
            comment=comment,
            source=source
        )
    except (ValueError, IndexError):
        return None


def _parse_nonbonded_line(line: str, source: str) -> Optional[NonbondedParameter]:
    """
    Parse a nonbonded (VDW) parameter line.

    Format: ATOMTYPE  RADIUS  EPSILON  [COMMENT]
    Example:   C   1.9080  0.1094
    """
    parts = line.split()
    if len(parts) < 3:
        return None

    try:
        atom_type = parts[0].strip()
        radius = float(parts[1])      # Rmin/2
        well_depth = float(parts[2])  # epsilon

        comment = ' '.join(parts[3:]).strip() if len(parts) > 3 else ""

        return NonbondedParameter(
            atom_type=atom_type,
            radius=radius,
            well_depth=well_depth,
            comment=comment,
            source=source
        )
    except (ValueError, IndexError):
        return None


def _parse_bond_line(line: str, source: str) -> Optional[BondParameter]:
    """
    Parse a bond parameter line.

    Fixed-width format: FORMAT(A2,1X,A2,2F10.2)
    - Chars 0-1: First atom type
    - Char 2: Dash
    - Chars 3-4: Second atom type
    - Rest: Force constant and equilibrium length

    Example: C -N   490.0    1.335
    """
    try:
        if len(line) < 15:
            return None

        type1 = line[0:2].strip()
        type2 = line[3:5].strip()

        if not type1 or not type2:
            return None

        # Parse numeric values from rest of line
        values_str = line[5:].strip()
        values = values_str.split()

        if len(values) < 2:
            return None

        force_constant = float(values[0])
        eq_length = float(values[1])

        return BondParameter(
            type1=type1,
            type2=type2,
            force_constant=force_constant,
            eq_length=eq_length,
            source=source
        )
    except (ValueError, IndexError):
        return None


def _parse_angle_line(line: str, source: str) -> Optional[AngleParameter]:
    """
    Parse an angle parameter line.

    Fixed-width format: FORMAT(A2,1X,A2,1X,A2,2F10.2)
    - Chars 0-1: First atom type
    - Char 2: Dash
    - Chars 3-4: Second atom type (central)
    - Char 5: Dash
    - Chars 6-7: Third atom type
    - Rest: Force constant and equilibrium angle

    Example: C -C -O     80.0      120.00
    """
    try:
        if len(line) < 18:
            return None

        type1 = line[0:2].strip()
        type2 = line[3:5].strip()
        type3 = line[6:8].strip()

        if not type1 or not type2 or not type3:
            return None

        # Parse numeric values
        values_str = line[8:].strip()
        values = values_str.split()

        if len(values) < 2:
            return None

        force_constant = float(values[0])
        eq_angle = float(values[1])

        return AngleParameter(
            type1=type1,
            type2=type2,
            type3=type3,
            force_constant=force_constant,
            eq_angle=eq_angle,
            source=source
        )
    except (ValueError, IndexError):
        return None


def _parse_dihedral_line(line: str, source: str) -> Optional[DihedralParameter]:
    """
    Parse a dihedral parameter line.

    Fixed-width format: FORMAT(A2,1X,A2,1X,A2,1X,A2,I4,3F15.2)
    - Chars 0-1, 3-4, 6-7, 9-10: Four atom types
    - Chars 11-14: IDIVF (divider)
    - Rest: PK (barrier), PHASE, PN (periodicity)

    Example: C -N -CX-C    1    0.00          0.0            -4.
    """
    try:
        if len(line) < 15:
            return None

        type1 = line[0:2].strip()
        type2 = line[3:5].strip()
        type3 = line[6:8].strip()
        type4 = line[9:11].strip()

        if not type1 or not type2 or not type3 or not type4:
            return None

        # Parse numeric values
        values_str = line[11:].strip()
        values = values_str.split()

        if len(values) < 4:
            return None

        divider = int(values[0])
        barrier = float(values[1])
        phase = float(values[2])
        periodicity = float(values[3])

        return DihedralParameter(
            type1=type1,
            type2=type2,
            type3=type3,
            type4=type4,
            divider=divider,
            barrier=barrier,
            phase=phase,
            periodicity=periodicity,
            source=source
        )
    except (ValueError, IndexError):
        return None


def _parse_improper_line(line: str, source: str) -> Optional[ImproperParameter]:
    """
    Parse an improper dihedral parameter line.

    Format similar to dihedral but different interpretation.
    The first atom is typically the central atom.

    Example: X -X -C -O          10.5         180.          2.
    """
    try:
        if len(line) < 15:
            return None

        type1 = line[0:2].strip()
        type2 = line[3:5].strip()
        type3 = line[6:8].strip()
        type4 = line[9:11].strip()

        if not type1 or not type2 or not type3 or not type4:
            return None

        # Parse numeric values
        values_str = line[11:].strip()
        values = values_str.split()

        if len(values) < 3:
            return None

        barrier = float(values[0])
        phase = float(values[1])
        periodicity = int(float(values[2]))

        return ImproperParameter(
            type1=type1,
            type2=type2,
            type3=type3,
            type4=type4,
            barrier=barrier,
            phase=phase,
            periodicity=periodicity,
            source=source
        )
    except (ValueError, IndexError):
        return None


# =============================================================================
# File Path Resolution Utilities
# =============================================================================

def find_amber_file(filename: str, amberhome: Path, search_type: str = 'lib') -> Optional[Path]:
    """
    Find an AMBER data file given its name.

    Args:
        filename: Name of file to find (e.g., 'amino12.lib', 'frcmod.ff14SB')
        amberhome: Path to AMBER installation
        search_type: Type of file ('lib', 'frcmod', 'parm', 'leaprc')

    Returns:
        Full path to file, or None if not found
    """
    # Common search directories
    search_dirs = {
        'lib': [
            amberhome / 'dat' / 'leap' / 'lib',
            amberhome / 'dat' / 'leap' / 'parm',
        ],
        'frcmod': [
            amberhome / 'dat' / 'leap' / 'parm',
        ],
        'parm': [
            amberhome / 'dat' / 'leap' / 'parm',
        ],
        'leaprc': [
            amberhome / 'dat' / 'leap' / 'cmd',
        ],
        'mol2': [
            amberhome / 'dat' / 'pymsmt',
            amberhome / 'dat' / 'leap' / 'lib',
        ],
    }

    dirs = search_dirs.get(search_type, [])

    for search_dir in dirs:
        # Try exact filename
        path = search_dir / filename
        if path.exists():
            return path

        # Try with common prefixes/suffixes removed
        base = Path(filename).name
        path = search_dir / base
        if path.exists():
            return path

    return None


def resolve_leaprc_path(name: str, amberhome: Path) -> Optional[Path]:
    """
    Resolve a leaprc file name to its full path.

    Handles:
    - Full names: 'leaprc.protein.ff14SB'
    - Relative paths: 'oldff/leaprc.ff10'
    - Short names: 'ff14SB' -> 'leaprc.protein.ff14SB'

    Args:
        name: Leaprc file name or path
        amberhome: Path to AMBER installation

    Returns:
        Full path to leaprc file, or None if not found
    """
    leap_cmd_dir = amberhome / 'dat' / 'leap' / 'cmd'

    # Try exact path first
    if name.startswith('leaprc') or '/' in name:
        path = leap_cmd_dir / name
        if path.exists():
            return path

    # Try common patterns
    patterns = [
        f'leaprc.protein.{name}',
        f'leaprc.{name}',
        f'oldff/leaprc.{name}',
    ]

    for pattern in patterns:
        path = leap_cmd_dir / pattern
        if path.exists():
            return path

    return None
