"""
AMBER .lib File Parser

Parses AMBER library files to extract atom types and charges for residues.
Used by ONIOM atom typer to load custom forcefield parameters.
"""

import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


class AmberLibParser:
    """Parser for AMBER .lib (OFF format) library files"""

    def __init__(self, lib_file_path: str):
        """
        Initialize parser with .lib file path

        Args:
            lib_file_path: Path to AMBER .lib file
        """
        self.lib_file_path = Path(lib_file_path)
        self.residues: Dict[str, Dict[str, any]] = {}

    def parse(self) -> Dict[str, Dict]:
        """
        Parse the .lib file and extract all residue definitions

        Returns:
            Dictionary mapping residue names to their atom data:
            {
                "HCR": {
                    "atoms": [
                        {"name": "FE", "type": "M8", "charge": 0.1234, ...},
                        ...
                    ],
                    "connectivity": [...],
                    "positions": [...]
                }
            }
        """
        if not self.lib_file_path.exists():
            raise FileNotFoundError(f".lib file not found: {self.lib_file_path}")

        with open(self.lib_file_path, 'r') as f:
            lines = f.readlines()

        # First pass: Read the index array to get residue names
        in_index_array = False
        for i, line in enumerate(lines):
            if line.startswith('!!index array str'):
                in_index_array = True
                continue
            elif in_index_array:
                # Read residue names from the index array
                if line.strip().startswith('"'):
                    resname = line.strip().strip('"')
                    if resname:
                        self.residues[resname] = {
                            "atoms": [],
                            "connectivity": [],
                            "positions": [],
                            "connection_points": {"n_terminus": 0, "c_terminus": 0}
                        }
                elif line.startswith('!'):
                    # End of index array
                    in_index_array = False

        # Second pass: Parse atom data, connectivity, and residueconnect for each residue
        current_residue = None
        in_atoms_section = False
        in_positions_section = False

        # NEW: Parse connectivity and residueconnect sections
        self._parse_connectivity_sections(lines)
        self._parse_residueconnect_sections(lines)

        for line in lines:
            # Detect atoms section
            # Check if we're entering a new residue's atom section
            if line.startswith('!entry.') and '.unit.atoms table' in line:
                # Extract residue name from line like "!entry.DHO.unit.atoms table"
                match = re.match(r'!entry\.([^.]+)\.unit\.atoms', line)
                if match:
                    current_residue = match.group(1)
                    if current_residue in self.residues:
                        in_atoms_section = True
                continue

            # Parse atom data
            elif in_atoms_section and current_residue:
                if line.startswith('!entry'):
                    in_atoms_section = False
                    continue

                # Parse atom line: "NAME" "TYPE" typex resx flags seq elmnt charge
                match = re.match(r'\s*"([^"]+)"\s+"([^"]+)"\s+\d+\s+\d+\s+\d+\s+\d+\s+(-?\d+)\s+([-\d.]+)', line)
                if match:
                    atom_name = match.group(1)
                    atom_type = match.group(2)
                    element_num = int(match.group(3))
                    charge = float(match.group(4))

                    self.residues[current_residue]["atoms"].append({
                        "name": atom_name,
                        "type": atom_type,
                        "charge": charge,
                        "element_num": element_num
                    })

            # Detect positions section
            elif line.startswith('!entry.') and '.unit.positions table' in line:
                # Extract residue name from line like "!entry.DHO.unit.positions table"
                match = re.match(r'!entry\.([^.]+)\.unit\.positions', line)
                if match:
                    current_residue = match.group(1)
                    if current_residue in self.residues:
                        in_positions_section = True
                continue

            # Parse position data
            elif in_positions_section and current_residue:
                if line.startswith('!entry'):
                    in_positions_section = False
                    continue

                # Parse position line: x y z
                match = re.match(r'\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', line)
                if match:
                    x = float(match.group(1))
                    y = float(match.group(2))
                    z = float(match.group(3))
                    self.residues[current_residue]["positions"].append((x, y, z))

        logger.debug(f"Parsed {len(self.residues)} residues from {self.lib_file_path.name}")
        return self.residues

    def get_residue_atoms(self, resname: str) -> List[Dict]:
        """
        Get atom list for a specific residue

        Args:
            resname: Residue name (e.g., "HCR", "DHR")

        Returns:
            List of atom dictionaries with name, type, and charge
        """
        if resname not in self.residues:
            if not self.residues:
                # Try parsing if not done yet
                self.parse()

            if resname not in self.residues:
                logger.debug(f"Residue {resname} not found in {self.lib_file_path.name}")
                return []

        return self.residues[resname].get("atoms", [])

    def get_atom_type_and_charge(self, resname: str, atom_name: str) -> Optional[Tuple[str, float]]:
        """
        Get atom type and charge for a specific atom in a residue

        Args:
            resname: Residue name
            atom_name: Atom name

        Returns:
            Tuple of (atom_type, charge) or None if not found
        """
        atoms = self.get_residue_atoms(resname)

        for atom in atoms:
            if atom["name"] == atom_name:
                return (atom["type"], atom["charge"])

        return None

    def _parse_connectivity_sections(self, lines: List[str]):
        """
        Parse !entry.RES.unit.connectivity table sections

        Example input:
            !entry.ALA.unit.connectivity table  int atom1x  int atom2x  int flags
             1 2 1
             1 3 1
             3 4 1

        Stores as list of tuples: [(1, 2), (1, 3), (3, 4), ...]
        """
        current_residue = None
        in_section = False

        for line in lines:
            if '!entry.' in line and '.unit.connectivity table' in line:
                # Extract residue name: !entry.ALA.unit.connectivity
                match = re.match(r'!entry\.([^.]+)\.unit\.connectivity', line)
                if match:
                    current_residue = match.group(1)
                    if current_residue in self.residues:
                        in_section = True
                continue

            if in_section:
                if line.startswith('!entry'):
                    # End of section
                    in_section = False
                    current_residue = None
                    continue

                # Parse bond line: " 1 2 1" -> (1, 2)
                # Format: atom1_idx atom2_idx flags
                match = re.match(r'\s*(\d+)\s+(\d+)\s+\d+', line)
                if match and current_residue in self.residues:
                    atom1 = int(match.group(1))
                    atom2 = int(match.group(2))
                    self.residues[current_residue]["connectivity"].append((atom1, atom2))

    def _parse_residueconnect_sections(self, lines: List[str]):
        """
        Parse !entry.RES.unit.residueconnect table sections

        Example input:
            !entry.ALA.unit.residueconnect table  int c1x  int c2x  int c3x  int c4x  int c5x  int c6x
             1 9 0 0 0 0

        First value (c1x) is N-terminus connection point
        Second value (c2x) is C-terminus connection point
        0 = no connection point

        Stores as: {"n_terminus": 1, "c_terminus": 9}
        """
        current_residue = None
        in_section = False

        for line in lines:
            if '!entry.' in line and '.unit.residueconnect table' in line:
                # Extract residue name: !entry.ALA.unit.residueconnect
                match = re.match(r'!entry\.([^.]+)\.unit\.residueconnect', line)
                if match:
                    current_residue = match.group(1)
                    if current_residue in self.residues:
                        in_section = True
                continue

            if in_section:
                if line.startswith('!entry'):
                    # End of section
                    in_section = False
                    current_residue = None
                    continue

                # Parse connection points: " 1 9 0 0 0 0"
                # Format: c1x c2x c3x c4x c5x c6x (we only care about first two)
                match = re.match(r'\s*(\d+)\s+(\d+)', line)
                if match and current_residue in self.residues:
                    n_term = int(match.group(1))
                    c_term = int(match.group(2))
                    self.residues[current_residue]["connection_points"] = {
                        "n_terminus": n_term,
                        "c_terminus": c_term
                    }

    def get_connectivity_template(self, resname: str) -> List[Tuple[int, int]]:
        """
        Get connectivity template for a residue

        Args:
            resname: Residue name (e.g., "ALA", "NALA", "CALA")

        Returns:
            List of (atom1_idx, atom2_idx) tuples with 1-based indices
        """
        if resname not in self.residues:
            if not self.residues:
                # Try parsing if not done yet
                self.parse()

            if resname not in self.residues:
                logger.debug(f"Residue {resname} not found in {self.lib_file_path.name}")
                return []

        return self.residues[resname].get("connectivity", [])

    def get_connection_points(self, resname: str) -> Dict[str, int]:
        """
        Get N/C-terminal connection points for a residue

        Args:
            resname: Residue name (e.g., "ALA", "NALA", "CALA")

        Returns:
            Dict with "n_terminus" and "c_terminus" atom indices (0 = no connection)

        Example:
            For ALA: {"n_terminus": 1, "c_terminus": 9}  # N and C atoms
            For NALA: {"n_terminus": 0, "c_terminus": 11}  # Only C connects
            For CALA: {"n_terminus": 1, "c_terminus": 0}  # Only N connects
        """
        if resname not in self.residues:
            if not self.residues:
                # Try parsing if not done yet
                self.parse()

            if resname not in self.residues:
                logger.debug(f"Residue {resname} not found in {self.lib_file_path.name}")
                return {"n_terminus": 0, "c_terminus": 0}

        return self.residues[resname].get("connection_points", {"n_terminus": 0, "c_terminus": 0})


class AmberPrepParser:
    """
    Parser for AMBER .prep (internal coordinate) files.

    Extracts residue names, atom names, atom types, and partial charges.
    Provides the same interface as AmberLibParser for interchangeability.
    """

    def __init__(self, prep_file_path: str):
        self.prep_file_path = Path(prep_file_path)
        self.residues: Dict[str, Dict[str, any]] = {}

    def parse(self) -> Dict[str, Dict]:
        """Parse the .prep file and extract residue definitions."""
        if not self.prep_file_path.exists():
            raise FileNotFoundError(f".prep file not found: {self.prep_file_path}")

        with open(self.prep_file_path, 'r') as f:
            content = f.read()

        # A .prep file can contain multiple residues separated by DONE/STOP
        # Split into blocks between header and DONE
        blocks = re.split(r'\bSTOP\b', content)

        for block in blocks:
            block = block.strip()
            if not block:
                continue
            self._parse_block(block)

        logger.debug(f"Parsed {len(self.residues)} residues from {self.prep_file_path.name}")
        return self.residues

    def _parse_block(self, block: str):
        """Parse a single residue block from a .prep file."""
        lines = block.split('\n')

        # Find the residue name line: "RES   INT  0" pattern
        resname = None
        atom_start = None

        for i, line in enumerate(lines):
            stripped = line.strip()
            # Skip empty lines and remark lines
            if not stripped:
                continue

            # Look for residue declaration: "ORE   INT  0" or similar
            match = re.match(r'^(\w{1,4})\s+INT\s+\d+', stripped)
            if match:
                resname = match.group(1)
                # Atoms start after the "CORRECT..." line and the charge line
                # Find them by looking for DUMM or numbered atom lines
                for j in range(i + 1, len(lines)):
                    atom_match = re.match(r'\s*\d+\s+\w+', lines[j].strip())
                    if atom_match:
                        atom_start = j
                        break
                break

        if not resname or atom_start is None:
            return

        # Parse atom lines
        atoms = []
        for i in range(atom_start, len(lines)):
            line = lines[i].strip()
            if not line or line.startswith('LOOP') or line.startswith('IMPROPER') or line.startswith('DONE'):
                break

            # Format: INDEX ATOM_NAME ATOM_TYPE STATUS P1 P2 P3 DIST ANGLE DIHEDRAL CHARGE
            match = re.match(
                r'\s*(\d+)\s+(\S+)\s+(\S+)\s+[MESB3]\s+'
                r'[-\d]+\s+[-\d]+\s+[-\d]+\s+'
                r'[\d.]+\s+[-\d.]+\s+[-\d.]+\s+'
                r'([-\d.]+)',
                line
            )
            if match:
                atom_name = match.group(2)
                atom_type = match.group(3)
                charge = float(match.group(4))

                # Skip dummy atoms
                if atom_name == 'DUMM' or atom_type == 'DU':
                    continue

                atoms.append({
                    "name": atom_name,
                    "type": atom_type,
                    "charge": charge,
                    "element_num": -1  # Not available in prep format
                })

        if atoms:
            self.residues[resname] = {
                "atoms": atoms,
                "connectivity": [],
                "positions": [],
                "connection_points": {"n_terminus": 0, "c_terminus": 0}
            }

    def get_residue_atoms(self, resname: str) -> List[Dict]:
        """Get atom list for a specific residue."""
        if resname not in self.residues:
            if not self.residues:
                self.parse()
            if resname not in self.residues:
                return []
        return self.residues[resname].get("atoms", [])

    def get_atom_type_and_charge(self, resname: str, atom_name: str) -> Optional[Tuple[str, float]]:
        """Get atom type and charge for a specific atom in a residue."""
        atoms = self.get_residue_atoms(resname)
        for atom in atoms:
            if atom["name"] == atom_name:
                return (atom["type"], atom["charge"])
        return None

    def get_connectivity_template(self, resname: str) -> List[Tuple[int, int]]:
        """Get connectivity template (not available from prep format)."""
        return []

    def get_connection_points(self, resname: str) -> Dict[str, int]:
        """Get connection points (not available from prep format)."""
        return {"n_terminus": 0, "c_terminus": 0}
