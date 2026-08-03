"""
ESP Data Extractor

Extracts electrostatic potential (ESP) data from Gaussian output files.
Generates .esp files in format required by AmberTools resp program.

Based on MCPB's get_esp_from_gau() implementation.
"""

from pathlib import Path
from typing import List, Tuple, Optional
import logging

# Bohr to Angstrom conversion
B_TO_A = 0.529177249


class ESPDataExtractor:
    """
    Extract electrostatic potential (ESP) data from Gaussian output.

    Generates .esp file in format required by AmberTools resp program.
    """

    def __init__(self, logger=None):
        """
        Initialize ESP data extractor.

        Args:
            logger: Optional logger instance
        """
        self.logger = logger or logging.getLogger(__name__)

    def extract_esp_data(self,
                        gaussian_log: str,
                        output_file: str = "large_resp.esp",
                        output_dir: Optional[str] = None) -> str:
        """
        Extract ESP data from Gaussian log file.

        Searches for two sections:
        1. "Electrostatic Properties Using The SCF Density" (coordinates)
        2. "Electrostatic Properties (Atomic Units)" (ESP values)

        Args:
            gaussian_log: Path to Gaussian .log file with ESP calculation
            output_file: Name for output .esp file
            output_dir: Directory for output (defaults to gaussian_log directory)

        Returns:
            Path to generated .esp file

        Raises:
            ValueError: If ESP data not found in log file
            FileNotFoundError: If gaussian_log doesn't exist
        """
        log_path = Path(gaussian_log)
        if not log_path.exists():
            raise FileNotFoundError(f"Gaussian log file not found: {gaussian_log}")

        # Determine output path
        if output_dir is None:
            output_dir = log_path.parent
        output_path = Path(output_dir) / output_file

        self.logger.debug(f"Extracting ESP data from {log_path.name}")

        # Parse ESP sections
        atom_coords, esp_coords, esp_values = self._parse_esp_sections(str(log_path))

        # Write ESP file
        self._write_esp_file(output_path, atom_coords, esp_coords, esp_values)

        self.logger.debug(f"ESP file written: {output_path.name}")
        return str(output_path)

    def _parse_esp_sections(self, log_file: str) -> Tuple[List, List, List]:
        """
        Parse ESP sections from Gaussian log.

        Section 1: "Electrostatic Properties Using The SCF Density"
        - Contains coordinates of atomic centers and ESP fit centers
        - Format: "      Atomic Center   N is at   x.xxx   y.yyy   z.zzz"
        - Coordinates in columns 32-42, 42-52, 52-62 (in Bohr)

        Section 2: "Electrostatic Properties (Atomic Units)"
        - Contains ESP values at each point
        - Starts 6 lines after section header
        - Format: "      Atom    N  ... ESP=  value" or "      Fit     N  ... ESP=  value"

        Args:
            log_file: Path to Gaussian log file

        Returns:
            Tuple of (atom_coords, esp_coords, esp_values)
            - atom_coords: List of (x, y, z) for atomic centers (Angstrom)
            - esp_coords: List of (x, y, z) for ESP fit points (Angstrom)
            - esp_values: List of ESP values at fit points (atomic units)

        Raises:
            ValueError: If required sections not found or data mismatch
        """
        # ================================================================
        # Step 1: Find section boundaries
        # ================================================================
        # For optimization jobs, there may be multiple iterations.
        # We need the LAST coordinate section and its matching ESP section.
        coord_section_lines = []
        esp_section_lines = []

        with open(log_file, 'r') as f:
            for ln, line in enumerate(f, 1):
                if 'Electrostatic Properties Using The SCF Density' in line:
                    coord_section_lines.append(ln)
                elif 'Electrostatic Properties (Atomic Units)' in line:
                    esp_section_lines.append(ln)

        if not coord_section_lines:
            raise ValueError(
                "No 'Electrostatic Properties Using The SCF Density' section found. "
                "Ensure Gaussian job finished normally and used Pop(MK,ReadRadii) keyword."
            )

        if not esp_section_lines:
            raise ValueError(
                "No 'Electrostatic Properties (Atomic Units)' section found. "
                "Ensure Gaussian job finished normally."
            )

        # Find the LAST coordinate section
        last_coord_start = coord_section_lines[-1]

        # Find the FIRST ESP section AFTER the last coordinate section
        # This ensures we get the matching pair from the same iteration
        last_esp_header = None
        for esp_line in esp_section_lines:
            if esp_line > last_coord_start:
                last_esp_header = esp_line
                break

        if last_esp_header is None:
            # Fallback: use the last ESP section (shouldn't happen with valid output)
            last_esp_header = esp_section_lines[-1]

        esp_data_start = last_esp_header + 6  # ESP values start 6 lines after header

        self.logger.debug(f"  Found {len(coord_section_lines)} coordinate section(s), using last at line {last_coord_start}")
        self.logger.debug(f"  Found {len(esp_section_lines)} ESP section(s), using one at line {last_esp_header}")

        # ================================================================
        # Step 2: Extract coordinates (atomic centers and ESP fit points)
        # ================================================================
        # Only extract from the last coordinate section, bounded by the ESP header
        atom_coords = []
        esp_coords = []

        with open(log_file, 'r') as f:
            for ln, line in enumerate(f, 1):
                # Only scan between coord section start and ESP section header
                if last_coord_start < ln < last_esp_header:
                    if '      Atomic Center' in line:
                        # Format: "      Atomic Center   1 is at   x.xxx   y.yyy   z.zzz"
                        # Coords in columns 32-42, 42-52, 52-62 (Bohr)
                        x = float(line[32:42]) / B_TO_A
                        y = float(line[42:52]) / B_TO_A
                        z = float(line[52:62]) / B_TO_A
                        atom_coords.append((x, y, z))
                    elif '     ESP Fit Center' in line:
                        # Same format as atomic centers
                        # Note: **** overflow in index field doesn't affect coords (fixed columns)
                        x = float(line[32:42]) / B_TO_A
                        y = float(line[42:52]) / B_TO_A
                        z = float(line[52:62]) / B_TO_A
                        esp_coords.append((x, y, z))
                elif ln >= last_esp_header:
                    # Past the ESP header, stop scanning for coordinates
                    break

        # ================================================================
        # Step 3: Extract ESP values
        # ================================================================
        atom_esp_values = []
        fit_esp_values = []

        with open(log_file, 'r') as f:
            for ln, line in enumerate(f, 1):
                if ln >= esp_data_start:
                    # Format: "    1 Atom     -1.049950" or "    1 Fit      -0.123456"
                    # Also handles **** overflow: " **** Atom    -1.049950" or " **** Fit     -0.123456"
                    stripped = line.strip()
                    if not stripped:
                        continue

                    # Check for Atom/Fit lines (index can be number or ****)
                    # Skip header lines like "Atom Element Radius" (no leading number or ****)
                    if ' Atom ' in line and (stripped[0].isdigit() or stripped.startswith('*')):
                        parts = stripped.split()
                        esp_value = float(parts[-1])
                        atom_esp_values.append(esp_value)
                    elif ' Fit ' in line and (stripped[0].isdigit() or stripped.startswith('*')):
                        parts = stripped.split()
                        esp_value = float(parts[-1])
                        fit_esp_values.append(esp_value)

        # ================================================================
        # Step 4: Validation
        # ================================================================
        if len(atom_coords) != len(atom_esp_values):
            raise ValueError(
                f"Coordinate/ESP mismatch for atoms: "
                f"{len(atom_coords)} coords vs {len(atom_esp_values)} ESP values"
            )

        if len(esp_coords) != len(fit_esp_values):
            raise ValueError(
                f"Coordinate/ESP mismatch for fit points: "
                f"{len(esp_coords)} coords vs {len(fit_esp_values)} ESP values"
            )

        self.logger.debug(f"  Extracted {len(atom_coords)} atomic centers")
        self.logger.debug(f"  Extracted {len(esp_coords)} ESP fit points")

        return atom_coords, esp_coords, fit_esp_values

    def _write_esp_file(self,
                       output_path: Path,
                       atom_coords: List[Tuple[float, float, float]],
                       esp_coords: List[Tuple[float, float, float]],
                       esp_values: List[float]) -> None:
        """
        Write .esp file in RESP format.

        Format (based on MCPB pymsmt/mol/gauio.py lines 676-683):
        Line 1: natoms  nesppts  0
        Lines 2-(natoms+1): 16 spaces + atomic coords (scientific notation, Angstrom)
        Lines (natoms+2)-end: ESP value + ESP point coords (scientific notation)

        Args:
            output_path: Output file path
            atom_coords: List of (x, y, z) tuples for atoms (Angstrom)
            esp_coords: List of (x, y, z) tuples for ESP points (Angstrom)
            esp_values: List of ESP values at each fit point (atomic units)
        """
        with open(output_path, 'w') as f:
            # Header: natoms nesppts 0
            f.write(f"{len(atom_coords):5d}{len(esp_coords):5d}{0:5d}\n")

            # Atomic centers (no ESP values)
            for x, y, z in atom_coords:
                f.write(f"{' ':16s} {x:15.7E} {y:15.7E} {z:15.7E}\n")

            # ESP fit points with values
            for (x, y, z), esp_val in zip(esp_coords, esp_values):
                f.write(f"{esp_val:16.7E} {x:15.7E} {y:15.7E} {z:15.7E}\n")
