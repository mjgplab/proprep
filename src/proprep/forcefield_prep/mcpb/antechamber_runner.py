"""
Antechamber Runner

Automates antechamber and parmchk2 workflow for ligand parameterization.
Extracts ligand atoms to PDB, runs GAFF typing and AM1-BCC charge calculation,
generates frcmod for missing parameters, and parses results.

© 2024 ProPrep Developer. All rights reserved.
"""

from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
from pathlib import Path
import subprocess
import logging
import re
import tempfile
import shutil

# Import RedoxSite types
import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from structure_prep.comprehensive_redox_detector import RedoxSiteAtom
from proprep.utils.prompts import prompt_with_context


@dataclass
class LigandParameters:
    """Results from antechamber parameterization"""
    atom_types: Dict[Tuple[float, float, float], str]      # coords -> GAFF atom type
    elements: Dict[Tuple[float, float, float], str]        # coords -> element symbol
    charges: Dict[Tuple[float, float, float], float]       # coords -> AM1-BCC charge
    frcmod_file: Optional[Path]                            # Path to generated frcmod
    mol2_file: Optional[Path]                              # Path to generated mol2
    success: bool                                          # Whether parameterization succeeded
    error_message: str = ""                                # Error details if failed


class AntechamberRunner:
    """
    Automate antechamber + parmchk2 workflow for ligand groups.

    Workflow:
    1. Write ligand atoms to temporary PDB file
    2. Run antechamber: PDB -> MOL2 with GAFF types and AM1-BCC charges
    3. Run parmchk2: Generate frcmod for missing parameters
    4. Parse MOL2 to extract atom types, elements, and charges
    5. Return LigandParameters object
    """

    def __init__(self, gaff_version: str = 'gaff2', logger=None, processor=None):
        """
        Initialize antechamber runner.

        Args:
            gaff_version: GAFF version to use ('gaff' or 'gaff2')
            logger: Optional logger instance
        """
        self.gaff_version = gaff_version
        self.logger = logger or logging.getLogger(__name__)
        self.processor = processor

    def parameterize_ligand_group(self,
                                  atoms: List[RedoxSiteAtom],
                                  group_id: str,
                                  formal_charge: int,
                                  output_dir: Path) -> LigandParameters:
        """
        Run complete antechamber workflow for a ligand group.

        Args:
            atoms: List of RedoxSiteAtom objects in this ligand group
            group_id: Unique identifier for this group (e.g., "ligand_1")
            formal_charge: Net formal charge for antechamber
            output_dir: Directory to save output files

        Returns:
            LigandParameters object with results or error information
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info(
            f"Parameterizing {group_id}: {len(atoms)} atoms, charge={formal_charge:+d}"
        )

        try:
            # Step 1: Write PDB file
            pdb_file = output_dir / f"{group_id}.pdb"
            self._write_ligand_pdb(atoms, pdb_file)
            self.logger.debug(f"Wrote ligand PDB: {pdb_file}")

            # Step 2: Run antechamber
            mol2_file = output_dir / f"{group_id}.mol2"
            success = self._run_antechamber(
                input_pdb=pdb_file,
                output_mol2=mol2_file,
                formal_charge=formal_charge
            )

            if not success:
                return LigandParameters(
                    atom_types={},
                    elements={},
                    charges={},
                    frcmod_file=None,
                    mol2_file=None,
                    success=False,
                    error_message="Antechamber failed (see logs)"
                )

            # Step 3: Run parmchk2
            frcmod_file = output_dir / f"{group_id}.frcmod"
            success = self._run_parmchk2(
                mol2_file=mol2_file,
                frcmod_file=frcmod_file
            )

            if not success:
                self.logger.warning(f"Parmchk2 failed for {group_id}, but continuing...")

            # Step 4: Parse MOL2 file
            atom_types, elements, charges = self._parse_mol2(mol2_file, atoms)

            self.logger.info(
                f"✓ Successfully parameterized {group_id}: "
                f"{len(atom_types)} atoms assigned GAFF types"
            )

            return LigandParameters(
                atom_types=atom_types,
                elements=elements,
                charges=charges,
                frcmod_file=frcmod_file if frcmod_file.exists() else None,
                mol2_file=mol2_file,
                success=True
            )

        except Exception as e:
            self.logger.error(f"Error parameterizing {group_id}: {e}")
            return LigandParameters(
                atom_types={},
                elements={},
                charges={},
                frcmod_file=None,
                mol2_file=None,
                success=False,
                error_message=str(e)
            )

    def _write_ligand_pdb(self, atoms: List[RedoxSiteAtom], output_file: Path):
        """
        Write ligand atoms to PDB file.

        Args:
            atoms: List of RedoxSiteAtom objects
            output_file: Path to output PDB file
        """
        with open(output_file, 'w') as f:
            for i, atom in enumerate(atoms, start=1):
                # PDB ATOM record format
                # ATOM   serial name resName chainID resSeq   x      y      z     occ   bfac element
                line = (
                    f"ATOM  {i:5d}  {atom.atom_name:<4s}{atom.resname:>3s} "
                    f"{atom.chain:1s}{atom.resid:4d}    "
                    f"{atom.coords[0]:8.3f}{atom.coords[1]:8.3f}{atom.coords[2]:8.3f}"
                    f"{atom.occupancy or 1.0:6.2f}{atom.bfactor or 0.0:6.2f}          "
                    f"{atom.element:>2s}\n"
                )
                f.write(line)
            f.write("END\n")

    def _run_antechamber(self, input_pdb: Path, output_mol2: Path,
                        formal_charge: int) -> bool:
        """
        Run antechamber to generate MOL2 file with GAFF types and AM1-BCC charges.

        Args:
            input_pdb: Input PDB file
            output_mol2: Output MOL2 file
            formal_charge: Net formal charge

        Returns:
            True if successful, False otherwise
        """
        cmd = [
            'antechamber',
            '-i', str(input_pdb),
            '-fi', 'pdb',
            '-o', str(output_mol2),
            '-fo', 'mol2',
            '-c', 'bcc',  # AM1-BCC charges
            '-at', self.gaff_version,
            '-nc', str(formal_charge),
            '-pf', 'y',  # Remove intermediate files
            '-dr', 'no'  # Don't check bond types
        ]

        self.logger.debug(f"Running: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout
            )

            if result.returncode != 0:
                self.logger.error(f"Antechamber failed with return code {result.returncode}")
                self.logger.error(f"STDOUT: {result.stdout}")
                self.logger.error(f"STDERR: {result.stderr}")
                return False

            if not output_mol2.exists():
                self.logger.error(f"Antechamber did not create output file: {output_mol2}")
                return False

            return True

        except subprocess.TimeoutExpired:
            self.logger.error("Antechamber timed out after 5 minutes")
            return False
        except FileNotFoundError:
            self.logger.error("Antechamber not found. Is AMBERHOME set correctly?")
            return False
        except Exception as e:
            self.logger.error(f"Error running antechamber: {e}")
            return False

    def _run_parmchk2(self, mol2_file: Path, frcmod_file: Path) -> bool:
        """
        Run parmchk2 to generate frcmod file for missing parameters.

        Args:
            mol2_file: Input MOL2 file from antechamber
            frcmod_file: Output frcmod file

        Returns:
            True if successful, False otherwise
        """
        cmd = [
            'parmchk2',
            '-i', str(mol2_file),
            '-f', 'mol2',
            '-o', str(frcmod_file),
            '-a', 'Y',  # Print all parameters
            '-s', self.gaff_version
        ]

        self.logger.debug(f"Running: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=60
            )

            if result.returncode != 0:
                self.logger.warning(f"Parmchk2 failed with return code {result.returncode}")
                self.logger.warning(f"STDERR: {result.stderr}")
                return False

            if not frcmod_file.exists():
                self.logger.warning(f"Parmchk2 did not create output file: {frcmod_file}")
                return False

            return True

        except subprocess.TimeoutExpired:
            self.logger.warning("Parmchk2 timed out after 60 seconds")
            return False
        except FileNotFoundError:
            self.logger.warning("Parmchk2 not found. Is AMBERHOME set correctly?")
            return False
        except Exception as e:
            self.logger.warning(f"Error running parmchk2: {e}")
            return False

    def _parse_mol2(self, mol2_file: Path,
                   original_atoms: List[RedoxSiteAtom]) -> Tuple[Dict, Dict, Dict]:
        """
        Parse MOL2 file to extract atom types, elements, and charges.

        Matches atoms by coordinates to original RedoxSiteAtom objects.

        Args:
            mol2_file: Path to MOL2 file
            original_atoms: Original RedoxSiteAtom objects (for coordinate matching)

        Returns:
            (atom_types, elements, charges) dicts mapping coords -> value
        """
        atom_types = {}
        elements = {}
        charges = {}

        # Read MOL2 file
        with open(mol2_file, 'r') as f:
            lines = f.readlines()

        # Find @<TRIPOS>ATOM section
        in_atom_section = False
        for line in lines:
            if line.startswith('@<TRIPOS>ATOM'):
                in_atom_section = True
                continue
            elif line.startswith('@<TRIPOS>'):
                in_atom_section = False
                continue

            if in_atom_section and line.strip():
                # MOL2 ATOM format:
                # atom_id atom_name x y z atom_type [subst_id] [subst_name] [charge] [status_bit]
                parts = line.split()
                if len(parts) < 9:
                    continue

                try:
                    atom_name = parts[1]
                    x = float(parts[2])
                    y = float(parts[3])
                    z = float(parts[4])
                    gaff_type = parts[5]
                    charge = float(parts[8])

                    # Match to original atom by coordinates (within tolerance)
                    coords = self._match_coords_to_original((x, y, z), original_atoms)

                    if coords:
                        # Extract element from GAFF type (first 1-2 chars)
                        element = self._extract_element_from_gaff(gaff_type)

                        atom_types[coords] = gaff_type
                        elements[coords] = element
                        charges[coords] = charge

                except (ValueError, IndexError) as e:
                    self.logger.warning(f"Failed to parse MOL2 line: {line.strip()} ({e})")
                    continue

        return atom_types, elements, charges

    def _match_coords_to_original(self, coords: Tuple[float, float, float],
                                  original_atoms: List[RedoxSiteAtom],
                                  tolerance: float = 0.01) -> Optional[Tuple[float, float, float]]:
        """
        Match MOL2 coordinates to original RedoxSiteAtom coordinates.

        Args:
            coords: Coordinates from MOL2 file (x, y, z)
            original_atoms: Original RedoxSiteAtom objects
            tolerance: Coordinate matching tolerance in Angstroms

        Returns:
            Original coordinates (from RedoxSiteAtom) or None if no match
        """
        for atom in original_atoms:
            distance = (
                (coords[0] - atom.coords[0])**2 +
                (coords[1] - atom.coords[1])**2 +
                (coords[2] - atom.coords[2])**2
            )**0.5

            if distance < tolerance:
                return atom.coords

        self.logger.warning(f"No matching original atom found for coords {coords}")
        return None

    def _extract_element_from_gaff(self, gaff_type: str) -> str:
        """
        Extract element symbol from GAFF atom type.

        Examples:
            'c3' -> 'C'
            'n4' -> 'N'
            'oh' -> 'O'
            's6' -> 'S'
            'br' -> 'Br'
            'cl' -> 'Cl'

        Args:
            gaff_type: GAFF atom type

        Returns:
            Element symbol (capitalized)
        """
        # GAFF types start with element symbol (1-2 lowercase letters)
        # followed by optional numbers

        # Special cases for 2-letter elements
        two_letter = gaff_type[:2].lower()
        if two_letter in ['br', 'cl']:
            return two_letter.capitalize()

        # Otherwise, first letter is element
        return gaff_type[0].upper()

    def parse_existing_files(self,
                            mol2_file: Path,
                            frcmod_file: Optional[Path],
                            atoms: List[RedoxSiteAtom],
                            group_id: str) -> LigandParameters:
        """
        Parse existing mol2 and frcmod files (e.g., from small molecule parameterizer).

        Args:
            mol2_file: Path to existing mol2 file
            frcmod_file: Path to existing frcmod file (optional)
            atoms: List of RedoxSiteAtom objects in this ligand group
            group_id: Group identifier

        Returns:
            LigandParameters with data from existing files
        """
        mol2_file = Path(mol2_file)
        if not mol2_file.exists():
            self.logger.error(f"MOL2 file not found: {mol2_file}")
            return LigandParameters(
                atom_types={},
                elements={},
                charges={},
                frcmod_file=None,
                mol2_file=None,
                success=False,
                error_message=f"MOL2 file not found: {mol2_file}"
            )

        if frcmod_file:
            frcmod_file = Path(frcmod_file)
            if not frcmod_file.exists():
                self.logger.warning(f"Frcmod file not found: {frcmod_file}, continuing without it")
                frcmod_file = None

        self.logger.info(f"Parsing existing files for {group_id}: {mol2_file.name}")

        try:
            # Parse MOL2 file
            atom_types, elements, charges = self._parse_mol2(mol2_file, atoms)

            self.logger.info(
                f"✓ Successfully parsed {group_id}: "
                f"{len(atom_types)} atoms from existing mol2"
            )

            return LigandParameters(
                atom_types=atom_types,
                elements=elements,
                charges=charges,
                frcmod_file=frcmod_file,
                mol2_file=mol2_file,
                success=True
            )

        except Exception as e:
            self.logger.error(f"Error parsing existing files for {group_id}: {e}")
            return LigandParameters(
                atom_types={},
                elements={},
                charges={},
                frcmod_file=None,
                mol2_file=None,
                success=False,
                error_message=str(e)
            )

    def manual_entry_fallback(self, atoms: List[RedoxSiteAtom],
                             group_id: str) -> LigandParameters:
        """
        Fallback for manual entry when antechamber fails.

        Prompts user to manually enter atom types and charges.

        Args:
            atoms: List of RedoxSiteAtom objects
            group_id: Group identifier

        Returns:
            LigandParameters with manually entered data
        """
        print(f"\n⚠ Antechamber failed for {group_id}. Manual entry required.")
        print(f"\nAtoms in {group_id}:")

        atom_types = {}
        elements = {}
        charges = {}

        for i, atom in enumerate(atoms, start=1):
            print(f"\n  {i}. {atom.element} {atom.atom_name} at {atom.coords}")

            # Prompt for atom type
            while True:
                atom_type = prompt_with_context(
                    self.processor, "GAFF atom type",
                    module="Antechamber Runner",
                    description=f"GAFF atom type for {atom.element} {atom.atom_name}",
                ).strip()
                if atom_type:
                    break
                print("     Please enter a valid atom type.")

            # Prompt for charge
            while True:
                try:
                    charge = float(prompt_with_context(
                        self.processor, "Partial charge",
                        module="Antechamber Runner",
                        description=f"Partial charge for {atom.element} {atom.atom_name}",
                    ).strip())
                    break
                except ValueError:
                    print("     Please enter a valid number.")

            atom_types[atom.coords] = atom_type
            elements[atom.coords] = atom.element
            charges[atom.coords] = charge

        return LigandParameters(
            atom_types=atom_types,
            elements=elements,
            charges=charges,
            frcmod_file=None,
            mol2_file=None,
            success=True,
            error_message="Manual entry (antechamber unavailable)"
        )
