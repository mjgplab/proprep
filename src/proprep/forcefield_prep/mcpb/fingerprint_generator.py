"""
Fingerprint File Generator

Generates MCPB-compatible fingerprint files from atom type assignments.
"""

from typing import Dict, Tuple, Optional
import logging

from proprep.structure_prep.comprehensive_redox_detector import RedoxSite
from .atom_typer import AtomTypeAssignment


class FingerprintGenerator:
    """
    Generates MCPB-compatible fingerprint files.

    Fingerprint format:
        ResID-ResName-AtomName  AtomID  OriginalType -> RenamedType
        ...
        LINK MetalAtomID-MetalName LigandAtomID-LigandName
        ...
    """

    def __init__(self, logger=None):
        """
        Initialize fingerprint generator.

        Args:
            logger: Optional logger instance
        """
        self.logger = logger or logging.getLogger(__name__)

    def write_fingerprint(self,
                         site: RedoxSite,
                         type_assignments: Dict[Tuple, AtomTypeAssignment],
                         output_file: str,
                         include_links: bool = True):
        """
        Write fingerprint file from atom type assignments.

        Args:
            site: RedoxSite object
            type_assignments: Dictionary mapping coordinates to AtomTypeAssignment
            output_file: Path to output fingerprint file
            include_links: Include LINK records for metal-ligand bonds

        Raises:
            ValueError: If atom IDs cannot be determined
        """
        self.logger.info(f"Writing fingerprint file: {output_file}")

        # Build coordinate to atom ID mapping
        coord_to_id = self._build_atom_id_map(site)

        with open(output_file, 'w') as f:
            # Part 1: Write atom type transformation lines
            self._write_atom_lines(f, site, type_assignments, coord_to_id)

            # Part 2: Write LINK records for metal-ligand bonds
            if include_links:
                self._write_link_records(f, site, coord_to_id)

        self.logger.info(f"Fingerprint file written: {output_file}")

    def _build_atom_id_map(self, site: RedoxSite) -> Dict[Tuple, int]:
        """
        Build mapping from coordinates to PDB serial numbers.

        Uses the original PDB atom serial numbers (not sequential numbering).
        This matches MCPB convention where fingerprint atom IDs are PDB serial numbers.

        Args:
            site: RedoxSite object

        Returns:
            Dictionary mapping coordinates to PDB serial numbers
        """
        coord_to_id = {}
        for atom in site.atoms:
            # Use PDB serial number if available, otherwise fall back to sequential
            if hasattr(atom, 'properties') and 'serial_number' in atom.properties:
                atom_id = atom.properties['serial_number']
            else:
                # Fallback for atoms without serial numbers (e.g., generated atoms)
                # Use a sequential ID starting from 10000 to avoid conflicts
                atom_id = 10000 + len(coord_to_id)
                self.logger.warning(
                    f"Atom {atom.atom_name} in {atom.resname}:{atom.resid} "
                    f"has no serial number, using fallback ID {atom_id}"
                )

            coord_to_id[atom.coords] = atom_id

        return coord_to_id

    def _write_atom_lines(self,
                         f,
                         site: RedoxSite,
                         type_assignments: Dict[Tuple, AtomTypeAssignment],
                         coord_to_id: Dict[Tuple, int]):
        """
        Write atom type transformation lines to fingerprint file.

        Format: ResID-ResName-AtomName  AtomID  OriginalType -> RenamedType

        Args:
            f: File handle
            site: RedoxSite object
            type_assignments: Atom type assignments
            coord_to_id: Coordinate to atom ID mapping
        """
        # Process atoms in order they appear in site.atoms
        for atom in site.atoms:
            if atom.coords not in type_assignments:
                self.logger.warning(
                    f"Atom {atom.resname}:{atom.resid}-{atom.atom_name} "
                    f"not found in type assignments, skipping"
                )
                continue

            assignment = type_assignments[atom.coords]
            atom_id = coord_to_id[atom.coords]

            # Format line: "ResID-ResName-AtomName  AtomID  OrigType -> RenamedType"
            # Ensure atom types are padded to 2 characters
            # Support both dict and AtomTypeAssignment dataclass
            if hasattr(assignment, 'original_type'):
                # AtomTypeAssignment dataclass - use attribute access
                orig_type = self._pad_atom_type(assignment.original_type)
                renamed_type = self._pad_atom_type(assignment.renamed_type)
                resid = assignment.resid
                resname = assignment.resname
                atom_name = assignment.atom_name
            else:
                # Dict - use subscript access
                orig_type = self._pad_atom_type(assignment['original_type'])
                renamed_type = self._pad_atom_type(assignment['renamed_type'])
                resid = assignment['resid']
                resname = assignment['resname']
                atom_name = assignment['atom_name']

            line = (
                f"{resid}-{resname}-{atom_name} "
                f"{atom_id} "
                f"{orig_type} -> {renamed_type}\n"
            )
            f.write(line)

    def _write_link_records(self,
                           f,
                           site: RedoxSite,
                           coord_to_id: Dict[Tuple, int]):
        """
        Write LINK records for metal-ligand bonds.

        Format: LINK AtomID1-AtomName1 AtomID2-AtomName2

        Args:
            f: File handle
            site: RedoxSite object
            coord_to_id: Coordinate to atom ID mapping
        """
        # Find all coordinate bonds (metal-ligand)
        link_count = 0

        for bond in site.bonds:
            # A restrained contact (nonbonded water held by an MD restraint) emits
            # no LINK record, so MCPB creates no metal-ligand bonded term for it.
            if bond.chemical_type == 'coordinate' and getattr(bond, 'treatment', 'bonded') == 'bonded':
                # Get atom IDs
                atom1_id = coord_to_id.get(bond.atom1_coords)
                atom2_id = coord_to_id.get(bond.atom2_coords)

                if atom1_id is None or atom2_id is None:
                    self.logger.warning(
                        f"LINK bond atoms not found in coordinate map, skipping"
                    )
                    continue

                # Get atom names from bond info
                atom1_name = bond.atom1_residue_info.get('atom_name', 'UNK')
                atom2_name = bond.atom2_residue_info.get('atom_name', 'UNK')

                # Write LINK record
                line = f"LINK {atom1_id}-{atom1_name} {atom2_id}-{atom2_name}\n"
                f.write(line)
                link_count += 1

        self.logger.info(f"  Wrote {link_count} LINK records")

    def _pad_atom_type(self, atom_type: str) -> str:
        """
        Pad atom type to 2 characters with space if needed.

        AMBER atom types are typically 2 characters, padded with space.

        Args:
            atom_type: Atom type string (1-2 characters)

        Returns:
            2-character padded atom type

        Examples:
            >>> _pad_atom_type('N')
            'N '
            >>> _pad_atom_type('CT')
            'CT'
        """
        if len(atom_type) == 1:
            return atom_type + ' '
        elif len(atom_type) == 2:
            return atom_type
        else:
            # Keep as is if longer (shouldn't happen in standard usage)
            return atom_type[:2]

    def write_large_fingerprint(self,
                               site: RedoxSite,
                               output_file: str):
        """
        Write large model fingerprint (simplified format for RESP fitting).

        Format (no atom types, no LINK records):
            ResID-ResName-AtomName
            ...

        Args:
            site: RedoxSite object
            output_file: Path to output fingerprint file
        """
        self.logger.info(f"Writing large model fingerprint: {output_file}")

        with open(output_file, 'w') as f:
            for atom in site.atoms:
                line = f"{atom.resid}-{atom.resname}-{atom.atom_name}\n"
                f.write(line)

        self.logger.info(f"Large model fingerprint written: {output_file}")

    def validate_fingerprint(self,
                            fingerprint_file: str) -> Dict:
        """
        Validate a fingerprint file and return statistics.

        Args:
            fingerprint_file: Path to fingerprint file

        Returns:
            Dictionary with validation results and statistics
        """
        self.logger.info(f"Validating fingerprint file: {fingerprint_file}")

        atom_lines = 0
        link_lines = 0
        renamed_atoms = 0
        errors = []

        with open(fingerprint_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue

                if line.startswith('LINK'):
                    # Validate LINK format
                    parts = line.split()
                    if len(parts) != 3:
                        errors.append(f"Line {line_num}: Invalid LINK format: {line}")
                    else:
                        link_lines += 1
                else:
                    # Validate atom line format
                    parts = line.split()
                    if len(parts) < 5 or parts[-2] != '->':
                        errors.append(f"Line {line_num}: Invalid atom line format: {line}")
                    else:
                        atom_lines += 1
                        # Check if renamed
                        orig_type = parts[-3]
                        renamed_type = parts[-1]
                        if orig_type != renamed_type:
                            renamed_atoms += 1

        result = {
            'valid': len(errors) == 0,
            'atom_lines': atom_lines,
            'link_lines': link_lines,
            'renamed_atoms': renamed_atoms,
            'errors': errors
        }

        if result['valid']:
            self.logger.info(
                f"  Validation passed: {atom_lines} atoms, {link_lines} LINKs, "
                f"{renamed_atoms} renamed"
            )
        else:
            self.logger.error(f"  Validation failed with {len(errors)} errors")
            for error in errors[:5]:  # Show first 5 errors
                self.logger.error(f"    {error}")

        return result

    def compare_fingerprints(self,
                            fingerprint1: str,
                            fingerprint2: str) -> Dict:
        """
        Compare two fingerprint files for differences.

        Useful for debugging or verifying regenerated fingerprints.

        Args:
            fingerprint1: Path to first fingerprint file
            fingerprint2: Path to second fingerprint file

        Returns:
            Dictionary with comparison results
        """
        # Read both files
        lines1 = []
        lines2 = []

        with open(fingerprint1, 'r') as f:
            lines1 = [line.strip() for line in f if line.strip()]

        with open(fingerprint2, 'r') as f:
            lines2 = [line.strip() for line in f if line.strip()]

        # Find differences
        differences = []
        max_len = max(len(lines1), len(lines2))

        for i in range(max_len):
            line1 = lines1[i] if i < len(lines1) else None
            line2 = lines2[i] if i < len(lines2) else None

            if line1 != line2:
                differences.append({
                    'line_number': i + 1,
                    'file1': line1,
                    'file2': line2
                })

        result = {
            'identical': len(differences) == 0,
            'lines_file1': len(lines1),
            'lines_file2': len(lines2),
            'differences': differences
        }

        if result['identical']:
            self.logger.info("Fingerprints are identical")
        else:
            self.logger.warning(f"Found {len(differences)} differences")

        return result
