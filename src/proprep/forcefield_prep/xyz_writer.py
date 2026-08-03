"""
XYZ File Writer for QM Calculations

Generates QM-software independent XYZ files from model definitions.
XYZ files contain atom coordinates for geometry optimization and RESP fitting.

Format:
  Line 1: Number of atoms
  Line 2: Comment line (total charge + description)
  Lines 3+: Element X Y Z
"""

from pathlib import Path
from typing import List, Dict, Tuple, Optional
import numpy as np

from Bio.PDB import PDBParser
from rich.console import Console

from .model_builder import ModelBuilder, ModelResidue, CapType
from proprep.structure_prep.comprehensive_redox_detector import RedoxSite


# MCPB capping group atom definitions
# Based on MCPB gene_model_files.py write_ace, write_nme, write_gly functions

ACE_ATOMS = {
    # ACE: CH3-CO-
    # Keeps: C, O, CH3 (from CA), HH31, HH32, HH33 (from HA, CB, N)
    'C': 'C',
    'O': 'O',
    'CH3': 'C',  # Renamed from CA
    'HH31': 'H',
    'HH32': 'H',
    'HH33': 'H'
}

NME_ATOMS = {
    # NME: -NH-CH3
    # Keeps: N, H, CH3 (from CA), HH31, HH32, HH33 (from C, HA, CB)
    'N': 'N',
    'H': 'H',
    'CH3': 'C',  # Renamed from CA
    'HH31': 'H',
    'HH32': 'H',
    'HH33': 'H'
}

ANT_ATOMS = {
    # ANT: CH3NH3+ (protonated methylamine)
    # Keeps: N, H1, H2, H3, CH3, HH31, HH32, HH33
    'N': 'N',
    'H1': 'H',
    'H2': 'H',
    'H3': 'H',
    'CH3': 'C',
    'HH31': 'H',
    'HH32': 'H',
    'HH33': 'H'
}

ACT_ATOMS = {
    # ACT: CH3CO2- (acetate)
    # Keeps: C, O, OXT, CH3, HH31, HH32, HH33
    'C': 'C',
    'O': 'O',
    'OXT': 'O',
    'CH3': 'C',
    'HH31': 'H',
    'HH32': 'H',
    'HH33': 'H'
}


class XYZWriter:
    """
    Writes XYZ files for QM calculations.

    XYZ format is software-independent and contains:
    - Atom count
    - Total charge + description
    - Element coordinates (from PDB)
    """

    def __init__(self, pdb_file: str, console: Optional[Console] = None):
        """
        Initialize XYZ writer.

        Args:
            pdb_file: Path to PDB file
            console: Rich console for output
        """
        self.pdb_file = Path(pdb_file)
        self.console = console or Console()

        # Parse PDB structure
        parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure('structure', str(self.pdb_file))

        # Build residue map
        self.residue_map = {}  # (chain, resid) -> BioPython Residue
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    het_flag, resseq, icode = residue.get_id()
                    key = (chain.id, resseq)
                    self.residue_map[key] = residue

    def write_xyz(self, model: ModelBuilder, atom_type_assignments: Dict,
                  output_file: str, description: str = "QM model"):
        """
        Write XYZ file from model.

        Also creates an atom ID mapping file (.mapping) that maps XYZ indices to original PDB atom IDs.
        This is required for the Seminario method to correctly map fingerprint atom IDs to XYZ coordinates.

        Args:
            model: ModelBuilder with residues to include
            atom_type_assignments: Dict mapping coords -> {'charge': float, ...}
            output_file: Path to output XYZ file
            description: Description for comment line
        """
        xyz_atoms = []  # List of (element, x, y, z, charge)
        atom_id_mapping = []  # List of (xyz_index, pdb_atom_id, resid, resname, atom_name)
        total_charge = 0.0
        xyz_index = 0

        # Process each residue in model
        for model_residue in model.model_residues:
            residue = self.residue_map.get((model_residue.chain, model_residue.resid))
            if residue is None:
                continue

            if model_residue.cap_type != CapType.NONE:
                # Capping group - use special handling
                cap_atoms = self._get_capping_group_atoms(model_residue, atom_type_assignments)
                # Convert to (element, x, y, z, charge) for XYZ writing
                for atom_name, element, x, y, z, charge in cap_atoms:
                    xyz_atoms.append((element, x, y, z, charge))
                    total_charge += charge

                    # For capping groups, find matching atom in PDB to get serial number
                    # Capping groups reuse atom positions from original residue
                    for atom in residue.get_atoms():
                        atom_coords = tuple(atom.get_coord())
                        if abs(atom_coords[0] - x) < 0.01 and abs(atom_coords[1] - y) < 0.01 and abs(atom_coords[2] - z) < 0.01:
                            atom_id_mapping.append((xyz_index, atom.serial_number, model_residue.resid, model_residue.resname, atom_name))
                            break
                    else:
                        # Generated atom (H on caps) - use a placeholder ID
                        atom_id_mapping.append((xyz_index, -1, model_residue.resid, model_residue.resname, atom_name))

                    xyz_index += 1
            else:
                # CapType.NONE - apply sidechain-only transformation
                # Removes backbone, converts CA→CH3, generates 3 H atoms
                sidechain_atoms = self._extract_sidechain_atoms(model_residue, atom_type_assignments)
                # Convert to (element, x, y, z, charge) for XYZ writing
                for atom_name, element, x, y, z, charge in sidechain_atoms:
                    xyz_atoms.append((element, x, y, z, charge))
                    total_charge += charge

                    # Map XYZ index to original PDB atom ID
                    # For generated H atoms (H1, H2, H3), use placeholder
                    if atom_name in ['H1', 'H2', 'H3']:
                        atom_id_mapping.append((xyz_index, -1, model_residue.resid, model_residue.resname, atom_name))
                    else:
                        # Find original atom in PDB to get serial number
                        for atom in residue.get_atoms():
                            if atom.get_id() == atom_name or (atom_name == 'CH3' and atom.get_id() == 'CA'):
                                atom_id_mapping.append((xyz_index, atom.serial_number, model_residue.resid, model_residue.resname, atom_name))
                                break
                        else:
                            atom_id_mapping.append((xyz_index, -1, model_residue.resid, model_residue.resname, atom_name))

                    xyz_index += 1

        # Write XYZ file
        self._write_xyz_file(xyz_atoms, total_charge, description, output_file)

        # Write atom ID mapping file
        mapping_file = str(output_file).replace('.xyz', '.mapping')
        self._write_mapping_file(atom_id_mapping, mapping_file)

    def _get_capping_group_atoms(self, model_residue: ModelResidue,
                                 atom_type_assignments: Dict) -> List[Tuple[str, str, float, float, float, float]]:
        """
        Get atoms for a capping group with properly generated coordinates.

        Based on MCPB gene_model_files.py:130-230 (write_ace, write_nme, etc.)

        Args:
            model_residue: ModelResidue with cap_type set
            atom_type_assignments: Atom type assignments for charges

        Returns:
            List of (atom_name, element, x, y, z, charge) tuples
        """
        residue = self.residue_map.get((model_residue.chain, model_residue.resid))
        if residue is None:
            return []

        atoms = []
        cap_type = model_residue.cap_type

        # CH bond length for hydrogen generation
        CH_BOND_LENGTH = 1.09  # Angstroms

        if cap_type == CapType.ACE:
            # ACE: CH3-CO-
            # Keeps: C, O, CH3 (from CA), HH31-33 (from HA/CB/N)
            ca_coord = None
            c_coord = None
            o_coord = None

            # Get coordinates
            for atom in residue.get_atoms():
                if atom.get_id() == 'CA':
                    ca_coord = atom.get_coord()
                elif atom.get_id() == 'C':
                    c_coord = atom.get_coord()
                elif atom.get_id() == 'O':
                    o_coord = atom.get_coord()

            if ca_coord is not None:
                # Add CH3 (renamed from CA)
                charge = self._get_charge_from_coords(tuple(ca_coord), atom_type_assignments)
                atoms.append(('CH3', 'C', ca_coord[0], ca_coord[1], ca_coord[2], charge))

                # Generate H atoms from HA, CB, N positions if they exist
                h_sources = ['HA', 'CB', 'N', 'HA2', 'HA3']
                h_names = ['HH31', 'HH32', 'HH33']
                h_count = 0
                for atom in residue.get_atoms():
                    if atom.get_id() in h_sources and h_count < 3:
                        # Generate H at correct bond length
                        atom_coord = atom.get_coord()
                        vec = atom_coord - ca_coord
                        dist = (vec[0]**2 + vec[1]**2 + vec[2]**2)**0.5
                        if dist > 0.01:
                            h_coord = ca_coord + (CH_BOND_LENGTH / dist) * vec
                            atoms.append((h_names[h_count], 'H', h_coord[0], h_coord[1], h_coord[2], 0.0))
                            h_count += 1

            # Add C and O
            if c_coord is not None:
                charge = self._get_charge_from_coords(tuple(c_coord), atom_type_assignments)
                atoms.append(('C', 'C', c_coord[0], c_coord[1], c_coord[2], charge))

            if o_coord is not None:
                charge = self._get_charge_from_coords(tuple(o_coord), atom_type_assignments)
                atoms.append(('O', 'O', o_coord[0], o_coord[1], o_coord[2], charge))

        elif cap_type == CapType.NME:
            # NME: -NH-CH3
            # Keeps: N, H, CH3 (from CA), HH31-33 (from C/HA/CB)
            ca_coord = None
            n_coord = None
            h_coord = None

            for atom in residue.get_atoms():
                if atom.get_id() == 'CA':
                    ca_coord = atom.get_coord()
                elif atom.get_id() == 'N':
                    n_coord = atom.get_coord()
                elif atom.get_id() in ['H', 'HN']:
                    h_coord = atom.get_coord()

            # Add N
            if n_coord is not None:
                charge = self._get_charge_from_coords(tuple(n_coord), atom_type_assignments)
                atoms.append(('N', 'N', n_coord[0], n_coord[1], n_coord[2], charge))

            # Add H on nitrogen
            if h_coord is not None:
                charge = self._get_charge_from_coords(tuple(h_coord), atom_type_assignments)
                atoms.append(('H', 'H', h_coord[0], h_coord[1], h_coord[2], charge))

            # Add CH3 and hydrogens
            if ca_coord is not None:
                charge = self._get_charge_from_coords(tuple(ca_coord), atom_type_assignments)
                atoms.append(('CH3', 'C', ca_coord[0], ca_coord[1], ca_coord[2], charge))

                # Generate H atoms on methyl
                h_sources = ['C', 'HA', 'CB', 'HA2', 'HA3']
                h_names = ['HH31', 'HH32', 'HH33']
                h_count = 0
                for atom in residue.get_atoms():
                    if atom.get_id() in h_sources and h_count < 3:
                        atom_coord = atom.get_coord()
                        vec = atom_coord - ca_coord
                        dist = (vec[0]**2 + vec[1]**2 + vec[2]**2)**0.5
                        if dist > 0.01:
                            h_coord = ca_coord + (CH_BOND_LENGTH / dist) * vec
                            atoms.append((h_names[h_count], 'H', h_coord[0], h_coord[1], h_coord[2], 0.0))
                            h_count += 1

        elif cap_type == CapType.GLY:
            # GLY: Keep all backbone atoms
            # Glycine should have: N, H, CA, HA2, HA3, C, O
            # Handle case where source residue is not GLY (e.g., has HA instead of HA2/HA3)
            # Special case: PRO has no H on N (secondary amine in ring) - must generate it
            ca_coord = None
            n_coord = None
            c_coord = None
            has_n_hydrogen = False
            ha_atoms = []  # Collect H-alpha atoms
            NH_BOND_LENGTH = 1.01  # Angstroms

            for atom in residue.get_atoms():
                atom_id = atom.get_id()
                if atom_id == 'CA':
                    ca_coord = atom.get_coord()
                    element = atom.element
                    charge = self._get_charge_from_coords(tuple(ca_coord), atom_type_assignments)
                    atoms.append(('CA', element, ca_coord[0], ca_coord[1], ca_coord[2], charge))
                elif atom_id == 'N':
                    n_coord = atom.get_coord()
                    element = atom.element
                    coords = atom.get_coord()
                    charge = self._get_charge_from_coords(tuple(coords), atom_type_assignments)
                    atoms.append(('N', element, coords[0], coords[1], coords[2], charge))
                elif atom_id == 'C':
                    c_coord = atom.get_coord()
                    element = atom.element
                    coords = atom.get_coord()
                    charge = self._get_charge_from_coords(tuple(coords), atom_type_assignments)
                    atoms.append(('C', element, coords[0], coords[1], coords[2], charge))
                elif atom_id == 'O':
                    element = atom.element
                    coords = atom.get_coord()
                    charge = self._get_charge_from_coords(tuple(coords), atom_type_assignments)
                    atoms.append(('O', element, coords[0], coords[1], coords[2], charge))
                elif atom_id in ['H', 'HN']:
                    has_n_hydrogen = True
                    element = atom.element
                    coords = atom.get_coord()
                    charge = self._get_charge_from_coords(tuple(coords), atom_type_assignments)
                    atoms.append(('H', element, coords[0], coords[1], coords[2], charge))
                elif atom_id in ['HA', 'HA2', 'HA3']:
                    # Collect H-alpha atoms
                    ha_atoms.append((atom_id, atom.element, atom.get_coord()))

            # Check if we need to generate N-H hydrogen (e.g., PRO→GLY conversion)
            if n_coord is not None and not has_n_hydrogen:
                # Generate H on nitrogen with proper trigonal planar geometry
                # Amide N is sp² hybridized - H should be in peptide plane (C(i-1)-N-CA)
                # at ~120° bond angle, where C(i-1) is the carbonyl C from PRECEDING residue
                import numpy as np

                # Get the preceding residue's carbonyl carbon
                prev_residue = self.residue_map.get((model_residue.chain, model_residue.resid - 1))
                prev_c_coord = None
                if prev_residue:
                    for atom in prev_residue.get_atoms():
                        if atom.get_id() == 'C':
                            prev_c_coord = atom.get_coord()
                            break

                if ca_coord is not None and prev_c_coord is not None:
                    # Use C from previous residue and CA from current to define peptide plane
                    # H should be in this plane, making ~120° with both C-N and CA-N

                    # Vectors from N to its substituents
                    v_nc = prev_c_coord - n_coord  # N to carbonyl C (residue i-1)
                    v_nca = ca_coord - n_coord      # N to CA (residue i)

                    # Normalize
                    v_nc = v_nc / np.linalg.norm(v_nc)
                    v_nca = v_nca / np.linalg.norm(v_nca)

                    # Bisector of C-N-CA angle (opposite direction is where H should go)
                    bisector = v_nc + v_nca
                    bisector = bisector / np.linalg.norm(bisector)

                    # H goes opposite to bisector for trigonal planar geometry
                    h_direction = -bisector

                    # Place H at N-H bond length
                    h_coord = n_coord + NH_BOND_LENGTH * h_direction
                    atoms.append(('H', 'H', h_coord[0], h_coord[1], h_coord[2], 0.0))

                elif ca_coord is not None:
                    # Fallback: only CA available, place H opposite to CA
                    vec = n_coord - ca_coord
                    dist = (vec[0]**2 + vec[1]**2 + vec[2]**2)**0.5
                    if dist > 0.01:
                        h_coord = n_coord + (NH_BOND_LENGTH / dist) * vec
                        atoms.append(('H', 'H', h_coord[0], h_coord[1], h_coord[2], 0.0))

            # Process H-alpha atoms: ensure we have HA2 and HA3
            if ha_atoms:
                # If we have HA2 and/or HA3, use them as-is
                has_ha2 = any(name == 'HA2' for name, _, _ in ha_atoms)
                has_ha3 = any(name == 'HA3' for name, _, _ in ha_atoms)

                for atom_id, element, coords in ha_atoms:
                    charge = self._get_charge_from_coords(tuple(coords), atom_type_assignments)

                    if atom_id == 'HA':
                        # Rename HA to HA2 for glycine if HA2 doesn't exist
                        if not has_ha2:
                            atoms.append(('HA2', element, coords[0], coords[1], coords[2], charge))
                            has_ha2 = True
                        elif not has_ha3:
                            atoms.append(('HA3', element, coords[0], coords[1], coords[2], charge))
                            has_ha3 = True
                    else:
                        # HA2 or HA3 - use as-is
                        atoms.append((atom_id, element, coords[0], coords[1], coords[2], charge))

                # If we only have one H-alpha but GLY needs two, generate the second one
                # by placing a hydrogen where CB would be (correct tetrahedral geometry)
                if ca_coord is not None and (not has_ha2 or not has_ha3):
                    # Find CB - this is where the second H should go for glycine
                    for atom in residue.get_atoms():
                        if atom.get_id() == 'CB':
                            cb_coord = atom.get_coord()
                            vec = cb_coord - ca_coord
                            dist = (vec[0]**2 + vec[1]**2 + vec[2]**2)**0.5
                            if dist > 0.01:
                                # Place hydrogen at correct C-H bond length along CA→CB vector
                                h_coord = ca_coord + (CH_BOND_LENGTH / dist) * vec
                                missing_name = 'HA3' if has_ha2 else 'HA2'
                                atoms.append((missing_name, 'H', h_coord[0], h_coord[1], h_coord[2], 0.0))
                                break

        elif cap_type == CapType.KNH:
            # KNH: Keep N, H, and sidechain
            # Remove: C, O, HA
            # Transform: CA → CH3
            # Generate: 2 H atoms in directions of C and HA
            # Based on MCPB write_sc_knh() at gene_model_files.py:1010-1084

            ca_coord = None
            h_sources = {}  # For H generation

            # First pass: find CA and collect H source positions
            for atom in residue.get_atoms():
                atom_id = atom.get_id()
                if atom_id == 'CA':
                    ca_coord = atom.get_coord()
                elif atom_id in ['C', 'HA']:
                    h_sources[atom_id] = atom.get_coord()

            # Second pass: add atoms
            for atom in residue.get_atoms():
                atom_id = atom.get_id()
                coords = atom.get_coord()
                element = atom.element
                charge = self._get_charge_from_coords(tuple(coords), atom_type_assignments)

                if atom_id == 'CA':
                    # Rename CA to CH3
                    atoms.append(('CH3', 'C', coords[0], coords[1], coords[2], charge))
                elif atom_id in ['N', 'H', 'HN']:
                    # Keep N-H group
                    if atom_id == 'HN':
                        atom_id = 'H'  # Normalize to H
                    atoms.append((atom_id, element, coords[0], coords[1], coords[2], charge))
                elif atom_id not in ['C', 'O', 'OXT', 'HA']:
                    # Keep sidechain (CB and beyond)
                    atoms.append((atom_id, element, coords[0], coords[1], coords[2], charge))

            # Generate 2 H atoms in directions of C and HA
            if ca_coord is not None:
                h_names = ['H1', 'H2']
                h_count = 0
                for source_name in ['C', 'HA']:
                    if source_name in h_sources and h_count < 2:
                        source_coord = h_sources[source_name]
                        vec = source_coord - ca_coord
                        dist = (vec[0]**2 + vec[1]**2 + vec[2]**2)**0.5
                        if dist > 0.01:
                            h_coord = ca_coord + (CH_BOND_LENGTH / dist) * vec
                            atoms.append((h_names[h_count], 'H', h_coord[0], h_coord[1], h_coord[2], 0.0))
                            h_count += 1

        elif cap_type == CapType.KCO:
            # KCO: Keep C, O, and sidechain
            # Remove: N, H, HA
            # Transform: CA → CH3
            # Generate: 2 H atoms in directions of N and HA
            # Based on MCPB write_sc_kco() at gene_model_files.py:1092-1160

            ca_coord = None
            h_sources = {}  # For H generation

            # First pass: find CA and collect H source positions
            for atom in residue.get_atoms():
                atom_id = atom.get_id()
                if atom_id == 'CA':
                    ca_coord = atom.get_coord()
                elif atom_id in ['N', 'HA']:
                    h_sources[atom_id] = atom.get_coord()

            # Second pass: add atoms
            for atom in residue.get_atoms():
                atom_id = atom.get_id()
                coords = atom.get_coord()
                element = atom.element
                charge = self._get_charge_from_coords(tuple(coords), atom_type_assignments)

                if atom_id == 'CA':
                    # Rename CA to CH3
                    atoms.append(('CH3', 'C', coords[0], coords[1], coords[2], charge))
                elif atom_id in ['C', 'O']:
                    # Keep C=O group
                    atoms.append((atom_id, element, coords[0], coords[1], coords[2], charge))
                elif atom_id not in ['N', 'H', 'HN', 'HA', 'H1', 'H2', 'H3', 'HN1', 'HN2', 'HN3']:
                    # Keep sidechain (CB and beyond)
                    atoms.append((atom_id, element, coords[0], coords[1], coords[2], charge))

            # Generate 2 H atoms in directions of N and HA
            if ca_coord is not None:
                h_names = ['H1', 'H2']
                h_count = 0
                for source_name in ['N', 'HA']:
                    if source_name in h_sources and h_count < 2:
                        source_coord = h_sources[source_name]
                        vec = source_coord - ca_coord
                        dist = (vec[0]**2 + vec[1]**2 + vec[2]**2)**0.5
                        if dist > 0.01:
                            h_coord = ca_coord + (CH_BOND_LENGTH / dist) * vec
                            atoms.append((h_names[h_count], 'H', h_coord[0], h_coord[1], h_coord[2], 0.0))
                            h_count += 1

        return atoms

    def _extract_sidechain_atoms(self, model_residue: ModelResidue,
                                 atom_type_assignments: Dict) -> List[Tuple[str, str, float, float, float, float]]:
        """
        Extract sidechain-only atoms with CA→CH3 transformation.

        Based on MCPB write_sc() at gene_model_files.py:935-1002

        Transformation for sidechain-coordinating residues:
        1. Remove: N, C, O, H, HN, OXT, HA, H1, H2, H3, HN1, HN2, HN3
        2. Keep: CB and all sidechain atoms (beyond CB)
        3. Transform: CA → CH3 (rename, keep coordinates)
        4. Generate: 3 H atoms in directions of original N, C, HA

        Args:
            model_residue: ModelResidue with CapType.NONE
            atom_type_assignments: Atom type assignments for charges

        Returns:
            List of (atom_name, element, x, y, z, charge) tuples
        """
        residue = self.residue_map.get((model_residue.chain, model_residue.resid))
        if residue is None:
            return []

        atoms = []
        CH_BOND_LENGTH = 1.09  # Angstroms

        # Backbone atoms to remove
        BACKBONE_ATOMS = {'N', 'C', 'O', 'H', 'HN', 'OXT', 'HA', 'H1', 'H2', 'H3',
                          'HN1', 'HN2', 'HN3'}

        # Find CA coordinates (needed for H generation)
        ca_coord = None
        for atom in residue.get_atoms():
            if atom.get_id() == 'CA':
                ca_coord = atom.get_coord()
                break

        if ca_coord is None:
            # No CA atom - just keep all atoms (e.g., non-standard residue)
            for atom in residue.get_atoms():
                element = atom.element
                coords = atom.get_coord()
                charge = self._get_charge_from_coords(tuple(coords), atom_type_assignments)
                atoms.append((atom.get_id(), element, coords[0], coords[1], coords[2], charge))
            return atoms

        # Collect coordinates for H generation (N, C, HA positions)
        h_sources = {}  # atom_name -> coordinates
        for atom in residue.get_atoms():
            if atom.get_id() in ['N', 'C', 'HA']:
                h_sources[atom.get_id()] = atom.get_coord()

        # Add CA as CH3
        charge = self._get_charge_from_coords(tuple(ca_coord), atom_type_assignments)
        atoms.append(('CH3', 'C', ca_coord[0], ca_coord[1], ca_coord[2], charge))

        # Generate 3 H atoms in directions of N, C, HA
        h_names = ['H1', 'H2', 'H3']
        h_count = 0
        for source_name in ['N', 'C', 'HA']:
            if source_name in h_sources and h_count < 3:
                source_coord = h_sources[source_name]
                vec = source_coord - ca_coord
                dist = (vec[0]**2 + vec[1]**2 + vec[2]**2)**0.5
                if dist > 0.01:
                    # Place H at correct C-H bond length along CA→source vector
                    h_coord = ca_coord + (CH_BOND_LENGTH / dist) * vec
                    atoms.append((h_names[h_count], 'H', h_coord[0], h_coord[1], h_coord[2], 0.0))
                    h_count += 1

        # Add all sidechain atoms (not in backbone set, not CA)
        for atom in residue.get_atoms():
            atom_name = atom.get_id()
            if atom_name not in BACKBONE_ATOMS and atom_name != 'CA':
                element = atom.element
                coords = atom.get_coord()
                charge = self._get_charge_from_coords(tuple(coords), atom_type_assignments)
                atoms.append((atom_name, element, coords[0], coords[1], coords[2], charge))

        return atoms

    def _get_charge_from_coords(self, coords: Tuple[float, float, float],
                                atom_type_assignments: Dict) -> float:
        """
        Get charge for atom from assignments based on coordinates.

        Args:
            coords: Atom coordinates as tuple
            atom_type_assignments: Dict mapping coords -> {'charge': float, ...}

        Returns:
            Partial charge (0.0 if not found)
        """
        if coords in atom_type_assignments:
            return atom_type_assignments[coords].get('charge', 0.0)
        return 0.0

    def _write_xyz_file(self, atoms: List[Tuple[str, float, float, float, float]],
                       total_charge: float, description: str, output_file: str):
        """
        Write atoms to XYZ file.

        Args:
            atoms: List of (element, x, y, z, charge) tuples
            total_charge: Total charge of system
            description: Description for comment line
            output_file: Path to output file
        """
        with open(output_file, 'w') as f:
            # Line 1: Number of atoms
            f.write(f"{len(atoms)}\n")

            # Line 2: Comment (total charge + description)
            f.write(f"Total charge: {total_charge:.2f} | {description}\n")

            # Lines 3+: Element X Y Z
            for element, x, y, z, charge in atoms:
                f.write(f"{element:2s}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")

        self.console.print(f"[green]✓ Wrote XYZ file: {Path(output_file).name} ({len(atoms)} atoms, charge={total_charge:.2f})[/green]")

    def _write_mapping_file(self, atom_id_mapping: List[Tuple[int, int, int, str, str]], output_file: str):
        """
        Write atom ID mapping file.

        Maps XYZ file indices (0-based) to original PDB atom IDs (serial numbers).
        This is required for Seminario method to map fingerprint atom IDs to XYZ coordinates.

        Args:
            atom_id_mapping: List of (xyz_index, pdb_atom_id, resid, resname, atom_name) tuples
            output_file: Path to output mapping file
        """
        with open(output_file, 'w') as f:
            f.write("# XYZ_index  PDB_atom_id  ResID  ResName  AtomName\n")
            for xyz_idx, pdb_id, resid, resname, atom_name in atom_id_mapping:
                f.write(f"{xyz_idx:6d}  {pdb_id:6d}  {resid:5d}  {resname:5s}  {atom_name:5s}\n")

        self.console.print(f"[green]✓ Wrote atom ID mapping: {Path(output_file).name}[/green]")

    def write_xyz_from_residues(self, residues: List[ModelResidue],
                                atom_type_assignments: Dict,
                                output_file: str, description: str = "QM model"):
        """
        Simplified interface: write XYZ directly from residue list.

        Args:
            residues: List of ModelResidue objects
            atom_type_assignments: Dict mapping coords -> {'charge': float, ...}
            output_file: Path to output XYZ file
            description: Description for comment line
        """
        self.console.print(f"\n[grey50]Writing atoms to XYZ file: {Path(output_file).name}[/grey50]")

        xyz_atoms = []
        total_charge = 0.0

        for model_residue in residues:
            residue = self.residue_map.get((model_residue.chain, model_residue.resid))
            if residue is None:
                self.console.print(f"  [yellow]Warning: Residue {model_residue.chain}:{model_residue.resid} not found in PDB[/yellow]")
                continue

            atom_count = 0
            for atom in residue.get_atoms():
                element = atom.element
                coords = atom.get_coord()
                x, y, z = coords

                # Get charge from atom type assignments
                atom_coords_tuple = tuple(coords)
                charge = 0.0
                if atom_coords_tuple in atom_type_assignments:
                    charge = atom_type_assignments[atom_coords_tuple].get('charge', 0.0)

                xyz_atoms.append((element, x, y, z, charge))
                total_charge += charge
                atom_count += 1

            cap_str = f" [{model_residue.cap_type.value}]" if model_residue.cap_type.value != "NONE" else ""
            self.console.print(f"  • {model_residue.chain}:{model_residue.resid} {model_residue.resname:3s}{cap_str} - {atom_count} atoms")

        self._write_xyz_file(xyz_atoms, total_charge, description, output_file)
