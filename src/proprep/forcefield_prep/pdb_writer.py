"""
PDB File Writer for MCPB QM Models

Generates PDB files for small and large QM models following MCPB conventions.
PDB format preserves atom serial numbers, enabling direct matching with fingerprint files.

Format:
  ATOM records: Standard PDB format with serial numbers
  CONECT records: Metal-ligand coordination bonds
  END: Terminator

This replaces the XYZ+mapping file approach with direct PDB serial number matching.
"""

from pathlib import Path
from typing import List, Dict, Tuple, Optional, Set, TYPE_CHECKING
import logging

from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain, Residue, Atom, Select
from Bio.PDB.Atom import DisorderedAtom
from rich.console import Console

from .model_builder import ModelBuilder, ModelResidue, CapType
from proprep.structure_prep.comprehensive_redox_detector import RedoxSite

if TYPE_CHECKING:
    from .forcefield_data import ForceFieldData


class PDBWriter:
    """
    Writes PDB files for MCPB QM calculations.

    Features:
    - Preserves original PDB serial numbers for fingerprint matching
    - Handles capping groups (ACE, NME, GLY, KNH, KCO)
    - Writes CONECT records for metal-ligand bonds
    - Compatible with MCPB workflow

    Charge Assignment Strategy:
    ============================

    The charge assignment follows a hierarchical strategy to ensure charge conservation
    and proper treatment of different atom types:

    1. **Metal Ions**:
       - Use user-provided formal charge (e.g., Cu2+ = +2.0, Fe3+ = +3.0)
       - Stored in MetalConfig and written to atom_type_assignments

    2. **Capping Groups (ACE, NME)**:
       - Always use force field library charges (e.g., ff14SB, ff19SB)
       - These atoms are generated/renamed and not in atom_type_assignments
       - Example: ACE CH3 = -0.3662, HH31-33 = +0.1123, C = +0.5972, O = -0.5679

    3. **GLY Backbone Caps**:
       - Primary: Look up in atom_type_assignments (if GLY is part of metal site)
       - Fallback: Use force field charges (if GLY is a spacer in large model)
       - This dual approach handles both metal-coordinating and non-coordinating GLY

    4. **Sidechain Truncation (CA→CH3 + H1/H2/H3)**:
       - Applies to: Sidechain-only residues (CapType.NONE)
       - Generated H atoms inherit charges to preserve residue neutrality:
         * H1 (toward N): charge = charge(N) + charge(H)  [replaces N-H group]
         * H2 (toward C): charge = charge(C) + charge(O)  [replaces C=O group]
         * H3 (toward HA): charge = charge(HA)            [replaces HA]
       - Charge conservation: total charge of truncated residue equals original residue

    5. **KNH Capping (Keep N+H, remove C+O, CA→CH3 + H1/H2)**:
       - H1 (toward C): charge = charge(C) + charge(O)    [replaces C=O group]
       - H2 (toward HA): charge = charge(HA)              [replaces HA]

    6. **KCO Capping (Keep C+O, remove N+H, CA→CH3 + H1/H2)**:
       - H1 (toward N): charge = charge(N) + charge(H)    [replaces N-H group]
       - H2 (toward HA): charge = charge(HA)              [replaces HA]

    7. **Real Atoms (sidechains, backbone of metal-coordinating residues)**:
       - Use atom_type_assignments (populated during MCPB Step 1 atom typing)
       - Contains AMBER atom types and charges for metal site and coordinating residues

    Charge Conservation Principle:
    ==============================

    When truncating residues (removing backbone, keeping sidechain), generated H atoms
    on the CH3 group receive the sum of charges from the atoms they replace. This ensures
    the truncated residue has the same total charge as the original residue.

    Example - Neutral HIE residue:
        Original: N(-0.42) + H(+0.27) + CA(-0.06) + HA(+0.14) + sidechain(+0.04) + C(+0.60) + O(-0.57) = 0.00
        Truncated: CH3(-0.06) + sidechain(+0.04) + H1(-0.15) + H2(+0.03) + H3(+0.14) = 0.00

        Where:
        - H1 charge = N + H = -0.42 + 0.27 = -0.15
        - H2 charge = C + O = +0.60 - 0.57 = +0.03
        - H3 charge = HA = +0.14

    This approach maintains charge neutrality for neutral residues and preserves
    the correct total charge for charged residues (e.g., CYM thiolate ≈ -1.0).
    """

    def __init__(self, pdb_file: str, ff_data: Optional['ForceFieldData'] = None,
                 console: Optional[Console] = None):
        """
        Initialize PDB writer.

        Args:
            pdb_file: Path to original PDB file
            ff_data: ForceFieldData for charge/type lookups (unified FF data)
            console: Rich console for output
        """
        self.pdb_file = Path(pdb_file)
        self.ff_data = ff_data
        self.console = console or Console()
        self.logger = logging.getLogger(__name__)

        # Coordinate-keyed assignments, re-keyed to a rounded tuple. Rebuilt
        # when a different dict arrives. See _coord_key.
        self._assignment_index = None
        self._assignment_index_src = None

        # Parse PDB structure
        parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure('structure', str(self.pdb_file))

        # Build residue map: (chain, resid) -> BioPython Residue
        self.residue_map = {}
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    het_flag, resseq, icode = residue.get_id()
                    key = (chain.id, resseq)
                    self.residue_map[key] = residue

        # Build atom map: coords -> BioPython Atom (for serial number lookup)
        self.atom_map = {}
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    for atom in residue.get_atoms():
                        coords = tuple(round(x, 3) for x in atom.coord)
                        self.atom_map[coords] = atom

        # Track serial numbers for new atoms (caps)
        self.next_serial = 1
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    for atom in residue.get_atoms():
                        if atom.serial_number >= self.next_serial:
                            self.next_serial = atom.serial_number + 1

    def write_pdb(self, model: ModelBuilder, atom_type_assignments: Dict,
                  output_file: str, title: str = "QM model",
                  redox_site: Optional[RedoxSite] = None):
        """
        Write PDB file from model.

        Args:
            model: ModelBuilder with residues to include
            atom_type_assignments: Dict mapping coords -> {'charge': float, ...}
            output_file: Path to output PDB file
            title: Title for HEADER record
            redox_site: Source of metal-ligand coordinate bonds for CONECT records.
                Within-residue bonds are derived by distance heuristic.
        """
        atoms_data = []  # List of (serial, name, resname, chain, resid, x, y, z, element, charge)
        total_charge = 0.0

        # Process each residue in model
        for model_residue in model.model_residues:
            residue = self.residue_map.get((model_residue.chain, model_residue.resid))
            if residue is None:
                continue

            if model_residue.cap_type != CapType.NONE:
                # Capping group - determine correct residue name
                if model_residue.cap_type == CapType.ACE:
                    cap_resname = 'ACE'
                elif model_residue.cap_type == CapType.NME:
                    cap_resname = 'NME'
                elif model_residue.cap_type == CapType.GLY:
                    cap_resname = 'GLY'
                else:
                    # FULL, KNH, KCO - keep original name
                    cap_resname = model_residue.resname

                cap_atoms = self._get_capping_group_atoms(model_residue, atom_type_assignments)
                for atom_name, element, x, y, z, charge, serial in cap_atoms:
                    atoms_data.append((
                        serial, atom_name, cap_resname,  # Use cap residue name
                        model_residue.chain, model_residue.resid,
                        x, y, z, element, charge
                    ))
                    total_charge += charge
            else:
                # Sidechain-only (CapType.NONE) - keep original residue name
                sidechain_atoms = self._extract_sidechain_atoms(model_residue, atom_type_assignments)
                for atom_name, element, x, y, z, charge, serial in sidechain_atoms:
                    atoms_data.append((
                        serial, atom_name, model_residue.resname,  # Original name for sidechain
                        model_residue.chain, model_residue.resid,
                        x, y, z, element, charge
                    ))
                    total_charge += charge

        # Sort atoms by serial number before writing
        atoms_data.sort(key=lambda x: x[0])  # Sort by first element (serial)

        # Print charge breakdown table
        self._print_charge_breakdown(atoms_data, total_charge)

        # Derive CONECT records: within-residue bonds by distance, plus
        # metal-ligand coordinate bonds from the redox_site if provided.
        # RESP equivalence detection (resp_input_generator._find_ch_groups)
        # relies on CONECT to find CH2/CH3 hydrogens that should share
        # charge; without it those equivalences silently miss.
        conect_bonds = self._derive_conect_bonds(atoms_data, redox_site)

        self._write_pdb_file(atoms_data, total_charge, title, output_file, conect_bonds)

        self.console.print(f"[green]✓ Wrote PDB file: {Path(output_file).name} ({len(atoms_data)} atoms, charge={total_charge:.2f})[/green]")

    def _get_capping_group_atoms(self, model_residue: ModelResidue,
                                  atom_type_assignments: Dict) -> List[Tuple[str, str, float, float, float, float, int]]:
        """
        Get atoms for a capping group with properly generated coordinates.

        Applies charge assignment strategies #2, #3, #5, #6 depending on cap type:
        - ACE/NME: Force field charges (Strategy #2)
        - GLY: atom_type_assignments with force field fallback (Strategy #3)
        - KNH/KCO: Charge conservation for generated H atoms (Strategy #5/#6)

        Returns:
            List of (atom_name, element, x, y, z, charge, serial_number) tuples
        """
        residue = self.residue_map.get((model_residue.chain, model_residue.resid))
        if residue is None:
            return []

        atoms = []
        cap_type = model_residue.cap_type
        CH_BOND_LENGTH = 1.09  # Angstroms
        NH_BOND_LENGTH = 1.01  # Angstroms

        if cap_type == CapType.ACE:
            # ACE: CH3-CO-
            # Strategy #2: All atoms use force field charges (not in atom_type_assignments)
            # Keeps: C, O, CH3 (from CA), HH31-33 (from HA/CB/N)
            ca_coord = None
            c_coord = None
            o_coord = None

            for atom in residue.get_atoms():
                if atom.get_id() == 'CA':
                    ca_coord = atom.get_coord()
                elif atom.get_id() == 'C':
                    c_coord = atom.get_coord()
                elif atom.get_id() == 'O':
                    o_coord = atom.get_coord()

            if ca_coord is not None:
                # Add CH3 (renamed from CA) - use ACE charge from FF
                charge = self._get_charge_from_ff('ACE', 'CH3')
                serial = self._get_serial_from_coords(tuple(ca_coord))
                atoms.append(('CH3', 'C', ca_coord[0], ca_coord[1], ca_coord[2], charge, serial))

                # Generate H atoms - use ACE charges from FF
                h_sources = ['HA', 'CB', 'N', 'HA2', 'HA3']
                h_names = ['HH31', 'HH32', 'HH33']
                h_count = 0
                for atom in residue.get_atoms():
                    if atom.get_id() in h_sources and h_count < 3:
                        atom_coord = atom.get_coord()
                        vec = atom_coord - ca_coord
                        dist = (vec[0]**2 + vec[1]**2 + vec[2]**2)**0.5
                        if dist > 0.01:
                            h_coord = ca_coord + (CH_BOND_LENGTH / dist) * vec
                            # Use source atom's serial number (MCPB behavior)
                            serial = atom.serial_number
                            charge = self._get_charge_from_ff('ACE', h_names[h_count])
                            atoms.append((h_names[h_count], 'H', h_coord[0], h_coord[1], h_coord[2], charge, serial))
                            h_count += 1

            # Add C and O - use ACE charges from FF
            if c_coord is not None:
                charge = self._get_charge_from_ff('ACE', 'C')
                serial = self._get_serial_from_coords(tuple(c_coord))
                atoms.append(('C', 'C', c_coord[0], c_coord[1], c_coord[2], charge, serial))

            if o_coord is not None:
                charge = self._get_charge_from_ff('ACE', 'O')
                serial = self._get_serial_from_coords(tuple(o_coord))
                atoms.append(('O', 'O', o_coord[0], o_coord[1], o_coord[2], charge, serial))

        elif cap_type == CapType.NME:
            # NME: -NH-CH3
            # Strategy #2: All atoms use force field charges (not in atom_type_assignments)
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

            # Add N - use NME charge from FF
            if n_coord is not None:
                charge = self._get_charge_from_ff('NME', 'N')
                serial = self._get_serial_from_coords(tuple(n_coord))
                atoms.append(('N', 'N', n_coord[0], n_coord[1], n_coord[2], charge, serial))

            # Add H on nitrogen - use NME charge from FF
            if h_coord is not None:
                charge = self._get_charge_from_ff('NME', 'H')
                serial = self._get_serial_from_coords(tuple(h_coord))
                atoms.append(('H', 'H', h_coord[0], h_coord[1], h_coord[2], charge, serial))

            # Add CH3 and hydrogens - use NME charges from FF
            if ca_coord is not None:
                charge = self._get_charge_from_ff('NME', 'CH3')
                serial = self._get_serial_from_coords(tuple(ca_coord))
                atoms.append(('CH3', 'C', ca_coord[0], ca_coord[1], ca_coord[2], charge, serial))

                # Generate H atoms on methyl - use NME charges from FF
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
                            # Use source atom's serial number (MCPB behavior)
                            serial = atom.serial_number
                            charge = self._get_charge_from_ff('NME', h_names[h_count])
                            atoms.append((h_names[h_count], 'H', h_coord[0], h_coord[1], h_coord[2], charge, serial))
                            h_count += 1

        elif cap_type == CapType.GLY:
            # GLY: Keep all backbone atoms (N, H, CA, HA2, HA3, C, O)
            # Strategy #3: Primary=atom_type_assignments, Fallback=force field
            import numpy as np

            ca_coord = None
            n_coord = None
            c_coord = None
            has_n_hydrogen = False
            ha_atoms = []

            for atom in residue.get_atoms():
                atom_id = atom.get_id()
                coords = atom.get_coord()
                charge = self._get_charge_from_coords(tuple(coords), atom_type_assignments)
                # Strategy #3: Fall back to force field if not in atom_type_assignments
                # (handles GLY spacers in large model that aren't part of metal site)
                if charge == 0.0 and self.ff_data is not None:
                    charge = self._get_charge_from_ff('GLY', atom_id if atom_id != 'HN' else 'H')
                serial = self._get_serial_from_coords(tuple(coords))

                if atom_id == 'CA':
                    ca_coord = coords
                    atoms.append(('CA', atom.element, coords[0], coords[1], coords[2], charge, serial))
                elif atom_id == 'N':
                    n_coord = coords
                    atoms.append(('N', atom.element, coords[0], coords[1], coords[2], charge, serial))
                elif atom_id == 'C':
                    c_coord = coords
                    atoms.append(('C', atom.element, coords[0], coords[1], coords[2], charge, serial))
                elif atom_id == 'O':
                    atoms.append(('O', atom.element, coords[0], coords[1], coords[2], charge, serial))
                elif atom_id in ['H', 'HN']:
                    has_n_hydrogen = True
                    atoms.append(('H', atom.element, coords[0], coords[1], coords[2], charge, serial))
                elif atom_id in ['HA', 'HA2', 'HA3']:
                    ha_atoms.append((atom_id, atom.element, coords))

            # Generate N-H hydrogen if missing (e.g., PRO→GLY conversion)
            if n_coord is not None and not has_n_hydrogen:
                prev_residue = self.residue_map.get((model_residue.chain, model_residue.resid - 1))
                prev_c_coord = None
                if prev_residue:
                    for atom in prev_residue.get_atoms():
                        if atom.get_id() == 'C':
                            prev_c_coord = atom.get_coord()
                            break

                charge = self._get_charge_from_ff('GLY', 'H')
                if ca_coord is not None and prev_c_coord is not None:
                    v_nc = prev_c_coord - n_coord
                    v_nca = ca_coord - n_coord
                    v_nc = v_nc / np.linalg.norm(v_nc)
                    v_nca = v_nca / np.linalg.norm(v_nca)
                    bisector = v_nc + v_nca
                    bisector = bisector / np.linalg.norm(bisector)
                    h_direction = -bisector
                    h_coord = n_coord + NH_BOND_LENGTH * h_direction
                    serial = self._get_next_serial()
                    atoms.append(('H', 'H', h_coord[0], h_coord[1], h_coord[2], charge, serial))
                elif ca_coord is not None:
                    vec = n_coord - ca_coord
                    dist = np.linalg.norm(vec)
                    if dist > 0.01:
                        h_coord = n_coord + (NH_BOND_LENGTH / dist) * vec
                        serial = self._get_next_serial()
                        atoms.append(('H', 'H', h_coord[0], h_coord[1], h_coord[2], charge, serial))

            # Process H-alpha atoms
            if ha_atoms:
                has_ha2 = any(name == 'HA2' for name, _, _ in ha_atoms)
                has_ha3 = any(name == 'HA3' for name, _, _ in ha_atoms)

                for atom_id, element, coords in ha_atoms:
                    charge = self._get_charge_from_coords(tuple(coords), atom_type_assignments)
                    # Fall back to force field if not in atom_type_assignments
                    if charge == 0.0 and self.ff_data is not None:
                        # For HA→HA2/HA3 renaming, use HA2 charge from FF
                        if atom_id == 'HA':
                            charge = self._get_charge_from_ff('GLY', 'HA2')
                        else:
                            charge = self._get_charge_from_ff('GLY', atom_id)
                    serial = self._get_serial_from_coords(tuple(coords))

                    if atom_id == 'HA':
                        if not has_ha2:
                            atoms.append(('HA2', element, coords[0], coords[1], coords[2], charge, serial))
                            has_ha2 = True
                        elif not has_ha3:
                            atoms.append(('HA3', element, coords[0], coords[1], coords[2], charge, serial))
                            has_ha3 = True
                    else:
                        atoms.append((atom_id, element, coords[0], coords[1], coords[2], charge, serial))

                # Generate missing H-alpha if needed (e.g., PHE→GLY conversion)
                if ca_coord is not None and (not has_ha2 or not has_ha3):
                    for atom in residue.get_atoms():
                        if atom.get_id() == 'CB':
                            cb_coord = atom.get_coord()
                            vec = cb_coord - ca_coord
                            dist = np.linalg.norm(vec)
                            if dist > 0.01:
                                h_coord = ca_coord + (CH_BOND_LENGTH / dist) * vec
                                missing_name = 'HA3' if has_ha2 else 'HA2'
                                # Use CB atom's serial number (MCPB convention for generated H atoms)
                                serial = atom.serial_number
                                charge = self._get_charge_from_ff('GLY', missing_name)
                                atoms.append((missing_name, 'H', h_coord[0], h_coord[1], h_coord[2], charge, serial))
                                break

        elif cap_type == CapType.FULL:
            # FULL: Keep all atoms (backbone + sidechain) - for gap residues in connected chains
            # Strategy: Use atom_type_assignments for charges, fall back to force field
            for atom in residue.get_atoms():
                atom_id = atom.get_id()
                coords = atom.get_coord()
                element = atom.element
                charge = self._get_charge_from_coords(tuple(coords), atom_type_assignments)
                # Fall back to force field if not in atom_type_assignments
                if charge == 0.0 and self.ff_data is not None:
                    charge = self._get_charge_from_ff(model_residue.resname, atom_id)
                serial = self._get_serial_from_coords(tuple(coords))
                atoms.append((atom_id, element, coords[0], coords[1], coords[2], charge, serial))

        elif cap_type == CapType.KNH:
            # KNH: Keep N, H, sidechain; CA→CH3; remove C, O
            # Strategy #5: Generated H atoms get sum of replaced atom charges
            ca_coord = None
            h_sources = {}
            backbone_atoms = {}

            for atom in residue.get_atoms():
                atom_id = atom.get_id()
                coords = tuple(atom.get_coord())
                if atom_id == 'CA':
                    ca_coord = atom.get_coord()
                elif atom_id in ['C', 'HA']:
                    h_sources[atom_id] = atom.get_coord()
                if atom_id in ['C', 'O', 'HA']:
                    backbone_atoms[atom_id] = coords

            for atom in residue.get_atoms():
                atom_id = atom.get_id()
                coords = atom.get_coord()
                element = atom.element
                charge = self._get_charge_from_coords(tuple(coords), atom_type_assignments)
                serial = self._get_serial_from_coords(tuple(coords))

                if atom_id == 'CA':
                    atoms.append(('CH3', 'C', coords[0], coords[1], coords[2], charge, serial))
                elif atom_id in ['N', 'H', 'HN']:
                    if atom_id == 'HN':
                        atom_id = 'H'
                    atoms.append((atom_id, element, coords[0], coords[1], coords[2], charge, serial))
                elif atom_id not in ['C', 'O', 'OXT', 'HA']:
                    atoms.append((atom_id, element, coords[0], coords[1], coords[2], charge, serial))

            # Generate H atoms with charges that preserve residue neutrality
            # Strategy #5 (KNH): H replaces groups: C+O, HA (charge conservation)
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

                            # Calculate charge: sum of replaced atoms
                            if source_name == 'C':
                                # H replaces C=O group
                                charge = 0.0
                                if 'C' in backbone_atoms:
                                    charge += self._get_charge_from_coords(backbone_atoms['C'], atom_type_assignments)
                                if 'O' in backbone_atoms:
                                    charge += self._get_charge_from_coords(backbone_atoms['O'], atom_type_assignments)
                            else:  # HA
                                # H replaces HA
                                charge = self._get_charge_from_coords(backbone_atoms.get('HA', tuple(source_coord)), atom_type_assignments)

                            # Use source atom's serial number
                            serial = None
                            for atom in residue.get_atoms():
                                if atom.get_id() == source_name:
                                    serial = atom.serial_number
                                    break
                            if serial is None:
                                serial = self._get_next_serial()
                            atoms.append((h_names[h_count], 'H', h_coord[0], h_coord[1], h_coord[2], charge, serial))
                            h_count += 1

        elif cap_type == CapType.KCO:
            # KCO: Keep C, O, sidechain; CA→CH3; remove N, H
            # Strategy #6: Generated H atoms get sum of replaced atom charges
            ca_coord = None
            h_sources = {}
            backbone_atoms = {}

            for atom in residue.get_atoms():
                atom_id = atom.get_id()
                coords = tuple(atom.get_coord())
                if atom_id == 'CA':
                    ca_coord = atom.get_coord()
                elif atom_id in ['N', 'HA']:
                    h_sources[atom_id] = atom.get_coord()
                if atom_id in ['N', 'H', 'HN', 'HA']:
                    backbone_atoms[atom_id] = coords

            for atom in residue.get_atoms():
                atom_id = atom.get_id()
                coords = atom.get_coord()
                element = atom.element
                charge = self._get_charge_from_coords(tuple(coords), atom_type_assignments)
                serial = self._get_serial_from_coords(tuple(coords))

                if atom_id == 'CA':
                    atoms.append(('CH3', 'C', coords[0], coords[1], coords[2], charge, serial))
                elif atom_id in ['C', 'O']:
                    atoms.append((atom_id, element, coords[0], coords[1], coords[2], charge, serial))
                elif atom_id not in ['N', 'H', 'HN', 'HA', 'H1', 'H2', 'H3', 'HN1', 'HN2', 'HN3']:
                    atoms.append((atom_id, element, coords[0], coords[1], coords[2], charge, serial))

            # Generate H atoms with charges that preserve residue neutrality
            # Strategy #6 (KCO): H replaces groups: N+H, HA (charge conservation)
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

                            # Calculate charge: sum of replaced atoms
                            if source_name == 'N':
                                # H replaces N+H group
                                charge = 0.0
                                if 'N' in backbone_atoms:
                                    charge += self._get_charge_from_coords(backbone_atoms['N'], atom_type_assignments)
                                if 'H' in backbone_atoms:
                                    charge += self._get_charge_from_coords(backbone_atoms['H'], atom_type_assignments)
                                elif 'HN' in backbone_atoms:
                                    charge += self._get_charge_from_coords(backbone_atoms['HN'], atom_type_assignments)
                            else:  # HA
                                # H replaces HA
                                charge = self._get_charge_from_coords(backbone_atoms.get('HA', tuple(source_coord)), atom_type_assignments)

                            # Use source atom's serial number
                            serial = None
                            for atom in residue.get_atoms():
                                if atom.get_id() == source_name:
                                    serial = atom.serial_number
                                    break
                            if serial is None:
                                serial = self._get_next_serial()
                            atoms.append((h_names[h_count], 'H', h_coord[0], h_coord[1], h_coord[2], charge, serial))
                            h_count += 1

        return atoms

    def _extract_sidechain_atoms(self, model_residue: ModelResidue,
                                 atom_type_assignments: Dict) -> List[Tuple[str, str, float, float, float, float, int]]:
        """
        Extract sidechain-only atoms with CA→CH3 transformation.

        Strategy #4: Generated H atoms (H1/H2/H3) get sum of replaced backbone atom charges
        to maintain charge neutrality.

        Returns:
            List of (atom_name, element, x, y, z, charge, serial_number) tuples
        """
        residue = self.residue_map.get((model_residue.chain, model_residue.resid))
        if residue is None:
            return []

        atoms = []
        CH_BOND_LENGTH = 1.09
        BACKBONE_ATOMS = {'N', 'C', 'O', 'H', 'HN', 'OXT', 'HA', 'H1', 'H2', 'H3',
                          'HN1', 'HN2', 'HN3'}

        # Find CA
        ca_coord = None
        for atom in residue.get_atoms():
            if atom.get_id() == 'CA':
                ca_coord = atom.get_coord()
                break

        if ca_coord is None:
            # No CA - keep all atoms
            for atom in residue.get_atoms():
                coords = atom.get_coord()
                charge = self._get_charge_from_coords(tuple(coords), atom_type_assignments)
                serial = self._get_serial_from_coords(tuple(coords))
                atoms.append((atom.get_id(), atom.element, coords[0], coords[1], coords[2], charge, serial))
            return atoms

        # Collect H source positions and atoms for charge calculation
        h_sources = {}
        backbone_atoms = {}  # Store coords for N, H, C, O, HA
        for atom in residue.get_atoms():
            atom_id = atom.get_id()
            coords = tuple(atom.get_coord())
            if atom_id in ['N', 'C', 'HA']:
                h_sources[atom_id] = atom.get_coord()
            if atom_id in ['N', 'H', 'HN', 'C', 'O', 'HA']:
                backbone_atoms[atom_id] = coords

        # Add CA as CH3
        charge = self._get_charge_from_coords(tuple(ca_coord), atom_type_assignments)
        serial = self._get_serial_from_coords(tuple(ca_coord))
        atoms.append(('CH3', 'C', ca_coord[0], ca_coord[1], ca_coord[2], charge, serial))

        # Generate 3 H atoms with charges that preserve residue neutrality
        # Strategy #4: H replaces groups: N+H, C+O, HA (charge conservation)
        h_names = ['H1', 'H2', 'H3']
        h_count = 0
        for source_name in ['N', 'C', 'HA']:
            if source_name in h_sources and h_count < 3:
                source_coord = h_sources[source_name]
                vec = source_coord - ca_coord
                dist = (vec[0]**2 + vec[1]**2 + vec[2]**2)**0.5
                if dist > 0.01:
                    h_coord = ca_coord + (CH_BOND_LENGTH / dist) * vec

                    # Calculate charge: sum of replaced atoms (maintains residue charge)
                    if source_name == 'N':
                        # H replaces N+H group
                        charge = 0.0
                        if 'N' in backbone_atoms:
                            charge += self._get_charge_from_coords(backbone_atoms['N'], atom_type_assignments)
                        if 'H' in backbone_atoms:
                            charge += self._get_charge_from_coords(backbone_atoms['H'], atom_type_assignments)
                        elif 'HN' in backbone_atoms:
                            charge += self._get_charge_from_coords(backbone_atoms['HN'], atom_type_assignments)
                    elif source_name == 'C':
                        # H replaces C=O group
                        charge = 0.0
                        if 'C' in backbone_atoms:
                            charge += self._get_charge_from_coords(backbone_atoms['C'], atom_type_assignments)
                        if 'O' in backbone_atoms:
                            charge += self._get_charge_from_coords(backbone_atoms['O'], atom_type_assignments)
                    else:  # HA
                        # H replaces HA
                        charge = self._get_charge_from_coords(backbone_atoms.get('HA', source_coord), atom_type_assignments)

                    # Use source atom's serial number
                    serial = None
                    for atom in residue.get_atoms():
                        if atom.get_id() == source_name:
                            serial = atom.serial_number
                            break
                    if serial is None:
                        serial = self._get_next_serial()

                    atoms.append((h_names[h_count], 'H', h_coord[0], h_coord[1], h_coord[2], charge, serial))
                    h_count += 1

        # Add sidechain atoms
        for atom in residue.get_atoms():
            atom_name = atom.get_id()
            if atom_name not in BACKBONE_ATOMS and atom_name != 'CA':
                coords = atom.get_coord()
                charge = self._get_charge_from_coords(tuple(coords), atom_type_assignments)
                serial = self._get_serial_from_coords(tuple(coords))
                atoms.append((atom_name, atom.element, coords[0], coords[1], coords[2], charge, serial))

        return atoms

    @staticmethod
    def _coord_key(coords, ndigits: int = 3):
        """A coordinate tuple that compares equal across float widths.

        BioPython hands back float32; assignments that have been through JSON
        (a resumed session, a saved atom_type_assignments.json) come back as
        float64. The two are not equal and do not hash alike, yet numpy prints
        the shortest repr that round-trips in float32, so both display as
        ``(-41.416, -16.455, -50.088)`` while ``==`` is False.

        That made every metal-site charge lookup miss silently: the writer
        summed 0.0 for all 44 site atoms, the small model's suggested QM charge
        came out 0 instead of -2, and nothing in the output looked wrong --
        a table of 0.0000 is what a withheld cluster atom legitimately shows.

        3 decimals is the precision a PDB file carries, so rounding there is
        lossless for coordinates that came from one.
        """
        try:
            return tuple(round(float(c), ndigits) for c in coords)
        except (TypeError, ValueError):
            return None

    def _assignment_by_coord_key(self, atom_type_assignments: Dict, coords):
        """Look up an assignment by rounded coordinate, indexing on first use."""
        if self._assignment_index_src is not atom_type_assignments:
            index = {}
            for key, value in atom_type_assignments.items():
                normalized = self._coord_key(key)
                if normalized is not None:
                    index[normalized] = value
            self._assignment_index = index
            self._assignment_index_src = atom_type_assignments

        normalized = self._coord_key(coords)
        if normalized is None:
            return None
        return self._assignment_index.get(normalized)

    def _get_charge_from_coords(self, coords: Tuple[float, float, float],
                                atom_type_assignments: Dict) -> float:
        """
        Get charge for atom from atom_type_assignments.

        Strategy #7: Used for real atoms that were atom-typed during MCPB Step 1.
        Returns 0.0 if atom not in assignments (caller should fall back to force field).

        Falls back to a rounded-coordinate match when the exact tuple misses;
        see _coord_key for why an exact match is not reliable here.
        """
        assignment = atom_type_assignments.get(coords)
        if assignment is None:
            assignment = self._assignment_by_coord_key(atom_type_assignments, coords)
        if assignment is not None:
            # Support both dataclass and dict. Guard None explicitly: a withheld
            # cluster atom (e.g. an Fe-S bridging sulfide, removed from the prmtop
            # with the rest of its cluster) has an assignment whose 'charge' key
            # exists but is None — .get('charge', 0.0) would return that None and
            # crash the total-charge sum. Treat it as 0.0 like every other
            # unknown charge in this writer.
            if hasattr(assignment, 'charge'):
                return assignment.charge if assignment.charge is not None else 0.0
            elif isinstance(assignment, dict):
                charge = assignment.get('charge')
                return charge if charge is not None else 0.0
        return 0.0

    def _get_charge_from_ff(self, resname: str, atom_name: str) -> float:
        """
        Get charge for atom from force field data.

        Strategy #2 (ACE/NME), #3 (GLY fallback): Used for capping groups and
        spacer residues not in atom_type_assignments.

        Uses unified ForceFieldData which includes:
        - lib file data (amino12.lib, conste.lib, etc.)
        - mol2 file data (parmXX.mol2)
        - Residue aliases (e.g., HEM → HEH)

        Args:
            resname: Residue name (e.g., 'ACE', 'NME', 'GLY', 'HEM', 'HEH')
            atom_name: Atom name (e.g., 'HH31', 'C', 'O', 'FE')

        Returns:
            Charge from force field, or 0.0 if not found
        """
        if self.ff_data is None:
            return 0.0

        # ForceFieldData handles alias lookup internally
        charge = self.ff_data.get_atom_charge(resname, atom_name)
        if charge is not None:
            return charge

        return 0.0

    def _get_serial_from_coords(self, coords: Tuple[float, float, float]) -> int:
        """Get PDB serial number from coordinates, or assign new serial."""
        if coords in self.atom_map:
            return self.atom_map[coords].serial_number
        return self._get_next_serial()

    def _get_next_serial(self) -> int:
        """Get next available serial number for generated atoms."""
        serial = self.next_serial
        self.next_serial += 1
        return serial

    def _print_charge_breakdown(self, atoms_data: List[Tuple], total_charge: float):
        """
        Print detailed charge breakdown table.

        Args:
            atoms_data: List of (serial, name, resname, chain, resid, x, y, z, element, charge)
            total_charge: Total charge of system
        """
        from rich.table import Table

        self.console.print("\n[bold cyan]Charge Breakdown:[/bold cyan]")

        table = Table(show_header=True, header_style="bold magenta")
        table.add_column("Serial", justify="right", style="cyan")
        table.add_column("Atom", justify="left", style="green")
        table.add_column("Res", justify="left", style="yellow")
        table.add_column("ResID", justify="right", style="yellow")
        table.add_column("Charge", justify="right", style="white")
        table.add_column("Running Sum", justify="right", style="blue")

        running_sum = 0.0
        for serial, name, resname, chain, resid, x, y, z, element, charge in atoms_data:
            running_sum += charge
            table.add_row(
                f"{serial}",
                f"{name}",
                f"{resname}",
                f"{resid}",
                f"{charge:>7.4f}",
                f"{running_sum:>7.4f}"
            )

        self.console.print(table)
        self.console.print(f"\n[bold]Total Charge: {total_charge:.4f}[/bold]\n")

    def _derive_conect_bonds(self, atoms_data: List[Tuple],
                              redox_site: Optional[RedoxSite]) -> List[Tuple[int, int]]:
        """
        Derive bond list for CONECT records.

        Within-residue bonds are detected by distance (covers C-H, C-C, C-N,
        C-O, C-S, N-H, O-H, S-H single bonds up to ~1.95 A). Cross-residue
        bonds are taken from redox_site.bonds (peptide bonds and metal-ligand
        coordinate bonds get re-derived as distance-based pairs if not
        present, so RESP equivalence at least sees within-residue connectivity).

        Returns:
            Sorted list of unique (serial1, serial2) tuples with serial1 < serial2.
        """
        bonds: Set[Tuple[int, int]] = set()

        # Group atoms by (chain, resid). Heavy-atom-to-H bond cutoff is ~1.3 A;
        # C-S can reach 1.81 A; use 1.95 A to capture C-S and miss non-bonded contacts.
        residue_groups: Dict[Tuple[str, int], List[Tuple]] = {}
        for entry in atoms_data:
            serial, _name, _resname, chain, resid, _x, _y, _z, _element, _charge = entry
            residue_groups.setdefault((chain, resid), []).append(entry)

        cutoff_sq = 1.95 * 1.95
        for atoms in residue_groups.values():
            for i, a in enumerate(atoms):
                sa, _, _, _, _, ax, ay, az, _, _ = a
                for j in range(i + 1, len(atoms)):
                    sb, _, _, _, _, bx, by, bz, _, _ = atoms[j]
                    d2 = (ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2
                    if d2 < cutoff_sq:
                        pair = (sa, sb) if sa < sb else (sb, sa)
                        bonds.add(pair)

        # Add metal-ligand bonds from the redox_site (longer than 1.95 A so
        # not caught by distance scan above).
        if redox_site is not None:
            coord_to_serial: Dict[Tuple[float, float, float], int] = {}
            for entry in atoms_data:
                serial, _name, _resname, _chain, _resid, x, y, z, _element, _charge = entry
                coord_to_serial[(round(x, 3), round(y, 3), round(z, 3))] = serial

            for bond in getattr(redox_site, "bonds", []):
                k1 = tuple(round(c, 3) for c in bond.atom1_coords)
                k2 = tuple(round(c, 3) for c in bond.atom2_coords)
                s1 = coord_to_serial.get(k1)
                s2 = coord_to_serial.get(k2)
                if s1 is not None and s2 is not None and s1 != s2:
                    pair = (s1, s2) if s1 < s2 else (s2, s1)
                    bonds.add(pair)

        return sorted(bonds)

    def _write_pdb_file(self, atoms_data: List[Tuple],
                       total_charge: float, title: str, output_file: str,
                       conect_bonds: Optional[List[Tuple[int, int]]] = None):
        """
        Write atoms to PDB file.

        Args:
            atoms_data: List of (serial, name, resname, chain, resid, x, y, z, element, charge)
            total_charge: Total charge of system
            title: Title for HEADER
            output_file: Path to output file
            conect_bonds: Optional list of (serial1, serial2) bond pairs to emit
                as CONECT records before the END marker.
        """
        with open(output_file, 'w') as f:
            # HEADER
            f.write(f"HEADER    {title:<60}\n")
            f.write(f"REMARK    Total charge: {total_charge:.2f}\n")
            f.write(f"REMARK    Generated by ProPrep MCPB module\n")

            # ATOM records
            for serial, name, resname, chain, resid, x, y, z, element, charge in atoms_data:
                # PDB format: ATOM serial name resname chain resid x y z occ bfac element
                atom_line = (
                    f"ATOM  {serial:5d} {name:4s} {resname:3s} {chain:1s}{resid:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:2s}\n"
                )
                f.write(atom_line)

            # CONECT records: group all partners of each atom on one line.
            # RESP's _find_ch_groups reads CONECT to detect CH2/CH3 equivalence.
            if conect_bonds:
                partners: Dict[int, List[int]] = {}
                for s1, s2 in conect_bonds:
                    partners.setdefault(s1, []).append(s2)
                    partners.setdefault(s2, []).append(s1)
                for s in sorted(partners):
                    others = sorted(set(partners[s]))
                    # PDB CONECT format: 'CONECT' + 5d serial + up to 4 5d partners per record.
                    for chunk_start in range(0, len(others), 4):
                        chunk = others[chunk_start:chunk_start + 4]
                        record = "CONECT" + f"{s:5d}" + "".join(f"{p:5d}" for p in chunk) + "\n"
                        f.write(record)

            f.write("END\n")
