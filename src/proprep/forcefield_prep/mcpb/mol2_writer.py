"""
Mol2 Writer

Writes AMBER-compatible mol2 files with atom types and RESP charges.
Uses pre-computed bond topology from Step 2A.

Based on MCPB's print_mol2f() implementation.
"""

from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
import json
import logging


@dataclass
class Mol2Atom:
    """Atom data for mol2 file"""
    atom_id: int
    atom_name: str
    x: float
    y: float
    z: float
    atom_type: str
    residue_id: int
    residue_name: str
    charge: float
    element: str


class Mol2Writer:
    """
    Write AMBER-compatible mol2 files with atom types and RESP charges.

    Combines:
    - Coordinates from large.pdb
    - Atom types from fingerprint file
    - Charges from RESP fitting
    - Bonds from saved topology (Step 2A)
    """

    def __init__(self, logger=None):
        """
        Initialize mol2 writer.

        Args:
            logger: Optional logger instance
        """
        self.logger = logger or logging.getLogger(__name__)

    def write_mol2_files(self,
                        output_dir: str,
                        pdb_file: str,
                        fingerprint_file: str,
                        resp_charges: List[float],
                        bond_topology_file: str,
                        metal_site_resids: set = None) -> Dict[str, str]:
        """
        Write mol2 files for each residue.

        Args:
            output_dir: Directory for output mol2 files
            pdb_file: Path to large.pdb with coordinates
            fingerprint_file: Path to standard.fingerprint with atom types
            resp_charges: List of RESP charges (in PDB atom order)
            bond_topology_file: Path to bond_topology.json from Step 2A
            metal_site_resids: Set of (chain, resid) tuples for metal-site residues.
                              Only these get mol2 files. If None, all residues are written.

        Returns:
            Dict mapping residue_name → mol2_file_path

        Raises:
            ValueError: If charge count doesn't match atom count
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Parse PDB for coordinates and residue info
        pdb_atoms = self._parse_pdb_atoms(pdb_file)

        # Validate charge count
        if len(resp_charges) != len(pdb_atoms):
            raise ValueError(
                f"Charge count ({len(resp_charges)}) doesn't match atom count ({len(pdb_atoms)}). "
                "Check resp2.chg file."
            )

        # Parse fingerprint for atom types
        atom_types = self._parse_fingerprint(fingerprint_file)

        # Load bond topology
        bonds_all = self._load_bond_topology(bond_topology_file)

        # Group atoms by residue
        residue_groups = self._group_atoms_by_residue(pdb_atoms)

        self.logger.debug(f"Writing mol2 files for {len(residue_groups)} residues...")

        # Write mol2 files only for metal-site residues (from RedoxSite definition).
        # Scaffolding residues (ACE, NME, gap fillers) exist only in the large model
        # for RESP fitting; their charges are fitting artifacts that get discarded.
        mol2_files = {}
        for residue_key, residue_atoms in residue_groups.items():
            chain, resid, resname = residue_key

            if metal_site_resids and (chain, resid) not in metal_site_resids:
                self.logger.debug(f"  Skipping scaffolding residue {resname}{resid} (not in metal site)")
                continue

            # Generate mol2 filename
            mol2_name = f"{resname}{resid}.mol2"
            mol2_path = output_dir / mol2_name

            # Build mol2 atoms for this residue
            mol2_atoms, serial_to_mol2_id = self._build_mol2_atoms(
                residue_atoms, atom_types, resp_charges, resname
            )

            # Build this residue's connectivity. The Step-2 bond_topology comes
            # from the truncated small model and is INCOMPLETE for residues that
            # ligate through their backbone (sidechain left unbonded), so derive
            # the full intra-residue covalent graph by distance from the (full,
            # large-model) mol2 atoms and union in whatever the topology recorded.
            derived = self._derive_intraresidue_bonds(residue_atoms, serial_to_mol2_id)
            topo = self._filter_bonds_for_residue(
                bonds_all, residue_atoms, serial_to_mol2_id
            )
            seen = set()
            bonds_for_mol2 = []
            for pair in list(derived) + list(topo):
                key = tuple(sorted(pair))
                if key not in seen:
                    seen.add(key)
                    bonds_for_mol2.append(pair)

            # Write mol2 file
            self._write_mol2_format(mol2_path, resname, mol2_atoms, bonds_for_mol2)

            mol2_files[f"{resname}{resid}"] = str(mol2_path)
            self.logger.debug(f"  Written {mol2_name}: {len(mol2_atoms)} atoms, {len(bonds_for_mol2)} bonds")

        return mol2_files

    def _parse_pdb_atoms(self, pdb_file: str) -> List[Dict]:
        """
        Parse atom data from PDB file.

        Args:
            pdb_file: Path to PDB file

        Returns:
            List of dicts with atom info
        """
        atoms = []
        pdb_position = 0  # 0-based index matching RESP charge ordering

        with open(pdb_file) as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    serial = int(line[6:11].strip())
                    atom_name = line[12:16].strip()
                    resname = line[17:20].strip()
                    chain = line[21] if len(line) > 21 else ' '
                    resid = int(line[22:26].strip())
                    element = line[76:78].strip() if len(line) > 77 else atom_name[0]
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())

                    atoms.append({
                        'serial': serial,
                        'atom_name': atom_name,
                        'resname': resname,
                        'chain': chain,
                        'resid': resid,
                        'element': element,
                        'x': x,
                        'y': y,
                        'z': z,
                        'pdb_position': pdb_position
                    })
                    pdb_position += 1

        return atoms

    def _parse_fingerprint(self, fingerprint_file: str) -> Dict[int, str]:
        """
        Parse atom types from fingerprint file.

        Format: ResID-ResName-AtomName  AtomID  OrigType -> RenamedType

        Args:
            fingerprint_file: Path to .fingerprint file

        Returns:
            Dict mapping PDB serial → atom_type (renamed)
        """
        atom_types = {}

        with open(fingerprint_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('LINK'):
                    continue

                parts = line.split()
                if len(parts) >= 5 and '->' in line:
                    try:
                        atom_id = int(parts[1])
                        renamed_type = parts[4].strip()
                        atom_types[atom_id] = renamed_type
                    except (ValueError, IndexError):
                        pass

        return atom_types

    def _load_bond_topology(self, topology_file: str) -> List[Tuple[int, int]]:
        """
        Load bond topology from JSON.

        Args:
            topology_file: Path to bond_topology.json

        Returns:
            List of (serial1, serial2) bond tuples (PDB serial numbers)
        """
        with open(topology_file) as f:
            data = json.load(f)

        return [tuple(b) for b in data["bonds"]]

    def _group_atoms_by_residue(self, pdb_atoms: List[Dict]) -> Dict[Tuple, List[Dict]]:
        """
        Group atoms by residue.

        Args:
            pdb_atoms: List of PDB atom dicts

        Returns:
            Dict mapping (chain, resid, resname) → list of atoms
        """
        residue_groups = {}

        for atom in pdb_atoms:
            key = (atom['chain'], atom['resid'], atom['resname'])
            if key not in residue_groups:
                residue_groups[key] = []
            residue_groups[key].append(atom)

        return residue_groups

    def _build_mol2_atoms(self,
                         residue_atoms: List[Dict],
                         atom_types: Dict[int, str],
                         resp_charges: List[float],
                         resname: str) -> Tuple[List[Mol2Atom], Dict[int, int]]:
        """
        Build Mol2Atom objects for a residue.

        Args:
            residue_atoms: List of atom dicts for this residue
            atom_types: Dict mapping PDB serial → atom type
            resp_charges: Full list of RESP charges (ordered by PDB atom position)
            resname: Residue name

        Returns:
            Tuple of (mol2_atoms list, serial_to_mol2_id mapping)
        """
        mol2_atoms = []
        serial_to_mol2_id = {}

        # Sort atoms by serial number for consistent ordering
        sorted_atoms = sorted(residue_atoms, key=lambda a: a['serial'])

        for mol2_id, atom in enumerate(sorted_atoms, 1):
            serial = atom['serial']
            atom_name = atom['atom_name']

            # Get atom type from fingerprint
            atom_type = atom_types.get(serial, 'XX')

            # Get charge from RESP using PDB position index (not serial number)
            # resp_charges are ordered by atom position in PDB, not by serial
            charge_idx = atom.get('pdb_position', None)
            if charge_idx is not None and charge_idx < len(resp_charges):
                charge = resp_charges[charge_idx]
            elif serial <= len(resp_charges):
                charge = resp_charges[serial - 1]
            else:
                charge = 0.0

            mol2_atom = Mol2Atom(
                atom_id=mol2_id,
                atom_name=atom_name,
                x=atom['x'],
                y=atom['y'],
                z=atom['z'],
                atom_type=atom_type,
                residue_id=1,  # All atoms in one residue for mol2
                residue_name=resname,
                charge=charge,
                element=atom['element']
            )
            mol2_atoms.append(mol2_atom)
            serial_to_mol2_id[serial] = mol2_id

        return mol2_atoms, serial_to_mol2_id

    def _filter_bonds_for_residue(self,
                                  bonds_all: List[Tuple[int, int]],
                                  residue_atoms: List[Dict],
                                  serial_to_mol2_id: Dict[int, int]) -> List[Tuple[int, int]]:
        """
        Filter bonds to only those within this residue.

        Args:
            bonds_all: All bonds (PDB serial numbers)
            residue_atoms: List of atom dicts for this residue
            serial_to_mol2_id: Mapping from PDB serial → mol2 ID

        Returns:
            List of (mol2_id1, mol2_id2) bond tuples for mol2 file
        """
        # Get set of serials in this residue
        residue_serials = {atom['serial'] for atom in residue_atoms}

        bonds_mol2 = []
        for serial1, serial2 in bonds_all:
            # Both atoms must be in this residue
            if serial1 in residue_serials and serial2 in residue_serials:
                # Convert to mol2 IDs
                mol2_id1 = serial_to_mol2_id.get(serial1)
                mol2_id2 = serial_to_mol2_id.get(serial2)

                if mol2_id1 and mol2_id2:
                    bonds_mol2.append((mol2_id1, mol2_id2))

        return bonds_mol2

    # Covalent radii (Cordero 2008, Å) for intra-residue bond perception.
    # A default covers anything unlisted.
    _COVALENT_RADII = {
        'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57,
        'P': 1.07, 'S': 1.05, 'CL': 1.02, 'BR': 1.20, 'I': 1.39,
        'SE': 1.20, 'B': 0.84,
    }
    _COVALENT_DEFAULT = 0.77
    _COVALENT_TOLERANCE = 0.45

    def _derive_intraresidue_bonds(self,
                                   residue_atoms: List[Dict],
                                   serial_to_mol2_id: Dict[int, int]) -> List[Tuple[int, int]]:
        """
        Derive a residue's COMPLETE intra-residue covalent connectivity by distance.

        The Step-2 bond_topology is built from the truncated SMALL QM model, so a
        residue that ligates the metal through its BACKBONE (e.g. an Ile whose
        carbonyl O binds the metal) appears there with backbone atoms only and its
        entire sidechain is left unbonded in the mol2 -> OFF lib -> topology
        ("Unbonded Hydrogen ..."). The mol2 atoms, however, come from the FULL
        large model, so perceive bonds directly from those coordinates.

        Uses covalent radii + a tolerance rather than a flat cutoff. This is both
        more accurate AND naturally excludes spurious H-H bonds (geminal hydrogens
        sit ~1.78 A apart, well beyond 0.31 + 0.31 + 0.45 = 1.07 A), which a flat
        ~1.9 A cutoff would wrongly bond. Metal-ligand and peptide bonds are
        inter-residue and are added by tleap separately, so they are omitted here.
        """
        atoms = list(residue_atoms)
        bonds = []
        for i, a in enumerate(atoms):
            id_a = serial_to_mol2_id.get(a['serial'])
            if id_a is None:
                continue
            ra = self._COVALENT_RADII.get((a.get('element') or '').strip().upper(), self._COVALENT_DEFAULT)
            for j in range(i + 1, len(atoms)):
                b = atoms[j]
                id_b = serial_to_mol2_id.get(b['serial'])
                if id_b is None:
                    continue
                rb = self._COVALENT_RADII.get((b.get('element') or '').strip().upper(), self._COVALENT_DEFAULT)
                cutoff = ra + rb + self._COVALENT_TOLERANCE
                d2 = (a['x'] - b['x']) ** 2 + (a['y'] - b['y']) ** 2 + (a['z'] - b['z']) ** 2
                if d2 < cutoff * cutoff:
                    bonds.append((id_a, id_b) if id_a < id_b else (id_b, id_a))
        return bonds

    def _write_mol2_format(self,
                          output_path: Path,
                          residue_name: str,
                          atoms: List[Mol2Atom],
                          bonds: List[Tuple[int, int]]) -> None:
        """
        Write mol2 file in AMBER format.

        Format (based on MCPB print_mol2f):
        @<TRIPOS>MOLECULE
        {residue_name}
        {num_atoms} {num_bonds}     1     0     0
        SMALL
        RESP Charge

        @<TRIPOS>ATOM
        {atom_id:7d} {atom_name:<4s}    {x:10.4f}{y:10.4f}{z:10.4f} {atom_type:<4s} {res_id:6d} {res_name:<4s} {charge:12.6f}
        ...

        @<TRIPOS>BOND
        {bond_id:6d}{atom1:5d}{atom2:5d}  1
        ...

        @<TRIPOS>SUBSTRUCTURE
             1 {residue_name}        1 TEMP              0 ****  ****    0 ROOT

        Args:
            output_path: Path for output file
            residue_name: Residue name
            atoms: List of Mol2Atom objects
            bonds: List of bond tuples
        """
        with open(output_path, 'w') as f:
            # MOLECULE section
            f.write("@<TRIPOS>MOLECULE\n")
            f.write(f"{residue_name}\n")
            f.write(f"{len(atoms):5d}{len(bonds):6d}     1     0     0\n")
            f.write("SMALL\n")
            f.write("RESP Charge\n")
            f.write(" \n")
            f.write(" \n")

            # ATOM section
            f.write("@<TRIPOS>ATOM\n")
            for atom in atoms:
                f.write(
                    f"{atom.atom_id:7d} {atom.atom_name:<4s}    "
                    f"{atom.x:10.4f}{atom.y:10.4f}{atom.z:10.4f} "
                    f"{atom.atom_type:<4s} {atom.residue_id:6d} {atom.residue_name:<4s} "
                    f"{atom.charge:12.6f}\n"
                )

            # BOND section
            f.write("@<TRIPOS>BOND\n")
            for bond_id, (atom1, atom2) in enumerate(bonds, 1):
                f.write(f"{bond_id:6d}{atom1:5d}{atom2:5d}  1\n")  # All single bonds

            # SUBSTRUCTURE section
            f.write("@<TRIPOS>SUBSTRUCTURE\n")
            f.write(f"     1 {residue_name}        1 TEMP              0 ****  ****    0 ROOT\n")
