"""
RESP Input Generator

Generates resp input files for two-stage RESP fitting.
Implements backbone restraints and custom charge constraints.

Based on MCPB's gene_resp_input_file() implementation.
"""

from pathlib import Path
from typing import List, Dict, Optional, Tuple, Set
from dataclasses import dataclass
import logging


@dataclass
class RESPRestraint:
    """Configuration for RESP charge restraints"""
    atom_index: int      # 1-based atom index
    charge: float        # Fixed charge value
    restraint_type: str  # 'fixed' or 'restrained'


class RESPInputGenerator:
    """
    Generate resp input files for two-stage RESP fitting.

    Implements:
    - Stage 1: qwt=0.0005 (weak restraint), all atoms free
    - Stage 2: qwt=0.001 (strong restraint), equivalences enforced
    - Backbone restraints via chgmod parameter
    - Custom fixed atom charges
    - Automatic equivalence detection for CH2/CH3 groups
    """

    # Backbone atoms for chgmod restraints
    BACKBONE_ATOMS = {
        0: [],  # No restraints
        1: ['CA', 'N', 'C', 'O'],
        2: ['CA', 'H', 'HA', 'N', 'C', 'O'],
        3: ['CA', 'H', 'HA', 'N', 'C', 'O', 'CB']
    }

    def __init__(self, logger=None):
        """
        Initialize RESP input generator.

        Args:
            logger: Optional logger instance
        """
        self.logger = logger or logging.getLogger(__name__)

    def generate_resp_inputs(self,
                           pdb_file: str,
                           esp_file: str,
                           output_dir: str,
                           chgmod: int = 1,
                           custom_restraints: Optional[List[RESPRestraint]] = None,
                           total_charge: int = 0,
                           ff_collector = None,
                           type_assignments: Optional[Dict] = None,
                           metal_site_resids: set = None,
                           cross_residue_equivalence_groups: Optional[List[List[Tuple[str, int]]]] = None) -> Tuple[str, str]:
        """
        Generate resp1.in and resp2.in files.

        Args:
            pdb_file: Path to large.pdb file (for atom identification)
            esp_file: Path to .esp file from ESP extraction
            output_dir: Directory for output files
            chgmod: Backbone restraint level (0-3)
            custom_restraints: Optional list of custom atom restraints
            total_charge: Total molecular charge
            ff_collector: ForcefieldDataCollector from Step 1 (for backbone charges)
            type_assignments: Dict of coord -> assignment dict (backup for charges from prmtop)

        Returns:
            Tuple of (respin1_path, respin2_path)

        Raises:
            ValueError: If chgmod not in [0, 1, 2, 3]
        """
        if chgmod not in [0, 1, 2, 3]:
            raise ValueError(f"chgmod must be 0, 1, 2, or 3 (got {chgmod})")

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.logger.debug("Generating RESP input files...")

        # Parse PDB to get atom information
        atoms_info = self._parse_pdb_atoms(pdb_file)
        n_atoms = len(atoms_info)

        # Identify scaffolding vs metal-site residues
        # Metal-site residues (from RedoxSite): get free RESP fitting + backbone restraints
        # Scaffolding (ACE, NME, gap fillers): constrained to sum to 0, no mol2 output
        residue_atoms_by_key = {}
        for idx, atom in enumerate(atoms_info, 1):
            key = (atom['chain'], atom['resid'])
            if key not in residue_atoms_by_key:
                residue_atoms_by_key[key] = {'resname': atom['resname'], 'indices': []}
            residue_atoms_by_key[key]['indices'].append(idx)

        # Build group constraints for scaffolding residues (sum to 0)
        residue_groups = []
        scaffolding_resids = set()
        for (chain, resid), info in residue_atoms_by_key.items():
            is_metal_site = metal_site_resids and (chain, resid) in metal_site_resids
            if not is_metal_site:
                residue_groups.append(info['indices'])
                scaffolding_resids.add((chain, resid))

        if residue_groups:
            self.logger.debug(f"  Scaffolding constraints: {len(residue_groups)} residues (sum to 0)")
        if metal_site_resids:
            self.logger.debug(f"  Metal-site residues (free fitting): {len(metal_site_resids)}")

        # Identify frozen atoms for capping groups (MCPB convention)
        # ACE: C, O frozen (-99 in stage 2)
        # NME: N, H frozen (-99 in stage 2)
        frozen_cap_atoms = set()  # 1-based indices
        for idx, atom in enumerate(atoms_info, 1):
            if atom['resname'] == 'ACE' and atom['atom_name'] in ('C', 'O'):
                frozen_cap_atoms.add(idx)
            elif atom['resname'] == 'NME' and atom['atom_name'] in ('N', 'H'):
                frozen_cap_atoms.add(idx)

        # Identify backbone atoms for restraints (Part 4 — only for metal-site residues)
        backbone_restraints = {}
        if chgmod > 0:
            # Scaffolding residues get group-sum constraints instead of per-atom
            # backbone restraints, so exclude them at lookup time (this also
            # suppresses charge-lookup warnings for cap atoms like NME:N).
            backbone_restraints = self._get_backbone_restraints(
                atoms_info, chgmod, ff_collector, type_assignments,
                scaffolding_resids=scaffolding_resids
            )
            self.logger.debug(f"  Backbone restraints: {len(backbone_restraints)} atoms (chgmod={chgmod})")

        # Add custom restraints
        all_restraints = backbone_restraints.copy()
        if custom_restraints:
            for restraint in custom_restraints:
                all_restraints[restraint.atom_index] = restraint.charge
            self.logger.debug(f"  Custom restraints: {len(custom_restraints)} atoms")

        # Detect equivalences for stage 2
        equivalences = self._detect_equivalences(pdb_file, atoms_info)
        self.logger.debug(f"  Auto-detected equivalences: {len(equivalences)} groups")

        # Layer user-specified cross-residue equivalences on top, then merge
        # transitively-connected groups via union-find so that within-residue
        # and cross-residue equivalences compose correctly.
        if cross_residue_equivalence_groups:
            cross_groups = self._expand_cross_residue_equivalences(
                atoms_info, cross_residue_equivalence_groups
            )
            self.logger.debug(f"  Cross-residue equivalences: {len(cross_groups)} groups")
            equivalences = self._merge_equivalence_groups(equivalences + cross_groups)
            self.logger.debug(f"  Total after merging: {len(equivalences)} groups")

        # Write resp1.in (stage 1)
        respin1_path = output_dir / "resp1.in"
        self._write_respin1(respin1_path, atoms_info, total_charge,
                           all_restraints, residue_groups)
        self.logger.debug(f"  Written: {respin1_path.name}")

        # Write resp2.in (stage 2)
        respin2_path = output_dir / "resp2.in"
        self._write_respin2(respin2_path, atoms_info, total_charge,
                           all_restraints, equivalences, residue_groups,
                           frozen_cap_atoms)
        self.logger.debug(f"  Written: {respin2_path.name}")

        return str(respin1_path), str(respin2_path)

    def _parse_pdb_atoms(self, pdb_file: str) -> List[Dict]:
        """
        Parse atom information from PDB file.

        Args:
            pdb_file: Path to PDB file

        Returns:
            List of dicts with atom info (serial, name, element, resname, resid, chain)
        """
        atoms = []
        with open(pdb_file) as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    serial = int(line[6:11].strip())
                    atom_name = line[12:16].strip()
                    resname = line[17:20].strip()
                    chain = line[21] if len(line) > 21 else ' '
                    resid = int(line[22:26].strip())
                    element = line[76:78].strip() if len(line) > 77 else atom_name[0]

                    atoms.append({
                        'serial': serial,
                        'atom_name': atom_name,
                        'resname': resname,
                        'chain': chain,
                        'resid': resid,
                        'element': element
                    })
        return atoms

    def _get_atomic_number(self, element: str) -> int:
        """
        Get atomic number from element symbol.

        Args:
            element: Element symbol (H, C, N, O, etc.)

        Returns:
            Atomic number

        Raises:
            ValueError: If element symbol is not recognized
        """
        # Complete periodic table (elements 1-118)
        atomic_numbers = {
            # Period 1
            'H': 1, 'He': 2,
            # Period 2
            'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
            # Period 3
            'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,
            # Period 4
            'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26,
            'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34,
            'Br': 35, 'Kr': 36,
            # Period 5
            'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44,
            'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52,
            'I': 53, 'Xe': 54,
            # Period 6
            'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62,
            'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70,
            'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78,
            'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86,
            # Period 7
            'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94,
            'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102,
            'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109,
            'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116,
            'Ts': 117, 'Og': 118
        }

        # Normalize element symbol (handle both 'C' and 'c', 'CA' and 'Ca')
        element_clean = element.strip()
        if len(element_clean) == 1:
            element_key = element_clean.upper()
        else:
            element_key = element_clean[0].upper() + element_clean[1:].lower()

        if element_key not in atomic_numbers:
            raise ValueError(
                f"Unknown element symbol: '{element}'. "
                f"Cannot determine atomic number for RESP input. "
                f"Please check PDB file element column."
            )

        return atomic_numbers[element_key]

    def _get_backbone_restraints(self,
                                atoms_info: List[Dict],
                                chgmod: int,
                                ff_collector,
                                type_assignments: Optional[Dict] = None,
                                scaffolding_resids: Optional[set] = None) -> Dict[int, float]:
        """
        Get backbone atom restraints from force field or type_assignments.

        Args:
            atoms_info: List of atom info dicts (from PDB parsing)
            chgmod: Backbone restraint level (1-3)
            ff_collector: ForcefieldDataCollector from Step 1 (may be None)
            type_assignments: Dict of coord -> assignment dict (backup from prmtop)
            scaffolding_resids: (chain, resid) tuples that get a group-sum
                constraint instead of per-atom backbone restraints. Their
                atoms are skipped here — including cap residues like NME/ACE
                whose backbone-named atoms (N, C, O) have no library charge
                to look up, which would otherwise emit spurious warnings.

        Returns:
            Dict mapping atom_index (1-based) → charge
        """
        backbone_atom_names = self.BACKBONE_ATOMS[chgmod]
        scaffolding_resids = scaffolding_resids or set()
        restraints = {}

        # Build lookup for type_assignments by (resname, atom_name) -> charge
        # This is used when ff_collector is None (prmtop-based workflow)
        charge_lookup = {}
        if type_assignments:
            for coords, assignment in type_assignments.items():
                if isinstance(assignment, dict):
                    key = (assignment.get('resname', ''), assignment.get('atom_name', ''))
                    charge = assignment.get('charge')
                    if charge is not None:
                        charge_lookup[key] = charge

        for idx, atom_info in enumerate(atoms_info, 1):
            atom_name = atom_info['atom_name']
            resname = atom_info['resname']

            # Scaffolding residues (caps, standard protein) get a group-sum
            # constraint instead, so their backbone atoms are never restrained
            # individually. Skip them here to avoid a fruitless library lookup
            # (and the "Could not get charge for NME:N" warning it produces).
            if (atom_info['chain'], atom_info['resid']) in scaffolding_resids:
                continue

            if atom_name in backbone_atom_names:
                charge = None

                # Try ff_collector first (legacy path)
                if ff_collector is not None:
                    try:
                        atom_type_def = ff_collector.get_atom_type(resname, atom_name)
                        if atom_type_def and hasattr(atom_type_def, 'atom_charge') and atom_type_def.atom_charge is not None:
                            charge = atom_type_def.atom_charge
                    except Exception as e:
                        self.logger.debug(f"ff_collector lookup failed for {resname}:{atom_name}: {e}")

                # Fall back to type_assignments (prmtop-based workflow)
                if charge is None and charge_lookup:
                    charge = charge_lookup.get((resname, atom_name))

                if charge is not None:
                    restraints[idx] = charge
                else:
                    self.logger.warning(f"Could not get charge for backbone atom {resname}:{atom_name}")

        return restraints

    def _detect_equivalences(self, pdb_file: str, atoms_info: List[Dict]) -> List[List[int]]:
        """
        Auto-detect chemically equivalent atoms for stage 2.

        Detects:
        - Methyl hydrogens (H1, H2, H3 on same carbon)
        - Methylene hydrogens (HA, HB on same carbon)
        - Equivalent groups in ACE, NME, GLY capping residues

        Based on MCPB's get_equal_atoms() implementation.

        Args:
            pdb_file: Path to PDB file
            atoms_info: List of atom info dicts

        Returns:
            List of equivalent atom groups (1-based RESP-input indices, NOT PDB serials)
        """
        # Build connectivity from CONECT records (in PDB-serial numbering).
        raw_bonds = self._parse_conect_records(pdb_file)

        # CONECT records reference PDB serial numbers, but the rest of the
        # pipeline (resp.in column-2, equivalence groups returned here)
        # uses 1-based RESP-input indices that match atoms_info order.
        # Build a serial->RESP-index map and translate the bond list once.
        serial_to_resp = {info['serial']: idx
                          for idx, info in enumerate(atoms_info, 1)}
        bonds = []
        for s1, s2 in raw_bonds:
            r1 = serial_to_resp.get(s1)
            r2 = serial_to_resp.get(s2)
            if r1 is not None and r2 is not None:
                bonds.append((r1, r2))

        equivalences = []

        # Build residue-based atom groups
        residue_atoms = {}
        for idx, atom_info in enumerate(atoms_info, 1):
            key = (atom_info['chain'], atom_info['resid'], atom_info['resname'])
            if key not in residue_atoms:
                residue_atoms[key] = []
            residue_atoms[key].append((idx, atom_info))

        # For each residue, find equivalent atoms
        for (chain, resid, resname), res_atoms in residue_atoms.items():
            # Get bonds within this residue (RESP-input indices on both sides).
            res_atom_indices = {idx for idx, _ in res_atoms}
            res_bonds = [(a, b) for a, b in bonds if a in res_atom_indices and b in res_atom_indices]

            # Find CH2 and CH3 groups
            equiv_groups = self._find_ch_groups(res_atoms, res_bonds, atoms_info)
            equivalences.extend(equiv_groups)

            # Special handling for ACE, NME, GLY capping groups
            if resname in ['ACE', 'NME', 'GLY']:
                special_equiv = self._get_capping_equivalences(res_atoms, resname, atoms_info)
                equivalences.extend(special_equiv)

        return equivalences

    def _expand_cross_residue_equivalences(
        self,
        atoms_info: List[Dict],
        groups: List[List[Tuple[str, int]]],
    ) -> List[List[int]]:
        """
        Expand user-specified cross-residue groups into per-atom-name equivalence groups.

        For each group (list of (chain, resid) tuples), every atom name that
        appears in ALL group members becomes a single equivalence group
        containing the matching RESP-input index from each member.

        Args:
            atoms_info: Parsed PDB atoms (each dict has chain, resid, atom_name, ...).
            groups: User-selected residue groups, each as a list of (chain, resid)
                tuples. Same-resname is the caller's responsibility to enforce.

        Returns:
            List of equivalence groups (1-based RESP-input indices), one per
            atom name that's present in all members of a group.
        """
        # Build (chain, resid) -> {atom_name: resp_idx}
        residue_atom_map: Dict[Tuple[str, int], Dict[str, int]] = {}
        for idx, info in enumerate(atoms_info, 1):
            key = (info['chain'], info['resid'])
            residue_atom_map.setdefault(key, {})[info['atom_name']] = idx

        expanded: List[List[int]] = []
        for group in groups:
            if len(group) < 2:
                continue

            members = [residue_atom_map.get(res_key, {}) for res_key in group]
            if any(not m for m in members):
                self.logger.debug(
                    f"Cross-residue group {group} has missing residue(s); skipping"
                )
                continue

            # Atom names common to every group member
            common_names = set(members[0].keys())
            for m in members[1:]:
                common_names &= set(m.keys())

            for atom_name in sorted(common_names):
                indices = [m[atom_name] for m in members]
                if len(set(indices)) >= 2:
                    expanded.append(sorted(set(indices)))

        return expanded

    def _merge_equivalence_groups(self, groups: List[List[int]]) -> List[List[int]]:
        """
        Merge transitively-connected equivalence groups via union-find.

        Within-residue and cross-residue equivalences can overlap (e.g.,
        HB2≡HB3 within each Cys, and HB2 across Cys equivalent, and HB3
        across Cys equivalent → all 8 HB atoms ultimately one group).
        Without merging, duplicate column-2 anchor pointers would conflict
        in the resp.in file. With merging, each connected component
        collapses to one anchor.
        """
        parent: Dict[int, int] = {}

        def find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(a: int, b: int) -> None:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb

        for group in groups:
            for atom in group:
                parent.setdefault(atom, atom)
            for atom in group[1:]:
                union(group[0], atom)

        components: Dict[int, List[int]] = {}
        for atom in parent:
            root = find(atom)
            components.setdefault(root, []).append(atom)

        return [sorted(component) for component in components.values()
                if len(component) >= 2]

    def _parse_conect_records(self, pdb_file: str) -> List[Tuple[int, int]]:
        """
        Parse CONECT records from PDB file to get bond information.

        Args:
            pdb_file: Path to PDB file

        Returns:
            List of (atom1_idx, atom2_idx) tuples (1-based serial numbers)
        """
        bonds = []
        with open(pdb_file) as f:
            for line in f:
                if line.startswith('CONECT'):
                    parts = line[6:].split()
                    if len(parts) >= 2:
                        try:
                            atom1 = int(parts[0])
                            for i in range(1, len(parts)):
                                atom2 = int(parts[i])
                                # Store bond once (avoid duplicates)
                                bond = tuple(sorted([atom1, atom2]))
                                if bond not in bonds:
                                    bonds.append(bond)
                        except ValueError:
                            pass
        return bonds

    def _find_ch_groups(self, res_atoms: List[Tuple], res_bonds: List[Tuple],
                       atoms_info: List[Dict]) -> List[List[int]]:
        """
        Find CH2 and CH3 groups in a residue.

        Args:
            res_atoms: List of (index, atom_info) tuples for residue
            res_bonds: List of bonds within residue
            atoms_info: Full atom info list (0-indexed for lookup)

        Returns:
            List of equivalent atom groups
        """
        equiv_groups = []
        idx_to_info = {idx: info for idx, info in res_atoms}

        # Find all hydrogens bonded to carbons
        for carbon_idx, carbon_info in res_atoms:
            if carbon_info['element'] != 'C':
                continue

            # Find hydrogens bonded to this carbon
            bonded_h = []
            for bond in res_bonds:
                if carbon_idx in bond:
                    other_idx = bond[1] if bond[0] == carbon_idx else bond[0]
                    if other_idx in idx_to_info:
                        other_info = idx_to_info[other_idx]
                        if other_info['element'] == 'H':
                            bonded_h.append(other_idx)

            # If 2 or more hydrogens on same carbon, they're equivalent
            if len(bonded_h) >= 2:
                equiv_groups.append(sorted(bonded_h))

        return equiv_groups

    def _get_capping_equivalences(self, res_atoms: List[Tuple], resname: str,
                                  atoms_info: List[Dict]) -> List[List[int]]:
        """
        Get special equivalences for capping groups (ACE, NME, GLY).

        Based on MCPB's get_equ_for_ang() implementation.

        Args:
            res_atoms: List of (index, atom_info) tuples for residue
            resname: Residue name
            atoms_info: Full atom info list

        Returns:
            List of equivalent atom groups
        """
        equiv_groups = []
        name_to_idx = {info['atom_name']: idx for idx, info in res_atoms}

        if resname == 'ACE':
            # HH31, HH32, HH33 are equivalent
            h_atoms = []
            for h_name in ['HH31', 'HH32', 'HH33']:
                if h_name in name_to_idx:
                    h_atoms.append(name_to_idx[h_name])
            if len(h_atoms) >= 2:
                equiv_groups.append(sorted(h_atoms))

        elif resname == 'NME':
            # HH31, HH32, HH33 are equivalent
            h_atoms = []
            for h_name in ['HH31', 'HH32', 'HH33']:
                if h_name in name_to_idx:
                    h_atoms.append(name_to_idx[h_name])
            if len(h_atoms) >= 2:
                equiv_groups.append(sorted(h_atoms))

        elif resname == 'GLY':
            # HA2, HA3 are equivalent
            if 'HA2' in name_to_idx and 'HA3' in name_to_idx:
                equiv_groups.append(sorted([name_to_idx['HA2'], name_to_idx['HA3']]))

        return equiv_groups

    def _write_respin1(self,
                      output_path: Path,
                      atoms_info: List[Dict],
                      total_charge: int,
                      restraints: Dict[int, float],
                      residue_groups: List[List[int]] = None) -> None:
        """
        Write stage 1 resp input (weak restraint).

        Format (based on MCPB resp_fitting.py lines 339-387):
        &cntrl
         nmol = 1,
         ihfree = 1,
         ioutopt = 1,
         qwt = 0.00050,
        &end
            1.0
        Resp charges for organic molecule
         {total_charge}  {n_atoms}
        {atomic_number}  0
        ...
        [restraint lines]

        Args:
            output_path: Path for resp1.in file
            atoms_info: List of atom info dicts
            total_charge: Total molecular charge
            restraints: Dict of atom_index → fixed_charge
        """
        n_atoms = len(atoms_info)

        with open(output_path, 'w') as f:
            # Header
            f.write("Resp charges for organic molecule\n")
            f.write(" \n")
            f.write(" &cntrl\n")
            f.write(" \n")
            f.write(" nmol = 1,\n")
            f.write(" ihfree = 1,\n")
            f.write(" ioutopt = 1,\n")
            f.write(" qwt = 0.00050,\n")
            f.write(" \n")
            f.write(" &end\n")
            f.write("    1.0\n")
            f.write("Resp charges for organic molecule\n")
            f.write(f"{total_charge:5d} {n_atoms:4d}\n")

            # Atom section (atomic numbers, all free in stage 1)
            for atom_info in atoms_info:
                atomic_num = self._get_atomic_number(atom_info['element'])
                f.write(f"{atomic_num:5d}    0\n")

            # Constraints section (no blank line before — blank line terminates)
            # Part 3: Residue group constraints (capping groups sum to 0)
            if residue_groups:
                for group_atoms in residue_groups:
                    n_grp = len(group_atoms)
                    f.write(f"{n_grp:5d}{'0.00000':>10s}\n")
                    for idx, atom_idx in enumerate(group_atoms):
                        f.write(f"{1:5d}{atom_idx:5d}")
                        if (idx + 1) % 8 == 0:
                            f.write("\n")
                    if len(group_atoms) % 8 != 0:
                        f.write("\n")

            # Part 4: Backbone charge restraints (metal-site residues only)
            if restraints:
                for atom_idx, charge in sorted(restraints.items()):
                    f.write(f"{1:5d}{charge:10.5f}\n")
                    f.write(f"{1:5d}{atom_idx:5d}\n")

            # Termination (blank lines to end sections 7 and 8)
            f.write("\n\n\n\n")

    def _write_respin2(self,
                      output_path: Path,
                      atoms_info: List[Dict],
                      total_charge: int,
                      restraints: Dict[int, float],
                      equivalences: List[List[int]],
                      residue_groups: List[List[int]] = None,
                      frozen_cap_atoms: set = None) -> None:
        """
        Write stage 2 resp input (strong restraint + equivalences).

        Args:
            output_path: Path for resp2.in file
            atoms_info: List of atom info dicts
            total_charge: Total molecular charge
            restraints: Dict of atom_index → fixed_charge
            equivalences: List of atom groups to constrain equal charges
            residue_groups: List of atom index groups for capping residues (sum to 0)
            frozen_cap_atoms: Set of 1-based atom indices to freeze (-99) in capping groups
        """
        n_atoms = len(atoms_info)

        with open(output_path, 'w') as f:
            # Header
            f.write("Resp charges for organic molecule\n")
            f.write(" \n")
            f.write(" &cntrl\n")
            f.write(" \n")
            f.write(" nmol = 1,\n")
            f.write(" ihfree = 1,\n")
            f.write(" ioutopt = 1,\n")
            f.write(" iqopt = 2,\n")  # Use charges from stage 1
            f.write(" qwt = 0.00100,\n")  # Stronger restraint
            f.write(" \n")
            f.write(" &end\n")
            f.write("    1.0\n")
            f.write("Resp charges for organic molecule\n")
            f.write(f"{total_charge:5d} {n_atoms:4d}\n")

            # Build equivalence map (atom_idx → reference_idx)
            equiv_map = {}
            for group in equivalences:
                if len(group) > 1:
                    ref_idx = group[0]  # First atom is reference
                    for atom_idx in group[1:]:
                        equiv_map[atom_idx] = ref_idx

            # Build restraint fixed map (atoms with restraints get -99)
            restraint_atoms = set(restraints.keys())

            # Atom section with equivalences
            frozen_set = frozen_cap_atoms or set()

            for idx in range(1, n_atoms + 1):
                atomic_num = self._get_atomic_number(atoms_info[idx-1]['element'])

                # Determine equivalence flag
                if idx in restraint_atoms:
                    # Restrained backbone atoms are fixed (-99)
                    equiv_flag = -99
                elif idx in frozen_set:
                    # Frozen capping group atoms (ACE: C,O; NME: N,H) — MCPB convention
                    equiv_flag = -99
                elif idx in equiv_map:
                    # Equivalent to another atom
                    equiv_flag = equiv_map[idx]
                else:
                    # Free atom (0)
                    equiv_flag = 0

                f.write(f"{atomic_num:5d}{equiv_flag:5d}\n")

            # Constraints section (no blank line before — blank line terminates)
            # Part 3: Residue group constraints (capping groups sum to 0)
            if residue_groups:
                for group_atoms in residue_groups:
                    n_grp = len(group_atoms)
                    f.write(f"{n_grp:5d}{'0.00000':>10s}\n")
                    for idx, atom_idx in enumerate(group_atoms):
                        f.write(f"{1:5d}{atom_idx:5d}")
                        if (idx + 1) % 8 == 0:
                            f.write("\n")
                    if len(group_atoms) % 8 != 0:
                        f.write("\n")

            # Part 4: Backbone charge restraints (metal-site residues only)
            if restraints:
                for atom_idx, charge in sorted(restraints.items()):
                    f.write(f"{1:5d}{charge:10.5f}\n")
                    f.write(f"{1:5d}{atom_idx:5d}\n")

            # Termination (blank lines to end sections 7 and 8)
            f.write("\n\n\n\n")
