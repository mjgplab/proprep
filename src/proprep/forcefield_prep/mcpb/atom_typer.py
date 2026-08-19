"""
MCPB Atom Typer

Assigns AMBER atom types to RedoxSite atoms following MCPB conventions.
Implements systematic 2-character atom type renaming with unlimited capacity.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional, TYPE_CHECKING
from enum import Enum
import logging
import math

from proprep.structure_prep.comprehensive_redox_detector import RedoxSite, RedoxSiteAtom, METALS
from .terminal_classifier import TerminalClassifier, TerminalType


# Amber names a hydrogen for the atom it is bonded to, not for hydrogen
# generally. 'H' is specifically the amide/amine hydrogen and carries
# r* = 0.6000 / eps = 0.0157; 'HO', the hydroxyl hydrogen, is 0.0000 / 0.0000
# by convention. Identical in parm10.dat and parm19.dat, so this is not tied to
# one protein force field, and GAFF mirrors it lowercase (ho/hn/hs/hc).
#
#   H   1.008  0.161  H bonded to nitrogen atoms
#   HO  1.008  0.135  hydroxyl group
#   HS  1.008  0.135  hydrogen bonded to sulphur
#
# Typing an inferred hydrogen from its element alone gave a metal-bound hydroxo
# proton a 0.6 A vdW sphere that the hydroxyl convention deliberately removes.
H_TYPE_BY_NEIGHBOR = {'O': 'HO', 'S': 'HS', 'N': 'H', 'C': 'HC'}


def hydrogen_type_from_neighbors(coords, neighbors, default: str = 'H') -> str:
    """Amber type for a hydrogen, from the heavy atom it is bonded to.

    Args:
        coords: the hydrogen's (x, y, z)
        neighbors: iterable of ``(coords, element)`` for OTHER atoms in the same
            residue. Same-residue filtering is the caller's job; hydrogens in it
            are ignored here.
        default: type when no bonded heavy atom is found.

    Bond test is the covalent-radius one the cluster bond perception uses, so
    the typer and the bond perceiver cannot disagree about what is bonded.
    """
    from proprep.forcefield_prep.mcpb.mol2_writer import Mol2Writer

    radii, fallback = Mol2Writer._COVALENT_RADII, Mol2Writer._COVALENT_DEFAULT
    r_h = radii.get('H', fallback)

    best, best_distance = None, None
    for other_coords, element in neighbors:
        element = (element or '').strip().upper()
        if element in ('H', ''):
            continue
        if other_coords == coords:
            continue
        cutoff = r_h + radii.get(element, fallback) + Mol2Writer._COVALENT_TOLERANCE
        distance = math.dist(coords, other_coords)
        if distance < cutoff and (best_distance is None or distance < best_distance):
            best, best_distance = element, distance

    return H_TYPE_BY_NEIGHBOR.get(best, default)

if TYPE_CHECKING:
    from proprep.forcefield_prep.forcefield_data import ForceFieldData


class AtomTypeRenamingMode(Enum):
    """MCPB atom type renaming modes"""
    NO_RENAME = 0      # Mode 1n: no renaming
    METAL_ONLY = 1     # Mode 1m: only rename metals
    FULL_AUTO = 2      # Mode 1a: rename metals + ligands (default)


@dataclass
class AtomTypeAssignment:
    """Result of atom type assignment for one atom"""
    # Atom identification
    coords: Tuple[float, float, float]
    chain: str
    resname: str
    resid: int
    atom_name: str
    element: str

    # Atom type progression
    original_type: str       # From force field library
    renamed_type: str        # After MCPB renaming

    # Classification
    is_center: bool          # Is this a redox center?
    is_metal_ligand: bool    # Is this bonded to metal?
    terminal_type: TerminalType

    # Library source
    library_source: str      # "ff19SB", "custom_mol2", "inferred"
    charge: Optional[float] = None


class MCPBAtomTyper:
    """
    Assigns AMBER atom types to RedoxSite atoms following MCPB conventions.

    Features:
    - Multi-method terminal residue detection
    - Force field library lookup
    - Systematic 2-character atom type renaming (unlimited capacity)
    - User-configurable naming schemes
    """

    def __init__(self,
                 ff_data: 'ForceFieldData',
                 custom_mol2_files: Optional[List[str]] = None,
                 renaming_mode: AtomTypeRenamingMode = AtomTypeRenamingMode.FULL_AUTO,
                 logger=None):
        """
        Initialize atom typer.

        Args:
            ff_data: ForceFieldData with atom type definitions loaded
            custom_mol2_files: Additional mol2 files for non-standard residues
            renaming_mode: How to rename metal site atoms
            logger: Optional logger instance
        """
        self.ff_data = ff_data
        self.renaming_mode = renaming_mode
        self.logger = logger or logging.getLogger(__name__)

        # Load custom mol2 files into ff_data
        if custom_mol2_files:
            for mol2_file in custom_mol2_files:
                self.ff_data.load_mol2_file(Path(mol2_file))

        # Atom type assignments (keyed by coordinates)
        self.type_assignments: Dict[Tuple[float, float, float], AtomTypeAssignment] = {}

        # Metal-ligand tracking for renaming
        self.metal_atoms: List[Tuple[float, float, float]] = []
        self.ligand_atoms: List[Tuple[float, float, float]] = []

        # Terminal classifications
        self.terminal_map: Dict[Tuple[str, int, str], TerminalType] = {}

    def assign_atom_types(self,
                         site: RedoxSite,
                         pdb_file: str,
                         require_terminal_confirmation: bool = True,
                         user_config: Optional[Dict] = None) -> Dict[Tuple, AtomTypeAssignment]:
        """
        Assign atom types to all atoms in a RedoxSite.

        Args:
            site: RedoxSite object to process
            pdb_file: Path to original PDB file
            require_terminal_confirmation: Prompt user to confirm terminal residues
            user_config: Optional user configuration for renaming (see _apply_renaming)

        Returns:
            Dictionary mapping coordinates to type assignments
        """
        self.logger.info("Starting atom type assignment...")
        self.type_assignments.clear()

        # Step 1: Classify residue terminals
        self.logger.info("Step 1: Classifying terminal residues...")
        terminal_classifier = TerminalClassifier(site, pdb_file, logger=self.logger)
        self.terminal_map = terminal_classifier.classify_all_residues(
            require_confirmation=require_terminal_confirmation
        )

        # Step 2: Identify metals and ligands
        self.logger.info("Step 2: Identifying metal centers and ligands...")
        self._identify_metal_ligand_atoms(site)
        self.logger.info(f"  Found {len(self.metal_atoms)} metal atoms")
        self.logger.info(f"  Found {len(self.ligand_atoms)} ligand atoms")

        # Step 3: Assign original types from force field
        self.logger.info("Step 3: Assigning original atom types from force field...")
        for atom in site.atoms:
            assignment = self._assign_atom_type(atom, site)
            self.type_assignments[atom.coords] = assignment

        # Count inferred types (potential issues)
        inferred_count = sum(1 for a in self.type_assignments.values()
                            if a.library_source == 'inferred')
        if inferred_count > 0:
            self.logger.warning(
                f"  {inferred_count} atoms have inferred types "
                f"(consider providing custom mol2 files)"
            )

        # Step 4: Apply renaming based on mode
        if self.renaming_mode != AtomTypeRenamingMode.NO_RENAME:
            self.logger.info(f"Step 4: Applying atom type renaming (mode: {self.renaming_mode.name})...")
            self._apply_renaming(user_config)

        self.logger.info(f"Atom type assignment complete: {len(self.type_assignments)} atoms")
        return self.type_assignments

    def _identify_metal_ligand_atoms(self, site: RedoxSite):
        """
        Identify which atoms are metals and which are metal ligands.
        Uses site.centers and site.bonds information.

        Metals come from two center types:
        - METAL_ION: the center itself is the metal (e.g., isolated Zn)
        - ORGANOMETALLIC_COFACTOR: contains embedded metal atoms (e.g., Fe in heme)
        """
        self.metal_atoms.clear()
        self.ligand_atoms.clear()

        # Get metal atom coordinates from centers
        for center in site.centers:
            if center.center_type.value == 'metal_ion':
                # Isolated metal - center coords ARE the metal
                self.metal_atoms.append(center.coords)
            elif center.center_type.value == 'organometallic_cofactor':
                # Embedded metal - find metal atom(s) in this residue
                for atom in site.atoms:
                    if (atom.chain == center.chain and
                        atom.resid == center.resid and
                        atom.resname == center.resname and
                        atom.element.upper() in METALS):
                        if atom.coords not in self.metal_atoms:
                            self.metal_atoms.append(atom.coords)

        # Get ligand atoms (non-metal atoms bonded to metals)
        for bond in site.bonds:
            if bond.chemical_type == 'coordinate':
                # One atom is metal, other is ligand
                atom1_metal = bond.atom1_coords in self.metal_atoms
                atom2_metal = bond.atom2_coords in self.metal_atoms

                if atom1_metal and not atom2_metal:
                    if bond.atom2_coords not in self.ligand_atoms:
                        self.ligand_atoms.append(bond.atom2_coords)
                elif atom2_metal and not atom1_metal:
                    if bond.atom1_coords not in self.ligand_atoms:
                        self.ligand_atoms.append(bond.atom1_coords)

    def _assign_atom_type(self, atom: RedoxSiteAtom, site: RedoxSite) -> AtomTypeAssignment:
        """
        Assign original atom type from force field library.

        Args:
            atom: RedoxSiteAtom to assign type to
            site: Parent RedoxSite

        Returns:
            AtomTypeAssignment with original type and metadata
        """
        res_key = (atom.chain, atom.resid, atom.insertion_code)
        terminal_type = self.terminal_map.get(res_key, TerminalType.INTERNAL)

        # Build residue name with terminal prefix if needed
        if terminal_type == TerminalType.NTERM:
            ff_resname = f"N{atom.resname}"
        elif terminal_type == TerminalType.CTERM:
            ff_resname = f"C{atom.resname}"
        else:
            ff_resname = atom.resname

        # Lookup in force field data
        original_type = self.ff_data.get_atom_type(ff_resname, atom.atom_name)

        if original_type:
            charge = self.ff_data.get_atom_charge(ff_resname, atom.atom_name)
            library_source = self.ff_data.get_source_file(ff_resname, atom.atom_name) or "force_field"
        else:
            # Fallback: infer from element
            original_type = self._infer_type_from_element(atom.element, atom, site)
            charge = None
            library_source = "inferred"
            self.logger.debug(
                f"No library entry for {ff_resname}-{atom.atom_name}, inferred type: {original_type}"
            )

        # Check classifications
        is_center = atom.coords in self.metal_atoms
        is_metal_ligand = atom.coords in self.ligand_atoms

        return AtomTypeAssignment(
            coords=atom.coords,
            chain=atom.chain,
            resname=atom.resname,
            resid=atom.resid,
            atom_name=atom.atom_name,
            element=atom.element,
            original_type=original_type,
            renamed_type=original_type,  # Will be updated in renaming step
            is_center=is_center,
            is_metal_ligand=is_metal_ligand,
            terminal_type=terminal_type,
            library_source=library_source,
            charge=charge
        )

    # Amber names a hydrogen for the atom it is bonded to, not for hydrogen
    # generally. 'H' is specifically the amide/amine hydrogen and carries
    # r* = 0.6000 / eps = 0.0157; 'HO', the hydroxyl hydrogen, is 0.0000 /
    # 0.0000 by convention. Identical in parm10.dat and parm19.dat, so this is
    # not tied to one protein force field.
    #
    #   H   1.008  0.161  H bonded to nitrogen atoms
    #   HO  1.008  0.135  hydroxyl group
    #   HS  1.008  0.135  hydrogen bonded to sulphur
    #
    # Typing every inferred hydrogen 'H' gave a metal-bound hydroxo proton a
    # 0.6 A vdW sphere that the hydroxyl convention deliberately removes.
    _H_TYPE_BY_NEIGHBOR = {'O': 'HO', 'S': 'HS', 'N': 'H ', 'C': 'HC'}

    def _bonded_heavy_neighbor(self, atom, site) -> Optional[str]:
        """Element of the heavy atom this one is bonded to, within its residue.

        Same covalent-radius test the cluster bond perception uses, so the two
        cannot disagree about what is bonded to what.
        """
        from proprep.forcefield_prep.mcpb.mol2_writer import Mol2Writer

        if site is None or atom is None:
            return None

        radii, default = Mol2Writer._COVALENT_RADII, Mol2Writer._COVALENT_DEFAULT
        r_h = radii.get((atom.element or '').strip().upper(), default)

        best, best_distance = None, None
        for other in getattr(site, 'atoms', None) or []:
            if other is atom or other.coords == atom.coords:
                continue
            if (other.chain, other.resid) != (atom.chain, atom.resid):
                continue
            element = (other.element or '').strip().upper()
            if element in ('H', ''):
                continue
            cutoff = r_h + radii.get(element, default) + Mol2Writer._COVALENT_TOLERANCE
            distance = math.dist(atom.coords, other.coords)
            if distance < cutoff and (best_distance is None or distance < best_distance):
                best, best_distance = element, distance
        return best

    def _infer_type_from_element(self, element: str, atom=None, site=None) -> str:
        """
        Fallback: infer atom type from element symbol.

        A hydrogen is typed from the atom it is bonded to when that can be
        determined; see _H_TYPE_BY_NEIGHBOR. Everything else is by element.

        Args:
            element: Element symbol (e.g., 'C', 'N', 'Zn')
            atom: The atom being typed, for hydrogen neighbour lookup
            site: Parent RedoxSite, for hydrogen neighbour lookup

        Returns:
            Generic 2-character atom type
        """
        if (element or '').strip().upper() == 'H' and atom is not None and site is not None:
            same_residue = [
                (a.coords, a.element) for a in (getattr(site, 'atoms', None) or [])
                if (a.chain, a.resid) == (atom.chain, atom.resid)
            ]
            inferred = hydrogen_type_from_neighbors(atom.coords, same_residue)
            return inferred if len(inferred) == 2 else inferred + ' '
        # Pad to 2 characters (AMBER requirement)
        element_types = {
            'C': 'CT', 'N': 'N ', 'O': 'O ', 'S': 'S ',
            'H': 'H ', 'P': 'P ', 'F': 'F ', 'CL': 'Cl',
            'BR': 'Br', 'I': 'I ',
            # Metals
            'ZN': 'Zn', 'FE': 'Fe', 'CU': 'Cu', 'MG': 'Mg',
            'CA': 'Ca', 'MN': 'Mn', 'NI': 'Ni', 'CO': 'Co',
            'NA': 'Na', 'K': 'K ',
        }
        inferred = element_types.get(element.upper(), 'XX')
        # Ensure 2 characters
        if len(inferred) == 1:
            inferred = inferred + ' '
        return inferred

    def _apply_renaming(self, user_config: Optional[Dict] = None):
        """
        Apply systematic renaming to atom type assignments.

        Args:
            user_config: Optional user configuration:
                {
                    'metal_prefix': 'M',     # Custom prefix for metals
                    'ligand_prefix': 'L',    # Custom prefix for ligands
                    'rename_metals': True,   # Rename metals?
                    'rename_ligands': True,  # Rename ligands?
                    'custom_names': {        # Manual overrides
                        (x, y, z): 'ZN',
                        (x2, y2, z2): 'S1'
                    }
                }
        """
        # Parse user config
        config = user_config or {}
        metal_prefix = config.get('metal_prefix', 'M')
        ligand_prefix = config.get('ligand_prefix', 'L')
        rename_metals = config.get('rename_metals', True)
        rename_ligands = config.get('rename_ligands', True)
        custom_names = config.get('custom_names', {})

        # Validate prefixes are 1 character
        if len(metal_prefix) != 1 or len(ligand_prefix) != 1:
            raise ValueError("Prefixes must be exactly 1 character for 2-char atom types")

        if self.renaming_mode == AtomTypeRenamingMode.NO_RENAME:
            # No renaming at all
            return

        elif self.renaming_mode == AtomTypeRenamingMode.METAL_ONLY:
            # Rename only metals with charge-based names
            if rename_metals:
                for metal_coord in self.metal_atoms:
                    if metal_coord in custom_names:
                        # Use custom name
                        self.type_assignments[metal_coord].renamed_type = custom_names[metal_coord]
                    else:
                        # Use element + charge format (e.g., "Zn2+")
                        assignment = self.type_assignments[metal_coord]
                        element = assignment.element
                        charge = self._get_metal_charge(assignment)
                        if charge > 0:
                            assignment.renamed_type = f"{element}{charge}+"
                        elif charge < 0:
                            assignment.renamed_type = f"{element}{abs(charge)}-"

        elif self.renaming_mode == AtomTypeRenamingMode.FULL_AUTO:
            # Systematic renaming for metals and ligands
            if rename_metals:
                for i, metal_coord in enumerate(self.metal_atoms):
                    if metal_coord in custom_names:
                        self.type_assignments[metal_coord].renamed_type = custom_names[metal_coord]
                    else:
                        # Generate systematic 2-char name
                        self.type_assignments[metal_coord].renamed_type = \
                            self._generate_systematic_name(metal_prefix, i)

            if rename_ligands:
                for i, ligand_coord in enumerate(self.ligand_atoms):
                    if ligand_coord in custom_names:
                        self.type_assignments[ligand_coord].renamed_type = custom_names[ligand_coord]
                    else:
                        # Generate systematic 2-char name
                        self.type_assignments[ligand_coord].renamed_type = \
                            self._generate_systematic_name(ligand_prefix, i)

    def _generate_systematic_name(self, prefix: str, index: int) -> str:
        """
        Generate systematic 2-character atom type name.

        Available names per prefix (2-char limit):
        - Digits: M1-M9 (9 names)
        - Letters: MA-MZ (26 names)
        Total: 35 names per prefix

        For additional capacity, can use X, Y, Z prefixes:
        - M: M1-M9, MA-MZ (35)
        - X: X1-X9, XA-XZ (35)
        - Y: Y1-Y9, YA-YZ (35)
        - Z: Z1-Z9, ZA-ZZ (35)
        Total: 140 unique 2-character names

        Args:
            prefix: Single character prefix ('M', 'L', 'X', 'Y', 'Z')
            index: 0-based index

        Returns:
            2-character atom type name

        Examples:
            >>> _generate_systematic_name('M', 0)
            'M1'
            >>> _generate_systematic_name('M', 8)
            'M9'
            >>> _generate_systematic_name('M', 9)
            'MA'
            >>> _generate_systematic_name('M', 34)
            'MZ'
            >>> _generate_systematic_name('L', 15)
            'LG'
        """
        if index < 9:
            # First 9: use digits 1-9
            return f"{prefix}{index + 1}"
        elif index < 35:
            # Next 26: use letters A-Z
            letter_index = index - 9
            letter = chr(ord('A') + letter_index)
            return f"{prefix}{letter}"
        else:
            raise ValueError(
                f"Index {index} exceeds 2-character naming capacity (0-34) for prefix '{prefix}'. "
                f"Consider using additional prefixes (X, Y, Z) or switching to 3-character names."
            )

    def _get_metal_charge(self, assignment: AtomTypeAssignment) -> int:
        """
        Get metal ion charge.

        Args:
            assignment: AtomTypeAssignment for metal atom

        Returns:
            Integer charge (e.g., 2 for Zn2+)
        """
        if assignment.charge is not None:
            return int(round(assignment.charge))

        # If charge not available, ASK USER (no hard-coding!)
        print(f"\nCharge unknown for {assignment.element} at {assignment.resname}:{assignment.resid}")
        while True:
            try:
                charge_str = input(f"  Enter charge for {assignment.element} (e.g., 2 for 2+): ")
                charge = int(charge_str)
                return charge
            except ValueError:
                print("  Invalid input. Please enter an integer charge.")

    def get_validation_report(self) -> Dict:
        """
        Generate validation report for assigned atom types.

        Returns:
            Dictionary with validation statistics and warnings
        """
        if not self.type_assignments:
            return {'error': 'No atom types assigned yet'}

        total_atoms = len(self.type_assignments)
        inferred_atoms = [a for a in self.type_assignments.values()
                         if a.library_source == 'inferred']
        metal_atoms = [a for a in self.type_assignments.values() if a.is_center]
        ligand_atoms = [a for a in self.type_assignments.values() if a.is_metal_ligand]
        renamed_atoms = [a for a in self.type_assignments.values()
                        if a.original_type != a.renamed_type]

        report = {
            'total_atoms': total_atoms,
            'metal_atoms': len(metal_atoms),
            'ligand_atoms': len(ligand_atoms),
            'renamed_atoms': len(renamed_atoms),
            'inferred_atoms': len(inferred_atoms),
            'loaded_leaprcs': list(self.ff_data.loaded_leaprcs),
            'renaming_mode': self.renaming_mode.name,
            'warnings': []
        }

        # Generate warnings
        if inferred_atoms:
            report['warnings'].append(
                f"{len(inferred_atoms)} atoms have inferred types (no library entry). "
                f"Consider providing custom mol2 files."
            )
            report['inferred_details'] = [
                f"{a.resname}-{a.atom_name}: {a.element} → {a.original_type}"
                for a in inferred_atoms[:10]  # Limit to first 10
            ]

        return report
