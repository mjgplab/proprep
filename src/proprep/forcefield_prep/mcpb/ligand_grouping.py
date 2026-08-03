"""
Ligand Grouping Interface

Interactive table-based interface for users to define ligand grouping
for antechamber parameterization. No hard-coded assumptions about
specific clusters - completely user-controlled.

© 2024 ProPrep Developer. All rights reserved.
"""

from typing import List, Dict, Tuple, Set
from dataclasses import dataclass
from collections import defaultdict
import logging

# Import RedoxSite types - adjust path as needed
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from structure_prep.comprehensive_redox_detector import RedoxSite, RedoxSiteAtom, RedoxSiteBond
from proprep.utils.prompts import prompt_with_context


@dataclass
class LigandGroup:
    """A group of ligand atoms to be parameterized together"""
    group_id: str                    # e.g., "ligand_1"
    group_number: int                # User-specified number (1, 2, 3, ...)
    atoms: List[RedoxSiteAtom]       # Atoms in this group
    formal_charge: int = 0           # Net charge for antechamber
    parameterization_method: str = 'antechamber'  # 'antechamber', 'existing', or 'manual'
    existing_mol2: str = None        # Path to existing mol2 file (if method='existing')
    existing_frcmod: str = None      # Path to existing frcmod file (if method='existing')

    def __repr__(self):
        atom_summary = f"{len(self.atoms)} atoms"
        if self.atoms:
            elements = [a.element for a in self.atoms]
            elem_counts = {}
            for e in elements:
                elem_counts[e] = elem_counts.get(e, 0) + 1
            elem_str = ", ".join(f"{count}{elem}" for elem, count in sorted(elem_counts.items()))
            atom_summary = f"{len(self.atoms)} atoms ({elem_str})"
        return f"LigandGroup(id={self.group_id}, method={self.parameterization_method}, {atom_summary})"


class LigandGroupingInterface:
    """
    Table-based interface for user to control ligand grouping.

    Features:
    - Displays coordinating atoms grouped by residue
    - User-friendly syntax: "1:HEM 2:46 3:47-50"
    - No hard-coded assumptions about cluster types
    - Prompt for formal charge per group
    """

    def __init__(self, console=None, logger=None, processor=None):
        """
        Initialize ligand grouping interface.

        Args:
            console: Rich console for formatted output (optional)
            logger: Logger instance (optional)
            processor: PDBProcessor for session-recording context (optional)
        """
        self.console = console
        self.logger = logger or logging.getLogger(__name__)
        self.processor = processor

    def identify_ligand_atoms(self, redox_site: RedoxSite) -> List[RedoxSiteAtom]:
        """
        Identify all atoms that are ligands (bonded to metals via coordinate bonds).

        Uses RedoxSite.bonds to find atoms bonded to metal centers with
        chemical_type='coordinate'.

        Args:
            redox_site: RedoxSite object

        Returns:
            List of ligand atoms (RedoxSiteAtom objects)
        """
        # Get metal center coordinates
        metal_coords = set(center.coords for center in redox_site.centers)

        ligand_coords = set()

        # Find all coordinate bonds
        for bond in redox_site.bonds:
            if bond.chemical_type == 'coordinate':
                # One end is metal, other is ligand
                atom1_is_metal = bond.atom1_coords in metal_coords
                atom2_is_metal = bond.atom2_coords in metal_coords

                if atom1_is_metal and not atom2_is_metal:
                    ligand_coords.add(bond.atom2_coords)
                elif atom2_is_metal and not atom1_is_metal:
                    ligand_coords.add(bond.atom1_coords)
                # If both are metals, skip (metal-metal bond)

        # Get RedoxSiteAtom objects for ligand coordinates
        ligand_atoms = [
            atom for atom in redox_site.atoms
            if atom.coords in ligand_coords
        ]

        return ligand_atoms

    def group_by_residue(self, atoms: List[RedoxSiteAtom]) -> Dict[Tuple[str, str, int], List[RedoxSiteAtom]]:
        """
        Group atoms by residue for organized display.

        Args:
            atoms: List of RedoxSiteAtom objects

        Returns:
            Dict mapping (chain, resname, resid) -> list of atoms
        """
        by_residue = defaultdict(list)

        for atom in atoms:
            key = (atom.chain, atom.resname, atom.resid)
            by_residue[key].append(atom)

        # Sort atoms within each residue by atom name
        for key in by_residue:
            by_residue[key].sort(key=lambda a: a.atom_name)

        return dict(by_residue)

    def display_coordinating_atoms(self, redox_site: RedoxSite) -> Dict[int, RedoxSiteAtom]:
        """
        Display table of coordinating atoms grouped by residue.

        Args:
            redox_site: RedoxSite object

        Returns:
            Mapping of index -> atom for parsing user input
        """
        ligand_atoms = self.identify_ligand_atoms(redox_site)

        if not ligand_atoms:
            self._print("[yellow]No coordinating ligand atoms detected.[/yellow]")
            return {}

        return self.display_atoms(ligand_atoms)

    def display_atoms(self, atoms: List[RedoxSiteAtom]) -> Dict[int, RedoxSiteAtom]:
        """
        Display table of specific atoms grouped by residue.

        Args:
            atoms: List of RedoxSiteAtom objects to display

        Returns:
            Mapping of index -> atom for parsing user input
        """
        if not atoms:
            self._print("[yellow]No atoms to display.[/yellow]")
            return {}

        by_residue = self.group_by_residue(atoms)

        self._print("\n" + "[cyan]" + "═" * 60 + "[/cyan]")
        self._print("[bold cyan]Coordinating Atoms Detected:[/bold cyan]")
        self._print("[cyan]" + "═" * 60 + "[/cyan]")

        idx_to_atom = {}
        current_idx = 1

        for (chain, resname, resid), atoms in sorted(by_residue.items()):
            self._print(f"\n[yellow]Residue: {resname} (Chain {chain}, Resid {resid})[/yellow]")
            self._print(f"  {'Index':<8} {'Atom':<8} {'Element':<10} Coords")
            self._print("  " + "─" * 58)

            for atom in atoms:
                coords_str = f"({atom.coords[0]:.1f}, {atom.coords[1]:.1f}, {atom.coords[2]:.1f})"
                self._print(f"  {current_idx:<8} {atom.atom_name:<8} {atom.element:<10} {coords_str}")
                idx_to_atom[current_idx] = atom
                current_idx += 1

        self._print("[cyan]" + "═" * 60 + "[/cyan]\n")

        return idx_to_atom

    def prompt_user_grouping(self, idx_to_atom: Dict[int, RedoxSiteAtom],
                           redox_site: RedoxSite) -> List[LigandGroup]:
        """
        Interactive prompt for ligand grouping.

        Displays syntax help and prompts user to define groups using
        the syntax: "1:HEM 2:46 3:47-50,52"

        Args:
            idx_to_atom: Mapping of index -> atom from display_coordinating_atoms()
            redox_site: RedoxSite object

        Returns:
            List of user-defined LigandGroup objects with formal charges
        """
        self._print_syntax_help()

        while True:
            user_input = prompt_with_context(
                self.processor, "Enter ligand groups",
                module="Ligand Grouping",
                description="Define ligand atom groups",
            ).strip()

            if not user_input:
                self._print("[yellow]⚠ No groups defined. Please enter at least one group.[/yellow]")
                continue

            try:
                # Parse grouping
                groups = self.parse_grouping_input(user_input, idx_to_atom, redox_site)

                # Show summary
                self._print(f"\n[green]✓ Defined {len(groups)} ligand group(s):[/green]")
                for group in groups:
                    atom_summary = self._get_atom_summary(group.atoms)
                    self._print(f"  Group {group.group_number}: {atom_summary}")

                # Prompt for parameterization method and details per group
                self._print("\n[cyan]Configure parameterization for each ligand group:[/cyan]")
                for group in groups:
                    atom_summary = self._get_atom_summary(group.atoms)
                    self._print(f"\n[yellow]Group {group.group_number}: {atom_summary}[/yellow]")

                    # Ask for parameterization method
                    self._print("  Options:")
                    self._print("    1. Run antechamber (GAFF + AM1-BCC)")
                    self._print("    2. Use existing mol2/frcmod files")
                    self._print("    3. Manual entry")

                    while True:
                        method_input = prompt_with_context(
                            self.processor, "Choose method (1/2/3)",
                            module="Ligand Grouping",
                            description="Parameterization method for ligand group",
                            options_map={
                                "1": "Run antechamber (GAFF + AM1-BCC)",
                                "2": "Use existing mol2/frcmod files",
                                "3": "Manual entry",
                            },
                        ).strip()
                        if method_input in ['1', '2', '3']:
                            break
                        self._print(f"[red]Invalid choice: {method_input}. Enter 1, 2, or 3.[/red]")

                    if method_input == '1':
                        # Antechamber
                        group.parameterization_method = 'antechamber'
                        while True:
                            charge_input = prompt_with_context(
                                self.processor, "Formal charge",
                                module="Ligand Grouping",
                                description="Formal charge for ligand group",
                            ).strip()
                            try:
                                group.formal_charge = int(charge_input)
                                break
                            except ValueError:
                                self._print(f"[red]Invalid charge: {charge_input}. Please enter an integer.[/red]")

                    elif method_input == '2':
                        # Existing files
                        group.parameterization_method = 'existing'
                        group.existing_mol2 = prompt_with_context(
                            self.processor, "Path to mol2 file",
                            module="Ligand Grouping",
                            description="Path to existing mol2 file",
                        ).strip()
                        group.existing_frcmod = prompt_with_context(
                            self.processor, "Path to frcmod file (optional, press Enter to skip)",
                            default="",
                            module="Ligand Grouping",
                            description="Path to existing frcmod file (optional)",
                        ).strip()
                        if not group.existing_frcmod:
                            group.existing_frcmod = None

                    elif method_input == '3':
                        # Manual entry
                        group.parameterization_method = 'manual'

                # Confirm
                self._print("\n[bold cyan]Final configuration:[/bold cyan]")
                for group in groups:
                    atom_summary = self._get_atom_summary(group.atoms)
                    if group.parameterization_method == 'antechamber':
                        self._print(f"  Group {group.group_number}: {atom_summary}, antechamber (charge={group.formal_charge:+d})")
                    elif group.parameterization_method == 'existing':
                        self._print(f"  Group {group.group_number}: {atom_summary}, existing files")
                        self._print(f"    mol2: {group.existing_mol2}")
                        if group.existing_frcmod:
                            self._print(f"    frcmod: {group.existing_frcmod}")
                    elif group.parameterization_method == 'manual':
                        self._print(f"  Group {group.group_number}: {atom_summary}, manual entry")

                confirm = prompt_with_context(
                    self.processor, "Confirm grouping? (y/n)",
                    module="Ligand Grouping",
                    description="Confirm ligand grouping configuration",
                ).strip().lower()
                if confirm == 'y':
                    return groups
                else:
                    self._print("[yellow]Let's try again...[/yellow]")

            except Exception as e:
                self._print(f"[red]Error: {e}[/red]")
                self._print("[yellow]Please try again with correct syntax.[/yellow]")

    def parse_grouping_input(self, user_input: str,
                            idx_to_atom: Dict[int, RedoxSiteAtom],
                            redox_site: RedoxSite) -> List[LigandGroup]:
        """
        Parse user grouping syntax into LigandGroup objects.

        Syntax examples:
            "1:HEM"              -> Group 1 with all HEM atoms
            "1:1-45"             -> Group 1 with atoms 1-45
            "1:5,7,9"            -> Group 1 with atoms 5, 7, 9
            "1:1-10,15,20-25"    -> Group 1 with atoms 1-10, 15, and 20-25
            "1:HEM 2:46 3:47-50" -> 3 groups

        Args:
            user_input: User input string
            idx_to_atom: Mapping of index -> atom
            redox_site: RedoxSite object

        Returns:
            List of LigandGroup objects (without charges yet)
        """
        groups = []

        # Split by spaces to get each group definition
        for group_def in user_input.strip().split():
            if ':' not in group_def:
                raise ValueError(
                    f"Invalid syntax: '{group_def}'. Expected format 'N:selector' (e.g., '1:HEM' or '2:5-10')"
                )

            # Parse group number and selector
            parts = group_def.split(':', 1)
            if len(parts) != 2:
                raise ValueError(f"Invalid format: '{group_def}'. Use 'N:selector'")

            try:
                group_num = int(parts[0])
            except ValueError:
                raise ValueError(f"Group number must be an integer, got: '{parts[0]}'")

            selector = parts[1]

            # Parse selector to get atoms
            atoms = self._parse_selector(selector, idx_to_atom)

            if not atoms:
                raise ValueError(f"No atoms selected for group {group_num} with selector '{selector}'")

            # Create LigandGroup (charge will be set later)
            groups.append(LigandGroup(
                group_id=f"ligand_{group_num}",
                group_number=group_num,
                atoms=atoms,
                formal_charge=0  # Will be set by user later
            ))

        # Check for duplicate group numbers
        group_numbers = [g.group_number for g in groups]
        if len(group_numbers) != len(set(group_numbers)):
            raise ValueError("Duplicate group numbers detected. Each group must have a unique number.")

        return groups

    def _parse_selector(self, selector: str,
                       idx_to_atom: Dict[int, RedoxSiteAtom]) -> List[RedoxSiteAtom]:
        """
        Parse selector portion of group definition.

        Examples:
            "HEM"      -> All atoms with resname HEM
            "1-45"     -> Atoms at indices 1 through 45
            "5,7,9"    -> Atoms at indices 5, 7, 9
            "1-10,15,20-25" -> Atoms 1-10, 15, and 20-25

        Args:
            selector: Selector string
            idx_to_atom: Mapping of index -> atom

        Returns:
            List of selected atoms
        """
        # Check if selector is a residue name (alphabetic, uppercase)
        if selector.replace('_', '').replace('-', '').isalpha():
            return self._get_atoms_by_resname(selector.upper(), idx_to_atom)

        # Otherwise, parse as indices
        indices = self._parse_indices(selector)

        atoms = []
        for idx in indices:
            if idx not in idx_to_atom:
                raise ValueError(f"Invalid atom index: {idx}. Valid range: 1-{len(idx_to_atom)}")
            atoms.append(idx_to_atom[idx])

        return atoms

    def _parse_indices(self, selector: str) -> List[int]:
        """
        Parse index selector like "1-45", "5,7,9", or "1-10,15,20-25".

        Args:
            selector: Index selector string

        Returns:
            List of indices
        """
        indices = []

        for part in selector.split(','):
            part = part.strip()

            if '-' in part:
                # Range: "1-45"
                range_parts = part.split('-')
                if len(range_parts) != 2:
                    raise ValueError(f"Invalid range format: '{part}'. Expected 'start-end'")

                try:
                    start = int(range_parts[0])
                    end = int(range_parts[1])
                except ValueError:
                    raise ValueError(f"Invalid range: '{part}'. Start and end must be integers")

                if start > end:
                    raise ValueError(f"Invalid range: '{part}'. Start must be ≤ end")

                indices.extend(range(start, end + 1))

            else:
                # Single index: "5"
                try:
                    indices.append(int(part))
                except ValueError:
                    raise ValueError(f"Invalid index: '{part}'. Must be an integer")

        return indices

    def _get_atoms_by_resname(self, resname: str,
                              idx_to_atom: Dict[int, RedoxSiteAtom]) -> List[RedoxSiteAtom]:
        """
        Get all atoms with specified residue name.

        Args:
            resname: Residue name (uppercase)
            idx_to_atom: Mapping of index -> atom

        Returns:
            List of atoms with matching residue name
        """
        atoms = [atom for atom in idx_to_atom.values() if atom.resname == resname]

        if not atoms:
            available_resnames = sorted(set(atom.resname for atom in idx_to_atom.values()))
            raise ValueError(
                f"No atoms found with residue name '{resname}'. "
                f"Available residue names: {', '.join(available_resnames)}"
            )

        return atoms

    def _get_atom_summary(self, atoms: List[RedoxSiteAtom]) -> str:
        """Get concise summary of atoms for display."""
        if not atoms:
            return "Formula: (empty)"

        elements = [a.element for a in atoms]
        elem_counts = {}
        for e in elements:
            elem_counts[e] = elem_counts.get(e, 0) + 1

        elem_str = "".join(f"{elem}{count}" for elem, count in sorted(elem_counts.items()))
        return f"Formula: {elem_str}"

    def _print_syntax_help(self):
        """Print syntax help for user."""
        self._print("\n[bold cyan]Define Ligand Groups for Antechamber Parameterization:[/bold cyan]")
        self._print("\n[yellow]Syntax:[/yellow]")
        self._print("  • 'N:RESNAME'     = Group N with all atoms of residue RESNAME")
        self._print("  • 'N:X-Y'         = Group N with atoms at indices X through Y")
        self._print("  • 'N:X,Y,Z'       = Group N with specific atom indices")
        self._print("  • 'N:X-Y,Z,A-B'   = Group N with ranges and individual indices")
        self._print("\n[yellow]Examples:[/yellow]")
        self._print("  Heme (all atoms):          [green]1:HEM[/green]")
        self._print("  Fe4S4 (4 separate):        [green]1:46 2:47 3:48 4:49[/green]")
        self._print("  Mixed grouping:            [green]1:1-45 2:46,47 3:SF4[/green]")

    def _print(self, message: str):
        """Print with Rich console if available, otherwise plain print."""
        if self.console:
            self.console.print(message)
        else:
            # Strip Rich markup for plain print
            import re
            plain_msg = re.sub(r'\[/?\w+.*?\]', '', message)
            print(plain_msg)

    def auto_group_single_atoms(self, idx_to_atom: Dict[int, RedoxSiteAtom]) -> List[LigandGroup]:
        """
        Auto-generate groups with one atom per group (conservative fallback).

        Args:
            idx_to_atom: Mapping of index -> atom

        Returns:
            List of LigandGroup objects (one atom each)
        """
        groups = []
        for idx, atom in sorted(idx_to_atom.items()):
            groups.append(LigandGroup(
                group_id=f"ligand_{idx}",
                group_number=idx,
                atoms=[atom],
                formal_charge=0  # User will need to set manually
            ))
        return groups
