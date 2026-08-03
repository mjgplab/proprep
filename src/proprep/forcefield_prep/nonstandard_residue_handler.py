"""
Non-Standard Residue Handler

Interactive prompts for handling residues not found in force field libraries.
Allows users to provide mol2 files or manually input atom type data.

TODO: Add session recording context to prompts
----------------------------------------------
The Prompt.ask calls at lines ~337 and ~397 (in handle_unknown_residue and
handle_no_matches methods) need session context for proper session replay.

To fix:
1. Add 'processor' parameter to NonStandardResidueHandler.__init__
2. Update callers to pass processor (likely in forcefield_parameterizer.py)
3. Import prompt_with_context from proprep.application.menu_commands
4. Replace prompt_with_context(None, "Select option", ...) with prompt_with_context() calls
5. Build options_map dynamically since choices are generated at runtime
   based on matches found, antechamber availability, etc.

The prompts are for:
- Selecting from matched residue mappings
- Choosing antechamber generation vs mol2 input vs manual input vs skip
"""

from pathlib import Path
from typing import List, Set, Dict, Tuple, Optional
from dataclasses import dataclass
from rich.console import Console
from proprep.utils.prompts import prompt_with_context, confirm_with_context, int_prompt_with_context, float_prompt_with_context
from rich.panel import Panel
from rich.table import Table

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .forcefield_data import ForceFieldData


@dataclass
class ResidueMatch:
    """Match between PDB residue and force field residue."""
    ff_residue_name: str
    match_type: str  # "perfect" or "partial"
    pdb_atom_count: int
    ff_atom_count: int
    matched_atoms: Set[str]
    ff_extra_atoms: Set[str]


class NonStandardResidueHandler:
    """
    Handles interactive input for non-standard residues.

    Prompts users to provide either mol2 files or manual atom type data
    for residues not found in the force field library.
    """

    def __init__(self, ff_data: 'ForceFieldData', console: Console = None,
                 redox_site=None, pdb_file: str = None, processor=None):
        """
        Initialize handler.

        Args:
            ff_data: ForceFieldData instance to add definitions to
            console: Rich console for output
            redox_site: RedoxSite object (optional) for checking metal renaming
            pdb_file: Path to PDB file (optional) for extracting residues
            processor: Optional ProPrep processor for session recording context
        """
        self.ff_data = ff_data
        self.console = console or Console()
        self.redox_site = redox_site
        self.pdb_file = pdb_file
        self.processor = processor

    def handle_missing_residues(self, missing_residues: Set[str],
                                pdb_residue_atoms: Optional[Dict[str, Set[str]]] = None,
                                interactive: bool = True) -> bool:
        """
        Handle missing residues interactively.

        Args:
            missing_residues: Set of residue names not in force field
            pdb_residue_atoms: Dict mapping residue_name -> set of atom names in PDB
            interactive: If True, prompt user; if False, skip missing residues

        Returns:
            True if all missing residues were handled, False otherwise
        """
        if not missing_residues:
            return True

        if not interactive:
            self.console.print(f"[yellow]Warning: {len(missing_residues)} residues not in force field[/yellow]")
            self.console.print(f"[yellow]Missing: {', '.join(sorted(missing_residues))}[/yellow]")
            self.console.print("[yellow]Run in interactive mode to provide atom type data[/yellow]")
            return False

        # Show missing residues
        self.console.print("\n" + "="*70)
        self.console.print(f"[bold yellow]⚠️  {len(missing_residues)} Non-Standard Residue(s) Detected[/bold yellow]")
        self.console.print("="*70)

        # Show which force fields were selected
        if len(self.ff_data.loaded_leaprcs) == 1:
            self.console.print(f"The following residues are not in {self.ff_data.loaded_leaprcs[0]}:")
        else:
            self.console.print(f"The following residues are not in the selected force fields:")
            self.console.print(f"[grey50]  ({', '.join(self.ff_data.loaded_leaprcs)})[/grey50]")

        for res_name in sorted(missing_residues):
            self.console.print(f"  • {res_name}")
        self.console.print()

        # Ask if user wants to provide data
        if not confirm_with_context(
            self.processor,
            "Would you like to provide atom type data for these residues?",
            default=True,
            module="Non-Standard Residue Handler",
            description="Provide atom type data for non-standard residues",
        ):
            self.console.print("[yellow]Skipping non-standard residues[/yellow]")
            return False

        # Handle each residue
        for res_name in sorted(missing_residues):
            pdb_atoms = pdb_residue_atoms.get(res_name) if pdb_residue_atoms else None
            if not self._handle_single_residue(res_name, pdb_atoms):
                self.console.print(f"[red]Failed to get data for {res_name}[/red]")
                return False

        return True

    def _find_matching_residues(self, residue_name: str, pdb_atoms: Set[str]) -> List[ResidueMatch]:
        """
        Find force field residues that match the PDB residue based on atom names.

        Matching criteria:
        - Perfect match: All PDB atoms in FF, all FF atoms in PDB (bidirectional)
        - Partial match: All PDB atoms in FF, but FF has additional atoms

        Args:
            residue_name: PDB residue name
            pdb_atoms: Set of atom names in PDB residue

        Returns:
            List of ResidueMatch objects, sorted by quality (perfect first, then by atom count)
        """
        matches = []

        # Get all residue names in force field
        all_ff_residues = set()
        for (ff_res_name, _) in self.ff_data.atom_definitions.keys():
            all_ff_residues.add(ff_res_name)

        # Check each force field residue
        for ff_res_name in all_ff_residues:
            # Get all atoms for this FF residue
            ff_atoms_list = self.ff_data.get_residue_atoms(ff_res_name)
            ff_atoms = {atom.atom_name for atom in ff_atoms_list}

            # Check if all PDB atoms are in FF residue
            if pdb_atoms.issubset(ff_atoms):
                matched_atoms = pdb_atoms & ff_atoms
                ff_extra_atoms = ff_atoms - pdb_atoms

                # Determine match type
                if len(ff_extra_atoms) == 0:
                    match_type = "perfect"
                else:
                    match_type = "partial"

                matches.append(ResidueMatch(
                    ff_residue_name=ff_res_name,
                    match_type=match_type,
                    pdb_atom_count=len(pdb_atoms),
                    ff_atom_count=len(ff_atoms),
                    matched_atoms=matched_atoms,
                    ff_extra_atoms=ff_extra_atoms
                ))

        # Sort: perfect matches first, then by smallest number of extra atoms
        matches.sort(key=lambda m: (0 if m.match_type == "perfect" else 1, len(m.ff_extra_atoms)))

        return matches

    def _check_if_will_rename(self, residue_name: str) -> Optional[Dict[str, int]]:
        """
        Check if residue atoms will be systematically renamed (M1-M9, L1-L9).

        Args:
            residue_name: Residue name to check

        Returns:
            Dict with total_count, metal_count, ligand_count, other_count, or None if no renaming
        """
        if not self.redox_site:
            return None

        # Get atoms from this residue in RedoxSite
        residue_atoms_coords = [
            atom.coords for atom in self.redox_site.atoms
            if atom.resname == residue_name
        ]

        if not residue_atoms_coords:
            return None

        total_count = len(residue_atoms_coords)

        # Count metals (atoms that are centers)
        metal_coords = {center.coords for center in self.redox_site.centers}
        metal_count = sum(1 for coord in residue_atoms_coords if coord in metal_coords)

        # Count ligands (atoms bonded to metals via coordinate bonds)
        ligand_coords = set()
        for bond in self.redox_site.bonds:
            if bond.chemical_type == 'coordinate':
                # Check which end is the metal
                if bond.atom1_coords in metal_coords:
                    ligand_coords.add(bond.atom2_coords)
                elif bond.atom2_coords in metal_coords:
                    ligand_coords.add(bond.atom1_coords)

        ligand_count = sum(1 for coord in residue_atoms_coords if coord in ligand_coords)

        # Count atoms that won't be renamed
        other_count = total_count - metal_count - ligand_count

        if metal_count > 0 or ligand_count > 0:
            return {
                'total_count': total_count,
                'metal_count': metal_count,
                'ligand_count': ligand_count,
                'other_count': other_count
            }

        return None

    def _residue_contains_metal(self, residue_name: str) -> bool:
        """
        Check if residue contains metal centers from RedoxSite.

        Args:
            residue_name: Residue name to check

        Returns:
            True if residue contains at least one metal center
        """
        if not self.redox_site:
            return False

        metal_coords = {center.coords for center in self.redox_site.centers}

        for atom in self.redox_site.atoms:
            if atom.resname == residue_name and atom.coords in metal_coords:
                return True

        return False

    def _handle_single_residue(self, residue_name: str, pdb_atoms: Optional[Set[str]] = None) -> bool:
        """
        Handle a single non-standard residue.

        Args:
            residue_name: Residue name
            pdb_atoms: Set of atom names in PDB residue (for matching)

        Returns:
            True if successful
        """
        self.console.print(f"\n[bold cyan]Handling residue: {residue_name}[/bold cyan]")

        # Check if atoms will be renamed (metal centers or ligands)
        will_rename = self._check_if_will_rename(residue_name)
        if will_rename:
            self.console.print(f"\n[grey50]ℹ️  {residue_name} is part of metal site ({will_rename['total_count']} atoms total):[/grey50]")

            if will_rename['metal_count'] > 0:
                self.console.print(f"[grey50]  • {will_rename['metal_count']} atom(s) → M1-M9 (metal centers)[/grey50]")

            if will_rename['ligand_count'] > 0:
                self.console.print(f"[grey50]  • {will_rename['ligand_count']} atom(s) → L1-L9 (ligands)[/grey50]")

            if will_rename['other_count'] > 0:
                self.console.print(f"[grey50]  • {will_rename['other_count']} atom(s) → No renaming (not directly coordinating)[/grey50]")

            self.console.print(f"[grey50]  • Renamed atom types/charges will be replaced in Step 3 (RESP)[/grey50]\n")

        # Try to find matching residues based on atom names
        matches = []
        if pdb_atoms:
            matches = self._find_matching_residues(residue_name, pdb_atoms)

        # Check if contains metals (for antechamber option)
        contains_metal = self._residue_contains_metal(residue_name)

        # Show matching suggestions if found
        if matches:
            self.console.print(f"\n⚠️  {residue_name} not found in force field")
            self.console.print("🔍 Found matching force field residues based on atom names:\n")

            # Create match table
            match_table = Table(show_header=True, header_style="bold magenta")
            match_table.add_column("Option", style="yellow", no_wrap=True)
            match_table.add_column("FF Residue", style="cyan", no_wrap=True)
            match_table.add_column("Match Quality", style="green")
            match_table.add_column("Atoms", style="grey50")

            option_num = 1
            for match in matches[:5]:  # Show top 5 matches
                if match.match_type == "perfect":
                    quality = f"Perfect match ({match.pdb_atom_count}/{match.ff_atom_count} atoms)"
                else:
                    extra_count = len(match.ff_extra_atoms)
                    quality = f"Partial match ({match.pdb_atom_count}/{match.ff_atom_count} atoms) - FF has {extra_count} additional"

                atoms_str = ", ".join(sorted(match.matched_atoms)[:8])
                if len(match.matched_atoms) > 8:
                    atoms_str += "..."

                match_table.add_row(str(option_num), match.ff_residue_name, quality, atoms_str)
                option_num += 1

            self.console.print(match_table)
            self.console.print()

            # Build choice list
            choices = [str(i+1) for i in range(len(matches[:5]))]
            next_option = len(matches[:5]) + 1

            self.console.print(f"[yellow]Choose option:[/yellow]")
            for i, match in enumerate(matches[:5], 1):
                self.console.print(f"  {i}. Map to {match.ff_residue_name} ({match.match_type} match)")

            # Antechamber option (disabled for metals)
            if not contains_metal:
                self.console.print(f"  {next_option}. Auto-generate with reduce + antechamber")
                antechamber_option = next_option
                choices.append(str(next_option))
                next_option += 1
            else:
                self.console.print(f"  [grey50]{next_option}. Auto-generate with reduce + antechamber (disabled - contains metal)[/grey50]")
                antechamber_option = None
                next_option += 1

            self.console.print(f"  {next_option}. Provide mol2 file")
            mol2_option = next_option
            choices.append(str(next_option))
            next_option += 1

            self.console.print(f"  {next_option}. Manually input atom data")
            manual_option = next_option
            choices.append(str(next_option))
            next_option += 1

            self.console.print(f"  {next_option}. Skip this residue")
            skip_option = next_option
            choices.append(str(next_option))

            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=choices,
                default="1",
                module="Non-Standard Residue Handler",
                description=f"Select handling option for residue {residue_name}",
            )

            choice_num = int(choice)

            # Handle mapping choices
            if choice_num <= len(matches[:5]):
                selected_match = matches[choice_num - 1]
                return self._handle_residue_mapping(residue_name, selected_match)

            # Handle other options
            if antechamber_option and choice_num == antechamber_option:
                return self._handle_antechamber_generation(residue_name, pdb_atoms)
            elif choice_num == mol2_option:
                return self._handle_mol2_input(residue_name)
            elif choice_num == manual_option:
                return self._handle_manual_input(residue_name, pdb_atoms)
            elif choice_num == skip_option:
                self.console.print(f"[yellow]Skipping {residue_name}[/yellow]")
                return False
            else:
                self.console.print(f"[yellow]Skipping {residue_name}[/yellow]")
                return False

        else:
            # No matches found - show original options
            self.console.print("\n[yellow]Choose how to provide atom type data:[/yellow]")

            option_num = 1
            choices = []

            # Antechamber option (disabled for metals)
            if not contains_metal:
                self.console.print(f"  {option_num}. Auto-generate with reduce + antechamber")
                antechamber_option = option_num
                choices.append(str(option_num))
                option_num += 1
            else:
                self.console.print(f"  [grey50]{option_num}. Auto-generate with reduce + antechamber (disabled - contains metal)[/grey50]")
                antechamber_option = None
                option_num += 1

            self.console.print(f"  {option_num}. Provide mol2 file")
            mol2_option = option_num
            choices.append(str(option_num))
            option_num += 1

            self.console.print(f"  {option_num}. Manually input atom data")
            manual_option = option_num
            choices.append(str(option_num))
            option_num += 1

            self.console.print(f"  {option_num}. Skip this residue")
            skip_option = option_num
            choices.append(str(option_num))

            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=choices,
                default="1" if not contains_metal else str(mol2_option),
                module="Non-Standard Residue Handler",
                description=f"Select handling method for unknown residue {residue_name}",
            )

            choice_num = int(choice)

            if antechamber_option and choice_num == antechamber_option:
                return self._handle_antechamber_generation(residue_name, pdb_atoms)
            elif choice_num == mol2_option:
                return self._handle_mol2_input(residue_name)
            elif choice_num == manual_option:
                return self._handle_manual_input(residue_name, pdb_atoms)
            elif choice_num == skip_option:
                self.console.print(f"[yellow]Skipping {residue_name}[/yellow]")
                return False
            else:
                self.console.print(f"[yellow]Skipping {residue_name}[/yellow]")
                return False

    def _handle_residue_mapping(self, source_residue: str, match: ResidueMatch) -> bool:
        """
        Map a PDB residue to a force field residue.

        Args:
            source_residue: PDB residue name
            match: ResidueMatch object with target FF residue

        Returns:
            True if successful
        """
        target_residue = match.ff_residue_name

        self.console.print(f"\n[cyan]Mapping {source_residue} → {target_residue}[/cyan]")

        # Get all atoms from target residue
        target_atoms = self.ff_data.get_residue_atoms(target_residue)

        if not target_atoms:
            self.console.print(f"[red]Error: No atoms found for {target_residue}[/red]")
            return False

        # Copy all atom definitions to source residue name
        atom_count = 0
        for atom_def in target_atoms:
            self.ff_data.add_manual_atom(
                source_residue,
                atom_def.atom_name,
                atom_def.atom_type,
                atom_def.charge
            )
            atom_count += 1

        self.console.print(f"[green]✓ Mapped {atom_count} atoms from {target_residue} to {source_residue}[/green]")

        # Show what was mapped
        if match.ff_extra_atoms:
            self.console.print(f"[grey50]Note: {target_residue} has {len(match.ff_extra_atoms)} additional atoms not in PDB:[/grey50]")
            self.console.print(f"[grey50]  {', '.join(sorted(match.ff_extra_atoms))}[/grey50]")

        return True

    def _handle_mol2_input(self, residue_name: str) -> bool:
        """
        Handle mol2 file input.

        Args:
            residue_name: Residue name

        Returns:
            True if successful
        """
        self.console.print(f"\n[cyan]Please provide a mol2 file for {residue_name}[/cyan]")
        self.console.print("[grey50]The mol2 file should contain atom names, types, and charges[/grey50]")

        # List mol2 files in current working directory
        import os
        cwd = Path.cwd()
        mol2_files = sorted(cwd.glob("*.mol2"))

        if mol2_files:
            self.console.print(f"\n[grey50]Found {len(mol2_files)} mol2 file(s) in {cwd.name}/:[/grey50]")
            for mol2_file in mol2_files:
                self.console.print(f"[grey50]  • {mol2_file.name}[/grey50]")
            self.console.print()

        while True:
            mol2_path_str = prompt_with_context(
                self.processor,
                "Path to mol2 file (or filename if in current directory)",
                module="Non-Standard Residue Handler",
                description="mol2 file path for non-standard residue",
            )

            if not mol2_path_str:
                if confirm_with_context(
                    self.processor,
                    "No file provided. Try again?",
                    default=True,
                    module="Non-Standard Residue Handler",
                    description="Retry mol2 file path entry",
                ):
                    continue
                else:
                    return False

            mol2_path = Path(mol2_path_str).expanduser()

            # If path doesn't exist, try as filename in current directory
            if not mol2_path.exists() and not mol2_path.is_absolute():
                mol2_path = cwd / mol2_path_str

            if not mol2_path.exists():
                self.console.print(f"[red]File not found: {mol2_path}[/red]")
                if confirm_with_context(
                    self.processor,
                    "Try again?",
                    default=True,
                    module="Non-Standard Residue Handler",
                    description="Retry mol2 file entry",
                ):
                    continue
                else:
                    return False

            # Try to parse the mol2 file
            atom_count = self.ff_data.add_mol2_file(mol2_path, residue_name)

            if atom_count > 0:
                return True
            else:
                self.console.print("[red]Failed to parse mol2 file[/red]")
                if confirm_with_context(
                    self.processor,
                    "Try another file?",
                    default=True,
                    module="Non-Standard Residue Handler",
                    description="Try another mol2 file",
                ):
                    continue
                else:
                    return False

    def _handle_manual_input(self, residue_name: str, pdb_atoms: Optional[Set[str]] = None) -> bool:
        """
        Handle manual atom data input with interactive table.

        Args:
            residue_name: Residue name
            pdb_atoms: Set of atom names found in PDB (optional)

        Returns:
            True if successful
        """
        from rich.live import Live
        from proprep.utils.prompts import float_prompt_with_context

        if not pdb_atoms:
            self.console.print(f"[yellow]No atoms found in PDB for {residue_name}[/yellow]")
            return False

        self.console.print(f"\n[cyan]Manual atom type input for {residue_name}[/cyan]")
        self.console.print(f"[grey50]Enter atom type and charge for each atom.[/grey50]\n")

        # Create data structure for atoms
        atom_list = sorted(pdb_atoms)
        atom_data = {}
        for atom_name in atom_list:
            atom_data[atom_name] = {
                'type': None,
                'charge': None,
                'status': 'pending'
            }

        def create_table():
            """Create the current state table."""
            table = Table(
                title=f"Atom Parameters for {residue_name}",
                show_header=True,
                header_style="bold cyan"
            )
            table.add_column("Atom Name", style="yellow", no_wrap=True)
            table.add_column("Atom Type", style="green")
            table.add_column("Charge", style="magenta", justify="right")
            table.add_column("Status", style="grey50")

            for atom_name in atom_list:
                data = atom_data[atom_name]
                atom_type = data['type'] if data['type'] else '?'
                charge = f"{data['charge']:.4f}" if data['charge'] is not None else '?'

                if data['status'] == 'complete':
                    status = '✓'
                elif data['status'] == 'current':
                    status = '→'
                else:
                    status = '⏸'

                table.add_row(atom_name, atom_type, charge, status)

            return table

        # Show initial table
        self.console.print(create_table())
        self.console.print()

        # Loop through each atom
        for atom_name in atom_list:
            atom_data[atom_name]['status'] = 'current'

            # Prompt for type
            atom_type = prompt_with_context(
                self.processor,
                f"[cyan]{atom_name}[/cyan] - Atom type",
                module="Non-Standard Residue Handler",
                description=f"Atom type for {atom_name}",
            )
            if not atom_type:
                self.console.print("[yellow]Skipping atom[/yellow]")
                atom_data[atom_name]['status'] = 'pending'
                continue

            # Prompt for charge
            try:
                charge = float_prompt_with_context(
                    self.processor,
                    f"[cyan]{atom_name}[/cyan] - Partial charge",
                    module="Non-Standard Residue Handler",
                    description=f"Partial charge for {atom_name}",
                )
            except Exception:
                self.console.print("[yellow]Skipping atom[/yellow]")
                atom_data[atom_name]['status'] = 'pending'
                continue

            # Update data
            atom_data[atom_name]['type'] = atom_type
            atom_data[atom_name]['charge'] = charge
            atom_data[atom_name]['status'] = 'complete'

            # Show updated table
            self.console.print()
            self.console.print(create_table())
            self.console.print()

        # Show final summary
        self.console.print("\n[bold]Final Summary:[/bold]")
        final_table = create_table()
        self.console.print(final_table)

        # Count completed atoms
        completed = sum(1 for data in atom_data.values() if data['status'] == 'complete')

        if completed == 0:
            self.console.print(f"\n[yellow]No atoms completed for {residue_name}[/yellow]")
            return False

        self.console.print(f"\n[cyan]Completed {completed}/{len(atom_list)} atoms[/cyan]")

        # Confirm
        if not confirm_with_context(
            self.processor,
            "Accept these assignments?",
            default=True,
            module="Non-Standard Residue Handler",
            description="Accept manual atom type/charge assignments",
        ):
            if confirm_with_context(
                self.processor,
                "Try again?",
                default=True,
                module="Non-Standard Residue Handler",
                description="Retry manual atom assignments",
            ):
                return self._handle_manual_input(residue_name, pdb_atoms)
            return False

        # Add to force field collector
        atoms_added = 0
        for atom_name, data in atom_data.items():
            if data['status'] == 'complete':
                self.ff_data.add_manual_atom(
                    residue_name,
                    atom_name,
                    data['type'],
                    data['charge']
                )
                atoms_added += 1

        if atoms_added > 0:
            self.console.print(f"\n[green]✓ Added {atoms_added} atoms for {residue_name}[/green]")
            return True
        else:
            self.console.print(f"[yellow]No atoms added for {residue_name}[/yellow]")
            return False

    def _extract_residue_from_pdb(self, residue_name: str, output_pdb: Path) -> bool:
        """
        Extract a specific residue from PDB file using grep/awk.

        Args:
            residue_name: Residue name to extract
            output_pdb: Output PDB file path

        Returns:
            True if successful
        """
        if not self.pdb_file:
            self.console.print("[red]Error: No PDB file path provided[/red]")
            return False

        import subprocess

        try:
            # Use grep to extract ATOM/HETATM lines for this residue
            # PDB format: columns 18-20 contain residue name
            cmd = f"grep -E '^(ATOM|HETATM)' {self.pdb_file} | awk '{{if (substr($0,18,3) ~ /{residue_name}/) print}}' > {output_pdb}"

            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                timeout=30
            )

            if result.returncode != 0:
                self.console.print(f"[red]Error extracting residue: {result.stderr}[/red]")
                return False

            # Check if file has content
            if not output_pdb.exists() or output_pdb.stat().st_size == 0:
                self.console.print(f"[red]Error: No atoms found for residue {residue_name}[/red]")
                return False

            return True

        except subprocess.TimeoutExpired:
            self.console.print("[red]Error: PDB extraction timed out[/red]")
            return False
        except Exception as e:
            self.console.print(f"[red]Error extracting residue: {str(e)}[/red]")
            return False

    def _handle_antechamber_generation(self, residue_name: str, pdb_atoms: Optional[Set[str]] = None) -> bool:
        """
        Generate mol2 using reduce + antechamber workflow.

        Args:
            residue_name: Residue name
            pdb_atoms: Set of atom names in PDB (optional)

        Returns:
            True if successful
        """
        import subprocess
        from proprep.utils.prompts import int_prompt_with_context

        self.console.print(f"\n[bold cyan]Auto-generating parameters for {residue_name}[/bold cyan]")
        self.console.print("[grey50]This will use antechamber to assign atom types and BCC charges[/grey50]\n")

        # Ask if user wants to use reduce
        use_reduce = confirm_with_context(
            self.processor,
            "Add hydrogens with reduce first?",
            default=True,
            module="Non-Standard Residue Handler",
            description="Use AmberTools 'reduce' to add hydrogens first",
        )

        # Get charge from user
        charge = int_prompt_with_context(
            self.processor,
            f"Net charge for {residue_name}",
            default=0,
            module="Non-Standard Residue Handler",
            description=f"Net charge for {residue_name}",
        )

        # Create working directory
        work_dir = Path(f"antechamber_{residue_name}")
        work_dir.mkdir(exist_ok=True)

        self.console.print(f"\n[grey50]Working directory: {work_dir.absolute()}[/grey50]\n")

        try:
            # Step 1: Extract residue from PDB
            step_num = 1
            total_steps = 3 if use_reduce else 2

            self.console.print(f"[cyan]{step_num}/{total_steps}[/cyan] Extracting {residue_name} from PDB...")
            residue_pdb = work_dir / f"{residue_name}.pdb"

            if not self._extract_residue_from_pdb(residue_name, residue_pdb):
                raise RuntimeError(f"Failed to extract {residue_name} from PDB")

            self.console.print(f"[green]✓ Extracted to {residue_pdb.name}[/green]")

            # Determine input file for antechamber
            antechamber_input = residue_pdb

            # Step 2 (optional): Add hydrogens with reduce
            if use_reduce:
                step_num += 1
                self.console.print(f"\n[cyan]{step_num}/{total_steps}[/cyan] Adding hydrogens with reduce...")
                residue_h_pdb = work_dir / f"{residue_name}_H.pdb"

                reduce_cmd = f"reduce {residue_pdb} > {residue_h_pdb}"
                self.console.print(f"[grey50]Command: {reduce_cmd}[/grey50]")

                result = subprocess.run(
                    reduce_cmd,
                    shell=True,
                    capture_output=True,
                    text=True,
                    timeout=120
                )

                if result.returncode != 0 or not residue_h_pdb.exists():
                    raise RuntimeError(f"reduce failed: {result.stderr}")

                self.console.print(f"[green]✓ Added hydrogens to {residue_h_pdb.name}[/green]")

                # Give user option to inspect/edit the file
                self.console.print(f"\n[yellow]Protonation complete. File location:[/yellow]")
                self.console.print(f"[yellow]  {residue_h_pdb.absolute()}[/yellow]")

                if confirm_with_context(
                    self.processor,
                    "\nWould you like to inspect/edit the file before running antechamber?",
                    default=False,
                    module="Non-Standard Residue Handler",
                    description="Inspect/edit PDB file before antechamber",
                ):
                    self.console.print(f"\n[cyan]Please edit the file:[/cyan] {residue_h_pdb.absolute()}")
                    self.console.print(f"[grey50]You may need to fix protonation assignments.[/grey50]")

                    if not confirm_with_context(
                        self.processor,
                        "\nPress Enter when done editing and ready to continue",
                        default=True,
                        module="Non-Standard Residue Handler",
                        description="Continue after editing PDB",
                    ):
                        self.console.print("[yellow]Aborting antechamber generation[/yellow]")
                        return False

                antechamber_input = residue_h_pdb

            # Step 3: Generate mol2 with antechamber
            step_num += 1
            self.console.print(f"\n[cyan]{step_num}/{total_steps}[/cyan] Generating mol2 with antechamber (BCC charges)...")
            mol2_file = work_dir / f"{residue_name}.mol2"

            antechamber_cmd = (
                f"antechamber -fi pdb -fo mol2 "
                f"-i {antechamber_input} -o {mol2_file} "
                f"-c bcc -nc {charge} -pf y"
            )
            self.console.print(f"[grey50]Command: {antechamber_cmd}[/grey50]")

            result = subprocess.run(
                antechamber_cmd,
                shell=True,
                capture_output=True,
                text=True,
                timeout=300
            )

            if result.returncode != 0 or not mol2_file.exists():
                self.console.print(f"[red]Error output:[/red]")
                self.console.print(f"[grey50]{result.stderr}[/grey50]")
                raise RuntimeError(f"antechamber failed")

            self.console.print(f"[green]✓ Generated {mol2_file.name}[/green]")

            # Step 4: Load the mol2 into force field collector
            self.console.print(f"\n[cyan]Loading parameters...[/cyan]")
            atom_count = self.ff_data.add_mol2_file(mol2_file, residue_name)

            if atom_count > 0:
                self.console.print(f"[green]✓ Successfully loaded {atom_count} atoms[/green]")
                self.console.print(f"[grey50]Intermediate files saved in: {work_dir.absolute()}[/grey50]\n")
                return True
            else:
                raise RuntimeError("No atoms found in generated mol2 file")

        except subprocess.TimeoutExpired:
            self.console.print(f"[red]❌ Command timed out[/red]")
            self.console.print(f"[yellow]Intermediate files saved in: {work_dir.absolute()}[/yellow]\n")

            if confirm_with_context(
                self.processor,
                "Try another method?",
                default=True,
                module="Non-Standard Residue Handler",
                description="Try a different handling method for residue",
            ):
                return self._handle_single_residue(residue_name, pdb_atoms)
            return False

        except Exception as e:
            self.console.print(f"[red]❌ Antechamber generation failed: {str(e)}[/red]")
            self.console.print(f"[yellow]Intermediate files saved in: {work_dir.absolute()}[/yellow]\n")

            if confirm_with_context(
                self.processor,
                "Try another method?",
                default=True,
                module="Non-Standard Residue Handler",
                description="Try a different handling method for residue",
            ):
                return self._handle_single_residue(residue_name, pdb_atoms)
            return False

    def show_atom_summary(self, residue_name: str):
        """
        Show a summary of atoms added for a residue.

        Args:
            residue_name: Residue name
        """
        atoms = self.ff_data.get_residue_atoms(residue_name)

        if not atoms:
            self.console.print(f"[yellow]No atoms found for {residue_name}[/yellow]")
            return

        table = Table(title=f"Atoms for {residue_name}")
        table.add_column("Atom Name", style="cyan")
        table.add_column("Atom Type", style="green")
        table.add_column("Charge", style="yellow", justify="right")

        for atom in atoms:
            table.add_row(atom.atom_name, atom.atom_type, f"{atom.charge:.6f}")

        self.console.print(table)
