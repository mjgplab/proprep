"""
Amino Acid Mutator Module

Allows users to specify and apply amino acid mutations to PDB structures.
Displays sequence with numbering for easy mutation specification.
Integrates with MODELLER for structural modeling of mutations.
"""

import os
from typing import Dict, Any, Optional, List, Tuple
from pathlib import Path
import tempfile
import re
import logging

from Bio.PDB import Structure, Chain, Residue
from Bio.PDB.PDBIO import PDBIO

from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from proprep.utils.module_registry import ProcessingModule, register_module
from proprep.utils.prompts import prompt_with_context, confirm_with_context
from proprep.structure_prep.amino_acid_mutator_commands import (
    DisplaySequenceCommand,
    SpecifyMutationsCommand,
)

# Try to import MODELLER for structural modeling
# Note: Mutations will be applied by Structure Completeness module using MODELLER
# Suppress C-level error output from _modeller.mod_start().
try:
    import ctypes as _ctypes
    _devnull = os.open(os.devnull, os.O_WRONLY)
    _old_stdout, _old_stderr = os.dup(1), os.dup(2)
    try:
        os.dup2(_devnull, 1)
        os.dup2(_devnull, 2)
        import modeller
        from modeller.automodel import AutoModel
        HAS_MODELLER = True
    finally:
        _ctypes.CDLL(None).fflush(None)
        os.dup2(_old_stdout, 1)
        os.dup2(_old_stderr, 2)
        os.close(_devnull)
        os.close(_old_stdout)
        os.close(_old_stderr)
except Exception:
    HAS_MODELLER = False


@register_module
class AminoAcidMutator(ProcessingModule):
    """Amino Acid Mutator module for MPSA."""
    
    NAME = "Amino Acid Mutator"
    CATEGORY = "structure_prep"
    DESCRIPTION = "Mutate residues in proteins"
    VERSION = "1.0.0"
    REQUIRES = ["PDB Loader"]
    PRIORITY = 3  # After PDB Filter, before Structure Completeness
    
    # Standard amino acid residues
    STD_AMINO_ACIDS = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        "MSE": "M",  # Selenomethionine
    }
    
    # Reverse mapping for display
    AA_THREE_LETTER = {v: k for k, v in STD_AMINO_ACIDS.items() if k != "MSE"}
    AA_THREE_LETTER["M"] = "MET"  # Prefer MET over MSE for display
    
    def __init__(self):
        super().__init__()
        self.mutations = []  # List of (chain_id, res_num, from_aa, to_aa) for standard mutations
        self.nonstandard_mutations = []  # List of (chain_id, res_num, from_aa, to_aa, atoms_to_keep)
        self.sequence_data = {}  # Store sequence information per chain
        
    def set_processor(self, processor):
        """Set the processor reference."""
        self.processor = processor
        
    @property 
    def console(self):
        """Get console from processor if available."""
        if hasattr(self, 'processor') and self.processor and hasattr(self.processor, 'console'):
            return self.processor.console
        else:
            return Console()
    
    def _get_structure_from_workspace(self, workspace: Dict[str, Any], silent: bool = False):
        """Get structure object from workspace using structure selector.

        Args:
            workspace: Current workspace
            silent: If True, suppress console output

        Returns:
            BioPython Structure object, or None if no structure available
        """
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(workspace, self.console, processor=self.processor)

        # Use priority-based structure object selection
        structure = selector.get_structure_object(silent=silent)

        return structure

    def _get_pdb_file_from_workspace(self, workspace: Dict[str, Any]) -> str:
        """Get PDB file from workspace using structure selector.

        Args:
            workspace: Current workspace

        Returns:
            Path to PDB file, or None if no structure available
        """
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(workspace, self.console, self.processor)

        # Check if interactive mode (has processor with workspace)
        if self.processor and hasattr(self.processor, 'workspace'):
            # Interactive mode - let user select
            pdb_file = selector.get_structure(interactive=True)
        else:
            # Automatic mode - use priority selection
            pdb_file = selector.get_structure(interactive=False)

        return pdb_file
        
    def get_menu_options(self) -> Dict[str, str]:
        """Get available menu options."""
        return {
            "display": "Display protein sequences with numbering",
            "specify": "Specify amino acid mutations",
            "status": "Show current mutation status"
        }

    def get_enhanced_menu_options(self, workspace):
        """
        Get menu options with enhanced status information.

        Args:
            workspace: Current workspace

        Returns:
            List of MenuOption objects with status
        """
        from proprep.utils.enhanced_menu import MenuOption, OptionStatus

        options = []

        # Check workspace state
        has_sequence_data = bool(self.sequence_data)  # Sequences have been displayed
        stored_mutations = workspace.get("pending_mutations", [])
        stored_nonstandard = workspace.get("pending_nonstandard_mutations", [])
        has_mutations = len(stored_mutations) > 0 or len(stored_nonstandard) > 0

        # Option 1: Display sequences - always available or completed
        options.append(MenuOption(
            key="1",
            description="Display protein sequences with numbering",
            status=OptionStatus.COMPLETED if has_sequence_data else OptionStatus.AVAILABLE
        ))

        # Option 2: Specify mutations - requires sequences to be displayed first
        if has_sequence_data and has_mutations:
            status = OptionStatus.COMPLETED
            dep_text = ""
        elif has_sequence_data:
            status = OptionStatus.READY
            dep_text = ""
        else:
            status = OptionStatus.BLOCKED
            dep_text = "[Need to display sequences first] ○"

        options.append(MenuOption(
            key="2",
            description="Specify amino acid mutations",
            status=status,
            dependency_text=dep_text
        ))

        # Option 3: Show status - requires mutations to be specified
        if has_mutations:
            status = OptionStatus.READY
            dep_text = ""
        else:
            status = OptionStatus.BLOCKED
            dep_text = "[Need to specify mutations first] ○"

        options.append(MenuOption(
            key="3",
            description="Show current mutation status",
            status=status,
            dependency_text=dep_text
        ))

        return options

    def get_menu_suggestion(self, workspace):
        """
        Get a suggestion for the next recommended action.

        Args:
            workspace: Current workspace

        Returns:
            Suggestion text or None
        """
        has_sequence_data = bool(self.sequence_data)
        stored_mutations = workspace.get("pending_mutations", [])
        stored_nonstandard = workspace.get("pending_nonstandard_mutations", [])
        has_mutations = len(stored_mutations) > 0 or len(stored_nonstandard) > 0

        if not has_sequence_data:
            return "Start by displaying sequences (option 1) to see residue numbering for mutation specification"
        elif not has_mutations:
            return "Specify mutations (option 2) using format: CHAIN:POSITION:FROM->TO (e.g., A:123:ALA->GLY)"
        else:
            num_mutations = len(stored_mutations) + len(stored_nonstandard)
            return f"✓ {num_mutations} mutation(s) ready. View status with option 3, or press [m] to return to the main menu, then open the Structure Fixer module"
    
    def get_workspace_requirements(self) -> List[str]:
        """List what this module needs from workspace - needs at least one structure loaded."""
        return [
            "rcsb_pdb_file | local_pdb_file | alphafold_pdb_file | alphafill_pdb_file | alphafold_homolog_pdb_file | aligned_target_pdb_file | aligned_ref_pdb_file"
        ]

    def get_workspace_outputs(self) -> List[str]:
        """Get workspace outputs"""
        return ["pending_mutations", "pending_nonstandard_mutations"]

    def can_process(self, workspace) -> bool:
        """Check if module can process current workspace."""
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(workspace, self.console)
        status = selector.get_structure_status()
        return status.get("has_any", False)
        
    def handle_menu_option(self, option: str) -> bool:
        """Handle menu option selection using command pattern."""
        if option == "display":
            command = DisplaySequenceCommand(self.processor)
            return command.execute()
        elif option == "specify":
            command = SpecifyMutationsCommand(self.processor)
            return command.execute()
        elif option == "status":
            return self._show_mutation_status()
        return False
    
    def _display_sequences(self) -> bool:
        """Display protein sequences with numbering for easy mutation specification."""
        workspace = self.processor.workspace

        # Check for filtered structure first - use it automatically if available
        # Use StructureSelector with priority for filtered structure
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(workspace, self.console, processor=self.processor)

        # First try filtered_structure specifically
        structure = selector.get_structure_by_key("filtered_structure", require_exists=False)
        structure_type = "filtered"

        if not structure:
            # No filtered structure - select interactively using priority order
            result = selector.get_structure_object(interactive=True, return_key=True)
            if result:
                structure, structure_key = result
                structure_type = "selected"
            else:
                structure = None

        if not structure:
            self.console.print("[red]No structure loaded. Please load a PDB file first.[/red]")
            return False

        self.console.print("\n[bold cyan]Protein Sequences with Residue Numbers[/bold cyan]")
        self.console.print("Use these numbers to specify mutations")
        self.console.print("Example: [yellow]A:123:ALA->GLY[/yellow] (Chain A, position 123, Alanine to Glycine)")

        if structure_type == "filtered":
            self.console.print("[grey50]Showing filtered structure sequences[/grey50]")
        else:
            self.console.print("[grey50]Showing selected structure sequences[/grey50]")
        
        self.sequence_data = {}
        
        for model in structure:
            for chain in model:
                chain_id = chain.id
                
                # Extract protein residues only
                protein_residues = []
                for residue in chain:
                    if residue.id[0] == " " and residue.resname in self.STD_AMINO_ACIDS:
                        protein_residues.append(residue)
                
                if not protein_residues:
                    continue
                
                # Store sequence data
                self.sequence_data[chain_id] = []
                
                self.console.print(f"\n[bold]Chain {chain_id}:[/bold]")
                
                # Display sequence in blocks of 60 residues (matches Structure Completeness)
                BLOCK_SIZE = 60
                
                for start_idx in range(0, len(protein_residues), BLOCK_SIZE):
                    end_idx = min(start_idx + BLOCK_SIZE, len(protein_residues))
                    block_residues = protein_residues[start_idx:end_idx]
                    
                    # Store residue info
                    for residue in block_residues:
                        res_num = residue.id[1]
                        res_name = residue.resname
                        aa_code = self.STD_AMINO_ACIDS[res_name]
                        
                        self.sequence_data[chain_id].append({
                            'position': res_num,
                            'residue_name': res_name,
                            'aa_code': aa_code
                        })
                    
                    # Get start and end position numbers for this block
                    start_pos = block_residues[0].id[1]
                    end_pos = block_residues[-1].id[1]
                    
                    # Create the sequence line with position prefix first to get exact length
                    sequence = "".join(self.STD_AMINO_ACIDS[res.resname] for res in block_residues)
                    seq_line = f"{start_pos:>4}-{end_pos:<3} {sequence}"
                    
                    # Calculate the actual prefix length (everything before the sequence)
                    prefix_length = len(seq_line) - len(sequence)
                    
                    # Create the position ruler with dynamic prefix spacing
                    ruler = " " * prefix_length  # Match the actual prefix length
                    for i in range(10, len(block_residues) + 1, 10):
                        # Add the number right-aligned at the 10th position
                        number = str(i)
                        # Calculate current position in ruler
                        current_pos = len(ruler) - prefix_length
                        # We want the number to end at position i
                        spaces_needed = i - current_pos - len(number)
                        ruler += " " * spaces_needed + number
                    
                    # Pad ruler to match sequence length if needed
                    while len(ruler) - prefix_length < len(block_residues):
                        ruler += " "
                    
                    # Print the block
                    self.console.print(f"[bright_blue]{ruler}[/bright_blue]")
                    self.console.print(f"[bold white]{seq_line}[/bold white]")
                    self.console.print()  # Empty line between blocks
        
        if not self.sequence_data:
            self.console.print("[yellow]No protein chains found in the structure.[/yellow]")
            return False
        
        # Show summary
        total_residues = sum(len(chain_data) for chain_data in self.sequence_data.values())
        self.console.print(f"[green]Found {len(self.sequence_data)} protein chain(s) with {total_residues} total residues[/green]")
        
        return True
    
    def _specify_mutations(self) -> bool:
        """Allow user to specify mutations interactively."""
        if not self.sequence_data:
            self.console.print("[yellow]Please display sequences first to see residue numbering.[/yellow]")
            return False
        
        # Check for existing mutations in workspace
        workspace = self.processor.workspace
        existing_mutations = workspace.get("pending_mutations", [])
        
        if existing_mutations:
            self.console.print(f"\n[yellow]Found {len(existing_mutations)} existing mutations in workspace:[/yellow]")
            for mutation in existing_mutations:
                chain_id, res_num, from_aa, to_aa = mutation
                self.console.print(f"  • Chain {chain_id}, Position {res_num}, {from_aa}->{to_aa}")

            # Halo the existing-mutation sites so the user sees spatially
            # what they're about to clear (or keep). Same label as the
            # working pending halo so the visual state is continuous —
            # if the user keeps them, they remain halo'd; if they clear,
            # the cleanup below unhighlights.
            try:
                from proprep.structure_prep.viewer_coordinator import viewer as _viewer
                clauses = [f"(:{m[0]} and {m[1]})" for m in existing_mutations]
                if clauses:
                    _viewer.highlight(
                        " or ".join(clauses),
                        style="halo",
                        color="#e31a1c",
                        label="mutator_pending",
                    )
            except Exception:
                pass

            if confirm_with_context(
                processor=self.processor,
                prompt="Clear existing mutations before specifying new ones?",
                default=True,
                module="Amino Acid Mutator",
                description="Clear existing mutations"
            ):
                workspace.set("pending_mutations", [])
                self.console.print("[green]Existing mutations cleared.[/green]")
                try:
                    from proprep.structure_prep.viewer_coordinator import viewer as _viewer
                    _viewer.unhighlight("mutator_pending")
                except Exception:
                    pass
        
        self.console.print("\n[bold cyan]Specify Amino Acid Mutations[/bold cyan]")
        self.console.print("\n[bold]Format:[/bold] CHAIN:POSITION:FROM->TO")

        # Standard mutation examples
        self.console.print("\n[bold]Standard mutations (MODELLER will model the change):[/bold]")
        self.console.print("  [yellow]A:123:ALA->GLY[/yellow]   - Alanine to Glycine at chain A position 123")
        self.console.print("  [yellow]B:456:PHE->TYR[/yellow]   - Phenylalanine to Tyrosine at chain B position 456")

        # Non-standard mutation examples
        self.console.print("\n[bold]Non-standard mutations (for forcefield-defined residues):[/bold]")
        self.console.print("  [yellow]A:154:PHE->CNF[/yellow]   - Phenylalanine to cyanophenylalanine")
        self.console.print("  [cyan]→ You'll select which atoms to keep (e.g., backbone)[/cyan]")
        self.console.print("  [cyan]→ TLEaP will build missing atoms during force field setup[/cyan]")

        # Commands
        self.console.print("\n[bold]Commands:[/bold]")
        self.console.print("  [green]done[/green]    - Finish entering mutations")
        self.console.print("  [green]clear[/green]   - Clear all pending mutations")
        self.console.print("  [green]help[/green]    - Show additional help")

        # Show available standard amino acids
        aa_list = sorted(list(self.AA_THREE_LETTER.values()))
        self.console.print(f"\n[grey50]Standard amino acids: {', '.join(aa_list)}[/grey50]")
        self.console.print("[grey50]Any other 3-letter code will be treated as non-standard[/grey50]")

        while True:
            mutation_input = prompt_with_context(
                processor=self.processor,
                prompt="\nEnter mutation",
                module="Amino Acid Mutator",
                description="Enter mutation specification",
                options_map={"done": "Finish entering mutations", "help": "Show help", "custom": "Mutation (e.g., A:123:ALA->GLY)"}
            ).strip()
            
            if mutation_input.lower() == 'done':
                break
            elif mutation_input.lower() == 'help':
                self._show_mutation_help()
                continue
            elif mutation_input.lower() == 'clear':
                self.mutations.clear()
                self.console.print("[green]All mutations cleared.[/green]")
                continue
            
            # Parse mutation
            success, mutation = self._parse_mutation(mutation_input)
            if success:
                # Check if this is a non-standard mutation (starts with "NONSTANDARD" marker)
                if isinstance(mutation, tuple) and len(mutation) > 4 and mutation[0] == "NONSTANDARD":
                    # Non-standard mutation: ("NONSTANDARD", chain_id, res_num, from_aa, to_aa, atoms_to_keep)
                    _, chain_id, res_num, from_aa, to_aa, atoms_to_keep = mutation
                    self.nonstandard_mutations.append((chain_id, res_num, from_aa, to_aa, atoms_to_keep))
                    self.console.print(f"[green]✓ Added non-standard mutation: Chain {chain_id}, Position {res_num}, {from_aa}->{to_aa}[/green]")
                    self.console.print(f"[grey50]  Keeping {len(atoms_to_keep)} atoms: {', '.join(atoms_to_keep)}[/grey50]")
                else:
                    # Standard mutation: (chain_id, res_num, from_aa, to_aa)
                    chain_id, res_num, from_aa, to_aa = mutation
                    self.mutations.append(mutation)
                    self.console.print(f"[green]✓ Added standard mutation: Chain {chain_id}, Position {res_num}, {from_aa}->{to_aa}[/green]")
                # Halo every pending mutation site so the user sees the
                # accumulating set as they type. Single stable label so
                # the rep grows in place rather than spawning per-call
                # overlays.
                try:
                    from proprep.structure_prep.viewer_coordinator import viewer as _viewer
                    sites = []
                    for m in self.mutations:
                        sites.append((m[0], m[1]))  # (chain, resid)
                    for m in self.nonstandard_mutations:
                        sites.append((m[0], m[1]))
                    if sites:
                        clauses = [f"(:{c} and {r})" for c, r in sites]
                        _viewer.highlight(
                            " or ".join(clauses),
                            style="halo",
                            color="#e31a1c",
                            label="mutator_pending",
                        )
                except Exception:
                    pass
            else:
                self.console.print(f"[red]Invalid mutation format: {mutation}[/red]")

        # Store mutations in workspace
        total_mutations = len(self.mutations) + len(self.nonstandard_mutations)

        if total_mutations > 0:
            self.console.print(f"\n[green]Total mutations specified: {total_mutations}[/green]")
            if self.mutations:
                self.console.print(f"  Standard: {len(self.mutations)}")
                self._show_mutation_summary()
            if self.nonstandard_mutations:
                self.console.print(f"  Non-standard: {len(self.nonstandard_mutations)}")
                self._show_nonstandard_mutation_summary()

            # Automatically store mutations in workspace
            try:
                workspace = self.processor.workspace
                workspace.set("pending_mutations", self.mutations.copy())
                workspace.set("pending_nonstandard_mutations", self.nonstandard_mutations.copy())

                self.console.print(f"\n[green]✓ Automatically stored mutations in workspace![/green]")
                self.console.print("[cyan]Mutations will be applied when you run Structure Fixer → Process structure[/cyan]")
                self.console.print("[yellow]Note: Non-standard mutations will be applied after MODELLER (or without it)[/yellow]")

            except Exception as e:
                self.console.print(f"\n[red]Error storing mutations: {str(e)}[/red]")
                return False
        else:
            self.console.print("[yellow]No mutations specified.[/yellow]")

        return total_mutations > 0
    
    def _parse_mutation(self, mutation_str: str) -> Tuple[bool, Any]:
        """Parse a mutation string into components.

        Handles both standard and non-standard residue mutations.
        Returns: (success, mutation_data or error_message)
        """
        # Expected format: CHAIN:POSITION:FROM->TO
        # Allow longer residue names for non-standard residues
        pattern = r'^([A-Za-z0-9]):(\d+):([A-Z0-9]{3,})->([A-Z0-9]{3,})$'
        match = re.match(pattern, mutation_str.upper())

        if not match:
            return False, "Invalid format. Use CHAIN:POSITION:FROM->TO (e.g., A:123:PHE->CNF)"

        chain_id, pos_str, from_aa, to_aa = match.groups()
        res_num = int(pos_str)

        # Validate chain exists
        if chain_id not in self.sequence_data:
            return False, f"Chain {chain_id} not found in structure"

        # Validate position exists
        positions = [res['position'] for res in self.sequence_data[chain_id]]
        if res_num not in positions:
            return False, f"Position {res_num} not found in chain {chain_id}"

        # Validate from_aa matches current residue
        current_res = next(res for res in self.sequence_data[chain_id] if res['position'] == res_num)
        if current_res['residue_name'] != from_aa:
            return False, f"Position {res_num} is {current_res['residue_name']}, not {from_aa}"

        # Check for duplicate mutations
        for existing in self.mutations:
            if existing[0] == chain_id and existing[1] == res_num:
                return False, f"Mutation already specified for chain {chain_id} position {res_num}"
        for existing in self.nonstandard_mutations:
            if existing[0] == chain_id and existing[1] == res_num:
                return False, f"Mutation already specified for chain {chain_id} position {res_num}"

        # Determine if this is a standard or non-standard mutation
        is_from_standard = from_aa in self.AA_THREE_LETTER.values()
        is_to_standard = to_aa in self.AA_THREE_LETTER.values()

        if not is_from_standard:
            return False, f"Unknown source residue: {from_aa}. FROM residue must be a standard amino acid."

        if is_to_standard:
            # Standard mutation - validated
            return True, (chain_id, res_num, from_aa, to_aa)
        else:
            # Non-standard mutation - need atom selection
            self.console.print(f"\n[yellow]⚠ Detected non-standard residue target: {to_aa}[/yellow]")
            self.console.print("[cyan]Non-standard mutations require atom selection (atoms to keep from original residue)[/cyan]")

            # Get structure and residue object for atom selection
            workspace = self.processor.workspace
            structure = self._get_structure_from_workspace(workspace, silent=True)
            if not structure:
                return False, "Cannot access structure for atom selection"

            # Get the specific residue
            target_residue = self._get_residue_from_structure(structure, chain_id, res_num)
            if not target_residue:
                return False, f"Cannot find residue {chain_id}:{res_num} in structure"

            # Interactive atom selection
            atoms_to_keep = self._select_atoms_to_keep(target_residue, from_aa, to_aa)
            if atoms_to_keep is None:
                return False, "Atom selection cancelled"

            # Filter to only atoms that actually exist in the residue
            existing_atom_names = {atom.name for atom in target_residue.get_atoms()}
            atoms_to_keep_filtered = [atom for atom in atoms_to_keep if atom in existing_atom_names]

            return True, ("NONSTANDARD", chain_id, res_num, from_aa, to_aa, atoms_to_keep_filtered)
    
    def _show_mutation_help(self):
        """Show detailed help for mutation specification."""
        help_text = """
[bold]Standard Mutations:[/bold]
MODELLER will model the complete structural change, including all atoms.
Best for: Conservative mutations, structure refinement

Example workflow:
  1. Specify: A:123:ALA->GLY
  2. Run Structure Fixer → MODELLER applies mutation
  3. Done! Complete structure with all atoms

[bold]Non-Standard Mutations:[/bold]
For residues defined in your force field but not standard amino acids.
You keep shared atoms, TLEaP builds the rest.
Best for: Modified residues (phospho-, methylated, unnatural AAs)

Example workflow:
  1. Specify: A:154:PHE->CNF (cyanophenylalanine)
  2. Select atoms to keep: "backbone" or "bb+cb"
  3. Run Structure Fixer → Atoms removed, residue renamed
  4. Run TLEaP Generator → Missing atoms built from forcefield

[bold]Important Notes:[/bold]
- Use sequence display (option 1) to see exact residue numbers
- FROM residue must match what's currently at that position
- All residue codes must be 3-letter format (ALA, PHE, CNF, etc.)
- Chain IDs are case-sensitive (A ≠ a)
- Multiple mutations allowed - enter them one at a time

[bold]Common Non-Standard Residues:[/bold]
  CNF - para-cyanophenylalanine
  SEP - phosphoserine
  TPO - phosphothreonine
  PTR - phosphotyrosine
  MLY - N-dimethyl-lysine
  (Requires residue definition in your forcefield)
        """
        self.console.print(Panel(help_text, title="Mutation Help", expand=False))
    
    def _show_mutation_summary(self):
        """Show summary of specified mutations."""
        if not self.mutations:
            return
        
        table = Table(title="Specified Mutations")
        table.add_column("Chain", style="cyan")
        table.add_column("Position", style="magenta")
        table.add_column("From", style="red")
        table.add_column("To", style="green")
        table.add_column("Change Type", style="yellow")
        
        for chain_id, res_num, from_aa, to_aa in self.mutations:
            # Classify mutation type
            change_type = self._classify_mutation(from_aa, to_aa)
            table.add_row(chain_id, str(res_num), from_aa, to_aa, change_type)
        
        self.console.print(table)

    def _show_nonstandard_mutation_summary(self):
        """Show summary of specified non-standard mutations."""
        if not self.nonstandard_mutations:
            return

        table = Table(title="Non-Standard Mutations")
        table.add_column("Chain", style="cyan")
        table.add_column("Position", style="magenta")
        table.add_column("From", style="red")
        table.add_column("To", style="green")
        table.add_column("Atoms to Keep", style="yellow")

        for chain_id, res_num, from_aa, to_aa, atoms_to_keep in self.nonstandard_mutations:
            # Show all atoms, no truncation
            atoms_str = ", ".join(atoms_to_keep)
            table.add_row(chain_id, str(res_num), from_aa, to_aa, atoms_str)

        self.console.print(table)

    def _classify_mutation(self, from_aa: str, to_aa: str) -> str:
        """Classify the type of amino acid change."""
        # More detailed classification based on properties
        hydrophobic = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}
        polar_uncharged = {"SER", "THR", "ASN", "GLN", "TYR", "CYS"}
        positive_charged = {"ARG", "LYS", "HIS"}  # Basic
        negative_charged = {"ASP", "GLU"}  # Acidic
        small = {"GLY", "ALA"}
        
        def get_detailed_class(aa):
            if aa in positive_charged:
                return "positive"
            elif aa in negative_charged:
                return "negative"
            elif aa in hydrophobic:
                return "hydrophobic"
            elif aa in polar_uncharged:
                return "polar"
            elif aa == "GLY":
                return "flexible"
            else:
                return "other"
        
        from_class = get_detailed_class(from_aa)
        to_class = get_detailed_class(to_aa)
        
        # Special cases for very conservative changes
        if from_aa == to_aa:
            return "identical"
        elif (from_aa in small and to_aa in small):
            return "conservative"
        elif (from_aa in {"LEU", "ILE"} and to_aa in {"LEU", "ILE"}):
            return "conservative"
        elif (from_aa in {"ASP", "ASN"} and to_aa in {"ASP", "ASN"}):
            return "conservative"
        elif (from_aa in {"GLU", "GLN"} and to_aa in {"GLU", "GLN"}):
            return "conservative"
        elif (from_aa in {"SER", "THR"} and to_aa in {"SER", "THR"}):
            return "conservative"
        # Charge reversal is always radical
        elif (from_class == "positive" and to_class == "negative") or \
             (from_class == "negative" and to_class == "positive"):
            return "radical"
        # Same general class
        elif from_class == to_class:
            return "moderate"
        else:
            return f"{from_class}→{to_class}"
    
    def _apply_mutations(self) -> bool:
        """Store mutations in workspace for Structure Completeness to apply."""
        if not self.mutations:
            self.console.print("[yellow]No mutations specified. Use 'specify' option first.[/yellow]")
            return False
        
        workspace = self.processor.workspace
        
        self.console.print(f"\n[bold cyan]Storing {len(self.mutations)} mutations for processing...[/bold cyan]")
        self._show_mutation_summary()

        # Re-fire the pending halo right before the confirm prompt so the
        # user is confirming against the visible structural set, not just
        # the table. Hook A keeps this fresh during specify, but a stale
        # path here (e.g. user came directly to apply) needs the refresh.
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer
            sites = [(m[0], m[1]) for m in self.mutations]
            sites.extend((m[0], m[1]) for m in self.nonstandard_mutations)
            if sites:
                clauses = [f"(:{c} and {r})" for c, r in sites]
                _viewer.highlight(
                    " or ".join(clauses),
                    style="halo",
                    color="#e31a1c",
                    label="mutator_pending",
                )
        except Exception:
            pass

        confirm = confirm_with_context(
            processor=self.processor,
            prompt="Store these mutations for application during structure modeling?",
            default=True,
            module="Amino Acid Mutator",
            description="Confirm storing mutations for modeling"
        )
        if not confirm:
            return False
        
        try:
            # Store mutations in workspace for Structure Completeness module
            workspace.set("pending_mutations", self.mutations.copy())            

            self.console.print(f"[green]Successfully stored {len(self.mutations)} mutations![/green]")
            self.console.print("[cyan]Mutations will be applied when Structure Completeness runs MODELLER[/cyan]")
            
            # Keep mutations visible until they're actually applied
            # They'll be cleared by Structure Completeness after successful application
            
            return True
            
        except Exception as e:
            self.console.print(f"[red]Error storing mutations: {str(e)}[/red]")
            return False

    
    def _clear_mutations(self) -> bool:
        """Clear all pending mutations."""
        if not self.mutations:
            self.console.print("[yellow]No mutations to clear.[/yellow]")
            return False

        count = len(self.mutations)
        self.mutations.clear()
        self.console.print(f"[green]Cleared {count} pending mutations.[/green]")
        # Drop the pending halo since the working set is now empty.
        # Non-standard mutations may still be present — re-halo them so
        # the rep tracks the actual remaining state rather than going
        # blank for both kinds.
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer
            if self.nonstandard_mutations:
                clauses = [f"(:{m[0]} and {m[1]})" for m in self.nonstandard_mutations]
                _viewer.highlight(
                    " or ".join(clauses),
                    style="halo",
                    color="#e31a1c",
                    label="mutator_pending",
                )
            else:
                _viewer.unhighlight("mutator_pending")
        except Exception:
            pass
        return True
    
    def _show_mutation_status(self) -> bool:
        """Show current mutation status."""
        self.console.print(f"\n[bold cyan]Mutation Status[/bold cyan]")

        workspace = self.processor.workspace

        # Show pending standard mutations
        if self.mutations:
            self.console.print(f"Pending standard mutations (not yet stored): {len(self.mutations)}")
            self._show_mutation_summary()

        # Show pending non-standard mutations
        if self.nonstandard_mutations:
            self.console.print(f"\nPending non-standard mutations (not yet stored): {len(self.nonstandard_mutations)}")
            self._show_nonstandard_mutation_summary()

        # Show stored standard mutations
        stored_mutations = workspace.get("pending_mutations", [])
        if stored_mutations:
            self.console.print(f"\n[yellow]Stored standard mutations (awaiting Structure Fixer): {len(stored_mutations)}[/yellow]")
            from rich.table import Table
            table = Table(title="Stored Standard Mutations")
            table.add_column("Chain", style="cyan")
            table.add_column("Position", style="magenta")
            table.add_column("From", style="red")
            table.add_column("To", style="green")
            table.add_column("Change Type", style="yellow")

            for chain_id, res_num, from_aa, to_aa in stored_mutations:
                change_type = self._classify_mutation(from_aa, to_aa)
                table.add_row(chain_id, str(res_num), from_aa, to_aa, change_type)

            self.console.print(table)

        # Show stored non-standard mutations
        stored_nonstandard = workspace.get("pending_nonstandard_mutations", [])
        if stored_nonstandard:
            self.console.print(f"\n[yellow]Stored non-standard mutations (awaiting Structure Fixer): {len(stored_nonstandard)}[/yellow]")
            from rich.table import Table
            table = Table(title="Stored Non-Standard Mutations")
            table.add_column("Chain", style="cyan")
            table.add_column("Position", style="magenta")
            table.add_column("From", style="red")
            table.add_column("To", style="green")
            table.add_column("Atoms to Keep", style="yellow")

            for chain_id, res_num, from_aa, to_aa, atoms_to_keep in stored_nonstandard:
                # Show all atoms, no truncation
                atoms_str = ", ".join(atoms_to_keep)
                table.add_row(chain_id, str(res_num), from_aa, to_aa, atoms_str)

            self.console.print(table)

        # Show applied mutations
        applied_mutations = workspace.get("mutations_applied", [])
        if applied_mutations:
            self.console.print(f"\n[green]Successfully applied: {len(applied_mutations)} mutations[/green]")

        # Check if nothing to show
        if (not self.mutations and not self.nonstandard_mutations and
            not stored_mutations and not stored_nonstandard and not applied_mutations):
            self.console.print("[grey50]No mutations specified, stored, or applied[/grey50]")

        # Halo every mutation site visible from this status view —
        # pending working set + stored workspace set, both standard and
        # non-standard. Single label so the rep is "everything the
        # mutator currently knows about". Applied mutations aren't
        # halo'd separately since the structure has changed by then.
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer
            sites = set()
            for m in self.mutations:
                sites.add((m[0], m[1]))
            for m in self.nonstandard_mutations:
                sites.add((m[0], m[1]))
            for m in stored_mutations:
                sites.add((m[0], m[1]))
            for m in stored_nonstandard:
                sites.add((m[0], m[1]))
            if sites:
                clauses = [f"(:{c} and {r})" for c, r in sorted(sites)]
                _viewer.highlight(
                    " or ".join(clauses),
                    style="halo",
                    color="#e31a1c",
                    label="mutator_pending",
                )
            else:
                _viewer.unhighlight("mutator_pending")
        except Exception:
            pass

        return True

    def _get_residue_from_structure(self, structure: Structure, chain_id: str, res_num: int) -> Optional[Residue]:
        """Get a specific residue from the structure."""
        for model in structure:
            if chain_id in model:
                chain = model[chain_id]
                for residue in chain:
                    if residue.id[1] == res_num and residue.id[0] == " ":
                        return residue
        return None

    def _select_atoms_to_keep(self, residue: Residue, from_aa: str, to_aa: str) -> Optional[List[str]]:
        """Interactive atom selection for non-standard mutations.

        Args:
            residue: Bio.PDB.Residue object
            from_aa: Source residue name (3-letter code)
            to_aa: Target residue name (3-letter code)

        Returns:
            List of atom names to keep, or None if cancelled
        """
        # Get all atoms in the residue
        atoms = list(residue.get_atoms())

        if not atoms:
            self.console.print("[red]No atoms found in residue[/red]")
            return None

        self.console.print(f"\n[bold cyan]Atom Selection for {from_aa} → {to_aa} Mutation[/bold cyan]")
        self.console.print(f"Residue: {residue.get_parent().id}:{residue.id[1]}")
        self.console.print(f"\n[yellow]Select which atoms to KEEP (shared between {from_aa} and {to_aa}):[/yellow]")
        self.console.print(f"[grey50]TLEaP will build the missing atoms for {to_aa} during force field setup[/grey50]\n")

        # Phase 1: zoom attention onto the residue being modified — red
        # ball+stick on the target so the user sees what's about to be
        # mutated before picking which atoms to retain. Replaced by the
        # kept-atoms rep (Phase 2) once they pick a preset / custom set.
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer
            chain_id_target = residue.get_parent().id
            res_num_target = residue.id[1]
            _viewer.highlight(
                f":{chain_id_target} and {res_num_target}",
                style="ball+stick",
                color="#e31a1c",
                label="mutator_target",
            )
        except Exception:
            pass

        # Display atoms in a table
        table = Table(title=f"Atoms in {from_aa} {residue.id[1]}")
        table.add_column("Index", style="cyan")
        table.add_column("Atom Name", style="yellow")
        table.add_column("Element", style="green")
        table.add_column("Coords (x, y, z)", style="grey50")

        for idx, atom in enumerate(atoms, 1):
            coords = f"({atom.coord[0]:.2f}, {atom.coord[1]:.2f}, {atom.coord[2]:.2f})"
            table.add_row(str(idx), atom.name, atom.element, coords)

        self.console.print(table)

        # Offer preset options for common cases
        self.console.print("\n[bold]Preset Options:[/bold]")
        self.console.print("  [cyan]backbone[/cyan]  - Keep backbone atoms only (N, CA, C, O)")
        self.console.print("  [cyan]bb+cb[/cyan]     - Keep backbone + CB (beta carbon)")
        self.console.print("  [cyan]all[/cyan]       - Keep all atoms (for testing)")
        self.console.print("  [cyan]custom[/cyan]    - Manually select specific atoms")
        self.console.print("  [cyan]cancel[/cyan]    - Cancel this mutation\n")

        choice = prompt_with_context(
            processor=self.processor,
            prompt="Select preset or 'custom'",
            module="Amino Acid Mutator",
            description="Atom selection method",
            options_map={
                "backbone": "Keep backbone atoms only",
                "bb+cb": "Keep backbone + CB",
                "all": "Keep all atoms",
                "custom": "Custom selection",
                "cancel": "Cancel mutation"
            }
        ).strip().lower()

        if choice == "cancel":
            # Drop both the target and any kept-atoms halo since the user
            # backed out of this mutation.
            try:
                from proprep.structure_prep.viewer_coordinator import viewer as _viewer
                _viewer.unhighlight("mutator_target")
                _viewer.unhighlight("mutator_kept_atoms")
            except Exception:
                pass
            return None
        elif choice == "backbone":
            kept = ["N", "CA", "C", "O", "H"]  # Include H for proper peptide bond
        elif choice == "bb+cb":
            kept = ["N", "CA", "C", "O", "H", "CB"]
        elif choice == "all":
            kept = [atom.name for atom in atoms]
        elif choice == "custom":
            kept = self._custom_atom_selection(atoms)
            if kept is None:
                return None
        else:
            self.console.print(f"[red]Invalid choice: {choice}[/red]")
            return None

        # Phase 2: halo only the kept atoms in green so the user sees
        # which atoms survive the mutation. Replaces the red target rep
        # so the visual transitions from "this whole residue is the
        # target" to "this subset is what stays".
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer
            chain_id_target = residue.get_parent().id
            res_num_target = residue.id[1]
            atom_clause = " or ".join(f".{n}" for n in kept)
            _viewer.highlight(
                f":{chain_id_target} and {res_num_target} and ({atom_clause})",
                style="ball+stick",
                color="#33a02c",
                label="mutator_kept_atoms",
            )
        except Exception:
            pass

        return kept

    def _custom_atom_selection(self, atoms: List) -> Optional[List[str]]:
        """Custom atom selection by index or name.

        Args:
            atoms: List of Bio.PDB.Atom objects

        Returns:
            List of atom names to keep, or None if cancelled
        """
        self.console.print("\n[bold]Custom Atom Selection[/bold]")
        self.console.print("Enter atom indices (comma-separated) or names (comma-separated)")
        self.console.print("Examples: [yellow]1,2,3,4[/yellow] or [yellow]N,CA,C,O,CB[/yellow]\n")

        selection_input = prompt_with_context(
            processor=self.processor,
            prompt="Enter selection",
            module="Amino Acid Mutator",
            description="Custom atom selection"
        ).strip()

        if not selection_input:
            return None

        selected_atoms = []

        # Try parsing as indices first
        try:
            indices = [int(x.strip()) for x in selection_input.split(",")]
            for idx in indices:
                if 1 <= idx <= len(atoms):
                    selected_atoms.append(atoms[idx - 1].name)
                else:
                    self.console.print(f"[red]Invalid index: {idx}[/red]")
                    return None
            return selected_atoms
        except ValueError:
            # Not indices, try as atom names
            names = [x.strip().upper() for x in selection_input.split(",")]
            atom_names = {atom.name for atom in atoms}

            for name in names:
                if name in atom_names:
                    selected_atoms.append(name)
                else:
                    self.console.print(f"[red]Atom not found: {name}[/red]")
                    return None

            return selected_atoms

