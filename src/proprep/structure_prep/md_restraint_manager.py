"""
MD Restraint Manager for any PDB structure
Interactive console interface for setting distance, angle, and torsion restraints
"""

import numpy as np
import logging
from typing import List, Tuple, Dict, Optional, Any, Set
from dataclasses import dataclass, field
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich import box
from rich.text import Text

# ProPrep module imports
from proprep.utils.module_registry import ProcessingModule, register_module
from proprep.utils.prompts import prompt_with_context, confirm_with_context, float_prompt_with_context, int_prompt_with_context
from .md_restraint_commands import (
    ConfigureRestraintsCommand,
    DisplayRestraintsCommand,
    ExportDisangCommand,
    ImportRestraintsCommand,
    GenerateRestraintMaskCommand,
)

logger = logging.getLogger(__name__)


@dataclass
class StructureAtom:
    """Atom from any PDB structure (generalized from RedoxSiteAtom)"""
    chain: str
    resname: str  
    resid: int
    atom_name: str
    coords: Tuple[float, float, float]
    element: str
    altloc: str = ""
    insertion_code: str = ""
    occupancy: Optional[float] = None
    bfactor: Optional[float] = None

@dataclass
class MDRestraint:
    """Individual restraint for MD simulations"""
    restraint_type: str  # 'distance', 'angle', 'torsion'
    atom_coords: List[Tuple[float, float, float]]  # Coordinates of atoms involved
    atom_info: List[Dict[str, Any]]  # Atom identification info
    
    # AMBER flat-bottom potential parameters
    r1: float = 0.0      # Lower bound (parabolic)
    r2: float = 0.0      # Start of flat region
    r3: float = 0.0      # End of flat region  
    r4: float = 999.0    # Upper bound (linear/parabolic)
    rk2: float = 10.0    # Force constant for lower parabola
    rk3: float = 10.0    # Force constant for upper parabola
    
    # Simple parabolic well (alternative to flat-bottom)
    r0: Optional[float] = None  # Target value
    k0: Optional[float] = None  # Force constant
    
    # Optional parameters
    ifvari: int = 0      # Time-varying restraints
    description: str = ""
    active: bool = True
    
    # Computed current values
    current_value: Optional[float] = None
    amber_index: List[int] = field(default_factory=list)  # 1-based atom indices

class RestraintManager:
    """Manager for MD restraints on any PDB structure"""
    
    def __init__(self, console: Console = None, processor=None):
        self.console = console or Console()
        self.restraints: List[MDRestraint] = []
        self.processor = processor  # Reference to main processor for workspace access
        self.current_structure = None
        self.current_structure_source = None
        
    def add_restraints_to_structure(self, workspace=None):
        """Main entry point for adding restraints to any PDB structure"""
        if workspace:
            # Get the best available structure from workspace
            structure_info = self._get_priority_structure(workspace)
            if not structure_info:
                self.console.print("[red]✗ No PDB structure found in workspace[/red]")
                return workspace

            structure_file, structure_source = structure_info
            self.current_structure_source = structure_source

            # Load the structure
            from Bio.PDB import PDBParser
            parser = PDBParser(QUIET=True)
            self.current_structure = parser.get_structure("current", structure_file)

        elif not self.current_structure:
            self.console.print("[red]✗ No structure available for restraint configuration[/red]")
            return workspace if workspace else None
            
        self.console.print(f"\n[bold cyan]═══ MD RESTRAINT CONFIGURATION ═══[/bold cyan]")
        self.console.print(f"Structure source: [bold]{self.current_structure_source}[/bold]")

        # Extract all atoms from structure
        all_atoms = self._extract_structure_atoms()
        self.console.print(f"Available atoms: [bold]{len(all_atoms)}[/bold]")

        # Check for RedoxSites in workspace
        redox_sites = self.get_from_workspace("detected_redox_sites", []) if self.processor else []

        # Display brief summary
        if redox_sites:
            total_bonds = sum(len(site.bonds) for site in redox_sites)
            self.console.print(f"Redox sites detected: [bold]{len(redox_sites)}[/bold] sites with [bold]{total_bonds}[/bold] bonds")

        # Count chains
        chains = set(atom.chain for atom in all_atoms)
        self.console.print(f"Structure chains: [bold]{len(chains)}[/bold] chains ({', '.join(sorted(chains))})")

        # Educational workflow message
        self.console.print("\n[bold cyan]DISANG Restraint Configuration Process:[/bold cyan]")
        self.console.print("  [bold]Step 1:[/bold] Review redox site bonds → Optional: Auto-generate distance restraints from detected bonds")
        self.console.print("  [bold]Step 2:[/bold] Configuration menu → Add custom restraints (distance/angle/torsion), edit, or delete")
        self.console.print("  [bold]Step 3:[/bold] Finish → Select 'done' to save to workspace and export DISANG file")
        self.console.print("\n[cyan]Residue selection modes (available when adding restraints):[/cyan]")
        self.console.print("  • [yellow](s)[/yellow] Sequence view - Browse protein sequences with redox site highlighting")
        self.console.print("  • [yellow](n)[/yellow] Non-standard residues - Numbered list of cofactors, metals, modified residues")
        self.console.print("  • [yellow](r)[/yellow] Redox sites - Select entire sites by ID (ideal for 93 sites!)")
        self.console.print("  • [yellow](m)[/yellow] Manual entry - Direct chain:resid input")
        self.console.print("\n[grey50]'Done' automatically exports DISANG file and saves all data to workspace[/grey50]")

        # Note: Detailed sequence view moved to selection modes

        # Pre-menu RedoxSite bond candidate processing
        if redox_sites:
            total_bonds = sum(len(site.bonds) for site in redox_sites)
            if total_bonds > 0:
                self.console.print()  # Blank line for spacing
                if confirm_with_context(
                    processor=self.processor,
                    prompt="Review RedoxSite bond candidates for distance restraints?",
                    default=False,
                    module="MD Restraint Manager",
                    description="Review RedoxSite bond candidates"
                ):
                    self._process_redox_site_candidates(redox_sites)

        while True:
            action = prompt_with_context(
                processor=self.processor,
                prompt="\n[bold]Restraint Configuration Menu[/bold]\n"
                "[green]add[/green] (a) - Add new restraint\n"
                "[yellow]list[/yellow] (l) - Show current restraints\n"
                "[blue]edit[/blue] (e) - Edit existing restraint\n"
                "[red]delete[/red] (d) - Delete restraint\n"
                "[white]done[/white] - Finish and export restraints\n"
                "Choose action",
                choices=["add", "a", "list", "l", "edit", "e", "delete", "d", "done"],
                default="add",
                module="MD Restraint Manager",
                description="Select restraint action",
                options_map={
                    "add": "Add new restraint", "a": "Add new restraint",
                    "list": "Show current restraints", "l": "Show current restraints",
                    "edit": "Edit existing restraint", "e": "Edit existing restraint",
                    "delete": "Delete restraint", "d": "Delete restraint",
                    "done": "Finish and export restraints"
                }
            )

            # Normalize single-letter inputs
            action_map = {"a": "add", "l": "list", "e": "edit", "d": "delete"}
            action = action_map.get(action, action)

            if action == "add":
                self._add_restraint_interactive(all_atoms, redox_sites)
            elif action == "list":
                self._display_restraints()
            elif action == "edit":
                self._edit_restraint(all_atoms)
            elif action == "delete":
                self._delete_restraint()
            elif action == "done":
                break

        # Save restraints data to workspace if processor available
        if self.processor and workspace:
            self._update_workspace_obj(workspace, "md_restraints", self.restraints)
            self._update_workspace_obj(workspace, "restraint_count", len(self.restraints))
            self._update_workspace_obj(workspace, "restraint_structure_source", self.current_structure_source)

        # Auto-export DISANG file if restraints were defined
        if len(self.restraints) > 0:
            self.console.print(f"\n[cyan]Exporting {len(self.restraints)} restraint(s) to DISANG file...[/cyan]")
            self._export_disang_file()
        else:
            self.console.print(f"\n[yellow]No restraints defined - nothing to export[/yellow]")

        self.console.print(f"\n[green]✓ Restraint configuration completed: {len(self.restraints)} restraints defined[/green]")
        return workspace if workspace else None
        
    # Mapping from workspace key to source name for backward compatibility
    _KEY_TO_SOURCE_NAME = {
        "transformed_pdb_file": "transformed",
        "protonation_pdb_file": "protonation-updated",
        "structure_with_prot_resnames": "protonation-updated",
        "repaired_pdb_file": "repaired",
        "filtered_pdb_file": "filtered",
        "topology_extracted_pdb": "topology-extracted",
        "local_pdb_file": "local",
        "rcsb_pdb_file": "rcsb",
        "alphafold_pdb_file": "alphafold",
        "alphafill_pdb_file": "alphafill",
        "alphafold_homolog_pdb_file": "alphafold-homolog",
        "aligned_target_pdb_file": "aligned-target",
        "aligned_ref_pdb_file": "aligned-ref",
    }

    def _get_priority_structure(self, workspace):
        """
        Get the highest priority structure file from workspace.

        IMPORTANT: For MD Restraint Manager, transformed structures ALWAYS take precedence
        because residue IDs change during transformation but not during protonation updates.

        Uses StructureSelector with custom priority_override to ensure transformed
        structures are checked first, followed by protonation, repaired, etc.
        """
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(workspace, self.console)
        result = selector.get_structure(
            priority_override=[
                "topology_extracted_pdb",    # FIRST - authoritative post-tLEaP source
                "transformed_pdb_file",      # Transformed (also has tLEaP numbering)
                "protonation_pdb_file",
                "structure_with_prot_resnames",  # Legacy support
                "repaired_pdb_file",
                "filtered_pdb_file",
                "hstripped_pdb_file",
                "local_pdb_file",
                "rcsb_pdb_file",
                "alphafold_pdb_file",
                "alphafill_pdb_file",
            ],
            interactive=True,
            return_key=True,
            silent=True,  # We'll display our own messages
        )

        if result is None:
            return None

        structure_file, workspace_key = result

        # Map workspace key to display name for backward compatibility
        source_name = self._KEY_TO_SOURCE_NAME.get(workspace_key, "selected")

        return structure_file, source_name
        
    def _extract_structure_atoms(self):
        """Extract all atoms from current structure"""
        atoms = []
        
        for model in self.current_structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        atoms.append(StructureAtom(
                            chain=chain.id,
                            resname=residue.resname,
                            resid=residue.id[1],
                            atom_name=atom.name,
                            coords=tuple(atom.coord),
                            element=atom.element,
                            insertion_code=residue.id[2] if residue.id[2] != ' ' else '',
                            occupancy=atom.occupancy,
                            bfactor=atom.bfactor
                        ))
        return atoms
        
    def _display_structure_overview(self, atoms, redox_sites=None):
        """Display structure overview with sequence view and RedoxSite coloring"""
        from rich.text import Text
        
        # Standard amino acid residues (including protonation/redox states)
        STD_AMINO_ACIDS = {
            "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
            "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
            "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
            "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
            "MSE": "M",  # Selenomethionine
            # Protonation states
            "ASH": "D",  # Protonated aspartate
            "GLH": "E",  # Protonated glutamate
            "HIE": "H",  # Epsilon-protonated histidine
            "HID": "H",  # Delta-protonated histidine
            "HIP": "H",  # Doubly-protonated histidine
            "HSD": "H",  # Delta-protonated histidine (alt naming)
            "HSE": "H",  # Epsilon-protonated histidine (alt naming)
            "HSP": "H",  # Doubly-protonated histidine (alt naming)
            # Cysteine states
            "CYM": "C",  # Deprotonated cysteine (thiolate)
            "CYX": "C",  # Cysteine in disulfide bridge
            # Lysine states
            "LYN": "K",  # Deprotonated lysine
        }
        
        # Create RedoxSite residue mapping if available
        site_residue_map = {}
        site_colors = []
        if redox_sites:
            # Expanded 256-color palette with ~60 visually distinct bright colors
            # Organized by color families for better visual distinction
            color_indices = [
                # Reds/Magentas
                196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206,
                # Oranges/Yellows
                214, 215, 216, 220, 221, 226, 227, 228,
                # Greens
                46, 47, 48, 49, 50, 51, 82, 83, 84, 85, 86, 118, 119, 120, 121, 122, 154, 155, 156,
                # Cyans/Blues
                51, 45, 44, 43, 42, 39, 33, 27, 21, 57, 63, 69, 75, 81, 87, 93,
                # Purples
                129, 135, 141, 147, 165, 171, 177, 183, 189, 195,
                # Additional bright colors
                226, 190, 154, 118, 82, 46
            ]

            for site_idx, site in enumerate(redox_sites):
                color_idx = color_indices[site_idx % len(color_indices)]
                site_colors.append(f"color({color_idx})")

                # Map residues to this site's color
                for atom in site.atoms:
                    res_key = f"{atom.chain}:{atom.resid}"
                    if res_key not in site_residue_map:
                        site_residue_map[res_key] = (site_idx, f"color({color_idx})")
        
        # Group atoms by chain and residue
        chain_sequences = {}
        for atom in atoms:
            chain_id = atom.chain
            if chain_id not in chain_sequences:
                chain_sequences[chain_id] = {}
            
            res_key = atom.resid
            if res_key not in chain_sequences[chain_id]:
                chain_sequences[chain_id][res_key] = {
                    'resname': atom.resname,
                    'atoms': []
                }
            chain_sequences[chain_id][res_key]['atoms'].append(atom)
        
        self.console.print("\n[bold cyan]═══ PROTEIN SEQUENCES FOR RESTRAINT SELECTION ═══[/bold cyan]")
        self.console.print("Workflow: [bold]1)[/bold] Specify residues → [bold]2)[/bold] Select atoms from those residues")
        self.console.print("Syntax: [yellow]A:123 B:45[/yellow] (distance), [yellow]A:123 A:124 A:125[/yellow] (angle), [yellow]A:123 A:124 A:125 A:126[/yellow] (torsion)")

        if redox_sites:
            site_count = len(redox_sites)

            # Threshold for showing full legend vs summary
            LEGEND_THRESHOLD = 10

            if site_count <= LEGEND_THRESHOLD:
                # Show full legend for small number of sites
                self.console.print(f"\n[cyan]RedoxSite coloring:[/cyan] {site_count} sites detected")
                for i, site in enumerate(redox_sites):
                    color = site_colors[i] if i < len(site_colors) else "white"
                    self.console.print(f"  Site {i+1}: ", style=color, end="")
                    self.console.print(f"{site.site_id}", style=color)
            else:
                # Show summary only for many sites
                from rich.panel import Panel
                info_text = (
                    f"[cyan]RedoxSite coloring:[/cyan] {site_count} sites detected (residues highlighted in sequence)\n\n"
                    f"[grey50]When >10 redox sites are present, the detailed site legend is hidden to keep the display compact.\n"
                    f"Redox site residues are still color-coded in the sequence view below.\n"
                    f"Use the redox site selection mode to select residues by site.[/grey50]"
                )
                panel = Panel(info_text, border_style="cyan", expand=False)
                self.console.print()
                self.console.print(panel)
        
        total_residues = 0
        
        # Create dynamic mapping for non-standard residues to lowercase letters
        nonstandard_mapping = {}
        lowercase_letters = 'abcdefghijklmnopqrstuvwxyz'
        letter_index = 0
        
        # First pass: identify all unique non-standard residues across all chains
        all_nonstandard_residues = set()
        for chain_id in chain_sequences:
            for res_num, res_data in chain_sequences[chain_id].items():
                resname = res_data['resname']
                if resname not in STD_AMINO_ACIDS:
                    all_nonstandard_residues.add(resname)
        
        # Assign letters to non-standard residues in alphabetical order for consistency
        for resname in sorted(all_nonstandard_residues):
            if letter_index < len(lowercase_letters):
                nonstandard_mapping[resname] = lowercase_letters[letter_index]
                letter_index += 1
            else:
                # Fallback to 'z' if we run out of letters (very unlikely)
                nonstandard_mapping[resname] = 'z'

        # Display each chain
        for chain_id in sorted(chain_sequences.keys()):
            residues = chain_sequences[chain_id]
            if not residues:
                continue
                
            # Include ALL residues (both standard and non-standard) and sort by position
            all_residues = []
            for res_num in sorted(residues.keys()):
                res_data = residues[res_num]
                resname = res_data['resname']
                all_residues.append((res_num, resname))
            
            if not all_residues:
                continue
                
            total_residues += len(all_residues)
            self.console.print(f"\n[bold]Chain {chain_id}:[/bold]")
            
            # Display sequence in blocks of 60 residues
            BLOCK_SIZE = 60
            
            for start_idx in range(0, len(all_residues), BLOCK_SIZE):
                end_idx = min(start_idx + BLOCK_SIZE, len(all_residues))
                block_residues = all_residues[start_idx:end_idx]
                
                # Get start and end position numbers for this block
                start_pos = block_residues[0][0]
                end_pos = block_residues[-1][0]
                
                # Create sequence with RedoxSite coloring
                sequence_text = Text()
                for res_num, resname in block_residues:
                    # Get display character (uppercase for standard, lowercase for non-standard)
                    if resname in STD_AMINO_ACIDS:
                        display_char = STD_AMINO_ACIDS[resname]
                    else:
                        display_char = nonstandard_mapping.get(resname, 'x')  # fallback to 'x'
                    
                    res_key = f"{chain_id}:{res_num}"
                    
                    # Apply RedoxSite coloring if available
                    if res_key in site_residue_map:
                        _, color = site_residue_map[res_key]
                        sequence_text.append(display_char, style=color)
                    else:
                        # Use different default colors for standard vs non-standard
                        if resname in STD_AMINO_ACIDS:
                            sequence_text.append(display_char, style="white")
                        else:
                            sequence_text.append(display_char, style="bright_yellow")
                
                # Create position line
                seq_line = Text()
                seq_line.append(f"{start_pos:>4}-{end_pos:<3} ", style="bright_blue")
                seq_line.append(sequence_text)
                
                # Create ruler
                prefix_length = 9  # Length of "1234-567 " format
                ruler = Text(" " * prefix_length, style="bright_blue")
                
                for i in range(10, len(block_residues) + 1, 10):
                    number = str(i)
                    current_pos = len(ruler.plain) - prefix_length
                    spaces_needed = i - current_pos - len(number)
                    ruler.append(" " * spaces_needed + number, style="bright_blue")
                
                # Print the block
                self.console.print(ruler)
                self.console.print(seq_line)
                self.console.print()  # Empty line between blocks
        
        self.console.print(f"[green]Found {len(chain_sequences)} chain(s) with {total_residues} residues[/green]")
        
        # Display legend for non-standard residues if any exist
        if nonstandard_mapping:
            self.console.print(f"\n[bold cyan]Non-standard residue legend:[/bold cyan]")
            legend_items = []
            for resname in sorted(nonstandard_mapping.keys()):
                letter = nonstandard_mapping[resname]
                legend_items.append(f"[bright_yellow]{letter}[/bright_yellow]={resname}")

            # Display legend in rows of 6 items
            for i in range(0, len(legend_items), 6):
                row_items = legend_items[i:i+6]
                self.console.print("  " + "  ".join(row_items))

        # Note: Non-protein residues list removed - users can access via (n) Non-standard residues mode
        # when adding restraints

        return chain_sequences

    def _get_residue_selection_multimode(self, all_atoms, redox_sites=None):
        """
        Multi-mode interface for residue selection.

        Returns:
            str: Residue specification in format "A:123 B:45" etc.
        """
        self.console.print("\n[bold cyan]How would you like to select residues for restraints?[/bold cyan]\n")

        mode_choice = prompt_with_context(
            processor=self.processor,
            prompt="  (s) Sequence view - Browse protein sequences with redox site highlighting\n"
                   "  (n) Non-standard residues - View numbered list of cofactors, metals, modified residues\n"
                   "  (r) Redox sites - Select entire redox sites by ID\n"
                   "  (m) Manual entry - Directly enter chain:resid specifications\n"
                   "Selection mode",
            choices=["s", "n", "r", "m"],
            default="m",
            module="MD Restraint Manager",
            description="Choose residue selection mode",
            options_map={
                "s": "Sequence view",
                "n": "Non-standard residues list",
                "r": "Redox site selection",
                "m": "Manual entry"
            }
        )

        if mode_choice == "s":
            # Display sequence view
            self._display_structure_overview(all_atoms, redox_sites)
            # Then get manual input
            return self._get_manual_residue_input()
        elif mode_choice == "n":
            return self._select_from_nonstandard_residues(all_atoms)
        elif mode_choice == "r":
            return self._select_from_redox_sites(redox_sites, all_atoms)
        else:  # mode_choice == "m"
            return self._get_manual_residue_input()

    def _get_manual_residue_input(self):
        """Get residue specification via manual text input."""
        residue_input = prompt_with_context(
            processor=self.processor,
            prompt="Enter residue specifications directly",
            module="MD Restraint Manager",
            description="Enter residue specifications (e.g., A:123 B:45)"
        )
        return residue_input

    def _select_from_nonstandard_residues(self, all_atoms):
        """Show numbered list of non-standard residues for selection."""
        from rich.table import Table
        from rich import box

        # Standard amino acids (same as in _display_structure_overview)
        STD_AMINO_ACIDS = {
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
            "MSE", "ASH", "GLH", "HIE", "HID", "HIP", "HSD", "HSE", "HSP",
            "CYM", "CYX", "LYN"
        }

        # Collect non-standard residues
        nonstandard_residues = []
        seen = set()

        for atom in all_atoms:
            if atom.resname not in STD_AMINO_ACIDS:
                res_key = (atom.chain, atom.resid, atom.resname)
                if res_key not in seen:
                    seen.add(res_key)
                    nonstandard_residues.append({
                        'chain': atom.chain,
                        'resid': atom.resid,
                        'resname': atom.resname,
                        'chain_resid': f"{atom.chain}:{atom.resid}"
                    })

        if not nonstandard_residues:
            self.console.print("[yellow]No non-standard residues found. Using manual entry mode.[/yellow]")
            return self._get_manual_residue_input()

        # Sort by chain then resid
        nonstandard_residues.sort(key=lambda x: (x['chain'], x['resid']))

        # Display numbered table
        self.console.print("\n[bold cyan]Non-standard Residues Available for Selection:[/bold cyan]\n")

        table = Table(box=box.MINIMAL_DOUBLE_HEAD, show_header=True, header_style="bold cyan")
        table.add_column("#", style="yellow", justify="right", width=5)
        table.add_column("Chain:Resid", style="white", width=15)
        table.add_column("Residue", style="cyan", width=10)

        for idx, res_data in enumerate(nonstandard_residues, 1):
            table.add_row(
                str(idx),
                res_data['chain_resid'],
                res_data['resname']
            )

        self.console.print(table)
        self.console.print(f"\n[grey50]Total: {len(nonstandard_residues)} non-standard residues[/grey50]")

        # Get user selection
        self.console.print("\n[cyan]Enter selection:[/cyan] Numbers (1,2,15), ranges (1-10), or 'all'")
        selection = prompt_with_context(
            self.processor,
            "Selection",
            module="MD Restraint Manager",
            description="Select non-standard residues (numbers, ranges, or 'all')",
        )

        # Parse selection
        selected_indices = self._parse_number_selection(selection, len(nonstandard_residues))

        # Convert to chain:resid format
        selected_specs = []
        for idx in selected_indices:
            res_data = nonstandard_residues[idx - 1]  # Convert to 0-based
            selected_specs.append(res_data['chain_resid'])

        result = ' '.join(selected_specs)
        self.console.print(f"\n[green]Selected: {result}[/green]")
        return result

    def _select_from_redox_sites(self, redox_sites, all_atoms):
        """Show redox sites for selection, with optional filtering."""
        from rich.table import Table
        from rich import box

        if not redox_sites:
            self.console.print("[yellow]No redox sites detected. Using manual entry mode.[/yellow]")
            return self._get_manual_residue_input()

        self.console.print("\n[bold cyan]Redox Sites Available ({} total):[/bold cyan]\n".format(len(redox_sites)))

        # Build site info
        site_info_list = []
        for idx, site in enumerate(redox_sites, 1):
            # Count unique residues in site
            unique_residues = set()
            chains = set()
            for atom in site.atoms:
                unique_residues.add((atom.chain, atom.resid))
                chains.add(atom.chain)

            site_info_list.append({
                'number': idx,
                'site_id': site.site_id,
                'residue_count': len(unique_residues),
                'chains': ','.join(sorted(chains)),
                'residues': sorted(unique_residues),
                'site_obj': site
            })

        # Display table
        table = Table(box=box.MINIMAL_DOUBLE_HEAD, show_header=True, header_style="bold cyan")
        table.add_column("#", style="yellow", justify="right", width=5)
        table.add_column("Site ID", style="white", width=30)
        table.add_column("Chain(s)", style="cyan", width=10)
        table.add_column("Residues", style="green", justify="right", width=10)

        for site_info in site_info_list:
            table.add_row(
                str(site_info['number']),
                site_info['site_id'],
                site_info['chains'],
                str(site_info['residue_count'])
            )

        self.console.print(table)

        # Hook 6: highlight every redox site listed in the table,
        # palette-coloured by row index — the same numbering the user
        # types in below. Lets the user pick spatially rather than by
        # reading site IDs off a table. Single per-site label so a
        # filter step (which re-fires this method) replaces cleanly.
        try:
            viewer = self._viewer_or_none()
            if viewer is not None:
                for lbl in getattr(self, "_redox_site_labels", None) or []:
                    viewer.unhighlight(lbl)
                applied: List[str] = []
                for site_info in site_info_list:
                    idx = site_info['number']
                    pairs = site_info['residues']
                    if not pairs:
                        continue
                    clauses = [f"(:{c} and {r})" for c, r in sorted(pairs)]
                    label = f"{self._VIEWER_LABEL_PREFIX}site_{idx}"
                    viewer.highlight(
                        " or ".join(clauses),
                        style="ball+stick",
                        color=f"palette:{idx}",
                        label=label,
                    )
                    applied.append(label)
                self._redox_site_labels = applied
        except Exception as exc:
            logger.debug("redox-site picker hook silenced: %s", exc)

        # Get user selection
        self.console.print("\n[cyan]Enter selection:[/cyan] Numbers (1,2,15), ranges (1-10), 'all', or 'filter'")
        selection = prompt_with_context(
            self.processor,
            "Selection",
            module="MD Restraint Manager",
            description="Select redox sites (numbers, ranges, 'all', or 'filter')",
        )

        if selection.lower() == 'filter':
            self.console.print("\n[cyan]Filter options:[/cyan]")
            self.console.print("  [yellow]chain:X[/yellow] - Show only sites in chain X")
            self.console.print("  [yellow]type:heme[/yellow] - Show only sites with 'heme' in ID")
            filter_input = prompt_with_context(
                self.processor,
                "Filter",
                module="MD Restraint Manager",
                description="Redox site filter (chain:X or type:STRING)",
            )

            # Apply filter and show filtered list
            filtered_sites = self._filter_redox_sites(site_info_list, filter_input)
            if not filtered_sites:
                self.console.print("[yellow]No sites match filter. Returning to site list.[/yellow]")
                return self._select_from_redox_sites(redox_sites, all_atoms)

            # Re-display filtered sites and get selection
            return self._select_from_filtered_sites(filtered_sites)

        # Parse numeric selection
        selected_indices = self._parse_number_selection(selection, len(site_info_list))

        # Collect all residues from selected sites
        all_selected_residues = set()
        for idx in selected_indices:
            site_info = site_info_list[idx - 1]
            all_selected_residues.update(site_info['residues'])

        # Convert to chain:resid format
        selected_specs = [f"{chain}:{resid}" for chain, resid in sorted(all_selected_residues)]
        result = ' '.join(selected_specs)

        self.console.print(f"\n[green]Selected {len(all_selected_residues)} residues from {len(selected_indices)} site(s)[/green]")
        return result

    def _parse_number_selection(self, selection, max_number):
        """
        Parse selection string like "1,2,5-10,15" into list of indices.

        Returns:
            List of selected indices (1-based)
        """
        if selection.lower() == 'all':
            return list(range(1, max_number + 1))

        selected = set()
        parts = selection.replace(' ', '').split(',')

        for part in parts:
            if '-' in part:
                # Range like "5-10"
                try:
                    start, end = part.split('-')
                    start_num = int(start)
                    end_num = int(end)
                    for num in range(start_num, end_num + 1):
                        if 1 <= num <= max_number:
                            selected.add(num)
                except ValueError:
                    self.console.print(f"[yellow]Warning: Invalid range '{part}' - skipping[/yellow]")
            else:
                # Single number
                try:
                    num = int(part)
                    if 1 <= num <= max_number:
                        selected.add(num)
                    else:
                        self.console.print(f"[yellow]Warning: Number {num} out of range (1-{max_number}) - skipping[/yellow]")
                except ValueError:
                    self.console.print(f"[yellow]Warning: Invalid number '{part}' - skipping[/yellow]")

        return sorted(selected)

    def _filter_redox_sites(self, site_info_list, filter_input):
        """Filter redox sites by chain or type."""
        filter_input = filter_input.lower().strip()

        if filter_input.startswith('chain:'):
            chain_filter = filter_input.split(':', 1)[1].strip().upper()
            return [s for s in site_info_list if chain_filter in s['chains']]
        elif filter_input.startswith('type:'):
            type_filter = filter_input.split(':', 1)[1].strip().lower()
            return [s for s in site_info_list if type_filter in s['site_id'].lower()]
        else:
            self.console.print(f"[yellow]Unknown filter format: {filter_input}[/yellow]")
            return site_info_list

    def _select_from_filtered_sites(self, filtered_sites):
        """Display filtered sites and get selection."""
        from rich.table import Table
        from rich import box

        self.console.print(f"\n[cyan]Filtered Sites ({len(filtered_sites)}):[/cyan]\n")

        table = Table(box=box.MINIMAL_DOUBLE_HEAD, show_header=True, header_style="bold cyan")
        table.add_column("#", style="yellow", justify="right", width=5)
        table.add_column("Site ID", style="white", width=30)
        table.add_column("Chain(s)", style="cyan", width=10)
        table.add_column("Residues", style="green", justify="right", width=10)

        for site_info in filtered_sites:
            table.add_row(
                str(site_info['number']),
                site_info['site_id'],
                site_info['chains'],
                str(site_info['residue_count'])
            )

        self.console.print(table)

        # Get selection from filtered list
        self.console.print("\n[cyan]Enter selection:[/cyan] Numbers (1,2,15), ranges (1-10), or 'all'")
        selection = prompt_with_context(
            self.processor,
            "Selection",
            module="MD Restraint Manager",
            description="Select from filtered redox sites (numbers, ranges, or 'all')",
        )

        # Parse selection (numbers refer to original site numbers)
        selected_indices = self._parse_number_selection(selection, len(filtered_sites))

        # Collect residues
        all_selected_residues = set()
        for idx in selected_indices:
            site_info = filtered_sites[idx - 1]
            all_selected_residues.update(site_info['residues'])

        selected_specs = [f"{chain}:{resid}" for chain, resid in sorted(all_selected_residues)]
        result = ' '.join(selected_specs)

        self.console.print(f"\n[green]Selected {len(all_selected_residues)} residues from {len(selected_indices)} site(s)[/green]")
        return result

    def _add_restraint_interactive(self, all_atoms, redox_sites=None):
        """Interactive restraint addition with RedoxSite context"""
        self.console.print("\n[bold cyan]═══ RESTRAINT SPECIFICATION ═══[/bold cyan]")
        self.console.print("Smart syntax: 2 residues=distance, 3=angle, 4=torsion")
        self.console.print("Examples:")
        self.console.print("  [yellow]A:123 B:45[/yellow] → distance restraint")
        self.console.print("  [yellow]A:123 A:124 A:125[/yellow] → angle restraint")
        self.console.print("  [yellow]A:123 A:124 A:125 A:126[/yellow] → torsion restraint")

        # Get residue selection via multi-mode interface
        residue_input = self._get_residue_selection_multimode(all_atoms, redox_sites)
        
        # Smart parsing: determine restraint type from residue count
        try:
            parts = residue_input.strip().split()
            if len(parts) == 2:
                restraint_type = "distance"
            elif len(parts) == 3:
                restraint_type = "angle"
            elif len(parts) == 4:
                restraint_type = "torsion"
            else:
                self.console.print(f"[red]Invalid: {len(parts)} residues. Need 2, 3, or 4 residues.[/red]")
                return
            
            required_atoms = len(parts)
            
            # Show what type was detected
            type_names = {"distance": "distance", "angle": "angle", "torsion": "torsion"}
            self.console.print(f"[green]→ Detected {type_names[restraint_type]} restraint[/green]")
            
            # Parse residue specification
            residue_specs = self._parse_residue_specification(residue_input, required_atoms)
            if not residue_specs:
                self.console.print("[red]Invalid residue specification[/red]")
                return

        except Exception as e:
            self.console.print(f"[red]Error parsing specification: {e}[/red]")
            return

        # Hook 2: highlight the picked residues so the user can verify
        # the spec ("yes, those are the right two residues") before
        # naming specific atoms.
        self._highlight_picked_residues(residue_specs)

        # Show atoms from selected residues and let user pick specific atoms
        selected_atoms = self._select_atoms_from_residues(all_atoms, residue_specs, restraint_type)
        if not selected_atoms or len(selected_atoms) != required_atoms:
            self.console.print("[red]Atom selection cancelled or incomplete[/red]")
            return

        # Hook 3: narrow the highlight from whole residues down to the
        # specific atoms the user named.
        self._highlight_picked_atoms(selected_atoms)

        # Compute current geometric value and validate
        current_value = self._compute_geometric_value(restraint_type, selected_atoms)
        if current_value is None:
            self.console.print("[red]Failed to compute current geometric value[/red]")
            return
            
        # Distance validation
        if restraint_type == "distance" and current_value > 15.0:
            if not confirm_with_context(
                processor=self.processor,
                prompt=f"[yellow]Warning: Distance is {current_value:.2f}Å (>15Å). Continue?[/yellow]",
                default=False,
                module="MD Restraint Manager",
                description="Confirm large distance"
            ):
                return
        
        # Display current value prominently
        value_units = {"distance": "Å", "angle": "°", "torsion": "°"}[restraint_type]
        self.console.print(f"\n[bold green]Current {restraint_type}: {current_value:.3f} {value_units}[/bold green]")
        
        # Configure restraint parameters
        restraint = self._configure_restraint_parameters(
            restraint_type, selected_atoms, current_value
        )
        
        if restraint:
            # Generate AMBER indices
            restraint.amber_index = self._get_amber_indices(all_atoms, selected_atoms)
            # Stash atom info on the restraint so Hook 4 / Hook 8 can
            # rebuild the NGL selection without re-resolving the StructureAtom list.
            restraint.atom_info = [
                {"chain": a.chain, "resid": a.resid, "atom_name": a.atom_name}
                for a in selected_atoms
            ]
            self.restraints.append(restraint)
            self.console.print(f"[green]✓ {restraint_type.capitalize()} restraint added successfully[/green]")
            # Hook 4: persistent bond overlay for distance restraints
            # (the per-atom highlight from Hook 3 is transient — it
            # gets cleared on the next iteration; the bond overlay
            # accumulates so the viewer shows every restraint added so far).
            self._draw_restraint_overlay(restraint)
    
    def _parse_residue_specification(self, residue_input, required_count):
        """Parse residue specification like 'A:123 B:45' or 'A:123 A:124 A:125'"""
        residue_specs = []
        parts = residue_input.strip().split()
        
        # Note: required_count validation happens in caller now
        for part in parts:
            if ':' not in part:
                raise ValueError(f"Invalid format '{part}'. Use 'CHAIN:RESID'")
            
            chain, resid_str = part.split(':', 1)
            try:
                resid = int(resid_str)
            except ValueError:
                raise ValueError(f"Invalid residue number '{resid_str}'")
            
            residue_specs.append((chain, resid))
        
        return residue_specs
    
    def _select_atoms_from_residues(self, all_atoms, residue_specs, restraint_type):
        """Show focused atom table for selected residues and let user pick atoms"""
        # Group atoms by residue
        residue_atoms = {}
        for atom in all_atoms:
            res_key = (atom.chain, atom.resid)
            if res_key not in residue_atoms:
                residue_atoms[res_key] = []
            residue_atoms[res_key].append(atom)
        
        # Collect atoms from specified residues
        selected_residue_atoms = []
        for chain, resid in residue_specs:
            res_key = (chain, resid)
            if res_key not in residue_atoms:
                self.console.print(f"[red]Residue {chain}:{resid} not found[/red]")
                return None
            selected_residue_atoms.extend(residue_atoms[res_key])
        
        # Display side-by-side columnar atom table
        restraint_names = {"distance": "Distance", "angle": "Angle", "torsion": "Torsion"}
        table = Table(title=f"Atoms for {restraint_names[restraint_type]} Restraint", box=box.ROUNDED)
        
        # Create columns for each residue (in user-specified order)
        residue_atom_lists = []
        for chain, resid in residue_specs:
            res_atoms = []
            res_name = None
            for atom in selected_residue_atoms:
                if atom.chain == chain and atom.resid == resid:
                    res_atoms.append(atom.atom_name)
                    res_name = atom.resname
            
            # Sort atoms in a logical order (backbone first, then sidechain)
            backbone_order = ['N', 'CA', 'C', 'O', 'OXT']
            backbone_atoms = [a for a in res_atoms if a in backbone_order]
            sidechain_atoms = [a for a in res_atoms if a not in backbone_order]
            ordered_atoms = []
            for bb_atom in backbone_order:
                if bb_atom in backbone_atoms:
                    ordered_atoms.append(bb_atom)
            ordered_atoms.extend(sorted(sidechain_atoms))
            
            residue_atom_lists.append((f"{chain}:{resid}({res_name})", ordered_atoms))
            table.add_column(f"{chain}:{resid}({res_name})", style="green", width=15)
        
        # Find the maximum number of atoms in any residue for row count
        max_atoms = max(len(atoms) for _, atoms in residue_atom_lists)
        
        # Fill table rows
        for row_idx in range(max_atoms):
            row_data = []
            for _, atoms in residue_atom_lists:
                if row_idx < len(atoms):
                    row_data.append(atoms[row_idx])
                else:
                    row_data.append("")  # Empty cell if this residue has fewer atoms
            table.add_row(*row_data)
        
        self.console.print("\n", table)
        
        # Get atom selection from user
        required_count = len(residue_specs)
        if restraint_type == "distance":
            self.console.print("Select 2 atoms: [yellow]atom1-atom2[/yellow] (e.g., CA-CB)")
        elif restraint_type == "angle":
            self.console.print("Select 3 atoms: [yellow]atom1-atom2-atom3[/yellow] (e.g., N-CA-C)")
        else:  # torsion
            self.console.print("Select 4 atoms: [yellow]atom1-atom2-atom3-atom4[/yellow] (e.g., N-CA-CB-CG)")
        
        atom_selection = prompt_with_context(
            processor=self.processor,
            prompt="Atoms",
            module="MD Restraint Manager",
            description="Select atoms (e.g., N-CA-C)"
        )
        
        # Parse atom selection
        try:
            atom_names = [name.strip() for name in atom_selection.split('-')]
            if len(atom_names) != required_count:
                self.console.print(f"[red]Expected {required_count} atoms, got {len(atom_names)}[/red]")
                return None
            
            # Find atoms by name in the same order as residue specs
            selected_atoms = []
            for i, (chain, resid) in enumerate(residue_specs):
                atom_name = atom_names[i]
                
                # Find atom in this residue
                found_atom = None
                for atom in selected_residue_atoms:
                    if atom.chain == chain and atom.resid == resid and atom.atom_name == atom_name:
                        found_atom = atom
                        break
                
                if not found_atom:
                    self.console.print(f"[red]Atom '{atom_name}' not found in {chain}:{resid}[/red]")
                    return None
                
                selected_atoms.append(found_atom)
            
            return selected_atoms
            
        except Exception as e:
            self.console.print(f"[red]Error parsing atom selection: {e}[/red]")
            return None
    
    def _process_redox_site_candidates(self, redox_sites):
        """Process RedoxSite bonds as distance restraint candidates"""
        self.console.print("\n[bold cyan]═══ REDOX SITE BOND CANDIDATES ═══[/bold cyan]")
        
        # Extract and organize bonds by site
        all_bond_candidates = []
        bond_index = 1
        
        for site in redox_sites:
            if not site.bonds:
                continue
                
            self.console.print(f"\n[bold]Site: {site.site_id} ({len(site.bonds)} bonds)[/bold]")
            
            # Create table for this site's bonds
            site_table = Table(box=box.ROUNDED)
            site_table.add_column("#", style="cyan", width=4)
            site_table.add_column("Bond", style="green", width=40)
            site_table.add_column("Distance", style="yellow", width=10)
            site_table.add_column("Type", style="blue", width=12)
            site_table.add_column("Treatment", width=12)

            for bond in site.bonds:
                # Find atom information from coordinates
                atom1_info = self._find_atom_by_coords(site.atoms, bond.atom1_coords)
                atom2_info = self._find_atom_by_coords(site.atoms, bond.atom2_coords)
                
                if atom1_info and atom2_info:
                    bond_desc = f"{atom1_info['id']} - {atom2_info['id']}"

                    treatment = getattr(bond, 'treatment', 'bonded')
                    treat_cell = ("[yellow]restrained[/yellow]"
                                  if treatment == 'restrained'
                                  else "[grey50]bonded[/grey50]")
                    site_table.add_row(
                        str(bond_index),
                        bond_desc,
                        f"{bond.distance:.2f}Å",
                        bond.chemical_type,
                        treat_cell
                    )
                    
                    # Store candidate information
                    all_bond_candidates.append({
                        'index': bond_index,
                        'site_id': site.site_id,
                        'bond': bond,
                        'atom1': atom1_info,
                        'atom2': atom2_info,
                        'description': bond_desc
                    })
                    
                    bond_index += 1
            
            self.console.print(site_table)
        
        if not all_bond_candidates:
            self.console.print("[yellow]No valid bond candidates found[/yellow]")
            return
        
        # Bonds the user marked "restrained" during Redox Site Sync are
        # pre-selected: those metal-ligand contacts are held ONLY by a distance
        # restraint (no bonded term), so they need a restraint by construction.
        restrained_default = ",".join(
            str(c['index']) for c in all_bond_candidates
            if getattr(c['bond'], 'treatment', 'bonded') == 'restrained'
        )

        # Get user selection
        self.console.print(f"\n[bold]Select bonds to convert to restraints[/bold]")
        self.console.print("Examples: [yellow]1,3,5[/yellow] or [yellow]1-10[/yellow] or [yellow]all[/yellow]")
        if restrained_default:
            self.console.print(
                f"[grey50]Restrained metal-ligand contacts are pre-selected "
                f"({restrained_default}); press Enter to accept.[/grey50]"
            )

        selection = prompt_with_context(
            processor=self.processor,
            prompt="Bond selection",
            default=restrained_default,
            module="MD Restraint Manager",
            description="Select bonds (e.g., 1,3,5 or all)"
        )
        if not selection.strip():
            self.console.print("[yellow]No bonds selected[/yellow]")
            return
        
        # Parse selection
        selected_indices = self._parse_bond_selection(selection, len(all_bond_candidates))
        if not selected_indices:
            self.console.print("[red]Invalid selection[/red]")
            return
        
        # Create distance restraints from selected bonds
        created_count = self._create_restraints_from_bonds(all_bond_candidates, selected_indices)
        
        if created_count > 0:
            self.console.print(f"\n[green]✓ Created {created_count} distance restraints from RedoxSite bonds[/green]")
            self.console.print("[grey50]Use 'edit' option to customize specific restraint parameters[/grey50]")
        
    def _find_atom_by_coords(self, atoms, target_coords):
        """Find atom information by coordinates"""
        for atom in atoms:
            if (abs(atom.coords[0] - target_coords[0]) < 0.01 and
                abs(atom.coords[1] - target_coords[1]) < 0.01 and 
                abs(atom.coords[2] - target_coords[2]) < 0.01):
                return {
                    'atom': atom,
                    'id': f"{atom.chain}:{atom.resid}({atom.resname}):{atom.atom_name}"
                }
        return None
    
    def _parse_bond_selection(self, selection, max_index):
        """Parse bond selection string into list of indices"""
        selection = selection.strip().lower()
        
        if selection == "all":
            return list(range(1, max_index + 1))
        
        indices = []
        try:
            parts = selection.split(',')
            for part in parts:
                part = part.strip()
                if '-' in part:
                    # Range like "1-5"
                    start, end = part.split('-', 1)
                    start, end = int(start.strip()), int(end.strip())
                    indices.extend(range(start, end + 1))
                else:
                    # Single number
                    indices.append(int(part))
            
            # Filter valid indices
            valid_indices = [i for i in indices if 1 <= i <= max_index]
            return sorted(set(valid_indices))  # Remove duplicates and sort
            
        except ValueError:
            return []
    
    def _create_restraints_from_bonds(self, bond_candidates, selected_indices):
        """Create distance restraints from selected bond candidates"""
        created_count = 0
        
        for index in selected_indices:
            candidate = next((c for c in bond_candidates if c['index'] == index), None)
            if not candidate:
                continue
            
            bond = candidate['bond']
            atom1 = candidate['atom1']['atom']
            atom2 = candidate['atom2']['atom']
            
            # Create restraint with current distance as target and reasonable defaults
            restraint = MDRestraint(
                restraint_type="distance",
                atom_coords=[atom1.coords, atom2.coords],
                atom_info=[
                    {
                        'chain': atom1.chain,
                        'resname': atom1.resname,
                        'resid': atom1.resid,
                        'atom_name': atom1.atom_name,
                        'element': atom1.element
                    },
                    {
                        'chain': atom2.chain,
                        'resname': atom2.resname,
                        'resid': atom2.resid,
                        'atom_name': atom2.atom_name,
                        'element': atom2.element
                    }
                ],
                # Use current distance as target with reasonable flat-bottom parameters
                r2=bond.distance - 0.2,  # Start of flat region
                r3=bond.distance + 0.2,  # End of flat region  
                r1=bond.distance - 0.5,  # Lower bound
                r4=bond.distance + 0.5,  # Upper bound
                rk2=10.0,  # Force constant for lower parabola
                rk3=10.0,  # Force constant for upper parabola
                current_value=bond.distance,
                description=f"RedoxSite bond: {candidate['description']} ({bond.chemical_type})",
                active=True
            )
            
            self.restraints.append(restraint)
            created_count += 1
        
        return created_count
    
    def _select_atoms_from_redox_sites(self, redox_sites, required_count):
        """Select atoms from detected redox sites"""
        selected_atoms = []
        
        while len(selected_atoms) < required_count:
            self.console.print(f"\n[yellow]Select atom {len(selected_atoms) + 1} of {required_count} from redox sites[/yellow]")
            
            # Display redox sites
            site_table = Table(title="Available Redox Sites", box=box.MINIMAL_DOUBLE_HEAD)
            site_table.add_column("#", style="cyan", width=3)
            site_table.add_column("Site ID", style="green", width=20)
            site_table.add_column("Centers", style="yellow", width=8)
            site_table.add_column("Total Atoms", style="white", width=8)
            site_table.add_column("Description", style="blue")
            
            for i, site in enumerate(redox_sites, 1):
                centers_str = f"{len(site.centers)} centers"
                description = f"{', '.join(set(atom.resname for atom in site.atoms))}"
                site_table.add_row(str(i), site.site_id, centers_str, str(len(site.atoms)), description)
            
            self.console.print(site_table)
            
            # Select site
            try:
                site_choice = int_prompt_with_context(
                    processor=self.processor,
                    prompt=f"Select redox site (1-{len(redox_sites)})",
                    default=1,
                    module="MD Restraint Manager",
                    description="Select redox site"
                )
                if 1 <= site_choice <= len(redox_sites):
                    selected_site = redox_sites[site_choice - 1]
                    
                    # Convert RedoxSiteAtoms to StructureAtoms for consistency
                    site_atoms = []
                    for atom in selected_site.atoms:
                        site_atoms.append(StructureAtom(
                            chain=atom.chain,
                            resname=atom.resname,
                            resid=atom.resid,
                            atom_name=atom.atom_name,
                            coords=atom.coords,
                            element=atom.element,
                            insertion_code=atom.insertion_code,
                            occupancy=atom.occupancy,
                            bfactor=atom.bfactor
                        ))
                    
                    # Now select specific atom from this site
                    atom_table = Table(title=f"Atoms in {selected_site.site_id}", box=box.MINIMAL_DOUBLE_HEAD)
                    atom_table.add_column("#", style="cyan", width=4)
                    atom_table.add_column("Residue", style="green", width=12)
                    atom_table.add_column("Atom", style="white", width=8)
                    atom_table.add_column("Element", style="yellow", width=3)
                    
                    available_atoms = []
                    for i, atom in enumerate(site_atoms, 1):
                        if atom not in selected_atoms:
                            atom_table.add_row(
                                str(len(available_atoms) + 1),
                                f"{atom.resname} {atom.chain}:{atom.resid}",
                                atom.atom_name,
                                atom.element
                            )
                            available_atoms.append(atom)
                    
                    self.console.print(atom_table)
                    
                    if available_atoms:
                        atom_choice = int_prompt_with_context(
                            processor=self.processor,
                            prompt=f"Select atom (1-{len(available_atoms)})",
                            default=1,
                            module="MD Restraint Manager",
                            description="Select atom from site"
                        )
                        if 1 <= atom_choice <= len(available_atoms):
                            selected_atom = available_atoms[atom_choice - 1]
                            selected_atoms.append(selected_atom)
                            self.console.print(f"[green]✓ Added {selected_atom.resname} {selected_atom.chain}:{selected_atom.resid} @{selected_atom.atom_name}[/green]")
                        else:
                            self.console.print("[red]✗ Invalid atom number[/red]")
                    else:
                        self.console.print("[yellow]⚠ No more atoms available in this site[/yellow]")
                else:
                    self.console.print("[red]✗ Invalid site number[/red]")
            except Exception as e:
                self.console.print(f"[red]✗ Error: {str(e)}[/red]")
                
        return selected_atoms
    
    def _specify_atoms_by_identifier(self, all_atoms, required_count):
        """Specify atoms using [chain]:[resname]:[resid]:[atom_name] format"""
        selected_atoms = []
        
        self.console.print(f"\n[bold cyan]Specify atoms using format: [chain]:[resname]:[resid]:[atom_name][/bold cyan]")
        self.console.print("[grey50]Example: A:CYS:142:SG or B:HIS:158:NE2[/grey50]")
        
        while len(selected_atoms) < required_count:
            self.console.print(f"\n[yellow]Specify atom {len(selected_atoms) + 1} of {required_count}[/yellow]")
            
            # Show currently selected atoms
            if selected_atoms:
                self.console.print("[grey50]Currently selected:[/grey50]")
                for i, atom in enumerate(selected_atoms, 1):
                    self.console.print(f"  {i}. {atom.chain}:{atom.resname}:{atom.resid}:{atom.atom_name}")
            
            atom_spec = prompt_with_context(
                processor=self.processor,
                prompt="Enter atom specification",
                module="MD Restraint Manager",
                description=f"Specify atom {len(selected_atoms) + 1} of {required_count}"
            )
            
            try:
                # Parse the specification
                parts = atom_spec.strip().split(':')
                if len(parts) != 4:
                    self.console.print("[red]✗ Invalid format. Use chain:resname:resid:atom_name[/red]")
                    continue
                
                chain, resname, resid_str, atom_name = parts
                resid = int(resid_str)
                
                # Find matching atom
                matching_atom = None
                for atom in all_atoms:
                    if (atom.chain == chain and 
                        atom.resname == resname and 
                        atom.resid == resid and 
                        atom.atom_name == atom_name):
                        matching_atom = atom
                        break
                
                if matching_atom:
                    if matching_atom not in selected_atoms:
                        selected_atoms.append(matching_atom)
                        self.console.print(f"[green]✓ Found and added {chain}:{resname}:{resid}:{atom_name}[/green]")
                    else:
                        self.console.print("[yellow]⚠ Atom already selected[/yellow]")
                else:
                    self.console.print(f"[red]✗ Atom {atom_spec} not found in structure[/red]")
                    
                    # Suggest similar atoms
                    suggestions = []
                    for atom in all_atoms[:50]:  # Check first 50 atoms for suggestions
                        if (atom.chain == chain and atom.resname == resname and atom.resid == resid):
                            suggestions.append(atom.atom_name)
                    
                    if suggestions:
                        self.console.print(f"[grey50]Available atoms in {chain}:{resname}:{resid}: {', '.join(set(suggestions))}[/grey50]")
                        
            except ValueError:
                self.console.print("[red]✗ Invalid residue number. Must be an integer[/red]")
            except Exception as e:
                self.console.print(f"[red]✗ Error parsing specification: {str(e)}[/red]")
                
        return selected_atoms
    
    def _select_atoms_interactive(self, all_atoms, required_count):
        """Interactive atom selection from any structure"""
        selected_atoms = []
        
        while len(selected_atoms) < required_count:
            remaining = required_count - len(selected_atoms)
            self.console.print(f"\n[yellow]Select atom {len(selected_atoms) + 1} of {required_count}[/yellow]")
            
            # Option to filter atoms by residue first for large structures
            if len(all_atoms) > 100:
                filter_choice = prompt_with_context(
                    processor=self.processor,
                    prompt="Filter atoms by residue first?",
                    choices=["yes", "no"],
                    default="yes",
                    module="MD Restraint Manager",
                    description="Filter atoms by residue",
                    options_map={
                        "yes": "Filter by residue",
                        "no": "Show all atoms"
                    }
                )
                
                if filter_choice == "yes":
                    filtered_atoms = self._filter_atoms_by_residue(all_atoms)
                    if filtered_atoms:
                        display_atoms = filtered_atoms
                    else:
                        display_atoms = all_atoms
                else:
                    display_atoms = all_atoms
            else:
                display_atoms = all_atoms
            
            # Display available atoms in a table
            atom_table = Table(title="Available Atoms", box=box.MINIMAL_DOUBLE_HEAD)
            atom_table.add_column("#", style="cyan", width=4)
            atom_table.add_column("Residue", style="green", width=12)
            atom_table.add_column("Atom", style="white", width=8)
            atom_table.add_column("Element", style="yellow", width=3)
            atom_table.add_column("Coordinates", style="blue")
            
            # Show up to 50 atoms at a time
            display_count = min(50, len(display_atoms))
            for i, atom in enumerate(display_atoms[:display_count], 1):
                if atom not in selected_atoms:  # Don't show already selected atoms
                    coords_str = f"({atom.coords[0]:.1f}, {atom.coords[1]:.1f}, {atom.coords[2]:.1f})"
                    atom_table.add_row(
                        str(i),
                        f"{atom.resname} {atom.chain}:{atom.resid}",
                        atom.atom_name,
                        atom.element,
                        coords_str
                    )
            
            if len(display_atoms) > 50:
                atom_table.add_row("...", f"({len(display_atoms)-50} more)", "...", "...", "...")
            
            self.console.print(atom_table)
            
            # Show currently selected atoms
            if selected_atoms:
                self.console.print("\n[grey50]Currently selected:[/grey50]")
                for i, atom in enumerate(selected_atoms, 1):
                    self.console.print(f"  {i}. {atom.resname} {atom.chain}:{atom.resid} @{atom.atom_name}")
            
            # Get user selection
            try:
                max_choice = min(50, len([a for a in display_atoms if a not in selected_atoms]))
                if max_choice == 0:
                    self.console.print("[red]✗ No more atoms available to select[/red]")
                    break
                    
                choice = int_prompt_with_context(
                    processor=self.processor,
                    prompt=f"Select atom (1-{max_choice})",
                    default=1,
                    module="MD Restraint Manager",
                    description=f"Select atom {len(selected_atoms) + 1}"
                )
                if 1 <= choice <= max_choice:
                    available_atoms = [a for a in display_atoms[:50] if a not in selected_atoms]
                    if choice <= len(available_atoms):
                        atom = available_atoms[choice - 1]
                        selected_atoms.append(atom)
                        self.console.print(f"[green]✓ Added {atom.resname} {atom.chain}:{atom.resid} @{atom.atom_name}[/green]")
                    else:
                        self.console.print("[red]✗ Invalid atom number[/red]")
                else:
                    self.console.print("[red]✗ Invalid atom number[/red]")
            except Exception as e:
                self.console.print(f"[red]✗ Error: {str(e)}[/red]")
                
        return selected_atoms
        
    def _filter_atoms_by_residue(self, all_atoms):
        """Filter atoms by selecting specific residues first"""
        # Group atoms by residue
        residue_groups = {}
        for atom in all_atoms:
            res_key = f"{atom.resname} {atom.chain}:{atom.resid}"
            if res_key not in residue_groups:
                residue_groups[res_key] = []
            residue_groups[res_key].append(atom)
        
        # Display residue options
        residue_table = Table(title="Available Residues", box=box.MINIMAL_DOUBLE_HEAD)
        residue_table.add_column("#", style="cyan", width=4)
        residue_table.add_column("Residue", style="green", width=15)
        residue_table.add_column("Atoms", style="white", width=8)
        residue_table.add_column("Elements", style="yellow")
        
        residue_list = list(residue_groups.keys())
        for i, res_key in enumerate(residue_list[:20], 1):  # Show first 20
            atoms = residue_groups[res_key]
            elements = sorted(set(atom.element for atom in atoms))
            residue_table.add_row(str(i), res_key, str(len(atoms)), ", ".join(elements))
        
        if len(residue_list) > 20:
            residue_table.add_row("...", f"({len(residue_list)-20} more)", "...", "...")
            
        self.console.print(residue_table)
        
        try:
            choice = int_prompt_with_context(
                processor=self.processor,
                prompt=f"Select residue (1-{min(20, len(residue_list))})",
                module="MD Restraint Manager",
                description="Select residue for atom filtering"
            )
            if 1 <= choice <= min(20, len(residue_list)):
                selected_residue = residue_list[choice - 1]
                return residue_groups[selected_residue]
            else:
                return None
        except:
            return None
    
    def _compute_geometric_value(self, restraint_type, atoms):
        """Compute current geometric value (distance, angle, or torsion)"""
        coords = [np.array(atom.coords) for atom in atoms]
        
        if restraint_type == "distance":
            return self._compute_distance(coords[0], coords[1])
        elif restraint_type == "angle":
            return self._compute_angle(coords[0], coords[1], coords[2])
        elif restraint_type == "torsion":
            return self._compute_torsion(coords[0], coords[1], coords[2], coords[3])
        else:
            raise ValueError(f"Unknown restraint type: {restraint_type}")
    
    def _compute_distance(self, coord1, coord2):
        """Compute distance between two atoms"""
        return np.linalg.norm(coord2 - coord1)
    
    def _compute_angle(self, coord1, coord2, coord3):
        """Compute angle between three atoms (1-2-3, where 2 is vertex)"""
        vec1 = coord1 - coord2  # Vector from vertex to atom 1
        vec2 = coord3 - coord2  # Vector from vertex to atom 3
        
        # Normalize vectors
        vec1_norm = vec1 / np.linalg.norm(vec1)
        vec2_norm = vec2 / np.linalg.norm(vec2)
        
        # Calculate angle using dot product
        cos_angle = np.clip(np.dot(vec1_norm, vec2_norm), -1.0, 1.0)
        angle_rad = np.arccos(cos_angle)
        angle_deg = np.degrees(angle_rad)
        
        return angle_deg
    
    def _compute_torsion(self, coord1, coord2, coord3, coord4):
        """Compute torsion/dihedral angle between four atoms"""
        # Vectors along bonds
        b1 = coord2 - coord1
        b2 = coord3 - coord2  
        b3 = coord4 - coord3
        
        # Normal vectors to planes
        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)
        
        # Normalize normal vectors
        n1_norm = n1 / np.linalg.norm(n1)
        n2_norm = n2 / np.linalg.norm(n2)
        
        # Calculate torsion angle
        cos_torsion = np.clip(np.dot(n1_norm, n2_norm), -1.0, 1.0)
        torsion_rad = np.arccos(cos_torsion)
        
        # Determine sign using scalar triple product
        if np.dot(np.cross(n1_norm, n2_norm), b2 / np.linalg.norm(b2)) < 0:
            torsion_rad = -torsion_rad
            
        torsion_deg = np.degrees(torsion_rad)
        
        return torsion_deg
    
    def _configure_restraint_parameters(self, restraint_type, selected_atoms, current_value):
        """Configure restraint parameters interactively"""
        # Create restraint object with atom info
        atom_info = []
        atom_coords = []
        
        for atom in selected_atoms:
            atom_info.append({
                'chain': atom.chain,
                'resname': atom.resname, 
                'resid': atom.resid,
                'atom_name': atom.atom_name,
                'element': atom.element
            })
            atom_coords.append(atom.coords)
        
        restraint = MDRestraint(
            restraint_type=restraint_type,
            atom_coords=atom_coords,
            atom_info=atom_info,
            current_value=current_value
        )
        
        # Choose parameter configuration method
        param_method = prompt_with_context(
            processor=self.processor,
            prompt="\n[bold]Parameter configuration method[/bold]\n"
            "[green]simple[/green] (s) - Simple parabolic well (target value + force constant)\n"
            "[yellow]flat[/yellow] (f) - Full flat-bottom potential (r1, r2, r3, r4, rk2, rk3)\n"
            "Choose method",
            choices=["simple", "s", "flat", "f"],
            default="simple",
            module="MD Restraint Manager",
            description="Select parameter configuration method",
            options_map={
                "simple": "Simple parabolic well",
                "s": "Simple parabolic well",
                "flat": "Full flat-bottom potential",
                "f": "Full flat-bottom potential"
            }
        )

        if param_method in ["simple", "s"]:
            self._configure_simple_parameters(restraint, current_value)
        else:
            self._configure_advanced_parameters(restraint, current_value)

        # Optional description
        description = prompt_with_context(
            processor=self.processor,
            prompt="\nOptional description for this restraint",
            default="",
            module="MD Restraint Manager",
            description="Enter restraint description"
        )
        restraint.description = description
        
        return restraint
    
    def _configure_simple_parameters(self, restraint, current_value):
        """Configure simple parabolic well parameters"""
        units = {"distance": "Å", "angle": "°", "torsion": "°"}[restraint.restraint_type]
        
        self.console.print(f"\n[bold]Simple parabolic well configuration[/bold]")
        self.console.print(f"Energy = k0 × (value - r0)²")
        
        # Target value (default to current)
        restraint.r0 = float_prompt_with_context(
            processor=self.processor,
            prompt=f"Target value r0 ({units})",
            default=round(current_value, 2),
            module="MD Restraint Manager",
            description=f"Set target {restraint.restraint_type} value"
        )

        # Force constant
        default_k0 = {"distance": 10.0, "angle": 50.0, "torsion": 20.0}[restraint.restraint_type]
        restraint.k0 = float_prompt_with_context(
            processor=self.processor,
            prompt=f"Force constant k0 (kcal/mol-{units}²)",
            default=default_k0,
            module="MD Restraint Manager",
            description="Set force constant"
        )
        
        self.console.print(f"[green]✓ Simple restraint configured: target={restraint.r0:.2f} {units}, k={restraint.k0:.1f}[/green]")
    
    def _configure_advanced_parameters(self, restraint, current_value):
        """Configure advanced flat-bottom potential parameters"""
        units = {"distance": "Å", "angle": "°", "torsion": "°"}[restraint.restraint_type]
        
        self.console.print(f"\n[bold]Flat-bottom potential configuration[/bold]")
        self.console.print("Energy profile:")
        self.console.print("  r < r2:  parabolic with force rk2")
        self.console.print("  r2 ≤ r ≤ r3:  zero energy (allowed region)")  
        self.console.print("  r > r3:  parabolic with force rk3")
        
        # Suggest reasonable defaults based on current value
        tolerance = {"distance": 0.2, "angle": 10.0, "torsion": 20.0}[restraint.restraint_type]
        
        restraint.r2 = float_prompt_with_context(
            processor=self.processor,
            prompt=f"Minimum allowed value r2 ({units})",
            default=round(current_value - tolerance, 2),
            module="MD Restraint Manager",
            description="Set minimum allowed value (r2)"
        )

        restraint.r3 = float_prompt_with_context(
            processor=self.processor,
            prompt=f"Maximum allowed value r3 ({units})",
            default=round(current_value + tolerance, 2),
            module="MD Restraint Manager",
            description="Set maximum allowed value (r3)"
        )

        restraint.r1 = float_prompt_with_context(
            processor=self.processor,
            prompt=f"Lower bound r1 ({units})",
            default=0.0 if restraint.restraint_type == "distance" else round(restraint.r2 - tolerance, 2),
            module="MD Restraint Manager",
            description="Set lower bound (r1)"
        )

        restraint.r4 = float_prompt_with_context(
            processor=self.processor,
            prompt=f"Upper bound r4 ({units})",
            default=999.0 if restraint.restraint_type == "distance" else round(restraint.r3 + tolerance * 2, 2),
            module="MD Restraint Manager",
            description="Set upper bound (r4)"
        )

        # Force constants
        default_k = {"distance": 10.0, "angle": 50.0, "torsion": 20.0}[restraint.restraint_type]

        restraint.rk2 = float_prompt_with_context(
            processor=self.processor,
            prompt=f"Lower force constant rk2 (kcal/mol-{units}²)",
            default=default_k,
            module="MD Restraint Manager",
            description="Set lower force constant (rk2)"
        )

        restraint.rk3 = float_prompt_with_context(
            processor=self.processor,
            prompt=f"Upper force constant rk3 (kcal/mol-{units}²)",
            default=default_k,
            module="MD Restraint Manager",
            description="Set upper force constant (rk3)"
        )
        
        self.console.print(f"[green]✓ Advanced restraint configured: allowed range [{restraint.r2:.2f}, {restraint.r3:.2f}] {units}[/green]")
    
    def _display_restraints(self):
        """Display current restraints in a table"""
        if not self.restraints:
            self.console.print("\n[yellow]No restraints configured[/yellow]")
            return
            
        restraint_table = Table(title=f"Current Restraints ({len(self.restraints)})", box=box.ROUNDED)
        restraint_table.add_column("#", style="cyan", width=3)
        restraint_table.add_column("Type", style="green", width=8)
        restraint_table.add_column("Atoms", style="white", width=25)
        restraint_table.add_column("Current", style="yellow", width=10)
        restraint_table.add_column("Parameters", style="blue", width=35)
        restraint_table.add_column("Status", style="magenta", width=8)
        
        for i, restraint in enumerate(self.restraints, 1):
            # Format atom string
            atoms_str = " → ".join([
                f"{info['resname']}{info['resid']}@{info['atom_name']}"
                for info in restraint.atom_info
            ])
            
            # Format current value
            units = {"distance": "Å", "angle": "°", "torsion": "°"}[restraint.restraint_type]
            current_str = f"{restraint.current_value:.2f} {units}" if restraint.current_value else "N/A"
            
            # Format parameters with restraint type and force constants
            if restraint.r0 is not None:
                # Simple parabolic restraint
                params_str = f"Simple: r0={restraint.r0:.2f}, k={restraint.k0:.1f}"
            else:
                # Flat-bottom restraint - show full range and force constants
                params_str = f"Flat: [{restraint.r2:.1f}-{restraint.r3:.1f}], k={restraint.rk2:.1f}/{restraint.rk3:.1f}"
            
            status_str = "Active" if restraint.active else "Inactive"
            
            restraint_table.add_row(
                str(i),
                restraint.restraint_type.capitalize(),
                atoms_str,
                current_str,
                params_str,
                status_str
            )

        self.console.print("\n", restraint_table)

        # Hook 8: re-render every active restraint in the viewer so
        # the user can compare what's set up against the spatial
        # picture. Distance restraints become yellow lines (one per
        # active restraint), angle/torsion restraints highlight their
        # constituent atoms in palette colours so the user can pick
        # them out. Stable per-restraint labels so this view replaces
        # cleanly on subsequent calls.
        try:
            viewer = self._viewer_or_none()
            if viewer is None:
                return
            # Drop any prior display reps before redrawing.
            for lbl in getattr(self, "_display_restraint_labels", None) or []:
                viewer.unhighlight(lbl)
            applied: List[str] = []
            for i, restraint in enumerate(self.restraints, 1):
                if not restraint.active:
                    continue
                atoms = getattr(restraint, "atom_info", None) or []
                rtype = getattr(restraint, "restraint_type", "")
                if rtype == "distance" and len(atoms) >= 2:
                    sel_a = self._atom_info_to_selection(atoms[0])
                    sel_b = self._atom_info_to_selection(atoms[1])
                    if sel_a and sel_b:
                        label = f"{self._VIEWER_LABEL_PREFIX}displ_bond_{i}"
                        viewer.show_bonds(
                            [(sel_a, sel_b)],
                            label=label,
                            color="#ffff00",
                            show_labels=True,
                        )
                        applied.append(label)
                else:
                    # Angle/torsion: highlight each atom palette-coloured by
                    # its position in the restraint so the geometry is
                    # readable.
                    for j, atom in enumerate(atoms, 1):
                        sel = self._atom_info_to_selection(atom)
                        if not sel:
                            continue
                        label = f"{self._VIEWER_LABEL_PREFIX}displ_atom_{i}_{j}"
                        viewer.highlight(
                            sel,
                            style="ball+stick",
                            color=f"palette:{j}",
                            label=label,
                        )
                        applied.append(label)
            self._display_restraint_labels = applied
        except Exception as exc:
            logger.debug("display-restraints viewer hook silenced: %s", exc)
    
    def _get_amber_indices(self, all_atoms, selected_atoms):
        """Get AMBER 1-based indices for selected atoms"""
        indices = []
        for selected_atom in selected_atoms:
            # Find atom index in all_atoms
            for i, atom in enumerate(all_atoms):
                if (abs(atom.coords[0] - selected_atom.coords[0]) < 1e-6 and
                    abs(atom.coords[1] - selected_atom.coords[1]) < 1e-6 and
                    abs(atom.coords[2] - selected_atom.coords[2]) < 1e-6):
                    indices.append(i + 1)  # AMBER uses 1-based indexing
                    break
        return indices
    
    def _edit_restraint(self, all_atoms):
        """Edit existing restraint(s) with multi-selection support"""
        if not self.restraints:
            self.console.print("\n[yellow]No restraints to edit[/yellow]")
            return
            
        self._display_restraints()
        
        self.console.print(f"\n[bold]Select restraints to edit[/bold]")
        self.console.print("Examples: [yellow]1,3,5[/yellow] or [yellow]1-10[/yellow] or [yellow]all[/yellow] or [yellow]none[/yellow]")
        
        selection = prompt_with_context(
            processor=self.processor,
            prompt="Restraint selection",
            default="1",
            module="MD Restraint Manager",
            description="Select restraints to edit"
        )
        
        # Handle "none" option to cancel editing
        if selection.strip().lower() == "none":
            self.console.print("[yellow]Edit cancelled[/yellow]")
            return
        
        # Parse selection using existing method (reuse bond selection logic)
        selected_indices = self._parse_bond_selection(selection, len(self.restraints))
        if not selected_indices:
            self.console.print("[red]Invalid selection[/red]")
            return
        
        self.console.print(f"[green]Editing {len(selected_indices)} restraint(s)[/green]")
        
        # Edit each selected restraint
        for i, restraint_index in enumerate(selected_indices, 1):
            if 1 <= restraint_index <= len(self.restraints):
                restraint = self.restraints[restraint_index - 1]
                
                self.console.print(f"\n[bold cyan]Editing restraint {restraint_index} of {len(selected_indices)} ({i}/{len(selected_indices)})[/bold cyan]")
                self.console.print(f"Current: {restraint.description}")
                
                # Ask if user wants to edit this one (for batch operations)
                if len(selected_indices) > 1:
                    if not confirm_with_context(
                        processor=self.processor,
                        prompt=f"Edit this restraint?",
                        default=True,
                        module="MD Restraint Manager",
                        description="Confirm edit restraint"
                    ):
                        continue

                edit_choice = prompt_with_context(
                    processor=self.processor,
                    prompt="\n[bold]What to edit?[/bold] [green]parameters[/green] (p), [yellow]description[/yellow] (d), [blue]status[/blue] (s)",
                    choices=["parameters", "p", "description", "d", "status", "s"],
                    default="parameters",
                    module="MD Restraint Manager",
                    description="Select what to edit",
                    options_map={
                        "parameters": "Edit parameters",
                        "p": "Edit parameters",
                        "description": "Edit description",
                        "d": "Edit description",
                        "status": "Edit status",
                        "s": "Edit status"
                    }
                )
                
                # Normalize single-letter inputs
                edit_map = {"p": "parameters", "d": "description", "s": "status"}
                edit_choice = edit_map.get(edit_choice, edit_choice)
                
                if edit_choice == "parameters":
                    # Re-configure parameters
                    current_value = restraint.current_value
                    param_method = prompt_with_context(
                        processor=self.processor,
                        prompt="Parameter type: [green]simple[/green] (s) or [yellow]flat[/yellow] (f)",
                        choices=["simple", "s", "flat", "f"],
                        default="simple" if restraint.r0 is not None else "flat",
                        module="MD Restraint Manager",
                        description="Select parameter type",
                        options_map={
                            "simple": "Simple parabolic well",
                            "s": "Simple parabolic well",
                            "flat": "Flat-bottom potential",
                            "f": "Flat-bottom potential"
                        }
                    )

                    if param_method in ["simple", "s"]:
                        self._configure_simple_parameters(restraint, current_value)
                    else:
                        restraint.r0 = None  # Clear simple parameters
                        restraint.k0 = None
                        self._configure_advanced_parameters(restraint, current_value)

                elif edit_choice == "description":
                    restraint.description = prompt_with_context(
                        processor=self.processor,
                        prompt="New description",
                        default=restraint.description,
                        module="MD Restraint Manager",
                        description="Update restraint description"
                    )

                elif edit_choice == "status":
                    restraint.active = confirm_with_context(
                        processor=self.processor,
                        prompt="Activate this restraint?",
                        default=restraint.active,
                        module="MD Restraint Manager",
                        description="Set restraint active status"
                    )
                
                self.console.print("[green]✓ Restraint updated[/green]")
            else:
                self.console.print(f"[red]✗ Invalid restraint number: {restraint_index}[/red]")
        
        self.console.print(f"[green]✓ Completed editing {len(selected_indices)} restraint(s)[/green]")
    
    def _delete_restraint(self):
        """Delete existing restraint"""
        if not self.restraints:
            self.console.print("\n[yellow]No restraints to delete[/yellow]")
            return
            
        self._display_restraints()
        
        try:
            choice = int_prompt_with_context(
                processor=self.processor,
                prompt=f"Select restraint to delete (1-{len(self.restraints)})",
                default=1,
                module="MD Restraint Manager",
                description="Select restraint to delete"
            )

            if 1 <= choice <= len(self.restraints):
                restraint = self.restraints[choice - 1]

                if confirm_with_context(
                    processor=self.processor,
                    prompt=f"Delete {restraint.restraint_type} restraint?",
                    default=False,
                    module="MD Restraint Manager",
                    description="Confirm delete restraint"
                ):
                    self.restraints.pop(choice - 1)
                    self.console.print("[green]✓ Restraint deleted[/green]")
            else:
                self.console.print("[red]✗ Invalid restraint number[/red]")
                
        except Exception as e:
            self.console.print(f"[red]✗ Error deleting restraint: {str(e)}[/red]")
    
    def _export_disang_file(self):
        """Export restraints to AMBER DISANG format"""
        if not self.restraints:
            self.console.print("\n[yellow]No restraints to export[/yellow]")
            return
        
        active_restraints = [r for r in self.restraints if r.active]
        if not active_restraints:
            self.console.print("\n[yellow]No active restraints to export[/yellow]")
            return
            
        filename = prompt_with_context(
            processor=self.processor,
            prompt="DISANG filename",
            default="restraints.disang",
            module="MD Restraint Manager",
            description="Enter DISANG filename"
        )
        
        try:
            with open(filename, 'w') as f:
                f.write("# AMBER restraint file generated by MD Restraint Manager\n")
                f.write(f"# {len(active_restraints)} restraints defined\n\n")
                
                for i, restraint in enumerate(active_restraints, 1):
                    f.write(f"# Restraint {i}: {restraint.restraint_type}\n")
                    
                    # Write atom info as comments
                    for j, info in enumerate(restraint.atom_info, 1):
                        f.write(f"# Atom {j}: {info['resname']} {info['chain']}:{info['resid']} @{info['atom_name']}\n")
                    
                    if restraint.current_value:
                        units = {"distance": "Å", "angle": "°", "torsion": "°"}[restraint.restraint_type]
                        f.write(f"# Current value: {restraint.current_value:.3f} {units}\n")
                    
                    if restraint.description:
                        f.write(f"# Description: {restraint.description}\n")
                    
                    # Write AMBER namelist
                    f.write("&rst\n")
                    f.write(f"  iat={','.join(map(str, restraint.amber_index))},\n")
                    
                    if restraint.r0 is not None:
                        # Simple parabolic well
                        f.write(f"  r0={restraint.r0:.3f}, k0={restraint.k0:.3f},\n")
                    else:
                        # Flat-bottom potential
                        f.write(f"  r1={restraint.r1:.3f}, r2={restraint.r2:.3f}, ")
                        f.write(f"r3={restraint.r3:.3f}, r4={restraint.r4:.3f},\n")
                        f.write(f"  rk2={restraint.rk2:.3f}, rk3={restraint.rk3:.3f},\n")
                    
                    if restraint.ifvari > 0:
                        f.write(f"  ifvari={restraint.ifvari},\n")
                    
                    f.write("/\n\n")
            
            # Convert to absolute path
            import os
            abs_filename = os.path.abspath(filename)
            
            self.console.print(f"[green]✓ DISANG file exported: {filename}[/green]")
            self.console.print(f"  {len(active_restraints)} active restraints written")
            
            # Save file path to workspace following structure_completeness pattern
            if self.processor:
                self.update_workspace("disang_file", abs_filename)
                self.update_workspace("disang_export_results", {
                    "file_path": abs_filename,
                    "active_restraints": len(active_restraints),
                    "total_restraints": len(self.restraints),
                    "export_successful": True,
                    "source": "restraint_manager"
                })
                self.console.print(f"[grey50]Saved DISANG file path to workspace: {abs_filename}[/grey50]")
            
        except Exception as e:
            self.console.print(f"[red]✗ Error exporting DISANG file: {str(e)}[/red]")
            
            # Save error info to workspace
            if self.processor:
                self.update_workspace("disang_export_results", {
                    "file_path": None,
                    "export_successful": False,
                    "error_message": str(e),
                    "source": "restraint_manager"
                })

    # ========================================================================
    # tLEaP Numbering and Restraint Mask Generation Helper Methods
    # ========================================================================

    def _build_tleap_numbering_map(self, structure_file: str) -> Tuple[Dict, Dict]:
        """
        Parse structure and create mapping from (chain, resid) → tleap_resid.
        tLEaP numbers residues consecutively across all chains.

        Returns:
            chain_info: {chain: {'start': int, 'end': int, 'count': int, 'ranges': [(start, end), ...]}}
            tleap_map: {(chain, resid): tleap_resid}
        """
        from collections import OrderedDict

        chain_info = OrderedDict()
        tleap_map = {}

        # Parse PDB to get residues by chain
        # Note: PDB residue numbers can wrap (max 9999), so we track by sequential
        # appearance rather than unique residue numbers
        # Water residue names to exclude from restraint display
        WATER_RESNAMES = {'WAT', 'HOH', 'TIP', 'TIP3', 'TIP4', 'SPC', 'SOL', 'T3P', 'T4P', 'T4E', 'TP3', 'TP4', 'TP5'}

        chain_residues = OrderedDict()  # {chain: [(resid, is_water), ...]}
        chain_last_resid = {}  # Track last residue number per chain to detect changes

        try:
            with open(structure_file, 'r') as f:
                for line in f:
                    if not line.startswith(('ATOM  ', 'HETATM')):
                        continue

                    chain = line[21:22].strip()
                    if not chain:
                        chain = 'A'  # Default chain if missing

                    resname = line[17:20].strip().upper()
                    is_water = resname in WATER_RESNAMES

                    try:
                        resid = int(line[22:26].strip())
                    except ValueError:
                        continue

                    if chain not in chain_residues:
                        chain_residues[chain] = []
                        chain_last_resid[chain] = None

                    # Add residue when number changes (handles both normal and wrapped numbering)
                    if chain_last_resid[chain] != resid:
                        chain_residues[chain].append((resid, is_water))
                        chain_last_resid[chain] = resid

        except Exception as e:
            logger.error(f"Error parsing structure for tLEaP numbering: {e}")
            return {}, {}

        # Build tLEaP consecutive numbering
        tleap_counter = 1

        for chain in chain_residues:
            residue_tuples = chain_residues[chain]  # List of (resid, is_water) tuples
            if not residue_tuples:
                continue

            # Separate water and non-water residues
            resids = [r[0] for r in residue_tuples]
            non_water_indices = [i for i, (resid, is_water) in enumerate(residue_tuples) if not is_water]
            water_count = sum(1 for _, is_water in residue_tuples if is_water)
            non_water_count = len(residue_tuples) - water_count

            count = len(residue_tuples)

            # Detect if PDB residue numbers wrap (more residues than max resid difference)
            # PDB format supports max 4-digit residue numbers (0-9999)
            min_resid = min(resids)
            max_resid = max(resids)
            resid_span = max_resid - min_resid + 1
            numbers_wrap = count > resid_span

            # Build ranges for non-water residues only (what users care about for restraints)
            ranges = []
            if non_water_count > 0:
                non_water_resids = [resids[i] for i in non_water_indices]
                if numbers_wrap:
                    # Residue numbers wrap - just note first/last non-water PDB numbers
                    ranges.append((non_water_resids[0], non_water_resids[-1]))
                else:
                    # Normal case - detect gaps to identify separate ranges
                    current_range_start = non_water_resids[0]
                    current_range_end = non_water_resids[0]

                    for i in range(1, len(non_water_resids)):
                        gap = non_water_resids[i] - non_water_resids[i-1]
                        if gap > 1 or gap < 0:
                            ranges.append((current_range_start, current_range_end))
                            current_range_start = non_water_resids[i]
                            current_range_end = non_water_resids[i]
                        else:
                            current_range_end = non_water_resids[i]

                    ranges.append((current_range_start, current_range_end))

            # Map each PDB residue to tLEaP number and track non-water tLEaP range
            # Note: For wrapped numbering, same PDB resid maps to multiple tLEaP numbers
            # We store the FIRST occurrence for each (chain, resid) pair
            tleap_start = tleap_counter
            non_water_tleap_start = None
            non_water_tleap_end = None

            for i, (resid, is_water) in enumerate(residue_tuples):
                if (chain, resid) not in tleap_map:
                    tleap_map[(chain, resid)] = tleap_counter
                if not is_water:
                    if non_water_tleap_start is None:
                        non_water_tleap_start = tleap_counter
                    non_water_tleap_end = tleap_counter
                tleap_counter += 1

            tleap_end = tleap_counter - 1

            chain_info[chain] = {
                'start': tleap_start,
                'end': tleap_end,
                'count': count,
                'non_water_count': non_water_count,
                'water_count': water_count,
                'non_water_tleap_start': non_water_tleap_start or tleap_start,
                'non_water_tleap_end': non_water_tleap_end or tleap_end,
                'ranges': ranges,  # List of (start, end) tuples for non-water only
                'numbers_wrap': numbers_wrap  # Flag if PDB numbers cycle
            }

        return chain_info, tleap_map

    def _display_tleap_numbering_info(self, chain_info: Dict, is_transformed: bool = False, structure_file: str = "", structure_source: str = ""):
        """Display educational message about tLEaP numbering with chain composition."""
        from rich.panel import Panel
        from pathlib import Path

        # Build info text
        info_lines = []

        # Show which structure is being used
        if structure_file:
            structure_name = Path(structure_file).name
            info_lines.append(f"[bold cyan]Structure:[/bold cyan] {structure_name}")
            info_lines.append(f"[bold cyan]Source:[/bold cyan] {structure_source}\n")

        if is_transformed:
            info_lines.append("[bold]AMBER Consecutive Residue Numbering[/bold]\n")
            info_lines.append("Residues are numbered consecutively across all chains,")
            info_lines.append("matching the numbering used by AMBER for restraint masks.\n")
            info_lines.append("[bold cyan]Structure composition:[/bold cyan]\n")
        else:
            info_lines.append("[bold]AMBER tLEaP Residue Numbering[/bold]\n")
            info_lines.append("AMBER's tLEaP renumbers residues consecutively across all chains,")
            info_lines.append("ignoring PDB chain IDs. This tool generates restraint masks")
            info_lines.append("compatible with tLEaP numbering.\n")
            info_lines.append("[bold cyan]Your structure composition:[/bold cyan]\n")

        total_non_water = 0
        total_water = 0
        any_numbers_wrap = False

        for chain, info in chain_info.items():
            ranges = info.get('ranges', [])  # Non-water ranges only
            numbers_wrap = info.get('numbers_wrap', False)
            non_water_count = info.get('non_water_count', info['count'])
            water_count = info.get('water_count', 0)
            non_water_tleap_start = info.get('non_water_tleap_start', info['start'])
            non_water_tleap_end = info.get('non_water_tleap_end', info['end'])

            if numbers_wrap:
                any_numbers_wrap = True

            # Use non-water tLEaP range for restraint-relevant display
            tleap_range = f"{non_water_tleap_start}-{non_water_tleap_end}" if non_water_tleap_start != non_water_tleap_end else str(non_water_tleap_start)

            if is_transformed:
                # For transformed structures, show actual ranges intelligently
                if len(ranges) == 1:
                    # Single contiguous range
                    start, end = ranges[0]
                    if start == end:
                        range_str = str(start)
                    else:
                        range_str = f"{start}-{end}"
                    info_lines.append(f"  Chain {chain}: residues {range_str:>15s} ({non_water_count:>4d} res)")
                elif len(ranges) > 1:
                    # Multiple ranges (e.g., protein + cofactors)
                    range_strs = []
                    for start, end in ranges:
                        if start == end:
                            range_strs.append(str(start))
                        else:
                            range_strs.append(f"{start}-{end}")
                    combined_ranges = ", ".join(range_strs)
                    info_lines.append(f"  Chain {chain}: residues {combined_ranges:>30s} ({non_water_count:>4d} res)")
                else:
                    # No non-water residues in this chain (water-only chain)
                    pass
            elif numbers_wrap:
                # PDB residue numbers wrap (large system) - just show tLEaP numbering for non-water
                if non_water_count > 0:
                    info_lines.append(f"  Chain {chain}: {non_water_count:>6d} residues  →  tLEaP residues {tleap_range}")
            else:
                # For non-transformed structures, show PDB → tLEaP mapping
                if len(ranges) == 1:
                    start, end = ranges[0]
                    if start == end:
                        pdb_range = str(start)
                    else:
                        pdb_range = f"{start}-{end}"
                    info_lines.append(f"  Chain {chain}: residues {pdb_range:>15s} ({non_water_count:>4d} res)  "
                                    f"→  tLEaP residues {tleap_range}")
                elif len(ranges) > 1:
                    # Multiple ranges in PDB
                    range_strs = []
                    for start, end in ranges:
                        if start == end:
                            range_strs.append(str(start))
                        else:
                            range_strs.append(f"{start}-{end}")
                    pdb_range = ", ".join(range_strs)
                    info_lines.append(f"  Chain {chain}: residues {pdb_range:>15s} ({non_water_count:>4d} res)  "
                                    f"→  tLEaP residues {tleap_range}")
                # If no ranges, skip (water-only chain)

            total_non_water += non_water_count
            total_water += water_count

        # Display total with water count shown separately
        if total_water > 0:
            info_lines.append(f"\n[bold]Total: {len(chain_info)} chains, {total_non_water} residues[/bold] [grey50](+ {total_water} water)[/grey50]")
        else:
            info_lines.append(f"\n[bold]Total: {len(chain_info)} chains, {total_non_water} residues[/bold]")

        if any_numbers_wrap:
            info_lines.append("\n[grey50]Note: PDB residue numbers cycle (>9999 residues). Use tLEaP numbering for restraints.[/grey50]")

        if is_transformed:
            info_lines.append("\n[grey50]Note: Residue numbers can be used directly in restraint masks[/grey50]")
        else:
            info_lines.append("\n[grey50]Note: Restraint masks will use tLEaP consecutive numbering[/grey50]")

        panel = Panel(
            "\n".join(info_lines),
            title="[bold cyan]Structure Numbering Information[/bold cyan]",
            border_style="cyan",
            expand=False
        )

        self.console.print()
        self.console.print(panel)

    def _configure_general_restraints(self, chain_info: Dict, tleap_map: Dict) -> Optional[Dict]:
        """
        Configure general structural restraints (not redox-specific).

        Args:
            chain_info: Chain composition information
            tleap_map: Mapping from (chain, resid) to tleap_resid

        Returns:
            Dict with keys: 'atom_selection', 'chains', 'tleap_ranges', 'pdb_format', 'tleap_format'
            or None if user cancels
        """
        from rich.panel import Panel

        self.console.print("\n[bold cyan]═══ General Structural Restraints ═══[/bold cyan]\n")

        # Step 1: Choose atom types
        self.console.print("Choose atom types to restrain:")
        self.console.print("  [cyan](b)[/cyan] Backbone atoms (@CA,C,O,N)")
        self.console.print("  [cyan](h)[/cyan] Heavy atoms only (!@H=)")
        self.console.print("  [cyan](c)[/cyan] Alpha carbons only (@CA)")
        self.console.print("  [cyan](u)[/cyan] Custom atom selection (you specify)")

        atom_choice = prompt_with_context(
            processor=self.processor,
            prompt="Select atom type",
            choices=["b", "h", "c", "u"],
            default="c",
            module="MD Restraint Manager",
            description="Select atom types for restraints",
            options_map={
                "b": "Backbone atoms (@CA,C,O,N)",
                "h": "Heavy atoms (!@H=)",
                "c": "Alpha carbons (@CA)",
                "u": "Custom atom selection"
            }
        )

        if atom_choice == "b":
            atom_selection = "@CA,C,O,N"
            atom_description = "Backbone atoms"
        elif atom_choice == "h":
            atom_selection = "!@H="
            atom_description = "Heavy atoms"
        elif atom_choice == "c":
            atom_selection = "@CA"
            atom_description = "Alpha carbons"
        else:  # custom
            self.console.print("\n[grey50]Enter the complete AMBER mask syntax. This will be used directly without modification.[/grey50]")
            atom_selection = prompt_with_context(
                self.processor,
                "Enter complete restraint mask (AMBER syntax)",
                module="MD Restraint Manager",
                description="Custom AMBER restraint mask",
            )
            atom_description = f"Custom: {atom_selection}"
            # Custom mask is used as-is - skip chain/range processing
            return {
                'atom_selection': atom_selection,
                'atom_description': atom_description,
                'chains': [],
                'pdb_ranges': [],
                'tleap_ranges': '',
                'pdb_format': atom_selection,
                'tleap_format': atom_selection,  # Use directly without modification
                'is_custom_mask': True  # Flag to skip further processing
            }

        # Step 2: Choose chains/ranges
        self.console.print(f"\n[bold]Apply to which chains?[/bold]\n")

        # Display available chains (show non-water residue counts for restraint relevance)
        chain_display = []
        for chain, info in chain_info.items():
            non_water_count = info.get('non_water_count', info['count'])
            if non_water_count > 0:  # Only show chains with non-water residues
                chain_display.append(f"[{chain}] Chain {chain} ({non_water_count} residues)")

        self.console.print("Available chains:")
        for line in chain_display:
            self.console.print(f"  {line}")
        self.console.print("  [All] All chains")

        self.console.print("\n[grey50]Enter chain IDs (e.g., 'A,B,C'), 'all', or PDB ranges (e.g., 'A:10-50,B:10-50')[/grey50]")

        chain_input = prompt_with_context(
            self.processor,
            "Chain selection",
            default="all",
            module="MD Restraint Manager",
            description="Chain selection for restraint mask",
        )

        # Parse chain input
        selected_chains = []
        pdb_ranges = []

        if chain_input.lower() == "all":
            selected_chains = list(chain_info.keys())
            # Use full ranges for all chains
            for chain in selected_chains:
                info = chain_info[chain]
                ranges = info['ranges']
                # If single range, use it directly
                if len(ranges) == 1:
                    pdb_ranges.append(f"{chain}:{ranges[0][0]}-{ranges[0][1]}")
                else:
                    # Multiple ranges - add each separately
                    for start, end in ranges:
                        pdb_ranges.append(f"{chain}:{start}-{end}")

        elif ':' in chain_input:
            # Range format: A:10-50,B:10-50
            pdb_ranges = [r.strip() for r in chain_input.split(',')]
            # Extract chains
            for range_spec in pdb_ranges:
                chain = range_spec.split(':')[0].strip()
                if chain in chain_info:
                    selected_chains.append(chain)

        else:
            # Simple chain list: A,B,C
            selected_chains = [c.strip() for c in chain_input.split(',')]
            # Use full ranges
            for chain in selected_chains:
                if chain in chain_info:
                    info = chain_info[chain]
                    ranges = info['ranges']
                    # If single range, use it directly
                    if len(ranges) == 1:
                        pdb_ranges.append(f"{chain}:{ranges[0][0]}-{ranges[0][1]}")
                    else:
                        # Multiple ranges - add each separately
                        for start, end in ranges:
                            pdb_ranges.append(f"{chain}:{start}-{end}")

        # Hook 5: highlight the residue ranges the user selected so
        # they can see what the general mask covers before continuing
        # to the redox-site step / preview. NGL supports range syntax
        # (`:A and 10-50`) directly. Single label per chain — re-firing
        # the prompt replaces the prior selection.
        try:
            viewer = self._viewer_or_none()
            if viewer is not None:
                # Drop any prior general-mask labels (track on self).
                for lbl in getattr(self, "_general_mask_labels", None) or []:
                    viewer.unhighlight(lbl)
                applied: List[str] = []
                for idx, range_spec in enumerate(pdb_ranges, 1):
                    if ':' not in range_spec:
                        continue
                    chain, rng = range_spec.split(':', 1)
                    chain = chain.strip()
                    rng = rng.strip()
                    if not chain or not rng:
                        continue
                    label = f"{self._VIEWER_LABEL_PREFIX}genmask_{chain}_{idx}"
                    viewer.highlight(
                        f":{chain} and {rng}",
                        style="ball+stick",
                        color=f"palette:{idx}",
                        label=label,
                    )
                    applied.append(label)
                self._general_mask_labels = applied
        except Exception as exc:
            logger.debug("general mask viewer hook silenced: %s", exc)

        # Convert chain:range format to AMBER mask format
        tleap_ranges = self._convert_pdb_ranges_to_tleap(pdb_ranges, tleap_map)

        # Build mask components
        pdb_format = f"{atom_selection} in chains {','.join(selected_chains)}"
        tleap_format = f"{atom_selection}&{tleap_ranges}"

        return {
            'atom_selection': atom_selection,
            'atom_description': atom_description,
            'chains': selected_chains,
            'pdb_ranges': pdb_ranges,
            'tleap_ranges': tleap_ranges,
            'pdb_format': pdb_format,
            'tleap_format': tleap_format,
        }

    def _convert_pdb_ranges_to_tleap(self, pdb_ranges: List[str], tleap_map: Dict) -> str:
        """
        Convert PDB ranges like ['A:10-50', 'B:10-50'] to tLEaP format ':10-50,260-300'.

        Args:
            pdb_ranges: List of PDB range specifications
            tleap_map: Mapping from (chain, resid) to tleap_resid

        Returns:
            tLEaP range string
        """
        tleap_residues = []

        for range_spec in pdb_ranges:
            if ':' not in range_spec:
                continue

            chain, res_range = range_spec.split(':', 1)
            chain = chain.strip()
            res_range = res_range.strip()

            if '-' in res_range:
                # Range format: 10-50
                start_str, end_str = res_range.split('-')
                try:
                    start_pdb = int(start_str.strip())
                    end_pdb = int(end_str.strip())

                    # Convert to tLEaP
                    start_tleap = tleap_map.get((chain, start_pdb))
                    end_tleap = tleap_map.get((chain, end_pdb))

                    if start_tleap and end_tleap:
                        tleap_residues.append(f"{start_tleap}-{end_tleap}")

                except ValueError:
                    logger.warning(f"Invalid range format: {range_spec}")
                    continue

            else:
                # Single residue
                try:
                    resid = int(res_range)
                    tleap_resid = tleap_map.get((chain, resid))
                    if tleap_resid:
                        tleap_residues.append(str(tleap_resid))
                except ValueError:
                    continue

        # Collapse contiguous/overlapping ranges to avoid overflowing
        # AMBER's mask parser stack (limited evaluation stack size).
        # E.g. "1-198,199-209,210-407" → "1-407"
        tleap_residues = self._collapse_ranges(tleap_residues)

        return ':' + ','.join(tleap_residues) if tleap_residues else ''

    @staticmethod
    def _collapse_ranges(ranges: List[str]) -> List[str]:
        """Merge contiguous or overlapping residue ranges.

        Converts ['1-198', '199-209', '210-407'] → ['1-407'].
        Single residues like '5' are treated as '5-5'.
        """
        # Parse each range into (start, end) integers
        parsed = []
        for r in ranges:
            if '-' in r:
                parts = r.split('-', 1)
                parsed.append((int(parts[0]), int(parts[1])))
            else:
                val = int(r)
                parsed.append((val, val))

        if not parsed:
            return []

        # Sort by start then end
        parsed.sort()

        # Merge contiguous/overlapping intervals
        merged = [parsed[0]]
        for start, end in parsed[1:]:
            prev_start, prev_end = merged[-1]
            if start <= prev_end + 1:
                # Contiguous or overlapping — extend
                merged[-1] = (prev_start, max(prev_end, end))
            else:
                merged.append((start, end))

        # Format back to strings
        result = []
        for start, end in merged:
            if start == end:
                result.append(str(start))
            else:
                result.append(f"{start}-{end}")
        return result

    def _select_atoms_interactive(self, atom_info_map: Dict, bonded_atoms: set) -> set:
        """
        Display table of unique atoms and allow user to select which to restrain.

        Args:
            atom_info_map: Dict mapping atom_name -> {count, elements}
            bonded_atoms: Set of atom names that are involved in redox site bonds

        Returns:
            Set of selected atom names
        """
        from rich.table import Table
        from rich import box

        self.console.print("\n[bold cyan]Atom Selection[/bold cyan]\n")

        # Create table showing all unique atoms
        table = Table(
            title="Atoms in Redox Sites",
            box=box.MINIMAL_DOUBLE_HEAD,
            show_header=True,
            header_style="bold cyan"
        )

        table.add_column("Atom Name", style="white", justify="left")
        table.add_column("Element(s)", style="cyan", justify="left")
        table.add_column("Count", style="yellow", justify="right")
        table.add_column("In Bonds", style="green", justify="center")

        # Sort atoms: bonded first, then by name
        sorted_atoms = sorted(
            atom_info_map.keys(),
            key=lambda x: (x not in bonded_atoms, x)
        )

        for atom_name in sorted_atoms:
            info = atom_info_map[atom_name]
            elements = ', '.join(sorted(info['elements']))
            count = str(info['count'])
            in_bonds = "✓" if atom_name in bonded_atoms else ""

            # Highlight bonded atoms
            if atom_name in bonded_atoms:
                table.add_row(
                    f"[bold]{atom_name}[/bold]",
                    f"[bold]{elements}[/bold]",
                    count,
                    in_bonds
                )
            else:
                table.add_row(atom_name, elements, count, in_bonds)

        self.console.print(table)

        # Educational message
        info_panel = Panel(
            "[cyan]The 'In Bonds' column indicates atoms involved in redox site bonds.[/cyan]\n"
            "These are typically coordination bonds or structural bonds within the site.\n\n"
            "Enter atom names to restrain (comma-separated), or 'all' for all atoms.",
            border_style="cyan",
            expand=False
        )
        self.console.print()
        self.console.print(info_panel)

        # Get user selection
        self.console.print()
        selection = prompt_with_context(
            self.processor,
            "Select atoms to restrain",
            default="all",
            module="MD Restraint Manager",
            description="Select atoms to restrain",
        )

        # Parse selection
        if selection.lower().strip() == 'all':
            return set(atom_info_map.keys())

        # Parse comma-separated list
        selected = set()
        for atom_name in selection.split(','):
            atom_name = atom_name.strip()
            if atom_name in atom_info_map:
                selected.add(atom_name)
            else:
                self.console.print(f"[yellow]Warning: '{atom_name}' not found in redox sites - skipping[/yellow]")

        return selected

    def _configure_redox_restraints(self, redox_sites, chain_info: Dict, tleap_map: Dict, is_transformed: bool = False):
        """
        Configure redox site restraints (either entire residues or specific atoms).

        Args:
            redox_sites: List of RedoxSite objects
            chain_info: Chain composition info from _build_tleap_numbering_map
            tleap_map: Mapping from (chain, resid) to tLEaP numbering
            is_transformed: Whether structure is already in tLEaP format

        Returns:
            Dict with keys: 'method', 'site_count', 'pdb_residues', 'tleap_residues',
                           'tleap_format', 'atom_names'
            or None if user cancels
        """
        from rich.panel import Panel

        self.console.print("\n[bold cyan]═══ Redox Site Restraints ═══[/bold cyan]\n")

        # Educational message
        panel = Panel(
            f"Detected [bold]{len(redox_sites)}[/bold] redox sites across your structure.\n"
            "You can restrain either entire residues or specific atoms.",
            title="[cyan]Redox Site Information[/cyan]",
            border_style="cyan",
            expand=False
        )
        self.console.print(panel)

        # Choose method
        self.console.print("\nChoose redox restraint method:")
        self.console.print("  [cyan](r)[/cyan] Entire residues - Restrain all atoms in redox site residues")
        self.console.print("  [cyan](a)[/cyan] Specific atoms - Select specific atom types from redox sites")

        method_choice = prompt_with_context(
            processor=self.processor,
            prompt="Select method",
            choices=["r", "a"],
            default="r",
            module="MD Restraint Manager",
            description="Select redox restraint method",
            options_map={
                "r": "Entire residues",
                "a": "Specific atoms only"
            }
        )

        method = "entire_residues" if method_choice == "r" else "specific_atoms"

        # Analyze redox sites to get residues and atoms
        pdb_residues = set()
        atom_info_map = {}  # Map atom_name -> {count, in_bonds, elements}
        bonded_atoms = set()  # Track which atom names are in bonds

        for site in redox_sites:
            if hasattr(site, 'atoms'):
                for atom in site.atoms:
                    pdb_residues.add((atom.chain, atom.resid))

                    # Track atom info
                    if atom.atom_name not in atom_info_map:
                        atom_info_map[atom.atom_name] = {
                            'count': 0,
                            'elements': set()
                        }
                    atom_info_map[atom.atom_name]['count'] += 1
                    atom_info_map[atom.atom_name]['elements'].add(atom.element)

            # Track which atoms are involved in bonds
            if hasattr(site, 'bonds'):
                for bond in site.bonds:
                    # Get atom names from coordinates
                    atom1_info = site.coord_to_pdb.get(bond.atom1_coords)
                    atom2_info = site.coord_to_pdb.get(bond.atom2_coords)
                    if atom1_info:
                        bonded_atoms.add(atom1_info['atom_name'])
                    if atom2_info:
                        bonded_atoms.add(atom2_info['atom_name'])

        # If specific atoms method, show selection table
        selected_atom_names = set()
        if method == "specific_atoms":
            selected_atom_names = self._select_atoms_interactive(atom_info_map, bonded_atoms)
            if not selected_atom_names:
                self.console.print("[yellow]No atoms selected - cancelling redox restraints[/yellow]")
                return None

        # Convert to tLEaP numbering (or use directly if already transformed)
        tleap_residues = []
        pdb_residue_strs = []

        for chain, resid in sorted(pdb_residues):
            if is_transformed:
                # Already in tLEaP format - use resid directly
                tleap_residues.append(resid)
                pdb_residue_strs.append(f"{chain}:{resid}")
            else:
                # Non-transformed - map to tLEaP numbering
                tleap_resid = tleap_map.get((chain, resid))
                if tleap_resid:
                    tleap_residues.append(tleap_resid)
                    pdb_residue_strs.append(f"{chain}:{resid}")

        # Build mask format
        if method == "entire_residues":
            # Residue-based mask
            tleap_format = ':' + ','.join(map(str, sorted(tleap_residues)))
        else:
            # Atom-based mask, SCOPED to the redox-site residues. A bare
            # "@CA,CB,..." name mask matches every same-named atom in the whole
            # system (CA/CB exist in every amino acid), so without the residue
            # qualifier the restraint leaks far beyond the redox sites. AMBER
            # intersection syntax "(:residues)&(@names)" keeps only the selected
            # names that also fall within the site residues.
            atom_part = '@' + ','.join(sorted(selected_atom_names))
            if tleap_residues:
                res_part = ':' + ','.join(map(str, sorted(set(tleap_residues))))
                tleap_format = f'({res_part})&({atom_part})'
            else:
                tleap_format = atom_part

        # Display summary with clear explanations
        self.console.print(f"\n[cyan]Restraint Configuration Summary:[/cyan]")

        if method == "entire_residues":
            self.console.print(f"  Residues to restrain: {', '.join(pdb_residue_strs[:10])}" +
                             (f" ... ({len(pdb_residue_strs)} total)" if len(pdb_residue_strs) > 10 else ""),
                             highlight=False)
            self.console.print(f"    [grey50](Chain:Residue# notation from structure)[/grey50]")
        else:
            # Atom-based restraints
            self.console.print(f"  Atoms to restrain: {', '.join(sorted(selected_atom_names))}", highlight=False)
            self.console.print(f"  Found in residues: {', '.join(pdb_residue_strs[:10])}" +
                             (f" ... ({len(pdb_residue_strs)} total)" if len(pdb_residue_strs) > 10 else ""),
                             highlight=False)
            self.console.print(f"    [grey50](Chain:Residue# notation from structure)[/grey50]")

        # Show the full mask (no truncation) so the user sees exactly what
        # goes into the MD input; highlight=False so AMBER syntax (numbers,
        # colons, '@') isn't recoloured by Rich's auto-highlighter.
        self.console.print(f"  AMBER restraint mask: {tleap_format}", highlight=False)
        self.console.print(f"    [grey50](Syntax used in MD input files)[/grey50]")

        return {
            'method': method,
            'site_count': len(redox_sites),
            'pdb_residues': pdb_residue_strs,
            'tleap_residues': tleap_residues,
            'tleap_format': tleap_format,
            'atom_names': list(selected_atom_names) if method == "specific_atoms" else [],
        }

    def _preview_and_confirm_mask(self, general_config: Optional[Dict],
                                  redox_config: Optional[Dict],
                                  chain_info: Dict) -> Tuple[bool, str]:
        """
        Preview combined restraint mask and get user confirmation.

        Returns:
            (accepted: bool, final_mask: str)
        """
        from rich.panel import Panel
        from rich.table import Table

        self.console.print("\n[bold cyan]═══ Restraint Mask Preview ═══[/bold cyan]\n")

        # Build preview content
        preview_lines = []
        mask_components = []

        if general_config:
            if general_config.get('is_custom_mask', False):
                # Custom mask - simplified display
                preview_lines.append("[bold]Custom Restraint Mask:[/bold]")
                preview_lines.append(f"  {general_config['atom_description']}")
                preview_lines.append(f"  Mask: {general_config['tleap_format']}\n")
            else:
                preview_lines.append("[bold]General Restraints:[/bold]")
                preview_lines.append(f"  Atom selection: {general_config['atom_description']}")
                preview_lines.append(f"  Chains: {', '.join(general_config['chains'])}")
                preview_lines.append(f"  tLEaP format: {general_config['tleap_format']}\n")
            mask_components.append(general_config['tleap_format'])

        if redox_config:
            preview_lines.append("[bold]Redox Site Restraints:[/bold]")
            preview_lines.append(f"  Method: {redox_config['method'].replace('_', ' ').title()}")
            preview_lines.append(f"  Redox sites: {redox_config['site_count']} sites")
            preview_lines.append(f"  Residues: {len(redox_config['pdb_residues'])} residues")
            preview_lines.append(f"  tLEaP format: {redox_config['tleap_format'][:60]}" +
                               ("..." if len(redox_config['tleap_format']) > 60 else "") + "\n")
            mask_components.append(redox_config['tleap_format'])

        # Combine masks
        if len(mask_components) > 1:
            combined_mask = '|'.join(mask_components)
        elif len(mask_components) == 1:
            combined_mask = mask_components[0]
        else:
            combined_mask = ""

        # Check if this is a custom mask (user provided complete mask syntax)
        is_custom_mask = general_config and general_config.get('is_custom_mask', False)

        # Ask user about water exclusion (skip for custom masks - user already specified everything)
        if combined_mask:
            if is_custom_mask:
                # Custom mask - use as-is without modification
                final_mask = combined_mask
                preview_lines.append("[grey50]Custom mask used directly without modification[/grey50]\n")
            else:
                exclude_water = confirm_with_context(
                    processor=self.processor,
                    prompt="[yellow]Exclude water molecules from restraints?[/yellow]",
                    default=True,
                    module="MD Restraint Manager",
                    description="Decide whether to exclude water from restraint mask"
                )

                if exclude_water:
                    # Apply water exclusion to entire mask (use parentheses for multiple components)
                    if len(mask_components) > 1:
                        final_mask = f"({combined_mask})&!:WAT"
                    else:
                        final_mask = f"{combined_mask}&!:WAT"
                    preview_lines.append("[grey50]Water molecules will be excluded from restraints[/grey50]\n")
                else:
                    final_mask = combined_mask
                    preview_lines.append("[grey50]Water molecules will NOT be excluded from restraints[/grey50]\n")
        else:
            final_mask = ""

        preview_lines.append("[bold yellow]FINAL COMBINED MASK (saved to workspace):[/bold yellow]")
        preview_lines.append(f"  [bold]{final_mask}[/bold]\n")

        # Estimate atom count (simplified - just show that we're restraining atoms)
        preview_lines.append("[grey50]This mask is compatible with tLEaP/AMBER MD simulations.[/grey50]")

        panel = Panel(
            "\n".join(preview_lines),
            title="[bold cyan]Restraint Mask Preview[/bold cyan]",
            border_style="cyan",
            expand=False
        )

        self.console.print(panel)

        # Hook 7: highlight the union of residues covered by the
        # final mask so the user can verify the spatial extent
        # before accepting. We use the chain:resid spec lists from
        # general_config and redox_config (parsed pre-tLEaP) rather
        # than re-translating the tLEaP mask back. Two distinct
        # colours so the user can tell which restraint type covers
        # which residue. Cleanup happens at every return path.
        try:
            viewer = self._viewer_or_none()
            if viewer is not None:
                # Drop earlier per-chain general_mask labels so they
                # don't conflict with the unified preview rep below.
                for lbl in getattr(self, "_general_mask_labels", None) or []:
                    viewer.unhighlight(lbl)
                self._general_mask_labels = []

                viewer.unhighlight(f"{self._VIEWER_LABEL_PREFIX}preview_general")
                viewer.unhighlight(f"{self._VIEWER_LABEL_PREFIX}preview_redox")

                if general_config and not general_config.get('is_custom_mask', False):
                    pairs: set = set()
                    for spec in general_config.get('pdb_ranges') or []:
                        if ':' not in spec:
                            continue
                        chain, rng = spec.split(':', 1)
                        chain = chain.strip()
                        rng = rng.strip()
                        if '-' in rng:
                            try:
                                start, end = (int(x) for x in rng.split('-', 1))
                            except ValueError:
                                continue
                            for r in range(start, end + 1):
                                pairs.add((chain, r))
                        else:
                            try:
                                pairs.add((chain, int(rng)))
                            except ValueError:
                                continue
                    if pairs:
                        clauses = [f"(:{c} and {r})" for c, r in sorted(pairs)]
                        viewer.highlight(
                            " or ".join(clauses),
                            style="ball+stick",
                            color="#1f78b4",
                            label=f"{self._VIEWER_LABEL_PREFIX}preview_general",
                        )

                if redox_config:
                    pairs = set()
                    for spec in redox_config.get('pdb_residues') or []:
                        if ':' not in spec:
                            continue
                        chain, resid = spec.split(':', 1)
                        try:
                            pairs.add((chain.strip(), int(resid.strip())))
                        except ValueError:
                            continue
                    if pairs:
                        clauses = [f"(:{c} and {r})" for c, r in sorted(pairs)]
                        viewer.highlight(
                            " or ".join(clauses),
                            style="ball+stick",
                            color="#e31a1c",
                            label=f"{self._VIEWER_LABEL_PREFIX}preview_redox",
                        )
        except Exception as exc:
            logger.debug("mask preview hook silenced: %s", exc)

        # Get confirmation
        self.console.print("\nOptions:")
        self.console.print("  [cyan](a)[/cyan] Accept this mask")
        self.console.print("  [cyan](e)[/cyan] Edit - start over")
        self.console.print("  [cyan](c)[/cyan] Cancel")

        choice = prompt_with_context(
            processor=self.processor,
            prompt="Choose action",
            choices=["a", "e", "c"],
            default="a",
            module="MD Restraint Manager",
            description="Accept, edit, or cancel restraint mask",
            options_map={
                "a": "Accept mask",
                "e": "Edit - start over",
                "c": "Cancel"
            }
        )

        # Drop the preview reps regardless of the user's choice — the
        # mask is either accepted (preview no longer needed), or the
        # workflow restarts/cancels (next step shouldn't inherit
        # stale highlights).
        viewer = self._viewer_or_none()
        if viewer is not None:
            viewer.unhighlight(f"{self._VIEWER_LABEL_PREFIX}preview_general")
            viewer.unhighlight(f"{self._VIEWER_LABEL_PREFIX}preview_redox")

        if choice == "a":
            return True, final_mask
        elif choice == "e":
            return False, ""  # Signal to restart
        else:
            return None, ""  # Signal to cancel

    def generate_restraint_mask(self, workspace=None, interactive=True) -> Optional[Dict[str, Any]]:
        """
        Generate AMBER restraint mask with support for general and redox-specific restraints.

        This new implementation supports:
        - General structural restraints (by chains, atom types)
        - Redox site-specific restraints (by residues or atoms)
        - Combination of both
        - tLEaP consecutive numbering for multi-chain structures

        Args:
            workspace: Workspace containing structure and optional redox sites
            interactive: Whether to prompt user for options

        Returns:
            Dict containing restraint mask and metadata, or None if failed
        """
        from datetime import datetime
        from rich.panel import Panel

        try:
            # Get structure (required)
            structure_info = self._get_priority_structure(workspace) if workspace else None
            if not structure_info:
                self.console.print("[yellow]No structure available for restraint mask generation.[/yellow]")
                return None

            # Unpack structure info (should be tuple of 2 values)
            try:
                structure_file, structure_source = structure_info
            except (TypeError, ValueError) as e:
                logger.error(f"Invalid structure_info format: {structure_info}, error: {e}")
                self.console.print(f"[red]Error: Invalid structure information format[/red]")
                return None

            # Detect if structure already has AMBER-consecutive numbering
            is_transformed = (structure_source in ("transformed", "topology-extracted"))

            # Build tLEaP numbering map
            chain_info, tleap_map = self._build_tleap_numbering_map(structure_file)
            if not chain_info:
                self.console.print("[red]Error: Could not parse structure for tLEaP numbering[/red]")
                return None

            # Display tLEaP numbering information
            if interactive:
                self._display_tleap_numbering_info(chain_info, is_transformed, structure_file, structure_source)

            # Get redox sites (optional)
            redox_sites = self.get_from_workspace("detected_redox_sites", []) if self.processor else []
            if workspace:
                redox_sites = self._get_from_workspace_obj(workspace, "detected_redox_sites", [])

            # Main workflow loop (allow user to restart if they choose "edit")
            while True:
                if not interactive:
                    # Non-interactive mode: default to general restraints only
                    general_config = None
                    redox_config = None
                    break

                # Step 1: Choose restraint strategy
                self.console.print("\n[bold cyan]═══ Restraint Mask Generation ═══[/bold cyan]\n")

                if redox_sites:
                    # Redox sites available - show all options
                    self.console.print("Choose base restraint strategy:")
                    self.console.print("  [cyan](a)[/cyan] General structural restraints only")
                    self.console.print("  [cyan](b)[/cyan] Redox site restraints only")
                    self.console.print("  [cyan](c)[/cyan] Combination: General + Redox site restraints")
                    strategy_choices = ["a", "b", "c"]
                else:
                    # No redox sites - only show general option
                    self.console.print("[yellow]Note: No redox sites detected in workspace[/yellow]")
                    self.console.print("\nOnly general structural restraints are available:")
                    self.console.print("  [cyan](a)[/cyan] General structural restraints")
                    strategy_choices = ["a"]

                strategy = prompt_with_context(
                    processor=self.processor,
                    prompt="Select strategy",
                    choices=strategy_choices,
                    default="a",
                    module="MD Restraint Manager",
                    description="Select restraint strategy",
                    options_map={
                        "a": "General restraints only",
                        "b": "Redox site restraints only",
                        "c": "Combination (general + redox)"
                    }
                )

                # Step 2: Configure restraints based on strategy
                general_config = None
                redox_config = None

                if strategy in ["a", "c"]:
                    # Configure general restraints
                    general_config = self._configure_general_restraints(chain_info, tleap_map)

                if strategy in ["b", "c"]:
                    # Configure redox restraints
                    if redox_sites:
                        redox_config = self._configure_redox_restraints(redox_sites, chain_info, tleap_map, is_transformed)
                    else:
                        self.console.print("[yellow]No redox sites available for restraint configuration[/yellow]")

                # Step 3: Preview and confirm
                accepted, final_mask = self._preview_and_confirm_mask(general_config, redox_config, chain_info)

                if accepted is True:
                    # User accepted - break out of loop
                    break
                elif accepted is False:
                    # User wants to edit - restart loop
                    self.console.print("\n[yellow]Restarting restraint configuration...[/yellow]")
                    continue
                else:
                    # User cancelled
                    self.console.print("[yellow]Restraint mask generation cancelled[/yellow]")
                    return None

            # Build comprehensive metadata for workspace
            metadata = {
                "mask": final_mask,
                "mask_format": "tleap_consecutive",
                "timestamp": datetime.now().isoformat(),
                "structure_source": structure_file,
            }

            # Add general restraints info
            if general_config:
                metadata["general_restraints"] = {
                    "enabled": True,
                    "atom_selection": general_config['atom_selection'],
                    "atom_description": general_config['atom_description'],
                    "chains_selected": general_config['chains'],
                    "pdb_ranges": general_config['pdb_ranges'],
                    "tleap_range": general_config['tleap_ranges'],
                    "tleap_mask_component": general_config['tleap_format']
                }
            else:
                metadata["general_restraints"] = None

            # Add redox restraints info
            if redox_config:
                metadata["redox_restraints"] = {
                    "enabled": True,
                    "method": redox_config['method'],
                    "site_count": redox_config['site_count'],
                    "pdb_residues": redox_config['pdb_residues'],
                    "tleap_residues": redox_config['tleap_residues'],
                    "tleap_mask_component": redox_config['tleap_format'],
                    "atom_names": redox_config.get('atom_names', [])
                }
            else:
                metadata["redox_restraints"] = None

            # Add tLEaP numbering map for reference
            metadata["tleap_numbering"] = {}
            for chain, info in chain_info.items():
                metadata["tleap_numbering"][chain] = {
                    "start": info['start'],
                    "end": info['end'],
                    "count": info['count'],
                    "ranges": info['ranges']  # List of (start, end) tuples for PDB residue ranges
                }

            # Save to workspace - CRITICAL!
            if self.processor:
                self.update_workspace("redox_restraint_mask", final_mask)
                self.update_workspace("redox_restraint_info", metadata)
                self.update_workspace("restraint_structure_source", structure_file)
                self.update_workspace("restraint_mask_generated", True)
            elif workspace:
                self._update_workspace_obj(workspace, "redox_restraint_mask", final_mask)
                self._update_workspace_obj(workspace, "redox_restraint_info", metadata)
                self._update_workspace_obj(workspace, "restraint_structure_source", structure_file)
                self._update_workspace_obj(workspace, "restraint_mask_generated", True)

            # Display success message with educational context
            if interactive:
                success_panel = Panel(
                    "[bold green]✓ Restraint mask saved to workspace[/bold green]\n\n"
                    "The mask is saved as [cyan]'redox_restraint_mask'[/cyan] and will be\n"
                    "automatically used by:\n"
                    "  • MD Manager (minimization/equilibration)\n\n"
                    "You can also manually use it in AMBER input files:\n"
                    f"  [grey50]restraintmask='{final_mask}'[/grey50]",
                    title="[bold green]Success[/bold green]",
                    border_style="green",
                    expand=False
                )
                self.console.print()
                self.console.print(success_panel)

            result = {
                "restraint_mask": final_mask,
                "info": metadata
            }

            return result
            
        except Exception as e:
            import traceback
            logger.error(f"Error generating restraint mask: {str(e)}")
            logger.error(f"Traceback: {traceback.format_exc()}")
            self.console.print(f"[red]Error generating restraint mask: {str(e)}[/red]")
            self.console.print(f"[grey50]{traceback.format_exc()}[/grey50]")
            return None

    def _analyze_redox_sites(self, redox_sites):
        """Analyze RedoxSite objects to identify atoms for restraints."""
        for site in redox_sites:
            try:
                # Add redox centers (the primary redox-active atoms)
                for center in site.centers:
                    self.redox_atoms.add((center.chain, str(center.resid), center.atom_name))
                
                # Add other site atoms (atoms that are part of the redox site)
                center_keys = {(c.chain, str(c.resid), c.atom_name) for c in site.centers}
                for atom in site.atoms:
                    atom_key = (atom.chain, str(atom.resid), atom.atom_name)
                    if atom_key not in center_keys:
                        self.redox_atoms.add(atom_key)
                        
                logger.debug(f"Analyzed redox site {site.site_id}: found {len(self.redox_atoms)} atoms")
                
            except Exception as e:
                logger.warning(f"Error analyzing redox site: {e}")
                continue

    def _analyze_redox_sites_full_residues(self, redox_sites):
        """Analyze RedoxSite objects to identify all atoms in redox site residues."""
        for site in redox_sites:
            try:
                # Add all residues from all atoms in the site
                for atom in site.atoms:
                    self.redox_residues.add((atom.chain, str(atom.resid)))

                logger.debug(f"Analyzed redox site {site.site_id}: found {len(self.redox_residues)} residues")

            except Exception as e:
                logger.warning(f"Error analyzing redox site: {e}")
                continue

        # Collect all atoms from these residues
        self._collect_atoms_from_residues()

    def _collect_atoms_from_residues(self):
        """Collect all atoms from the identified redox residues."""
        for line in self.structure_lines:
            if not line.startswith(('ATOM  ', 'HETATM')):
                continue
                
            try:
                # Parse PDB line format
                chain = line[21:22].strip()
                resid = line[22:26].strip()
                atom_name = line[12:16].strip()
                
                # Check if this atom is in our redox residues
                res_key = (chain, resid)
                if res_key in self.redox_residues:
                    self.redox_atoms.add((chain, resid, atom_name))
                    logger.debug(f"Added atom: {atom_name} from residue {chain}:{resid}")
                
            except Exception as e:
                logger.debug(f"Error parsing structure line: {e}")
                continue

    def _parse_structure_for_redox_atoms(self) -> Set[str]:
        """Parse structure to find actual atom names for redox atoms."""
        found_atoms = set()
        
        for line in self.structure_lines:
            if not line.startswith(('ATOM  ', 'HETATM')):
                continue
                
            try:
                # Parse PDB line format
                chain = line[21:22].strip()
                resid = line[22:26].strip()
                atom_name = line[12:16].strip()
                
                # Check if this atom is in our redox atoms set
                atom_key = (chain, resid, atom_name)
                if atom_key in self.redox_atoms:
                    found_atoms.add(atom_name)
                    logger.debug(f"Found redox atom: {atom_name} in {chain}:{resid}")
                
            except Exception as e:
                logger.debug(f"Error parsing structure line: {e}")
                continue
        
        return found_atoms

    def _create_restraint_mask(self, redox_atom_names: Set[str]) -> str:
        """Create AMBER restraint mask including backbone and redox atoms."""
        # Start with backbone atoms (excluding water)
        backbone_mask = "@CA,C,O,N&!:WAT"
        
        if not redox_atom_names:
            return backbone_mask
        
        # Sort atom names for consistent output
        sorted_atoms = sorted(redox_atom_names)
        redox_mask = "@" + ",".join(sorted_atoms)
        
        # Combine backbone and redox atom masks
        combined_mask = f"{backbone_mask}|{redox_mask}"
        
        return combined_mask

    def _create_residue_based_mask(self, interactive: bool = True, redox_sites=None) -> str:
        """Create AMBER restraint mask based on entire residues."""
        from rich.table import Table
        from proprep.utils.prompts import confirm_with_context
        
        # Get residue information with names, grouped by site
        residue_info = self._get_residue_info_with_names()
        
        if interactive and residue_info:
            # Group residues by RedoxSite
            residues_by_site = {}
            if redox_sites:
                for site in redox_sites:
                    site_residues = []
                    site_res_keys = set()
                    # Get all residues from this site
                    for atom in site.atoms:
                        res_key = (atom.chain, str(atom.resid))
                        site_res_keys.add(res_key)
                    
                    # Find matching residues in our residue_info
                    for res_info in residue_info:
                        res_key = (res_info['chain'], res_info['resid'])
                        if res_key in site_res_keys:
                            site_residues.append(res_info)
                    
                    if site_residues:
                        residues_by_site[site.site_id] = site_residues
            
            # Display residues grouped by site
            total_residues = sum(len(site_residues) for site_residues in residues_by_site.values())
            self.console.print(f"\n[bold]Found {total_residues} redox residues for restraints:[/bold]")
            
            table = Table(title="Redox Residues for Restraints (Grouped by Site)")
            table.add_column("#", style="cyan", width=4)
            table.add_column("Res Name", style="magenta", width=8)
            table.add_column("Chain:ResID", style="green")
            table.add_column("Include", style="yellow", width=8)
            
            row_num = 1
            first_site = True
            
            for site_id, site_residues in residues_by_site.items():
                # Add horizontal divider between sites (except before first site)
                if not first_site:
                    table.add_row("", "─" * 8, "─" * 11, "─" * 8, style="grey50")
                
                # Add site header row
                site_header = f"[bold cyan]Site: {site_id}[/bold cyan]"
                table.add_row("", site_header, "", "", style="bold cyan")
                
                # Add residues for this site
                for res_info in site_residues:
                    table.add_row(
                        str(row_num),
                        res_info['resname'],
                        f"{res_info['chain']}:{res_info['resid']}",
                        "✓"
                    )
                    row_num += 1
                
                first_site = False
            
            self.console.print(table)
            
            # Allow removal of entire residues
            if confirm_with_context(
                processor=self.processor,
                prompt="\nWould you like to edit this list?",
                default=False,
                module="MD Restraint Manager",
                description="Edit residue restraint list"
            ):
                residue_info = self._interactive_edit_residues(residue_info)
        
        # Create mask from selected residues
        if not residue_info:
            return "@CA,C,O,N&!:WAT"
        
        # Build residue-based mask
        residue_specs = []
        for res_info in residue_info:
            residue_specs.append(f":{res_info['resid']}")
        
        # Combine with backbone
        backbone_mask = "@CA,C,O,N&!:WAT"
        residue_mask = "|".join(residue_specs)
        
        combined_mask = f"{backbone_mask}|{residue_mask}"
        return combined_mask

    def _get_residue_info_with_names(self) -> List[Dict[str, str]]:
        """Get residue information including residue names from PDB structure."""
        residue_info = []
        seen_residues = set()
        
        # Parse structure to get residue names
        for line in self.structure_lines:
            if not line.startswith(('ATOM  ', 'HETATM')):
                continue
                
            try:
                # Parse PDB line format  
                resname = line[17:20].strip()
                chain = line[21:22].strip()
                resid = line[22:26].strip()
                
                # Check if this is one of our redox residues
                res_key = (chain, resid)
                if res_key in self.redox_residues and res_key not in seen_residues:
                    residue_info.append({
                        'resname': resname,
                        'chain': chain,
                        'resid': resid
                    })
                    seen_residues.add(res_key)
                    
            except Exception as e:
                logger.debug(f"Error parsing structure line for residue names: {e}")
                continue
        
        # Sort by resname then chain then resid
        return sorted(residue_info, key=lambda x: (x['resname'], x['chain'], int(x['resid']) if x['resid'].isdigit() else x['resid']))

    def _interactive_edit_atoms(self, atom_names: Set[str], redox_sites=None) -> Set[str]:
        """Allow user to interactively edit the list of atoms to include."""
        from rich.table import Table

        if not atom_names:
            return atom_names
        
        # Group atoms by residue and get residue names  
        atoms_by_residue = self._group_atoms_by_residue_with_names(atom_names)
        total_atom_count = sum(len(atoms) for atoms in atoms_by_residue.values())
        
        self.console.print(f"\n[bold]Found {total_atom_count} redox atoms in {len(atoms_by_residue)} residues for restraints:[/bold]")
        
        # Group by RedoxSite if available
        if redox_sites:
            atoms_by_site = self._group_atoms_by_site(atoms_by_residue, redox_sites)
        else:
            atoms_by_site = {"Unknown Site": atoms_by_residue}
        
        # Create a flat list for indexing while preserving visual grouping
        display_atoms = []  # List of (index, atom_name, residue_info)
        index = 1
        
        # Display atoms grouped by site then residue in a table
        table = Table(title="Redox Atoms for Restraints (Grouped by Site)")
        table.add_column("#", style="cyan", width=4)
        table.add_column("Residue", style="magenta", width=15)
        table.add_column("Atom Name", style="green")
        table.add_column("Include", style="yellow", width=8)
        
        first_site = True
        for site_id, site_residues in atoms_by_site.items():
            # Add horizontal divider between sites (except before first site)
            if not first_site:
                table.add_row("", "─" * 15, "─" * 11, "─" * 8, style="grey50")
            
            # Add site header row
            site_header = f"[bold cyan]Site: {site_id}[/bold cyan]"
            table.add_row("", site_header, "", "", style="bold cyan")
            
            # Add residues and atoms for this site
            for (chain, resid, resname), residue_atoms in site_residues.items():
                residue_label = f"{resname} {chain}:{resid}"
                
                # Sort atoms within each residue alphabetically
                sorted_residue_atoms = sorted(residue_atoms)
                
                for i, atom in enumerate(sorted_residue_atoms):
                    # Only show residue label on first atom of each residue
                    display_residue = residue_label if i == 0 else ""
                    table.add_row(str(index), display_residue, atom, "✓")
                    display_atoms.append((index, atom, residue_label))
                    index += 1
            
            first_site = False
        
        self.console.print(table)

        if not confirm_with_context(
            processor=self.processor,
            prompt="\nWould you like to edit this list?",
            default=False,
            module="MD Restraint Manager",
            description="Edit atom restraint list"
        ):
            return atom_names
        
        # Allow user to remove atoms
        edited_atoms = set(atom_names)
        
        self.console.print("\n[bold]Edit restraint atoms:[/bold]")
        self.console.print("Enter atom numbers to REMOVE (e.g., '1,3,5' or '1-3,5')")
        self.console.print("Press Enter to keep all atoms")
        
        remove_input = prompt_with_context(
            self.processor, "Atoms to remove",
            module="MD Restraint Manager",
            description="Atoms to remove from restraint",
        ).strip()
        
        if remove_input:
            try:
                # Parse removal indices
                indices_to_remove = self._parse_atom_indices(remove_input, len(display_atoms))
                
                for idx in sorted(indices_to_remove, reverse=True):
                    if 1 <= idx <= len(display_atoms):
                        _, atom_to_remove, residue_info = display_atoms[idx - 1]
                        edited_atoms.discard(atom_to_remove)
                        self.console.print(f"[yellow]Removed: {atom_to_remove} (from {residue_info})[/yellow]")
                
                self.console.print(f"\n[green]Final restraint atoms ({len(edited_atoms)}):[/green] {', '.join(sorted(edited_atoms))}")
                
            except Exception as e:
                self.console.print(f"[red]Error parsing input: {e}[/red]")
                self.console.print("[yellow]Keeping all atoms[/yellow]")
                edited_atoms = atom_names
        
        return edited_atoms

    def _group_atoms_by_residue(self, atom_names: Set[str]) -> Dict[Tuple[str, str], List[str]]:
        """Group atoms by their residue (chain, resid)."""
        atoms_by_residue = {}
        
        # Map atom names back to their residues using redox_atoms
        for chain, resid, atom_name in self.redox_atoms:
            if atom_name in atom_names:
                res_key = (chain, resid)
                if res_key not in atoms_by_residue:
                    atoms_by_residue[res_key] = []
                atoms_by_residue[res_key].append(atom_name)
        
        # Sort residues by chain then resid for consistent display
        return dict(sorted(atoms_by_residue.items(), key=lambda x: (x[0][0], int(x[0][1]) if x[0][1].isdigit() else x[0][1])))

    def _group_atoms_by_residue_with_names(self, atom_names: Set[str]) -> Dict[Tuple[str, str, str], List[str]]:
        """Group atoms by their residue (chain, resid, resname)."""
        atoms_by_residue = {}
        
        # Parse structure to get residue names
        residue_names = {}
        for line in self.structure_lines:
            if not line.startswith(('ATOM  ', 'HETATM')):
                continue
            try:
                resname = line[17:20].strip()
                chain = line[21:22].strip()
                resid = line[22:26].strip()
                residue_names[(chain, resid)] = resname
            except:
                continue
        
        # Map atom names back to their residues using redox_atoms
        for chain, resid, atom_name in self.redox_atoms:
            if atom_name in atom_names:
                resname = residue_names.get((chain, resid), "UNK")
                res_key = (chain, resid, resname)
                if res_key not in atoms_by_residue:
                    atoms_by_residue[res_key] = []
                atoms_by_residue[res_key].append(atom_name)
        
        # Sort residues by resname, then chain, then resid for consistent display
        return dict(sorted(atoms_by_residue.items(), key=lambda x: (x[0][2], x[0][0], int(x[0][1]) if x[0][1].isdigit() else x[0][1])))

    def _group_atoms_by_site(self, atoms_by_residue: Dict[Tuple[str, str, str], List[str]], redox_sites) -> Dict[str, Dict[Tuple[str, str, str], List[str]]]:
        """Group residues by RedoxSite."""
        atoms_by_site = {}
        
        for site in redox_sites:
            site_residues = {}
            site_res_keys = set()
            
            # Get all residues from this site
            for atom in site.atoms:
                res_key = (atom.chain, str(atom.resid))
                site_res_keys.add(res_key)
            
            # Find matching residues in our atoms_by_residue
            for (chain, resid, resname), residue_atoms in atoms_by_residue.items():
                res_key = (chain, resid)
                if res_key in site_res_keys:
                    site_residues[(chain, resid, resname)] = residue_atoms
            
            if site_residues:
                atoms_by_site[site.site_id] = site_residues
        
        return atoms_by_site

    def _parse_atom_indices(self, input_str: str, max_index: int) -> Set[int]:
        """Parse atom indices from user input (supports ranges like '1-3,5,7-9')."""
        indices = set()
        
        for part in input_str.split(','):
            part = part.strip()
            if '-' in part:
                # Range like '1-3'
                try:
                    start, end = map(int, part.split('-', 1))
                    indices.update(range(start, end + 1))
                except ValueError:
                    continue
            else:
                # Single number
                try:
                    indices.add(int(part))
                except ValueError:
                    continue
        
        # Filter valid indices
        return {idx for idx in indices if 1 <= idx <= max_index}

    def _interactive_edit_residues(self, residue_info: List[Dict[str, str]]) -> List[Dict[str, str]]:
        """Allow user to remove entire residues from restraints."""
        self.console.print("\n[bold]Edit residue list:[/bold]")
        self.console.print("Enter residue numbers to REMOVE (e.g., '1,3,5' or '1-3,5')")
        self.console.print("Press Enter to keep all residues")
        
        remove_input = prompt_with_context(
            self.processor, "Residues to remove",
            module="MD Restraint Manager",
            description="Residues to remove from restraint",
        ).strip()
        
        if remove_input:
            try:
                # Parse removal indices
                indices_to_remove = self._parse_atom_indices(remove_input, len(residue_info))
                
                # Remove residues (in reverse order to maintain indices)
                for idx in sorted(indices_to_remove, reverse=True):
                    if 1 <= idx <= len(residue_info):
                        removed = residue_info.pop(idx - 1)
                        self.console.print(f"[yellow]Removed: {removed['resname']} {removed['chain']}:{removed['resid']}[/yellow]")
                
                self.console.print(f"\n[green]Final residues ({len(residue_info)}):[/green]")
                for res in residue_info:
                    self.console.print(f"  {res['resname']} {res['chain']}:{res['resid']}")
                    
            except Exception as e:
                self.console.print(f"[red]Error parsing input: {e}[/red]")
                self.console.print("[yellow]Keeping all residues[/yellow]")
        
        return residue_info

    def _get_redox_residues_info(self) -> List[Dict[str, Any]]:
        """Get information about redox residues."""
        redox_residues = []
        
        # Group redox atoms by residue
        residue_atoms = {}
        for chain, resid, atom_name in self.redox_atoms:
            res_key = (chain, resid)
            if res_key not in residue_atoms:
                residue_atoms[res_key] = []
            residue_atoms[res_key].append(atom_name)
        
        # Create residue info
        for (chain, resid), atoms in residue_atoms.items():
            redox_residues.append({
                "chain": chain,
                "resid": resid,
                "atom_count": len(atoms),
                "atoms": sorted(atoms)
            })
        
        return redox_residues

    def _count_atoms_in_selected_residues(self):
        """Count atoms in the selected redox residues for residue-based restraints."""
        atom_names = set()
        
        # Get the redox residues that were selected
        redox_residues = self._get_redox_residues_info()
        
        try:
            # Count atoms in each selected redox residue
            for residue_info in redox_residues:
                chain_id = residue_info["chain"]
                resid = residue_info["resid"]
                
                # Add atom identifiers from the structure
                for line in self.structure_lines:
                    if not line.startswith(('ATOM  ', 'HETATM')):
                        continue
                        
                    try:
                        # Parse PDB line format
                        chain = line[21:22].strip()
                        res_id = line[22:26].strip()
                        atom_name = line[12:16].strip()
                        
                        if chain == chain_id and res_id == resid:
                            atom_names.add(f"{chain_id}:{resid}:{atom_name}")
                            
                    except Exception as e:
                        continue
                        
        except Exception as e:
            self.console.print(f"[yellow]Warning: Could not count atoms in residues: {e}[/yellow]")
            
        return atom_names

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Workspace management methods (following structure_completeness pattern)
    
    def _update_workspace_obj(self, workspace, key, value):
        """Update workspace object (for module-level processing)"""
        if hasattr(workspace, 'set'):
            old_value = workspace.get(key)
            workspace.set(key, value)
            
            debug_enabled = workspace.get("debug", False)
            if debug_enabled:
                self._debug_value(key, value, old_value, "updated")
        elif isinstance(workspace, dict):
            workspace[key] = value
    
    def _get_from_workspace_obj(self, workspace, key, default=None):
        """Get from workspace object (for module-level processing)"""
        if hasattr(workspace, 'get'):
            return workspace.get(key, default)
        elif isinstance(workspace, dict):
            return workspace.get(key, default)
        return default
    
    def update_workspace(self, key, value):
        """
        Helper method to update the processor's workspace
        
        Args:
            key: Key to update
            value: New value
        """
        if self.processor and hasattr(self.processor, 'workspace'):
            old_value = self.processor.workspace.get(key)
            self.processor.workspace.set(key, value)
            
            # Debug output if enabled
            debug_enabled = self.processor.workspace.get("debug", False)
            if debug_enabled:
                self._debug_value(key, value, old_value, "updated")
    
    def get_from_workspace(self, key, default=None):
        """
        Helper method to get values from the processor's workspace
        
        Args:
            key: Key to retrieve
            default: Default value if key not found
            
        Returns:
            The value for the given key or default
        """
        if self.processor and hasattr(self.processor, 'workspace'):
            return self.processor.workspace.get(key, default)
        return default
    
    def _debug_value(self, key, value, old_value, action):
        """
        Debug output for value updates (following structure_completeness pattern)
        
        Args:
            key: Key being updated
            value: New value
            old_value: Previous value
            action: Action being performed
        """
        if key == "debug":
            return  # Don't debug the debug flag itself

        value_type = type(value).__name__
        if hasattr(value, "__len__") and not isinstance(value, (str, bytes)):
            value_info = f"{value_type}({len(value)})"
        else:
            value_info = value_type

        if old_value is None:
            self.console.print(
                f"[yellow]DEBUG: RestraintManager {action} new value for '{key}': {value_info}[/yellow]"
            )
        else:
            old_type = type(old_value).__name__
            if hasattr(old_value, "__len__") and not isinstance(old_value, (str, bytes)):
                old_info = f"{old_type}({len(old_value)})"
            else:
                old_info = old_type
            self.console.print(
                f"[yellow]DEBUG: RestraintManager {action} '{key}': {old_info} -> {value_info}[/yellow]"
            )

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Module interface methods
    
    def get_workspace_requirements(self) -> List[str]:
        """Get workspace requirements for this module"""
        return []  # No hard requirements - will find best available structure
    
    def can_process(self, workspace) -> bool:
        """Check if workspace has required data for processing"""
        return self._get_priority_structure(workspace) is not None

    def process(self, workspace):
        """Process workspace - this is a RestraintManager method, not for module interface"""
        # Note: This method exists for backward compatibility but shouldn't be called directly
        # The MDRestraintModule.process() should be used as the module entry point
        return self.add_restraints_to_structure(workspace)

    # =========================================================================
    # VIEWER HOOKS — best-effort visualisation of restraint picks
    # =========================================================================
    # All hooks log silently on failure; the configurator must complete even
    # if the coordinator can't reach the viewer. Stable per-hook label
    # prefixes let later hooks clear what earlier hooks drew without
    # disturbing other modules' overlays.

    _VIEWER_LABEL_PREFIX = "mdrestraint_"

    def _viewer_or_none(self):
        """Return the coordinator viewer or None on import failure."""
        try:
            from proprep.structure_prep.viewer_coordinator import (
                viewer as _viewer,
            )
            return _viewer
        except Exception as exc:
            logger.debug("md restraint viewer import silenced: %s", exc)
            return None

    def _clear_pick_labels(self) -> None:
        """Drop labels left behind by Hooks 2/3/4 between iterations.
        Tracked on ``self`` rather than per-call so the cleanup is
        possible from any later hook without threading state through.
        """
        viewer = self._viewer_or_none()
        if viewer is None:
            return
        for label in getattr(self, "_pick_labels", None) or []:
            viewer.unhighlight(label)
        self._pick_labels = []

    def _highlight_picked_residues(self, residue_specs) -> None:
        """Hook 2: highlight the 2/3/4 residues the user typed in for
        a DISANG restraint, palette-coloured by their position in the
        spec (residue 1 = palette:1, etc.). Lets the user verify their
        spec before they choose specific atoms in the next prompt.
        """
        viewer = self._viewer_or_none()
        if viewer is None:
            return
        self._clear_pick_labels()
        applied = []
        for idx, (chain, resid) in enumerate(residue_specs, 1):
            label = f"{self._VIEWER_LABEL_PREFIX}res_{idx}"
            try:
                viewer.highlight(
                    f":{chain} and {resid}",
                    style="ball+stick",
                    color=f"palette:{idx}",
                    label=label,
                )
                applied.append(label)
            except Exception as exc:
                logger.debug("residue highlight silenced: %s", exc)
        self._pick_labels = applied

    def _highlight_picked_atoms(self, selected_atoms) -> None:
        """Hook 3: narrow the highlight from whole residues down to
        just the picked atoms once the user has named them. Same
        palette indexing as Hook 2 so colours line up between the two
        views.
        """
        viewer = self._viewer_or_none()
        if viewer is None:
            return
        # Drop the whole-residue reps so we don't show two rep layers
        # for the same residue.
        self._clear_pick_labels()
        applied = []
        for idx, atom in enumerate(selected_atoms, 1):
            chain = getattr(atom, "chain", None)
            resid = getattr(atom, "resid", None)
            atom_name = getattr(atom, "atom_name", None)
            if not chain or resid is None or not atom_name:
                continue
            label = f"{self._VIEWER_LABEL_PREFIX}atom_{idx}"
            try:
                viewer.highlight(
                    f":{chain} and {resid} and .{atom_name}",
                    style="ball+stick",
                    color=f"palette:{idx}",
                    label=label,
                )
                applied.append(label)
            except Exception as exc:
                logger.debug("atom highlight silenced: %s", exc)
        self._pick_labels = applied

    def _draw_restraint_overlay(self, restraint) -> None:
        """Hook 4: after a restraint is added, draw its geometric
        signature. Distance restraints become a yellow line via
        ``show_bonds``; angle/torsion keep the per-atom highlights
        from Hook 3 (a polyline visualiser would be useful but the
        coordinator's bond primitive only supports atom pairs, and
        chaining bonds would clutter the viewer for short angle
        restraints).
        """
        viewer = self._viewer_or_none()
        if viewer is None:
            return
        if getattr(restraint, "restraint_type", None) != "distance":
            return
        atoms = getattr(restraint, "atom_info", None) or []
        if len(atoms) < 2:
            return
        try:
            sel_a = self._atom_info_to_selection(atoms[0])
            sel_b = self._atom_info_to_selection(atoms[1])
            if not sel_a or not sel_b:
                return
            label = (
                f"{self._VIEWER_LABEL_PREFIX}bond_"
                f"{len(getattr(self, 'restraints', []) or [])}"
            )
            viewer.show_bonds(
                [(sel_a, sel_b)],
                label=label,
                color="#ffff00",
                show_labels=True,
            )
            existing = getattr(self, "_restraint_bond_labels", None) or []
            existing.append(label)
            self._restraint_bond_labels = existing
        except Exception as exc:
            logger.debug("restraint bond overlay silenced: %s", exc)

    @staticmethod
    def _atom_info_to_selection(atom_info) -> Optional[str]:
        """Build an NGL atom selector from a restraint's atom_info dict
        (or StructureAtom). Returns None if the required keys are
        missing.
        """
        if hasattr(atom_info, "chain"):
            chain = getattr(atom_info, "chain", None)
            resid = getattr(atom_info, "resid", None)
            atom_name = getattr(atom_info, "atom_name", None)
        elif isinstance(atom_info, dict):
            chain = atom_info.get("chain") or atom_info.get("chain_id")
            resid = atom_info.get("resid") or atom_info.get("resnum")
            atom_name = atom_info.get("atom_name") or atom_info.get("name")
        else:
            return None
        if not chain or resid is None or not atom_name:
            return None
        return f":{chain} and {resid} and .{atom_name}"

    def _highlight_residue_set(
        self,
        residue_pairs,
        *,
        label: str,
        color: str = "#1f78b4",
    ) -> None:
        """Generic helper used by Hooks 5 and 7 to halo a set of
        ``(chain, resid)`` pairs under a single label. Idempotent —
        re-firing replaces the prior selection.
        """
        viewer = self._viewer_or_none()
        if viewer is None:
            return
        viewer.unhighlight(label)
        clauses = [
            f"(:{c} and {r})"
            for c, r in sorted({(c, r) for c, r in residue_pairs if c and r is not None})
        ]
        if not clauses:
            return
        try:
            viewer.highlight(
                " or ".join(clauses),
                style="ball+stick",
                color=color,
                label=label,
            )
        except Exception as exc:
            logger.debug("residue-set highlight silenced (%s): %s", label, exc)

    def _clear_label_set(self, labels) -> None:
        viewer = self._viewer_or_none()
        if viewer is None or not labels:
            return
        for lbl in labels:
            viewer.unhighlight(lbl)


# Standalone integration functions
def add_restraints_to_redox_site(site, structure, console=None, processor=None):
    """
    Add MD restraints to a RedoxSite object (legacy compatibility)
    Call this after site detection and refinement
    
    Args:
        site: RedoxSite object
        structure: Bio.PDB Structure object
        console: Rich console for output
        processor: Main processor object for workspace access
    """
    restraint_manager = RestraintManager(console, processor)
    restraint_manager.current_structure = structure
    restraint_manager.current_structure_source = "redox_site"
    
    # Convert RedoxSite atoms to StructureAtoms
    structure_atoms = []
    for atom in site.atoms:
        structure_atoms.append(StructureAtom(
            chain=atom.chain,
            resname=atom.resname,
            resid=atom.resid,
            atom_name=atom.atom_name,
            coords=atom.coords,
            element=atom.element,
            insertion_code=atom.insertion_code,
            occupancy=atom.occupancy,
            bfactor=atom.bfactor
        ))
    
    restraint_manager._display_structure_overview(structure_atoms)
    
    # Run interactive restraint configuration
    while True:
        action = prompt_with_context(
            processor=restraint_manager.processor,
            prompt="\n[bold]Restraint Configuration Menu[/bold]\n"
            "[green]add[/green] - Add new restraint\n"
            "[yellow]list[/yellow] - Show current restraints\n"
            "[cyan]export[/cyan] - Export DISANG file\n"
            "[white]done[/white] - Finish restraint configuration\n"
            "Choose action",
            choices=["add", "list", "export", "done"],
            default="add",
            module="MD Restraint Manager",
            description="Select restraint action (standalone mode)",
            options_map={
                "add": "Add new restraint",
                "list": "Show current restraints",
                "export": "Export DISANG file",
                "done": "Finish restraint configuration"
            }
        )
        
        if action == "add":
            restraint_manager._add_restraint_interactive(structure_atoms, [])
        elif action == "list":
            restraint_manager._display_restraints()
        elif action == "export":
            restraint_manager._export_disang_file()
        elif action == "done":
            break
    
    # Attach restraints to site
    site.restraints = restraint_manager.restraints
    return site


@register_module
class MDRestraintModule(ProcessingModule):
    """MD Restraint Manager module for ProPrep structure preparation."""
    
    NAME = "MD Restraint Manager"
    DESCRIPTION = "Configure distance, angle, and torsion restraints for MD simulations"
    VERSION = "1.0.0"
    CATEGORY = "structure_prep"
    REQUIRES = ["PDB Loader"]
    PRIORITY = 6  # After redox site preparation
    
    def __init__(self):
        super().__init__()
        self.restraint_manager = None
        
    def initialize(self):
        """Initialize the module"""
        if hasattr(self, 'processor') and self.processor:
            self.restraint_manager = RestraintManager(
                console=self.processor.console,
                processor=self.processor
            )
        else:
            self.restraint_manager = RestraintManager()
    
    @property
    def console(self):
        """Get console from processor if available."""
        if hasattr(self, 'processor') and self.processor and hasattr(self.processor, 'console'):
            return self.processor.console
        else:
            from rich.console import Console
            return Console()
    
    def get_workspace_requirements(self) -> List[str]:
        """Get workspace requirements"""
        return ["structure"]  # Requires any PDB structure to be loaded

    def get_workspace_outputs(self) -> List[str]:
        """Get workspace outputs"""
        return [
            "disang_file",
            "disang_export_results",
            "redox_restraint_mask",
            "redox_restraint_info",
            "restraint_structure_source",
            "restraint_mask_generated",
            "restraints",
        ]

    def can_process(self, workspace) -> bool:
        """Check if module can process current workspace"""
        # Check if any structure is available from priority list
        priority_keys = [
            "protonation_pdb_file",          # protonation-updated
            "structure_with_prot_resnames",  # protonation-updated (legacy)
            "transformed_pdb_file",          # transformed
            "repaired_pdb_file",             # repaired
            "filtered_pdb_file",             # filtered
            "rcsb_pdb_file",                 # RCSB PDB
            "local_pdb_file",                # local PDB
            "alphafold_pdb_file",            # AlphaFold
            "alphafold_homolog_pdb_file",    # AlphaFold homolog
            "aligned_target_pdb_file",       # aligned target
            "aligned_ref_pdb_file",          # aligned reference
        ]

        # Check if any of the structure sources are available
        for key in priority_keys:
            # Check both attribute-style and dict-style access
            value = None
            if hasattr(workspace, key):
                value = getattr(workspace, key, None)
            elif hasattr(workspace, 'get'):
                value = workspace.get(key, None)
                
            if value is not None:
                # For file paths, verify they exist
                if key.endswith('_file') or key == 'pdb_file':
                    from pathlib import Path
                    if isinstance(value, str) and Path(value).exists():
                        return True
                else:
                    # For structure objects or other non-file values
                    return True
        
        return False

    def process(self, workspace):
        """
        Main module entry point for workspace processing.

        Presents an interactive menu with both DISANG and restraintmask options.
        """
        if not self.restraint_manager:
            self.initialize()

        from proprep.utils.enhanced_menu import EnhancedMenuDisplay
        from proprep.application.menu_commands import prompt_with_context

        console = self.processor.console if self.processor else None
        if not console:
            # Fallback if no console available
            return self.restraint_manager.add_restraints_to_structure(workspace)

        menu_display = EnhancedMenuDisplay(console)

        # Interactive menu loop
        while True:
            console.print("\n[bold cyan]═══ MD RESTRAINT CONFIGURATION ═══[/bold cyan]")

            # Get and display structure source
            structure_info = self.restraint_manager._get_priority_structure(workspace)
            if structure_info:
                structure_file, source_name = structure_info
                console.print(f"[grey50]Structure source: {source_name}[/grey50]\n")

                # Show the structure being restrained so the viewer is in
                # sync with what every option below will be picking
                # atoms from. Idempotent for the same path — only swaps
                # if a different priority structure has appeared.
                try:
                    from proprep.structure_prep.viewer_coordinator import (
                        viewer as _viewer,
                    )
                    _viewer.show_structure(structure_file)
                except Exception as exc:
                    logger.debug("MD restraint baseline silenced: %s", exc)

            # Display menu with enhanced options
            enhanced_options = self.get_enhanced_menu_options(workspace)
            for option in enhanced_options:
                menu_display.print_option(option)

            # Get suggestion
            suggestion = self.get_menu_suggestion(workspace)
            if suggestion:
                console.print(f"\n[cyan]→ Suggestion: {suggestion}[/cyan]")

            # Build choices list
            choices = ["1", "2", "3", "4", "5", "6", "7", "8", "d", "x"]

            # Build options map for context
            options_map = {
                "1": "Import DISANG restraints from file",
                "2": "Configure DISANG restraints interactively",
                "3": "Display current DISANG restraints",
                "4": "Export restraints to DISANG file",
                "5": "Generate restraint mask for MD minimization",
                "6": "Configure GROUP specification",
                "7": "Display current GROUP specification",
                "8": "Export GROUP specification to file",
                "d": "Done - Return to MD Manager",
                "x": "Exit ProPrep"
            }

            choice = prompt_with_context(
                processor=self.processor,
                prompt="\nChoose option (or 'd' when done, 'x' to exit)",
                choices=choices,
                default="2",
                module="MD Restraint Manager",
                description="Select restraint configuration option",
                options_map=options_map
            )

            # Handle choice
            if choice == "d":
                # Done - return to caller
                console.print("\n[green]Restraint configuration complete[/green]")
                return workspace
            elif choice == "x":
                # Exit program
                console.print("[yellow]Exiting ProPrep...[/yellow]")
                return None
            elif choice == "1":
                self.handle_menu_option("import")
            elif choice == "2":
                self.handle_menu_option("configure")
            elif choice == "3":
                self.handle_menu_option("display")
            elif choice == "4":
                self.handle_menu_option("export")
            elif choice == "5":
                self.handle_menu_option("mask")
            elif choice == "6":
                self.handle_menu_option("group_configure")
            elif choice == "7":
                self.handle_menu_option("group_display")
            elif choice == "8":
                self.handle_menu_option("group_export")

    def get_menu_options(self) -> Dict[str, str]:
        """Get available menu options"""
        return {
            "import": "Import DISANG restraints from file",
            "configure": "Configure DISANG restraints interactively",
            "display": "Display current DISANG restraints",
            "export": "Export restraints to DISANG file",
            "mask": "Generate restraint mask for MD minimization",
            "group_configure": "Configure GROUP specification",
            "group_display": "Display current GROUP specification",
            "group_export": "Export GROUP specification to file",
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
        restraints = workspace.get('md_restraints', [])
        has_restraints = len(restraints) > 0
        has_mask = workspace.get('redox_restraint_mask') is not None

        # DISANG Restraints Section Header
        options.append(MenuOption(
            key="",
            description="DISANG Restraints (Distance/Angle/Torsion)",
            is_separator=True
        ))

        # Option 1: Import DISANG - always available
        options.append(MenuOption(
            key="1",
            description="Import DISANG restraints from file",
            status=OptionStatus.AVAILABLE
        ))

        # Option 2: Configure DISANG - needs a loaded structure; ● once restraints exist
        if has_restraints:
            disang_status, disang_dep = OptionStatus.COMPLETED, ""
        elif self.can_process(workspace):
            disang_status, disang_dep = OptionStatus.AVAILABLE, ""
        else:
            disang_status = OptionStatus.BLOCKED
            disang_dep = self.availability_note(workspace) or "Load a structure first"
        options.append(MenuOption(
            key="2",
            description="Configure DISANG restraints interactively",
            status=disang_status,
            dependency_text=disang_dep,
        ))

        # Option 3: Display DISANG - requires restraints
        if has_restraints:
            status = OptionStatus.READY
            dep_text = ""
        else:
            status = OptionStatus.BLOCKED
            dep_text = "[Need to configure or import restraints first] ○"

        options.append(MenuOption(
            key="3",
            description="Display current DISANG restraints",
            status=status,
            dependency_text=dep_text
        ))

        # Option 4: Export DISANG - requires restraints
        if has_restraints:
            status = OptionStatus.READY
            dep_text = ""
        else:
            status = OptionStatus.BLOCKED
            dep_text = "[Need to configure or import restraints first] ○"

        options.append(MenuOption(
            key="4",
            description="Export restraints to DISANG file",
            status=status,
            dependency_text=dep_text
        ))

        # Positional Restraints Section Header
        options.append(MenuOption(
            key="",
            description="Positional Restraints",
            is_separator=True
        ))

        # Option 5: Generate restraintmask - needs a loaded structure
        if has_mask:
            mask_status, mask_dep = OptionStatus.COMPLETED, ""
        elif self.can_process(workspace):
            mask_status, mask_dep = OptionStatus.AVAILABLE, ""
        else:
            mask_status = OptionStatus.BLOCKED
            mask_dep = self.availability_note(workspace) or "Load a structure first"

        options.append(MenuOption(
            key="5",
            description="Generate restraint mask (simple positional restraints)",
            status=mask_status,
            dependency_text=mask_dep,
        ))

        # Check for GROUP specification in workspace
        has_group = workspace.get('group_restraints') is not None

        # Option 6: Configure GROUP specification - always available
        if has_group:
            status = OptionStatus.COMPLETED
        else:
            status = OptionStatus.AVAILABLE

        options.append(MenuOption(
            key="6",
            description="Configure GROUP specification",
            status=status
        ))

        # Option 7: Display GROUP - requires GROUP configuration
        if has_group:
            status = OptionStatus.READY
            dep_text = ""
        else:
            status = OptionStatus.BLOCKED
            dep_text = "[Need to configure GROUP first] ○"

        options.append(MenuOption(
            key="7",
            description="Display current GROUP specification",
            status=status,
            dependency_text=dep_text
        ))

        # Option 8: Export GROUP - requires GROUP configuration
        if has_group:
            status = OptionStatus.READY
            dep_text = ""
        else:
            status = OptionStatus.BLOCKED
            dep_text = "[Need to configure GROUP first] ○"

        options.append(MenuOption(
            key="8",
            description="Export GROUP specification to file",
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
        restraints = workspace.get('md_restraints', [])
        has_restraints = len(restraints) > 0
        has_mask = workspace.get('redox_restraint_mask') is not None

        # Check what's been done
        has_group = workspace.get('group_restraints') is not None

        if not has_restraints and not has_mask and not has_group:
            return "Import DISANG restraints (option 1), configure new ones (option 2), generate restraintmask (option 5), or configure GROUP specification (option 6)"
        elif has_restraints and not has_mask:
            return f"✓ {len(restraints)} DISANG restraint(s) configured. Display (option 3), export (option 4), or generate mask (option 5)"
        elif not has_restraints and has_mask:
            return "✓ Positional restraint mask generated. Import DISANG restraints (option 1) or configure new ones (option 2) if needed"
        else:
            return f"✓ {len(restraints)} DISANG restraint(s) configured and mask generated. Display (option 3), export (option 4), or continue to next module"
    
    def handle_menu_option(self, option: str) -> bool:
        """Handle menu option selection using command pattern"""
        if not self.restraint_manager:
            self.initialize()
            
        if option == "configure":
            command = ConfigureRestraintsCommand(self.processor)
            return command.execute()
        elif option == "display":
            command = DisplayRestraintsCommand(self.processor)
            return command.execute()
        elif option == "import":
            command = ImportRestraintsCommand(self.processor)
            return command.execute()
        elif option == "export":
            command = ExportDisangCommand(self.processor)
            return command.execute()
        elif option == "mask":
            command = GenerateRestraintMaskCommand(self.processor)
            return command.execute()
        elif option == "group_configure":
            return self._configure_group_restraints()
        elif option == "group_display":
            return self._display_group_restraints()
        elif option == "group_export":
            return self._export_group_restraints()

        return False
    
    def _configure_restraints(self) -> bool:
        """Configure restraints using the RestraintManager"""
        try:
            workspace = self.processor.workspace
            result = self.restraint_manager.add_restraints_to_structure(workspace)
            return result is not None
        except Exception as e:
            logger.error(f"Error configuring restraints: {e}")
            self.console.print(f"[red]Error: {e}[/red]")
            return False
    
    def _display_restraints(self) -> bool:
        """Display current restraints"""
        try:
            if not self.restraint_manager:
                self.initialize()
            
            # Get restraints from workspace using the correct access pattern
            workspace = self.processor.workspace
            restraints = workspace.get('md_restraints', []) if workspace else []
            
            if restraints:
                # Set restraints in manager and display them
                self.restraint_manager.restraints = restraints
                self.restraint_manager._display_restraints()
            else:
                self.console.print("[yellow]No restraints configured yet[/yellow]")
            return True
        except Exception as e:
            logger.error(f"Error displaying restraints: {e}")
            self.console.print(f"[red]Error: {e}[/red]")
            return False
    
    def _import_restraints(self) -> bool:
        """Import restraints from DISANG file"""
        try:
            if not self.restraint_manager:
                self.initialize()
            
            # Find available DISANG files
            import glob
            import os

            disang_files = []
            for pattern in ["*.disang", "*restraint*.disang", "*_restraints.disang"]:
                disang_files.extend(glob.glob(pattern))
            
            # Remove duplicates and sort
            disang_files = sorted(list(set(disang_files)))
            
            if disang_files:
                self.console.print(f"\n[cyan]Available DISANG files in current directory:[/cyan]")
                for i, filename in enumerate(disang_files, 1):
                    # Get file size and modification time for info
                    stat = os.stat(filename)
                    size_kb = stat.st_size / 1024
                    from datetime import datetime
                    mod_time = datetime.fromtimestamp(stat.st_mtime).strftime('%Y-%m-%d %H:%M')
                    self.console.print(f"  {i}. {filename} ({size_kb:.1f} KB, modified {mod_time})")
                
                # Prompt with suggestion
                default_file = disang_files[0]
                self.console.print(f"\n[grey50]Suggested: {default_file}[/grey50]")
            else:
                self.console.print("[yellow]No DISANG files found in current directory[/yellow]")
                default_file = "restraints.disang"
            
            filename = prompt_with_context(
                processor=self.processor,
                prompt="Enter DISANG filename to import",
                default=default_file if disang_files else "restraints.disang",
                module="MD Restraint Manager",
                description="Enter DISANG import filename"
            )
            
            if not os.path.exists(filename):
                self.console.print(f"[red]✗ File not found: {filename}[/red]")
                return False
            
            # Parse and import restraints
            imported_restraints = self._parse_disang_file(filename)
            
            if imported_restraints:
                # Add to existing restraints
                self.restraint_manager.restraints.extend(imported_restraints)
                
                # Update workspace
                workspace_restraints = self.get_from_workspace("restraints", [])
                workspace_restraints.extend([r.__dict__ for r in imported_restraints])
                self.update_workspace("restraints", workspace_restraints)
                
                self.console.print(f"[green]✓ Successfully imported {len(imported_restraints)} restraints from {filename}[/green]")
                
                # Display summary
                types_count = {}
                for restraint in imported_restraints:
                    types_count[restraint.restraint_type] = types_count.get(restraint.restraint_type, 0) + 1
                
                for rtype, count in types_count.items():
                    self.console.print(f"  - {count} {rtype} restraint(s)")
                
                return True
            else:
                self.console.print(f"[yellow]No restraints found in {filename}[/yellow]")
                return False
                
        except Exception as e:
            self.console.print(f"[red]✗ Error importing restraints: {str(e)}[/red]")
            return False
    
    def _parse_disang_file(self, filename: str) -> List:
        """Parse DISANG file and return list of MDRestraint objects"""
        from dataclasses import field
        restraints = []
        
        try:
            with open(filename, 'r') as f:
                content = f.read()
            
            # Split into sections by &rst.../ blocks
            import re
            rst_blocks = re.findall(r'&rst\s*(.*?)\s*/', content, re.DOTALL | re.MULTILINE)
            
            # Find comment blocks before each &rst
            comment_blocks = re.split(r'&rst\s*.*?\s*/', content)[:-1]  # Remove last empty part
            
            for i, (rst_block, comment_block) in enumerate(zip(rst_blocks, comment_blocks)):
                try:
                    restraint = self._parse_single_restraint(rst_block, comment_block)
                    if restraint:
                        restraints.append(restraint)
                except Exception as e:
                    self.console.print(f"[yellow]Warning: Failed to parse restraint {i+1}: {str(e)}[/yellow]")
                    continue
            
            return restraints
            
        except Exception as e:
            self.console.print(f"[red]Error reading DISANG file: {str(e)}[/red]")
            return []
    
    def _parse_single_restraint(self, rst_block: str, comment_block: str):
        """Parse a single restraint from &rst block and associated comments"""
        import re
        from dataclasses import field
        
        # Extract restraint type from comments
        type_match = re.search(r'# Restraint \d+: (\w+)', comment_block)
        restraint_type = type_match.group(1) if type_match else "distance"
        
        # Extract atom information from comments
        atom_info = []
        atom_matches = re.findall(r'# Atom \d+: (\w+) ([A-Z]):([\w\d]+) @(\w+)', comment_block)
        for resname, chain, resid, atom_name in atom_matches:
            atom_info.append({
                'resname': resname,
                'chain': chain, 
                'resid': resid,
                'atom_name': atom_name,
                'element': 'C'  # Default, will be updated when structure is loaded
            })
        
        # Extract current value from comments
        current_value = None
        value_match = re.search(r'# Current value: ([\d.]+)', comment_block)
        if value_match:
            current_value = float(value_match.group(1))
        
        # Extract description from comments
        description = ""
        desc_match = re.search(r'# Description: (.+)', comment_block)
        if desc_match:
            description = desc_match.group(1).strip()
        
        # Parse parameters from &rst block
        params = {}
        
        # Extract iat parameter (amber indices)
        iat_match = re.search(r'iat\s*=\s*([\d,\s]+)', rst_block)
        amber_index = []
        if iat_match:
            # Parse comma-separated indices, handling trailing comma
            iat_str = iat_match.group(1).strip().rstrip(',')
            if iat_str:
                amber_index = [int(x.strip()) for x in iat_str.split(',') if x.strip()]
        
        # Extract other parameters
        param_patterns = {
            'r0': r'r0\s*=\s*([\d.]+)',
            'k0': r'k0\s*=\s*([\d.]+)', 
            'r1': r'r1\s*=\s*([\d.]+)',
            'r2': r'r2\s*=\s*([\d.]+)',
            'r3': r'r3\s*=\s*([\d.]+)',
            'r4': r'r4\s*=\s*([\d.]+)',
            'rk2': r'rk2\s*=\s*([\d.]+)',
            'rk3': r'rk3\s*=\s*([\d.]+)',
            'ifvari': r'ifvari\s*=\s*(\d+)'
        }
        
        for param, pattern in param_patterns.items():
            match = re.search(pattern, rst_block)
            if match:
                if param == 'ifvari':
                    params[param] = int(match.group(1))
                else:
                    params[param] = float(match.group(1))
        
        # Create restraint object
        restraint = MDRestraint(
            restraint_type=restraint_type,
            atom_coords=[],  # Will be filled when structure is loaded
            atom_info=atom_info,
            current_value=current_value,
            description=description,
            amber_index=amber_index,
            # Set parameters with defaults
            r0=params.get('r0'),
            k0=params.get('k0', 10.0),
            r1=params.get('r1', 0.0),
            r2=params.get('r2', 0.0), 
            r3=params.get('r3', 0.0),
            r4=params.get('r4', 0.0),
            rk2=params.get('rk2', 10.0),
            rk3=params.get('rk3', 10.0),
            ifvari=params.get('ifvari', 0),
            active=True
        )
        
        return restraint
    
    def _export_disang(self) -> bool:
        """Export restraints to DISANG file"""
        try:
            if not self.restraint_manager:
                self.initialize()
            
            # Get restraints from workspace using the correct access pattern
            workspace = self.processor.workspace
            restraints = workspace.get('md_restraints', []) if workspace else []
            
            if restraints:
                # Set restraints in manager and export them
                self.restraint_manager.restraints = restraints
                self.restraint_manager._export_disang_file()
                return True
            else:
                self.console.print("[yellow]No restraints to export. Configure restraints first.[/yellow]")
            return True
        except Exception as e:
            logger.error(f"Error exporting DISANG file: {e}")
            self.console.print(f"[red]Error: {e}[/red]")
            return False

    def _generate_restraint_mask(self) -> bool:
        """Generate restraint mask using the RestraintManager"""
        try:
            if not self.restraint_manager:
                self.initialize()

            workspace = self.processor.workspace
            result = self.restraint_manager.generate_restraint_mask(workspace, interactive=True)
            return result is not None
        except Exception as e:
            logger.error(f"Error generating restraint mask: {e}")
            self.console.print(f"[red]Error: {e}[/red]")
            return False

    def _configure_group_restraints(self) -> bool:
        """Configure GROUP specification restraints interactively"""
        try:
            workspace = self.processor.workspace

            # Introduction
            self.console.print(Panel.fit(
                "[bold cyan]GROUP Specification Configuration[/bold cyan]\n\n"
                "GROUP specifications allow multiple restraint groups with:\n"
                "  • Different force constants per group\n"
                "  • Atom filtering via FIND criteria\n"
                "  • Complex residue selections\n\n"
                "Example: Restrain CA atoms of terminal chains (weak) and metal\n"
                "sites (strong) with different force constants.\n\n"
                "Note: Cannot be combined with restraintmask in the same simulation step.",
                border_style="cyan"
            ))

            if not confirm_with_context(
                self.processor,
                "\nContinue with GROUP configuration?",
                default=True,
                module="MD Restraint Manager",
                description="Continue with GROUP-restraint configuration",
            ):
                return False

            # Ask number of groups
            num_groups = int_prompt_with_context(
                self.processor,
                "How many restraint groups to define?",
                default=1,
                module="MD Restraint Manager",
                description="Number of GROUP restraint groups (1-10)",
            )
            if num_groups < 1 or num_groups > 10:
                self.console.print("[red]Number of groups must be between 1 and 10[/red]")
                return False

            groups = []

            # Configure each group
            for i in range(num_groups):
                self.console.print(f"\n{'─' * 63}")
                self.console.print(f"Configuring Group {i+1} of {num_groups}")
                self.console.print(f"{'─' * 63}\n")

                # Group title
                title = prompt_with_context(
                    self.processor,
                    "Group title/description",
                    default=f"Group {i+1}",
                    module="MD Restraint Manager",
                    description=f"Title for GROUP restraint {i+1}",
                )

                # Force constant
                force_constant = float_prompt_with_context(
                    self.processor,
                    "Force constant (kcal/mol/Å²)",
                    default=10.0,
                    module="MD Restraint Manager",
                    description=f"Force constant for GROUP restraint {i+1}",
                )

                # FIND criteria
                self.console.print(f"\n{'-' * 45}")
                self.console.print("FIND Criteria (Atom Selection)")
                self.console.print(f"{'-' * 45}\n")

                self.console.print("FIND format: atom_name  atom_type  tree_name  residue_name")
                self.console.print("  • Use * as wildcard for any field")
                self.console.print("  • All 4 fields required\n")

                self.console.print("Field explanations:")
                self.console.print("  atom_name   : PDB atom name (e.g., CA, CB, FE, N)")
                self.console.print("  atom_type   : AMBER atom type (rarely used, usually *)")
                self.console.print("  tree_name   : Tree structure (M=main, S=side, B=both, usually *)")
                self.console.print("  residue_name: PDB residue name (e.g., ALA, HEM, CYS)\n")

                self.console.print("Common examples:")
                self.console.print("  CA * * *     → All CA (alpha carbon) atoms")
                self.console.print("  * * * HEM    → All atoms in HEM residues")
                self.console.print("  FE * * *     → All iron atoms")
                self.console.print("  N * M *      → All backbone N atoms (main chain)")
                self.console.print("  CB * S *     → All CB (beta carbon) atoms (side chain)\n")

                num_find = int_prompt_with_context(
                    self.processor,
                    "How many FIND criteria for this group?",
                    default=1,
                    module="MD Restraint Manager",
                    description=f"Number of FIND criteria for GROUP restraint {i+1}",
                )
                if num_find < 1 or num_find > 5:
                    self.console.print("[red]Number of FIND criteria must be between 1 and 5[/red]")
                    return False

                find_criteria = []
                for j in range(num_find):
                    criterion = prompt_with_context(
                        self.processor,
                        f"FIND criterion {j+1}",
                        module="MD Restraint Manager",
                        description=f"FIND criterion {j+1} for GROUP restraint {i+1}",
                    )
                    find_criteria.append(criterion)

                # Residue ranges
                self.console.print(f"\n{'-' * 45}")
                self.console.print("Residue Selection")
                self.console.print(f"{'-' * 45}\n")

                self.console.print("Enter residue ranges for this group:")
                self.console.print("  • Format: start-end (e.g., 1-198)")
                self.console.print("  • Separate multiple ranges with commas or spaces")
                self.console.print("  • Maximum 7 ranges per group\n")

                # Loop until we parse a valid set of ranges — one typo shouldn't
                # discard the group's title / force constant / FIND criteria.
                residue_ranges = []
                while True:
                    ranges_input = prompt_with_context(
                        self.processor,
                        "Residue ranges",
                        module="MD Restraint Manager",
                        description=f"Residue ranges for GROUP restraint {i+1}",
                    )
                    # Normalize Unicode dashes (en-dash, em-dash, minus sign) that
                    # shells, editors, or macOS text-substitution often produce in
                    # place of an ASCII hyphen.
                    normalized = (
                        ranges_input
                        .replace('\u2013', '-')  # en dash
                        .replace('\u2014', '-')  # em dash
                        .replace('\u2212', '-')  # minus sign
                    )
                    residue_ranges = []
                    try:
                        for range_str in normalized.replace(',', ' ').split():
                            range_str = range_str.strip()
                            if '-' in range_str:
                                start, end = range_str.split('-')
                                residue_ranges.append((int(start.strip()), int(end.strip())))
                            else:
                                res = int(range_str)
                                residue_ranges.append((res, res))
                        break
                    except ValueError:
                        self.console.print(
                            f"[red]Invalid range format: {ranges_input!r}[/red]  "
                            f"[grey50](expected: start-end, e.g. '1-198 793-990')[/grey50]"
                        )

                if len(residue_ranges) > 7:
                    self.console.print("[yellow]Warning: Maximum 7 ranges per group. Using first 7.[/yellow]")
                    residue_ranges = residue_ranges[:7]

                # Add group
                group = {
                    "title": title,
                    "force_constant": force_constant,
                    "find_criteria": find_criteria,
                    "residue_ranges": residue_ranges
                }
                groups.append(group)

                # Display confirmation
                self.console.print(f"\n✓ Group {i+1} configured:")
                self.console.print(f"  Title: {title}")
                self.console.print(f"  Force constant: {force_constant} kcal/mol/Å²")
                self.console.print(f"  FIND Criteria:")
                for criterion in find_criteria:
                    self.console.print(f"    - {criterion}")
                range_strs = [f"{start}-{end}" if start != end else str(start) for start, end in residue_ranges]
                self.console.print(f"  Residues: {', '.join(range_strs)}")

            # Preview
            self.console.print(f"\n{'─' * 63}")
            self.console.print("Preview: GROUP Specification Output")
            self.console.print(f"{'─' * 63}\n")

            self.console.print("This will be added to mdin files after the &cntrl namelist:\n")

            preview_lines = []
            for group in groups:
                preview_lines.append(group["title"])
                preview_lines.append(str(group["force_constant"]))
                if group["find_criteria"]:
                    preview_lines.append("FIND")
                    for criterion in group["find_criteria"]:
                        preview_lines.append(criterion)
                    preview_lines.append("SEARCH")

                # Format RES line
                res_line = "RES"
                for start, end in group["residue_ranges"]:
                    res_line += f" {start} {end}"
                preview_lines.append(res_line)
                preview_lines.append("END")

            preview_lines.append("END")

            for line in preview_lines:
                self.console.print(line)

            self.console.print("\nNote: ntr=1 will be automatically added to &cntrl namelist\n")

            # Confirm save
            if not confirm_with_context(
                self.processor,
                "Save this GROUP specification?",
                default=True,
                module="MD Restraint Manager",
                description="Save GROUP specification to workspace",
            ):
                return False

            # Save to workspace
            workspace['group_restraints'] = groups

            self.console.print("\n[green]✓ GROUP specification saved to workspace[/green]")
            self.console.print("  • Available in MD workflow generation")
            self.console.print("  • View with option 7")
            self.console.print("  • Export with option 8")

            return True

        except Exception as e:
            logger.error(f"Error configuring GROUP restraints: {e}")
            self.console.print(f"[red]Error: {e}[/red]")
            import traceback
            traceback.print_exc()
            return False

    def _display_group_restraints(self) -> bool:
        """Display current GROUP specification"""
        try:
            workspace = self.processor.workspace
            groups = workspace.get('group_restraints')

            if not groups:
                self.console.print("[yellow]No GROUP specification configured yet[/yellow]")
                return False

            self.console.print("\n" + "=" * 60)
            self.console.print("Current GROUP Specification")
            self.console.print("=" * 60 + "\n")

            for i, group in enumerate(groups, 1):
                self.console.print(f"Group {i}: {group['title']}")
                self.console.print(f"  Force constant: {group['force_constant']} kcal/mol/Å²")
                self.console.print("  FIND Criteria:")
                for criterion in group['find_criteria']:
                    self.console.print(f"    - {criterion}")
                range_strs = [f"{start}-{end}" if start != end else str(start)
                             for start, end in group['residue_ranges']]
                self.console.print(f"  Residues: {', '.join(range_strs)}")
                self.console.print()

            prompt_with_context(
                self.processor,
                "\nPress Enter to continue",
                default="",
                module="MD Restraint Manager",
                description="Pause after viewing GROUP spec",
            )
            return True

        except Exception as e:
            logger.error(f"Error displaying GROUP restraints: {e}")
            self.console.print(f"[red]Error: {e}[/red]")
            return False

    def _export_group_restraints(self) -> bool:
        """Export GROUP specification to file"""
        try:
            workspace = self.processor.workspace
            groups = workspace.get('group_restraints')

            if not groups:
                self.console.print("[yellow]No GROUP specification configured. Configure with option 6 first.[/yellow]")
                return False

            # Ask for filename
            filename = prompt_with_context(
                self.processor,
                "Export filename",
                default="group_restraints.txt",
                module="MD Restraint Manager",
                description="Output filename for GROUP specification export",
            )

            # Generate GROUP specification
            lines = []
            for group in groups:
                lines.append(group["title"])
                lines.append(str(group["force_constant"]))
                if group["find_criteria"]:
                    lines.append("FIND")
                    for criterion in group["find_criteria"]:
                        lines.append(criterion)
                    lines.append("SEARCH")

                # Format RES line
                res_line = "RES"
                for start, end in group["residue_ranges"]:
                    res_line += f" {start} {end}"
                lines.append(res_line)
                lines.append("END")

            lines.append("END")

            # Write to file
            with open(filename, 'w') as f:
                f.write('\n'.join(lines) + '\n')

            self.console.print(f"[green]✓ GROUP specification exported to {filename}[/green]")
            return True

        except Exception as e:
            logger.error(f"Error exporting GROUP specification: {e}")
            self.console.print(f"[red]Error: {e}[/red]")
            return False