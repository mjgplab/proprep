"""
Structure Alignment Module

Progressive superimposition of multiple PDB structures onto a reference structure
using specified residues. Integrates with RedoxSite objects for automatic residue
selection from detected redox-active sites.
"""

import os
import numpy as np
import logging
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Any
from dataclasses import dataclass

from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Superimposer, Structure, Residue
from Bio.PDB.Atom import Atom

# Optional TM-align support for structure-based alignment
try:
    import tmtools
    TMTOOLS_AVAILABLE = True
except ImportError:
    TMTOOLS_AVAILABLE = False
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich import box

from proprep.utils.module_registry import ProcessingModule, register_module
from proprep.utils.prompts import prompt_with_context, int_prompt_with_context, confirm_with_context
from proprep.structure_prep.structure_alignment_commands import (
    AlignStructuresCommand,
    AlignOnRedoxSitesCommand,
    ViewAlignmentResultsCommand,
)

logger = logging.getLogger(__name__)


@dataclass
class ResidueSpec:
    """Specification for a residue in alignment"""
    chain_id: str
    resname: str
    resid: int

    def __str__(self):
        return f"{self.chain_id}:{self.resname}:{self.resid}"


@dataclass
class ResidueMapping:
    """Mapping between reference residue and structure-specific residues"""
    reference: ResidueSpec
    mappings: Dict[int, ResidueSpec]  # structure_index -> ResidueSpec

    def get_spec_for_structure(self, struct_idx: int, is_reference: bool) -> ResidueSpec:
        """Get the appropriate residue spec for a given structure"""
        if is_reference:
            return self.reference
        return self.mappings.get(struct_idx, self.reference)


class AlignmentMode:
    """Defines which atoms to use for alignment"""
    FULL = "full"
    BACKBONE = "backbone"
    ALPHA_CARBON = "alpha_carbon"

    BACKBONE_ATOMS = {'N', 'CA', 'C', 'O'}

    @staticmethod
    def get_atoms(residue: Residue, mode: str) -> List[Atom]:
        """Extract atoms from residue based on alignment mode"""
        if mode == AlignmentMode.ALPHA_CARBON:
            try:
                ca = residue['CA']
                return [ca]
            except KeyError:
                return []
        elif mode == AlignmentMode.BACKBONE:
            atoms = []
            for atom_name in AlignmentMode.BACKBONE_ATOMS:
                try:
                    atom = residue[atom_name]
                    atoms.append(atom)
                except KeyError:
                    pass
            return atoms
        else:  # FULL
            return list(residue.get_atoms())


class MappingStrategy:
    """Defines the method for mapping residues between structures"""
    SEQUENCE_BASED = "sequence_based"
    STRUCTURE_BASED = "structure_based"
    MANUAL = "manual"


@register_module
class StructureAlignmentModule(ProcessingModule):
    """Module for progressive structure superimposition"""

    NAME = "Structure Aligner"
    DESCRIPTION = "Align structures by sequence or fold"
    VERSION = "1.0.0"
    CATEGORY = "analysis"
    PRIORITY = 10  # Run after structure prep modules

    def __init__(self):
        """Initialize the module"""
        super().__init__()

    def initialize(self):
        """Initialize module resources"""
        self.structures = []  # List of (filename, structure) tuples
        self.aligned_structures = {}  # Dict mapping index to aligned structure
        self.alignment_results = []  # List of alignment RMSD results
        self.reference_idx = None
        self.alignment_mode = None
        self.residue_mappings = []
        self.save_intermediates = False
        self.hetatm_to_add = []  # List of (ref_residue, target_idx_list) for HETATMs to add after alignment
        self.final_transformation_matrices = {}  # Dict mapping target_idx to final transformation matrix

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

    def get_menu_options(self) -> Dict[str, str]:
        """Get module menu options"""
        return {
            "align_structures": "Align multiple PDB structures",
            "align_on_redox_sites": "Align structures using redox site residues",
            "view_results": "View alignment results",
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
        detected_redox_sites = workspace.get("detected_redox_sites")
        alignment_results = workspace.get("alignment_results")
        has_redox_sites = detected_redox_sites is not None and len(detected_redox_sites) > 0

        # Option 1: Align structures - needs a loaded structure; ● once aligned
        if alignment_results:
            align_status, align_dep = OptionStatus.COMPLETED, ""
        elif self.can_process(workspace):
            align_status, align_dep = OptionStatus.AVAILABLE, ""
        else:
            align_status = OptionStatus.BLOCKED
            align_dep = self.availability_note(workspace) or "Load a structure first"
        options.append(MenuOption(
            key="1",
            description="Align multiple PDB structures",
            status=align_status,
            dependency_text=align_dep,
        ))

        # Option 2: Align on redox sites - requires detected redox sites; ● once aligned
        if alignment_results:
            status = OptionStatus.COMPLETED
            dep_text = ""
        elif has_redox_sites:
            status = OptionStatus.AVAILABLE
            dep_text = ""
        else:
            status = OptionStatus.BLOCKED
            dep_text = "[Detect redox sites first] ○"

        options.append(MenuOption(
            key="2",
            description="Align structures using redox site residues",
            status=status,
            dependency_text=dep_text
        ))

        # Option 3: View results - requires alignment to be done
        if alignment_results:
            status = OptionStatus.READY
            dep_text = ""
        else:
            status = OptionStatus.BLOCKED
            dep_text = "[Need to align structures first] ○"

        options.append(MenuOption(
            key="3",
            description="View alignment results",
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
        detected_redox_sites = workspace.get("detected_redox_sites")
        alignment_results = workspace.get("alignment_results")
        has_redox_sites = detected_redox_sites is not None and len(detected_redox_sites) > 0

        if not alignment_results:
            if not self.can_process(workspace):
                return f"{self.availability_note(workspace) or 'A structure is required'}. Load one via the Structure Loader."
            if has_redox_sites:
                return "Align structures manually (option 1) or use detected redox sites (option 2)"
            else:
                return "Start by aligning structures (option 1)"
        else:
            return "View alignment results with option 3, or press [m] to return to the main menu"

    def handle_menu_option(self, option: str) -> bool:
        """Handle menu option selection using command pattern"""
        try:
            if option == "align_structures":
                command = AlignStructuresCommand(self.processor)
                return command.execute_with_error_handling()
            elif option == "align_on_redox_sites":
                command = AlignOnRedoxSitesCommand(self.processor)
                return command.execute_with_error_handling()
            elif option == "view_results":
                command = ViewAlignmentResultsCommand(self.processor)
                return command.execute_with_error_handling()
        except Exception as e:
            import traceback
            logger.error(f"Error executing menu option '{option}': {e}")
            logger.error(traceback.format_exc())
            self.console.print(f"[red]Error: {e}[/red]")
            self.console.print(f"[grey50]{traceback.format_exc()}[/grey50]")

        return False

    def get_workspace_requirements(self) -> List[str]:
        """Get workspace requirements - need at least one structure loaded"""
        return ["rcsb_pdb_file | local_pdb_file | alphafold_pdb_file | alphafill_pdb_file | alphafold_homolog_pdb_file"]

    def get_workspace_outputs(self) -> List[str]:
        """Get workspace outputs - aligned structures"""
        return [
            "aligned_ref_pdb_file",      # Reference aligned structure file
            "aligned_ref_structure",     # Reference aligned structure object
            "aligned_target_pdb_file",   # Selected target aligned structure file
            "aligned_target_structure",  # Selected target aligned structure object
            "aligned_structures",        # Dict of all aligned structures (index -> structure)
            "alignment_results",         # List of alignment RMSD results
            "alignment_residues",        # Residue mappings used for alignment
        ]

    def can_process(self, workspace: Dict[str, Any]) -> bool:
        """Check if module can process current workspace.

        Uses StructureSelector to check for any available structure.
        """
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(workspace, self.console)
        status = selector.get_structure_status()
        return status.get("has_any", False)

    def process(self, workspace: Dict[str, Any]) -> Dict[str, Any]:
        """Process the workspace"""
        # This module is primarily interactive, no automatic processing
        return workspace

    def _align_structures_interactive(self, workspace: Dict[str, Any] = None) -> Dict[str, Any]:
        """Interactive structure alignment workflow"""
        if workspace is None:
            workspace = self.processor._get_workspace()

        self.console.print(Panel.fit(
            "[bold cyan]Structure Alignment Tool[/bold cyan]\n"
            "Progressive superimposition with cumulative residue addition",
            border_style="cyan"
        ))

        # Step 1: Load structures
        self.console.print("\n[bold]Step 1: Load Structures[/bold]")
        self.structures = self._load_structures_from_input(workspace)

        if len(self.structures) < 2:
            self.console.print("[red]Error: At least 2 structures required for alignment[/red]")
            return workspace

        # Step 2: Select reference
        self.console.print("\n[bold]Step 2: Select Reference Structure[/bold]")
        self.reference_idx = self._select_reference()

        # Step 3: Select residues for alignment
        self.console.print("\n[bold]Step 3: Select Residues for Alignment[/bold]")
        use_all_residues = confirm_with_context(
            processor=self.processor,
            prompt="Use all common residues for alignment?",
            default=True,
            module="Structure Alignment",
            description="Use all common residues vs. subset"
        )

        if use_all_residues:
            # Extract all residues from reference structure
            self.residue_mappings = self._extract_all_residues_from_reference()
            self.console.print(f"\n[green]Extracted {len(self.residue_mappings)} residue(s) from reference structure[/green]")
        else:
            # User specifies subset of residues
            self.console.print("\n[cyan]Specify subset of residues to use for alignment[/cyan]")
            self.residue_mappings = self._get_residue_subset_interactive()

        if not self.residue_mappings:
            self.console.print("[red]Error: No residues specified[/red]")
            return workspace

        # Step 4: Residue Mapping Strategy
        self.console.print("\n[bold]Step 4: Residue Mapping Strategy[/bold]")
        mapping_strategy = self._select_mapping_strategy()

        # Handle structure-based alignment separately (uses TM-align directly)
        if mapping_strategy == MappingStrategy.STRUCTURE_BASED:
            # Structure-based alignment handles its own superposition
            self._perform_structure_based_alignment()

            # Save aligned structures
            self.console.print("\n[bold]Saving Aligned Structures...[/bold]")
            self._save_results_to_workspace(workspace)

            return workspace

        # For sequence-based and manual strategies, continue with residue-based alignment
        # Step 5: Select alignment mode (which atoms within each residue)
        self.console.print("\n[bold]Step 5: Select Atom Mode[/bold]")
        self.alignment_mode = self._select_alignment_mode()

        if mapping_strategy == MappingStrategy.SEQUENCE_BASED:
            self.console.print("\n[cyan]Performing automatic sequence alignment...[/cyan]")

            # Separate protein residues from HETATM residues
            protein_residues, hetatm_residues = self._separate_protein_and_hetatm_residues(
                self.residue_mappings,
                self.structures[self.reference_idx][1]
            )

            if hetatm_residues:
                self.console.print(
                    f"\n[yellow]Note: {len(hetatm_residues)} non-protein residues excluded from sequence alignment[/yellow]"
                )

            if protein_residues:
                # Prepare target structures list
                _, ref_structure = self.structures[self.reference_idx]
                target_structures = [
                    (i, filename, struct)
                    for i, (filename, struct) in enumerate(self.structures)
                    if i != self.reference_idx
                ]

                # Create automatic mappings for protein residues
                protein_residues = self._create_automatic_residue_mappings(
                    protein_residues,
                    ref_structure,
                    target_structures
                )

                # Combine protein and HETATM residues
                self.residue_mappings = protein_residues + hetatm_residues

            # Offer to map HETATM residues interactively
            if hetatm_residues:
                self.console.print("\n[bold]HETATM Residue Mapping[/bold]")
                self.console.print(
                    "[grey50]Non-protein residues cannot be mapped by sequence alignment.\n"
                    "Interactive mapping uses spatial proximity after aligning protein residues.\n"
                    "For each ion/ligand/water, you'll select the corresponding residue in each target structure.\n\n"
                    "[cyan]Important:[/cyan] HETATMs in target structures will be transformed along with the protein,\n"
                    "even if not explicitly mapped. Reference HETATMs can be added to targets that lack them.[/grey50]"
                )
                map_hetatm = confirm_with_context(
                    processor=self.processor,
                    prompt="Interactively map non-protein residues using spatial proximity?",
                    default=False,
                    module="Structure Alignment - HETATM Mapping",
                    description="Interactive mapping of HETATM residues"
                )

                if map_hetatm:
                    # Perform preliminary protein alignment to get spatial positions
                    self.console.print("\n[cyan]Performing preliminary protein alignment for spatial reference...[/cyan]")

                    # Map HETATM residues interactively
                    hetatm_residues = self._map_hetatm_residues_interactive(
                        hetatm_residues,
                        protein_residues
                    )

                    # Update combined list
                    self.residue_mappings = protein_residues + hetatm_residues

        # For MANUAL strategy, residue_mappings already have the correct mappings (same numbering assumed)

        # Step 6: Output options
        self.console.print("\n[bold]Step 6: Output Options[/bold]")
        self.save_intermediates = confirm_with_context(
            processor=self.processor,
            prompt="Save intermediate alignment steps for debugging?",
            default=False,
            module="Structure Alignment",
            description="Save intermediate alignment steps"
        )

        # Perform alignment
        self.console.print("\n[bold]Performing Progressive Alignment...[/bold]")
        self._perform_progressive_alignment()

        # Add HETATM residues using final transformation matrices
        if self.hetatm_to_add:
            self.console.print("\n[bold]Adding HETATM residues to aligned structures...[/bold]")
            self._add_hetatm_residues_after_alignment()

        # Save results to workspace
        workspace = self._save_results_to_workspace(workspace)

        self.console.print("\n[bold green]✓ Structure alignment complete![/bold green]")

        return workspace

    def _align_on_redox_sites_interactive(self, workspace: Dict[str, Any] = None) -> Dict[str, Any]:
        """Align structures using redox site residues"""
        if workspace is None:
            workspace = self.processor._get_workspace()

        # Check for redox sites in workspace
        detected_redox_sites = workspace.get('detected_redox_sites', [])

        if not detected_redox_sites:
            self.console.print("[yellow]No redox sites detected in workspace.[/yellow]")
            self.console.print("[cyan]Please run redox site detection first (PDB Filter module).[/cyan]")
            return workspace

        self.console.print(Panel.fit(
            "[bold cyan]Redox Site-Based Alignment[/bold cyan]\n"
            "Align structures using redox-active site residues",
            border_style="cyan"
        ))

        # Display available redox sites
        self._display_redox_sites(detected_redox_sites)

        # Select which site(s) to use for alignment
        site_choices = self._select_redox_sites(detected_redox_sites)
        if site_choices is None:
            return workspace

        selected_sites = [detected_redox_sites[idx] for idx in site_choices]

        # Extract residues from all selected redox sites
        self.residue_mappings = self._extract_residues_from_sites(selected_sites)

        if len(site_choices) == 1:
            self.console.print(f"\n[green]Extracted {len(self.residue_mappings)} residue(s) from redox site {site_choices[0] + 1}[/green]")
        else:
            self.console.print(f"\n[green]Extracted {len(self.residue_mappings)} unique residue(s) from {len(site_choices)} redox sites[/green]")

        for mapping in self.residue_mappings:
            # Format residue with consistent coloring: chain:resname in white, resid in cyan
            ref = mapping.reference
            from rich.text import Text
            line = Text("  • ")
            line.append(f"{ref.chain_id}:{ref.resname}:")
            line.append(str(ref.resid), style="cyan")
            self.console.print(line)

        # Load structures
        self.console.print("\n[bold]Load PDB Structures for Alignment[/bold]")
        self.console.print("[grey50]Tip: When loading multi-chain structures, select all chains you want in the final output.[/grey50]")
        self.console.print("[grey50]     Alignment uses only matching residues, but transformation applies to all selected chains.[/grey50]")
        self.structures = self._load_structures_from_input(workspace)

        if len(self.structures) < 2:
            self.console.print("[red]Error: At least 2 structures required[/red]")
            return workspace

        # Select reference
        self.reference_idx = self._select_reference()

        # Ask about residue mapping strategy
        self.console.print("\n[bold]Residue Mapping Strategy[/bold]")
        mapping_strategy = self._select_mapping_strategy()

        # Handle structure-based alignment separately (uses TM-align directly)
        if mapping_strategy == MappingStrategy.STRUCTURE_BASED:
            self.console.print("\n[cyan]Note: Structure-based alignment will align the full structures using TM-align,[/cyan]")
            self.console.print("[cyan]ignoring the redox site residue selection for the alignment calculation.[/cyan]")

            # Structure-based alignment handles its own superposition
            self._perform_structure_based_alignment()

            # Save aligned structures
            self.console.print("\n[bold]Saving Aligned Structures...[/bold]")
            self._save_results_to_workspace(workspace)

            self.console.print("\n[bold green]✓ Structure-based alignment complete![/bold green]")
            return workspace

        if mapping_strategy == MappingStrategy.SEQUENCE_BASED:
            self.console.print("\n[cyan]Performing automatic sequence alignment...[/cyan]")

            # Separate protein residues from HETATM residues
            # Only protein residues can be mapped via sequence alignment
            protein_residues, hetatm_residues = self._separate_protein_and_hetatm_residues(
                self.residue_mappings,
                self.structures[self.reference_idx][1]
            )

            if hetatm_residues:
                self.console.print(
                    f"\n[yellow]Note: {len(hetatm_residues)} non-protein residues excluded from sequence alignment[/yellow]"
                )
                self.console.print("[grey50]Non-protein residues (ions, waters, ligands):[/grey50]")
                for mapping in hetatm_residues:
                    from rich.text import Text
                    ref = mapping.reference
                    line = Text("  • ", style="grey50")
                    line.append(f"{ref.chain_id}:{ref.resname}:", style="grey50")
                    line.append(str(ref.resid), style="cyan")
                    self.console.print(line)
                self.console.print(
                    "[cyan]These will only align if they exist at same position in target structures[/cyan]"
                )

            if not protein_residues:
                self.console.print("[yellow]Warning: No protein residues for sequence alignment![/yellow]")
                self.console.print("[yellow]Using original mappings (assumes same residue numbering)[/yellow]")
            else:
                # Prepare target structures list
                _, ref_structure = self.structures[self.reference_idx]
                target_structures = [
                    (i, filename, struct)
                    for i, (filename, struct) in enumerate(self.structures)
                    if i != self.reference_idx
                ]

                # Create automatic mappings for protein residues only
                protein_residues = self._create_automatic_residue_mappings(
                    protein_residues,
                    ref_structure,
                    target_structures
                )

                # Combine protein and HETATM residues
                self.residue_mappings = protein_residues + hetatm_residues

            # Show mapped residues
            self.console.print("\n[green]Protein Residue Mappings:[/green]")
            for mapping in protein_residues:
                from rich.text import Text
                ref = mapping.reference

                # Build the line using Text to avoid Rich markup interpretation
                line = Text("  ")
                line.append(f"{ref.chain_id}:{ref.resname}:")
                line.append(str(ref.resid), style="cyan")

                if mapping.mappings:
                    line.append(" ← ")
                    # Format target mappings
                    target_strs = []
                    for i, spec in mapping.mappings.items():
                        target_strs.append(f"{i}:{spec.chain_id}:{spec.resname}:{spec.resid}")
                    line.append(", ".join(target_strs))
                else:
                    line.append(" ")
                    line.append("(no targets mapped)", style="yellow")

                self.console.print(line)

            # Offer to map HETATM residues interactively
            if hetatm_residues:
                self.console.print("\n[bold]HETATM Residue Mapping[/bold]")
                self.console.print(
                    "[grey50]Non-protein residues cannot be mapped by sequence alignment.\n"
                    "Interactive mapping uses spatial proximity after aligning protein residues.\n"
                    "For each ion/water, you'll select the corresponding residue in each target structure.[/grey50]"
                )
                map_hetatm = confirm_with_context(
                    processor=self.processor,
                    prompt="Interactively map non-protein residues using spatial proximity?",
                    default=False,
                    module="Structure Alignment - HETATM Mapping",
                    description="Interactive mapping of HETATM residues"
                )

                if map_hetatm:
                    # Perform preliminary protein alignment to get spatial positions
                    self.console.print("\n[cyan]Performing preliminary protein alignment for spatial reference...[/cyan]")

                    # Map HETATM residues interactively
                    hetatm_residues = self._map_hetatm_residues_interactive(
                        hetatm_residues,
                        protein_residues
                    )

                    # Update combined list
                    self.residue_mappings = protein_residues + hetatm_residues

                    # Show HETATM mappings
                    self.console.print("\n[green]HETATM Residue Mappings:[/green]")
                    for mapping in hetatm_residues:
                        ref_str = str(mapping.reference)
                        if mapping.mappings:
                            target_strs = [f"{i}:{spec}" for i, spec in mapping.mappings.items()]
                            self.console.print(f"  {ref_str} ← {', '.join(target_strs)}")
                        else:
                            self.console.print(f"  {ref_str} [grey50](skipped)[/grey50]")

        # Select alignment mode (suggest backbone for redox sites)
        self.console.print("\n[bold]Select Alignment Mode[/bold]")
        self.console.print("[cyan]Recommendation: Use backbone mode for redox site alignment[/cyan]")
        self.alignment_mode = self._select_alignment_mode(default="backbone")

        # Output options
        self.save_intermediates = confirm_with_context(
            processor=self.processor,
            prompt="Save intermediate alignment steps?",
            default=False,
            module="Structure Alignment - Redox Sites",
            description="Save intermediate alignment steps"
        )

        # Perform alignment
        self.console.print("\n[bold]Performing Redox Site-Based Alignment...[/bold]")
        self._perform_progressive_alignment()

        # Add HETATM residues using final transformation matrices
        if self.hetatm_to_add:
            self.console.print("\n[bold]Adding HETATM residues to aligned structures...[/bold]")
            self._add_hetatm_residues_after_alignment()

        # Save results
        workspace = self._save_results_to_workspace(workspace)

        self.console.print("\n[bold green]✓ Redox site-based alignment complete![/bold green]")

        return workspace

    def _view_alignment_results(self, workspace: Dict[str, Any] = None) -> Dict[str, Any]:
        """Display alignment results"""
        if workspace is None:
            workspace = self.processor._get_workspace()

        if not self.alignment_results:
            self.console.print("[yellow]No alignment results available.[/yellow]")
            self.console.print("[cyan]Please run structure alignment first.[/cyan]")
            return workspace

        # Display results table
        self._display_alignment_results()

        return workspace

    def _select_workspace_structures(self, workspace: Dict[str, Any]) -> List[str]:
        """Select workspace structures (single or multiple) using structure selector"""
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(workspace, self.console, self.processor)
        available = selector.get_available_structures()

        if not available:
            return []

        # Use the multi-selection method from structure selector
        return selector._interactive_multi_selection(available)

    def _get_aligned_structure_files(self, workspace: Dict[str, Any]) -> List[str]:
        """Get list of aligned structure files from previous alignment run

        Returns:
            List of file paths to aligned structures that exist on disk
        """
        aligned_files = []

        # Check for aligned_pdb_file (single reference structure from last run)
        aligned_pdb = workspace.get('aligned_pdb_file')
        if aligned_pdb and os.path.exists(aligned_pdb):
            aligned_files.append(aligned_pdb)

        # Also check for other *_aligned.pdb files in current directory
        # These would be the target structures from the previous alignment
        current_dir = Path.cwd()
        for aligned_file in current_dir.glob("*_aligned.pdb"):
            aligned_file_str = str(aligned_file)
            # Don't duplicate if already added
            if aligned_file_str not in aligned_files and os.path.exists(aligned_file_str):
                aligned_files.append(aligned_file_str)

        return aligned_files

    def _load_structures_from_input(self, workspace: Dict[str, Any] = None) -> List[Tuple[str, Structure]]:
        """Load structures from user input with enhanced options"""
        pdb_files = []

        # Check if workspace has aligned structures from previous alignment
        if workspace:
            aligned_files = self._get_aligned_structure_files(workspace)
            if aligned_files:
                self.console.print(f"\n[green]Found {len(aligned_files)} aligned structure(s) from previous alignment:[/green]")
                for aligned_file in aligned_files:
                    self.console.print(f"[green]  • {aligned_file}[/green]")

                use_aligned = confirm_with_context(
                    processor=self.processor,
                    prompt="Include these aligned structures for further alignment?",
                    default=True,
                    module="Structure Alignment",
                    description="Include aligned structures from previous run"
                )

                if use_aligned:
                    pdb_files.extend(aligned_files)
                    self.console.print(f"[green]✓ Added {len(aligned_files)} aligned structure(s)[/green]")
            else:
                # No aligned structures, check for other workspace structures
                workspace_pdbs = self._select_workspace_structures(workspace)
                if workspace_pdbs:
                    # User already selected these, so add them directly
                    pdb_files.extend(workspace_pdbs)
                    self.console.print(f"[green]✓ Added {len(workspace_pdbs)} workspace structure(s)[/green]")

        # Interactive structure loading loop
        self.console.print("\n[bold]Add Structures for Alignment[/bold]")

        while True:
            self.console.print(f"\n[cyan]Current structures: {len(pdb_files)}[/cyan]")

            # Show available local structure files (PDB and CIF)
            current_dir = Path.cwd()
            structure_files_in_dir = []
            for pattern in ["*.pdb", "*.cif", "*.mmcif"]:
                structure_files_in_dir.extend([f for f in current_dir.glob(pattern) if str(f) not in pdb_files])
            structure_files_in_dir = sorted(set(structure_files_in_dir))

            if structure_files_in_dir:
                self.console.print("\n[grey50]Available local structure files:[/grey50]")
                for structure_file in structure_files_in_dir:
                    self.console.print(f"[grey50]  • {structure_file.name}[/grey50]")

            # Menu options
            self.console.print("\n[bold]Options:[/bold]")
            self.console.print("  1. Browse for local structure file (PDB/CIF)")
            self.console.print("  2. Download PDB from database")
            self.console.print("  3. Done (proceed with alignment)")

            choice = prompt_with_context(
                processor=self.processor,
                prompt="Select option",
                choices=["1", "2", "3"],
                default="3" if len(pdb_files) >= 2 else "1",
                module="Structure Alignment",
                description="Add structures",
                options_map={"1": "Browse local file", "2": "Download PDB", "3": "Done"}
            )

            if choice == "3":
                break
            elif choice == "1":
                # Add local files using file browser
                from proprep.structure_prep.pdb_loader import display_pdb_file_menu
                file_path = display_pdb_file_menu(
                    directory=".",
                    console=self.console,
                    processor=self.processor
                )
                if file_path and file_path not in pdb_files:
                    pdb_files.append(file_path)
                    self.console.print(f"[green]✓ Added: {file_path}[/green]")
                elif file_path:
                    self.console.print(f"[yellow]Already added: {file_path}[/yellow]")

            elif choice == "2":
                # Download from PDB
                pdb_id = prompt_with_context(
                    processor=self.processor,
                    prompt="Enter PDB ID to download (e.g., 1CLL)",
                    default="",
                    module="Structure Alignment",
                    description="Download PDB from database"
                )
                if pdb_id.strip():
                    downloaded_file = self._download_pdb(pdb_id.strip())
                    if downloaded_file:
                        pdb_files.append(downloaded_file)

        if len(pdb_files) < 2:
            self.console.print("[yellow]At least 2 structures required for alignment[/yellow]")
            return []

        # Load structures with model/chain selection
        return self._load_and_select_structures(pdb_files)

    def _get_structure_parser(self, file_path: str):
        """Get the appropriate parser based on file extension"""
        file_ext = Path(file_path).suffix.lower()
        if file_ext in ['.cif', '.mmcif']:
            return MMCIFParser(QUIET=True)
        else:  # .pdb or any other extension
            return PDBParser(QUIET=True)

    def _load_and_verify_structures(self, pdb_files: List[str]) -> List[Tuple[str, Structure]]:
        """Load PDB/CIF structures and verify atom counts"""
        structures = []

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=self.console
        ) as progress:
            for pdb_file in pdb_files:
                task = progress.add_task(f"Loading {pdb_file}...", total=None)

                if not os.path.exists(pdb_file):
                    self.console.print(f"[red]✗ {pdb_file} not found[/red]")
                    continue

                # Use appropriate parser based on file extension
                parser = self._get_structure_parser(pdb_file)
                structure = parser.get_structure(Path(pdb_file).stem, pdb_file)
                is_valid, file_count, parsed_count, alt_count, explanation = self._verify_atom_count(pdb_file, structure)

                if not is_valid:
                    if explanation:
                        self.console.print(
                            f"[yellow]ℹ {pdb_file}: {parsed_count} atoms parsed from {file_count} lines "
                            f"(difference due to {explanation})[/yellow]"
                        )
                    else:
                        self.console.print(
                            f"[red]✗ {pdb_file}: BioPython parsed {parsed_count} atoms "
                            f"but file contains {file_count} atoms[/red]"
                        )
                        self.console.print("[yellow]  Some atoms may not be recognized. Aborting.[/yellow]")
                        continue

                structures.append((pdb_file, structure))
                self.console.print(f"[green]✓ {pdb_file} ({parsed_count} atoms)[/green]")
                progress.remove_task(task)

        return structures

    def _verify_atom_count(self, pdb_file: str, structure: Structure) -> Tuple[bool, int, int]:
        """Verify BioPython parsed all atoms from structure file

        Returns:
            Tuple of (is_valid, file_atom_count, parsed_atom_count, alternate_location_count, explanation)
        """
        parsed_atom_count = sum(1 for _ in structure.get_atoms())

        # Check file format - CIF files have different structure than PDB
        file_ext = Path(pdb_file).suffix.lower()
        if file_ext in ['.cif', '.mmcif']:
            # For CIF files, we can't easily count atoms from the file
            # Just trust BioPython's parsing
            return True, parsed_atom_count, parsed_atom_count, 0, None

        # For PDB files, verify atom count from file
        file_atom_count = 0
        alternate_locations = {}  # Track atoms with alternate locations

        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith(('ATOM  ', 'HETATM')):
                    file_atom_count += 1
                    # Check for alternate location indicator (column 17, index 16)
                    if len(line) > 16:
                        altloc = line[16]
                        if altloc not in (' ', ''):
                            # Get atom identifier (serial number from columns 7-11)
                            atom_serial = line[6:11].strip()
                            atom_name = line[12:16].strip()
                            residue = line[17:26].strip()
                            atom_key = (residue, atom_name)

                            if atom_key not in alternate_locations:
                                alternate_locations[atom_key] = []
                            alternate_locations[atom_key].append(altloc)

        # Count how many extra lines are due to alternate locations
        alt_count = sum(len(locs) - 1 for locs in alternate_locations.values() if len(locs) > 1)

        # Check if mismatch is explained by alternate locations
        is_valid = file_atom_count == parsed_atom_count
        expected_difference = alt_count
        actual_difference = file_atom_count - parsed_atom_count

        explanation = None
        if not is_valid and actual_difference == expected_difference and alt_count > 0:
            # Mismatch is fully explained by alternate locations
            num_atoms_with_alts = sum(1 for locs in alternate_locations.values() if len(locs) > 1)
            explanation = f"alternate locations ({num_atoms_with_alts} atoms with multiple conformations)"
            is_valid = True  # This is actually valid - BioPython is handling it correctly

        return is_valid, file_atom_count, parsed_atom_count, alt_count, explanation

    def _select_reference(self) -> int:
        """Prompt user to select reference structure"""
        self.console.print("\n[grey50]Choose your reference structure carefully:[/grey50]")
        self.console.print("[grey50]  • Experimental structures (RCSB): Use if they contain HETATMs (ligands, ions, cofactors)[/grey50]")
        self.console.print("[grey50]    you want to include in alignment or transfer to other structures[/grey50]")
        self.console.print("[grey50]  • Predicted structures (AlphaFold): Use if you want complete sequence coverage[/grey50]")
        self.console.print("[grey50]    as the reference (no missing residues)[/grey50]\n")

        for i, (filename, _) in enumerate(self.structures):
            self.console.print(f"  [{i + 1}] {filename}")

        # Build options map (displayed 1-based to the user)
        options_map = {str(i + 1): filename for i, (filename, _) in enumerate(self.structures)}

        ref_idx = int_prompt_with_context(
            processor=self.processor,
            prompt="Enter reference structure index",
            default=1,
            show_default=True,
            module="Structure Alignment",
            description="Select reference structure",
            options_map=options_map
        )

        if ref_idx < 1 or ref_idx > len(self.structures):
            self.console.print("[red]Error: Invalid reference index[/red]")
            return 0

        # Convert back to 0-based for use as a list subscript
        ref_idx -= 1

        self.console.print(f"[green]Reference: {self.structures[ref_idx][0]}[/green]")
        self.console.print(
            "\n[grey50]Note: The transformation matrix will be calculated using the selected alignment residues,\n"
            "but will be applied to ALL atoms in each target structure (including HETATMs, waters, etc.).[/grey50]"
        )
        # Swap the viewer to the chosen reference so the user has visual
        # confirmation of which structure the alignment will be based on.
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer
            _viewer.show_structure(self.structures[ref_idx][0])
        except Exception:
            pass
        return ref_idx

    def _select_alignment_mode(self, default="backbone") -> str:
        """Prompt user to select alignment mode"""
        self.console.print("  [1] Full residue (all atoms)")
        self.console.print("  [2] Backbone only (N, CA, C, O)")
        self.console.print("  [3] Alpha carbon only (CA)")

        default_choice = 2 if default == "backbone" else 3 if default == "alpha_carbon" else 1

        mode_choice = int_prompt_with_context(
            processor=self.processor,
            prompt="Enter mode",
            default=default_choice,
            show_default=True,
            module="Structure Alignment",
            description="Select alignment mode",
            options_map={
                "1": "Full residue (all atoms)",
                "2": "Backbone only (N, CA, C, O)",
                "3": "Alpha carbon only (CA)"
            }
        )

        mode_map = {
            1: AlignmentMode.FULL,
            2: AlignmentMode.BACKBONE,
            3: AlignmentMode.ALPHA_CARBON
        }

        if mode_choice not in mode_map:
            self.console.print("[yellow]Invalid mode, using backbone[/yellow]")
            return AlignmentMode.BACKBONE

        mode = mode_map[mode_choice]
        self.console.print(f"[green]Mode: {mode}[/green]")
        return mode

    def _select_mapping_strategy(self) -> str:
        """Prompt user to select residue mapping strategy"""
        self.console.print("  [1] Sequence-based (automatic) - for homologous proteins with similar sequences")

        if TMTOOLS_AVAILABLE:
            self.console.print("  [2] Structure-based (TM-align) - for proteins with same fold but different sequences")
        else:
            self.console.print("  [grey50][2] Structure-based (TM-align) - NOT AVAILABLE (install: pip install tmtools)[/grey50]")

        self.console.print("  [3] Manual specification - for identical numbering or custom residue selection")

        strategy_choice = int_prompt_with_context(
            processor=self.processor,
            prompt="Select method",
            default=1,
            show_default=True,
            module="Structure Alignment",
            description="Select residue mapping strategy",
            options_map={
                "1": "Sequence-based (automatic)",
                "2": "Structure-based (TM-align)",
                "3": "Manual specification"
            }
        )

        strategy_map = {
            1: MappingStrategy.SEQUENCE_BASED,
            2: MappingStrategy.STRUCTURE_BASED,
            3: MappingStrategy.MANUAL
        }

        if strategy_choice not in strategy_map:
            self.console.print("[yellow]Invalid choice, using sequence-based[/yellow]")
            return MappingStrategy.SEQUENCE_BASED

        strategy = strategy_map[strategy_choice]

        # Check if TM-align was selected but not available
        if strategy == MappingStrategy.STRUCTURE_BASED and not TMTOOLS_AVAILABLE:
            self.console.print("\n[red]TM-align is not available.[/red]")
            self.console.print("[yellow]To enable structure-based alignment, install tmtools:[/yellow]")
            self.console.print("[cyan]  pip install tmtools[/cyan]")
            self.console.print("\n[yellow]Falling back to sequence-based alignment.[/yellow]")
            return MappingStrategy.SEQUENCE_BASED

        self.console.print(f"[green]Strategy: {strategy}[/green]")
        return strategy

    def _get_ca_coords_and_seq(self, structure: Structure, chain_ids: Optional[List[str]] = None) -> Tuple[np.ndarray, str, List[Tuple[str, int]]]:
        """Extract CA coordinates, sequence, and residue info from structure.

        Args:
            structure: BioPython Structure object
            chain_ids: Optional list of chain IDs to include. If None, includes all chains.

        Returns:
            Tuple of (coordinates array, sequence string, list of (chain_id, resid) tuples)
        """
        aa_map = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
        }

        coords = []
        seq = []
        residue_info = []

        for model in structure:
            for chain in model:
                if chain_ids is not None and chain.id not in chain_ids:
                    continue
                for residue in chain:
                    hetero_flag, resseq, icode = residue.id
                    if hetero_flag != ' ':  # Skip HETATM
                        continue
                    if residue.resname.strip() not in aa_map:
                        continue
                    if 'CA' not in residue:
                        continue

                    coords.append(residue['CA'].get_coord())
                    seq.append(aa_map.get(residue.resname.strip(), 'X'))
                    residue_info.append((chain.id, resseq))

        return np.array(coords), ''.join(seq), residue_info

    def _select_tmalign_chain_mode(self, ref_chains: List[str], target_chains: List[str]) -> Tuple[str, Optional[str]]:
        """Prompt user to select how to handle multi-chain TM-align alignment.

        Returns:
            Tuple of (mode, selected_chain) where mode is 'combined' or 'single'
        """
        self.console.print("\n[bold]Structure-based Alignment Options[/bold]")
        self.console.print(f"  Reference chains: {', '.join(ref_chains)}")
        self.console.print(f"  Target chains: {', '.join(target_chains)}")
        self.console.print("")
        self.console.print("  [1] Align all chains combined (single best-fit transformation)")
        self.console.print("  [2] Align on specific chain (use one chain's transformation)")

        choice = int_prompt_with_context(
            processor=self.processor,
            prompt="Select option",
            default=1,
            show_default=True,
            module="Structure Alignment - TM-align",
            description="Select chain alignment mode",
            options_map={
                "1": "Align all chains combined",
                "2": "Align on specific chain"
            }
        )

        if choice == 2:
            # Let user select which chain to align on
            common_chains = sorted(set(ref_chains) & set(target_chains))
            if not common_chains:
                self.console.print("[yellow]No common chains found, using combined alignment[/yellow]")
                return 'combined', None

            self.console.print("\n[bold]Select chain for alignment:[/bold]")
            for i, chain_id in enumerate(common_chains, 1):
                self.console.print(f"  [{i}] Chain {chain_id}")

            chain_choice = int_prompt_with_context(
                processor=self.processor,
                prompt="Select chain",
                default=1,
                show_default=True,
                module="Structure Alignment - TM-align",
                description="Select chain for alignment",
                options_map={str(i): f"Chain {c}" for i, c in enumerate(common_chains, 1)}
            )

            if 1 <= chain_choice <= len(common_chains):
                return 'single', common_chains[chain_choice - 1]

        return 'combined', None

    def _perform_structure_based_alignment(self) -> Dict[str, Any]:
        """Perform alignment using TM-align for structure-based superposition.

        Returns:
            Dictionary with alignment results including RMSD and TM-scores
        """
        import copy

        self.console.print("\n[cyan]Performing TM-align structure-based alignment...[/cyan]")

        _, ref_structure = self.structures[self.reference_idx]

        # Get chain IDs from reference
        ref_chains = []
        for model in ref_structure:
            for chain in model:
                # Check if chain has protein residues
                has_protein = any(
                    res.id[0] == ' ' and res.resname.strip() in {
                        'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE',
                        'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER',
                        'THR', 'VAL', 'TRP', 'TYR'
                    }
                    for res in chain
                )
                if has_protein:
                    ref_chains.append(chain.id)

        self.aligned_structures = {}
        self.aligned_structures[self.reference_idx] = ref_structure
        self.alignment_results = []

        # Process each target structure
        for i, (target_file, target_structure) in enumerate(self.structures):
            if i == self.reference_idx:
                continue

            self.console.print(f"\n[bold]Aligning: {Path(target_file).name}[/bold]")

            # Get target chains
            target_chains = []
            for model in target_structure:
                for chain in model:
                    has_protein = any(
                        res.id[0] == ' ' and res.resname.strip() in {
                            'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE',
                            'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER',
                            'THR', 'VAL', 'TRP', 'TYR'
                        }
                        for res in chain
                    )
                    if has_protein:
                        target_chains.append(chain.id)

            # Determine alignment mode for multi-chain structures
            if len(ref_chains) > 1 or len(target_chains) > 1:
                chain_mode, selected_chain = self._select_tmalign_chain_mode(ref_chains, target_chains)
            else:
                chain_mode = 'combined'
                selected_chain = None

            # Extract coordinates
            if chain_mode == 'single' and selected_chain:
                ref_coords, ref_seq, ref_info = self._get_ca_coords_and_seq(ref_structure, [selected_chain])
                target_coords, target_seq, target_info = self._get_ca_coords_and_seq(target_structure, [selected_chain])
                self.console.print(f"  Aligning on chain {selected_chain}")
            else:
                ref_coords, ref_seq, ref_info = self._get_ca_coords_and_seq(ref_structure)
                target_coords, target_seq, target_info = self._get_ca_coords_and_seq(target_structure)
                self.console.print(f"  Aligning all chains combined")

            self.console.print(f"  Reference: {len(ref_coords)} CA atoms")
            self.console.print(f"  Target: {len(target_coords)} CA atoms")

            if len(ref_coords) == 0 or len(target_coords) == 0:
                self.console.print("[red]Error: No CA atoms found for alignment[/red]")
                continue

            # Run TM-align
            result = tmtools.tm_align(target_coords, ref_coords, target_seq, ref_seq)

            # Report results
            self.console.print(f"\n  [green]TM-score (normalized by target): {result.tm_norm_chain1:.4f}[/green]")
            self.console.print(f"  [green]TM-score (normalized by reference): {result.tm_norm_chain2:.4f}[/green]")
            self.console.print(f"  [green]RMSD: {result.rmsd:.3f} Å[/green]")

            # Apply transformation to target structure
            aligned_structure = copy.deepcopy(target_structure)
            for atom in aligned_structure.get_atoms():
                coord = atom.get_coord()
                new_coord = np.dot(coord, result.u.T) + result.t
                atom.set_coord(new_coord)

            self.aligned_structures[i] = aligned_structure

            # Store results
            self.alignment_results.append({
                'structure': target_file,
                'tm_score_target': result.tm_norm_chain1,
                'tm_score_ref': result.tm_norm_chain2,
                'rmsd': result.rmsd,
                'aligned_length': len(result.seqM.replace(' ', '').replace('-', '')),
                'chain_mode': chain_mode,
                'selected_chain': selected_chain
            })

        return {
            'method': 'structure_based',
            'results': self.alignment_results
        }

    def _extract_all_residues_from_reference(self) -> List[ResidueMapping]:
        """Extract all residues from reference structure including HETATM"""
        _, ref_structure = self.structures[self.reference_idx]

        mappings = []
        for model in ref_structure:
            for chain in model:
                for residue in chain:
                    hetero_flag, resseq, icode = residue.id

                    # Skip water molecules
                    if residue.resname.strip() in ['HOH', 'WAT', 'H2O']:
                        continue

                    # Create residue specification (includes protein and HETATM)
                    resname = residue.resname.strip()
                    chain_id = chain.id
                    resid = resseq

                    res_spec = ResidueSpec(chain_id, resname, resid)
                    mappings.append(ResidueMapping(res_spec, {}))

        return mappings

    def _get_residue_subset_interactive(self) -> List[ResidueMapping]:
        """Get subset of residues via range or manual entry"""
        self.console.print("\n[bold]Residue Selection Options:[/bold]")
        self.console.print("  [1] Specify residue range (e.g., A:50-150)")
        self.console.print("  [2] Manual residue-by-residue entry")

        choice = prompt_with_context(
            processor=self.processor,
            prompt="Select option",
            choices=["1", "2"],
            default="1",
            module="Structure Alignment",
            description="Select residue subset method",
            options_map={"1": "Residue range", "2": "Manual entry"}
        )

        if choice == "1":
            return self._get_residue_range_interactive()
        else:
            return self._get_residue_mappings_interactive()

    def _get_residue_range_interactive(self) -> List[ResidueMapping]:
        """Get residue range from user (e.g., A:50-150)"""
        self.console.print("\n[grey50]Format: chain:start-end (e.g., A:50-150)[/grey50]")
        self.console.print("[grey50]You can specify multiple ranges separated by commas[/grey50]")
        self.console.print("[grey50]Example: A:50-150,B:20-80[/grey50]\n")

        range_input = prompt_with_context(
            processor=self.processor,
            prompt="Enter residue range(s)",
            default="",
            module="Structure Alignment",
            description="Enter residue range specification",
            options_map={"custom": "Range spec (e.g., A:50-150)"}
        )

        if not range_input:
            return []

        # Parse range specifications
        _, ref_structure = self.structures[self.reference_idx]
        mappings = []

        for range_spec in range_input.split(','):
            range_spec = range_spec.strip()
            try:
                # Parse format: chain:start-end
                if ':' not in range_spec or '-' not in range_spec:
                    raise ValueError("Invalid format. Use chain:start-end (e.g., A:50-150)")

                chain_part, range_part = range_spec.split(':', 1)
                chain_id = chain_part.strip()
                start_str, end_str = range_part.split('-', 1)
                start_resid = int(start_str.strip())
                end_resid = int(end_str.strip())

                # Extract residues in this range from reference structure
                for model in ref_structure:
                    for chain in model:
                        if chain.id != chain_id:
                            continue

                        for residue in chain:
                            hetero_flag, resseq, icode = residue.id
                            # Skip non-protein residues
                            if hetero_flag != ' ' and hetero_flag != 'H_MSE':
                                continue

                            # Check if in range
                            if start_resid <= resseq <= end_resid:
                                resname = residue.resname.strip()
                                res_spec = ResidueSpec(chain_id, resname, resseq)
                                mappings.append(ResidueMapping(res_spec, {}))

                self.console.print(f"[green]✓ Added {sum(1 for m in mappings if m.reference.chain_id == chain_id)} residues from {range_spec}[/green]")

            except Exception as e:
                self.console.print(f"[red]✗ Error parsing '{range_spec}': {e}[/red]")

        return mappings

    def _get_residue_mappings_interactive(self) -> List[ResidueMapping]:
        """Get residue mappings from user input"""
        self.console.print("[grey50]Format: chain:resname:resid[/grey50]")
        self.console.print("[grey50]   OR:  chain:resname:resid <- idx:chain:resname:resid, ...[/grey50]")
        self.console.print("[grey50]Enter empty line to finish[/grey50]\n")

        mappings = []
        while True:
            spec = prompt_with_context(
                processor=self.processor,
                prompt=f"Residue {len(mappings) + 1}",
                default="",
                module="Structure Alignment",
                description=f"Enter residue specification {len(mappings) + 1}",
                options_map={"custom": "Residue spec (e.g., A:GLY:123)"}
            )
            if not spec:
                break

            try:
                mapping = self._parse_residue_mapping(spec)
                mappings.append(mapping)
                self.console.print(f"[green]✓ Added: {mapping.reference}[/green]")
            except Exception as e:
                self.console.print(f"[red]✗ Error: {e}[/red]")

        return mappings

    def _parse_residue_mapping(self, mapping_str: str) -> ResidueMapping:
        """Parse residue mapping from string"""
        if '<-' in mapping_str:
            ref_part, mappings_part = mapping_str.split('<-', 1)
            ref_part = ref_part.strip()

            # Parse reference
            struct_idx, ref_spec = self._parse_residue_spec(ref_part)
            if struct_idx is not None:
                raise ValueError("Reference residue should not have structure index")

            # Parse mappings
            mappings = {}
            for mapping in mappings_part.split(','):
                mapping = mapping.strip()
                if mapping:
                    struct_idx, spec = self._parse_residue_spec(mapping)
                    if struct_idx is None:
                        raise ValueError(f"Explicit mapping must include structure index: {mapping}")
                    mappings[struct_idx] = spec

            return ResidueMapping(ref_spec, mappings)
        else:
            # Simple format - applies to all structures
            struct_idx, ref_spec = self._parse_residue_spec(mapping_str)
            if struct_idx is not None:
                raise ValueError(f"Simple format should not have structure index: {mapping_str}")
            return ResidueMapping(ref_spec, {})

    def _parse_residue_spec(self, spec: str) -> Tuple[Optional[int], ResidueSpec]:
        """Parse residue specification"""
        parts = spec.strip().split(':')

        if len(parts) == 3:
            # Simple format: chain:resname:resid
            chain_id, resname, resid = parts
            return None, ResidueSpec(chain_id, resname.upper(), int(resid))
        elif len(parts) == 4:
            # Explicit format: struct_idx:chain:resname:resid
            struct_idx, chain_id, resname, resid = parts
            return int(struct_idx), ResidueSpec(chain_id, resname.upper(), int(resid))
        else:
            raise ValueError(f"Invalid residue specification: {spec}")

    def _perform_progressive_alignment(self):
        """Perform progressive superimposition with cumulative residues"""
        self.aligned_structures = {}
        _, ref_structure = self.structures[self.reference_idx]
        self.aligned_structures[self.reference_idx] = ref_structure
        self.alignment_results = []

        # Track RMSD values for plotting
        rmsd_by_step = {}  # step -> {structure_idx: rmsd}

        # Progressive alignment: add one residue at a time with progress bar
        from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn, TimeRemainingColumn

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            TextColumn("RMSD: {task.fields[rmsd]:.3f} Å"),
            TimeRemainingColumn(),
            console=self.console
        ) as progress:
            task = progress.add_task(
                "[cyan]Aligning structures...",
                total=len(self.residue_mappings),
                rmsd=0.0
            )

            for step in range(1, len(self.residue_mappings) + 1):
                rmsd_by_step[step] = {}

                # Build structure list for alignment
                if step > 1:
                    structures_for_alignment = []
                    for i, (filename, orig_structure) in enumerate(self.structures):
                        if i in self.aligned_structures:
                            structures_for_alignment.append((filename, self.aligned_structures[i]))
                        else:
                            structures_for_alignment.append((filename, orig_structure))
                else:
                    structures_for_alignment = self.structures

                # Collect atoms for this step
                alignment_data = self._collect_alignment_atoms(
                    structures_for_alignment,
                    self.reference_idx,
                    self.residue_mappings,
                    self.alignment_mode,
                    step
                )

                if not alignment_data:
                    progress.console.print("[red]No valid alignment data for this step[/red]")
                    break

                # Track average RMSD for this step
                step_rmsds = []

                # Align each target structure
                for i, (target_file, target_structure) in enumerate(self.structures):
                    if i == self.reference_idx:
                        continue

                    # Get structure to align
                    if step > 1 and i in self.aligned_structures:
                        structure_to_align = self.aligned_structures[i]
                    else:
                        structure_to_align = target_structure

                    # Find alignment data for this structure
                    target_idx = i if i < self.reference_idx else i - 1
                    if target_idx >= len(alignment_data):
                        continue

                    ref_atoms, target_atoms = alignment_data[target_idx]

                    # Perform superimposition
                    super_imposer = Superimposer()
                    super_imposer.set_atoms(ref_atoms, target_atoms)

                    # Apply transformation
                    atoms_to_transform = list(structure_to_align.get_atoms())
                    super_imposer.apply(atoms_to_transform)

                    rmsd = super_imposer.rms
                    self.aligned_structures[i] = structure_to_align
                    step_rmsds.append(rmsd)

                    # Store RMSD for plotting
                    rmsd_by_step[step][i] = rmsd

                    # Store final transformation matrix (will be overwritten each step until the last)
                    if step == len(self.residue_mappings):
                        self.final_transformation_matrices[i] = super_imposer

                    # Record result
                    residue_desc = f"1-{step}" if step > 1 else "1"
                    self.alignment_results.append({
                        'step': step,
                        'structure': target_file,
                        'residues': residue_desc,
                        'rmsd': rmsd
                    })

                    # Save intermediate if requested
                    if self.save_intermediates and step < len(self.residue_mappings):
                        base_name = Path(target_file).stem + f"_step{step}"
                        intermediate_file = self._get_indexed_filename(base_name)
                        self._save_structure(structure_to_align, intermediate_file)

                # Update progress bar with average RMSD
                avg_rmsd = sum(step_rmsds) / len(step_rmsds) if step_rmsds else 0.0
                progress.update(task, advance=1, rmsd=avg_rmsd)

        # After alignment, display ASCII plot of RMSD vs step
        self.console.print("\n")
        self._plot_rmsd_progress(rmsd_by_step)

    def _plot_rmsd_progress(self, rmsd_by_step):
        """Create ASCII plot of RMSD vs alignment step"""
        if not rmsd_by_step:
            return

        # Get all structure indices and steps
        structure_indices = set()
        for step_data in rmsd_by_step.values():
            structure_indices.update(step_data.keys())

        if not structure_indices:
            return

        # Plot for each structure
        for struct_idx in sorted(structure_indices):
            filename = Path(self.structures[struct_idx][0]).stem

            # Collect RMSD values for this structure
            steps = []
            rmsds = []
            for step in sorted(rmsd_by_step.keys()):
                if struct_idx in rmsd_by_step[step]:
                    steps.append(step)
                    rmsds.append(rmsd_by_step[step][struct_idx])

            if not rmsds:
                continue

            # Create simple ASCII plot
            self.console.print(f"\n[bold cyan]RMSD Progress: {filename}[/bold cyan]")
            self.console.print("[grey50]  ● = data point, · = connecting line[/grey50]")
            self._plot_ascii_line(steps, rmsds)

    def _plot_ascii_line(self, x_values, y_values, width=60, height=15):
        """Create a simple ASCII line plot"""
        if not x_values or not y_values:
            return

        min_x, max_x = min(x_values), max(x_values)
        min_y, max_y = min(y_values), max(y_values)

        # Add padding to y-axis
        y_range = max_y - min_y
        if y_range == 0:
            y_range = 1
        min_y = max(0, min_y - y_range * 0.1)
        max_y = max_y + y_range * 0.1

        # Create grid
        grid = [[' ' for _ in range(width)] for _ in range(height)]

        # Plot points and connect with lines
        for i in range(len(x_values)):
            x_pos = int((x_values[i] - min_x) / (max_x - min_x) * (width - 1)) if max_x > min_x else 0
            y_pos = height - 1 - int((y_values[i] - min_y) / (max_y - min_y) * (height - 1))

            if 0 <= x_pos < width and 0 <= y_pos < height:
                grid[y_pos][x_pos] = '●'

                # Draw line to next point
                if i < len(x_values) - 1:
                    next_x = int((x_values[i+1] - min_x) / (max_x - min_x) * (width - 1)) if max_x > min_x else 0
                    next_y = height - 1 - int((y_values[i+1] - min_y) / (max_y - min_y) * (height - 1))

                    # Simple line drawing
                    steps_x = abs(next_x - x_pos)
                    steps_y = abs(next_y - y_pos)
                    steps = max(steps_x, steps_y)

                    for s in range(1, steps):
                        interp_x = x_pos + int((next_x - x_pos) * s / steps)
                        interp_y = y_pos + int((next_y - y_pos) * s / steps)
                        if 0 <= interp_x < width and 0 <= interp_y < height:
                            if grid[interp_y][interp_x] == ' ':
                                grid[interp_y][interp_x] = '·'

        # Print grid with axes
        self.console.print(f"  RMSD (Å)")
        for i, row in enumerate(grid):
            y_val = max_y - (i / (height - 1)) * (max_y - min_y)
            self.console.print(f"  {y_val:5.2f} │{''.join(row)}")

        # X-axis
        self.console.print(f"        └{'─' * width}")
        self.console.print(f"         {min_x:<{width//2}}{max_x:>{width//2}}")
        self.console.print(f"         {'Residues':^{width}}")

        # Final RMSD
        final_rmsd = y_values[-1]
        self.console.print(f"\n  [green]Final RMSD: {final_rmsd:.3f} Å[/green]")

    def _add_hetatm_residues_after_alignment(self):
        """Add HETATM residues to aligned structures using final transformation matrices"""
        from copy import deepcopy

        ref_structure = self.structures[self.reference_idx][1]

        for ref_residue, ref_spec, target_indices in self.hetatm_to_add:
            # Use Text to avoid Rich markup interpretation of colons
            from rich.text import Text
            message = Text("\nAdding ")
            message.append(f"{ref_spec.chain_id}:{ref_spec.resname}:{ref_spec.resid}")
            message.append(" to target structures...")
            self.console.print(message)

            for target_idx in target_indices:
                # Get the final transformation matrix for this target
                if target_idx not in self.final_transformation_matrices:
                    self.console.print(f"  [yellow]Warning: No transformation matrix for structure {target_idx}, skipping[/yellow]")
                    continue

                superimposer = self.final_transformation_matrices[target_idx]

                # Get the aligned target structure
                if target_idx not in self.aligned_structures:
                    self.console.print(f"  [yellow]Warning: No aligned structure for structure {target_idx}, skipping[/yellow]")
                    continue

                target_structure = self.aligned_structures[target_idx]
                target_file = self.structures[target_idx][0]

                # Add the transformed residue
                result = self._add_transformed_residue(
                    target_structure,
                    ref_residue,
                    ref_spec.resname,
                    superimposer
                )

                if result:
                    target_chain_id, new_resid = result
                    self.console.print(f"  [green]✓ Added {ref_spec.resname} to {Path(target_file).stem} at {target_chain_id}:{new_resid}[/green]")
                else:
                    self.console.print(f"  [red]✗ Failed to add {ref_spec.resname} to {Path(target_file).stem}[/red]")

    def _collect_alignment_atoms(
        self,
        structures: List[Tuple[str, Structure]],
        ref_idx: int,
        mappings: List[ResidueMapping],
        mode: str,
        cumulative_count: int
    ) -> List[Tuple[List[Atom], List[Atom]]]:
        """Collect atoms for alignment from reference and target structures"""
        ref_file, ref_structure = structures[ref_idx]
        alignment_data = []

        # Use only the first N mappings for progressive alignment
        active_mappings = mappings[:cumulative_count]

        for i, (target_file, target_structure) in enumerate(structures):
            if i == ref_idx:
                continue

            ref_atoms = []
            target_atoms = []

            for mapping in active_mappings:
                # Get residues
                ref_spec = mapping.reference
                target_spec = mapping.get_spec_for_structure(i, is_reference=False)

                ref_res = self._get_residue(ref_structure, ref_spec)
                target_res = self._get_residue(target_structure, target_spec)

                if ref_res is None:
                    self.console.print(f"[yellow]⚠ Warning: Residue {ref_spec} not found in reference, skipping[/yellow]")
                    continue

                if target_res is None:
                    self.console.print(f"[yellow]⚠ Warning: Residue {target_spec} not found in {target_file}, skipping[/yellow]")
                    continue

                # Check compatibility
                is_compatible, warning = self._check_residue_compatibility(ref_res, target_res, mode)
                if not is_compatible:
                    self.console.print(f"[red]✗ {warning}[/red]")
                    continue
                if warning:
                    self.console.print(f"[yellow]⚠ {warning}[/yellow]")

                # Extract atoms
                ref_res_atoms = AlignmentMode.get_atoms(ref_res, mode)
                target_res_atoms = AlignmentMode.get_atoms(target_res, mode)

                if not ref_res_atoms:
                    self.console.print(f"[yellow]⚠ Warning: No atoms found for {ref_spec} in reference, skipping[/yellow]")
                    continue

                if not target_res_atoms:
                    self.console.print(f"[yellow]⚠ Warning: No atoms found for {target_spec} in {target_file}, skipping[/yellow]")
                    continue

                if len(ref_res_atoms) != len(target_res_atoms):
                    self.console.print(
                        f"[yellow]⚠ Warning: Atom count mismatch for {ref_spec} vs {target_spec} "
                        f"({len(ref_res_atoms)} vs {len(target_res_atoms)}), skipping[/yellow]"
                    )
                    continue

                ref_atoms.extend(ref_res_atoms)
                target_atoms.extend(target_res_atoms)

            if ref_atoms and target_atoms:
                alignment_data.append((ref_atoms, target_atoms))
            else:
                self.console.print(f"[red]✗ No valid atoms for alignment in {target_file}[/red]")

        return alignment_data

    def _get_residue(self, structure: Structure, spec: ResidueSpec) -> Optional[Residue]:
        """Get residue from structure based on specification"""
        for model in structure:
            for chain in model:
                if chain.id == spec.chain_id:
                    for residue in chain:
                        # residue.id is (hetfield, resseq, icode)
                        if residue.id[1] == spec.resid:
                            return residue
        return None

    def _check_residue_compatibility(
        self,
        res1: Residue,
        res2: Residue,
        mode: str
    ) -> Tuple[bool, Optional[str]]:
        """Check if two residues are compatible for alignment"""
        if res1.resname != res2.resname:
            # Get residue IDs for debugging
            res1_id = res1.get_full_id()
            res2_id = res2.get_full_id()
            res1_chain = res1_id[2]
            res1_resid = res1_id[3][1]
            res2_chain = res2_id[2]
            res2_resid = res2_id[3][1]

            if mode == AlignmentMode.FULL:
                return False, f"Residue names don't match ({res1_chain}:{res1.resname}:{res1_resid} vs {res2_chain}:{res2.resname}:{res2_resid}), cannot use full residue mode"
            else:
                return True, f"Residue names don't match ({res1_chain}:{res1.resname}:{res1_resid} vs {res2_chain}:{res2.resname}:{res2_resid}), using {mode} mode"
        return True, None

    def _save_structure(self, structure: Structure, output_file: str):
        """Save structure to PDB file"""
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_file)

    def _get_indexed_filename(self, base_name: str, extension: str = ".pdb") -> str:
        """Get next available indexed filename to avoid overwriting

        Args:
            base_name: Base filename without extension (e.g., "1CLL_aligned")
            extension: File extension (default: ".pdb")

        Returns:
            Next available filename (e.g., "1CLL_aligned.pdb", "1CLL_aligned_2.pdb", etc.)
        """
        # Check if base filename exists
        filename = base_name + extension
        if not os.path.exists(filename):
            return filename

        # Find next available index
        index = 2
        while True:
            filename = f"{base_name}_{index}{extension}"
            if not os.path.exists(filename):
                return filename
            index += 1

    def _save_results_to_workspace(self, workspace: Dict[str, Any]) -> Dict[str, Any]:
        """Save alignment results to workspace"""
        # Save final aligned structures
        self.console.print("\n[bold]Saving Aligned Structures...[/bold]")
        reference_aligned_file = None
        target_aligned_files = []  # List of (index, filename, structure) for target structures

        for i, (original_file, _) in enumerate(self.structures):
            if i in self.aligned_structures:
                # Get stem and strip any existing _aligned suffix to avoid _aligned_aligned_aligned...
                stem = Path(original_file).stem
                if stem.endswith("_aligned"):
                    stem = stem[:-8]  # Remove "_aligned" suffix
                # Strip any trailing numbers from previous indexing (e.g., "1CLL_aligned_2" -> "1CLL")
                import re
                stem = re.sub(r'_aligned(_\d+)?$', '', stem)

                base_name = stem + "_aligned"
                output_file = self._get_indexed_filename(base_name)
                self._save_structure(self.aligned_structures[i], output_file)
                self.console.print(f"[green]✓ {output_file}[/green]")

                # Store reference aligned structure file path
                if i == self.reference_idx:
                    reference_aligned_file = output_file
                else:
                    # Store target aligned structure info
                    target_aligned_files.append((i, output_file, self.aligned_structures[i]))

        # Store in workspace
        workspace = self.update_workspace(workspace, 'aligned_structures', self.aligned_structures)
        workspace = self.update_workspace(workspace, 'alignment_results', self.alignment_results)
        workspace = self.update_workspace(workspace, 'alignment_residues', self.residue_mappings)

        # Save reference aligned structure with dedicated keys
        if reference_aligned_file:
            workspace = self.update_workspace(workspace, 'aligned_ref_pdb_file', str(Path(reference_aligned_file).absolute()))
            workspace = self.update_workspace(workspace, 'aligned_ref_structure', self.aligned_structures[self.reference_idx])
            self.console.print(f"\n[cyan]Reference aligned structure saved to workspace:[/cyan] aligned_ref_pdb_file")

        # If there are target structures, prompt user to select one for downstream processing
        if target_aligned_files:
            # Overlay reference + every aligned target in the viewer so the
            # user can visually compare them while picking the one they
            # want for downstream processing. The user's selection below
            # then triggers update_workspace('aligned_target_pdb_file'),
            # whose auto-launch hook collapses the overlay back to a
            # single-structure view of the chosen target.
            try:
                from proprep.structure_prep.viewer_coordinator import viewer as _viewer
                overlay_paths = []
                if reference_aligned_file:
                    overlay_paths.append(str(Path(reference_aligned_file).absolute()))
                for _, target_path, _ in target_aligned_files:
                    overlay_paths.append(str(Path(target_path).absolute()))
                if overlay_paths:
                    _viewer.show_structures(overlay_paths)
            except Exception:
                pass

            self.console.print(f"\n[bold]Select Target Structure for Downstream Processing[/bold]")
            self.console.print("[grey50]Choose which aligned target structure to make available for subsequent modules[/grey50]\n")

            # Display table of target structures
            from rich.table import Table
            table = Table()
            table.add_column("#", style="cyan")
            table.add_column("File", style="green")
            table.add_column("Original", style="yellow")

            for idx, (struct_idx, filename, _) in enumerate(target_aligned_files, 1):
                original_name = Path(self.structures[struct_idx][0]).name
                table.add_row(str(idx), filename, original_name)

            self.console.print(table)

            # Build options for prompt
            choices = [str(i) for i in range(1, len(target_aligned_files) + 1)] + ["n"]
            options_map = {str(i): target_aligned_files[i-1][1] for i in range(1, len(target_aligned_files) + 1)}
            options_map["n"] = "None (skip)"

            from .structure_loader import prompt_with_context
            choice = prompt_with_context(
                processor=self.processor,
                prompt="Select target structure (or 'n' to skip)",
                choices=choices,
                default="1" if len(target_aligned_files) == 1 else None,
                module="Structure Alignment",
                description="Select aligned target for downstream processing",
                options_map=options_map
            )

            if choice != "n":
                selected_idx = int(choice) - 1
                _, selected_file, selected_structure = target_aligned_files[selected_idx]

                # Save selected target to workspace
                workspace = self.update_workspace(workspace, 'aligned_target_pdb_file', str(Path(selected_file).absolute()))
                workspace = self.update_workspace(workspace, 'aligned_target_structure', selected_structure)
                self.console.print(f"\n[cyan]Target aligned structure saved to workspace:[/cyan] aligned_target_pdb_file")
                self.console.print(f"[green]✓ {selected_file} is now available for downstream modules[/green]")

        return workspace

    def _display_redox_sites(self, redox_sites: List):
        """Display available redox sites"""
        table = Table(title="Available Redox Sites")
        table.add_column("Index", style="cyan")
        table.add_column("Site ID", style="green")
        table.add_column("Type", style="yellow")
        table.add_column("Residues", style="magenta")

        for i, site in enumerate(redox_sites):
            site_id = getattr(site, 'site_id', f"Site_{i}")
            site_type = getattr(site, 'site_type', 'unknown')
            residue_count = len(getattr(site, 'residue_groups', {}))

            table.add_row(str(i + 1), site_id, site_type, str(residue_count))

        self.console.print(table)

    def _select_redox_sites(self, redox_sites: List) -> Optional[List[int]]:
        """Prompt user to select one or more redox sites

        Returns:
            List of selected site indices, or None if invalid
        """
        if not redox_sites:
            return None

        # Build options map for rich context
        options_map = {}
        for i, site in enumerate(redox_sites):
            site_id = getattr(site, 'site_id', f"Site_{i}")
            site_type = getattr(site, 'site_type', 'unknown')
            options_map[str(i + 1)] = f"{site_id} ({site_type})"

        self.console.print("\n[cyan]Enter site indices (e.g., '1' or '1,3,4' or '1-4')[/cyan]")

        response = prompt_with_context(
            processor=self.processor,
            prompt="Select redox site(s) for alignment",
            default="1",
            module="Structure Alignment - Redox Sites",
            description="Select one or more redox sites",
            options_map=options_map
        )

        # Parse the response (handle single number, comma-separated, or range)
        try:
            selected_indices = []

            # Handle ranges (e.g., "0-3")
            if '-' in response:
                parts = response.split('-')
                if len(parts) == 2:
                    start, end = int(parts[0].strip()), int(parts[1].strip())
                    selected_indices = list(range(start, end + 1))
                else:
                    raise ValueError("Invalid range format")
            # Handle comma-separated (e.g., "0,2,3")
            elif ',' in response:
                selected_indices = [int(x.strip()) for x in response.split(',')]
            # Handle single number (e.g., "0")
            else:
                selected_indices = [int(response.strip())]

            # Validate all indices (displayed 1-based to the user)
            for idx in selected_indices:
                if idx < 1 or idx > len(redox_sites):
                    self.console.print(f"[red]Invalid index {idx} (must be 1-{len(redox_sites)})[/red]")
                    return None

            # Convert back to 0-based for use as list subscripts downstream
            return [idx - 1 for idx in selected_indices]

        except ValueError as e:
            self.console.print(f"[red]Invalid input format. Please enter a number, comma-separated numbers, or a range (e.g., '1-4')[/red]")
            return None

    def _extract_residues_from_site(self, site) -> List[ResidueMapping]:
        """Extract residue specifications from RedoxSite object"""
        residue_mappings = []

        # Get unique residues from site's residue_groups
        if hasattr(site, 'residue_groups'):
            for (chain_id, resid, insertion_code), coords in site.residue_groups.items():
                if not coords:
                    continue

                # Get resname from coord_to_pdb mapping
                sample_coord = coords[0]
                if hasattr(site, 'coord_to_pdb') and sample_coord in site.coord_to_pdb:
                    resname = site.coord_to_pdb[sample_coord]['resname']
                    residue_spec = ResidueSpec(chain_id, resname, resid)
                    residue_mapping = ResidueMapping(residue_spec, {})
                    residue_mappings.append(residue_mapping)

        return residue_mappings

    def _extract_residues_from_sites(self, sites: List) -> List[ResidueMapping]:
        """Extract residue specifications from multiple RedoxSite objects

        Args:
            sites: List of RedoxSite objects

        Returns:
            Combined list of ResidueMapping objects from all sites
        """
        all_residue_mappings = []
        seen_residues = set()  # Track unique residues to avoid duplicates

        for site in sites:
            site_residues = self._extract_residues_from_site(site)

            # Add only unique residues (by chain_id and resid)
            for mapping in site_residues:
                residue_key = (mapping.reference.chain_id, mapping.reference.resid)
                if residue_key not in seen_residues:
                    seen_residues.add(residue_key)
                    all_residue_mappings.append(mapping)

        return all_residue_mappings

    def _map_hetatm_residues_interactive(
        self,
        hetatm_residues: List[ResidueMapping],
        protein_residues: List[ResidueMapping]
    ) -> List[ResidueMapping]:
        """Interactively map HETATM residues based on spatial proximity

        Args:
            hetatm_residues: List of HETATM residue mappings to process
            protein_residues: List of protein residue mappings (already aligned)

        Returns:
            Updated list of HETATM ResidueMapping objects (only those mapped to existing residues)

        Side effects:
            Populates self.hetatm_to_add with (ref_residue, target_indices) for HETATMs to add after alignment
        """
        import numpy as np

        # First, do a quick protein-based alignment to get transformed coordinates
        ref_structure = self.structures[self.reference_idx][1]
        target_structures = [
            (i, filename, struct)
            for i, (filename, struct) in enumerate(self.structures)
            if i != self.reference_idx
        ]

        # Perform alignment using protein residues
        aligned_target_structures = self._perform_preliminary_alignment(
            ref_structure,
            target_structures,
            protein_residues
        )

        # Now map each HETATM residue
        updated_hetatm = []  # Only HETATMs mapped to existing residues

        for hetatm_idx, hetatm_mapping in enumerate(hetatm_residues):
            ref_spec = hetatm_mapping.reference
            ref_residue = self._get_residue(ref_structure, ref_spec)

            if ref_residue is None:
                self.console.print(f"\n[yellow]Warning: {ref_spec} not found in reference, skipping[/yellow]")
                continue

            # Get reference residue center
            ref_center = self._get_residue_center(ref_residue)

            self.console.print(f"\n[bold cyan]Mapping {ref_spec}[/bold cyan]")

            # Track which targets will have this HETATM added (vs mapped)
            target_mappings = {}
            targets_to_add = []  # List of target indices where this HETATM will be added

            for target_idx, target_file, original_target, aligned_target, superimposer in aligned_target_structures:
                # Find all matching HETATM residues in target (same resname)
                candidates = self._find_hetatm_candidates(
                    aligned_target,
                    ref_spec.resname,
                    ref_center
                )

                # Build options table
                table = Table(show_header=True, box=box.SIMPLE)
                table.add_column("Option", style="cyan")
                table.add_column("Residue", style="green")
                table.add_column("Distance (Å)", style="yellow", justify="right")

                table.add_row("0", "[grey50]Skip (no match)[/grey50]", "-")

                options_map = {"0": "Skip"}

                # Add existing candidates
                if candidates:
                    self.console.print(f"\n  [bold]{target_file}:[/bold] Found {len(candidates)} candidate(s)")
                    for i, (chain_id, resid, distance) in enumerate(candidates, start=1):
                        table.add_row(
                            str(i),
                            f"{chain_id}:{ref_spec.resname}:{resid}",
                            f"{distance:.2f}"
                        )
                        options_map[str(i)] = f"{chain_id}:{ref_spec.resname}:{resid}"
                else:
                    self.console.print(f"\n  [bold]{target_file}:[/bold] No matching {ref_spec.resname} residues found")

                # Add option to create new residue at aligned position
                add_option_num = len(candidates) + 1
                table.add_row(
                    str(add_option_num),
                    f"[cyan]Add new {ref_spec.resname} at aligned position[/cyan]",
                    "-"
                )
                options_map[str(add_option_num)] = f"Add new {ref_spec.resname}"

                self.console.print(table)

                # Ask user to select
                choice = prompt_with_context(
                    processor=self.processor,
                    prompt=f"Select option for {ref_spec} in {target_file}",
                    default="0",
                    choices=[str(i) for i in range(add_option_num + 1)],
                    module="Structure Alignment - HETATM Mapping",
                    description=f"Map {ref_spec}",
                    options_map=options_map
                )

                choice_idx = int(choice)

                if choice_idx == add_option_num:
                    # User wants to add new residue - will be added after final alignment
                    targets_to_add.append(target_idx)
                    self.console.print(f"  [cyan]Will add new {ref_spec.resname} after alignment completes[/cyan]")
                elif choice_idx > 0 and choice_idx <= len(candidates):
                    # User selected an existing candidate - use for progressive alignment
                    selected_chain, selected_resid, selected_dist = candidates[choice_idx - 1]
                    target_spec = ResidueSpec(
                        selected_chain,
                        ref_spec.resname,
                        selected_resid
                    )
                    target_mappings[target_idx] = target_spec
                    self.console.print(f"  [green]✓ Mapped to {target_spec}[/green]")
                else:
                    self.console.print(f"  [grey50]Skipped[/grey50]")

            # If any targets were mapped to existing residues, add to alignment mappings
            if target_mappings:
                updated_mapping = ResidueMapping(ref_spec, target_mappings)
                updated_hetatm.append(updated_mapping)

            # If any targets need this HETATM added, track for post-alignment
            if targets_to_add:
                self.hetatm_to_add.append((ref_residue, ref_spec, targets_to_add))

        return updated_hetatm

    def _perform_preliminary_alignment(
        self,
        ref_structure: Structure,
        target_structures: List[Tuple[int, str, Structure]],
        protein_residues: List[ResidueMapping]
    ) -> List[Tuple[int, str, Structure, Structure, object]]:
        """Perform alignment using ALL protein residues from redox sites

        This creates the transformation matrix that will be used to calculate
        where HETATM residues would end up after alignment.

        Returns:
            List of (target_idx, filename, original_structure, aligned_structure, superimposer) tuples
        """
        from copy import deepcopy

        aligned_targets = []

        for target_idx, target_file, target_structure in target_structures:
            # Collect atoms for alignment using ALL protein residues
            ref_atoms = []
            target_atoms = []

            for mapping in protein_residues:  # Use ALL protein residues
                ref_spec = mapping.reference
                target_spec = mapping.get_spec_for_structure(target_idx, is_reference=False)

                ref_res = self._get_residue(ref_structure, ref_spec)
                target_res = self._get_residue(target_structure, target_spec)

                if ref_res and target_res:
                    # Use CA atoms for alignment (or all backbone atoms depending on mode)
                    if 'CA' in ref_res and 'CA' in target_res:
                        ref_atoms.append(ref_res['CA'])
                        target_atoms.append(target_res['CA'])

            if len(ref_atoms) >= 3:
                # Create a deep copy to avoid modifying the original
                target_copy = deepcopy(target_structure)

                # Perform alignment
                super_imposer = Superimposer()
                super_imposer.set_atoms(ref_atoms, target_atoms)

                # Apply transformation to all atoms in the copy
                atoms_to_transform = list(target_copy.get_atoms())
                super_imposer.apply(atoms_to_transform)

                self.console.print(f"  [grey50]Aligned {target_file} using {len(ref_atoms)} protein residues[/grey50]")

                aligned_targets.append((target_idx, target_file, target_structure, target_copy, super_imposer))
            else:
                self.console.print(f"[yellow]Warning: Not enough atoms for alignment of {target_file}[/yellow]")
                aligned_targets.append((target_idx, target_file, target_structure, target_structure, None))

        return aligned_targets

    def _find_hetatm_candidates(
        self,
        structure: Structure,
        resname: str,
        ref_center: np.ndarray
    ) -> List[Tuple[str, int, float]]:
        """Find HETATM residues of given type and calculate distances

        Args:
            structure: Target structure to search
            resname: Residue name to match (e.g., 'CA', 'HOH')
            ref_center: Reference residue center coordinates

        Returns:
            List of (chain_id, resid, distance) tuples, sorted by distance
        """
        import numpy as np

        candidates = []

        for model in structure:
            for chain in model:
                for residue in chain:
                    # Check if residue name matches
                    if residue.resname == resname:
                        # Calculate center of this residue
                        center = self._get_residue_center(residue)
                        # Calculate distance to reference
                        distance = np.linalg.norm(center - ref_center)

                        candidates.append((chain.id, residue.id[1], distance))
            break  # Only first model

        # Sort by distance
        candidates.sort(key=lambda x: x[2])

        return candidates

    def _get_residue_center(self, residue: Residue) -> np.ndarray:
        """Calculate geometric center of a residue

        Args:
            residue: BioPython Residue object

        Returns:
            NumPy array with (x, y, z) coordinates
        """
        import numpy as np

        coords = [atom.coord for atom in residue.get_atoms()]
        if not coords:
            return np.array([0.0, 0.0, 0.0])

        return np.mean(coords, axis=0)

    def _add_transformed_residue(
        self,
        target_structure: Structure,
        ref_residue: Residue,
        resname: str,
        superimposer
    ) -> Optional[Tuple[str, int]]:
        """Add a new residue to target structure with transformed coordinates

        Args:
            target_structure: Target structure to add residue to
            ref_residue: Reference residue to copy and transform
            resname: Residue name
            superimposer: Superimposer object with transformation matrix

        Returns:
            Tuple of (chain_id, residue_id) if successful, None otherwise
        """
        from Bio.PDB.Residue import Residue as BioResidue
        from Bio.PDB.Atom import Atom as BioAtom
        from copy import deepcopy
        import numpy as np

        # Find the next available residue ID and select target chain
        max_resid = 0
        target_chain = None

        for model in target_structure:
            chains = list(model.get_chains())

            if len(chains) == 0:
                self.console.print("[red]Error: No chains found in target structure[/red]")
                return None
            elif len(chains) == 1:
                # Only one chain, use it
                target_chain = chains[0]
            else:
                # Multiple chains, ask user
                self.console.print(f"\n[bold]Select chain for new {resname} residue:[/bold]")
                for i, chain in enumerate(chains, 1):
                    residue_count = len(list(chain))
                    self.console.print(f"  {i}. Chain {chain.id} ({residue_count} residues)")

                options_map = {}
                for i, chain in enumerate(chains, 1):
                    options_map[str(i)] = f"Chain {chain.id}"

                choice = prompt_with_context(
                    processor=self.processor,
                    prompt="Select chain",
                    default="1",
                    choices=[str(i) for i in range(1, len(chains) + 1)],
                    module="Structure Alignment - Add Residue",
                    description=f"Select chain for {resname}",
                    options_map=options_map
                )
                target_chain = chains[int(choice) - 1]

            # Find max residue ID across all chains
            for chain in model:
                for residue in chain:
                    resid = residue.id[1]
                    if resid > max_resid:
                        max_resid = resid
            break  # Only use first model

        # New residue ID is max + 1
        new_resid = max_resid + 1

        # Create a deep copy of the reference residue
        ref_residue_copy = deepcopy(ref_residue)

        # Apply transformation to all atoms in the copied residue
        atoms_to_transform = list(ref_residue_copy.get_atoms())
        superimposer.apply(atoms_to_transform)

        # Create a new residue with the new ID
        # BioPython residue ID is a tuple: (hetfield, resseq, icode)
        # Use 'H_' prefix for HETATM residues (ions, waters, ligands)
        hetfield = ref_residue.id[0] if ref_residue.id[0].strip() else ' '
        new_residue = BioResidue(
            (hetfield, new_resid, ' '),
            resname,
            ''  # segid
        )

        # Add transformed atoms to the new residue
        for atom in ref_residue_copy.get_atoms():
            # Create new atom with same properties but new coordinates
            new_atom = BioAtom(
                atom.name,
                atom.coord,  # Already transformed
                atom.bfactor,
                atom.occupancy,
                atom.altloc,
                atom.fullname,
                atom.serial_number,
                element=atom.element
            )
            new_residue.add(new_atom)

        # Add the new residue to the target chain
        try:
            target_chain.add(new_residue)
            return (target_chain.id, new_resid)
        except Exception as e:
            self.console.print(f"[red]Error adding residue: {e}[/red]")
            return None

    def _separate_protein_and_hetatm_residues(
        self,
        residue_mappings: List[ResidueMapping],
        reference_structure: Structure
    ) -> Tuple[List[ResidueMapping], List[ResidueMapping]]:
        """Separate protein residues from HETATM residues (ions, waters, ligands)

        Args:
            residue_mappings: List of residue mappings
            reference_structure: Reference structure to check residue types

        Returns:
            Tuple of (protein_residues, hetatm_residues)
        """
        # Try different import methods for BioPython compatibility
        try:
            from Bio.PDB.Polypeptide import is_aa
        except ImportError:
            from Bio.PDB.Polypeptide import is_aa

        protein_residues = []
        hetatm_residues = []

        for mapping in residue_mappings:
            # Find the residue in the reference structure
            ref_residue = self._get_residue(reference_structure, mapping.reference)

            if ref_residue is None:
                # Can't find it, assume it's HETATM to be safe
                hetatm_residues.append(mapping)
                continue

            # Check if it's a standard amino acid
            if is_aa(ref_residue, standard=True):
                protein_residues.append(mapping)
            else:
                hetatm_residues.append(mapping)

        return protein_residues, hetatm_residues

    def _extract_chain_sequence(self, structure: Structure, chain_id: str) -> Tuple[str, Dict[int, int]]:
        """Extract sequence from a chain and create position mapping

        Args:
            structure: BioPython Structure object
            chain_id: Chain identifier

        Returns:
            Tuple of (sequence_string, position_to_resid_mapping)
            position_to_resid_mapping maps alignment position (0-based) to actual residue number
        """
        # Try different import methods for BioPython compatibility
        try:
            from Bio.PDB.Polypeptide import three_to_one, is_aa
        except ImportError:
            # Newer BioPython versions
            from Bio.SeqUtils import seq1
            from Bio.PDB.Polypeptide import is_aa

            # Create wrapper for three_to_one using seq1
            def three_to_one(resname):
                return seq1(resname)

        sequence = []
        position_to_resid = {}
        position = 0

        for model in structure:
            for chain in model:
                if chain.id == chain_id:
                    for residue in chain:
                        # Only process standard amino acids
                        if is_aa(residue, standard=True):
                            try:
                                # Convert three-letter code to one-letter
                                one_letter = three_to_one(residue.resname)
                                sequence.append(one_letter)
                                position_to_resid[position] = residue.id[1]  # residue.id = (hetfield, resseq, icode)
                                position += 1
                            except (KeyError, ValueError):
                                # Non-standard residue, use 'X'
                                sequence.append('X')
                                position_to_resid[position] = residue.id[1]
                                position += 1
                    break  # Only process first model
            break

        return ''.join(sequence), position_to_resid

    def _create_automatic_residue_mappings(
        self,
        reference_mappings: List[ResidueMapping],
        reference_structure: Structure,
        target_structures: List[Tuple[int, str, Structure]]
    ) -> List[ResidueMapping]:
        """Automatically create residue mappings using sequence alignment

        Args:
            reference_mappings: List of reference residue specifications
            reference_structure: Reference Structure object
            target_structures: List of (index, filename, Structure) tuples for target structures

        Returns:
            List of ResidueMapping objects with automatic mappings filled in
        """
        from Bio.Align import PairwiseAligner

        # Group reference residues by chain
        residues_by_chain = {}
        for mapping in reference_mappings:
            chain_id = mapping.reference.chain_id
            if chain_id not in residues_by_chain:
                residues_by_chain[chain_id] = []
            residues_by_chain[chain_id].append(mapping)

        # Create new mappings with automatic target residue assignment
        new_mappings = []

        for chain_id, chain_mappings in residues_by_chain.items():
            # Extract reference sequence
            ref_seq, ref_pos_to_resid = self._extract_chain_sequence(reference_structure, chain_id)

            if not ref_seq:
                self.console.print(f"[yellow]Warning: Could not extract sequence for chain {chain_id} in reference[/yellow]")
                # Keep original mappings without target mappings
                new_mappings.extend(chain_mappings)
                continue

            self.console.print(f"\n[cyan]Chain {chain_id} reference sequence: {len(ref_seq)} residues[/cyan]")

            # For each target structure, align and map
            for target_idx, target_file, target_structure in target_structures:
                # Try to find matching chain in target
                target_chain_id = self._find_matching_chain(
                    target_structure, chain_id, ref_seq
                )

                if not target_chain_id:
                    self.console.print(
                        f"[yellow]Warning: No matching chain found for {chain_id} in {target_file}[/yellow]"
                    )
                    continue

                # Extract target sequence
                target_seq, target_pos_to_resid = self._extract_chain_sequence(
                    target_structure, target_chain_id
                )

                if not target_seq:
                    self.console.print(
                        f"[yellow]Warning: Could not extract sequence for chain {target_chain_id} in {target_file}[/yellow]"
                    )
                    continue

                self.console.print(
                    f"[cyan]  → {target_file} chain {target_chain_id}: {len(target_seq)} residues[/cyan]"
                )

                # Perform sequence alignment
                aligner = PairwiseAligner()
                aligner.mode = 'global'
                aligner.match_score = 2
                aligner.mismatch_score = -1
                aligner.open_gap_score = -2
                aligner.extend_gap_score = -0.5

                alignments = aligner.align(ref_seq, target_seq)
                best_alignment = alignments[0]

                # Create position mapping from alignment
                ref_aligned, target_aligned = best_alignment
                ref_align_pos = 0
                target_align_pos = 0
                alignment_map = {}  # ref_position -> target_position

                for ref_char, target_char in zip(ref_aligned, target_aligned):
                    if ref_char != '-' and target_char != '-':
                        # Both have residues at this position
                        alignment_map[ref_align_pos] = target_align_pos

                    if ref_char != '-':
                        ref_align_pos += 1
                    if target_char != '-':
                        target_align_pos += 1

                # Show alignment statistics
                identity = sum(1 for r, t in zip(ref_aligned, target_aligned) if r == t and r != '-')
                length = len(ref_aligned)
                self.console.print(
                    f"    Alignment: {identity}/{length} identical "
                    f"({100*identity/length:.1f}%), score: {best_alignment.score:.1f}"
                )

        # Now create mappings for each reference residue
        for mapping in reference_mappings:
            chain_id = mapping.reference.chain_id
            ref_resid = mapping.reference.resid

            # Extract reference sequence for this chain
            ref_seq, ref_pos_to_resid = self._extract_chain_sequence(reference_structure, chain_id)

            # Find position of this residue in reference sequence
            ref_position = None
            for pos, resid in ref_pos_to_resid.items():
                if resid == ref_resid:
                    ref_position = pos
                    break

            if ref_position is None:
                # Reference residue doesn't exist in structure - skip it
                continue

            # Create mapping dictionary for this residue
            target_mappings = {}

            for target_idx, target_file, target_structure in target_structures:
                # Find matching chain
                target_chain_id = self._find_matching_chain(target_structure, chain_id, ref_seq)
                if not target_chain_id:
                    continue

                # Get target sequence
                target_seq, target_pos_to_resid = self._extract_chain_sequence(
                    target_structure, target_chain_id
                )
                if not target_seq:
                    continue

                # Perform alignment
                aligner = PairwiseAligner()
                aligner.mode = 'global'
                aligner.match_score = 2
                aligner.mismatch_score = -1
                aligner.open_gap_score = -2
                aligner.extend_gap_score = -0.5
                alignments = aligner.align(ref_seq, target_seq)
                best_alignment = alignments[0]

                # Map position using alignment
                ref_aligned, target_aligned = best_alignment
                ref_align_pos = 0
                target_align_pos = 0

                for ref_char, target_char in zip(ref_aligned, target_aligned):
                    # Check if this is the reference position we're looking for
                    # IMPORTANT: Check ONLY when ref has a non-gap character
                    if ref_char != '-' and ref_align_pos == ref_position:
                        if target_char != '-':
                            # Found the corresponding target position
                            target_resid = target_pos_to_resid[target_align_pos]

                            # Verify this target residue actually exists in structure
                            target_residue = None
                            for model in target_structure:
                                for chain in model:
                                    if chain.id == target_chain_id:
                                        for residue in chain:
                                            if residue.id[1] == target_resid:
                                                target_residue = residue
                                                break
                                    if target_residue:
                                        break
                                if target_residue:
                                    break

                            if target_residue:
                                # Create mapping with actual target residue name
                                target_resname = target_residue.resname.strip()
                                target_spec = ResidueSpec(
                                    target_chain_id,
                                    target_resname,
                                    target_resid
                                )
                                target_mappings[target_idx] = target_spec
                        # If target has gap at this position, no mapping for this target
                        break

                    # Increment positions for non-gap characters
                    if ref_char != '-':
                        ref_align_pos += 1
                    if target_char != '-':
                        target_align_pos += 1

            # Only create mapping if at least one target has this residue
            if target_mappings:
                new_mapping = ResidueMapping(mapping.reference, target_mappings)
                new_mappings.append(new_mapping)

        return new_mappings

    def _find_matching_chain(
        self,
        structure: Structure,
        preferred_chain_id: str,
        reference_seq: str
    ) -> Optional[str]:
        """Find the best matching chain in target structure

        Args:
            structure: Target structure
            preferred_chain_id: Preferred chain ID (try this first)
            reference_seq: Reference sequence to compare

        Returns:
            Chain ID of best match, or None
        """
        from Bio.Align import PairwiseAligner

        # First try the same chain ID
        for model in structure:
            for chain in model:
                if chain.id == preferred_chain_id:
                    return preferred_chain_id
            break

        # If not found, try to find best sequence match
        best_chain_id = None
        best_score = 0

        aligner = PairwiseAligner()
        aligner.mode = 'global'

        for model in structure:
            for chain in model:
                target_seq, _ = self._extract_chain_sequence(structure, chain.id)
                if target_seq:
                    alignments = aligner.align(reference_seq, target_seq)
                    score = alignments[0].score
                    if score > best_score:
                        best_score = score
                        best_chain_id = chain.id
            break

        return best_chain_id

    def _download_pdb(self, pdb_id: str) -> Optional[str]:
        """Download PDB file from RCSB database

        Args:
            pdb_id: PDB identifier (e.g., '1CLL')

        Returns:
            Path to downloaded file, or None if download failed
        """
        import urllib.request
        from urllib.error import HTTPError

        pdb_id = pdb_id.upper()
        output_file = f"{pdb_id}.pdb"

        if os.path.exists(output_file):
            self.console.print(f"[yellow]File {output_file} already exists, using it[/yellow]")
            return output_file

        self.console.print(f"[cyan]Downloading {pdb_id} from RCSB PDB...[/cyan]")
        pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

        try:
            urllib.request.urlretrieve(pdb_url, output_file)
            self.console.print(f"[green]✓ Downloaded {output_file}[/green]")
            return output_file
        except HTTPError as e:
            if e.code == 404:
                self.console.print(f"[red]✗ PDB {pdb_id} not found (may be a large structure requiring mmCIF format)[/red]")
            else:
                self.console.print(f"[red]✗ Error downloading {pdb_id}: {e}[/red]")
            return None
        except Exception as e:
            self.console.print(f"[red]✗ Error downloading {pdb_id}: {e}[/red]")
            return None

    def _load_and_select_structures(self, pdb_files: List[str]) -> List[Tuple[str, Structure]]:
        """Load structures and allow model/chain selection for each

        Args:
            pdb_files: List of PDB/CIF file paths

        Returns:
            List of (filename, filtered_structure) tuples
        """
        structures = []

        for pdb_file in pdb_files:
            if not os.path.exists(pdb_file):
                self.console.print(f"[red]✗ {pdb_file} not found[/red]")
                continue

            # Load full structure using appropriate parser
            self.console.print(f"\n[bold]Loading {pdb_file}...[/bold]")
            try:
                parser = self._get_structure_parser(pdb_file)
                full_structure = parser.get_structure(Path(pdb_file).stem, pdb_file)
            except Exception as e:
                self.console.print(f"[red]✗ Error loading {pdb_file}: {e}[/red]")
                continue

            # Verify atom count
            is_valid, file_count, parsed_count, alt_count, explanation = self._verify_atom_count(pdb_file, full_structure)
            if not is_valid:
                if explanation:
                    self.console.print(
                        f"[cyan]ℹ {parsed_count} atoms parsed from {file_count} lines "
                        f"(difference due to {explanation})[/cyan]"
                    )
                else:
                    self.console.print(
                        f"[yellow]Warning: BioPython parsed {parsed_count} atoms "
                        f"but file contains {file_count} atoms[/yellow]"
                    )

            # Count models and chains
            num_models = len(list(full_structure.get_models()))
            all_chains = []
            for model in full_structure.get_models():
                all_chains.extend([chain.id for chain in model])
            num_chains = len(set(all_chains))  # Unique chains

            # Display structure information
            self.console.print(f"[cyan]Structure contains {num_models} model(s), {num_chains} chain(s), {parsed_count} atoms[/cyan]")

            # Ask if user wants to select specific model/chain
            use_full = confirm_with_context(
                processor=self.processor,
                prompt=f"Use entire structure from {pdb_file}?",
                default=True,
                module="Structure Alignment",
                description="Use full structure or select model/chain"
            )

            if use_full:
                structures.append((pdb_file, full_structure))
                self.console.print(f"[green]✓ Added full structure ({parsed_count} atoms)[/green]")
            else:
                # Model selection
                selected_structure = self._select_model_and_chains(pdb_file, full_structure)
                if selected_structure:
                    structures.append((pdb_file, selected_structure))
                    atom_count = sum(1 for _ in selected_structure.get_atoms())
                    self.console.print(f"[green]✓ Added selected portion ({atom_count} atoms)[/green]")

        return structures

    def _select_model_and_chains(self, filename: str, structure: Structure) -> Optional[Structure]:
        """Select specific model and chains from a structure

        Args:
            filename: PDB filename (for labeling)
            structure: Full BioPython Structure object

        Returns:
            New Structure with selected model/chains, or None if cancelled
        """
        from Bio.PDB.Structure import Structure as BioStructure
        from Bio.PDB.Model import Model as BioModel

        # Model selection
        models = list(structure)
        if len(models) > 1:
            self.console.print(f"\n[bold]Structure has {len(models)} models[/bold]")
            for i, model in enumerate(models):
                chain_count = len(list(model))
                self.console.print(f"  {i}. Model {model.id} ({chain_count} chains)")

            model_choice = int_prompt_with_context(
                processor=self.processor,
                prompt="Select model index",
                default=0,
                module="Structure Alignment",
                description="Select model"
            )
            selected_model = models[model_choice] if 0 <= model_choice < len(models) else models[0]
        else:
            selected_model = models[0]

        # Chain selection
        chains = list(selected_model)
        if len(chains) > 1:
            self.console.print(f"\n[bold]Model has {len(chains)} chains[/bold]")
            for i, chain in enumerate(chains, 1):
                residue_count = len(list(chain))
                self.console.print(f"  {i}. Chain {chain.id} ({residue_count} residues)")

            self.console.print(f"\n[yellow]Important: Select ALL chains you want in the final aligned structure.[/yellow]")
            self.console.print(f"[yellow]  • Alignment calculation: Only uses residues that match the reference structure[/yellow]")
            self.console.print(f"[yellow]  • Transformation applied: All atoms in all selected chains (as a rigid body)[/yellow]")
            self.console.print(f"[yellow]  • Final output: Contains all selected chains in their transformed positions[/yellow]")
            self.console.print(f"\n[cyan]Enter chain numbers to include (e.g., '1' or '1,2' or 'all'):[/cyan]")

            # Build options map for context
            options_map = {}
            for i, chain in enumerate(chains, 1):
                options_map[str(i)] = f"Chain {chain.id}"
            options_map["all"] = "All chains"

            chain_input = prompt_with_context(
                processor=self.processor,
                prompt="Select chains",
                default="all",
                module="Structure Alignment",
                description="Select chains to include",
                options_map=options_map
            )

            if chain_input.strip().lower() == "all":
                selected_chain_ids = [c.id for c in chains]
            else:
                # Parse numerical input (e.g., "1", "1,2", "1-3")
                selected_indices = []
                for part in chain_input.replace(',', ' ').split():
                    part = part.strip()
                    if '-' in part:
                        # Handle ranges like "1-3"
                        try:
                            start, end = part.split('-')
                            selected_indices.extend(range(int(start), int(end) + 1))
                        except ValueError:
                            self.console.print(f"[yellow]Warning: Invalid range '{part}', skipping[/yellow]")
                    else:
                        try:
                            selected_indices.append(int(part))
                        except ValueError:
                            self.console.print(f"[yellow]Warning: Invalid input '{part}', skipping[/yellow]")

                # Convert indices to chain IDs
                selected_chain_ids = []
                for idx in selected_indices:
                    if 1 <= idx <= len(chains):
                        selected_chain_ids.append(chains[idx - 1].id)
                    else:
                        self.console.print(f"[yellow]Warning: Chain number {idx} out of range, skipping[/yellow]")

                # If nothing valid selected, use all chains
                if not selected_chain_ids:
                    self.console.print("[yellow]No valid chains selected, using all chains[/yellow]")
                    selected_chain_ids = [c.id for c in chains]
        else:
            selected_chain_ids = [chains[0].id]

        # Build new structure with selected chains
        from copy import deepcopy

        new_structure = BioStructure(f"{Path(filename).stem}_selected")
        new_model = BioModel(selected_model.id)

        for chain in selected_model:
            if chain.id in selected_chain_ids:
                # Use deepcopy to avoid detaching from original structure
                new_model.add(deepcopy(chain))

        new_structure.add(new_model)
        return new_structure

    def _display_alignment_results(self):
        """Display alignment results table"""
        if not self.alignment_results:
            self.console.print("[yellow]No alignment results available[/yellow]")
            return

        table = Table(title="Alignment Results Summary", box=box.ROUNDED)
        table.add_column("Step", style="cyan")
        table.add_column("Structure", style="green")
        table.add_column("Residues", style="yellow")
        table.add_column("RMSD (Å)", style="magenta")

        for result in self.alignment_results:
            table.add_row(
                str(result['step']),
                Path(result['structure']).stem,
                result['residues'],
                f"{result['rmsd']:.3f}"
            )

        self.console.print(table)

    def cleanup(self):
        """Clean up module resources"""
        self.structures = []
        self.aligned_structures = {}
        self.alignment_results = []
        self.residue_mappings = []
