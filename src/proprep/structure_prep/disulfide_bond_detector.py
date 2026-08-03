"""
Disulfide Bond Module

Module for detecting and processing disulfide bonds in protein structures.
Integrates with the MPSA processor workflow using the command pattern.
"""

import logging
import os
from typing import Any, Dict, List

from rich.panel import Panel
from proprep.utils.prompts import prompt_with_context, confirm_with_context
from rich.status import Status
from rich.table import Table
from Bio.PDB import PDBParser

from proprep.utils.module_registry import ProcessingModule, register_module
from proprep.application.processor_command import ModuleActionCommand
from proprep.tleap_prep.bond_management_base import BondManagementBase
from .disulfide_bond_detector_worker import DisulfideBondDetector
from .disulfide_commands import (
    DetectDisulfideBondsCommand,
    ViewDisulfideBondsCommand,
    UpdateCysResiduesCommand,
    GenerateTLeapCommandsCommand,
)


logger = logging.getLogger(__name__)


@register_module
class DisulfideBondModule(ProcessingModule):
    """Module for detecting and processing disulfide bonds in protein structures"""

    NAME = "Disulfide Bond Detector"
    DESCRIPTION = "Detects disulfide bonds and updates CYS residues for MD simulations"
    VERSION = "1.0.0"
    CATEGORY = "analysis"
    REQUIRES = ["structure", "filtered_structure"]
    PRIORITY = 5

    def initialize(self):
        """Initialize module resources"""
        self.detector = DisulfideBondDetector()
        self.disulfide_bonds = []
        self.modified_structure = None
        self.output_file = None
        self._local_command_history = []

    def get_workspace_display(self, workspace_key, value, console):
        """
        Custom display method for Disulfide Bond Module data

        Args:
            workspace_key: The key in the workspace dictionary
            value: The value stored in the workspace
            console: The rich console object for output

        Returns:
            bool: True if the module handled the display, False otherwise
        """
        if workspace_key == "disulfide_bonds" and isinstance(value, list):
            if not value:
                console.print(
                    Panel(
                        "No disulfide bonds detected.",
                        title="Disulfide Bonds",
                        border_style="yellow",
                    )
                )
                return True

            table = Table(title="Detected Disulfide Bonds")
            table.add_column("ID", style="cyan")
            table.add_column("Chain 1", style="green")
            table.add_column("Residue 1", style="green")
            table.add_column("Chain 2", style="blue")
            table.add_column("Residue 2", style="blue")

            for i, (chain1, res_num1, chain2, res_num2) in enumerate(value, 1):
                table.add_row(str(i), chain1, str(res_num1), chain2, str(res_num2))

            console.print(Panel(table, title="Disulfide Bonds", border_style="green"))
            return True

        elif workspace_key == "disulfide_structure" and value is not None:
            chain_count = 0
            cys_count = 0
            cyx_count = 0

            try:
                for model in value:
                    for chain in model:
                        chain_count += 1

                        for residue in chain:
                            if residue.get_resname() == "CYS":
                                cys_count += 1
                            elif residue.get_resname() == "CYX":
                                cyx_count += 1
            except Exception as e:
                console.print(f"[yellow]Error processing structure: {str(e)}[/yellow]")

            console.print(
                Panel(
                    f"[bold]Structure with Updated Disulfide Bonds[/bold]\n"
                    f"Chains: {chain_count}\n"
                    f"CYS residues: {cys_count}\n"
                    f"CYX residues (disulfide-bonded): {cyx_count}\n"
                    f"Output file: {self.output_file}",
                    title="Disulfide-Modified Structure",
                    border_style="blue",
                )
            )
            return True

        elif workspace_key == "disulfide_tleap_commands" and isinstance(value, list):
            if not value:
                console.print(
                    Panel(
                        "No tLEaP commands generated.",
                        title="tLEaP Disulfide Commands",
                        border_style="yellow",
                    )
                )
                return True

            commands_text = "\n".join(value)
            console.print(
                Panel(
                    f"[bold]tLEaP Commands for Disulfide Bonds[/bold]\n\n"
                    f"{commands_text}",
                    title="tLEaP Disulfide Commands",
                    border_style="green",
                )
            )
            return True

        return False

    def get_menu_options(self) -> Dict[str, str]:
        """Get module menu options"""
        return {
            "detect": "Detect disulfide bonds",
            "view": "View detected disulfide bonds",
            "update": "Update CYS to CYX for disulfide bonds",
            "tleap": "Generate tLEaP commands for disulfide bonds",
        }

    def handle_menu_option(self, option: str) -> bool:
        """Handle a menu option selection"""
        try:
            if option == "detect":
                command = DetectDisulfideBondsCommand(
                    self.processor, interactive=True, use_ssbond_records=True
                )
                command.execute_with_error_handling()
                return True
            elif option == "view":
                command = ViewDisulfideBondsCommand(self.processor)
                command.execute_with_error_handling()
                return True
            elif option == "update":
                command = UpdateCysResiduesCommand(self.processor)
                command.execute_with_error_handling()
                return True
            elif option == "tleap":
                command = GenerateTLeapCommandsCommand(self.processor)
                command.execute_with_error_handling()
                return True
        except Exception as e:
            logger.error(f"Error executing menu option '{option}': {e}")
            self.processor.console.print(f"[red]Error: {e}[/red]")

        return False

    def detect_disulfide_bonds(
        self,
        workspace: Dict[str, Any] = None,
        interactive=True,
        use_ssbond_records=True,
    ) -> Dict[str, Any]:
        """Detect disulfide bonds in the structure."""
        if workspace is None:
            workspace = self.processor._get_workspace()

        structure, structure_source, file_path, using_file_path = (
            self._get_structure_info(workspace)
        )

        if structure is None and file_path is None:
            self.processor.console.print(
                "[yellow]No structure available. Load and filter a structure first.[/yellow]"
            )
            return workspace

        self.processor.console.print(
            f"[bold cyan]Using {structure_source} for disulfide bond detection[/bold cyan]"
        )

        input_file = file_path if using_file_path else workspace.get("pdb_file")
        
        # CRITICAL FIX: Don't apply chain filtering when using repaired structure
        # The repaired structure already contains only relevant chains and has different numbering
        if structure_source.startswith("repaired"):
            selected_chains = []  # No filtering - analyze all chains in repaired structure
            self.processor.console.print(
                "[grey50]Skipping chain filtering for repaired structure (different numbering)[/grey50]"
            )
        else:
            selected_chains = self._get_selected_chains(workspace)
            
        ssbond_records = self._get_ssbond_records(workspace)

        with Status("[bold green]Setting up disulfide bond detector..."):
            if using_file_path:
                self.detector.setup(
                    structure=None,
                    input_file=input_file,
                    selected_chains=selected_chains,
                )
            else:
                self.detector.setup(
                    structure=structure,
                    input_file=input_file,
                    selected_chains=selected_chains,
                )

        disulfide_bonds = self._detect_bonds(
            ssbond_records, input_file, selected_chains, using_file_path, interactive
        )

        workspace = self.update_workspace(workspace, "disulfide_bonds", disulfide_bonds)
        self._print_disulfide_summary(disulfide_bonds)

        if (
            disulfide_bonds
            and interactive
            and confirm_with_context(
                self.processor,
                "Update CYS residues to CYX for disulfide bonds?",
                default=True,
                module="Disulfide Bond Detector",
                description="Update CYS residues to CYX for detected disulfides",
            )
        ):
            workspace = self.update_cys_residues(workspace)

        return workspace

    def view_disulfide_bonds(self, workspace: Dict[str, Any] = None) -> Dict[str, Any]:
        """View detected disulfide bonds."""
        if workspace is None:
            workspace = self.processor._get_workspace()

        disulfide_bonds = workspace.get("disulfide_bonds")

        if not disulfide_bonds:
            self.processor.console.print(
                "[yellow]No disulfide bonds detected or loaded.[/yellow]"
            )
            if confirm_with_context(
                self.processor,
                "Run disulfide bond detection?",
                module="Disulfide Bond Detector",
                description="Run disulfide bond detection",
            ):
                workspace = self.detect_disulfide_bonds(workspace)
            return workspace

        table = Table(title="Detected Disulfide Bonds")
        table.add_column("ID", style="cyan")
        table.add_column("Chain 1", style="green")
        table.add_column("Residue 1", style="green")
        table.add_column("Chain 2", style="blue")
        table.add_column("Residue 2", style="blue")

        for i, (chain1, res_num1, chain2, res_num2) in enumerate(disulfide_bonds, 1):
            table.add_row(str(i), chain1, str(res_num1), chain2, str(res_num2))

        self.processor.console.print(table)
        return workspace

    def update_cys_residues(self, workspace: Dict[str, Any] = None) -> Dict[str, Any]:
        """Update CYS residues to CYX for disulfide bonds."""
        if workspace is None:
            workspace = self.processor._get_workspace()

        disulfide_bonds = workspace.get("disulfide_bonds")

        if not disulfide_bonds:
            self.processor.console.print(
                "[yellow]No disulfide bonds detected or loaded.[/yellow]"
            )
            if confirm_with_context(
                self.processor,
                "Run disulfide bond detection first?",
                module="Disulfide Bond Detector",
                description="Run disulfide bond detection before updating CYS",
            ):
                workspace = self.detect_disulfide_bonds(workspace)
            return workspace

        input_structure, input_file_path, using_file_path = self._get_input_source(
            workspace
        )

        if input_structure is None and input_file_path is None:
            self.processor.console.print(
                "[yellow]No structure available for updating CYS residues.[/yellow]"
            )
            return workspace

        output_file = prompt_with_context(
            self.processor,
            "Enter output filename for updated structure",
            default="structure_disulfide.pdb",
            module="Disulfide Bond Detector",
            description="Output filename for CYX-updated structure",
        )

        with Status("[bold green]Updating CYS residues to CYX..."):
            if using_file_path:
                self.detector.disulfide_bonds = disulfide_bonds
                saved_file = self.detector.update_cys_in_file(
                    input_file_path, output_file
                )
            else:
                self.detector.setup(
                    structure=input_structure,
                    selected_chains=self._get_selected_chains(workspace),
                )
                self.detector.disulfide_bonds = disulfide_bonds
                saved_file = self.detector.update_cys_residues(output_file=output_file)

        if saved_file:
            self.processor.console.print(
                f"[green]Updated structure saved to: {saved_file}[/green]"
            )
            workspace = self.update_workspace(
                workspace, "disulfide_output_file", saved_file
            )

            # Update the appropriate structure in workspace hierarchy
            # Priority: repaired > filtered > original
            # ALWAYS save as file path, not BioPython object
            
            # Check which structure we used and update it
            repair_results = workspace.get("repair_results", {})
            if repair_results.get("output_file") and os.path.exists(repair_results.get("output_file")):
                # We used repaired structure, update it
                workspace = self.update_workspace(
                    workspace, "repaired_structure", saved_file
                )
                # Also update repair_results to reflect the disulfide modification
                repair_results["output_file"] = saved_file
                repair_results["disulfide_modified"] = True
                workspace = self.update_workspace(
                    workspace, "repair_results", repair_results
                )
                self.processor.console.print(
                    "[cyan]Updated repaired structure with disulfide modifications[/cyan]"
                )
            elif workspace.get("repaired_structure") is not None:
                # Repaired structure exists but as wrong type, update it
                workspace = self.update_workspace(
                    workspace, "repaired_structure", saved_file
                )
                self.processor.console.print(
                    "[cyan]Updated repaired structure with disulfide modifications[/cyan]"
                )
            elif workspace.get("filtered_structure") is not None:
                # We used filtered structure, update it (but keep as file path)
                workspace = self.update_workspace(
                    workspace, "filtered_structure", saved_file
                )
                self.processor.console.print(
                    "[cyan]Updated filtered structure with disulfide modifications[/cyan]"
                )
            else:
                # We used original structure, update it
                workspace = self.update_workspace(
                    workspace, "structure", saved_file
                )
                self.processor.console.print(
                    "[cyan]Updated original structure with disulfide modifications[/cyan]"
                )
            
            # Also keep a reference as disulfide_structure for backwards compatibility
            workspace = self.update_workspace(
                workspace, "disulfide_structure", saved_file
            )
        else:
            self.processor.console.print(
                "[yellow]No structure was saved (possibly all CYS residues already named CYX)[/yellow]"
            )

        return workspace

    def generate_tleap_commands(
        self, workspace: Dict[str, Any] = None
    ) -> Dict[str, Any]:
        """Generate tLEaP commands for disulfide bonds with professional management interface."""
        if workspace is None:
            workspace = self.processor._get_workspace()

        disulfide_bonds = workspace.get("disulfide_bonds")

        if not disulfide_bonds:
            self.processor.console.print(
                "[yellow]No disulfide bonds detected or loaded.[/yellow]"
            )
            if confirm_with_context(
                self.processor,
                "Run disulfide bond detection first?",
                module="Disulfide Bond Detector",
                description="Run disulfide bond detection before generating tLEaP commands",
            ):
                workspace = self.detect_disulfide_bonds(workspace)
                disulfide_bonds = workspace.get("disulfide_bonds")
                if not disulfide_bonds:
                    return workspace

        # Initialize bond management system
        bond_manager = BondManagementBase(self.processor, "Disulfide Bond Detector")
        
        # Define categories for disulfide bonds
        categories = {
            "disulfide_bonds": "Inter-residue Disulfide Bonds", 
            "intra_residue_bonds": "Intra-residue Bonds"
        }
        bond_manager.initialize_bond_categories(categories)

        # Generate original disulfide bond commands
        detector = DisulfideBondDetector()
        detector.disulfide_bonds = disulfide_bonds
        commands = detector.get_tleap_commands()

        if not commands or not any(commands.values()):
            self.processor.console.print("[yellow]No tLEaP commands generated.[/yellow]")
            return workspace

        # Categorize disulfide bonds (inter-residue vs intra-residue)
        disulfide_commands = commands.get("disulfide_bonds", [])
        inter_residue_bonds = []
        intra_residue_bonds = []
        
        for cmd in disulfide_commands:
            if self._is_intra_residue_disulfide_command(cmd):
                intra_residue_bonds.append(cmd)
            else:
                inter_residue_bonds.append(cmd)

        # Add bonds to categories
        bond_manager.add_bonds_to_category("disulfide_bonds", inter_residue_bonds)
        bond_manager.add_bonds_to_category("intra_residue_bonds", intra_residue_bonds)
        
        # Note: intra_residue_bonds category is available for user to add custom bonds
        # The correct flow is: this generates bonds → stores in workspace → Topology Generator gathers them

        # Backup original bonds for reset functionality
        bond_manager.backup_original_bonds()

        # Show professional bond management interface using standard MenuCommand
        from proprep.tleap_prep.bond_management_commands import BondManagementMenuCommand
        
        self.processor.console.print("\n[bold cyan]Disulfide Bond tLEaP Command Generation[/bold cyan]")
        
        menu_command = BondManagementMenuCommand(
            self.processor, 
            bond_manager, 
            menu_category="Disulfide Bond Detector"
        )
        result = menu_command.execute()
        
        # Handle navigation and save results
        if result in ["back", "section", "main", "exit"]:
            # User chose navigation - return without saving
            return workspace
        elif result is True:
            # User chose to save and continue
            pdb_name = workspace.get("pdb_file", "unknown")
            if isinstance(pdb_name, str):
                pdb_name = pdb_name.split('/')[-1].replace('.pdb', '')
                
            output_file = prompt_with_context(
                self.processor,
                "Enter output filename",
                default=f"{pdb_name}_disulfide_bonds.tleap",
                module="Disulfide Bond Detector",
                description="Output filename for tLEaP commands file",
            )
            
            # Save to file
            if bond_manager.save_bonds_to_file(output_file, pdb_name):
                # Update workspace with final bond commands
                all_commands = []
                for category_bonds in bond_manager.bonds_by_category.values():
                    all_commands.extend(category_bonds)
                
                workspace = self.update_workspace(
                    workspace, "disulfide_tleap_commands", all_commands
                )

        return workspace

    def _is_intra_residue_disulfide_command(self, command: str) -> bool:
        """Check if a disulfide bond command is intra-residue."""
        # Parse command to extract residue identifiers
        # Format: bond mol.resid1.SG mol.resid2.SG
        import re
        pattern = r'bond\s+mol\.(\w+)\.SG\s+mol\.(\w+)\.SG'
        match = re.match(pattern, command.strip(), re.IGNORECASE)
        
        if match:
            resid1, resid2 = match.groups()
            return resid1 == resid2  # Same residue = intra-residue
            
        return False

    def _get_structure_info(self, workspace):
        """Get structure information in order of preference."""
        for key, source in [
            ("repaired_structure", "repaired structure"),
            ("filtered_structure", "filtered structure"),
            ("structure", "original structure"),
        ]:
            structure_obj = workspace.get(key)
            if structure_obj is not None:
                if isinstance(structure_obj, str) or hasattr(structure_obj, "as_posix"):
                    return None, f"{source} (file path)", str(structure_obj), True
                else:
                    return structure_obj, source, None, False

        return None, None, None, False

    def _get_selected_chains(self, workspace):
        """Get selected chains from filter selections."""
        filter_selections = workspace.get("filter_selections", {})
        return list(filter_selections.keys()) if filter_selections else []

    def _get_ssbond_records(self, workspace):
        """Get SSBOND records from metadata."""
        metadata = workspace.get("metadata")
        ssbond_records = []

        if metadata and hasattr(metadata, "ssbond_records") and metadata.ssbond_records:
            self.processor.console.print(
                f"[green]Found {len(metadata.ssbond_records)} SSBOND records in metadata[/green]"
            )
            for bond in metadata.ssbond_records:
                ssbond_records.append(
                    (bond["chain1"], bond["res1_num"], bond["chain2"], bond["res2_num"])
                )

        return ssbond_records

    def _detect_bonds(
        self, ssbond_records, input_file, selected_chains, using_file_path, interactive
    ):
        """Detect disulfide bonds using appropriate method."""
        if ssbond_records:
            self.processor.console.print(
                f"[green]Using {len(ssbond_records)} pre-parsed SSBOND records from metadata[/green]"
            )
            
            # Display the SSBOND records for debugging and user awareness
            self.processor.console.print("\n[grey50]SSBOND records from PDB file:[/grey50]")
            for i, (chain1, res1, chain2, res2) in enumerate(ssbond_records, 1):
                self.processor.console.print(
                    f"[grey50]  {i}. Chain {chain1}:{res1} ←→ Chain {chain2}:{res2}[/grey50]"
                )

            if selected_chains:
                filtered_bonds = [
                    bond
                    for bond in ssbond_records
                    if bond[0] in selected_chains and bond[2] in selected_chains
                ]
                if len(filtered_bonds) < len(ssbond_records):
                    self.processor.console.print(
                        f"[yellow]Filtered out {len(ssbond_records) - len(filtered_bonds)} bonds not in selected chains[/yellow]"
                    )
                ssbond_records = filtered_bonds

            # CRITICAL FIX: For repaired structures, always use structure-based verification
            # instead of file-based, because the file may have different chain/residue numbering
            workspace = self.processor._get_workspace() if hasattr(self.processor, '_get_workspace') else {}
            is_repaired_structure = workspace.get('repaired_structure') is not None
            
            if using_file_path and not is_repaired_structure:
                # Use file-based verification for original structures
                verified_bonds = self.detector.verify_ssbonds_in_file(
                    ssbond_records, input_file
                )
            else:
                # Use structure-based verification (with mapping) for repaired structures
                verified_bonds = self._verify_bonds_in_structure(ssbond_records)

            if len(verified_bonds) < len(ssbond_records):
                failed_count = len(ssbond_records) - len(verified_bonds)
                self.processor.console.print(
                    f"\n[yellow]Could only verify {len(verified_bonds)} out of {len(ssbond_records)} bonds in the structure[/yellow]"
                )
                if failed_count > 0 and interactive:
                    self.processor.console.print(
                        f"[cyan]💡 Tip: The unverified bonds may be due to chain/residue renumbering in processed structures.[/cyan]"
                    )
                    # Show which bonds could be verified vs which failed
                    if len(verified_bonds) > 0:
                        self.processor.console.print("\n[green]✓ Verified bonds:[/green]")
                        for chain1, res1, chain2, res2 in verified_bonds:
                            self.processor.console.print(f"[green]  Chain {chain1}:{res1} ←→ Chain {chain2}:{res2}[/green]")

            disulfide_bonds = verified_bonds
            if interactive:
                disulfide_bonds = self._review_disulfide_bonds(verified_bonds)
        else:
            self.processor.console.print(
                "[bold]Running disulfide bond detection...[/bold]"
            )

            if using_file_path:
                potential_bonds = self.detector.detect_bonds_from_file(
                    input_file, selected_chains
                )
            else:
                potential_bonds = self.detector._detect_bonds_by_distance()

            disulfide_bonds = (
                self._confirm_detected_bonds(potential_bonds)
                if interactive
                else potential_bonds
            )

        # Enhanced manual specification UX - MOVED OUTSIDE IF/ELSE
        # This ensures manual specification is offered regardless of detection method
        if interactive:
            if len(disulfide_bonds) == 0:
                self.processor.console.print(
                    "\n[bold yellow]🔍 No disulfide bonds found automatically.[/bold yellow]"
                )
                
                # Show available CYS residues to help user decide
                # Need to get structure for preview
                structure = self.detector.structure
                if using_file_path:
                    preview_cys = self.detector.get_cys_residues_from_file(input_file, selected_chains)
                else:
                    preview_cys = self.detector.get_cys_residues_from_structure(structure, selected_chains)
                
                if preview_cys:
                    self.processor.console.print(
                        f"[grey50]Found {len(preview_cys)} CYS/CYX residue(s) available for manual specification:[/grey50]"
                    )
                    for chain_id, res_id, res_type in preview_cys[:5]:  # Show first 5
                        self.processor.console.print(f"[grey50]  • Chain {chain_id} Res {res_id} ({res_type})[/grey50]")
                    if len(preview_cys) > 5:
                        self.processor.console.print(f"[grey50]  ... and {len(preview_cys) - 5} more[/grey50]")
                
                self.processor.console.print(
                    "[cyan]You can manually specify disulfide bonds if you know they exist.[/cyan]"
                )
                manual_prompt = "Would you like to manually specify disulfide bonds?"
                manual_default = True  # Default to Yes when no bonds found
            else:
                manual_prompt = f"\nFound {len(disulfide_bonds)} disulfide bond(s). Would you like to specify additional disulfide bonds manually?"
                manual_default = False  # Default to No when bonds already found
            
            if confirm_with_context(
                self.processor,
                manual_prompt,
                default=manual_default,
                module="Disulfide Bond Detector",
                description="Specify disulfide bonds manually",
            ):
                additional_bonds = self._prompt_for_manual_bonds(
                    input_file if using_file_path else None, selected_chains
                )
                if additional_bonds:
                    disulfide_bonds.extend(additional_bonds)
                    self.processor.console.print(
                        f"[green]✓ Added {len(additional_bonds)} manual disulfide bond(s)[/green]"
                    )

        self.detector.disulfide_bonds = disulfide_bonds
        return disulfide_bonds

    def _verify_bonds_in_structure(self, ssbond_records):
        """Verify SSBOND records against BioPython structure."""
        verified_bonds = []
        structure = self.detector.structure

        # Check if we're working with a repaired structure and need mapping
        # Use the proper workspace access method
        workspace = self.processor._get_workspace() if hasattr(self.processor, '_get_workspace') else {}
        residue_mapping = workspace.get('residue_mapping', {})
        filter_selections = workspace.get('filter_selections', {})
        
        # Determine if this is a repaired structure scenario
        using_repaired = workspace.get('repaired_structure') is not None
        
        if using_repaired and residue_mapping:
            self.processor.console.print(
                f"\n[grey50]Applying chain/residue mapping for repaired structure (MODELLER renumbering)[/grey50]"
            )
            # Show available chains in the repaired structure for debugging
            available_chains = []
            if structure:
                try:
                    available_chains = [chain.id for chain in structure[0]]
                except (TypeError, IndexError):
                    available_chains = []
            if available_chains:
                self.processor.console.print(f"[grey50]Available chains in repaired structure: {', '.join(available_chains)}[/grey50]")

        for chain1, res_num1, chain2, res_num2 in ssbond_records:
            # Apply mapping if available for repaired structure
            mapped_chain1, mapped_res1 = self._map_chain_residue(
                chain1, res_num1, residue_mapping, filter_selections, structure
            )
            mapped_chain2, mapped_res2 = self._map_chain_residue(
                chain2, res_num2, residue_mapping, filter_selections, structure
            )
            
            if using_repaired and residue_mapping and (
                (mapped_chain1 != chain1 or mapped_res1 != res_num1) or 
                (mapped_chain2 != chain2 or mapped_res2 != res_num2)
            ):
                self.processor.console.print(
                    f"[grey50]  Mapping bond {chain1}:{res_num1}-{chain2}:{res_num2} → {mapped_chain1}:{mapped_res1}-{mapped_chain2}:{mapped_res2}[/grey50]"
                )

            # For repaired structures, use file-based verification since BioPython
            # can't handle transformed structures with non-standard residue IDs
            if using_repaired and hasattr(self.detector, 'input_file') and self.detector.input_file:
                # Use file-based verification for transformed structures
                cys1_exists, cys2_exists = self._verify_cys_in_file(
                    self.detector.input_file, 
                    mapped_chain1, mapped_res1, 
                    mapped_chain2, mapped_res2
                )
            elif structure:
                # Use BioPython for normal structures
                cys1_exists = False
                cys2_exists = False
                for model in structure:
                    if mapped_chain1 in model and mapped_chain2 in model:
                        for residue in model[mapped_chain1]:
                            if residue.id[1] == mapped_res1 and residue.get_resname() in [
                                "CYS",
                                "CYX",
                            ]:
                                cys1_exists = True
                        for residue in model[mapped_chain2]:
                            if residue.id[1] == mapped_res2 and residue.get_resname() in [
                                "CYS",
                                "CYX",
                            ]:
                                cys2_exists = True
            else:
                cys1_exists = False
                cys2_exists = False

            if cys1_exists and cys2_exists:
                # Store the mapped coordinates for the verified bond
                verified_bonds.append((mapped_chain1, mapped_res1, mapped_chain2, mapped_res2))
            else:
                # More detailed failure message
                if using_repaired and residue_mapping and (
                    (mapped_chain1 != chain1 or mapped_res1 != res_num1) or 
                    (mapped_chain2 != chain2 or mapped_res2 != res_num2)
                ):
                    self.processor.console.print(
                        f"[yellow]✗ Bond {chain1}:{res_num1} ←→ {chain2}:{res_num2} not verified[/yellow]"
                    )
                    self.processor.console.print(
                        f"[yellow]  Attempted mapping: {chain1}:{res_num1} → {mapped_chain1}:{mapped_res1}, {chain2}:{res_num2} → {mapped_chain2}:{mapped_res2}[/yellow]"
                    )
                    # Check what residues are actually at those positions
                    if structure is not None:
                        for model in structure:
                            if mapped_chain1 in model:
                                for residue in model[mapped_chain1]:
                                    if residue.id[1] == mapped_res1:
                                        self.processor.console.print(
                                            f"[yellow]  Found at {mapped_chain1}:{mapped_res1}: {residue.get_resname()}[/yellow]"
                                        )
                                        break
                            if mapped_chain2 in model:
                                for residue in model[mapped_chain2]:
                                    if residue.id[1] == mapped_res2:
                                        self.processor.console.print(
                                            f"[yellow]  Found at {mapped_chain2}:{mapped_res2}: {residue.get_resname()}[/yellow]"
                                        )
                                        break
                else:
                    self.processor.console.print(
                        f"[yellow]✗ Bond {chain1}:{res_num1} ←→ {chain2}:{res_num2} not found in structure[/yellow]"
                    )

        return verified_bonds

    def _verify_cys_in_file(self, file_path, chain1, res1, chain2, res2):
        """
        Verify if CYS residues exist in the PDB file.
        
        Args:
            file_path: Path to PDB file
            chain1, res1: First CYS location
            chain2, res2: Second CYS location
            
        Returns:
            Tuple (cys1_exists, cys2_exists)
        """
        cys1_exists = False
        cys2_exists = False
        
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')):
                        chain_id = line[21]
                        try:
                            residue_id = int(line[22:26])
                            residue_name = line[17:20].strip()
                            
                            if chain_id == chain1 and residue_id == res1 and residue_name in ['CYS', 'CYX']:
                                cys1_exists = True
                            if chain_id == chain2 and residue_id == res2 and residue_name in ['CYS', 'CYX']:
                                cys2_exists = True
                                
                            if cys1_exists and cys2_exists:
                                return (True, True)
                        except ValueError:
                            continue
        except Exception as e:
            self.processor.console.print(f"[red]Error reading file for CYS verification: {e}[/red]")
            
        return (cys1_exists, cys2_exists)
    
    def _map_chain_residue(self, original_chain, original_res, residue_mapping, filter_selections, structure):
        """
        Map original chain/residue IDs to repaired structure numbering.
        
        Handles both chain renaming (C → A) and residue renumbering by MODELLER.
        
        Args:
            original_chain: Original chain ID (e.g., 'C')
            original_res: Original residue number (e.g., 112)
            residue_mapping: Mapping from structure completeness module
            filter_selections: Original filter selections
            structure: Current structure to check available chains
        
        Returns:
            Tuple of (mapped_chain, mapped_residue)
        """
        # Get available chains in the current structure
        available_chains = []
        if structure:
            try:
                available_chains = [chain.id for chain in structure[0]]
            except (TypeError, IndexError):
                available_chains = []
        
        # Step 1: Map the chain ID
        # MODELLER renames chains to A, B, C... sequentially
        if filter_selections and original_chain in filter_selections:
            # Chain was filtered - map to corresponding position in repaired structure
            filtered_chains = [k for k in filter_selections.keys() if k != "selected_model"]
            if filtered_chains and original_chain in filtered_chains:
                chain_index = filtered_chains.index(original_chain)
                # MODELLER uses A, B, C... sequentially
                # First filtered chain becomes A, second becomes B, etc.
                mapped_chain = chr(ord('A') + chain_index)
            else:
                # Default to A if not found
                mapped_chain = 'A'
        else:
            # No filter info - assume chain mapping if original not present
            # For single chain, MODELLER always outputs as chain A
            mapped_chain = 'A'
        
        # Step 2: Map the residue number using MODELLER mapping if available
        if residue_mapping and original_chain in residue_mapping:
            chain_mapping = residue_mapping[original_chain]
            original_to_new = chain_mapping.get('original_to_new', {})
            mapped_res = original_to_new.get(original_res, original_res)
        else:
            # No explicit mapping available - keep original residue number
            # This may work in some cases but is not reliable for MODELLER repairs
            mapped_res = original_res
            
        return mapped_chain, mapped_res

    def _review_disulfide_bonds(self, bonds):
        """Present discovered disulfide bonds to the user for review."""
        if not bonds:
            self.processor.console.print("[yellow]No disulfide bonds detected[/yellow]")
            return []

        table = Table(title="Detected Disulfide Bonds")
        table.add_column("ID", style="cyan")
        table.add_column("Chain 1", style="green")
        table.add_column("Residue 1", style="green")
        table.add_column("Chain 2", style="blue")
        table.add_column("Residue 2", style="blue")

        for i, (chain1, res_num1, chain2, res_num2) in enumerate(bonds, 1):
            table.add_row(str(i), chain1, str(res_num1), chain2, str(res_num2))

        self.processor.console.print(table)
        self.processor.console.print(
            "\n[bold]Please review each disulfide bond:[/bold]"
        )

        confirmed_bonds = []
        for i, (chain1, res_num1, chain2, res_num2) in enumerate(bonds, 1):
            self.processor.console.print(
                f"\nDisulfide {i}: Chain {chain1} Res {res_num1} - Chain {chain2} Res {res_num2}"
            )

            if confirm_with_context(
                self.processor,
                "Keep this disulfide bond?",
                default=True,
                module="Disulfide Bond Detector",
                description=f"Keep disulfide bond {i}",
            ):
                confirmed_bonds.append((chain1, res_num1, chain2, res_num2))
                self.processor.console.print("[green]Bond confirmed[/green]")
            else:
                self.processor.console.print("[yellow]Bond excluded[/yellow]")

        if confirmed_bonds:
            self.processor.console.print(
                f"\n[green]Confirmed {len(confirmed_bonds)} disulfide bond(s)[/green]"
            )
        else:
            self.processor.console.print(
                "[yellow]No disulfide bonds confirmed[/yellow]"
            )

        return confirmed_bonds

    def _confirm_detected_bonds(self, potential_bonds):
        """Interactive confirmation of detected bonds."""
        confirmed_bonds = []

        if potential_bonds:
            self.processor.console.print(
                "\n=== Confirming Detected Disulfide Bonds ==="
            )
            self.processor.console.print(
                f"Found {len(potential_bonds)} potential disulfide bond(s)"
            )

            for i, (chain1, res_num1, chain2, res_num2) in enumerate(potential_bonds):
                self.processor.console.print(
                    f"\nDisulfide {i+1}: Chain {chain1} Residue {res_num1} - Chain {chain2} Residue {res_num2}"
                )

                if confirm_with_context(
                    self.processor,
                    "Accept this disulfide bond?",
                    default=True,
                    module="Disulfide Bond Detector",
                    description=f"Accept potential disulfide bond {i+1}",
                ):
                    confirmed_bonds.append((chain1, res_num1, chain2, res_num2))
                    self.processor.console.print(
                        "[green]Disulfide bond accepted.[/green]"
                    )
                else:
                    self.processor.console.print(
                        "[yellow]Disulfide bond rejected.[/yellow]"
                    )
        else:
            self.processor.console.print(
                "[yellow]No disulfide bonds detected.[/yellow]"
            )

        return confirmed_bonds

    def _prompt_for_manual_bonds(self, file_path, selected_chains):
        """Allow manual specification of disulfide bonds."""
        additional_bonds = []
        self.processor.console.print("\n=== Manual Disulfide Bond Specification ===")

        if file_path:
            cys_residues = self.detector.get_cys_residues_from_file(
                file_path, selected_chains
            )
        else:
            cys_residues = self.detector.get_cys_residues_from_structure(
                self.detector.structure, selected_chains
            )

        if not cys_residues:
            self.processor.console.print(
                "[yellow]No CYS/CYX residues found in the structure.[/yellow]"
            )
            return additional_bonds

        self.processor.console.print("\nAvailable CYS/CYX residues:")
        for idx, (chain_id, res_id, res_type) in enumerate(cys_residues, 1):
            self.processor.console.print(
                f"[{idx}] Chain {chain_id} - Residue {res_id} ({res_type})"
            )

        idx_to_residue = {idx: residue for idx, residue in enumerate(cys_residues, 1)}

        self.processor.console.print("\nInput options:")
        self.processor.console.print(
            "1. Select CYS/CYX residues by index number from the list above"
        )
        self.processor.console.print("2. Enter chain ID and residue number manually", highlight=False)
        self.processor.console.print("3. Cancel (no disulfide bonds)", highlight=False)

        input_method = ""
        while input_method not in ["1", "2", "3"]:
            input_method = prompt_with_context(
                self.processor,
                "Choose input method [1/2/3]",
                default="1",
                module="Disulfide Bond Detector",
                description="Select method for manual disulfide bond entry",
                options_map={
                    "1": "Select CYS/CYX residues by index",
                    "2": "Enter chain ID and residue number manually",
                    "3": "Cancel",
                },
            )

        input_method = int(input_method)

        if input_method == 1:
            additional_bonds.extend(self._input_by_index(idx_to_residue))
        elif input_method == 2:
            additional_bonds.extend(
                self._input_by_chain_residue(cys_residues, selected_chains)
            )
        else:
            # input_method == 3: Cancel
            self.processor.console.print("[yellow]Cancelled disulfide bond specification[/yellow]")
            return []

        return additional_bonds

    def _input_by_index(self, idx_to_residue):
        """Input disulfide bonds by selecting indices."""
        bonds = []

        while True:
            self.processor.console.print(
                "\nEnter two index numbers for disulfide bond (or 'done' to finish):"
            )
            self.processor.console.print("Format: idx1 idx2 (e.g., 1 4)")
            bond_input = prompt_with_context(
                self.processor,
                "> ",
                module="Disulfide Bond Detector",
                description="Enter two CYS indices for manual disulfide bond (or 'done')",
            )

            if bond_input.lower() == "done":
                break

            try:
                parts = bond_input.split()
                if len(parts) != 2:
                    self.processor.console.print(
                        "[yellow]Invalid format. Please enter two index numbers.[/yellow]"
                    )
                    continue

                idx1, idx2 = int(parts[0]), int(parts[1])

                if idx1 not in idx_to_residue:
                    self.processor.console.print(
                        f"[yellow]Error: Index {idx1} is not in the list[/yellow]"
                    )
                    continue
                if idx2 not in idx_to_residue:
                    self.processor.console.print(
                        f"[yellow]Error: Index {idx2} is not in the list[/yellow]"
                    )
                    continue

                chain1, res_num1, _ = idx_to_residue[idx1]
                chain2, res_num2, _ = idx_to_residue[idx2]

                bonds.append((chain1, res_num1, chain2, res_num2))
                self.processor.console.print(
                    f"[green]Added disulfide bond: Chain {chain1} Res {res_num1} - Chain {chain2} Res {res_num2}[/green]"
                )

            except ValueError:
                self.processor.console.print(
                    "[yellow]Invalid input. Please enter numeric indices.[/yellow]"
                )
            except Exception as e:
                self.processor.console.print(f"[yellow]Error: {str(e)}[/yellow]")

        return bonds

    def _input_by_chain_residue(self, cys_residues, selected_chains):
        """Input disulfide bonds by entering chain and residue numbers."""
        bonds = []

        while True:
            self.processor.console.print(
                "\nEnter disulfide bond information (or 'done' to finish):"
            )
            self.processor.console.print(
                "Format: chain1 resnum1 chain2 resnum2 (e.g., A 42 B 156)"
            )
            bond_input = prompt_with_context(
                self.processor,
                "> ",
                module="Disulfide Bond Detector",
                description="Enter chain+resnum pair for manual disulfide bond (or 'done')",
            )

            if bond_input.lower() == "done":
                break

            parts = bond_input.split()
            if len(parts) != 4:
                self.processor.console.print(
                    "[yellow]Invalid format. Please use: chain1 resnum1 chain2 resnum2[/yellow]"
                )
                continue

            try:
                chain1, res_num1, chain2, res_num2 = parts
                res_num1 = int(res_num1)
                res_num2 = int(res_num2)

                if selected_chains and (
                    chain1 not in selected_chains or chain2 not in selected_chains
                ):
                    self.processor.console.print(
                        f"[yellow]Warning: Chain(s) not in selected chains {selected_chains}[/yellow]"
                    )
                    if not confirm_with_context(
                        self.processor,
                        "Add anyway?",
                        default=False,
                        module="Disulfide Bond Detector",
                        description="Add bond despite validation warning",
                    ):
                        continue

                cys1_found = (chain1, res_num1) in [(c, r) for c, r, _ in cys_residues]
                cys2_found = (chain2, res_num2) in [(c, r) for c, r, _ in cys_residues]

                if not cys1_found:
                    self.processor.console.print(
                        f"[yellow]Warning: No CYS/CYX residue found at Chain {chain1} Res {res_num1}[/yellow]"
                    )
                    if not confirm_with_context(
                        self.processor,
                        "Add anyway?",
                        default=False,
                        module="Disulfide Bond Detector",
                        description="Add bond despite validation warning",
                    ):
                        continue

                if not cys2_found:
                    self.processor.console.print(
                        f"[yellow]Warning: No CYS/CYX residue found at Chain {chain2} Res {res_num2}[/yellow]"
                    )
                    if not confirm_with_context(
                        self.processor,
                        "Add anyway?",
                        default=False,
                        module="Disulfide Bond Detector",
                        description="Add bond despite validation warning",
                    ):
                        continue

                bonds.append((chain1, res_num1, chain2, res_num2))
                self.processor.console.print(
                    f"[green]Added disulfide bond: Chain {chain1} Res {res_num1} - Chain {chain2} Res {res_num2}[/green]"
                )

            except ValueError:
                self.processor.console.print(
                    "[yellow]Invalid residue numbers. Please enter numeric residue IDs.[/yellow]"
                )
            except Exception as e:
                self.processor.console.print(f"[yellow]Error: {str(e)}[/yellow]")

        return bonds

    def _print_disulfide_summary(self, bonds):
        """Print a summary of detected disulfide bonds."""
        if not bonds:
            self.processor.console.print(
                "\n[yellow]No disulfide bonds detected or specified.[/yellow]"
            )
            return

        self.processor.console.print("\n=== Disulfide Bonds Summary ===")
        table = Table(title="Detected Disulfide Bonds")
        table.add_column("Bond", style="cyan")
        table.add_column("Chain 1", style="green")
        table.add_column("Residue 1", style="green")
        table.add_column("Chain 2", style="blue")
        table.add_column("Residue 2", style="blue")

        for i, (chain1, res_num1, chain2, res_num2) in enumerate(bonds, 1):
            table.add_row(str(i), chain1, str(res_num1), chain2, str(res_num2))

        self.processor.console.print(table)

    def _get_input_source(self, workspace):
        """Get the input source for structure updating."""
        for key in ["repaired_structure", "filtered_structure", "structure"]:
            structure_obj = workspace.get(key)
            if structure_obj is not None:
                if isinstance(structure_obj, str) or hasattr(structure_obj, "as_posix"):
                    return None, str(structure_obj), True
                else:
                    return structure_obj, None, False

        return None, None, False

    def _record_local_command(self, command):
        """Record command in local history for module-specific tracking"""
        self._local_command_history.append(command)

    def get_breadcrumb_for_command(self, command):
        """Return breadcrumb string for a command"""
        action_to_run = getattr(command, "_action_to_run", None)
        if action_to_run:
            action_map = {
                "detect_disulfide_bonds": "Detect",
                "view_disulfide_bonds": "View",
                "update_cys_residues": "Update",
                "generate_tleap_commands": "tLEaP",
            }
            action_name = action_map.get(action_to_run, action_to_run)
            return f"{self.NAME} > {action_name}"
        return self.NAME

    def get_workspace_requirements(self) -> List[str]:
        """Get workspace requirements"""
        return ["structure"]

    def get_workspace_outputs(self) -> List[str]:
        """Get workspace outputs"""
        return ["disulfide_bonds", "disulfide_structure", "disulfide_tleap_commands"]

    def can_process(self, workspace: Dict[str, Any]) -> bool:
        """Check if the module can process the current workspace"""
        return any(
            workspace.get(key) is not None
            for key in ["repaired_structure", "filtered_structure", "structure"]
        )

    def process(self, workspace):
        """Process the workspace"""
        if self.can_process(workspace):
            if workspace.get("detect_disulfides", False):
                workspace = self.detect_disulfide_bonds(
                    workspace, interactive=False, use_ssbond_records=True
                )

                if workspace.get("update_disulfides", False) and workspace.get(
                    "disulfide_bonds"
                ):
                    workspace = self.update_cys_residues(workspace)

                if workspace.get("generate_tleap", False) and workspace.get(
                    "disulfide_bonds"
                ):
                    workspace = self.generate_tleap_commands(workspace)

        return workspace

    def cleanup(self):
        """Clean up module resources"""
        self.detector = None
        self.disulfide_bonds = []
        self.modified_structure = None
        self.output_file = None
        self._local_command_history.clear()
