"""
AltLoc Selector Module

Identifies and processes alternate locations (conformations) in protein structures.
"""

import logging
from collections import defaultdict
from typing import Any, Dict, List

from rich.table import Table
from Bio.PDB import PDBParser

from proprep.utils.debug_utils import debug_workspace
from proprep.utils.module_registry import ProcessingModule, register_module
from proprep.utils.prompts import prompt_with_context, confirm_with_context, float_prompt_with_context, int_prompt_with_context
from proprep.structure_prep.altloc_selector_worker import AltLocSelector
from proprep.application.processor_command import ModuleActionCommand
from .altloc_selector_commands import (
    IdentifyAlternateLocationsCommand,
    SelectConformationsCommand,
    ProcessStructureCommand,
    ShowSelectionResultsCommand,
)


logger = logging.getLogger(__name__)


@register_module
class AltLocSelectorModule(ProcessingModule):
    """Module for selecting alternate conformations in structures"""

    NAME = "AltLoc Selector"
    DESCRIPTION = "Select alternate conformations in protein structures"
    VERSION = "1.0.0"
    CATEGORY = "modification"
    PRIORITY = 2

    def initialize(self):
        """Initialize module resources"""
        self.selector = AltLocSelector(processor=self.processor)
        self.selection_results = None
        self._local_command_history = []

    def get_menu_options(self) -> Dict[str, str]:
        """Get module menu options"""
        return {
            "identify": "Identify alternate locations",
            "select": "Select conformations",
            "process": "Process structure with selected conformations",
            "view_results": "View selection results",
        }

    def handle_menu_option(self, option: str) -> bool:
        """Handle a menu option selection"""
        try:
            if option == "identify":
                command = IdentifyAlternateLocationsCommand(self.processor)
                command.execute_with_error_handling()
                return True
            elif option == "select":
                command = SelectConformationsCommand(self.processor)
                command.execute_with_error_handling()
                return True
            elif option == "process":
                command = ProcessStructureCommand(self.processor)
                command.execute_with_error_handling()
                return True
            elif option == "view_results":
                command = ShowSelectionResultsCommand(self.processor)
                command.execute_with_error_handling()
                return True
        except Exception as e:
            logger.error(f"Error executing menu option '{option}': {e}")
            self.processor.console.print(f"[red]Error: {e}[/red]")

        return False

    def identify_alternate_locations(
        self, workspace: Dict[str, Any] = None
    ) -> Dict[str, Any]:
        """Identify alternate locations in the structure"""
        if workspace is None:
            workspace = self.processor._get_workspace()

        completeness_results = workspace.get("completeness_results", {})
        alternate_locations = completeness_results.get("alternate_locations", {})

        if not alternate_locations:
            self.processor.console.print(
                "[yellow]No alternate locations found in completeness results[/yellow]"
            )
            self.processor.console.print(
                "Analyzing structure for alternate locations..."
            )

            # Use StructureSelector for consistent structure access
            from proprep.utils.structure_selector import StructureSelector
            selector = StructureSelector(workspace, self.processor.console, processor=self.processor)

            structure = selector.get_structure_by_key(
                "repaired_structure",
                fallback_keys=["filtered_structure", "structure"],
                require_exists=False  # BioPython objects
            )
            if not structure:
                self.processor.console.print(
                    "[red]No structure available to analyze[/red]"
                )
                return workspace

            pdb_file = selector.get_structure(silent=True)
            if not pdb_file:
                self.processor.console.print(
                    "[red]No PDB file available for analysis[/red]"
                )
                return workspace

            selected_chains = []
            filter_selections = workspace.get("filter_selections", {})
            if filter_selections:
                selected_chains = list(filter_selections.keys())

            self.selector.setup(
                structure=structure,
                input_file=pdb_file,
                selected_chains=selected_chains,
            )

            alt_locs = self.selector.identify_alt_locs()

            if not alt_locs:
                self.processor.console.print(
                    "[yellow]No alternate locations found in the structure[/yellow]"
                )
            else:
                self.processor.console.print(
                    f"[green]Found {len(alt_locs)} residues with alternate locations[/green]"
                )
        else:
            self.processor.console.print(
                "[green]Alternate locations found in completeness results[/green]"
            )
            self._setup_selector_from_completeness(alternate_locations, workspace)

        return workspace

    def select_conformations(self, workspace: Dict[str, Any] = None) -> Dict[str, Any]:
        """Select which conformation to keep for each residue"""
        if workspace is None:
            workspace = self.processor._get_workspace()

        if not getattr(self.selector, "alt_residues", None):
            self.processor.console.print(
                "[yellow]No alternate locations identified. Run 'Identify alternate locations' first.[/yellow]"
            )
            return workspace

        auto_select = False
        if confirm_with_context(
            processor=self.processor,
            prompt="Automatically select highest occupancy conformation for all residues?",
            default=True,
            module="AltLoc Selector",
            description="Auto-select highest occupancy conformations"
        ):
            auto_select = True

        selections = self.selector.select_conformations(
            auto_select=auto_select, interactive=not auto_select
        )

        if selections:
            self.processor.console.print(
                f"[green]Selected conformations for {len(selections)} residues[/green]"
            )
            workspace = self.update_workspace(
                workspace, "altloc_selections", selections
            )
            self.selection_results = selections

        return workspace

    def process_structure(self, workspace: Dict[str, Any] = None) -> Dict[str, Any]:
        """Process the structure with selected conformations"""
        if workspace is None:
            workspace = self.processor._get_workspace()

        if not getattr(self.selector, "selections", None):
            self.processor.console.print(
                "[yellow]No conformations selected. Run 'Select conformations' first.[/yellow]"
            )
            return workspace

        method = self.selector.get_output_options()

        default_name = "altloc_processed.pdb"
        pdb_id = None
        metadata = getattr(self.processor, "metadata", None)
        if metadata:
            header_info = getattr(metadata, "header_info", None)
            if header_info:
                pdb_id = header_info.get("pdb_id")

        if pdb_id:
            output_file = f"{pdb_id}_{default_name}"
        else:
            output_file = default_name

        result_file = self.selector.process_structure(
            output_file=output_file, method=method
        )

        if result_file:
            self.processor.console.print(
                f"[green]Processed structure saved to {result_file}[/green]"
            )

            try:
                parser = PDBParser(QUIET=True)
                processed_structure = parser.get_structure("processed", result_file)

                workspace = self.update_workspace(
                    workspace, "processed_structure", processed_structure
                )

                method_name = (
                    "remove_alternates" if method == 1 else "normalize_occupancies"
                )
                workspace = self.update_workspace(
                    workspace,
                    "altloc_processing_results",
                    {
                        "output_file": result_file,
                        "method": method_name,
                    },
                )

            except Exception as e:
                self.processor.console.print(
                    f"[yellow]Saved file but could not load processed structure: {str(e)}[/yellow]"
                )

        return workspace

    def show_selection_results(
        self, workspace: Dict[str, Any] = None
    ) -> Dict[str, Any]:
        """Show selection results"""
        if workspace is None:
            workspace = self.processor._get_workspace()

        if not self.selection_results:
            self.processor.console.print(
                "[yellow]No selection results available[/yellow]"
            )
            return workspace

        self.processor.console.print("\n[bold]AltLoc Selection Results:[/bold]")

        table = Table(title="Selected Conformations")
        table.add_column("Chain", style="cyan")
        table.add_column("Residue", style="magenta")
        table.add_column("Selected AltLoc", style="green")

        by_chain = defaultdict(list)
        for res_id, alt_loc in self.selection_results.items():
            try:
                # Safely handle tuple unpacking 
                if isinstance(res_id, (tuple, list)) and len(res_id) >= 4:
                    chain_id, res_name, res_num, ins_code = res_id[0], res_id[1], res_id[2], res_id[3]
                elif isinstance(res_id, str):
                    # Handle string representation
                    try:
                        import ast
                        parsed_res_id = ast.literal_eval(res_id)
                        if len(parsed_res_id) >= 4:
                            chain_id, res_name, res_num, ins_code = parsed_res_id[0], parsed_res_id[1], parsed_res_id[2], parsed_res_id[3]
                        else:
                            self.processor.console.print(f"[yellow]Skipping invalid residue ID: {res_id}[/yellow]")
                            continue
                    except (ValueError, SyntaxError):
                        self.processor.console.print(f"[yellow]Skipping unparseable residue ID: {res_id}[/yellow]")
                        continue
                else:
                    self.processor.console.print(f"[yellow]Skipping invalid residue ID format: {res_id}[/yellow]")
                    continue
                    
                ins_str = f":{ins_code}" if ins_code else ""
                by_chain[chain_id].append((f"{res_name} {res_num}{ins_str}", alt_loc))
            except Exception as e:
                self.processor.console.print(f"[yellow]Error processing residue {res_id}: {str(e)}[/yellow]")
                continue

        for chain_id, residues in sorted(by_chain.items()):
            for i, (res_str, alt_loc) in enumerate(
                sorted(residues, key=lambda x: x[0])
            ):
                if i == 0:
                    table.add_row(chain_id, res_str, alt_loc)
                else:
                    table.add_row("", res_str, alt_loc)

        self.processor.console.print(table)

        if workspace.has("altloc_processing_results"):
            results = workspace.get("altloc_processing_results")
            self.processor.console.print("\n[bold]Processing Results:[/bold]")
            method_str = (
                "Remove alternate conformations"
                if results.get("method") == "remove_alternates"
                else "Normalize occupancies"
            )
            self.processor.console.print(f"Method: {method_str}")
            self.processor.console.print(f"Output file: {results.get('output_file')}")

        return workspace

    def _record_local_command(self, command):
        """Record command in local history for module-specific tracking"""
        self._local_command_history.append(command)

    def get_breadcrumb_for_command(self, command):
        """Return breadcrumb string for a command"""
        action_to_run = getattr(command, "_action_to_run", None)
        if action_to_run:
            action_map = {
                "identify_alternate_locations": "Identify",
                "select_conformations": "Select",
                "process_structure": "Process",
                "show_selection_results": "View Results",
            }
            action_name = action_map.get(action_to_run, action_to_run)
            return f"{self.NAME} > {action_name}"
        return self.NAME

    def _setup_selector_from_completeness(self, alternate_locations, workspace):
        """Setup the selector from completeness results"""
        # Use StructureSelector for consistent structure access
        from proprep.utils.structure_selector import StructureSelector
        selector = StructureSelector(workspace, self.processor.console, processor=self.processor)

        structure = selector.get_structure_by_key(
            "repaired_structure",
            fallback_keys=["filtered_structure", "structure"],
            require_exists=False  # BioPython objects
        )
        if not structure:
            self.processor.console.print("[red]No structure available to analyze[/red]")
            return

        pdb_file = selector.get_structure(silent=True)
        if not pdb_file:
            self.processor.console.print(
                "[red]No PDB file available for analysis[/red]"
            )
            return

        selected_chains = []
        filter_selections = workspace.get("filter_selections", {})
        if filter_selections:
            selected_chains = list(filter_selections.keys())

        self.selector.setup(
            structure=structure, input_file=pdb_file, selected_chains=selected_chains
        )

        alt_residues = defaultdict(set)

        for method, chain_results in alternate_locations.items():
            for chain_id, alternates in chain_results.items():
                for alt in alternates:
                    if (
                        "chain_id" in alt
                        and "residue_number" in alt
                        and "residue_name" in alt
                    ):
                        chain_id = alt["chain_id"]
                        res_name = alt["residue_name"]
                        res_num = str(alt["residue_number"])
                        ins_code = alt.get("insertion_code", "")

                        # Ensure we create a proper tuple, not a string representation
                        res_id = (str(chain_id), str(res_name), str(res_num), str(ins_code))

                        if "altlocs" in alt:
                            for loc in alt["altlocs"]:
                                alt_residues[res_id].add(str(loc))
                        elif "atom_name" in alt:
                            if selected_chains and chain_id not in selected_chains:
                                continue

                            # Extract altloc from detection details if available
                            details = alt.get("details", "")
                            if "Alternate locations:" in details:
                                # Parse altloc from details like "Alternate locations: A"
                                altloc_match = details.split("Alternate locations:")[1].strip().split()[0]
                                alt_residues[res_id].add(altloc_match)
                            elif "atom_name" in alt:
                                # Default to 'A' if we can't determine the altloc
                                alt_residues[res_id].add('A')

        if alt_residues:
            self.selector.alt_residues = alt_residues
            self.processor.console.print(
                f"[green]Found {len(alt_residues)} residues with alternate locations[/green]"
            )
        else:
            self.selector.identify_alt_locs()

    def get_workspace_requirements(self) -> List[str]:
        """Get workspace requirements"""
        return ["structure"]

    def get_workspace_outputs(self) -> List[str]:
        """Get workspace outputs"""
        return ["processed_structure", "altloc_processing_results"]

    def process(self, workspace: Dict[str, Any]) -> Dict[str, Any]:
        """Process the workspace"""
        debug_workspace(
            "AltLocSelector", "receiving", workspace, self.processor.console
        )

        auto_process = workspace.get("auto_process_altloc", False)

        if auto_process:
            # Use StructureSelector for consistent structure access
            from proprep.utils.structure_selector import StructureSelector
            selector = StructureSelector(workspace, self.processor.console, processor=self.processor)

            structure = selector.get_structure_by_key(
                "repaired_structure",
                fallback_keys=["filtered_structure", "structure"],
                require_exists=False  # BioPython objects
            )
            pdb_file = selector.get_structure(silent=True)

            if structure and pdb_file:
                selected_chains = []
                filter_selections = workspace.get("filter_selections", {})
                if filter_selections:
                    selected_chains = list(filter_selections.keys())

                self.selector.setup(
                    structure=structure,
                    input_file=pdb_file,
                    selected_chains=selected_chains,
                )

                alt_locs = self.selector.identify_alt_locs()

                if alt_locs:
                    selections = self.selector.select_conformations(
                        auto_select=True, interactive=False
                    )

                    default_name = "altloc_processed.pdb"
                    pdb_id = None
                    metadata = workspace.get("metadata")
                    if metadata:
                        header_info = getattr(metadata, "header_info", None)
                        if header_info:
                            pdb_id = header_info.get("pdb_id")

                    if pdb_id:
                        output_file = f"{pdb_id}_{default_name}"
                    else:
                        output_file = default_name

                    result_file = self.selector.process_structure(
                        output_file=output_file, method=1
                    )

                    if result_file:
                        try:
                            parser = PDBParser(QUIET=True)
                            processed_structure = parser.get_structure(
                                "processed", result_file
                            )

                            workspace = self.update_workspace(
                                workspace, "processed_structure", processed_structure
                            )
                            workspace = self.update_workspace(
                                workspace,
                                "altloc_processing_results",
                                {
                                    "output_file": result_file,
                                    "method": "remove_alternates",
                                    "selections": len(selections),
                                },
                            )

                        except Exception as e:
                            self.processor.console.print(
                                f"[yellow]Automatic processing error: {str(e)}[/yellow]"
                            )

        debug_workspace(
            "AltLocSelector", "returning", workspace, self.processor.console
        )
        return workspace

    def cleanup(self):
        """Clean up module resources"""
        self.selector = None
        self.selection_results = None
        self._local_command_history.clear()
