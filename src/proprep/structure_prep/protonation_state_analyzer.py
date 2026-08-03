"""
Titratable Residue Protonation State Analyzer

This module analyzes the protonation state of titratable residues at a specified pH,
excluding metal-coordinating residues. It integrates with the PDB Structure Processor
workflow to provide insights into protein protonation states for simulation preparation.

Features:
- Predicts protonation states of titratable residues using PROPKA or default pKa values
- Excludes metal-coordinating residues from analysis
- Handles PROPKA errors by identifying problematic residues
- Recommends residue names for MD simulations (AMBER naming conventions)
- Supports preparation for constant pH simulations
- Generates pH profiles for titration analysis
- Updates structure with user-selected residue names
"""

import logging
import os
from typing import Any, Dict, List

from rich.table import Table

from proprep.utils.module_registry import ProcessingModule, register_module, registry
from proprep.utils.prompts import prompt_with_context, confirm_with_context
from .protonation_worker import ProtonationStateAnalyzer, PKA_VALUES
from .protonation_commands import (
    AnalyzeProtonationStatesCommand,
    RunPropkaCommand,
    ViewProtonationResultsCommand,
    SetMDResidueNamesCommand,
    UpdateStructureCommand,
    SaveProtonationResultsCommand,
    GeneratePHProfileCommand,
    CalculateIdealCapacitanceCommand,
    GenerateCapacitanceProfileCommand,
)

# Setup logging
logger = logging.getLogger(__name__)


@register_module
class ProtonationStateModule(ProcessingModule):
    """Module for analyzing protonation states of titratable residues"""

    NAME = "Protonation State Analyzer"
    DESCRIPTION = "Predict pKas and assign protonation states"
    VERSION = "1.0.0"
    CATEGORY = "analysis"
    REQUIRES = ["structure", "filtered_structure"]
    PRIORITY = 6  # After metal detector

    def initialize(self):
        """Initialize module resources"""
        # Import Rich UI components
        self.analyzer = ProtonationStateAnalyzer(processor=self.processor)
        self.analysis_results = None
        self.md_residue_names = None

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Module Display

    def get_workspace_display(self, workspace_key, value, console):
        """
        Custom display method for Protonation State Analyzer module data

        Args:
            workspace_key: The key in the workspace dictionary
            value: The value stored in the workspace
            console: The rich console object for output

        Returns:
            bool: True if the module handled the display, False otherwise
        """
        from collections import defaultdict

        from rich.panel import Panel
        from rich.table import Table
        from rich.text import Text

        # Handle protonation results display
        if workspace_key == "protonation_results" and isinstance(value, dict):
            if not value:
                console.print(
                    Panel(
                        "No protonation analysis results available.",
                        title="Protonation Analysis",
                        border_style="yellow",
                    )
                )
                return True

            # Organize results by chain
            by_chain = defaultdict(list)
            for key, res in value.items():
                chain = res.get("chain", "?")
                by_chain[chain].append(res)

            # Create a summary table for all chains
            summary_table = Table(title="Protonation States Summary")
            summary_table.add_column("Chain", style="cyan")
            summary_table.add_column("Residue Type", style="green")
            summary_table.add_column("Count", style="blue")
            summary_table.add_column("Protonated", style="magenta")
            summary_table.add_column("Deprotonated", style="yellow")

            # Calculate summary statistics by chain and residue type
            for chain in sorted(by_chain.keys()):
                # Group by residue type within chain
                by_type = defaultdict(list)
                for res in by_chain[chain]:
                    res_type = res.get("type", "?")
                    by_type[res_type].append(res)

                # Add statistics for each residue type
                for res_type in sorted(by_type.keys()):
                    residues = by_type[res_type]
                    total = len(residues)
                    protonated = sum(
                        1 for res in residues if res.get("protonated", False)
                    )
                    deprotonated = total - protonated

                    summary_table.add_row(
                        chain, res_type, str(total), str(protonated), str(deprotonated)
                    )

            # Get the pH value from analyzer if available
            ph_value = "Unknown"
            if hasattr(self, "analyzer") and hasattr(self.analyzer, "pH"):
                ph_value = f"{self.analyzer.pH:.2f}"

            # Get net charge from workspace
            net_charge = "Unknown"
            net_charge_value = self.get_from_workspace("net_charge")
            if net_charge_value is not None:
                net_charge = f"{net_charge_value:.2f}"

            # Create detailed tables for each chain (limited to keep display manageable)
            chain_tables = []
            max_chains_to_show = 3
            max_residues_per_chain = 10

            for i, chain in enumerate(sorted(by_chain.keys())):
                # Only show first few chains
                if i >= max_chains_to_show:
                    break

                residues = by_chain[chain]
                # Sort by residue number
                sorted_residues = sorted(residues, key=lambda x: x.get("number", 0))

                # Limit number of residues to display
                display_residues = sorted_residues[:max_residues_per_chain]
                remaining = len(sorted_residues) - max_residues_per_chain

                # Create table for this chain
                chain_table = Table(title=f"Chain {chain} Titratable Residues")
                chain_table.add_column("Residue", style="cyan", justify="right")
                chain_table.add_column("pKa", style="green", justify="right")
                chain_table.add_column("Protonated", style="magenta", justify="right")
                chain_table.add_column("Probability", style="yellow", justify="right")
                chain_table.add_column("Charge", style="blue", justify="right")

                # Add each residue
                for res in display_residues:
                    res_type = res.get("type", "?")
                    res_num = res.get("number", "?")
                    pka = res.get("pKa", 0.0)
                    protonated = res.get("protonated", False)
                    probability = res.get("probability", 0.0)
                    charge = res.get("charge", 0.0)

                    chain_table.add_row(
                        f"{res_type} {res_num}",
                        f"{pka:.2f}",
                        str(protonated),
                        f"{probability:.3f}",
                        f"{charge:.1f}",
                    )

                # Add a row indicating there are more residues if needed
                if remaining > 0:
                    chain_table.add_row(
                        f"... and {remaining} more residues", "", "", "", ""
                    )

                chain_tables.append(chain_table)

            # Indicate if there are more chains
            more_chains_note = ""
            if len(by_chain) > max_chains_to_show:
                more_chains_note = f"\n[italic]... and {len(by_chain) - max_chains_to_show} more chains[/italic]"

            # Combine all content
            console.print(
                Panel(
                    f"[bold]Protonation Analysis at pH {ph_value}[/bold]\n"
                    f"Net Charge: {net_charge}\n\n"
                    f"{summary_table}\n\n"
                    f"[bold]Detailed Residue Information:[/bold]",
                    title="Protonation State Analysis",
                    border_style="blue",
                )
            )

            # Print chain tables
            for table in chain_tables:
                console.print(table)

            # Print more chains note if needed
            if more_chains_note:
                console.print(more_chains_note)

            return True

        # Handle MD residue names display
        elif workspace_key == "md_residue_names" and isinstance(value, dict):
            if not value:
                console.print(
                    Panel(
                        "No MD residue names defined.",
                        title="MD Residue Names",
                        border_style="yellow",
                    )
                )
                return True

            # Organize by chain and residue type
            by_chain = defaultdict(list)

            for key, md_name in value.items():
                # Parse the key to get chain and residue info
                parts = key.split("_")
                if len(parts) >= 3:
                    chain = parts[0]
                    try:
                        res_num = int(parts[1])
                        res_type = "_".join(
                            parts[2:]
                        )  # Join in case there are underscores in residue type
                    except ValueError:
                        # Handle case with non-standard key format
                        res_num = "?"
                        res_type = "_".join(parts[1:])

                    by_chain[chain].append((res_num, res_type, md_name))

            # Create a summary table of residue name changes
            summary_table = Table(title="MD Residue Name Changes")
            summary_table.add_column("Original Type", style="cyan")
            summary_table.add_column("MD Name", style="green")
            summary_table.add_column("Count", style="yellow")

            # Track changes by residue type
            changes_by_type = defaultdict(lambda: defaultdict(int))

            for chain, residues in by_chain.items():
                for _, res_type, md_name in residues:
                    changes_by_type[res_type][md_name] += 1

            # Add rows for each change type
            for res_type in sorted(changes_by_type.keys()):
                for md_name, count in sorted(changes_by_type[res_type].items()):
                    summary_table.add_row(res_type, md_name, str(count))

            # Create detailed tables for each chain
            chain_tables = []
            max_chains_to_show = 3
            max_residues_per_chain = 10

            for i, chain in enumerate(sorted(by_chain.keys())):
                # Only show first few chains
                if i >= max_chains_to_show:
                    break

                # Sort residues by number
                residues = sorted(by_chain[chain], key=lambda x: x[0])

                # Limit display
                display_residues = residues[:max_residues_per_chain]
                remaining = len(residues) - max_residues_per_chain

                # Create table
                chain_table = Table(title=f"Chain {chain} Residue Names")
                chain_table.add_column("Residue", style="cyan")
                chain_table.add_column("Original Type", style="yellow")
                chain_table.add_column("MD Name", style="green")

                for res_num, res_type, md_name in display_residues:
                    chain_table.add_row(str(res_num), res_type, md_name)

                # Add row for remaining residues if needed
                if remaining > 0:
                    chain_table.add_row(f"... and {remaining} more", "", "")

                chain_tables.append(chain_table)

            # Indicate if there are more chains
            more_chains_note = ""
            if len(by_chain) > max_chains_to_show:
                more_chains_note = f"\n[italic]... and {len(by_chain) - max_chains_to_show} more chains[/italic]"

            # Combine all content
            console.print(
                Panel(
                    f"[bold]MD Residue Names for Simulation[/bold]\n"
                    f"Total residues with assigned names: {len(value)}\n\n"
                    f"{summary_table}",
                    title="MD Residue Names",
                    border_style="green",
                )
            )

            # Print chain tables
            for table in chain_tables:
                console.print(table)

            # Print more chains note if needed
            if more_chains_note:
                console.print(more_chains_note)

            return True

        # Handle updated structure display
        elif workspace_key == "updated_structure" and value is not None:
            # Count chains and residues
            chain_count = 0
            residue_count = 0
            modified_residues = 0

            # Track residue types after modification
            residue_types = defaultdict(int)

            # Check for MD residue names to compare against
            md_names = self.get_from_workspace("md_residue_names", {})

            # Process the structure
            try:
                for model in value:
                    for chain in model:
                        chain_count += 1

                        for residue in chain:
                            residue_count += 1

                            # Count residue types
                            res_type = residue.resname
                            residue_types[res_type] += 1

                            # Check if this residue was modified
                            if md_names:
                                chain_id = chain.id
                                res_num = residue.id[1]

                                # Check if any keys in md_names match this residue
                                for key in md_names:
                                    if key.startswith(f"{chain_id}_{res_num}_"):
                                        # This is likely a modified residue
                                        modified_residues += 1
                                        break
            except Exception as e:
                # Handle any errors in processing the structure
                console.print(f"[yellow]Error processing structure: {str(e)}[/yellow]")

            # Create residue type table
            type_table = Table(title="Residue Types in Updated Structure")
            type_table.add_column("Residue Type", style="cyan")
            type_table.add_column("Count", style="green")

            # Add most common residue types (limit to 10)
            for res_type, count in sorted(
                residue_types.items(), key=lambda x: x[1], reverse=True
            )[:10]:
                type_table.add_row(res_type, str(count))

            # Check if there are more types
            remaining_types = len(residue_types) - 10
            if remaining_types > 0:
                type_table.add_row(f"... and {remaining_types} more types", "")

            # Create the panel
            console.print(
                Panel(
                    f"[bold]Updated Structure for MD Simulation[/bold]\n"
                    f"Chains: {chain_count}\n"
                    f"Residues: {residue_count}\n"
                    f"Modified Residues: {modified_residues}\n\n"
                    f"{type_table}",
                    title="Structure with Updated Residue Names",
                    border_style="blue",
                )
            )

            return True

        # Handle pH profile display
        elif workspace_key == "ph_profile" and isinstance(value, dict):
            if not value or "pH" not in value or "net_charge" not in value:
                console.print(
                    Panel(
                        "No pH profile data available.",
                        title="pH Profile",
                        border_style="yellow",
                    )
                )
                return True

            # Create a table with key pH points
            table = Table(title="pH Profile Summary")
            table.add_column("pH", style="cyan")
            table.add_column("Net Charge", style="green")

            # Get data points
            ph_values = value["pH"]
            charges = value["net_charge"]

            if not ph_values or not charges:
                console.print(
                    Panel(
                        "Incomplete pH profile data.",
                        title="pH Profile",
                        border_style="yellow",
                    )
                )
                return True

            # Show a few key points from the profile
            num_points = len(ph_values)
            points_to_show = min(num_points, 8)

            if points_to_show < num_points:
                # Show evenly distributed points
                step = num_points // points_to_show
                indices = [i * step for i in range(points_to_show)]
                # Always include the last point
                if indices[-1] != num_points - 1:
                    indices[-1] = num_points - 1
            else:
                # Show all points
                indices = range(num_points)

            # Add rows to the table
            for i in indices:
                table.add_row(f"{ph_values[i]:.1f}", f"{charges[i]:.2f}")

            # Create the panel
            console.print(
                Panel(
                    f"[bold]pH Titration Profile[/bold]\n"
                    f"pH range: {ph_values[0]:.1f} to {ph_values[-1]:.1f}\n"
                    f"Net charge range: {min(charges):.2f} to {max(charges):.2f}\n\n"
                    f"{table}\n\n"
                    f"[italic]Note: This is a summary of the full profile. Use plotting software for visualization.[/italic]",
                    title="pH Titration Profile",
                    border_style="green",
                )
            )

            return True

        # Handle PROPKA results display
        elif workspace_key == "propka_pka_values" and isinstance(value, dict):
            if not value:
                console.print(
                    Panel(
                        "No PROPKA pKa results available.",
                        title="PROPKA Results",
                        border_style="yellow",
                    )
                )
                return True

            # Organize by chain and residue type
            by_chain = defaultdict(list)

            for key, pka in value.items():
                # Parse the key to get chain and residue info
                parts = key.split("_")
                if len(parts) >= 3:
                    chain = parts[0]
                    try:
                        res_num = int(parts[1])
                        res_type = "_".join(
                            parts[2:]
                        )  # Join in case there are underscores in residue type
                    except ValueError:
                        # Handle case with non-standard key format
                        res_num = "?"
                        res_type = "_".join(parts[1:])

                    by_chain[chain].append((res_num, res_type, pka))

            # Create a summary table by residue type
            summary_table = Table(title="PROPKA pKa Summary by Residue Type")
            summary_table.add_column("Residue Type", style="cyan", justify="right")
            summary_table.add_column("Count", style="blue", justify="right")
            summary_table.add_column("Avg pKa", style="green", justify="right")
            summary_table.add_column("Min pKa", style="yellow", justify="right")
            summary_table.add_column("Max pKa", style="red", justify="right")

            # Calculate statistics by residue type
            by_type = defaultdict(list)
            for chain in by_chain.values():
                for _, res_type, pka in chain:
                    by_type[res_type].append(pka)

            # Add rows for each residue type
            for res_type, pka_values in sorted(by_type.items()):
                count = len(pka_values)
                avg_pka = sum(pka_values) / count if count > 0 else 0
                min_pka = min(pka_values) if count > 0 else 0
                max_pka = max(pka_values) if count > 0 else 0

                summary_table.add_row(
                    res_type,
                    str(count),
                    f"{avg_pka:.2f}",
                    f"{min_pka:.2f}",
                    f"{max_pka:.2f}",
                )

            # Create the panel
            console.print(
                Panel(
                    f"[bold]PROPKA pKa Prediction Results[/bold]\n"
                    f"Total residues analyzed: {len(value)}\n\n"
                    f"{summary_table}",
                    title="PROPKA Results",
                    border_style="blue",
                )
            )

            return True

        return False

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Workspace helper functions

    def get_workspace(self):
        """
        Helper method to get the appropriate workspace object

        Returns:
            The current workspace object
        """
        return self.processor.workspace

    def get_from_workspace(self, key, default=None):
        """
        Helper method to get values from the processor's workspace

        Args:
            key: Key to retrieve
            default: Default value if key not found

        Returns:
            The value for the given key or default
        """
        return self.processor.workspace.get(key, default)

    def get_from_workspace_obj(self, workspace, key, default=None):
        """
        Helper method to get values from a workspace object

        Args:
            workspace: Workspace object
            key: Key to retrieve
            default: Default value if key not found

        Returns:
            The value for the given key or default
        """
        return workspace.get(key, default)

    def update_workspace(self, key, value):
        """
        Helper method to update the processor's workspace

        Args:
            key: Key to update
            value: New value
        """
        old_value = self.processor.workspace.get(key)
        self.processor.workspace.set(key, value)
        debug_enabled = self.processor.workspace.get("debug", False)

        # Debug output if enabled
        if debug_enabled and hasattr(self.processor, "console"):
            self._debug_value(key, value, old_value, "updated")

    def update_workspace_obj(self, workspace, key, value):
        """
        Helper method to update a workspace object

        Args:
            workspace: Workspace object
            key: Key to update
            value: New value

        Returns:
            Updated workspace object
        """
        debug_enabled = workspace.get("debug", False)
        old_value = workspace.get(key)

        workspace.set(key, value)

        # Debug output if enabled
        if debug_enabled and hasattr(self.processor, "console"):
            self._debug_value(key, value, old_value, "updated")

        return workspace

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Input/output with workspace

    def get_workspace_requirements(self) -> List[str]:
        """Get workspace requirements - needs at least one structure loaded"""
        return [
            "rcsb_pdb_file | local_pdb_file | alphafold_pdb_file | alphafill_pdb_file | alphafold_homolog_pdb_file | aligned_target_pdb_file | aligned_ref_pdb_file"
        ]

    def get_workspace_outputs(self) -> List[str]:
        """Get workspace outputs"""
        return [
            "protonation_results",
            "md_residue_names",
            "net_charge",
            "microstate_protonation_results",
            "protonation_method",
            "textbook_pkas",
            "propka_pka_values",
            "pdb2pqr_output_file",
            "pdb2pqr_pqr_file",
            "protonation_pdb_file",
            "pdb2pqr_target_ph",
            "constant_ph_residues",
            "constant_ph_data",
            "ideal_capacitance",
            "specific_capacitance",
            "capacitance_per_group",
            "capacitance_profile",
        ]

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Module Menu

    def get_menu_options(self):
        """Return available menu options for this module"""
        return {
            "analyze": "Analyze protonation states (PROPKA)",
            "capacitance": "Calculate ideal charge capacitance",
            "set_names": "Set MD residue names",
            "view": "View analysis results",
            "export": "Export analysis results",
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

        # Check workspace state — either single-structure or microstate results
        has_analysis = (workspace.get("protonation_results") is not None
                        or workspace.get("microstate_protonation_results"))

        # Option 1: Analyze protonation states (PROPKA) - needs a loaded structure
        if has_analysis:
            prot_status, prot_dep = OptionStatus.COMPLETED, ""
        elif self.can_process(workspace):
            prot_status, prot_dep = OptionStatus.AVAILABLE, ""
        else:
            prot_status = OptionStatus.BLOCKED
            prot_dep = self.availability_note(workspace) or "Load a structure first"
        options.append(MenuOption(
            key="1",
            description="Analyze protonation states (PROPKA)",
            status=prot_status,
            dependency_text=prot_dep,
        ))

        # Option 2: Calculate capacitance - requires analysis; ● once the
        # capacitance result ("ideal_capacitance") has been written.
        if workspace.get("ideal_capacitance") is not None:
            status, dep_text = OptionStatus.COMPLETED, ""
        elif has_analysis:
            status, dep_text = OptionStatus.READY, ""
        else:
            status = OptionStatus.BLOCKED
            dep_text = "[Need to analyze protonation states first] ○"

        options.append(MenuOption(
            key="2",
            description="Calculate ideal charge capacitance",
            status=status,
            dependency_text=dep_text
        ))

        # Option 3: Set MD names - requires analysis; ● once the names
        # ("md_residue_names") have been written to the workspace.
        if workspace.get("md_residue_names") is not None:
            status, dep_text = OptionStatus.COMPLETED, ""
        elif has_analysis:
            status, dep_text = OptionStatus.READY, ""
        else:
            status = OptionStatus.BLOCKED
            dep_text = "[Need to analyze protonation states first] ○"

        options.append(MenuOption(
            key="3",
            description="Set MD residue names",
            status=status,
            dependency_text=dep_text
        ))

        # Option 4: View results - requires analysis
        if has_analysis:
            status = OptionStatus.READY
            dep_text = ""
        else:
            status = OptionStatus.BLOCKED
            dep_text = "[Need to analyze protonation states first] ○"

        options.append(MenuOption(
            key="4",
            description="View analysis results",
            status=status,
            dependency_text=dep_text
        ))

        # Option 5: Export results - requires analysis
        if has_analysis:
            status = OptionStatus.READY
            dep_text = ""
        else:
            status = OptionStatus.BLOCKED
            dep_text = "[Need to analyze protonation states first] ○"

        options.append(MenuOption(
            key="5",
            description="Export analysis results",
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
        has_analysis = (workspace.get("protonation_results") is not None
                        or workspace.get("microstate_protonation_results"))

        if not has_analysis:
            if not self.can_process(workspace):
                return f"{self.availability_note(workspace) or 'A structure is required'}. Load one via the Structure Loader."
            return "Start by analyzing protonation states (option 1) to identify titratable residues and their states"
        else:
            # Check what actions have been completed
            has_capacitance = workspace.get("ideal_capacitance") is not None
            has_md_names = workspace.get("md_residue_names") is not None

            if not has_capacitance and not has_md_names:
                return "Calculate capacitance (option 2), set MD names (option 3), view results (option 4), or export (option 5)"
            elif not has_capacitance:
                return "✓ MD names set. Calculate capacitance (option 2), view (option 4), or export (option 5)"
            elif not has_md_names:
                return "✓ Capacitance calculated. Set MD names (option 3), view (option 4), or export (option 5)"
            else:
                return "✓ Analysis complete. View (option 4), export (option 5), or press [m] to return to the main menu"
        
    def handle_menu_option(self, option: str) -> bool:
        """Handle a menu option selection"""
        try:
            if option == "analyze":
                # Combined workflow: Analyze + optional PROPKA + optional pH titration profile
                # Check if microstate results exist before analysis
                microstate_results_before = self.get_from_workspace("microstate_protonation_results", {})

                command = AnalyzeProtonationStatesCommand(
                    self.processor, interactive=True
                )
                success = command.execute()

                # Check if microstate results were created during this analysis
                microstate_results_after = self.get_from_workspace("microstate_protonation_results", {})
                microstates_analyzed = len(microstate_results_after) > len(microstate_results_before)

                if success and not microstates_analyzed:
                    # Only ask for pH profile if we analyzed the default structure
                    # (microstates already got asked individually)
                    if confirm_with_context(
                        processor=self.processor,
                        prompt="Would you like to generate a pH titration profile?",
                        default=False,
                        module="Protonation State Analyzer",
                        description="Generate pH titration profile"
                    ):
                        ph_command = GeneratePHProfileCommand(
                            self.processor, ph_range=(0, 14), step=0.5
                        )
                        ph_command.execute()

                return success

            elif option == "view":
                command = ViewProtonationResultsCommand(self.processor)
                return command.execute()
                
            elif option == "capacitance":
                # Combined workflow: Calculate capacitance + optional pH profile
                command = CalculateIdealCapacitanceCommand(self.processor)
                success = command.execute()
                
                if success:
                    # Ask if user wants to generate capacitance pH profile
                    if confirm_with_context(
                        processor=self.processor,
                        prompt="Would you like to generate a capacitance pH profile?",
                        default=False,
                        module="Protonation State Analyzer",
                        description="Generate capacitance pH profile"
                    ):
                        profile_command = GenerateCapacitanceProfileCommand(self.processor)
                        profile_command.execute()
                
                return success
                
            elif option == "set_names":
                # Combined workflow: Set MD names + update structure
                command = SetMDResidueNamesCommand(self.processor, constant_ph=False)
                success = command.execute()
                
                if success:
                    # Ask if user wants to update structure with MD names
                    if confirm_with_context(
                        processor=self.processor,
                        prompt="Would you like to update the structure with MD names?",
                        default=True,
                        module="Protonation State Analyzer",
                        description="Update structure with MD names"
                    ):
                        update_command = UpdateStructureCommand(self.processor)
                        update_command.execute()
                
                return success
                
            elif option == "export":
                command = SaveProtonationResultsCommand(
                    self.processor, output_format="text"
                )
                return command.execute()
                
            return False
        except Exception as e:
            self.processor.console.print(
                f"[red]Error executing command: {str(e)}[/red]"
            )
            logger.error(f"Error in handle_menu_option for {option}: {str(e)}")
            return False

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Public Action Methods (for command pattern)

    def analyze_protonation_states(self, interactive=True, **kwargs):
        """Public action method for analyzing protonation states"""
        try:
            if interactive:
                self._analyze_protonation_states()
            else:
                # Handle non-interactive mode if needed
                self._analyze_protonation_states()
            return True
        except Exception as e:
            self.processor.console.print(f"[red]Error analyzing protonation states: {str(e)}[/red]")
            return False

    def run_propka(self, **kwargs):
        """Public action method for running PROPKA"""
        self._run_propka()

    def view_results(self, **kwargs):
        """Public action method for viewing results"""
        self._view_results()

    def set_residue_names(self, constant_ph=False, **kwargs):
        """Public action method for setting residue names"""
        self._set_residue_names()

    def update_structure(self, **kwargs):
        """Public action method for updating structure"""
        self._update_structure()

    def save_results(self, output_format="text", **kwargs):
        """Public action method for saving results"""
        self._save_results()

    def generate_ph_profile(self, ph_range=(0, 14), step=0.5, **kwargs):
        """Public action method for generating pH profile"""
        self._generate_ph_profile()

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # pKa estimations and renaming capabilities

    def _check_and_select_microstates(self):
        """
        Check if microstate PDB files exist in workspace and allow user to select one or more.

        Returns:
            list of selected microstate file paths, or None if no microstates or user wants default structure
        """
        console = self.processor.console
        workspace = self.get_workspace()

        # Check for generated microstate files
        microstate_data = workspace.get("generated_microstate_pdbs", None)

        if not microstate_data or len(microstate_data) == 0:
            # No microstates available
            return None

        # Microstates exist - ask user if they want to analyze them
        console.print(f"\n[bold cyan]Detected {len(microstate_data)} generated redox microstate PDB files[/bold cyan]")

        use_microstates = confirm_with_context(
            processor=self.processor,
            prompt="Would you like to analyze microstate(s) instead of the default structure?",
            default=True,
            module="Protonation State Analyzer",
            description="Use microstates"
        )

        if not use_microstates:
            return None

        # Display available microstates
        from rich.table import Table
        from pathlib import Path
        table = Table(title="Available Redox Microstates")
        table.add_column("#", style="cyan", justify="right")
        table.add_column("File", style="green")
        table.add_column("Summary", style="yellow")

        for i, ms_data in enumerate(microstate_data, 1):
            filename = Path(ms_data['file']).name
            summary = self._summarize_microstate_description(ms_data.get('description', ''))
            table.add_row(str(i), filename, summary)

        console.print(table)

        # Ask user to select one or more microstates

        console.print("\n[bold]Selection Options:[/bold]")
        console.print("  • Single number: e.g., '3'")
        console.print("  • Multiple: e.g., '1,3,5'")
        console.print("  • Range: e.g., '1-4'")
        console.print("  • All: enter 'all'")

        selection = prompt_with_context(
            processor=self.processor,
            prompt="Select microstate(s) to analyze",
            default="all",
            module="Protonation State Analyzer",
            description="Microstate selection"
        )

        # Parse selection
        selected_indices = self._parse_microstate_selection(selection, len(microstate_data))

        if not selected_indices:
            console.print("[yellow]No valid microstates selected. Using default structure.[/yellow]")
            return None

        # Get the selected microstate file paths
        selected_files = [microstate_data[i]['file'] for i in selected_indices]

        console.print(f"\n[green]Selected {len(selected_files)} microstate(s) for analysis[/green]")

        return selected_files

    def _parse_microstate_selection(self, selection: str, max_count: int):
        """Parse user's microstate selection string into list of indices (0-based)"""
        selection = selection.strip().lower()

        if selection == 'all':
            return list(range(max_count))

        indices = []
        parts = selection.split(',')

        for part in parts:
            part = part.strip()
            if '-' in part:
                # Range
                try:
                    start, end = part.split('-')
                    start_idx = int(start.strip()) - 1  # Convert to 0-based
                    end_idx = int(end.strip()) - 1
                    indices.extend(range(start_idx, end_idx + 1))
                except ValueError:
                    continue
            else:
                # Single number
                try:
                    idx = int(part) - 1  # Convert to 0-based
                    indices.append(idx)
                except ValueError:
                    continue

        # Validate and deduplicate
        valid_indices = [i for i in indices if 0 <= i < max_count]
        return sorted(list(set(valid_indices)))

    def _summarize_microstate_description(self, description: str) -> str:
        """Summarize a microstate description into a compact string.

        Counts sites per redox state. If a state has only 1-2 sites,
        their IDs are listed for quick identification.

        Args:
            description: Comma-separated "site_id:redox_spin" pairs

        Returns:
            e.g. "61 oxidized, 1 reduced (site_5), 23 unchanged"
        """
        if not description:
            return "N/A"

        # Parse "site_1:oxidized_low_spin, site_2:reduced_low_spin, ..."
        state_sites = {}  # redox_state -> [site_id, ...]
        for entry in description.split(','):
            entry = entry.strip()
            if ':' not in entry:
                continue
            site_id, state_str = entry.split(':', 1)
            site_id = site_id.strip()
            # Extract redox state (before the underscore separating spin)
            # Format is "oxidized_low_spin" or "unchanged_unchanged"
            redox = state_str.strip().split('_')[0]
            state_sites.setdefault(redox, []).append(site_id)

        if not state_sites:
            return "N/A"

        parts = []
        for redox_state in ['oxidized', 'reduced', 'unchanged']:
            sites = state_sites.get(redox_state, [])
            if not sites:
                continue
            part = f"{len(sites)} {redox_state}"
            if len(sites) <= 2:
                part += f" ({', '.join(sites)})"
            parts.append(part)

        # Include any other states not in the standard list
        for redox_state, sites in state_sites.items():
            if redox_state not in ['oxidized', 'reduced', 'unchanged']:
                part = f"{len(sites)} {redox_state}"
                if len(sites) <= 2:
                    part += f" ({', '.join(sites)})"
                parts.append(part)

        return ", ".join(parts) if parts else "N/A"

    def _analyze_protonation_states(self):
        """Analyze protonation states of titratable residues"""
        # Step 1: Check if user wants to analyze microstate(s)
        selected_microstate_files = self._check_and_select_microstates()

        if selected_microstate_files:
            # For microstates, collect all parameters ONCE before the loop
            console = self.processor.console
            console.print("\n[bold cyan]Batch Microstate Analysis Settings[/bold cyan]")
            console.print("These settings will be applied to all selected microstates.\n")

            # Redox site exclusion
            redox_excluded_residues = self._get_redox_site_residues_for_exclusion()

            # pH
            pH_str = prompt_with_context(
                processor=self.processor,
                prompt="Enter pH value for analysis",
                default="7.0",
                module="Protonation State Analyzer",
                description="Enter pH value for analysis"
            )
            try:
                pH = float(pH_str)
            except ValueError:
                console.print("[red]Invalid pH value, using default of 7.0[/red]")
                pH = 7.0

            # Threshold
            threshold_str = prompt_with_context(
                processor=self.processor,
                prompt="Enter protonation probability threshold (0-1)",
                default="0.5",
                module="Protonation State Analyzer",
                description="Enter protonation threshold"
            )
            try:
                threshold = float(threshold_str)
                if threshold < 0 or threshold > 1:
                    console.print("[red]Invalid threshold, using default of 0.5[/red]")
                    threshold = 0.5
            except ValueError:
                console.print("[red]Invalid threshold, using default of 0.5[/red]")
                threshold = 0.5

            # PROPKA
            use_propka = confirm_with_context(
                processor=self.processor,
                prompt="Run PROPKA for improved pKa prediction?",
                default=True,
                module="Protonation State Analyzer",
                description="Run PROPKA for pKa prediction"
            )

            # pH titration profile
            generate_ph_profile = confirm_with_context(
                processor=self.processor,
                prompt="Generate pH titration profiles?",
                default=False,
                module="Protonation State Analyzer",
                description="Generate pH titration profile"
            )

            batch_params = {
                'pH': pH,
                'threshold': threshold,
                'use_propka': use_propka,
                'generate_ph_profile': generate_ph_profile,
                'visualize': False,  # Skip per-microstate visualization in batch
            }

            console.print(f"\n[green]Analyzing {len(selected_microstate_files)} microstates with: "
                         f"pH={pH}, threshold={threshold}, PROPKA={'yes' if use_propka else 'no'}[/green]")

            # Analyze each microstate with pre-collected parameters
            from pathlib import Path
            for i, microstate_file in enumerate(selected_microstate_files, 1):
                # Reset PROPKA retry counter for each microstate
                self.analyzer.propka_run_count = 0

                console.print(f"\n[bold cyan]{'='*60}[/bold cyan]")
                console.print(f"[bold]Analyzing microstate {i}/{len(selected_microstate_files)}: {Path(microstate_file).name}[/bold]")
                console.print(f"[bold cyan]{'='*60}[/bold cyan]\n")
                self._analyze_single_structure(
                    microstate_file, "microstate",
                    redox_excluded_residues=redox_excluded_residues,
                    batch_params=batch_params
                )
            return

        # Step 2: No microstates selected, check for processed structures with priority:
        # 1. repaired_pdb_file (from structure completeness/repair)
        # 2. filtered_pdb_file (from PDB filter)
        # 3. Ask user to select interactively
        # Note: Excludes transformed structures (protonating an already-transformed structure doesn't make sense)
        from proprep.utils.structure_selector import StructureSelector

        workspace = self.processor._get_workspace()

        # Source name mapping for display
        source_mapping = {
            "transformed_pdb_file": "transformed",
            "preprocessing_protein_input": "preprocessing protein",
            "repaired_pdb_file": "repaired",
            "filtered_pdb_file": "filtered",
            "hstripped_pdb_file": "H-stripped",
            "rcsb_pdb_file": "RCSB PDB",
            "local_pdb_file": "local",
            "alphafold_pdb_file": "AlphaFold",
            "alphafold_homolog_pdb_file": "AlphaFold homolog",
            "aligned_target_pdb_file": "aligned target",
            "aligned_ref_pdb_file": "aligned reference",
        }

        selector = StructureSelector(
            workspace, self.processor.console, processor=self.processor
        )

        result = selector.get_structure(
            priority_override=[
                "transformed_pdb_file",  # Highest priority: post-redoxsite transformation
                "preprocessing_protein_input",  # Extracted protein from preprocessing
                "repaired_pdb_file",
                "filtered_pdb_file",
                "hstripped_pdb_file",
                "rcsb_pdb_file",
                "local_pdb_file",
                "alphafold_pdb_file",
                "alphafold_homolog_pdb_file",
                "aligned_target_pdb_file",
                "aligned_ref_pdb_file",
            ],
            exclude_keys=["protonation_pdb_file", "structure_with_prot_resnames"],
            interactive=True,
            return_key=True,
            silent=True,
        )

        if result is None:
            input_file_path = None
            structure_source = None
        else:
            input_file_path, file_key = result
            structure_source = source_mapping.get(file_key, "selected")
            self.processor.console.print(f"[cyan]Using {structure_source} structure: {input_file_path}[/cyan]")

        # Check if we have a file path
        if input_file_path is None:
            self.processor.console.print(
                "[yellow]No structure available. Load a structure first.[/yellow]"
            )
            return

        # Analyze the single default structure
        # Don't ask for pH profile here - it will be asked in handle_menu_option
        self._analyze_single_structure(input_file_path, structure_source, ask_ph_profile=False)

    def _analyze_single_structure(self, input_file_path, structure_source, ask_ph_profile=True, redox_excluded_residues=None, batch_params=None):
        """Analyze protonation states for a single structure file

        Args:
            input_file_path: Path to the structure file
            structure_source: Source description (e.g., "transformed", "microstate")
            ask_ph_profile: Whether to ask about pH profile generation (default True)
            redox_excluded_residues: Pre-determined redox site exclusions (optional, for microstates)
            batch_params: Pre-collected parameters dict to skip interactive prompts.
                Keys: pH, threshold, use_propka, generate_ph_profile, visualize
        """

        console = self.processor.console
        input_structure = None
        using_file_path = True

        # Set batch mode flag to suppress per-structure interactive prompts
        self._batch_mode = batch_params is not None

        # Establish the analyze-time structure baseline in the viewer.
        # Skipped in batch mode (where the viewer would flip for every
        # microstate) and silently no-ops if the coordinator can't
        # launch (CLI without prior viewer).
        if not self._batch_mode and input_file_path:
            try:
                from proprep.structure_prep.viewer_coordinator import (
                    viewer as _viewer,
                )
                _viewer.show_structure(input_file_path)
            except Exception:
                pass

        # IMPORTANT: For microstate analysis, treat as "transformed" type (analyze all chains)
        # But remember if it's a microstate for later
        is_microstate = (structure_source == "microstate")
        if is_microstate:
            structure_source = "transformed"

        # RedoxSite integration - Check for redox sites and allow user selection for exclusion
        # If redox_excluded_residues is provided (microstate batch mode), use it
        # Otherwise, ask the user
        if redox_excluded_residues is None:
            redox_excluded_residues = self._get_redox_site_residues_for_exclusion()

        excluded_residues = redox_excluded_residues if redox_excluded_residues else set()

        if excluded_residues:
            self.processor.console.print(
                f"[green]Excluding {len(excluded_residues)} residues from protonation analysis[/green]"
            )
        else:
            self.processor.console.print(
                "[yellow]No residues excluded from protonation analysis[/yellow]"
            )

        # Get selected chains and residues from filter selections
        selected_chains = []
        selected_residues = {}

        # IMPORTANT: If using transformed structure or microstate, don't use filter_selections
        # because residue IDs have been renumbered. Analyze all residues instead.
        if structure_source == "transformed":
            # For transformed structures/microstates, analyze ALL chains/residues (excluding redox sites)
            # Get all chains from the structure
            if using_file_path:
                all_chains = set()
                with open(input_file_path, 'r') as f:
                    for line in f:
                        if line.startswith(("ATOM", "HETATM")):
                            all_chains.add(line[21])
                selected_chains = list(all_chains)
            else:
                # Get chains from structure object
                selected_chains = [chain.id for chain in input_structure[0]]

            self.processor.console.print(
                f"[cyan]Analyzing all residues in {structure_source} structure (chains: {', '.join(selected_chains)})[/cyan]"
            )
        else:
            # For non-transformed structures, use filter_selections if available
            filter_selections = self.get_from_workspace("filter_selections", {})

            if filter_selections:
                selected_chains = list(filter_selections.keys())

                # Check for cross-chain selected residues
                for chain_id, component_selections in filter_selections.items():
                    for comp_type, residues in component_selections.items():
                        for resid in residues:
                            # Find residue name based on input type
                            if using_file_path:
                                resname = self._find_residue_name_from_file(
                                    input_file_path, chain_id, resid
                                )
                            else:
                                resname = self._find_residue_name(
                                    input_structure, chain_id, resid
                                )

                            if resname and chain_id not in selected_chains:
                                if chain_id not in selected_residues:
                                    selected_residues[chain_id] = []
                                selected_residues[chain_id].append((resid, resname))
            else:
                # No filter_selections: analyze ALL chains (e.g., after structure repair without filtering)
                if using_file_path:
                    all_chains = set()
                    with open(input_file_path, 'r') as f:
                        for line in f:
                            if line.startswith(("ATOM", "HETATM")):
                                all_chains.add(line[21])
                    selected_chains = list(all_chains)
                else:
                    # Get chains from structure object
                    selected_chains = [chain.id for chain in input_structure[0]]

                self.processor.console.print(
                    f"[cyan]No filter selections found - analyzing all chains: {', '.join(selected_chains)}[/cyan]"
                )

        # Get analysis parameters (from batch_params or interactively)
        if batch_params:
            pH = batch_params['pH']
            threshold = batch_params['threshold']
            use_propka = batch_params['use_propka']
        else:
            # Ask for pH
            pH = prompt_with_context(
                processor=self.processor,
                prompt="Enter pH value for analysis",
                default="7.0",
                module="Protonation State Analyzer",
                description="Enter pH value for analysis"
            )
            try:
                pH = float(pH)
            except ValueError:
                self.processor.console.print(
                    "[red]Invalid pH value, using default of 7.0[/red]"
                )
                pH = 7.0

            # Ask about protonation threshold with explanation
            self.processor.console.print("\n[bold cyan]Protonation Probability Threshold[/bold cyan]")
            self.processor.console.print(
                "[grey50]Protonation states exist as equilibrium populations, not binary on/off switches.[/grey50]"
            )
            self.processor.console.print(
                f"[grey50]For example, at pH {pH:.1f}:[/grey50]"
            )
            self.processor.console.print(
                "[grey50]  • Residues far from their pKa (~99% one state) - threshold doesn't matter[/grey50]"
            )
            self.processor.console.print(
                "[grey50]  • Residues near their pKa (e.g., HIS) - significant populations of both states[/grey50]"
            )
            self.processor.console.print(
                "\n[grey50]The threshold determines the decision boundary for assigning discrete states:[/grey50]"
            )
            self.processor.console.print(
                "[grey50]  • 0.5 (default) - Assign the more likely state[/grey50]"
            )
            self.processor.console.print(
                "[grey50]  • >0.5 (conservative) - Only assign protonated when very confident[/grey50]"
            )
            self.processor.console.print(
                "[grey50]  • <0.5 (liberal) - Assign protonated even when less certain[/grey50]"
            )
            self.processor.console.print(
                "\n[yellow]Note: Standard MD requires discrete states. For residues near pKa, consider constant pH MD.[/yellow]\n"
            )

            threshold = prompt_with_context(
                processor=self.processor,
                prompt="Enter protonation probability threshold (0-1)",
                default="0.5",
                module="Protonation State Analyzer",
                description="Enter protonation threshold"
            )
            try:
                threshold = float(threshold)
                if threshold < 0 or threshold > 1:
                    self.processor.console.print(
                        "[red]Invalid threshold, using default of 0.5[/red]"
                    )
                    threshold = 0.5
            except ValueError:
                self.processor.console.print(
                    "[red]Invalid threshold, using default of 0.5[/red]"
                )
                threshold = 0.5

            # Ask if user wants to run PROPKA BEFORE setup
            use_propka = confirm_with_context(
                processor=self.processor,
                prompt="Run PROPKA for improved pKa prediction?",
                default=True,
                module="Protonation State Analyzer",
                description="Run PROPKA for pKa prediction"
            )

        # Setup the analyzer
        with self.processor.console.status(
            "[bold green]Setting up protonation state analyzer..."
        ):
            input_file = (
                input_file_path
                if using_file_path
                else self.get_from_workspace("pdb_file")
            )

            self.analyzer.setup(
                structure=None if using_file_path else input_structure,
                input_file=input_file,
                excluded_residues=excluded_residues,
                selected_chains=selected_chains,
                selected_residues=selected_residues,
                pH=pH,
                threshold=threshold,
            )

        # Store textbook pKa values for comparison (before PROPKA)
        textbook_pkas = {}
        if use_propka and hasattr(self.analyzer, 'titratable_residues'):
            for key, residue in self.analyzer.titratable_residues.items():
                textbook_pkas[key] = residue.get('pKa', 0.0)

        # Run PROPKA if requested
        if use_propka:
            with self.processor.console.status(
                "[bold green]Running PROPKA for pKa prediction..."
            ):
                propka_success = self.analyzer.run_propka()

            if not propka_success:
                console.print("[yellow]PROPKA failed. Using textbook pKa values.[/yellow]")
                use_propka = False  # Fall back to textbook values

        # Analyze protonation states (with PROPKA or textbook pKa values)
        with self.processor.console.status(
            "[bold green]Analyzing protonation states..."
        ):
            results, net_charge = self.analyzer.analyze_protonation_states()
            self.analysis_results = results

        # Display summary with method indicator and optional delta column
        method = "PROPKA" if use_propka else "Textbook"
        self._display_summary(results, net_charge, method=method, textbook_pkas=textbook_pkas if use_propka else None)

        # Store results in workspace - handle both microstate and default structure cases
        if is_microstate:
            # This is a microstate - store separately
            from pathlib import Path
            microstate_filename = Path(input_file_path).name

            # Get or create the microstate results dictionary
            microstate_results = self.get_from_workspace("microstate_protonation_results", {})

            # Store this microstate's results
            microstate_results[input_file_path] = {
                'filename': microstate_filename,
                'filepath': input_file_path,
                'results': results,
                'net_charge': net_charge,
                'pH': pH,
                'threshold': threshold,
                'method': method,
                'textbook_pkas': textbook_pkas if use_propka else {}
            }

            self.update_workspace("microstate_protonation_results", microstate_results)
            console.print(f"[grey50]Stored results for microstate: {microstate_filename}[/grey50]")
        else:
            # Default structure - store as before
            self.update_workspace("protonation_results", results)
            self.update_workspace("net_charge", net_charge)
            self.update_workspace("protonation_method", method)
            self.update_workspace("textbook_pkas", textbook_pkas if use_propka else {})

        # pH titration profile
        if batch_params:
            if batch_params.get('generate_ph_profile', False):
                self._generate_ph_profile()
        elif ask_ph_profile:
            if confirm_with_context(
                processor=self.processor,
                prompt="Would you like to generate a pH titration profile?",
                default=False,
                module="Protonation State Analyzer",
                description="Generate pH titration profile"
            ):
                self._generate_ph_profile()

    def _find_residue_name(self, structure, chain_id, resid):
        """Find the residue name for a given chain and residue ID"""
        for model in structure:
            if chain_id in model:
                for residue in model[chain_id]:
                    if residue.id[1] == resid:
                        return residue.resname
        return None

    def _find_residue_name_from_file(self, file_path, chain_id, resid):
        """Find residue name from a PDB file directly"""
        try:
            with open(file_path, "r") as f:
                for line in f:
                    if line.startswith(("ATOM", "HETATM")):
                        line_chain = line[21]
                        line_resid = int(line[22:26])

                        if line_chain == chain_id and line_resid == resid:
                            resname = line[17:20].strip()
                            return resname
            return None
        except Exception as e:
            self.processor.console.print(
                f"[yellow]Error reading PDB file: {str(e)}[/yellow]"
            )
            return None

    def _run_propka(self):
        """Run PROPKA to predict pKa values"""
        # Check if analyzer is initialized
        if (
            not hasattr(self.analyzer, "titratable_residues")
            or not self.analyzer.titratable_residues
        ):
            self.processor.console.print(
                "[yellow]Analyzer not initialized. Run analysis first.[/yellow]"
            )
            if confirm_with_context(
                processor=self.processor,
                prompt="Run analysis now?",
                module="Protonation State Analyzer",
                description="Run analysis now"
            ):
                self._analyze_protonation_states()
            else:
                return

        # Check for input file from analyzer
        input_file = getattr(self.analyzer, "input_file", None)
        if not input_file or not os.path.exists(input_file):
            self.processor.console.print(
                "[yellow]No valid input PDB file available for PROPKA.[/yellow]"
            )
            return

        # Get selected chains
        selected_chains = []
        filter_selections = self.get_from_workspace("filter_selections", {})
        if filter_selections:
            selected_chains = list(filter_selections.keys())

        # Check if we have current pH and threshold values in the analyzer
        if not hasattr(self.analyzer, "pH"):
            # Get pH value from workspace or use default
            pH = self.get_from_workspace("default_pH", 7.0)
            self.analyzer.pH = pH

        if not hasattr(self.analyzer, "threshold"):
            self.analyzer.threshold = 0.5

        # Run PROPKA
        with self.processor.console.status(
            "[bold green]Running PROPKA for pKa prediction..."
        ):
            success = self.analyzer.run_propka()

        if success:
            # Re-analyze with PROPKA results
            results, net_charge = self.analyzer.analyze_protonation_states()
            self.analysis_results = results

            # Display updated results
            self.processor.console.print(
                "[green]PROPKA completed successfully. Updated results:[/green]"
            )

            # Get textbook pKa values if available for delta display
            textbook_pkas = self.get_from_workspace("textbook_pkas", {})
            self._display_summary(results, net_charge, method="PROPKA", textbook_pkas=textbook_pkas)

            # Update workspace
            self.update_workspace("protonation_results", results)
            self.update_workspace("net_charge", net_charge)
            self.update_workspace("protonation_method", "PROPKA")

            # Store PROPKA results in workspace for future reference
            if hasattr(self.analyzer, "custom_pkas") and self.analyzer.custom_pkas:
                self.update_workspace("propka_pka_values", self.analyzer.custom_pkas)

            # Return success status for chaining
            return True
        else:
            self.processor.console.print(
                "[yellow]PROPKA failed. Using default pKa values.[/yellow]"
            )

            # If we have problematic residues, store them in workspace
            if (
                hasattr(self.analyzer, "problematic_residues")
                and self.analyzer.problematic_residues
            ):
                self.update_workspace(
                    "propka_problematic_residues", self.analyzer.problematic_residues
                )

            # Return failure status
            return False

    def _run_pdb2pqr_workflow(self):
        """
        Run pdb2pqr to add hydrogens with AMBER naming conventions.

        This workflow:
        1. Gets the best available structure (prefers transformed_pdb_file)
        2. Asks for target pH (default 7.0)
        3. Asks about constant pH titratable residues (AS4, GL4, HIP)
        4. Runs pdb2pqr with AMBER forcefield
        5. Post-processes for titratable residue renaming if needed
        6. Saves output to workspace
        """
        import subprocess
        import shutil
        from proprep.utils.structure_selector import StructureSelector

        console = self.processor.console
        workspace = self.get_workspace()

        # Get the best available structure (with key for display)
        structure_selector = StructureSelector(workspace, console)
        result = structure_selector.get_structure(silent=True, return_key=True)

        if not result:
            console.print("[red]No structure file available for pdb2pqr.[/red]")
            return False

        input_file, structure_key = result
        if not os.path.exists(input_file):
            console.print("[red]Structure file not found on disk.[/red]")
            return False

        # Show what structure we're using
        console.print(f"\n[bold]Running pdb2pqr on:[/bold] {os.path.basename(input_file)}")
        console.print(f"[grey50]Source: {structure_key}[/grey50]\n")

        # Get pH value
        default_ph = self.get_from_workspace("default_pH", 7.0)
        ph_input = prompt_with_context(
            processor=self.processor,
            prompt=f"Enter target pH for protonation",
            default=str(default_ph),
            module="Protonation State Analyzer",
            description="pdb2pqr target pH"
        )
        try:
            target_ph = float(ph_input)
        except ValueError:
            target_ph = default_ph
            console.print(f"[yellow]Invalid pH, using default: {target_ph}[/yellow]")

        # Ask about constant pH titratable residues
        use_titratable = confirm_with_context(
            processor=self.processor,
            prompt="Will you be running constant pH MD simulations?",
            default=False,
            module="Protonation State Analyzer",
            description="Enable constant pH MD residue naming (AS4, GL4, HIP)"
        )

        titratable_residues = {}
        if use_titratable:
            console.print("\n[bold]Titratable Residue Selection[/bold]")
            console.print("[grey50]Select which residues should be titratable (AS4, GL4, HIP).[/grey50]")
            console.print("[grey50]For pKa-shifted residues, they need protons added.[/grey50]")
            console.print("[grey50]tLEaP will add missing protons from lib files.[/grey50]\n")

            # Let user specify titratable residues by type
            for restype, (titname, desc) in [
                ("ASP", ("AS4", "Aspartate -> AS4 (titratable)")),
                ("GLU", ("GL4", "Glutamate -> GL4 (titratable)")),
                ("HIS", ("HIP", "Histidine -> HIP (titratable, doubly protonated)")),
            ]:
                if confirm_with_context(
                    processor=self.processor,
                    prompt=f"Make all {restype} residues titratable ({titname})?",
                    default=False,
                    module="Protonation State Analyzer",
                    description=desc
                ):
                    titratable_residues[restype] = titname

        # Get coordinating residues from RedoxSites (these should NOT be titrated)
        coordinating_residues = self._get_redoxsite_coordinating_residues()
        if coordinating_residues:
            console.print(f"\n[yellow]Found {len(coordinating_residues)} residues coordinating to RedoxSites[/yellow]")
            console.print("[grey50]These will be excluded from pKa-based protonation.[/grey50]")
            for chain, resid, resname in sorted(coordinating_residues):
                console.print(f"  [grey50]• {chain}:{resname}{resid}[/grey50]")

        # Setup output paths
        input_dir = os.path.dirname(input_file)
        input_base = os.path.splitext(os.path.basename(input_file))[0]

        pdb2pqr_dir = os.path.join(input_dir, "pdb2pqr_results")
        os.makedirs(pdb2pqr_dir, exist_ok=True)

        output_pqr = os.path.join(pdb2pqr_dir, f"{input_base}_pH{target_ph}.pqr")
        output_pdb = os.path.join(pdb2pqr_dir, f"{input_base}_pH{target_ph}_pdb2pqr.pdb")

        # Build pdb2pqr command
        # Note: We do NOT use --drop-water because coordinating waters must be preserved
        cmd = [
            "pdb2pqr",
            "--ff=AMBER",
            "--titration-state-method=propka",
            f"--with-ph={target_ph}",
            "--keep-chain",
            f"--pdb-output={output_pdb}",
        ]

        # If we have coordinating residues, get titratable residues and exclude them
        if coordinating_residues:
            titrate_only = self._get_titrate_only_residues(input_file, coordinating_residues)
            if titrate_only:
                # Format: "A:10,A:11,B:25"
                titrate_str = ",".join(f"{c}:{r}" for c, r, _ in titrate_only)
                cmd.append(f"--titrate_only={titrate_str}")
                console.print(f"[grey50]Titrating {len(titrate_only)} residues (excluding coordinating)[/grey50]")

        # Add input and output files
        cmd.extend([input_file, output_pqr])

        # Run pdb2pqr
        console.print(f"[bold green]Running pdb2pqr at pH {target_ph}...[/bold green]")
        logger.info(f"Running pdb2pqr: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout
            )

            if result.returncode != 0:
                console.print(f"[red]pdb2pqr failed with error:[/red]")
                console.print(f"[red]{result.stderr}[/red]")
                logger.error(f"pdb2pqr failed: {result.stderr}")
                return False

            # Check output files exist
            if not os.path.exists(output_pdb):
                console.print("[red]pdb2pqr did not create PDB output file.[/red]")
                return False

            console.print(f"[green]pdb2pqr completed successfully![/green]")

            # Post-process for titratable residues if needed
            if titratable_residues:
                console.print("\n[bold]Post-processing for titratable residues...[/bold]")
                final_pdb = self._postprocess_pdb2pqr_titratable(
                    output_pdb, titratable_residues, pdb2pqr_dir
                )
                if final_pdb:
                    output_pdb = final_pdb

            # Count atoms and hydrogens
            h_count = 0
            total_atoms = 0
            with open(output_pdb, 'r') as f:
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')):
                        total_atoms += 1
                        atom_name = line[12:16].strip()
                        if atom_name.startswith('H') or (len(atom_name) > 1 and atom_name[0].isdigit() and atom_name[1] == 'H'):
                            h_count += 1

            # Parse PQR file for net charge
            pqr_net_charge, pqr_atoms = self._parse_pqr_net_charge(output_pqr)

            console.print(f"\n[bold]Results:[/bold]")
            console.print(f"  Total atoms: {total_atoms}")
            console.print(f"  Hydrogen atoms added: {h_count}")
            if pqr_net_charge is not None:
                console.print(f"  Net charge: {pqr_net_charge:+.2f}")
            console.print(f"  Output PDB: {output_pdb}")
            console.print(f"  Output PQR: {output_pqr}")

            # Update workspace
            self.update_workspace("pdb2pqr_output_file", output_pdb)
            self.update_workspace("pdb2pqr_pqr_file", output_pqr)
            self.update_workspace("protonation_pdb_file", output_pdb)
            self.update_workspace("pdb2pqr_target_ph", target_ph)
            if pqr_net_charge is not None:
                self.update_workspace("net_charge", pqr_net_charge)
            if titratable_residues:
                self.update_workspace("constant_ph_residues", titratable_residues)

            console.print("\n[green]Structure with hydrogens saved to workspace.[/green]")
            console.print("[grey50]Use Topology Generator for accurate solvation estimation.[/grey50]")

            return True

        except subprocess.TimeoutExpired:
            console.print("[red]pdb2pqr timed out after 5 minutes.[/red]")
            return False
        except FileNotFoundError:
            console.print("[red]pdb2pqr command not found. Is it installed?[/red]")
            console.print("[grey50]Install with: pip install pdb2pqr[/grey50]")
            return False
        except Exception as e:
            console.print(f"[red]Error running pdb2pqr: {str(e)}[/red]")
            logger.error(f"pdb2pqr error: {str(e)}")
            return False

    def _postprocess_pdb2pqr_titratable(self, pdb_file, titratable_residues, output_dir):
        """
        Post-process pdb2pqr output to rename residues for constant pH MD.

        pdb2pqr outputs standard AMBER names (ASP, GLU, HID/HIE/HIP).
        For constant pH, we need to rename to titratable forms (AS4, GL4, HIP).
        tLEaP will add any missing protons from the lib files.

        Args:
            pdb_file: Path to pdb2pqr output PDB
            titratable_residues: Dict mapping residue type to titratable name
                e.g., {"ASP": "AS4", "GLU": "GL4", "HIS": "HIP"}
            output_dir: Directory for output

        Returns:
            Path to processed PDB file, or None on failure
        """
        try:
            basename = os.path.splitext(os.path.basename(pdb_file))[0]
            output_file = os.path.join(output_dir, f"{basename}_titratable.pdb")

            # Build residue name mapping
            # pdb2pqr uses AMBER names: ASP, GLU, HID, HIE, HIP
            rename_map = {}
            if "ASP" in titratable_residues:
                rename_map["ASP"] = titratable_residues["ASP"]  # ASP -> AS4
            if "GLU" in titratable_residues:
                rename_map["GLU"] = titratable_residues["GLU"]  # GLU -> GL4
            if "HIS" in titratable_residues:
                # All histidine forms -> HIP for constant pH
                for his_form in ["HID", "HIE", "HIS"]:
                    rename_map[his_form] = titratable_residues["HIS"]
                # HIP stays as HIP

            renamed_count = 0
            with open(pdb_file, 'r') as f_in, open(output_file, 'w') as f_out:
                for line in f_in:
                    if line.startswith(('ATOM', 'HETATM')):
                        resname = line[17:20].strip()
                        if resname in rename_map:
                            new_resname = rename_map[resname].ljust(3)
                            line = line[:17] + new_resname + line[20:]
                            renamed_count += 1
                    f_out.write(line)

            self.processor.console.print(
                f"[green]Renamed {renamed_count} atoms for constant pH MD.[/green]"
            )

            return output_file

        except Exception as e:
            logger.error(f"Error post-processing titratable residues: {str(e)}")
            self.processor.console.print(
                f"[yellow]Warning: Could not rename titratable residues: {str(e)}[/yellow]"
            )
            return None

    def _parse_pqr_net_charge(self, pqr_file):
        """
        Parse a PQR file and calculate the total net charge.

        PQR format: ATOM/HETATM record with charge in column 55-62 (or after coordinates).
        The charge field comes after x, y, z coordinates.

        Args:
            pqr_file: Path to PQR file from pdb2pqr

        Returns:
            Tuple of (net_charge: float, n_atoms: int) or (None, 0) on error
        """
        net_charge = 0.0
        n_atoms = 0

        try:
            with open(pqr_file, 'r') as f:
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')):
                        # PQR format varies but charge is typically after coordinates
                        # Format: ATOM serial name res chain resSeq x y z charge radius
                        parts = line.split()
                        if len(parts) >= 9:
                            try:
                                # Charge is second-to-last field, radius is last
                                charge = float(parts[-2])
                                net_charge += charge
                                n_atoms += 1
                            except (ValueError, IndexError):
                                # Try parsing from fixed columns if whitespace split fails
                                try:
                                    # Charge typically at columns 55-62
                                    charge = float(line[54:62].strip())
                                    net_charge += charge
                                    n_atoms += 1
                                except (ValueError, IndexError):
                                    pass

            return net_charge, n_atoms

        except Exception as e:
            logger.warning(f"Could not parse PQR file for net charge: {e}")
            return None, 0

    def _get_redoxsite_coordinating_residues(self):
        """
        Get residues that are part of RedoxSite objects (coordinating to metals).

        These residues should NOT have their protonation state changed by pdb2pqr
        because their protonation is determined by metal coordination.

        Returns:
            Set of (chain, resid, resname) tuples for coordinating residues
        """
        coordinating = set()

        detected_sites = self.get_from_workspace("detected_redox_sites", [])
        if not detected_sites:
            return coordinating

        for site in detected_sites:
            # Get residues from the site's residue_groups or atoms
            if hasattr(site, 'residue_groups') and site.residue_groups:
                for (chain, resid, insertion_code), atom_coords in site.residue_groups.items():
                    # Get resname from coord_to_pdb for any atom in this residue
                    if atom_coords and hasattr(site, 'coord_to_pdb'):
                        atom_info = site.coord_to_pdb.get(atom_coords[0])
                        if atom_info:
                            resname = atom_info.get('resname', '')
                            # Only include titratable residue types
                            if resname in ('ASP', 'GLU', 'HIS', 'HID', 'HIE', 'HIP', 'CYS', 'TYR', 'LYS'):
                                coordinating.add((chain, resid, resname))

            # Also check atoms directly
            if hasattr(site, 'atoms'):
                for atom in site.atoms:
                    if hasattr(atom, 'resname') and atom.resname in ('ASP', 'GLU', 'HIS', 'HID', 'HIE', 'HIP', 'CYS', 'TYR', 'LYS'):
                        coordinating.add((atom.chain, atom.resid, atom.resname))

        return coordinating

    def _get_titrate_only_residues(self, pdb_file, exclude_residues):
        """
        Get list of titratable residues in PDB, excluding specified residues.

        Scans the PDB for titratable residue types (ASP, GLU, HIS, etc.)
        and returns those not in the exclude list.

        Args:
            pdb_file: Path to input PDB file
            exclude_residues: Set of (chain, resid, resname) tuples to exclude

        Returns:
            List of (chain, resid, resname) tuples to titrate
        """
        titratable_types = {'ASP', 'GLU', 'HIS', 'HID', 'HIE', 'HIP', 'CYS', 'TYR', 'LYS'}

        # Scan PDB for titratable residues
        all_titratable = set()
        try:
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')):
                        chain = line[21].strip() or ' '
                        try:
                            resid = int(line[22:26].strip())
                        except ValueError:
                            continue
                        resname = line[17:20].strip()

                        if resname in titratable_types:
                            all_titratable.add((chain, resid, resname))

        except Exception as e:
            logger.error(f"Error scanning PDB for titratable residues: {e}")
            return []

        # Exclude coordinating residues
        exclude_set = {(c, r, n) for c, r, n in exclude_residues}

        # Also match by chain:resid only (resname might differ)
        exclude_by_id = {(c, r) for c, r, _ in exclude_residues}

        titrate_only = []
        for chain, resid, resname in sorted(all_titratable):
            if (chain, resid, resname) not in exclude_set and (chain, resid) not in exclude_by_id:
                titrate_only.append((chain, resid, resname))

        return titrate_only

    def _view_results(self):
        """View protonation analysis results"""
        console = self.processor.console

        # Check if we have microstate results
        microstate_results = self.get_from_workspace("microstate_protonation_results", {})

        if microstate_results and len(microstate_results) > 0:
            # We have microstate results - ask user what to view
            console.print(f"\n[bold cyan]Available Results:[/bold cyan]")
            console.print(f"  • Default structure results: {'Yes' if self.get_from_workspace('protonation_results') else 'No'}")
            console.print(f"  • Microstate results: {len(microstate_results)} microstate(s)")


            choices = ["1", "2"]
            options_map = {
                "1": "View default structure results",
                "2": "View microstate results"
            }

            choice = prompt_with_context(
                processor=self.processor,
                prompt="What would you like to view?",
                choices=choices,
                default="2",
                module="Protonation State Analyzer",
                description="Select results to view",
                options_map=options_map
            )

            if choice == "2":
                self._view_microstate_results()
                return

        # View default structure results
        if not self.analysis_results:
            results = self.get_from_workspace("protonation_results")
            if results:
                self.analysis_results = results
                net_charge = self.get_from_workspace("net_charge", 0.0)
                method = self.get_from_workspace("protonation_method", "Textbook")
                textbook_pkas = self.get_from_workspace("textbook_pkas", {})
                self._display_summary(results, net_charge, method=method, textbook_pkas=textbook_pkas)
            else:
                self.processor.console.print(
                    "[yellow]No analysis results available. Run analysis first.[/yellow]"
                )
                if confirm_with_context(
                    processor=self.processor,
                    prompt="Run analysis now?",
                    module="Protonation State Analyzer",
                    description="Run analysis now"
                ):
                    self._analyze_protonation_states()
                return
        else:
            # Display existing results
            net_charge = self.get_from_workspace("net_charge", 0.0)
            method = self.get_from_workspace("protonation_method", "Textbook")
            textbook_pkas = self.get_from_workspace("textbook_pkas", {})
            self._display_summary(self.analysis_results, net_charge, method=method, textbook_pkas=textbook_pkas)

    def _view_microstate_results(self):
        """View protonation analysis results for microstates"""
        console = self.processor.console
        microstate_results = self.get_from_workspace("microstate_protonation_results", {})

        if not microstate_results:
            console.print("[yellow]No microstate results available.[/yellow]")
            return

        # Display table of available microstates
        from rich.table import Table
        table = Table(title="Analyzed Microstates")
        table.add_column("#", style="cyan", justify="right")
        table.add_column("Filename", style="green")
        table.add_column("pH", style="yellow")
        table.add_column("Net Charge", style="magenta")

        microstate_list = list(microstate_results.items())
        for i, (filepath, data) in enumerate(microstate_list, 1):
            table.add_row(
                str(i),
                data['filename'],
                f"{data['pH']:.1f}",
                f"{data['net_charge']:.2f}"
            )

        console.print(table)

        # Ask user to select which ones to view

        console.print("\n[bold]Selection Options:[/bold]")
        console.print("  • Single number: e.g., '3'")
        console.print("  • Multiple: e.g., '1,3,5'")
        console.print("  • Range: e.g., '1-4'")
        console.print("  • All: enter 'all'")

        selection = prompt_with_context(
            processor=self.processor,
            prompt="Select microstate(s) to view",
            default="all",
            module="Protonation State Analyzer",
            description="Microstate selection for viewing"
        )

        # Parse selection
        selected_indices = self._parse_microstate_selection(selection, len(microstate_list))

        if not selected_indices:
            console.print("[yellow]No valid microstates selected.[/yellow]")
            return

        # Display results for each selected microstate
        for idx in selected_indices:
            filepath, data = microstate_list[idx]
            console.print(f"\n[bold cyan]{'='*60}[/bold cyan]")
            console.print(f"[bold]Microstate: {data['filename']}[/bold]")
            console.print(f"[bold cyan]{'='*60}[/bold cyan]\n")

            method = data.get('method', 'Textbook')
            textbook_pkas = data.get('textbook_pkas', {})
            self._display_summary(data['results'], data['net_charge'], method=method, textbook_pkas=textbook_pkas)

    def _select_analyzed_microstates(self):
        """
        Check for analyzed microstates and prompt user to select which ones to process.

        Returns:
            list of tuples (filepath, results_data) or None for default structure
        """
        console = self.processor.console
        workspace = self.get_workspace()

        # Check for microstate analysis results
        microstate_results = workspace.get("microstate_protonation_results", {})

        if not microstate_results or len(microstate_results) == 0:
            # No microstate results
            return None

        # Microstates analyzed - ask user if they want to process them
        console.print(f"\n[bold cyan]Detected analysis results for {len(microstate_results)} microstate(s)[/bold cyan]")

        use_microstates = confirm_with_context(
            processor=self.processor,
            prompt="Set MD residue names for microstate(s) instead of default structure?",
            default=True,
            module="Protonation State Analyzer",
            description="Process microstates"
        )

        if not use_microstates:
            return None

        # Display available microstates
        from rich.table import Table
        table = Table(title="Analyzed Microstates")
        table.add_column("#", style="cyan", justify="right")
        table.add_column("Filename", style="green")
        table.add_column("pH", style="yellow")
        table.add_column("Net Charge", style="magenta")

        microstate_list = list(microstate_results.items())
        for i, (filepath, data) in enumerate(microstate_list, 1):
            table.add_row(
                str(i),
                data['filename'],
                f"{data['pH']:.1f}",
                f"{data['net_charge']:.2f}"
            )

        console.print(table)

        # Ask user to select one or more microstates

        console.print("\n[bold]Selection Options:[/bold]")
        console.print("  • Single number: e.g., '3'")
        console.print("  • Multiple: e.g., '1,3,5'")
        console.print("  • Range: e.g., '1-4'")
        console.print("  • All: enter 'all'")

        selection = prompt_with_context(
            processor=self.processor,
            prompt="Select microstate(s) to process",
            default="all",
            module="Protonation State Analyzer",
            description="Microstate selection for MD names"
        )

        # Parse selection (reuse the helper method)
        selected_indices = self._parse_microstate_selection(selection, len(microstate_list))

        if not selected_indices:
            console.print("[yellow]No valid microstates selected.[/yellow]")
            return None

        # Get the selected microstate data
        selected_microstates = [(microstate_list[i][0], microstate_list[i][1]) for i in selected_indices]

        console.print(f"\n[green]Selected {len(selected_microstates)} microstate(s) for MD residue name assignment[/green]")

        return selected_microstates

    def _set_residue_names_for_microstates(self, selected_microstates):
        """Set residue names for microstates using a compact comparison grid.

        Shows a grid comparing default recommendations across all microstates,
        highlighting only residues whose recommendations differ. The user can
        accept all defaults or override specific cells before applying.
        """
        from pathlib import Path

        console = self.processor.console
        num_microstates = len(selected_microstates)

        # Ask once if this is for constant pH simulation
        constant_pH = confirm_with_context(
            processor=self.processor,
            prompt="Is this for constant pH simulation?",
            default=False,
            module="Protonation State Analyzer",
            description="Set for constant pH simulation"
        )

        # Standard fixed states are mutually exclusive with constant pH (which
        # wants titratable names), so only offer this when not doing CpHMD.
        standard_states = False
        if not constant_pH:
            standard_states = confirm_with_context(
                processor=self.processor,
                prompt="Set all residues to standard protonation states "
                       "(ASP/GLU deprotonated, CYS/LYS/ARG/TYR protonated, "
                       "HIS->HIE, heme PRN->PRD)?",
                default=False,
                module="Protonation State Analyzer",
                description="Use standard protonation states (ignore PROPKA)"
            )

        # Phase A & B: Build the grid of recommendations across all microstates
        with console.status("[bold green]Computing recommendations for all microstates..."):
            grid, uniform, varying = self._build_protonation_grid(
                selected_microstates, constant_pH, standard_states
            )

        if not grid:
            console.print("[yellow]No titratable residues found in analysis results.[/yellow]")
            return

        # Phase C: Display
        self._display_uniform_summary(uniform)

        # Build microstate labels (short sequential numbers)
        microstate_names = [f"{i:03d}" for i in range(1, num_microstates + 1)]

        if varying:
            console.print(f"\n[bold yellow]Residues that DIFFER across microstates: "
                          f"{len(varying)}[/bold yellow]")
            self._display_protonation_grid(varying, microstate_names)
            self._display_uniform_summary(uniform)
        else:
            console.print(f"\n[green]All {len(uniform)} residues have identical "
                          f"recommendations across all {num_microstates} microstates.[/green]")
            self._display_uniform_detail(uniform)

        # Phase D: Override prompt
        console.print("\n[bold]Override syntax:[/bold]")
        console.print("  [cyan]A:94:ASH[/cyan]              "
                      "- Set across ALL microstates")
        if varying:
            console.print("  [cyan]5:A:94:ASH[/cyan]            "
                          "- Set for microstate 5 only")
            console.print("  [cyan]3-7:A:94:ASH[/cyan]          "
                          "- Set for microstate range 3-7")

        user_input = prompt_with_context(
            processor=self.processor,
            prompt="Enter overrides (space-separated), or press Enter to accept all",
            default="",
            module="Protonation State Analyzer",
            description="Override protonation recommendations"
        ).strip()

        if user_input:
            overrides = self._parse_protonation_overrides(
                user_input, grid, num_microstates
            )
            # Apply overrides to the grid
            for (key, ms_idx), name in overrides.items():
                grid[key][ms_idx] = name
            if overrides:
                console.print(f"[green]Applied {len(overrides)} override(s)[/green]")

        # Phase E: Apply grid to PDB files
        from rich.progress import Progress, TextColumn, BarColumn, MofNCompleteColumn, TimeElapsedColumn

        with Progress(
            TextColumn("[bold green]{task.description}"),
            BarColumn(),
            MofNCompleteColumn(),
            TimeElapsedColumn(),
            console=console,
        ) as progress:
            task = progress.add_task("Updating microstates...", total=num_microstates)

            for i, (filepath, microstate_data) in enumerate(selected_microstates):
                # Extract this microstate's column from the grid
                md_residue_names = {key: names[i] for key, names in grid.items()}

                # Update the PDB file
                base_name = Path(filepath).stem
                output_file = f"{base_name}_md_names.pdb"

                saved_file = self._update_structure_in_file(
                    filepath, md_residue_names, output_file
                )

                if saved_file:
                    microstate_data['updated_filepath'] = saved_file
                    microstate_data['md_residue_names'] = md_residue_names
                else:
                    console.print(f"[red]Failed to update: {Path(filepath).name}[/red]")

                progress.advance(task)

        # Update workspace
        workspace = self.get_workspace()
        microstate_results = workspace.get("microstate_protonation_results", {})
        for filepath, microstate_data in selected_microstates:
            microstate_results[filepath] = microstate_data
        self.update_workspace("microstate_protonation_results", microstate_results)

        console.print(f"\n[bold green]✓ Updated {num_microstates} microstate(s) "
                      f"with {len(grid)} residue name assignments each[/bold green]")

    # --- Protonation Grid helpers ---

    # Numeric protonation level encoding for compact grid display
    _NAME_TO_SYMBOL = {
        'ASP': '0', 'ASH': '1', 'AS4': '3',
        'GLU': '0', 'GLH': '1', 'GL4': '3',
        'HIE': '1', 'HID': '-1', 'HIP': '2',
        'LYS': '1', 'LYN': '0',
        'CYS': '1', 'CYM': '0',
        'TYR': '1', 'TYM': '0',
        'PRN': '1', 'PRP': '1', 'PRD': '0',
    }

    def _build_protonation_grid(self, selected_microstates, constant_pH,
                                standard_states=False):
        """Build a grid of default recommendations across all microstates.

        Args:
            standard_states: If True, assign canonical fixed protonation states
                instead of PROPKA-based recommendations (see
                protonation_worker.get_default_recommendations).

        Returns:
            grid: Dict[key] -> List[str]  (one AMBER name per microstate)
            uniform: Dict[key] -> str     (residues same across all microstates)
            varying: Dict[key] -> List[str]  (residues that differ)
        """
        all_recommendations = []

        for filepath, ms_data in selected_microstates:
            # Swap worker results to this microstate's analysis
            original_results = self.analyzer.results
            self.analyzer.results = ms_data['results']

            recs = self.analyzer.get_default_recommendations(
                constant_ph_simulation=constant_pH,
                standard_states=standard_states
            )
            all_recommendations.append(recs)

            # Restore
            self.analyzer.results = original_results

        # Union of all residue keys
        all_keys = set()
        for recs in all_recommendations:
            all_keys.update(recs.keys())

        # Build grid and partition
        grid = {}
        uniform = {}
        varying = {}

        for key in sorted(all_keys):
            names = [recs.get(key, '?') for recs in all_recommendations]
            grid[key] = names
            if len(set(names)) == 1:
                uniform[key] = names[0]
            else:
                varying[key] = names

        return grid, uniform, varying

    def _display_uniform_summary(self, uniform):
        """Display a compact summary of residues with identical recommendations."""
        console = self.processor.console

        if not uniform:
            return

        # Group by AMBER name
        by_name = {}
        for key, name in uniform.items():
            by_name.setdefault(name, []).append(key)

        # Build compact summary string
        parts = []
        for name in sorted(by_name.keys()):
            parts.append(f"{name}: {len(by_name[name])}")

        console.print(f"\n[green]Residues with identical recommendations across all microstates: "
                      f"{len(uniform)}[/green]")
        console.print(f"  [grey50]{', '.join(parts)}[/grey50]")

    def _display_uniform_detail(self, uniform):
        """Display uniform residues grouped by AMBER name with residue identities.

        Used when no residues vary across microstates, so the user can still
        see exactly which residues have which assignment for overrides.
        """
        console = self.processor.console

        if not uniform:
            return

        # Group by AMBER name, with formatted residue labels
        by_name = {}
        for key, name in uniform.items():
            parts = key.split("_")
            if len(parts) >= 3:
                chain = parts[0]
                resid = parts[1]
                label = f"{chain}:{resid}"
            else:
                label = key
            by_name.setdefault(name, []).append(label)

        for name in sorted(by_name.keys()):
            residues = sorted(by_name[name])
            residue_str = ", ".join(residues)
            # Wrap long lines
            if len(residue_str) > 80:
                residue_str = residue_str[:77] + "..."
            console.print(f"  [cyan]{name:>3}[/cyan] ({len(residues):>3}): [grey50]{residue_str}[/grey50]")

    def _display_protonation_grid(self, varying, microstate_names):
        """Display varying residues as a compact grid with auto-pagination and auto-transpose.

        Args:
            varying: Dict[key] -> List[str] of AMBER names per microstate
            microstate_names: List[str] of short microstate labels (e.g. ['001', '002', ...])
        """
        import shutil

        console = self.processor.console
        n_residues = len(varying)
        n_microstates = len(microstate_names)

        if n_residues == 0:
            return

        # Decide orientation
        transposed = n_microstates > n_residues

        term_width = shutil.get_terminal_size().columns
        label_width = 16
        col_width = 4
        cols_per_page = max(1, (term_width - label_width) // col_width)

        # Build residue labels: "A:ASP  94"
        residue_labels = {}
        for key in varying:
            parts = key.split("_")
            if len(parts) >= 3:
                chain = parts[0]
                resid = parts[1]
                res_type = "_".join(parts[2:])
                residue_labels[key] = f"{chain}:{res_type}{resid}"
            else:
                residue_labels[key] = key

        if transposed:
            self._display_grid_transposed(varying, microstate_names, residue_labels, cols_per_page)
        else:
            self._display_grid_normal(varying, microstate_names, residue_labels, cols_per_page)

        # Legend
        console.print("\n[grey50]Legend: 0=deprotonated  1=protonated  "
                      "-1=HID(Nδ)  2=HIP(both)  3=titratable(CpHMD)[/grey50]")

    def _display_grid_normal(self, varying, microstate_names, residue_labels, cols_per_page):
        """Display grid with residues as rows, microstates as columns."""
        console = self.processor.console
        n_microstates = len(microstate_names)
        n_pages = (n_microstates + cols_per_page - 1) // cols_per_page

        keys = list(varying.keys())

        for page in range(n_pages):
            start = page * cols_per_page
            end = min(start + cols_per_page, n_microstates)
            page_ms_names = microstate_names[start:end]

            title = "Residues with varying recommendations across microstates"
            if n_pages > 1:
                title += f" (page {page + 1}/{n_pages}, microstates {start + 1}-{end})"

            table = Table(title=title)
            table.add_column("Residue", style="cyan", no_wrap=True)
            for ms_name in page_ms_names:
                table.add_column(ms_name, justify="center", width=3)

            for key in keys:
                names = varying[key]
                page_names = names[start:end]
                symbols = [self._NAME_TO_SYMBOL.get(n, n[:2]) for n in page_names]
                table.add_row(residue_labels[key], *symbols)

            console.print(table)

    def _display_grid_transposed(self, varying, microstate_names, residue_labels, cols_per_page):
        """Display grid with microstates as rows, residues as columns (when fewer residues)."""
        console = self.processor.console
        keys = list(varying.keys())
        n_residues = len(keys)
        n_pages = (n_residues + cols_per_page - 1) // cols_per_page

        for page in range(n_pages):
            start = page * cols_per_page
            end = min(start + cols_per_page, n_residues)
            page_keys = keys[start:end]

            title = "Residues with varying recommendations across microstates (transposed)"
            if n_pages > 1:
                title += f" (page {page + 1}/{n_pages})"

            table = Table(title=title)
            table.add_column("Microstate", style="cyan", no_wrap=True)
            for key in page_keys:
                table.add_column(residue_labels[key], justify="center", width=4)

            for ms_idx, ms_name in enumerate(microstate_names):
                symbols = []
                for key in page_keys:
                    name = varying[key][ms_idx]
                    symbols.append(self._NAME_TO_SYMBOL.get(name, name[:2]))
                table.add_row(ms_name, *symbols)

            console.print(table)

    def _parse_protonation_overrides(self, user_input, grid, num_microstates):
        """Parse override syntax into grid modifications.

        Input format: space-separated tokens
          A:94:ASH             -> all microstates, chain A residue 94 -> ASH
          5:A:94:ASH           -> microstate 5 only
          3-7:A:94:ASH         -> microstates 3-7

        Returns:
            Dict[(key, ms_idx), name]
        """
        console = self.processor.console
        overrides = {}

        for token in user_input.strip().split():
            parts = token.split(':')
            if len(parts) == 3:
                # chain:resid:name -> all microstates
                chain, resid, name = parts
                try:
                    key = self._find_residue_key(chain, int(resid), grid)
                except ValueError:
                    console.print(f"[yellow]Invalid residue ID '{resid}' in '{token}', skipping[/yellow]")
                    continue
                if key:
                    for ms_idx in range(num_microstates):
                        overrides[(key, ms_idx)] = name.upper()
                else:
                    console.print(f"[yellow]Residue {chain}:{resid} not found in grid, skipping[/yellow]")
            elif len(parts) == 4:
                # ms_spec:chain:resid:name
                ms_spec, chain, resid, name = parts
                try:
                    key = self._find_residue_key(chain, int(resid), grid)
                except ValueError:
                    console.print(f"[yellow]Invalid residue ID '{resid}' in '{token}', skipping[/yellow]")
                    continue
                if key:
                    for ms_idx in self._parse_ms_range(ms_spec, num_microstates):
                        overrides[(key, ms_idx)] = name.upper()
                else:
                    console.print(f"[yellow]Residue {chain}:{resid} not found in grid, skipping[/yellow]")
            else:
                console.print(f"[yellow]Invalid override syntax '{token}', expected "
                              f"chain:resid:name or ms:chain:resid:name[/yellow]")

        return overrides

    def _find_residue_key(self, chain, resid, grid):
        """Find the matching key in the grid for a (chain, resid) pair.

        Keys are in format 'A_94_HIS'. Returns the key or None.
        """
        for key in grid:
            parts = key.split("_")
            if len(parts) >= 3:
                key_chain = parts[0]
                try:
                    key_resid = int(parts[1])
                except ValueError:
                    continue
                if key_chain == chain and key_resid == resid:
                    return key
        return None

    def _parse_ms_range(self, ms_spec, num_microstates):
        """Parse microstate range spec (1-indexed) into 0-indexed indices.

        '5'   -> [4]
        '3-7' -> [2, 3, 4, 5, 6]
        """
        indices = []
        try:
            if '-' in ms_spec:
                start_str, end_str = ms_spec.split('-', 1)
                start = int(start_str) - 1  # Convert to 0-indexed
                end = int(end_str)  # end is inclusive in user syntax
                for i in range(max(0, start), min(end, num_microstates)):
                    indices.append(i)
            else:
                idx = int(ms_spec) - 1  # Convert to 0-indexed
                if 0 <= idx < num_microstates:
                    indices.append(idx)
        except ValueError:
            pass
        return indices

    def _set_residue_names(self):
        """Set residue names for MD simulation"""
        # Check if we have analyzed microstates to process
        selected_microstates = self._select_analyzed_microstates()

        if selected_microstates:
            # Process each selected microstate
            self._set_residue_names_for_microstates(selected_microstates)
            return

        # No microstates - process default structure
        # Check if results are available
        if not self.analysis_results:
            results = self.get_from_workspace("protonation_results")
            if results:
                self.analysis_results = results
            else:
                self.processor.console.print(
                    "[yellow]No analysis results available. Run analysis first.[/yellow]"
                )
                if confirm_with_context(
                    processor=self.processor,
                    prompt="Run analysis now?",
                    module="Protonation State Analyzer",
                    description="Run analysis now"
                ):
                    self._analyze_protonation_states()
                return

        # Ask if this is for constant pH simulation
        constant_pH = confirm_with_context(
            processor=self.processor,
            prompt="Is this for constant pH simulation?",
            default=False,
            module="Protonation State Analyzer",
            description="Set for constant pH simulation"
        )

        # Standard fixed states are mutually exclusive with constant pH (which
        # wants titratable names), so only offer this when not doing CpHMD.
        standard_states = False
        if not constant_pH:
            standard_states = confirm_with_context(
                processor=self.processor,
                prompt="Set all residues to standard protonation states "
                       "(ASP/GLU deprotonated, CYS/LYS/ARG/TYR protonated, "
                       "HIS->HIE, heme PRN->PRD)?",
                default=False,
                module="Protonation State Analyzer",
                description="Use standard protonation states (ignore PROPKA)"
            )

        if standard_states:
            # Canonical states ignore the per-residue analysis, so the
            # interactive review would have nothing to decide; assign directly.
            self.md_residue_names = self.analyzer.get_default_recommendations(
                standard_states=True
            )
            self.processor.console.print(
                "[cyan]Assigned standard protonation states "
                "(PROPKA recommendations ignored).[/cyan]"
            )
        else:
            # Set residue names - WITHOUT the status indicator since this is interactive
            self.md_residue_names = self.analyzer.set_md_residue_names(
                interactive=True, constant_ph_simulation=constant_pH
            )

        # Store in workspace
        self.update_workspace("md_residue_names", self.md_residue_names)
        
        # Store constant pH data if available
        if constant_pH and hasattr(self.analyzer, 'constant_ph_data'):
            self.update_workspace("constant_ph_data", self.analyzer.constant_ph_data)
            self.processor.console.print("[cyan]Constant pH titration data saved to workspace[/cyan]")

        self.processor.console.print(
            f"[green]Set names for {len(self.md_residue_names)} residues[/green]"
        )

        # Ask if user wants to update the structure
        if confirm_with_context(
            processor=self.processor,
            prompt="Update structure with these residue names?",
            default=True,
            module="Protonation State Analyzer",
            description="Update structure with residue names"
        ):
            self._update_structure()

    def _update_structure(self):
        """Update structure with MD residue names - uses direct file manipulation to preserve modified residues"""
        # Check if residue names are available
        if not self.md_residue_names:
            names = self.get_from_workspace("md_residue_names")
            if names:
                self.md_residue_names = names
            else:
                self.processor.console.print(
                    "[yellow]No residue names defined. Set residue names first.[/yellow]"
                )
                if confirm_with_context(
                    processor=self.processor,
                    prompt="Set residue names now?",
                    module="Protonation State Analyzer",
                    description="Set residue names now"
                ):
                    self._set_residue_names()
                return

        # The residue keys in md_residue_names use chain IDs and residue
        # numbers from the file that was analyzed.  We MUST apply the update
        # to that same file; otherwise chain/residue numbering may differ
        # (e.g. MODELLER global numbering vs. original per-chain numbering).
        analysis_file = getattr(self.analyzer, "input_file", None)

        priority_keys = [
            ("protonation_pdb_file", "protonation-updated"),
            ("transformed_pdb_file", "transformed"),
            ("repaired_pdb_file", "repaired"),
            ("filtered_pdb_file", "filtered"),
            ("original_pdb_file", "original"),
            ("pdb_file", "loaded")  # Fallback
        ]

        input_file_path = None
        structure_source = None

        # First choice: the file the analysis was actually run on
        if analysis_file and os.path.exists(str(analysis_file)):
            input_file_path = str(analysis_file)
            structure_source = "analyzed"
            self.processor.console.print(
                f"[cyan]Using analyzed structure for MD residue name update: {analysis_file}[/cyan]"
            )
        else:
            # Fallback: try each priority level
            for key, source_name in priority_keys:
                structure_file = self.get_from_workspace(key)
                if structure_file and os.path.exists(str(structure_file)):
                    input_file_path = str(structure_file)
                    structure_source = source_name
                    self.processor.console.print(
                        f"[cyan]Using {source_name} structure for MD residue name update: {structure_file}[/cyan]"
                    )
                    break
        
        if input_file_path is None:
            self.processor.console.print(
                "[yellow]No structure file available. Please ensure a structure has been loaded.[/yellow]"
            )
            return

        # IMPORTANT: Always use direct file manipulation to preserve modified residues
        # Never use BioPython after redox site preparation
        
        # Ask for output filename before starting the spinner
        output_file = prompt_with_context(
            processor=self.processor,
            prompt="Enter output filename for updated structure",
            default="structure_with_md_names.pdb",
            module="Protonation State Analyzer",
            description="Enter output filename"
        )
        
        with self.processor.console.status(
            "[bold green]Updating structure with MD residue names..."
        ):
            saved_file = self._update_structure_in_file(
                input_file_path, self.md_residue_names, output_file
            )
        
        if saved_file:
            self.processor.console.print(
                f"[green]Saved updated structure to {saved_file}[/green]"
            )

            # Store the file path in workspace (NOT a BioPython object)
            # Using standardized key name: protonation_pdb_file
            self.update_workspace("protonation_pdb_file", saved_file)

            # CRITICAL: Update protonation_results with the new residue names
            # so that CPIN generation gets the correct names (AS4, GL4, HIP)
            protonation_results = self.get_from_workspace("protonation_results", {})
            if protonation_results and self.md_residue_names:
                updated_count = 0
                for key in protonation_results.keys():
                    if key in self.md_residue_names:
                        old_name = protonation_results[key].get("type", "")
                        new_name = self.md_residue_names[key]
                        if old_name != new_name:
                            protonation_results[key]["type"] = new_name
                            updated_count += 1

                if updated_count > 0:
                    self.update_workspace("protonation_results", protonation_results)
                    self.processor.console.print(
                        f"[green]✓ Updated {updated_count} residue names in protonation_results[/green]"
                    )

            # Note to user about the update
            if structure_source == "repaired":
                self.processor.console.print(
                    "[green]Updated the repaired structure with MD residue names[/green]"
                )
            elif structure_source == "filtered":
                self.processor.console.print(
                    "[yellow]The filtered structure was updated because no repaired structure was found.[/yellow]"
                )

    def _update_structure_in_file(self, input_file, md_residue_names, output_file):
        """
        Update residue names in a PDB file directly.

        Args:
            input_file: Path to input PDB file
            md_residue_names: Dictionary mapping residue keys to MD residue names
            output_file: Path to save the updated structure

        Returns:
            Path to the saved file, or None if update failed
        """
        try:
            # Parse residue keys into a more efficient lookup format
            residue_map = {}
            for key, new_name in md_residue_names.items():
                parts = key.split("_")
                if len(parts) >= 3:
                    chain_id = parts[0]
                    try:
                        res_num = int(parts[1])
                        res_type = "_".join(
                            parts[2:]
                        )  # Join in case there are underscores in residue type

                        residue_map[(chain_id, res_num, res_type)] = new_name
                    except ValueError:
                        continue

            # Track which residues were modified
            modified_residues = []

            # Process the input file and create updated output
            with open(input_file, "r") as infile, open(output_file, "w") as outfile:
                for line in infile:
                    if line.startswith(("ATOM", "HETATM")):
                        chain_id = line[21]
                        try:
                            res_num = int(line[22:26])
                            res_type = line[17:20].strip()

                            # Check if this residue has a new name
                            if (chain_id, res_num, res_type) in residue_map:
                                new_name = residue_map[(chain_id, res_num, res_type)]

                                # Create updated line with new residue name
                                updated_line = line[:17] + f"{new_name:<3}" + line[20:]
                                outfile.write(updated_line)

                                # Track that we modified this residue (only once per residue)
                                if (
                                    chain_id,
                                    res_num,
                                    res_type,
                                    new_name,
                                ) not in modified_residues:
                                    modified_residues.append(
                                        (chain_id, res_num, res_type, new_name)
                                    )

                                continue
                        except ValueError:
                            pass

                    # Write original line if not modified
                    outfile.write(line)

            if modified_residues:
                self.processor.console.print(
                    f"[green]Modified {len(modified_residues)} residue names in PDB file[/green]"
                )

                # Log details of changes at debug level
                debug_enabled = self.get_from_workspace("debug", False)
                if debug_enabled:
                    for chain, resnum, old_type, new_type in modified_residues:
                        self.processor.console.print(
                            f"[blue]Changed Chain {chain}, Residue {resnum} from {old_type} to {new_type}[/blue]"
                        )
            else:
                self.processor.console.print(
                    "[yellow]No residues were modified in the PDB file[/yellow]"
                )

            return output_file

        except Exception as e:
            self.processor.console.print(
                f"[red]Error updating residue names in PDB file: {str(e)}[/red]"
            )
            return None

    def _save_results(self):
        """Save analysis results to file"""
        console = self.processor.console

        # Check if we have microstate results
        microstate_results = self.get_from_workspace("microstate_protonation_results", {})

        if microstate_results and len(microstate_results) > 0:
            # We have microstate results - ask user what to save
            console.print(f"\n[bold cyan]Available Results:[/bold cyan]")
            console.print(f"  • Default structure results: {'Yes' if self.get_from_workspace('protonation_results') else 'No'}")
            console.print(f"  • Microstate results: {len(microstate_results)} microstate(s)")


            choices = ["1", "2"]
            options_map = {
                "1": "Save default structure results",
                "2": "Save microstate results"
            }

            choice = prompt_with_context(
                processor=self.processor,
                prompt="What would you like to save?",
                choices=choices,
                default="2",
                module="Protonation State Analyzer",
                description="Select results to save",
                options_map=options_map
            )

            if choice == "2":
                self._save_microstate_results()
                return

        # Save default structure results
        if not self.analysis_results:
            results = self.get_from_workspace("protonation_results")
            if results:
                self.analysis_results = results
            else:
                self.processor.console.print(
                    "[yellow]No analysis results available. Run analysis first.[/yellow]"
                )
                if confirm_with_context(
                    processor=self.processor,
                    prompt="Run analysis now?",
                    module="Protonation State Analyzer",
                    description="Run analysis now"
                ):
                    self._analyze_protonation_states()
                return

        # Ask for output format
        format_choice = prompt_with_context(
            processor=self.processor,
            prompt="Select output format",
            choices=["text", "csv"],
            default="text",
            module="Protonation State Analyzer",
            description="Select output format",
            options_map={"text": "Text format", "csv": "CSV format"}
        )

        # Ask for output filename
        output_file = prompt_with_context(
            processor=self.processor,
            prompt="Enter output filename",
            default=f"protonation_analysis.{format_choice}",
            module="Protonation State Analyzer",
            description="Enter output filename"
        )

        # Save results
        with self.processor.console.status("[bold green]Saving analysis results..."):
            saved_file = self.analyzer.save_results(
                output_file, output_format=format_choice
            )

        if saved_file:
            self.processor.console.print(
                f"[green]Saved analysis results to {saved_file}[/green]"
            )

    def _save_microstate_results(self):
        """Save protonation analysis results for microstates"""
        console = self.processor.console
        microstate_results = self.get_from_workspace("microstate_protonation_results", {})

        if not microstate_results:
            console.print("[yellow]No microstate results available.[/yellow]")
            return

        # Display table of available microstates
        from rich.table import Table
        table = Table(title="Analyzed Microstates")
        table.add_column("#", style="cyan", justify="right")
        table.add_column("Filename", style="green")
        table.add_column("pH", style="yellow")
        table.add_column("Net Charge", style="magenta")

        microstate_list = list(microstate_results.items())
        for i, (filepath, data) in enumerate(microstate_list, 1):
            table.add_row(
                str(i),
                data['filename'],
                f"{data['pH']:.1f}",
                f"{data['net_charge']:.2f}"
            )

        console.print(table)

        # Ask user to select which ones to save

        console.print("\n[bold]Selection Options:[/bold]")
        console.print("  • Single number: e.g., '3'")
        console.print("  • Multiple: e.g., '1,3,5'")
        console.print("  • Range: e.g., '1-4'")
        console.print("  • All: enter 'all'")

        selection = prompt_with_context(
            processor=self.processor,
            prompt="Select microstate(s) to save",
            default="all",
            module="Protonation State Analyzer",
            description="Microstate selection for saving"
        )

        # Parse selection
        selected_indices = self._parse_microstate_selection(selection, len(microstate_list))

        if not selected_indices:
            console.print("[yellow]No valid microstates selected.[/yellow]")
            return

        # Ask for output format (once for all)
        format_choice = prompt_with_context(
            processor=self.processor,
            prompt="Select output format",
            choices=["text", "csv"],
            default="text",
            module="Protonation State Analyzer",
            description="Select output format",
            options_map={"text": "Text format", "csv": "CSV format"}
        )

        # Save results for each selected microstate
        saved_count = 0
        for idx in selected_indices:
            filepath, data = microstate_list[idx]

            # Generate output filename based on microstate filename
            from pathlib import Path
            base_name = Path(data['filename']).stem
            output_file = f"{base_name}_protonation.{format_choice}"

            console.print(f"\n[cyan]Saving results for {data['filename']}...[/cyan]")

            # Temporarily set the analysis results for this microstate
            original_results = self.analysis_results
            self.analysis_results = data['results']

            with console.status(f"[bold green]Saving {output_file}..."):
                saved_file = self.analyzer.save_results(
                    output_file, output_format=format_choice
                )

            # Restore original results
            self.analysis_results = original_results

            if saved_file:
                console.print(f"[green]✓ Saved: {saved_file}[/green]")
                saved_count += 1
            else:
                console.print(f"[red]✗ Failed to save: {output_file}[/red]")

        console.print(f"\n[bold green]✓ Saved {saved_count}/{len(selected_indices)} microstate result(s)[/bold green]")

    def _generate_ph_profile(self):
        """Generate pH titration profile"""
        # Check if analyzer is initialized
        if (
            not hasattr(self.analyzer, "titratable_residues")
            or not self.analyzer.titratable_residues
        ):
            self.processor.console.print(
                "[yellow]Analyzer not initialized. Run analysis first.[/yellow]"
            )
            if confirm_with_context(
                processor=self.processor,
                prompt="Run analysis now?",
                module="Protonation State Analyzer",
                description="Run analysis now"
            ):
                self._analyze_protonation_states()
            else:
                return

        # Ask for pH range
        min_pH = prompt_with_context(
            processor=self.processor,
            prompt="Enter minimum pH",
            default="0",
            module="Protonation State Analyzer",
            description="Enter minimum pH"
        )
        max_pH = prompt_with_context(
            processor=self.processor,
            prompt="Enter maximum pH",
            default="14",
            module="Protonation State Analyzer",
            description="Enter maximum pH"
        )
        step = prompt_with_context(
            processor=self.processor,
            prompt="Enter pH step size",
            default="0.5",
            module="Protonation State Analyzer",
            description="Enter pH step size"
        )

        try:
            min_pH = float(min_pH)
            max_pH = float(max_pH)
            step = float(step)
        except ValueError:
            self.processor.console.print(
                "[red]Invalid pH values, using defaults (0-14, step 0.5)[/red]"
            )
            min_pH = 0.0
            max_pH = 14.0
            step = 0.5

        # Generate profile
        with self.processor.console.status("[bold green]Generating pH profile..."):
            profile = self.analyzer.generate_ph_profile(
                pH_range=(min_pH, max_pH), step=step
            )

        # Display summary
        self._display_ph_profile(profile)

        # Ask if user wants to save the profile
        if confirm_with_context(
            processor=self.processor,
            prompt="Save pH profile to file?",
            default=True,
            module="Protonation State Analyzer",
            description="Save pH profile to file"
        ):
            output_file = prompt_with_context(
                processor=self.processor,
                prompt="Enter output filename",
                default="ph_profile.csv",
                module="Protonation State Analyzer",
                description="Enter pH profile filename"
            )
            saved_file = self.analyzer.save_ph_profile(profile, output_file)
            if saved_file:
                self.processor.console.print(
                    f"[green]Saved pH profile to {saved_file}[/green]"
                )

    def _display_summary(self, results, net_charge, method="Textbook", textbook_pkas=None):
        """Display a summary of protonation analysis results

        Args:
            results: Analysis results dictionary
            net_charge: Net charge of the structure
            method: "PROPKA" or "Textbook" - indicates which pKa method was used
            textbook_pkas: Optional dict of textbook pKa values (for delta calculation when using PROPKA)
        """
        # Title indicates which method was used
        table = Table(title=f"\nProtonation States at pH {self.analyzer.pH} ({method} pKa)")
        table.add_column("Index", style="grey50", justify="right", width=6)
        table.add_column("Chain", style="cyan", justify="right")
        table.add_column("Residue", style="green", justify="right")
        table.add_column("pKa", style="yellow", justify="right")

        # Add delta column if using PROPKA
        if method == "PROPKA" and textbook_pkas:
            table.add_column("ΔpKa", style="bright_yellow", justify="right")

        # Desolvation is the dominant pKa driver for a buried group, and it's the
        # one contribution with no named partner — so it's invisible in the
        # partner-oriented inspect view unless you drill in. Surface its
        # magnitude here (regular + RE, as PROPKA reports it) so a burial-driven
        # shift is legible straight from the summary. PROPKA-only.
        show_desolv = (method == "PROPKA")
        if show_desolv:
            table.add_column("Desolv", style="#ff00ff", justify="right")

        table.add_column("Protonated", style="magenta", justify="right")
        table.add_column("Probability", style="blue", justify="right")
        table.add_column("Charge", style="red", justify="right")

        # Organize results by chain
        by_chain = {}
        for key, res in results.items():
            chain = res["chain"]
            if chain not in by_chain:
                by_chain[chain] = []
            by_chain[chain].append((key, res))

        # Build indexed list for visualization (store for later use)
        indexed_residues = []
        global_index = 1

        # Add rows for each chain
        for chain, residues in sorted(by_chain.items()):
            # Sort residues by number
            sorted_residues = sorted(residues, key=lambda x: x[1]["number"])

            for key, res in sorted_residues:
                # Store indexed residue info for visualization
                delta_pka = None
                if method == "PROPKA" and textbook_pkas:
                    textbook_pka = textbook_pkas.get(key, res['pKa'])
                    delta_pka = res['pKa'] - textbook_pka

                indexed_residues.append({
                    'index': global_index,
                    'key': key,
                    'chain': res["chain"],
                    'resnum': res["number"],
                    'restype': res["type"],
                    'pKa': res["pKa"],
                    'delta_pka': delta_pka,
                    'protonated': res["protonated"],
                    'probability': res["probability"],
                    'charge': res["charge"],
                    'determinants': res.get("determinants", {})
                })

                row = [
                    str(global_index),
                    res["chain"],
                    f"{res['type']} {res['number']}",
                    f"{res['pKa']:.2f}",
                ]

                # Add delta if using PROPKA
                if method == "PROPKA" and textbook_pkas:
                    if abs(delta_pka) >= 0.01:  # Only show meaningful deltas
                        delta_str = f"{delta_pka:+.2f}"
                    else:
                        delta_str = "—"
                    row.append(delta_str)

                # Desolvation magnitude (regular + RE). Shown as a bare positive
                # number the way PROPKA prints it — it's a penalty on the charged
                # state, so it pushes an acid's pKa up and a base's pKa down; the
                # direction is already encoded in ΔpKa, so we don't sign it here.
                if show_desolv:
                    det = res.get("determinants") or {}
                    if det:
                        dv = (det.get("desolvation_regular", 0.0)
                              + det.get("desolvation_re", 0.0))
                        row.append(f"{dv:.2f}" if abs(dv) >= 0.01 else "—")
                    else:
                        row.append("—")

                row.extend([
                    str(res["protonated"]),
                    f"{res['probability']:.4f}",
                    f"{res['charge']:.2f}",
                ])

                table.add_row(*row)
                global_index += 1

        self.processor.console.print(table)
        self.processor.console.print(f"[bold]Net Charge:[/bold] {net_charge:.2f}")

        # Display information about excluded residues
        if self.analyzer.excluded_residues:
            # Group excluded residues by chain
            excluded_by_chain = {}
            for chain, resnum, restype in self.analyzer.excluded_residues:
                if chain not in excluded_by_chain:
                    excluded_by_chain[chain] = []
                excluded_by_chain[chain].append((resnum, restype))

            self.processor.console.print(
                "\n[bold yellow]Excluded residues (user-selected):[/bold yellow]"
            )
            for chain, residues in sorted(excluded_by_chain.items()):
                residue_str = ", ".join(
                    [f"{restype} {resnum}" for resnum, restype in sorted(residues)]
                )
                self.processor.console.print(f"  Chain {chain}: {residue_str}")

        # Offer visualization of titratable residue environments
        # Get structure file using priority-based selection (same as structure_selector)
        from proprep.utils.structure_selector import StructureSelector

        workspace = self.processor._get_workspace()
        structure_selector = StructureSelector(workspace, self.processor.console)
        structure_file = structure_selector.get_structure(silent=True)

        if indexed_residues and structure_file and not getattr(self, '_batch_mode', False):
            # One index prompt drives both the structure view and the PROPKA
            # determinant breakdown for the chosen residues (Enter = skip).
            self._inspect_titratable_residues(indexed_residues)

    def _inspect_titratable_residues(self, indexed_residues):
        """Pick residues once, then show each one's environment + PROPKA pKa
        determinants in the structure viewer.

        A single index prompt replaces the old visualize→inspect two-step:
        entering indices both opens the viewer and prints the determinant
        breakdown for those residues; pressing Enter skips. PROPKA explains each
        predicted pKa as a desolvation penalty plus named side-chain / backbone /
        Coulombic determinants. This surfaces that breakdown — including which
        contributions come from cofactors/ligands — so the user can see *why* a
        residue's pKa is shifted (e.g. a nearby FMN phosphate pulling an Asp
        down). Residues PROPKA didn't explain are still shown in the viewer.
        """
        from rich.table import Table

        self.processor.console.print(
            "\n[grey50]Tip: visualize a residue to see its 10Å environment and the "
            "PROPKA determinants driving its pKa shift.[/grey50]"
        )
        selection = prompt_with_context(
            processor=self.processor,
            prompt="Enter indices to visualize (e.g., 1,3,5-7), or Enter to skip",
            default="",
            module="Protonation State Analyzer",
            description="Select residues to visualize and inspect determinants",
        ).strip()
        if not selection:
            return
        try:
            selected_indices = self._parse_index_selection(selection, len(indexed_residues))
        except ValueError as e:
            self.processor.console.print(f"[red]Error: {e}[/red]")
            return

        selected = [indexed_residues[i - 1] for i in selected_indices]
        viewer_payload = []  # (center_res, [determinant dicts])

        for res_info in selected:
            det = res_info.get('determinants') or {}
            head = (f"{res_info['restype']} {res_info['chain']}:{res_info['resnum']}"
                    f"  (predicted pKa {res_info['pKa']:.2f}")
            if res_info.get('delta_pka') is not None:
                head += f", Δ {res_info['delta_pka']:+.2f} vs textbook"
            head += ")"
            self.processor.console.print(f"\n[bold #ff00ff]{head}[/bold #ff00ff]")

            dets = det.get('determinants', [])
            if not dets:
                self.processor.console.print(
                    "  [grey50]No PROPKA determinants (pKa near its model value, "
                    "or default value used).[/grey50]"
                )
                # Still show it in the viewer (environment only, no drivers).
                viewer_payload.append((res_info, []))
                continue

            desolv = det.get('desolvation_regular', 0.0) + det.get('desolvation_re', 0.0)
            self.processor.console.print(
                f"  [#ff00ff]Buried {det.get('buried', '?')}; desolvation "
                f"{desolv:+.2f} (raises pKa)[/#ff00ff]"
            )

            table = Table(show_header=True, header_style="bold")
            table.add_column("Contribution", justify="right")
            table.add_column("Type")
            table.add_column("Partner")
            table.add_column("Effect on pKa")
            # Largest-magnitude contributors first.
            for d in sorted(dets, key=lambda x: -abs(x['value'])):
                partner = f"{d['partner_resname']} {d['partner_id']} {d['partner_chain']}"
                if d['is_ligand']:
                    partner += "  [magenta](ligand/cofactor)[/magenta]"
                # Negative determinant lowers the pKa (stabilizes the charged
                # state); positive raises it. Colors match the viewer overlay.
                effect = ("[bold #15803d]lowers (↓)[/bold #15803d]" if d['value'] < 0
                          else "[bold #c2410c]raises (↑)[/bold #c2410c]")
                table.add_row(f"{d['value']:+.2f}", d['kind'], partner, effect)
            self.processor.console.print(table)
            viewer_payload.append((res_info, dets))

        if viewer_payload:
            self._launch_determinant_viewer(viewer_payload)

    def _launch_determinant_viewer(self, payload):
        """Show each inspected residue and its pKa determinants in the viewer.

        Per residue: the titratable residue itself as magenta ball+stick, its
        pKa-lowering partners green and pKa-raising partners orange (also
        ball+stick, so the drivers dominate), over a thin line-rendered 10Å
        neighborhood that gives rough solvent-accessibility context without the
        clutter of a full element-coloured environment. Ligand/cofactor partners
        are selected by resname so the whole cofactor is visible.

        Green/orange are used instead of blue/red so the overlay isn't confused
        with the default blue nitrogen / red oxygen element colors elsewhere in
        the structure.

        force=True on the show only — the user explicitly chose to inspect, so the
        viewer must open in CLI mode too; subsequent highlight calls then update
        whatever viewer is running, web-shell or CLI.
        """
        from proprep.utils.structure_selector import StructureSelector

        workspace = self.processor._get_workspace()
        structure_file = StructureSelector(
            workspace, self.processor.console).get_structure(silent=True)
        if not structure_file:
            self.processor.console.print("[yellow]No PDB file found in workspace[/yellow]")
            return
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer
        except Exception:
            self.processor.console.print("[yellow]Structure viewer not available[/yellow]")
            return

        _viewer.show_structure(structure_file, force=True)

        focused_yet = False
        for res_info, dets in payload:
            index = res_info['index']
            key = f"determinant_{index}"
            center = f":{res_info['chain']} and {res_info['resnum']}"
            _viewer.unhighlight(f"{key}_center")
            _viewer.unhighlight(f"{key}_env")
            _viewer.unhighlight(f"{key}_lower")
            _viewer.unhighlight(f"{key}_raise")

            # Center residue magenta; the first focused=True call also hides the
            # default reps and frames the view on the residue. Use hex, not the
            # named "magenta": NGL reads a bare color word as a (nonexistent)
            # color-scheme id and silently drops the representation, so the
            # target residue never appears — only hex is treated as a uniform
            # color (matches color_map's "#ff00ff" for magenta).
            _viewer.highlight(
                center, style="ball+stick", color="#ff00ff",
                label=f"{key}_center", focused=not focused_yet,
            )
            focused_yet = True

            # Thin line-rendered 10Å neighborhood for rough solvent-accessibility
            # context; the ball+stick drivers below are drawn on top and dominate.
            neighborhood = self._build_neighborhood_selection(
                [res_info], radius=10.0, exclude_center=True,
            )
            if neighborhood:
                _viewer.highlight(
                    neighborhood, style="line", color="element",
                    label=f"{key}_env",
                )

            lower_sel, raise_sel = [], []
            for d in dets:
                if d['is_ligand']:
                    sel = f"[{d['partner_resname']}] and :{d['partner_chain']}"
                else:
                    sel = f":{d['partner_chain']} and {d['partner_id']}"
                (lower_sel if d['value'] < 0 else raise_sel).append(sel)

            if lower_sel:
                _viewer.highlight(
                    " or ".join(f"({s})" for s in dict.fromkeys(lower_sel)),
                    style="ball+stick", color="#15803d", label=f"{key}_lower",
                )
            if raise_sel:
                _viewer.highlight(
                    " or ".join(f"({s})" for s in dict.fromkeys(raise_sel)),
                    style="ball+stick", color="#c2410c", label=f"{key}_raise",
                )

        self.processor.console.print(
            "\n[cyan]Highlighted:[/cyan] [#ff00ff]magenta = residue[/#ff00ff], "
            "[#15803d]green = lowers pKa[/#15803d], "
            "[#c2410c]orange = raises pKa[/#c2410c]; thin lines = 10Å neighborhood. "
            "Toggle reps in the viewer's manager.\n"
        )

    def _build_neighborhood_selection(self, residues, radius=10.0, exclude_center=False):
        """
        Build NGL selection string for residues and their neighborhoods.

        Args:
            residues: List of residue dictionaries with 'chain' and 'resnum'
            radius: Neighborhood radius in Angstroms (default: 10.0)
            exclude_center: If True, exclude the center residue from the selection

        Returns:
            NGL selection string or empty string on error
        """
        from Bio.PDB import PDBParser
        from Bio.PDB.NeighborSearch import NeighborSearch

        structure = self.analyzer.structure

        # If structure is not loaded but we have an input file, load it
        if not structure and self.analyzer.input_file:
            try:
                parser = PDBParser(QUIET=True)
                structure = parser.get_structure("protein", self.analyzer.input_file)
            except Exception as e:
                self.processor.console.print(f"[yellow]Could not load structure from {self.analyzer.input_file}: {e}[/yellow]")
                return ""

        if not structure:
            return ""

        try:
            model = list(structure.get_models())[0]
        except (IndexError, StopIteration):
            self.processor.console.print("[yellow]No model found in structure[/yellow]")
            return ""

        # Get all atoms for NeighborSearch
        all_atoms = list(model.get_atoms())
        if not all_atoms:
            self.processor.console.print("[yellow]No atoms found in structure[/yellow]")
            return ""

        ns = NeighborSearch(all_atoms)

        # Collect all residues to display (selected + neighbors)
        residues_to_show = set()
        center_residues = set()  # Track center residues for exclusion if needed

        for res_info in residues:
            chain = res_info['chain']
            resnum = res_info['resnum']
            restype = res_info['restype']

            # Track this as a center residue
            center_residues.add((chain, resnum))

            # Add the titratable residue itself (unless excluding)
            if not exclude_center:
                residues_to_show.add((chain, resnum))

            # Find center of this residue and its neighbors
            try:
                # Handle residue ID - may need to account for insertion codes
                target_residue = None
                for residue in model[chain].get_residues():
                    if residue.id[1] == resnum:
                        target_residue = residue
                        break

                if not target_residue:
                    self.processor.console.print(
                        f"[yellow]Warning: Could not find {restype} {chain}:{resnum} in structure[/yellow]"
                    )
                    continue

                # Find all atoms in this residue
                center_atoms = list(target_residue.get_atoms())

                # Find neighbors within radius for each atom in the residue
                for center_atom in center_atoms:
                    neighbors = ns.search(center_atom.coord, radius)
                    for neighbor_atom in neighbors:
                        neighbor_res = neighbor_atom.get_parent()
                        neighbor_chain = neighbor_res.get_parent().id
                        neighbor_resnum = neighbor_res.id[1]
                        neighbor_key = (neighbor_chain, neighbor_resnum)

                        # If excluding center, skip center residues
                        if exclude_center and neighbor_key in center_residues:
                            continue

                        residues_to_show.add(neighbor_key)

            except Exception as e:
                self.processor.console.print(
                    f"[yellow]Warning: Could not find neighbors for {restype} {chain}:{resnum}: {e}[/yellow]"
                )
                continue

        # Build NGL selection string
        if not residues_to_show:
            return ""

        selections = [f":{chain} and {resnum}" for chain, resnum in sorted(residues_to_show)]
        ngl_selection = " or ".join(selections)

        self.processor.console.print(
            f"[grey50]Showing {len(residues_to_show)} residues total (selected + neighbors)[/grey50]"
        )

        return ngl_selection

    def _calculate_isoelectric_point(self, profile):
        """
        Calculate the isoelectric point (pH where net charge = 0).
        
        Args:
            profile: Dictionary with 'pH' and 'net_charge' lists
        
        Returns:
            float or None: The isoelectric point, or None if not found
        """
        pH_values = profile['pH']
        charges = profile['net_charge']
        
        if not pH_values or not charges:
            return None
        
        # Find the point where charge crosses zero
        for i in range(len(charges) - 1):
            charge1 = charges[i]
            charge2 = charges[i + 1]
            
            # Check if zero is between these two points
            if (charge1 <= 0 <= charge2) or (charge2 <= 0 <= charge1):
                # Linear interpolation to find exact pH
                pH1 = pH_values[i]
                pH2 = pH_values[i + 1]
                
                if charge2 == charge1:  # Avoid division by zero
                    return pH1
                
                # Linear interpolation: pH = pH1 + (0 - charge1) * (pH2 - pH1) / (charge2 - charge1)
                isoelectric_pH = pH1 + (0 - charge1) * (pH2 - pH1) / (charge2 - charge1)
                return isoelectric_pH
        
        return None

    def _interpolate_value(self, profile, target_value, lookup_type):
        """
        Interpolate to find pH from charge or charge from pH.
        
        Args:
            profile: Dictionary with 'pH' and 'net_charge' lists
            target_value: The value to look up
            lookup_type: 'charge' to find pH from charge, 'pH' to find charge from pH
        
        Returns:
            float or None: The interpolated value, or None if not found
        """
        pH_values = profile['pH']
        charges = profile['net_charge']
        
        if not pH_values or not charges:
            return None
        
        if lookup_type == 'charge':
            # Find pH from charge
            x_values = charges
            y_values = pH_values
            target = target_value
        else:
            # Find charge from pH
            x_values = pH_values
            y_values = charges
            target = target_value
        
        # Check if target is within range
        min_x = min(x_values)
        max_x = max(x_values)
        
        if target < min_x or target > max_x:
            return None
        
        # Find the closest points for interpolation
        for i in range(len(x_values) - 1):
            x1 = x_values[i]
            x2 = x_values[i + 1]
            
            if (x1 <= target <= x2) or (x2 <= target <= x1):
                y1 = y_values[i]
                y2 = y_values[i + 1]
                
                if x2 == x1:  # Avoid division by zero
                    return y1
                
                # Linear interpolation
                result = y1 + (target - x1) * (y2 - y1) / (x2 - x1)
                return result
        
        return None

    def _interactive_lookup(self, profile):
        """
        Interactive lookup for pH/charge values.
        User can enter charge values with 'e' suffix or pH values without suffix.
        """
        self.processor.console.print("\n[bold]Interactive pH/Charge Lookup[/bold]")
        self.processor.console.print("[italic]Enter values to look up. Type 'done' to finish.[/italic]")
        self.processor.console.print("[italic]Examples: '7.0' (finds charge at pH 7.0), '-1.5e' (finds pH at charge -1.5)[/italic]")

        while True:
            try:
                user_input = prompt_with_context(
                    processor=self.processor,
                    prompt="\nEnter value (or 'done')",
                    module="Protonation State Analyzer",
                    description="Enter lookup value"
                ).strip()
                
                if user_input.lower() == 'done':
                    break
                
                # Parse the input
                if user_input.endswith('e') or user_input.endswith('E'):
                    # Looking up pH from charge
                    charge_str = user_input[:-1]
                    try:
                        target_charge = float(charge_str)
                        result_pH = self._interpolate_value(profile, target_charge, 'charge')
                        
                        if result_pH is not None:
                            self.processor.console.print(
                                f"[green]At charge {target_charge:.2f}: pH = {result_pH:.2f}[/green]"
                            )
                        else:
                            min_charge = min(profile['net_charge'])
                            max_charge = max(profile['net_charge'])
                            self.processor.console.print(
                                f"[yellow]Charge {target_charge:.2f} is outside the range "
                                f"[{min_charge:.2f}, {max_charge:.2f}][/yellow]"
                            )
                    except ValueError:
                        self.processor.console.print(
                            f"[red]Invalid charge value: '{charge_str}'[/red]"
                        )
                else:
                    # Looking up charge from pH
                    try:
                        target_pH = float(user_input)
                        result_charge = self._interpolate_value(profile, target_pH, 'pH')
                        
                        if result_charge is not None:
                            self.processor.console.print(
                                f"[green]At pH {target_pH:.2f}: charge = {result_charge:.2f}[/green]"
                            )
                        else:
                            min_pH = min(profile['pH'])
                            max_pH = max(profile['pH'])
                            self.processor.console.print(
                                f"[yellow]pH {target_pH:.2f} is outside the range "
                                f"[{min_pH:.2f}, {max_pH:.2f}][/yellow]"
                            )
                    except ValueError:
                        self.processor.console.print(
                            f"[red]Invalid pH value: '{user_input}'[/red]"
                        )
            
            except KeyboardInterrupt:
                self.processor.console.print("\n[yellow]Lookup cancelled.[/yellow]")
                break
            except Exception as e:
                self.processor.console.print(f"[red]Error: {str(e)}[/red]")
        
    def _interactive_capacitance_lookup(self, profile):
        """
        Interactive lookup for pH/capacitance values.
        User can enter capacitance values with 'e' suffix or pH values without suffix.
        """
        self.processor.console.print("\n[bold]Interactive pH/Capacitance Lookup[/bold]")
        self.processor.console.print("Enter values to look up. Type 'done' to finish.")
        self.processor.console.print("Examples: '7.0' (finds capacitance at pH 7.0), '0.5e' (finds pH at capacitance 0.5)")

        while True:
            try:
                user_input = prompt_with_context(
                    processor=self.processor,
                    prompt="\nEnter value (or 'done')",
                    module="Protonation State Analyzer",
                    description="Enter capacitance lookup value"
                ).strip()
                
                if user_input.lower() == 'done':
                    break
                
                try:
                    pH_values = profile['pH_values']
                    capacitances = profile['total_capacitance']
                    
                    if user_input.endswith('e'):
                        # Find pH at given capacitance
                        target_cap = float(user_input[:-1])
                        
                        if target_cap < min(capacitances) or target_cap > max(capacitances):
                            self.processor.console.print(
                                f"[yellow]Capacitance {target_cap:.4f} is outside the range "
                                f"[{min(capacitances):.4f}, {max(capacitances):.4f}][/yellow]"
                            )
                            continue
                        
                        # Find pH using linear interpolation
                        result_pH = self._interpolate_value(capacitances, pH_values, target_cap)
                        if result_pH is not None:
                            self.processor.console.print(
                                f"[green]At capacitance {target_cap:.4f}: pH = {result_pH:.2f}[/green]"
                            )
                        else:
                            self.processor.console.print(
                                f"[yellow]Could not find pH for capacitance {target_cap:.4f}[/yellow]"
                            )
                    else:
                        # Find capacitance at given pH
                        target_pH = float(user_input)
                        
                        if target_pH < min(pH_values) or target_pH > max(pH_values):
                            self.processor.console.print(
                                f"[yellow]pH {target_pH:.2f} is outside the range "
                                f"[{min(pH_values):.2f}, {max(pH_values):.2f}][/yellow]"
                            )
                            continue
                        
                        # Find capacitance using linear interpolation
                        result_cap = self._interpolate_value(pH_values, capacitances, target_pH)
                        if result_cap is not None:
                            self.processor.console.print(
                                f"[green]At pH {target_pH:.2f}: capacitance = {result_cap:.4f}[/green]"
                            )
                        else:
                            self.processor.console.print(
                                f"[yellow]Could not find capacitance for pH {target_pH:.2f}[/yellow]"
                            )
                            
                except ValueError:
                    self.processor.console.print(
                        f"[red]Invalid input: '{user_input}'. Use numbers for pH or numbers with 'e' for capacitance.[/red]"
                    )
            
            except KeyboardInterrupt:
                self.processor.console.print("\n[yellow]Lookup cancelled.[/yellow]")
                break
            except Exception as e:
                self.processor.console.print(f"[red]Error: {str(e)}[/red]")
                
    def _plot_ascii_titration_curve(self, profile, width=60, height=20, isoelectric_point=None):
        """
        Create an ASCII plot of the pH titration curve.

        Args:
            profile: Dictionary with 'pH' and 'net_charge' lists
            width: Width of the plot in characters
            height: Height of the plot in characters
            isoelectric_point: The isoelectric point (pI) value, or None if not found

        Returns:
            String containing the ASCII plot
        """
        import math
        
        pH_values = profile['pH']
        charges = profile['net_charge']
        
        if not pH_values or not charges:
            return "No data to plot"
        
        # Get data ranges
        min_pH = min(pH_values)
        max_pH = max(pH_values)
        min_charge = min(charges)
        max_charge = max(charges)
        
        # Add some padding to the charge range
        charge_range = max_charge - min_charge
        if charge_range == 0:
            charge_range = 1
        padding = charge_range * 0.1
        min_charge -= padding
        max_charge += padding
        
        # Create the plot grid
        plot_lines = []
        
        # Create each row from top to bottom
        for row in range(height):
            line = []
            
            # Calculate the charge value for this row (top = max_charge, bottom = min_charge)
            row_charge = max_charge - (row / (height - 1)) * (max_charge - min_charge)
            
            for col in range(width):
                # Calculate the pH value for this column
                col_pH = min_pH + (col / (width - 1)) * (max_pH - min_pH)
                
                # Find the closest data point
                closest_idx = min(range(len(pH_values)), 
                                key=lambda i: abs(pH_values[i] - col_pH))
                closest_charge = charges[closest_idx]
                
                # Determine what character to place
                char = ' '
                
                # Check if this position should have a data point
                charge_tolerance = (max_charge - min_charge) / height * 0.6
                if abs(closest_charge - row_charge) <= charge_tolerance:
                    char = '*'
                
                # Add axis characters
                elif col == 0:  # Y-axis
                    char = '│'
                elif row == height - 1:  # X-axis
                    char = '─'
                elif col == 0 and row == height - 1:  # Origin
                    char = '└'
                
                line.append(char)
            
            plot_lines.append(''.join(line))
        
        # Add labels and build the complete plot
        plot_result = []
        
        # Title
        plot_result.append(f"{'pH Titration Curve':^{width + 10}}")
        plot_result.append("")
        
        # Y-axis label and plot
        for i, line in enumerate(plot_lines):
            if i == 0:
                charge_label = f"{max_charge:6.1f}"
            elif i == height // 2:
                charge_label = f"{(max_charge + min_charge) / 2:6.1f}"
            elif i == height - 1:
                charge_label = f"{min_charge:6.1f}"
            else:
                charge_label = "      "
            
            plot_result.append(f"{charge_label} {line}")
        
        # X-axis labels
        x_labels = "       "
        for col in range(0, width, max(1, width // 8)):
            col_pH = min_pH + (col / (width - 1)) * (max_pH - min_pH)
            x_labels += f"{col_pH:5.1f}"
            if col < width - 5:
                x_labels += " " * (max(1, width // 8) - 5)
        
        plot_result.append(x_labels)
        plot_result.append(f"{'pH':^{width + 7}}")

        # Add isoelectric point information
        plot_result.append("")
        if isoelectric_point is not None:
            plot_result.append(f"Isoelectric point (pI): {isoelectric_point:.2f}")
        else:
            plot_result.append("Isoelectric point (pI): not found in pH range")

        return "\n".join(plot_result)

    def _plot_ascii_capacitance_curve(self, profile, width=60, height=20):
        """
        Create an ASCII plot of the capacitance vs pH curve.
        
        Args:
            profile: Dictionary with 'pH_values' and 'total_capacitance' lists
            width: Width of the plot in characters
            height: Height of the plot in characters
        
        Returns:
            String containing the ASCII plot
        """
        import math
        
        pH_values = profile['pH_values']
        capacitances = profile['total_capacitance']
        
        if not pH_values or not capacitances:
            return "No data to plot"
        
        # Get data ranges
        min_pH = min(pH_values)
        max_pH = max(pH_values)
        min_cap = min(capacitances)
        max_cap = max(capacitances)
        
        # Add some padding to the capacitance range
        cap_range = max_cap - min_cap
        if cap_range == 0:
            cap_range = 0.1
        padding = cap_range * 0.1
        min_cap = max(0, min_cap - padding)  # Capacitance shouldn't be negative
        max_cap += padding
        
        # Create the plot grid
        plot_lines = []
        
        # Create each row from top to bottom
        for row in range(height):
            line = []
            
            # Calculate the capacitance value for this row
            row_cap = max_cap - (row / (height - 1)) * (max_cap - min_cap)
            
            for col in range(width):
                # Calculate the pH value for this column
                col_pH = min_pH + (col / (width - 1)) * (max_pH - min_pH)
                
                # Find the closest data point
                closest_idx = min(range(len(pH_values)), 
                                key=lambda i: abs(pH_values[i] - col_pH))
                closest_cap = capacitances[closest_idx]
                
                # Determine what character to place
                char = ' '
                
                # Check if this position should have a data point
                cap_tolerance = (max_cap - min_cap) / height * 1.0
                if abs(closest_cap - row_cap) <= cap_tolerance:
                    char = '*'
                
                # Add axis characters
                elif col == 0:  # Y-axis
                    char = '│'
                elif row == height - 1:  # X-axis
                    char = '─'
                elif col == 0 and row == height - 1:  # Origin
                    char = '└'
                
                line.append(char)
            
            plot_lines.append(''.join(line))
        
        # Add labels and build the complete plot
        plot_result = []
        
        # Title
        plot_result.append(f"{'Ideal Charge Capacitance vs pH':^{width + 10}}")
        plot_result.append("")
        
        # Y-axis label and plot
        for i, line in enumerate(plot_lines):
            if i == 0:
                cap_label = f"{max_cap:6.3f}"
            elif i == height // 2:
                cap_label = f"{(max_cap + min_cap) / 2:6.3f}"
            elif i == height - 1:
                cap_label = f"{min_cap:6.3f}"
            else:
                cap_label = "      "
            
            plot_result.append(f"{cap_label} {line}")
        
        # X-axis labels
        x_labels = "       "
        for col in range(0, width, max(1, width // 8)):
            col_pH = min_pH + (col / (width - 1)) * (max_pH - min_pH)
            x_labels += f"{col_pH:5.1f}"
            if col < width - 5:
                x_labels += " " * (max(1, width // 8) - 5)
        
        plot_result.append(x_labels)
        plot_result.append(f"{'pH':^{width + 7}}")
        
        # Add peak information
        plot_result.append("")
        plot_result.append(f"Peak capacitance: {profile['max_capacitance']:.4f} at pH {profile['pH_max_capacitance']:.2f}")
        
        return "\n".join(plot_result)
        
    def _display_ph_profile(self, profile):
        """Display pH profile summary with ASCII plot, isoelectric point, and interactive lookup"""
        self.processor.console.print("[bold]pH Titration Profile:[/bold]")
        self.processor.console.print(
            f"pH range: {profile['pH'][0]} to {profile['pH'][-1]}"
        )
        
        # Calculate and display isoelectric point
        isoelectric_point = self._calculate_isoelectric_point(profile)
        if isoelectric_point is not None:
            self.processor.console.print(
                f"[bold cyan]Isoelectric Point (pI): {isoelectric_point:.2f}[/bold cyan]"
            )
        else:
            self.processor.console.print(
                "[yellow]Isoelectric point not found in the pH range[/yellow]"
            )

        # Display a few key points
        idx_2 = next((i for i, pH in enumerate(profile["pH"]) if pH >= 2), 0)
        idx_7 = next((i for i, pH in enumerate(profile["pH"]) if pH >= 7), 0)
        idx_12 = next((i for i, pH in enumerate(profile["pH"]) if pH >= 12), 0)

        table = Table(title="Selected pH Points")
        table.add_column("pH", style="cyan", justify="right")
        table.add_column("Net Charge", style="green", justify="right")

        # Add key pH points
        table.add_row(f"{profile['pH'][0]:.1f}", f"{profile['net_charge'][0]:.2f}")
        table.add_row(
            f"{profile['pH'][idx_2]:.1f}", f"{profile['net_charge'][idx_2]:.2f}"
        )
        table.add_row(
            f"{profile['pH'][idx_7]:.1f}", f"{profile['net_charge'][idx_7]:.2f}"
        )
        table.add_row(
            f"{profile['pH'][idx_12]:.1f}", f"{profile['net_charge'][idx_12]:.2f}"
        )
        table.add_row(f"{profile['pH'][-1]:.1f}", f"{profile['net_charge'][-1]:.2f}")

        self.processor.console.print(table)
        
        # Generate and display the ASCII plot
        self.processor.console.print("\n[bold]Titration Curve:[/bold]")
        ascii_plot = self._plot_ascii_titration_curve(profile, isoelectric_point=isoelectric_point)
        
        # Print the plot with a monospace font style
        from rich.panel import Panel
        self.processor.console.print(
            Panel(
                f"[font=monospace]{ascii_plot}[/font]",
                title="Net Charge vs pH",
                border_style="bright_blue",
                expand=False
            )
        )

        # Add interactive lookup
        if confirm_with_context(
            processor=self.processor,
            prompt="\nWould you like to do interactive pH/charge lookups?",
            default=False,
            module="Protonation State Analyzer",
            description="Interactive pH/charge lookup"
        ):
            self._interactive_lookup(profile)

        # Add note about visualization
        self.processor.console.print(
            "\n[italic]Save the profile to CSV for high-resolution visualization with external tools.[/italic]"
        )
        
    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Process

    def can_process(self, workspace: Dict[str, Any]) -> bool:
        """Check if the module can process the current workspace"""
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(workspace, self.processor.console)
        status = selector.get_structure_status()
        return status.get("has_any", False)

    def process(self, workspace: Dict[str, Any]) -> Dict[str, Any]:
        """Process the workspace"""
        workspace_obj = workspace

        if self.can_process(workspace_obj):
            # Check if auto-analysis is requested
            if self.get_from_workspace_obj(workspace_obj, "analyze_protonation", False):
                # Get required parameters
                structure = self.get_from_workspace_obj(
                    workspace_obj, "filtered_structure"
                )
                excluded_residues = self.get_from_workspace_obj(
                    workspace_obj, "excluded_residues", set()
                )
                selected_chains = list(
                    self.get_from_workspace_obj(
                        workspace_obj, "filter_selections", {}
                    ).keys()
                )
                pH = self.get_from_workspace_obj(workspace_obj, "default_pH", 7.0)

                # Setup analyzer
                self.analyzer.setup(
                    structure=structure,
                    input_file=self.get_from_workspace_obj(workspace_obj, "pdb_file"),
                    excluded_residues=excluded_residues,
                    selected_chains=selected_chains,
                    pH=pH,
                )

                # Run analysis
                results, net_charge = self.analyzer.analyze_protonation_states()
                self.analysis_results = results

                # Store results in workspace
                workspace_obj = self.update_workspace_obj(
                    workspace_obj, "protonation_results", results
                )
                workspace_obj = self.update_workspace_obj(
                    workspace_obj, "net_charge", net_charge
                )

                # Optionally run PROPKA
                if self.get_from_workspace_obj(workspace_obj, "run_propka", False):
                    success = self.analyzer.run_propka()
                    if success:
                        # Re-analyze with PROPKA results
                        results, net_charge = self.analyzer.analyze_protonation_states()
                        workspace_obj = self.update_workspace_obj(
                            workspace_obj, "protonation_results", results
                        )
                        workspace_obj = self.update_workspace_obj(
                            workspace_obj, "net_charge", net_charge
                        )

                # Set MD residue names if requested
                if self.get_from_workspace_obj(workspace_obj, "set_md_names", False):
                    constant_pH = self.get_from_workspace_obj(
                        workspace_obj, "constant_pH", False
                    )
                    self.md_residue_names = self.analyzer.recommend_md_residue_names(
                        constant_ph_simulation=constant_pH
                    )
                    workspace_obj = self.update_workspace_obj(
                        workspace_obj, "md_residue_names", self.md_residue_names
                    )

                # Update structure if requested
                if (
                    self.get_from_workspace_obj(
                        workspace_obj, "update_structure", False
                    )
                    and self.md_residue_names
                ):
                    import copy

                    structure_copy = copy.deepcopy(structure)
                    updated_structure = self.analyzer.apply_residue_names_to_structure(
                        structure_copy
                    )
                    # Using standardized key name: protonation_pdb_file
                    workspace_obj = self.update_workspace_obj(
                        workspace_obj, "protonation_pdb_file", updated_structure
                    )

        return workspace_obj


# =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
# new features

    def calculate_capacitance(self):
        """Calculate ideal charge capacitance at specified pH"""
        console = self.processor.console

        # Check if we have microstate results
        microstate_results = self.get_from_workspace("microstate_protonation_results", {})

        if microstate_results and len(microstate_results) > 0:
            # We have microstate results - ask user what to calculate
            console.print(f"\n[bold cyan]Available Results:[/bold cyan]")
            console.print(f"  • Default structure results: {'Yes' if self.get_from_workspace('protonation_results') else 'No'}")
            console.print(f"  • Microstate results: {len(microstate_results)} microstate(s)")


            choices = ["1", "2"]
            options_map = {
                "1": "Calculate for default structure",
                "2": "Calculate for microstate(s)"
            }

            choice = prompt_with_context(
                processor=self.processor,
                prompt="Calculate capacitance for which structure?",
                choices=choices,
                default="2",
                module="Protonation State Analyzer",
                description="Select structure for capacitance",
                options_map=options_map
            )

            if choice == "2":
                self._calculate_capacitance_for_microstates()
                return True

        # Calculate for default structure
        if not self.analyzer.titratable_residues:
            self.processor.console.print(
                "[yellow]No titratable residues analyzed. Run analysis first.[/yellow]"
            )
            return

        # Explain what capacitance calculation does before asking for pH
        self.processor.console.print("\n[bold cyan]Ideal Charge Capacitance[/bold cyan]")
        self.processor.console.print(
            "Capacitance measures how sensitively your protein's charge responds to pH changes."
        )
        self.processor.console.print("\n[bold]What you'll learn:[/bold]")
        self.processor.console.print("• [cyan]Total capacitance[/cyan] - Overall pH buffering capacity (scales with protein size)")
        self.processor.console.print("• [cyan]Specific capacitance[/cyan] - Average per residue (enables size-independent comparisons)")
        self.processor.console.print("  [grey50]→ Maximum per-residue value is 0.25 (occurs when pH = pKa for that residue)[/grey50]")
        self.processor.console.print("• [cyan]Per-residue breakdown[/cyan] - Which residue types contribute most to pH buffering")

        self.processor.console.print("\n[bold]Interpretation:[/bold]")
        self.processor.console.print("• High capacitance at a pH means the protein is charge-responsive there")
        self.processor.console.print("• Maximum capacitance occurs when pH ≈ pKa values (greatest sensitivity)")
        self.processor.console.print("• Useful for identifying critical pH ranges for your system\n")

        # Prompt for pH directly
        pH = float(prompt_with_context(
            processor=self.processor,
            prompt=f"Enter pH for capacitance calculation",
            default=f"{self.analyzer.pH:.1f}",
            module="Protonation State Analyzer",
            description="Enter pH for capacitance"
        ))

        # Calculate capacitance
        capacitance, cap_per_group, specific_capacitance, total_titratable = self.analyzer.calculate_ideal_capacitance(pH)

        # Display results
        self._display_capacitance_results(capacitance, cap_per_group, specific_capacitance, total_titratable, pH)

        # Store in workspace
        self.update_workspace("ideal_capacitance", capacitance)
        self.update_workspace("specific_capacitance", specific_capacitance)
        self.update_workspace("capacitance_per_group", cap_per_group)

        return True  # Success

    def _calculate_capacitance_for_microstates(self):
        """Calculate capacitance for selected microstates"""
        console = self.processor.console
        microstate_results = self.get_from_workspace("microstate_protonation_results", {})

        if not microstate_results:
            console.print("[yellow]No microstate results available.[/yellow]")
            return

        # Display table of available microstates
        from rich.table import Table
        table = Table(title="Analyzed Microstates")
        table.add_column("#", style="cyan", justify="right")
        table.add_column("Filename", style="green")
        table.add_column("pH", style="yellow")
        table.add_column("Net Charge", style="magenta")

        microstate_list = list(microstate_results.items())
        for i, (filepath, data) in enumerate(microstate_list, 1):
            table.add_row(
                str(i),
                data['filename'],
                f"{data['pH']:.1f}",
                f"{data['net_charge']:.2f}"
            )

        console.print(table)

        # Ask user to select which ones to calculate

        console.print("\n[bold]Selection Options:[/bold]")
        console.print("  • Single number: e.g., '3'")
        console.print("  • Multiple: e.g., '1,3,5'")
        console.print("  • Range: e.g., '1-4'")
        console.print("  • All: enter 'all'")

        selection = prompt_with_context(
            processor=self.processor,
            prompt="Select microstate(s) to calculate capacitance",
            default="all",
            module="Protonation State Analyzer",
            description="Microstate selection for capacitance"
        )

        # Parse selection
        selected_indices = self._parse_microstate_selection(selection, len(microstate_list))

        if not selected_indices:
            console.print("[yellow]No valid microstates selected.[/yellow]")
            return

        # Ask for pH once (applies to all)
        # Check if we can use the analysis pH
        first_microstate = microstate_list[selected_indices[0]][1]
        analysis_pH = first_microstate.get('pH', 7.0)

        use_analysis_pH = confirm_with_context(
            processor=self.processor,
            prompt=f"Calculate capacitance at analysis pH ({analysis_pH:.1f})?",
            default=True,
            module="Protonation State Analyzer",
            description="Use analysis pH for capacitance"
        )

        if use_analysis_pH:
            pH = analysis_pH
        else:
            pH = float(prompt_with_context(
                processor=self.processor,
                prompt="Enter pH for capacitance calculation",
                default="7.0",
                module="Protonation State Analyzer",
                description="Enter pH for capacitance"
            ))

        # Calculate capacitance for each selected microstate
        for idx in selected_indices:
            filepath, data = microstate_list[idx]

            console.print(f"\n[bold cyan]{'='*60}[/bold cyan]")
            console.print(f"[bold]Microstate: {data['filename']}[/bold]")
            console.print(f"[bold cyan]{'='*60}[/bold cyan]\n")

            # Temporarily set the analysis results for this microstate
            original_results = self.analysis_results
            self.analysis_results = data['results']

            # Calculate capacitance (analyzer needs to have the right titratable_residues)
            # We need to ensure the analyzer's internal state is correct
            with console.status("[bold green]Calculating capacitance..."):
                capacitance, cap_per_group, specific_capacitance, total_titratable = self.analyzer.calculate_ideal_capacitance(pH)

            # Display results
            self._display_capacitance_results(capacitance, cap_per_group, specific_capacitance, total_titratable, pH)

            # Store results back in microstate data
            if 'capacitance_data' not in data:
                data['capacitance_data'] = {}
            data['capacitance_data'][pH] = {
                'capacitance': capacitance,
                'cap_per_group': cap_per_group,
                'specific_capacitance': specific_capacitance,
                'total_titratable': total_titratable
            }

            # Restore original results
            self.analysis_results = original_results

        # Update workspace with the modified microstate data
        workspace = self.get_workspace()
        for idx in selected_indices:
            filepath, data = microstate_list[idx]
            microstate_results[filepath] = data
        self.update_workspace("microstate_protonation_results", microstate_results)

        console.print(f"\n[bold green]✓ Calculated capacitance for {len(selected_indices)} microstate(s)[/bold green]")
        
    def generate_capacitance_profile(self):
        """Generate capacitance profile over pH range"""
        if not self.analyzer.titratable_residues:
            self.processor.console.print(
                "[yellow]No titratable residues analyzed. Run analysis first.[/yellow]"
            )
            return

        # Prompt for parameters
        min_pH = float(prompt_with_context(
            processor=self.processor,
            prompt="Enter minimum pH",
            default="0",
            module="Protonation State Analyzer",
            description="Enter minimum pH for capacitance profile"
        ))
        max_pH = float(prompt_with_context(
            processor=self.processor,
            prompt="Enter maximum pH",
            default="14",
            module="Protonation State Analyzer",
            description="Enter maximum pH for capacitance profile"
        ))
        step = float(prompt_with_context(
            processor=self.processor,
            prompt="Enter pH step size",
            default="0.1",
            module="Protonation State Analyzer",
            description="Enter pH step size for capacitance profile"
        ))

        with self.processor.console.status(
            "[bold green]Generating capacitance profile..."
        ):
            profile = self.analyzer.generate_capacitance_profile((min_pH, max_pH), step)
        
        # Display summary
        self.processor.console.print(
            f"\n[bold cyan]Capacitance Profile Summary:[/bold cyan]"
        )
        self.processor.console.print(
            f"pH range: {min_pH:.1f} - {max_pH:.1f}"
        )
        self.processor.console.print(
            f"Maximum total capacitance: {profile['max_capacitance']:.4f} at pH {profile['pH_max_capacitance']:.2f}"
        )
        self.processor.console.print(
            f"Maximum specific capacitance: {profile['max_specific_capacitance']:.4f} at pH {profile['pH_max_specific_capacitance']:.2f}"
        )
        
        # Generate and display the ASCII plot
        self.processor.console.print("\n[bold]Capacitance Curve:[/bold]")
        ascii_plot = self._plot_ascii_capacitance_curve(profile)
        
        # Print the plot with a monospace font style
        from rich.panel import Panel
        self.processor.console.print(
            Panel(
                f"[font=monospace]{ascii_plot}[/font]",
                title="Ideal Charge Capacitance",
                border_style="bright_blue",
                expand=False
            )
        )
        
        # Store in workspace
        self.update_workspace("capacitance_profile", profile)

        # Add interactive lookup
        if confirm_with_context(
            processor=self.processor,
            prompt="\nWould you like to do interactive pH/capacitance lookups?",
            default=False,
            module="Protonation State Analyzer",
            description="Interactive pH/capacitance lookup"
        ):
            self._interactive_capacitance_lookup(profile)

        # Offer to save or plot
        if confirm_with_context(
            processor=self.processor,
            prompt="Would you like to save the capacitance profile?",
            module="Protonation State Analyzer",
            description="Save capacitance profile"
        ):
            self._save_capacitance_profile(profile)
        
        return True  # Success
        
    def _display_capacitance_results(self, capacitance, cap_per_group, specific_capacitance, total_titratable, pH):
        """Display capacitance results in a formatted table"""
        from rich.table import Table
        from rich.panel import Panel
        from rich.align import Align
        
        # Display the equation in a more mathematical style
        equation_lines = [
            "[bold cyan]Ideal Charge Capacitance Equation:[/bold cyan]",
            "",
            "       [bold]C[/bold][italic]ideal[/italic] = Σᵢ [10^(pH-pKᵢ)] / [1 + 10^(pH-pKᵢ)]²",
            "",
            "where:",
            "  i = individual titratable residue",
            "  pKᵢ = pKa for residue i (from PROPKA or default values)",
            f"  pH = [bold green]{pH:.1f}[/bold green] (current pH)",
        ]
        
        equation_panel = Panel(
            Align.center("\n".join(equation_lines)),
            title="Ideal Charge Capacitance Theory",
            border_style="blue",
            width=80
        )
        
        self.processor.console.print(equation_panel)
        self.processor.console.print()
        
        # Create summary table
        table = Table(title=f"Ideal Charge Capacitance at pH {pH:.1f}")
        table.add_column("Residue Type", style="cyan")
        table.add_column("Count", justify="right")
        table.add_column("pKa Range", justify="center")
        table.add_column("Avg pKa", justify="right")
        table.add_column("Total C", justify="right")
        table.add_column("C/Residue", justify="right")

        # Check if we're using custom pKa values (from PROPKA)
        using_propka = hasattr(self.analyzer, 'custom_pkas') and self.analyzer.custom_pkas
        
        # Add rows for each residue type
        for res_type, data in sorted(cap_per_group.items()):
            # Format pKa range
            if "min_pKa" in data and "max_pKa" in data and data["min_pKa"] != data["max_pKa"]:
                pka_range = f"{data['min_pKa']:.1f}-{data['max_pKa']:.1f}"
            else:
                pka_range = f"{data.get('avg_pKa', data.get('pKa', 0)):.1f}"
            
            # Format average pKa with std if available
            if "pKa_std" in data and data["pKa_std"] > 0.1:
                avg_pka_str = f"{data['avg_pKa']:.1f}±{data['pKa_std']:.1f}"
            else:
                avg_pka_str = f"{data.get('avg_pKa', data.get('pKa', 0)):.1f}"
            
            table.add_row(
                res_type,
                str(data["count"]),
                pka_range,
                avg_pka_str,
                f"{data['capacitance']:.4f}",
                f"{data['capacitance_per_residue']:.4f}"
            )
        
        # Add total row
        table.add_row(
            "[bold]TOTAL[/bold]",
            f"[bold]{total_titratable}[/bold]",
            "",
            "",
            f"[bold]{capacitance:.4f}[/bold]",
            "",
            style="bold green"
        )
        
        self.processor.console.print(table)
        
        # Combined summary and key insights
        self.processor.console.print("\n[bold]Summary & Key Insights:[/bold]")
        self.processor.console.print(f"• Total titratable residues: {total_titratable}")
        self.processor.console.print(f"• Total ideal charge capacitance: {capacitance:.4f}")
        self.processor.console.print(f"• Specific capacitance: {specific_capacitance:.4f} = {capacitance:.4f} / {total_titratable} residues")
        self.processor.console.print(f"  → At pH {pH:.1f}, charge-responsiveness is {self._interpret_capacitance_level(specific_capacitance)}")

        if using_propka:
            self.processor.console.print(f"• pKa source: [cyan]PROPKA predictions[/cyan] (environment-adjusted for this structure)")

            # Check for significant pKa shifts
            significant_shifts = []
            for res_type, data in cap_per_group.items():
                if "avg_pKa" in data and res_type in PKA_VALUES:
                    shift = abs(data["avg_pKa"] - PKA_VALUES[res_type])
                    if shift > 1.0:
                        significant_shifts.append(f"{res_type} (Δ={shift:.1f})")

            if significant_shifts:
                self.processor.console.print(f"  → Significant pKa shifts from textbook values: {', '.join(significant_shifts)}")
        else:
            self.processor.console.print(f"• pKa source: [cyan]Textbook values[/cyan] (Lehninger - not adjusted for structural environment)")

    def _interpret_capacitance_level(self, specific_capacitance):
        """Provide qualitative interpretation of specific capacitance value"""
        # Theoretical maximum for a single residue is 0.25 (at pH = pKa)
        # Typical ranges based on how many residues are near their pKa
        if specific_capacitance > 0.15:
            return "[bold green]very high[/bold green] (many residues near pKa)"
        elif specific_capacitance > 0.10:
            return "[bold cyan]high[/bold cyan] (good buffering capacity)"
        elif specific_capacitance > 0.05:
            return "[yellow]moderate[/yellow] (some buffering)"
        else:
            return "[yellow]low[/yellow] (most residues far from pKa)"

    def _save_capacitance_profile(self, profile):
        """Save capacitance profile to file"""
        import csv

        filename = prompt_with_context(
            processor=self.processor,
            prompt="Enter filename for capacitance profile",
            default="capacitance_profile.csv",
            module="Protonation State Analyzer",
            description="Enter capacitance profile filename"
        )
        
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            
            # Write header
            header = ["pH", "Total_Capacitance", "Specific_Capacitance"]  # UPDATE THIS
            for res_type in sorted(profile["group_capacitances"].keys()):
                header.append(f"{res_type}_Capacitance")
            writer.writerow(header)

            # Write data
            for i, pH in enumerate(profile["pH_values"]):
                row = [pH, profile["total_capacitance"][i], profile["specific_capacitance"][i]]  # UPDATE THIS
                for res_type in sorted(profile["group_capacitances"].keys()):
                    row.append(profile["group_capacitances"][res_type][i])
                writer.writerow(row)
                        
        self.processor.console.print(f"[green]Capacitance profile saved to {filename}[/green]")

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # RedoxSite Integration Methods
    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#

    def _get_redox_site_residues_for_exclusion(self):
        """
        Get RedoxSite residues from workspace and allow user selection for exclusion.
        
        Returns:
            set: Set of (chain, resid, icode) tuples for excluded residues
        """
        # Check for detected redox sites in workspace
        detected_redox_sites = self.get_from_workspace("detected_redox_sites", [])
        
        if not detected_redox_sites:
            return set()
            
        self.processor.console.print("\n[bold cyan]Redox Sites Detected[/bold cyan]")
        self.processor.console.print("The following redox-active sites were found in your structure:")
        
        # Extract and display residues for user selection
        residue_catalog = self._extract_residues_from_redox_sites(detected_redox_sites)
        
        if not residue_catalog:
            self.processor.console.print("[yellow]No residues found in redox sites[/yellow]")
            return set()
            
        # Display sites and get user selection
        selected_residues = self._display_redox_sites_for_exclusion(residue_catalog)
        
        return selected_residues

    def _extract_residues_from_redox_sites(self, redox_sites):
        """
        Extract residue information from RedoxSite objects.
        
        Args:
            redox_sites: List of RedoxSite objects
            
        Returns:
            dict: Organized residue information by site
        """
        residue_catalog = {}
        
        for site in redox_sites:
            site_id = getattr(site, 'site_id', 'Unknown')
            site_type = getattr(site, 'site_type', 'Unknown')
            
            # Use residue_groups for efficient residue extraction (already grouped by residue)
            residue_groups = getattr(site, 'residue_groups', {})
            
            site_residues = []
            
            if residue_groups:
                # residue_groups is Dict[Tuple[str, int, str], List[coords]]
                # where key is (chain, resid, insertion_code)
                for residue_key, coords_list in residue_groups.items():
                    chain_id, res_num, icode = residue_key
                    
                    # Get residue name from first atom in this residue
                    res_name = 'UNK'  # Default
                    if coords_list and hasattr(site, 'coord_to_pdb'):
                        first_coord = coords_list[0]
                        atom_info = site.coord_to_pdb.get(first_coord, {})
                        res_name = atom_info.get('resname', 'UNK')
                    
                    site_residues.append({
                        'chain_id': chain_id,
                        'res_num': res_num,
                        'res_name': res_name,
                        'icode': icode,
                        'residue_key': residue_key
                    })
            else:
                # Fallback: extract from atoms list if residue_groups is empty
                atoms = getattr(site, 'atoms', [])
                residue_set = set()  # To avoid duplicates
                
                for atom in atoms:
                    # Extract residue information from atom
                    chain_id = getattr(atom, 'chain', '')
                    res_num = getattr(atom, 'resid', None)
                    res_name = getattr(atom, 'resname', 'UNK')
                    icode = getattr(atom, 'insertion_code', '') or ''
                    
                    if res_num is not None:
                        residue_key = (chain_id, res_num, icode)
                        if residue_key not in residue_set:
                            residue_set.add(residue_key)
                            site_residues.append({
                                'chain_id': chain_id,
                                'res_num': res_num,
                                'res_name': res_name,
                                'icode': icode,
                                'residue_key': residue_key
                            })
            
            if site_residues:
                residue_catalog[site_id] = {
                    'site_type': site_type,
                    'residues': site_residues
                }
                
        return residue_catalog

    def _display_redox_sites_for_exclusion(self, residue_catalog):
        """
        Display redox site residues and handle user selection.
        
        Args:
            residue_catalog: Dictionary of site information
            
        Returns:
            set: Set of (chain, resid, icode) tuples for selected residues
        """
        from rich.table import Table
        from rich.prompt import Prompt
        
        # Create display table
        table = Table(title="Redox Site Residues Available for Exclusion")
        table.add_column("Index", style="cyan", width=8)
        table.add_column("Site ID", style="green", width=12)
        table.add_column("Site Type", style="yellow", width=15)
        table.add_column("Chain", style="blue", width=8)
        table.add_column("Residue", style="white", width=15)
        table.add_column("Res ID", style="white", width=10)
        
        # Build indexed list of all residues
        residue_index = []
        index = 1
        previous_site_id = None
        
        for site_id, site_info in residue_catalog.items():
            site_type = site_info['site_type']
            
            # Add separator row between sites
            if previous_site_id is not None and previous_site_id != site_id:
                table.add_section()  # This adds a horizontal divider
            
            for residue in site_info['residues']:
                chain_id = residue['chain_id']
                res_name = residue['res_name']
                res_num = residue['res_num']
                icode = residue['icode']
                
                # Format residue number with icode if present (no colon for empty icode)
                if icode:
                    res_display = f"{res_num}{icode}"
                else:
                    res_display = str(res_num)
                
                table.add_row(
                    str(index),
                    site_id,
                    site_type,
                    chain_id,
                    res_name,
                    res_display
                )
                
                residue_index.append(residue['residue_key'])
                index += 1
            
            previous_site_id = site_id
        
        # Display table
        self.processor.console.print(table)

        # Highlight the listed residues in the viewer, palette-coloured
        # by parent redox site so the user can think "exclude all of
        # heme A, but keep heme B's residues" before reading individual
        # indices off the table. Stable per-site labels; cleaned up
        # after the user enters their selection.
        try:
            from proprep.structure_prep.viewer_coordinator import (
                viewer as _viewer,
            )
            applied_exclusion_labels: List[str] = []
            for site_idx, (site_id, site_info) in enumerate(
                residue_catalog.items(), 1
            ):
                pairs = set()
                for residue in site_info['residues']:
                    chain_id = residue.get('chain_id')
                    res_num = residue.get('res_num')
                    if chain_id and res_num is not None:
                        pairs.add((chain_id, res_num))
                if not pairs:
                    continue
                clauses = [f"(:{c} and {r})" for c, r in sorted(pairs)]
                lbl = f"prot_redox_excl_{site_id}"
                _viewer.highlight(
                    " or ".join(clauses),
                    style="ball+stick",
                    color=f"palette:{site_idx}",
                    label=lbl,
                )
                applied_exclusion_labels.append(lbl)
            self._exclusion_viewer_labels = applied_exclusion_labels
        except Exception as exc:
            logger.debug("redox exclusion viewer hook silenced: %s", exc)
            self._exclusion_viewer_labels = []

        # Get user selection
        self.processor.console.print(
            "\n[bold]Instructions:[/bold] Enter comma-separated indices of residues to "
            "[red]EXCLUDE[/red] from protonation analysis."
        )
        self.processor.console.print(
            "These residues will be skipped when analyzing protonation states."
        )
        self.processor.console.print(
            "Enter 'all' to exclude all redox site residues, or press Enter to include all."
        )
        self.processor.console.print(
            "\n[yellow]Note:[/yellow] Residues coordinated to metal centers typically exhibit "
            "unusually shifted pKa values."
        )
        self.processor.console.print(
            "[grey50]Metal coordination stabilizes specific charge states, causing large pKa shifts "
            "(often >2 pH units).[/grey50]"
        )
        self.processor.console.print(
            "[grey50]If you include metal-coordinating residues, expect PROPKA to predict shifted pKa values.[/grey50]\n"
        )

        selection = prompt_with_context(
            processor=self.processor,
            prompt="\nEnter indices to exclude (e.g., 1,3,5-7 or 'all')",
            default="",
            module="Protonation State Analyzer",
            description="Enter redox residue indices to exclude"
        ).strip().lower()

        try:
            if not selection:
                self.processor.console.print("[green]No residues excluded - all redox site residues will be analyzed[/green]")
                return set()

            # Handle 'all' option
            if selection == "all":
                all_residues = set(residue_index)
                self.processor.console.print(
                    f"[green]Excluding all {len(all_residues)} redox site residues from protonation analysis[/green]"
                )
                return all_residues

            # Parse selection
            try:
                selected_indices = self._parse_index_selection(selection, len(residue_index))
                selected_residues = set()

                for idx in selected_indices:
                    selected_residues.add(residue_index[idx - 1])  # Convert to 0-based

                self.processor.console.print(
                    f"[green]Selected {len(selected_residues)} residues for exclusion from protonation analysis[/green]"
                )

                return selected_residues

            except ValueError as e:
                self.processor.console.print(f"[red]Error parsing selection: {e}[/red]")
                self.processor.console.print("[yellow]No residues will be excluded[/yellow]")
                return set()
        finally:
            # Clear the per-site exclusion highlights now that the
            # picker has resolved — they don't carry meaning beyond
            # the selection prompt.
            self._clear_exclusion_viewer()

    def _clear_exclusion_viewer(self) -> None:
        """Drop the per-site highlights set by the redox-exclusion picker."""
        labels = getattr(self, "_exclusion_viewer_labels", None) or []
        if not labels:
            return
        try:
            from proprep.structure_prep.viewer_coordinator import (
                viewer as _viewer,
            )
            for lbl in labels:
                _viewer.unhighlight(lbl)
        except Exception as exc:
            logger.debug("redox exclusion cleanup silenced: %s", exc)
        self._exclusion_viewer_labels = []

    def _parse_index_selection(self, selection, max_index):
        """
        Parse comma-separated index selection with range support.
        
        Args:
            selection: String like "1,3,5-7,10"
            max_index: Maximum valid index
            
        Returns:
            list: List of selected indices
        """
        indices = set()
        
        for part in selection.split(','):
            part = part.strip()
            
            if '-' in part:
                # Handle range like "5-7"
                try:
                    start, end = part.split('-', 1)
                    start_idx = int(start.strip())
                    end_idx = int(end.strip())
                    
                    if start_idx < 1 or end_idx > max_index or start_idx > end_idx:
                        raise ValueError(f"Invalid range {part} (valid range: 1-{max_index})")
                    
                    for i in range(start_idx, end_idx + 1):
                        indices.add(i)
                        
                except ValueError:
                    raise ValueError(f"Invalid range format: {part}")
            else:
                # Handle single index
                try:
                    idx = int(part)
                    if idx < 1 or idx > max_index:
                        raise ValueError(f"Index {idx} out of range (valid range: 1-{max_index})")
                    indices.add(idx)
                except ValueError:
                    raise ValueError(f"Invalid index: {part}")
        
        return sorted(list(indices))
        
