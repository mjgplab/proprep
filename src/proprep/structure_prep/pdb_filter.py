"""
PDB Filter Module

Module for filtering PDB structures by component type with interface analysis.
Integrates with the MPSA processor workflow using the command pattern.
"""

import copy
import json
import logging
from collections import defaultdict
from typing import Any, Dict, List

from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from proprep.utils.module_registry import ProcessingModule, register_module
from proprep.utils.prompts import prompt_with_context, confirm_with_context
from proprep.structure_prep.viewer_coordinator import _is_web_shell_mode
from proprep.application.processor_command import ModuleActionCommand
from .pdb_filter_worker import PDBFilterWorker, ComponentClassifier
from .pdb_filter_commands import (
    FilterPDBStructureCommand,
    ShowFilterStatusCommand,
    ExportFilterStatisticsCommand,
)

# Import for redox detection integration
try:
    from .comprehensive_redox_detector import ComprehensiveRedoxDetector
    REDOX_DETECTION_AVAILABLE = True
except ImportError:
    REDOX_DETECTION_AVAILABLE = False
    logger.warning("Redox detection module not available. Filtering will proceed without redox site analysis.")

logger = logging.getLogger(__name__)


@register_module
class PDBFilterModule(ProcessingModule):
    """Module for filtering PDB structures"""

    NAME = "PDB Filter"
    DESCRIPTION = "Filter PDB structures by component"
    VERSION = "1.0.0"
    CATEGORY = "preparation"
    REQUIRES = ["PDB Loader"]
    PRIORITY = 2

    def initialize(self):
        """Initialize module resources"""
        self.filter_worker = None
        self.filtered_structure = None
        self.filter_selections = {}
        self._local_command_history = []
        self.skip_redox_detection = False

    @property
    def console(self):
        """Get console from processor or create new one"""
        if self.processor and hasattr(self.processor, 'console'):
            return self.processor.console
        else:
            from rich.console import Console
            return Console()

    def get_workspace_display(self, workspace_key, value, console):
        """
        Custom display method for PDB Filter module data

        Args:
            workspace_key: The key in the workspace dictionary
            value: The value stored in the workspace
            console: The rich console object for output

        Returns:
            bool: True if the module handled the display, False otherwise
        """
        if workspace_key == "filtered_structure" and value is not None:
            structure_table = Table(title="Filtered Structure Contents")
            structure_table.add_column("Component", style="cyan")
            structure_table.add_column("Details", style="green")
            structure_table.add_column("Count", style="yellow")

            model_count = len(value)
            structure_table.add_row(
                "Models", "Total models in structure", str(model_count)
            )

            for model_idx, model in enumerate(value):
                chain_count = len(model)
                total_residues = sum(len(chain) for chain in model)
                total_atoms = sum(
                    sum(len(residue) for residue in chain) for chain in model
                )

                structure_table.add_row(
                    f"Model {model_idx}",
                    f"{chain_count} chains, {total_residues} residues, {total_atoms} atoms",
                    "",
                )

                for chain in model:
                    chain_id = chain.id
                    residue_count = len(chain)
                    atom_count = sum(len(residue) for residue in chain)

                    structure_table.add_row(
                        f"Chain {chain_id}",
                        f"{residue_count} residues, {atom_count} atoms",
                        "",
                    )

                    residue_types = defaultdict(int)
                    hetero_count = 0
                    water_count = 0

                    for residue in chain:
                        if residue.id[0] != " ":
                            hetero_count += 1
                            if residue.resname in ("HOH", "WAT"):
                                water_count += 1
                        residue_types[residue.resname] += 1

                    if hetero_count > 0:
                        structure_table.add_row(
                            "",
                            f"Hetero groups: {hetero_count} (including {water_count} waters)",
                            "",
                        )

                    if len(residue_types) <= 10:
                        for res_name, count in sorted(
                            residue_types.items(), key=lambda x: x[1], reverse=True
                        ):
                            if count > 5:
                                structure_table.add_row(
                                    "", f"Residue {res_name}", str(count)
                                )
                    else:
                        top_types = sorted(
                            residue_types.items(), key=lambda x: x[1], reverse=True
                        )[:5]
                        other_count = sum(
                            count
                            for _, count in sorted(
                                residue_types.items(), key=lambda x: x[1], reverse=True
                            )[5:]
                        )

                        for res_name, count in top_types:
                            structure_table.add_row(
                                "", f"Residue {res_name}", str(count)
                            )
                        structure_table.add_row(
                            "", "Other residue types", str(other_count)
                        )

            console.print(structure_table)
            return True

        elif workspace_key == "filter_selections" and value is not None:
            console.print("[bold]Filter Selections[/bold]")

            overview_table = Table(title="Selection Overview")
            overview_table.add_column("Chain", style="cyan")
            overview_table.add_column("Component Types", style="green")
            overview_table.add_column("Residue Count", style="yellow")

            for chain_id, components in value.items():
                total_residues = 0
                component_names = []

                for comp_type, residues in components.items():
                    display_type = ComponentClassifier.display_name(comp_type)
                    component_names.append(display_type)

                    if isinstance(residues, (list, set)):
                        total_residues += len(residues)

                overview_table.add_row(
                    chain_id, ", ".join(component_names), str(total_residues)
                )

            console.print(overview_table)

            for chain_id, components in value.items():
                chain_table = Table(title=f"Chain {chain_id} Components")
                chain_table.add_column("Component Type", style="cyan")
                chain_table.add_column("Residue Count", style="green")
                chain_table.add_column("Residue Numbers", style="yellow")

                for comp_type, residues in components.items():
                    display_type = ComponentClassifier.display_name(comp_type)

                    if isinstance(residues, (list, set)):
                        count = len(residues)
                        sorted_residues = sorted(list(residues))

                        if len(sorted_residues) > 10:
                            ranges = []
                            start = sorted_residues[0]
                            prev = start

                            for i, num in enumerate(sorted_residues[1:], 1):
                                if num > prev + 1:
                                    if start == prev:
                                        ranges.append(str(start))
                                    else:
                                        ranges.append(f"{start}-{prev}")
                                    start = num
                                prev = num

                            if start == sorted_residues[-1]:
                                ranges.append(str(start))
                            else:
                                ranges.append(f"{start}-{sorted_residues[-1]}")

                            if len(ranges) > 5:
                                display_ranges = (
                                    ", ".join(ranges[:3])
                                    + f", ... ({len(ranges)-3} more ranges)"
                                )
                            else:
                                display_ranges = ", ".join(ranges)

                            residue_display = display_ranges
                        else:
                            residue_display = ", ".join(map(str, sorted_residues))
                    else:
                        count = 0
                        residue_display = "None"

                    chain_table.add_row(display_type, str(count), residue_display)

                console.print(chain_table)

            return True

        elif workspace_key == "filtered_pdb_file" and value is not None:
            console.print(f"[bold]Filtered PDB File:[/bold] {value}")
            return True

        return False

    def get_menu_options(self) -> Dict[str, str]:
        """Get module menu options"""
        return {
            "filter": "Filter PDB structure",
            "show_filter": "Show filtering status",
            "export_stats": "Export filter statistics",
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

        # Check if structure has been filtered
        filtered_structure = workspace.get("filtered_structure")
        has_filtered = filtered_structure is not None

        # Option 1: Filter PDB structure - needs a loaded structure
        if has_filtered:
            filter_status, filter_dep = OptionStatus.COMPLETED, ""
        elif self.can_process(workspace):
            filter_status, filter_dep = OptionStatus.AVAILABLE, ""
        else:
            filter_status = OptionStatus.BLOCKED
            filter_dep = self.availability_note(workspace) or "Load a structure first"
        options.append(MenuOption(
            key="1",
            description="Filter PDB structure",
            status=filter_status,
            dependency_text=filter_dep,
        ))

        # Option 2: Show filtering status - requires filtering to be done
        options.append(MenuOption(
            key="2",
            description="Show filtering status",
            status=OptionStatus.READY if has_filtered else OptionStatus.BLOCKED,
            dependency_text="[Need to filter structure first] ○" if not has_filtered else ""
        ))

        # Option 3: Export filter statistics - requires filtering to be done
        options.append(MenuOption(
            key="3",
            description="Export filter selection",
            status=OptionStatus.READY if has_filtered else OptionStatus.BLOCKED,
            dependency_text="[Need to filter structure first] ○" if not has_filtered else ""
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
        filtered_structure = workspace.get("filtered_structure")

        if not filtered_structure:
            if not self.can_process(workspace):
                return f"{self.availability_note(workspace) or 'A structure is required'}. Load one via the Structure Loader."
            return "Start by filtering the PDB structure (option 1) to select relevant components"
        else:
            return "View filter status with option 2, export the filter selection with option 3, or press [m] to return to the main menu"

    def handle_menu_option(self, option: str) -> bool:
        """Handle a menu option selection"""
        try:
            if option == "filter":
                command = FilterPDBStructureCommand(self.processor, interactive=True)
                command.execute_with_error_handling()
                return True
            elif option == "show_filter":
                command = ShowFilterStatusCommand(self.processor)
                command.execute_with_error_handling()
                return True
            elif option == "export_stats":
                command = ExportFilterStatisticsCommand(self.processor)
                command.execute_with_error_handling()
                return True
        except Exception as e:
            import traceback
            logger.error(f"Error executing menu option '{option}': {e}")
            logger.error(f"Traceback:\n{traceback.format_exc()}")
            self.processor.console.print(f"[red]Error: {e}[/red]")
            self.processor.console.print(f"[grey50]{traceback.format_exc()}[/grey50]")

        return False

    def filter_pdb_structure(
        self, workspace: Dict[str, Any] = None, interactive=True,
        skip_redox_detection=False
    ) -> Dict[str, Any]:
        """Filter PDB structure by component type.

        Args:
            workspace: Current workspace dict (defaults to processor workspace)
            interactive: Enable interactive filtering mode
            skip_redox_detection: Skip redox site detection entirely.
                Used when launching from structure_preprocessor, which
                runs its own redox detection on the final prepared structure.
        """
        if workspace is None:
            workspace = self.processor._get_workspace()

        pdb_file = self._get_pdb_file_from_workspace(workspace)
        if not pdb_file:
            self.processor.console.print(
                "[yellow]No PDB file loaded. Please load a PDB file first.[/yellow]"
            )
            return workspace

        existing_structure = self._get_structure_from_workspace(workspace, pdb_file)

        self.filter_worker = PDBFilterWorker(pdb_file, existing_structure, processor=self.processor)
        self.skip_redox_detection = skip_redox_detection

        # Store H++ detection flag in workspace for downstream modules
        if self.filter_worker.is_hplusplus_structure:
            workspace = self.update_workspace(workspace, "is_hplusplus_structure", True)
            # Display informative message about H++ structure
            self.processor.console.print(
                "\n[bold cyan]═══ H++ Structure Detected ═══[/bold cyan]\n"
                "[cyan]• Residue names:[/cyan] AMBER forcefield conventions (HID/HIE/HIP, CYX, etc.)\n"
                "[cyan]• Hydrogens:[/cyan] Already added by H++ server\n"
                "[cyan]• Analysis mode:[/cyan] Backbone-based residue classification\n"
                "[cyan]• Interface detection:[/cyan] Distance-based (freeSASA skipped)\n"
            )

        if interactive:
            result = self._run_interactive_filter()
        else:
            result = self._run_automated_filter()
            
        # Handle result (both interactive and automated now return tuple or None)
        if result and len(result) == 2:
            filtered_structure, filtered_pdb_file = result
        else:
            filtered_structure = result
            filtered_pdb_file = None

        if filtered_structure:
            self.filtered_structure = filtered_structure
            filter_selections_copy = copy.deepcopy(self.filter_worker.filter_selections)

            workspace = self.update_workspace(
                workspace, "filtered_structure", filtered_structure
            )
            if filtered_pdb_file:
                workspace = self.update_workspace(
                    workspace, "filtered_pdb_file", filtered_pdb_file
                )
            workspace = self.update_workspace(
                workspace, "filter_selections", filter_selections_copy
            )

            self.processor.console.print(
                "[green]Structure filtering completed successfully.[/green]"
            )

            # Sync redox sites with the filtered structure
            removed = self._sync_redox_sites_with_filtered_structure(workspace, filtered_structure)
            if removed > 0:
                self.processor.console.print(
                    f"[cyan]Synced redox sites: removed {removed} atom{'s' if removed != 1 else ''} "
                    f"no longer in filtered structure[/cyan]"
                )

        return workspace

    def show_filter_status(self, workspace: Dict[str, Any] = None) -> Dict[str, Any]:
        """Show the current filter status."""
        if workspace is None:
            workspace = self.processor._get_workspace()

        filtered_structure = workspace.get("filtered_structure")
        if not filtered_structure:
            self.processor.console.print(
                "[yellow]No filtered structure available.[/yellow]"
            )
            return workspace

        table = Table(title="Filtered Structure Status")
        table.add_column("Component", style="cyan")
        table.add_column("Count", style="green")

        chains = len(list(filtered_structure[0]))
        table.add_row("Chains", str(chains))

        residue_counts = {}
        for chain in filtered_structure[0]:
            for residue in chain:
                res_type = residue.id[0]
                if res_type not in residue_counts:
                    residue_counts[res_type] = 0
                residue_counts[res_type] += 1

        for res_type, count in residue_counts.items():
            type_name = "Standard" if res_type == " " else res_type
            table.add_row(f"Residues ({type_name})", str(count))

        self.processor.console.print(table)

        if self.filter_worker:
            stats = self.filter_worker.get_filter_statistics()
            if stats:
                self._display_filter_summary(stats)

        return workspace

    def export_filter_statistics(
        self, workspace: Dict[str, Any] = None
    ) -> Dict[str, Any]:
        """Export filter statistics to a JSON file."""
        if workspace is None:
            workspace = self.processor._get_workspace()

        if not self.filter_worker:
            self.processor.console.print(
                "[yellow]No filter worker available. Run filtering first.[/yellow]"
            )
            return workspace

        stats = self.filter_worker.get_filter_statistics()
        if not stats:
            self.processor.console.print(
                "[yellow]No filter statistics available.[/yellow]"
            )
            return workspace

        output_file = prompt_with_context(
            processor=self.processor,
            prompt="Enter output filename for statistics",
            default="filter_statistics.json",
            module="PDB Filter",
            description="Enter filename for filter statistics export"
        )

        try:
            with open(output_file, "w") as f:
                json.dump(stats, f, indent=2)

            self.processor.console.print(
                f"[green]Filter statistics exported to: {output_file}[/green]"
            )
        except Exception as e:
            self.processor.console.print(
                f"[red]Error exporting statistics: {str(e)}[/red]"
            )

        return workspace

    def _run_interactive_filter(self):
        """Run interactive filtering with user input."""
        # Check if structure has multiple models
        model_count = len(self.filter_worker.structure)

        selected_model_idx = self._get_model_selection()
        selected_model = self.filter_worker.structure[selected_model_idx]

        # Display BSA and topology analysis first
        self._display_initial_chain_analysis(selected_model_idx)

        # Redox site detection integration - after BSA/topology, before chain selection
        # When skip_redox_detection=True (called from structure_preprocessor), we still
        # check for pre-existing sites in workspace but don't prompt to run detection.
        detected_redox_sites = None

        # Always check for existing redox sites in workspace (from Redox Site Detector module)
        workspace = self.processor._get_workspace()
        existing_sites = workspace.get("detected_redox_sites")

        if existing_sites:
            # Sites already detected, use them for warnings
            num_sites = len(existing_sites)
            self.processor.console.print(f"\n[green]✓ Using {num_sites} previously detected redox site(s)[/green]")
            self.processor.console.print("[grey50]Redox sites were already detected earlier in the workflow[/grey50]")
            detected_redox_sites = existing_sites
        elif not self.skip_redox_detection:
            # Only prompt to run detection if not skipped
            if REDOX_DETECTION_AVAILABLE and confirm_with_context(
                processor=self.processor,
                prompt="\n[bold cyan]Detect redox-active sites before filtering?[/bold cyan]",
                default=True,
                module="PDB Filter",
                description="Detect redox-active sites before filtering"
            ):
                # Run detection using already-selected structure
                detected_redox_sites, detector = self._detect_redox_sites(self.filter_worker.filename, self.filter_worker.structure)
                if detected_redox_sites:
                    # Store redox sites in workspace for downstream modules
                    workspace = self.update_workspace(workspace, 'detected_redox_sites', detected_redox_sites)

                    # Store transformer mappings if available
                    if detector and hasattr(detector, 'transformer_mappings') and detector.transformer_mappings:
                        workspace = self.update_workspace(workspace, 'redox_transformer_mappings', detector.transformer_mappings)
                        self.processor.console.print(f"\n[green]✓ Stored transformer mappings for {len(detector.transformer_mappings)} site types[/green]")

                    self.processor.console.print("\n[green]✓ Redox sites stored in workspace for downstream processing[/green]")

        selected_chain_ids = self._get_chain_selection(selected_model_idx, detected_redox_sites, show_analysis=False)

        # Halo the chains the user just chose to keep so the about-to-
        # follow per-component prompts have visible context for "this
        # is what you're filtering". Stable label so a re-pick replaces.
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer
            chain_clause = " or ".join(f":{c}" for c in selected_chain_ids)
            if chain_clause:
                _viewer.highlight(
                    chain_clause,
                    style="halo",
                    color="#1f78b4",
                    label="pdbfilter_kept_chains",
                )
        except Exception:
            pass

        filter_selections = {}

        #Only store model selection for multi-model structures
        if model_count > 1:
            filter_selections["selected_model"] = selected_model_idx

        self.processor.console.print(
            "\n[bold underline blue]Chain Composition Overview[/bold underline blue]"
        )
        
        # Analyze composition once for all chains with redox integration
        chain_compositions = {}
        for chain_id in selected_chain_ids:
            chain = selected_model[chain_id]
            chain_compositions[chain_id] = self.filter_worker.analyze_chain_composition(chain, detected_redox_sites)

            self.processor.console.print(f"\n[bold blue]Chain {chain_id} Composition[/bold blue]")
            self._display_chain_composition(chain_compositions[chain_id], chain_id)

        # Confirm before proceeding with filtering
        self.processor.console.print()
        if not confirm_with_context(
            processor=self.processor,
            prompt="Do you want to proceed with filtering?",
            default=True,
            module="PDB Filter",
            description="Proceed with filtering"
        ):
            return None

        for chain_id in selected_chain_ids:
            chain = selected_model[chain_id]
            chain_selections = {}
            chain_composition = chain_compositions[chain_id]  # Use pre-computed composition
            
            # Check if this chain is involved in any redox sites
            chain_redox_involvement = self._check_chain_redox_involvement(chain_id, detected_redox_sites) if detected_redox_sites else None

            for comp_type, residue_counts in chain_composition.items():
                display_type = ComponentClassifier.display_name(comp_type)

                self.processor.console.print(
                    f"\n[bold blue]Chain {chain_id} - {display_type}[/bold blue]"
                )

                # Halo the current (chain, component) being prompted about
                # so the user sees what they're keeping/discarding. Stable
                # label so each iteration replaces the previous component
                # rather than accumulating halos across the loop.
                try:
                    from proprep.structure_prep.viewer_coordinator import viewer as _viewer
                    _component_selectors = {
                        "amino_acid":     "protein",
                        "water":          "water",
                        "ion":            "ion",
                        "dna_base":       "nucleic",
                        "rna_base":       "nucleic",
                        "hetero":         "hetero and not water and not ion",
                        "small_molecule": "hetero and not water and not ion",
                    }
                    sel_clause = _component_selectors.get(comp_type, "hetero")
                    _viewer.highlight(
                        f":{chain_id} and ({sel_clause})",
                        style="halo",
                        color="#ff7f00",
                        label="pdbfilter_current_component",
                    )
                except Exception:
                    pass

                # Check for redox site involvement before offering choices
                redox_conflict = self._check_redox_conflict_for_component_integrated(comp_type, chain_composition)
                
                if redox_conflict:
                    self.processor.console.print(f"[yellow]⚠️  This component contains redox site residues:[/yellow]")
                    for site_info in redox_conflict:
                        self.processor.console.print(f"   • Site {site_info['site_id']}: {site_info['residues']}")
                    
                    # Enhanced options for components with redox involvement
                    if comp_type == "amino_acid":
                        # Amino acid specific options
                        self.processor.console.print("\\[r] Retain entire component")
                        self.processor.console.print("\\[c] Keep redox site residues")
                        self.processor.console.print("\\[f] Keep redox site + and - N flanking residues")
                        self.processor.console.print("\\[s] Select specific residues")
                        self.processor.console.print("\\[o] Override redox protection and discard (⚠️ breaks redox sites)")
                        retention_choice = prompt_with_context(
                            processor=self.processor,
                            prompt="Choose option",
                            choices=["r", "c", "f", "s", "o"],
                            default="r",
                            module="PDB Filter",
                            description=f"Chain {chain_id} - {display_type} retention option (with redox sites)",
                            options_map={
                                "r": "Retain entire component",
                                "c": "Keep redox site residues",
                                "f": "Keep redox site + flanking residues",
                                "s": "Select specific residues",
                                "o": "Override redox protection and discard"
                            }
                        )
                    else:
                        # Non-amino acid options
                        self.processor.console.print("\\[r] Retain entire component")
                        self.processor.console.print("\\[c] Keep redox site residues")
                        self.processor.console.print("\\[s] Select specific residues")
                        self.processor.console.print("\\[o] Override redox protection and discard (⚠️ breaks redox sites)")
                        retention_choice = prompt_with_context(
                            processor=self.processor,
                            prompt="Choose option",
                            choices=["r", "c", "s", "o"],
                            default="r",
                            module="PDB Filter",
                            description=f"Chain {chain_id} - {display_type} retention option (with redox sites)",
                            options_map={
                                "r": "Retain entire component",
                                "c": "Keep redox site residues",
                                "s": "Select specific residues",
                                "o": "Override redox protection and discard"
                            }
                        )
                else:
                    # Standard options for components without redox involvement
                    self.processor.console.print("\\[r] Retain entire component")
                    self.processor.console.print("\\[s] Select specific residues")
                    self.processor.console.print("\\[d] Discard entire component")
                    retention_choice = prompt_with_context(
                        processor=self.processor,
                        prompt="Choose option",
                        choices=["r", "s", "d"],
                        default="r",
                        module="PDB Filter",
                        description=f"Chain {chain_id} - {display_type} retention option",
                        options_map={
                            "r": "Retain entire component",
                            "s": "Select specific residues",
                            "d": "Discard entire component"
                        }
                    )

                if retention_choice in ["d", "o"]:
                    if redox_conflict:
                        self.processor.console.print("[red]⚠️  Warning: Discarding this component will break redox sites![/red]")
                    continue

                if retention_choice == "s":
                    selected_residues = self._filter_component_type(chain, comp_type)
                    chain_selections[comp_type] = selected_residues
                elif retention_choice == "c":
                    # Keep only redox site residues
                    redox_residues = self._get_redox_site_residues_for_component_integrated(comp_type, chain_composition)
                    chain_selections[comp_type] = redox_residues
                    self.processor.console.print(f"[green]Keeping {len(redox_residues)} redox site residues[/green]")
                elif retention_choice == "f":
                    # Keep redox site + flanking residues
                    flanking_count = int(prompt_with_context(
                        processor=self.processor,
                        prompt="Number of flanking residues on each side",
                        default="2",
                        module="PDB Filter",
                        description="Enter number of flanking residues"
                    ))
                    redox_residues = self._get_redox_site_residues_for_component_integrated(comp_type, chain_composition)
                    flanking_residues = self._get_flanking_residues(chain, redox_residues, flanking_count)
                    chain_selections[comp_type] = flanking_residues
                    self.processor.console.print(f"[green]Keeping {len(flanking_residues)} residues (redox site + flanking)[/green]")

                    # Show what flanking N actually grabs — orange ball+stick
                    # on the (redox + flanking) set so the user sees the
                    # extended retention before moving on. Distinct color
                    # from the redox-at-risk halo (red) so the meaning is
                    # 'this is what you'll keep' rather than 'this is at
                    # risk'.
                    try:
                        from proprep.structure_prep.viewer_coordinator import viewer as _viewer
                        clauses = [f"(:{chain_id} and {r})" for r in sorted(flanking_residues)]
                        if clauses:
                            _viewer.highlight(
                                " or ".join(clauses),
                                style="ball+stick",
                                color="#ff7f00",
                                label="pdbfilter_flanking_kept",
                            )
                    except Exception:
                        pass
                else:
                    # Default: retain all - use existing component classification
                    chain_selections[comp_type] = {
                        residue.id[1]
                        for residue in chain
                        if ComponentClassifier.classify_residue(
                            residue, self.filter_worker.ccd_parser
                        )
                        == comp_type
                    }

            filter_selections[chain_id] = chain_selections

        self.filter_worker.filter_selections = filter_selections
        return self._review_selections(
            selected_model_idx, selected_chain_ids, filter_selections
        )

    def _run_automated_filter(self):
        """Run automated filtering with default settings."""
        if len(self.filter_worker.structure) == 0:
            return None

        model = self.filter_worker.structure[0]
        chain_ids = [chain.id for chain in model]

        filter_selections = {}
        for chain_id in chain_ids:
            chain = model[chain_id]
            chain_selections = {}

            for residue in chain:
                comp_type = ComponentClassifier.classify_residue(
                    residue, self.filter_worker.ccd_parser
                )
                if comp_type not in chain_selections:
                    chain_selections[comp_type] = set()
                chain_selections[comp_type].add(residue.id[1])

            filter_selections[chain_id] = chain_selections

        self.filter_worker.filter_selections = filter_selections
        filtered_structure = self.filter_worker.apply_filters(0, chain_ids, filter_selections)
        
        if filtered_structure:
            # Automatically save the filtered structure with standardized filename
            output_filename = "filtered_structure.pdb"
            self.filter_worker.save_filtered_structure(filtered_structure, output_filename)
            if hasattr(self, 'processor') and self.processor:
                self.processor.console.print(
                    f"[green]Filtered structure automatically saved to: {output_filename}[/green]"
                )
            return filtered_structure, output_filename
        
        return None

    def _get_model_selection(self) -> int:
        """Prompt user to select a model."""
        models = self.filter_worker.get_available_models()

        if len(models) == 1:
            return models[0]

        self.processor.console.print(
            "\n[bold underline]Available Models[/bold underline]"
        )
        for i, model_idx in enumerate(models):
            self.processor.console.print(f"{i+1}. Model {model_idx}")

        while True:
            # Build options map for model selection
            model_options = {str(i + 1): f"Model {models[i]}" for i in range(len(models))}

            choice = prompt_with_context(
                processor=self.processor,
                prompt="Select model",
                choices=[str(i + 1) for i in range(len(models))],
                default="1",
                module="PDB Filter",
                description="Select PDB model",
                options_map=model_options
            )
            picked = models[int(choice) - 1]
            if picked != models[0]:
                self.processor.console.print(
                    "[yellow]Note: viewer currently shows only the first model; "
                    "multi-model preview is on the backlog.[/yellow]"
                )
            return picked

    def _display_initial_chain_analysis(self, model_idx: int):
        """Display BSA and topology analysis before redox detection."""
        chain_info = self.filter_worker.get_model_chain_info(model_idx)

        self.processor.console.print(
            "\n[bold underline blue]Chain Interface Analysis[/bold underline blue]"
        )
        self.processor.console.print(
            "Calculating buried surface area between chains. This may take a moment..."
        )

        # Display the heatmap (which includes topology)
        self._display_chain_interface_heatmap(chain_info)

    def _get_chain_selection(self, model_idx: int, detected_redox_sites: List = None, show_analysis: bool = True) -> List[str]:
        """Prompt user to select chains in the given model."""
        chain_info = self.filter_worker.get_model_chain_info(model_idx)

        # Only show analysis if requested (default behavior for backward compatibility)
        if show_analysis:
            self.processor.console.print(
                "\n[bold underline blue]Chain Interface Analysis[/bold underline blue]"
            )
            self.processor.console.print(
                "Calculating buried surface area between chains. This may take a moment..."
            )
            # Display the heatmap (which includes topology)
            self._display_chain_interface_heatmap(chain_info)

        self.processor.console.print(
            "\n[bold underline blue]Available Chains[/bold underline blue]"
        )
#       chain_list = list(chain_info.keys())
        chain_list = [key for key in chain_info.keys() if not key.startswith('_')]

        for i, chain_id in enumerate(chain_list):
            info = chain_info[chain_id]

            # Check redox site involvement
            redox_display = self._get_chain_redox_display(chain_id, detected_redox_sites)

            if info["interfaces"]:
                sorted_interfaces = sorted(
                    [
                        (c, info["interface_areas"].get(c, 0))
                        for c in info["interfaces"]
                    ],
                    key=lambda x: x[1],
                    reverse=True,
                )
                interface_text = ", ".join(
                    [f"{c} ({int(area)}Å²)" for c, area in sorted_interfaces]
                )
                interface_display = f"[green]Interfaces with: {interface_text}[/green]"
            else:
                interface_display = "[dark_orange3]No interfaces detected[/dark_orange3]"

            self.processor.console.print(
                f"{i+1}. Chain {chain_id} - {info['residue_count']} residues - {interface_display}{redox_display}"
            )

        # Display heatmap only if analysis wasn't shown earlier
        if show_analysis:
            self._display_chain_interface_heatmap(chain_info)

        # Offer to view structure before selecting chains. Under proprep-web
        # the 3D viewer is already docked and open, so the Y/N prompt is
        # redundant — passively show the current structure in that panel
        # instead of asking. (show_structure no-ops if it's already the
        # bound structure, so this is cheap.) In the plain terminal there is
        # no viewer yet, so keep the explicit opt-in prompt.
        self.processor.console.print()
        if _is_web_shell_mode():
            self._launch_structure_viewer(force=False)
        elif confirm_with_context(
            processor=self.processor,
            prompt="Would you like to view the structure in the interactive 3D viewer?",
            default=False,
            module="PDB Filter",
            description="View structure before chain selection"
        ):
            self._launch_structure_viewer()

        while True:
            choice = prompt_with_context(
                processor=self.processor,
                prompt="Select chains (comma-separated indices, or 'all')",
                default="all",
                module="PDB Filter",
                description="Select chains to include in filtered structure"
            )

            if choice.lower() == "all":
                return chain_list

            try:
                selected_indices = [int(x.strip()) for x in choice.split(",")]

                for idx in selected_indices:
                    if idx < 1 or idx > len(chain_list):
                        self.processor.console.print(
                            f"[bold red]Invalid chain index: {idx}. Must be between 1 and {len(chain_list)}[/bold red]"
                        )
                        raise ValueError(f"Chain index out of range: {idx}")

                selected_chain_ids = [chain_list[idx - 1] for idx in selected_indices]
                
                # Check if excluding chains would break redox sites
                if detected_redox_sites:
                    excluded_chains = [chain_id for chain_id in chain_list if chain_id not in selected_chain_ids]
                    redox_conflicts = self._check_chain_exclusion_conflicts(excluded_chains, detected_redox_sites)
                    
                    if redox_conflicts:
                        self.processor.console.print(f"\n[yellow]⚠️  Warning: Excluding chains will affect redox sites:[/yellow]")
                        for conflict in redox_conflicts:
                            self.processor.console.print(f"   • {conflict}")

                        # Surface the at-risk redox residues spatially so the
                        # conflict isn't just text — user can see which sites
                        # break before choosing keep / add-chains / back. Red
                        # ball+stick on every residue of every affected site
                        # (any site with a center or atom in an excluded chain).
                        try:
                            from proprep.structure_prep.viewer_coordinator import viewer as _viewer
                            at_risk_residues = set()
                            excluded_set = set(excluded_chains)
                            for site in detected_redox_sites:
                                # Collect residues of this site
                                site_residues = set()
                                site_chains = set()
                                for c in getattr(site, 'centers', []) or []:
                                    site_chains.add(c.chain)
                                    site_residues.add((c.chain, c.resid))
                                for a in getattr(site, 'atoms', []) or []:
                                    site_chains.add(a.chain)
                                    site_residues.add((a.chain, a.resid))
                                # Site is at risk if any of its chains is being excluded
                                if site_chains & excluded_set:
                                    at_risk_residues.update(site_residues)
                            if at_risk_residues:
                                clauses = [f"(:{c} and {r})" for c, r in sorted(at_risk_residues)]
                                _viewer.highlight(
                                    " or ".join(clauses),
                                    style="ball+stick",
                                    color="#e31a1c",
                                    label="pdbfilter_redox_at_risk",
                                )
                        except Exception:
                            pass

                        self.processor.console.print("\n[bold]Options:[/bold]")
                        self.processor.console.print("\\[k] Keep current selection (breaks redox sites)")
                        self.processor.console.print("\\[a] Add required chains to maintain redox sites")
                        self.processor.console.print("\\[b] Go back and reselect chains")

                        conflict_choice = prompt_with_context(
                            processor=self.processor,
                            prompt="Choose option",
                            choices=["k", "a", "b"],
                            default="a",
                            module="PDB Filter",
                            description="Handle chain exclusion conflict with redox sites",
                            options_map={
                                "k": "Keep current selection (breaks redox sites)",
                                "a": "Add required chains to maintain redox sites",
                                "b": "Go back and reselect chains"
                            }
                        )
                        
                        if conflict_choice == "b":
                            continue  # Go back to chain selection
                        elif conflict_choice == "a":
                            # Add required chains
                            required_chains = self._get_required_chains_for_redox_sites(selected_chain_ids, detected_redox_sites)
                            selected_chain_ids.extend(required_chains)
                            selected_chain_ids = list(set(selected_chain_ids))  # Remove duplicates
                            self.processor.console.print(f"[green]Added chains {', '.join(required_chains)} to maintain redox sites[/green]")
                        # If "k", proceed with current selection
                
                return selected_chain_ids

            except ValueError:
                self.processor.console.print(
                    f"[bold red]Invalid input: {choice}. Please enter comma-separated numbers or 'all'[/bold red]"
                )

    def _display_chain_composition(
        self, chain_composition: Dict[str, Dict[str, Any]], chain_id: str
    ):
        """Display a hierarchical view of chain composition with redox site information."""
        # Colors chosen to stay legible on BOTH white and dark backgrounds
        # (manuscript figures + normal terminals): blue/green/dark_orange3
        # over the original cyan/yellow, which wash out on white.
        # header_style is set explicitly so the column headers don't fall
        # back to bold-default-foreground, which goes invisible on white
        # (see _display_chain_interface_heatmap for the full rationale).
        table = Table(
            title=f"Chain {chain_id} Composition", show_lines=False,
            header_style="bold blue",
        )
        # Color is spent on semantics only: blue for the grouping column,
        # dark_orange3 for the attention column. The two DATA columns use
        # the default foreground so they stay maximally legible on BOTH
        # white and dark backgrounds (saturated magenta/green clear contrast
        # on one background but not the other).
        # Cap the width of this column so long CCD chemical names (e.g.
        # "PROTOPORPHYRIN IX CONTAINING FE") wrap to several short lines
        # instead of stretching the column far wider than the others.
        table.add_column("Component Type", style="bold blue", max_width=18, overflow="fold")
        table.add_column("Residue\nName", style="default", justify="center")
        table.add_column("Count", style="default", justify="center")
        table.add_column("Redox\nSites", style="dark_orange3", justify="center")

        standard_types = ["amino_acid", "dna_base", "rna_base", "water"]

        standard_components = {
            t: chain_composition.get(t, {})
            for t in standard_types
            if t in chain_composition
        }

        special_components = {
            t: c for t, c in chain_composition.items() if t not in standard_types
        }

        for comp_type, residues in standard_components.items():
            # Handle both old format (int) and new format (dict with count/redox_residues)
            if residues and isinstance(list(residues.values())[0], dict):
                sorted_residues = sorted(residues.items(), key=lambda x: x[1]["count"], reverse=True)
                type_total = sum(data["count"] for _, data in sorted_residues)
            else:
                # Fallback for old format
                sorted_residues = sorted(residues.items(), key=lambda x: x[1], reverse=True)
                type_total = sum(count for _, count in sorted_residues)

            display_type = ComponentClassifier.display_name(comp_type)

            table.add_row(
                display_type, "[Total]", str(type_total), ""
            )

            for residue, data in sorted_residues:
                if isinstance(data, dict):
                    count = data["count"]
                    redox_residues = data.get("redox_residues", {})
                    if redox_residues:
                        redox_info = f"{len(redox_residues)} redox"
                    else:
                        redox_info = ""
                else:
                    # Fallback for old format
                    count = data
                    redox_info = ""
                
                table.add_row("", residue, str(count), redox_info)

            if sorted_residues:
                table.add_row("", "", "", "")

        for comp_type, residues in sorted(special_components.items()):
            # Handle both old format (int) and new format (dict with count/redox_residues)
            if residues and isinstance(list(residues.values())[0], dict):
                sorted_residues = sorted(residues.items(), key=lambda x: x[1]["count"], reverse=True)
                type_total = sum(data["count"] for _, data in sorted_residues)
            else:
                # Fallback for old format
                sorted_residues = sorted(residues.items(), key=lambda x: x[1], reverse=True)
                type_total = sum(count for _, count in sorted_residues)

            table.add_row(comp_type, "[Total]", str(type_total), "")

            for residue, data in sorted_residues:
                if isinstance(data, dict):
                    count = data["count"]
                    redox_residues = data.get("redox_residues", {})
                    if redox_residues:
                        redox_info = f"{len(redox_residues)} redox"
                    else:
                        redox_info = ""
                else:
                    # Fallback for old format
                    count = data
                    redox_info = ""
                
                table.add_row("", residue, str(count), redox_info)

            if sorted_residues:
                table.add_row("", "", "", "")

        self.processor.console.print(table)

    def _filter_component_type(self, chain, comp_type: str):
        """Filter specific component type within a chain with enhanced water analysis."""
        residues = self.filter_worker.get_component_residues(chain, comp_type)
        
        # NEW: Enhanced water analysis
        if comp_type == "water":
            return self.filter_worker.filter_water_with_analysis(
                chain, residues, console=self.processor.console
            )
        
        # ORIGINAL: Standard filtering for all other components (EXACT same logic)
        display_type = ComponentClassifier.display_name(comp_type)

        from rich.table import Table
        from rich.prompt import Prompt
        
        table = Table(title=f"Residues in Chain {chain.id} - {display_type}")
        table.add_column("Select", style="cyan")
        table.add_column("Record Type", style="yellow")
        table.add_column("Residue Name", style="magenta")
        table.add_column("Residue Number", style="green")
        
        for residue in residues:
            record_type = "Standard" if residue.id[0] == " " else "HETATM"
            table.add_row(
                str(residues.index(residue) + 1),
                record_type,
                residue.resname,
                str(residue.id[1]),
            )
            
        self.processor.console.print(table)

        choice = prompt_with_context(
            processor=self.processor,
            prompt="Select residues (comma-separated indices, 'all', or 'none')",
            default="all",
            module="PDB Filter",
            description=f"Select specific {display_type} residues to keep"
        )
        
        if choice.lower() == "all":
            picked_ids = {residue.id[1] for residue in residues}
        elif choice.lower() == "none":
            picked_ids = set()
        else:
            # ORIGINAL: Keep the exact same selection parsing logic
            selected_indices = [int(x.strip()) - 1 for x in choice.split(",")]
            picked_ids = {residues[idx].id[1] for idx in selected_indices}

        # Halo the picked residues so the user gets a visual confirmation
        # of what their numeric selection actually corresponds to. Blue
        # ball+stick distinguishes from the redox-at-risk (red) and
        # flanking (orange) reps elsewhere in the workflow.
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer
            if picked_ids:
                clauses = [f"(:{chain.id} and {r})" for r in sorted(picked_ids)]
                _viewer.highlight(
                    " or ".join(clauses),
                    style="ball+stick",
                    color="#1f78b4",
                    label="pdbfilter_specific_selection",
                )
            else:
                _viewer.unhighlight("pdbfilter_specific_selection")
        except Exception:
            pass

        return picked_ids

    def _review_selections(
        self,
        model_idx: int,
        chain_ids: List[str],
        filter_selections: Dict[str, Dict[str, set]],
    ):
        """Review and confirm filter selections."""
        selections_text = f"Model: {model_idx}\n\n"
        
        # Handle the filter_selections properly, skipping the 'selected_model' key
        for chain_id, chain_filters in filter_selections.items():
            # Skip the 'selected_model' key which contains an integer, not a dict
            if chain_id == "selected_model":
                continue
                
            selections_text += f"Chain {chain_id}:\n"
            if not chain_filters:
                selections_text += "  - No filters applied\n"
            else:
                for comp_type, residues in chain_filters.items():
                    selections_text += (
                        f"  - {comp_type.capitalize()}: {len(residues)} residues\n"
                    )

        self.processor.console.print(Panel(selections_text, title="Filter Selections", expand=False))

        # Visualize the entire retention set across all chains, coloured
        # by chain via the palette. Lets the user spot-check the final
        # composition spatially before answering 'apply these filters?'.
        # One ball+stick rep per chain so each chain is its own colour.
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer
            # Clear any leftover per-component / per-pick reps so the
            # final review isn't visually polluted by older overlays.
            _viewer.unhighlight("pdbfilter_current_component")
            _viewer.unhighlight("pdbfilter_specific_selection")
            _viewer.unhighlight("pdbfilter_flanking_kept")
            for ch_idx, (chain_id, chain_filters) in enumerate(
                (k, v) for k, v in filter_selections.items() if k != "selected_model"
            ):
                if not chain_filters:
                    _viewer.unhighlight(f"pdbfilter_review_{chain_id}")
                    continue
                kept_resids = set()
                for residue_set in chain_filters.values():
                    kept_resids.update(residue_set)
                if not kept_resids:
                    _viewer.unhighlight(f"pdbfilter_review_{chain_id}")
                    continue
                clauses = [f"(:{chain_id} and {r})" for r in sorted(kept_resids)]
                _viewer.highlight(
                    " or ".join(clauses),
                    style="ball+stick",
                    color=f"palette:{ch_idx + 1}",
                    label=f"pdbfilter_review_{chain_id}",
                )
        except Exception:
            pass

        if confirm_with_context(
            processor=self.processor,
            prompt="Do you want to apply these filters?",
            default=True,
            module="PDB Filter",
            description="Apply filter selections to structure"
        ):
            # Create a copy without the selected_model key for apply_filters
            chain_filters = {k: v for k, v in filter_selections.items() if k != "selected_model"}
            filtered_structure = self.filter_worker.apply_filters(
                model_idx, chain_ids, chain_filters
            )
            self.filter_worker.filtered_structure = filtered_structure

            # Create serializable version of selections
            serializable_selections = {}
            for chain_id, chain_data in filter_selections.items():
                if chain_id == "selected_model":
                    # Store the model selection as-is (it's already an integer)
                    serializable_selections[chain_id] = chain_data
                else:
                    # Convert sets to lists for chain data
                    serializable_selections[chain_id] = {}
                    for comp_type, residue_set in chain_data.items():
                        serializable_selections[chain_id][comp_type] = list(residue_set)

            self.filter_worker.filter_selections = serializable_selections

            # Automatically save the filtered structure with standardized filename
            output_filename = "filtered_structure.pdb"
            self.filter_worker.save_filtered_structure(
                filtered_structure, output_filename
            )
            self.processor.console.print(
                f"[green]Filtered structure automatically saved to: {output_filename}[/green]"
            )

            return filtered_structure, output_filename

        return None

    def _display_chain_interface_heatmap(self, chain_info: Dict[str, Dict[str, Any]]):
        """Display a visual heatmap of chain interfaces."""
#       chains = sorted(chain_info.keys())
        chains = sorted([key for key in chain_info.keys() if not key.startswith('_')])

        # header_style="bold blue": the chain-ID column headers (and the
        # row labels below) otherwise use Rich's default header style —
        # bold in the terminal's DEFAULT foreground, no color — which goes
        # invisible on a white background (macOS Terminal's "bright colors
        # for bold" remaps it to bright white). An explicit blue stays
        # legible on both white and dark.
        table = Table(
            title="Chain Interface Map (Buried Surface Area in Å²)",
            header_style="bold blue",
        )
        table.add_column("")
        for chain_id in chains:
            table.add_column(f"{chain_id}", justify="center")

        for chain1 in chains:
            row = [f"[bold blue]{chain1}[/bold blue]"]
            for chain2 in chains:
                if chain1 == chain2:
                    row.append("[grey50]■[/grey50]")
                elif chain2 in chain_info[chain1]["interfaces"]:
                    bsa = int(chain_info[chain1]["interface_areas"].get(chain2, 0))

                    if bsa > 2000:
                        row.append(f"[bold red]{bsa}[/bold red]")
                    elif bsa > 1000:
                        row.append(f"[bold green]{bsa}[/bold green]")
                    elif bsa > 500:
                        row.append(f"[green]{bsa}[/green]")
                    elif bsa > 200:
                        row.append(f"[dark_orange3]{bsa}[/dark_orange3]")
                    else:
                        row.append(f"[grey50]{bsa}[/grey50]")
                else:
                    row.append("[grey50]0[/grey50]")

            table.add_row(*row)

        self.processor.console.print(table)
        self.processor.console.print("\n[bold blue]BSA Color Legend:[/bold blue]")
        self.processor.console.print("[grey50]< 200 Å²[/grey50]: Minimal contact")
        self.processor.console.print("[dark_orange3]200-500 Å²[/dark_orange3]: Small interface")
        self.processor.console.print("[green]500-1000 Å²[/green]: Medium interface")
        self.processor.console.print(
            "[bold green]1000-2000 Å²[/bold green]: Large interface"
        )
        self.processor.console.print(
            "[bold red]> 2000 Å²[/bold red]: Very large interface"
        )

        topology_info = chain_info.get('_topology')
        if topology_info:
            self._display_topology_from_info(topology_info, chain_info)
            
    def _display_topology_from_info(self, topology_info: Dict[str, any], chain_info: Dict[str, Any] = None):
        """Render a streamlined chain-topology summary.

        A compact, label-aligned panel (chain/interface counts,
        connectivity with component count, a concise structure descriptor,
        and a degree histogram) followed by a one-line detection-criteria
        caption and the per-chain neighbor list. Replaces the older verbose
        panel + ASCII-art block.

        Degrees, component count, and the neighbor list are all derived from
        ``chain_info[chain]["interfaces"]`` — the same adjacency the
        interface heatmap uses — so the panel and the list cannot disagree.
        """
        from rich.panel import Panel
        from rich.text import Text

        # Adjacency (chain -> sorted neighbors), skipping private keys.
        adjacency = {}
        if chain_info:
            for chain, data in chain_info.items():
                if chain.startswith('_'):
                    continue
                adjacency[chain] = sorted(data.get('interfaces', []) or [])

        # Concise structure descriptor.
        structure_names = {
            'monomeric': 'monomeric',
            'linear': 'linear',
            'cyclic': 'cyclic',
            'complex_cyclic': 'cyclic, branched',
            'star_branched': 'star-branched',
            'multi_branched': 'branched',
            'disconnected': 'disconnected',
        }
        topology_type = topology_info.get('topology_type', 'unknown')
        structure_str = structure_names.get(topology_type, topology_type)

        # Prefer the displayed adjacency for all counts so the panel and the
        # neighbor list below it stay consistent; fall back to topology_info.
        degrees_str = None
        if adjacency:
            num_chains = len(adjacency)
            num_interfaces = len({
                frozenset((c, n)) for c, ns in adjacency.items() for n in ns if n != c
            })
            components = self._count_components(adjacency)
            degrees_str = self._format_degree_histogram(adjacency)
        else:
            num_chains = topology_info.get('num_chains', 'N/A')
            num_interfaces = topology_info.get('num_interfaces', 'N/A')
            components = 1 if topology_info.get('is_connected') else None

        if components is None:
            connected_str = 'unknown'
        else:
            plural = '' if components == 1 else 's'
            yn = 'yes' if components == 1 else 'no'
            connected_str = f"{yn} ({components} component{plural})"

        rows = [
            ("Chains", str(num_chains)),
            ("Interfaces", str(num_interfaces)),
            ("Connected", connected_str),
            ("Structure", structure_str),
        ]
        if degrees_str:
            rows.append(("Degrees", degrees_str))

        # Aligned label column: blue labels + default-foreground values,
        # built as Text so Rich's number highlighter never recolors them.
        body = Text()
        for i, (label, value) in enumerate(rows):
            if i:
                body.append("\n")
            body.append(f"{label + ':':<12}", style="blue")
            body.append(value)

        self.processor.console.print()
        self.processor.console.print(
            Panel(body, title="Chain Topology", border_style="blue", expand=False)
        )

        # One-line detection-criteria caption (grey, non-highlighted).
        self.processor.console.print(
            "Interfaces: buried area ≥200 Å² or contact ≤4.5 Å",
            style="grey50", highlight=False,
        )

        # Per-chain neighbor list with buried surface area, sorted by BSA
        # descending so the strongest interface leads each row. This is the
        # screen-reader-friendly complement to the visual heatmap: same
        # data, linear form. markup=False so chain letters/arrows/brackets
        # render literally; highlight=False to avoid auto-recoloring. The
        # header carries a literal "[BSA, Å²]" so it is a Text (not markup).
        if adjacency:
            self.processor.console.print()
            self.processor.console.print(
                Text("Connectivity (chain ↔ neighbor [BSA, Å²]):", style="bold blue")
            )
            lines = []
            for chain in sorted(adjacency):
                areas = chain_info[chain].get('interface_areas', {}) or {}
                # Ties (equal BSA) broken alphabetically by neighbor.
                ordered = sorted(
                    adjacency[chain],
                    key=lambda n: (-int(areas.get(n, 0)), n),
                )
                parts = [f"{n} ({int(areas.get(n, 0))})" for n in ordered]
                rhs = ", ".join(parts) if parts else "(none)"
                lines.append(f"  {chain} ↔ {rhs}")
            self.processor.console.print(
                "\n".join(lines), markup=False, highlight=False
            )

    def _count_components(self, adjacency: Dict[str, list]) -> int:
        """Number of connected components in the chain-interface graph."""
        seen = set()
        components = 0
        for start in adjacency:
            if start in seen:
                continue
            components += 1
            stack = [start]
            while stack:
                node = stack.pop()
                if node in seen:
                    continue
                seen.add(node)
                for neighbor in adjacency.get(node, []):
                    if neighbor not in seen:
                        stack.append(neighbor)
        return components

    def _format_degree_histogram(self, adjacency: Dict[str, list]) -> str:
        """Summarize the degree distribution, e.g. '6×5, 6×3 (range 3–5)'."""
        from collections import Counter
        degrees = [len(neighbors) for neighbors in adjacency.values()]
        if not degrees:
            return "0"
        counts = Counter(degrees)
        parts = [f"{counts[d]}×{d}" for d in sorted(counts, reverse=True)]
        lo, hi = min(degrees), max(degrees)
        span = f"(range {lo}–{hi})" if lo != hi else f"(uniform {lo})"
        return f"{', '.join(parts)} {span}"
             
    def _display_filter_summary(self, stats: Dict[str, Any]):
        """Display a summary of the filtering process."""
        self.processor.console.print("\n[bold]Original Structure Summary[/bold]")
        self.processor.console.print(f"Models: {stats['original']['models']}")
        self.processor.console.print(f"Chains: {len(stats['original']['chains'])}")

        total_residues = sum(
            chain["residue_count"] for chain in stats["original"]["chains"].values()
        )
        self.processor.console.print(f"Total residues: {total_residues}")

        self.processor.console.print("\n[bold]Filtered Structure Summary[/bold]")
        self.processor.console.print(f"Models: {stats['filtered']['models']}")
        self.processor.console.print(f"Chains: {len(stats['filtered']['chains'])}")

        total_residues = sum(
            chain["residue_count"] for chain in stats["filtered"]["chains"].values()
        )
        self.processor.console.print(f"Total residues: {total_residues}")

        self.processor.console.print("\n[bold]Removed Components[/bold]")
        if stats["removed"].get("chains"):
            self.processor.console.print(
                f"Chains: {', '.join(stats['removed']['chains'])}"
            )

        if stats["removed"].get("residues"):
            self.processor.console.print("Residues:")
            for chain_id, count in stats["removed"]["residues"].items():
                self.processor.console.print(
                    f"  Chain {chain_id}: {count} residues removed"
                )

        if stats["removed"].get("hetero_groups"):
            self.processor.console.print(
                f"Hetero groups: {', '.join(stats['removed']['hetero_groups'])}"
            )

    def _sync_redox_sites_with_filtered_structure(self, workspace: Dict[str, Any],
                                                    filtered_structure) -> int:
        """
        Sync RedoxSite objects with the filtered structure.

        After filtering, atoms/residues may have been removed from the structure.
        This updates RedoxSite objects to only reference atoms that still exist
        in the filtered structure, keeping the two in sync.

        Uses residue-level identity matching (chain, resid, insertion_code) since
        PDB filtering removes whole residues. This is more robust than coordinate
        matching which can fail due to floating-point representation differences
        between independent BioPython parse calls.

        Args:
            workspace: Current workspace containing detected_redox_sites
            filtered_structure: The filtered BioPython Structure object

        Returns:
            int: Total number of atoms removed from all RedoxSites
        """
        detected_redox_sites = workspace.get('detected_redox_sites', [])
        if not detected_redox_sites or filtered_structure is None:
            logger.debug("SYNC: No redox sites or filtered_structure is None, skipping")
            return 0

        # Build set of residues present in the filtered structure
        # Key: (chain_id, resid, insertion_code)
        # BioPython uses ' ' for empty insertion codes; normalize to ''
        filtered_residues = set()
        for model in filtered_structure:
            for chain in model:
                for residue in chain:
                    res_id = residue.get_id()
                    icode = res_id[2].strip()
                    filtered_residues.add((chain.id, res_id[1], icode))

        total_removed = 0

        for site in detected_redox_sites:
            # Identify atoms whose residue is no longer in the filtered structure
            retained_atoms = []
            coords_to_remove = []

            for atom in site.atoms:
                icode = (getattr(atom, 'insertion_code', '') or '').strip()
                res_key = (atom.chain, atom.resid, icode)
                if res_key in filtered_residues:
                    retained_atoms.append(atom)
                else:
                    coords_to_remove.append(atom.coords)
                    total_removed += 1

            site.atoms = retained_atoms

            # Remove stale coordinates from coord_to_pdb mapping
            for coord in coords_to_remove:
                if coord in site.coord_to_pdb:
                    del site.coord_to_pdb[coord]

            # Update residue_groups to remove stale coordinates
            for residue_key, coord_list in site.residue_groups.items():
                site.residue_groups[residue_key] = [
                    c for c in coord_list if c not in coords_to_remove
                ]

            # Remove bonds involving removed atoms
            site.bonds = [
                bond for bond in site.bonds
                if bond.atom1_coords not in coords_to_remove
                and bond.atom2_coords not in coords_to_remove
            ]

            # Remove centers that are no longer in the filtered structure
            retained_centers = []
            for center in site.centers:
                icode = (getattr(center, 'insertion_code', '') or '').strip()
                res_key = (center.chain, center.resid, icode)
                if res_key in filtered_residues:
                    retained_centers.append(center)
            site.centers = retained_centers

        return total_removed

    def _detect_redox_sites(self, pdb_file: str, structure):
        """
        Detect redox-active sites in the structure, with option to import from JSON.

        Args:
            pdb_file: Path to PDB file (already selected by user)
            structure: BioPython Structure object (already loaded)

        Returns:
            Tuple of (detected_sites, detector) or (None, None) if detection fails
        """
        if not REDOX_DETECTION_AVAILABLE:
            self.processor.console.print("[red]Redox detection module not available[/red]")
            return None, None

        try:

            # Check if user wants to import from existing JSON first
            json_file = self._prompt_for_json_import()
            if json_file:
                detected_sites = self._import_redox_sites_from_json(json_file)
                return detected_sites, None  # No detector instance for imported sites

            # Otherwise proceed with normal detection using provided structure
            if not pdb_file or not structure:
                self.processor.console.print("[red]No structure available for redox detection[/red]")
                return None, None

            # Initialize detector with console from processor
            detector = ComprehensiveRedoxDetector(console=self.processor.console, processor=self.processor)

            # Set source PDB filename for export functionality
            detector.source_pdb_file = pdb_file
            
            # Run detection on the original structure
            self.processor.console.print("[cyan]Running comprehensive redox site detection...[/cyan]")
            detected_sites = detector.detect_redox_sites(
                structure=structure,  # Pass structure object directly
                selected_chains=None,  # Analyze all chains initially
                interactive=True  # Allow configuration prompts
            )
            
            if detected_sites:
                self.processor.console.print(f"[green]✓ Found {len(detected_sites)} redox-active site(s)[/green]")
            else:
                self.processor.console.print("[yellow]No redox-active sites detected[/yellow]")
                
            return detected_sites, detector
            
        except Exception as e:
            self.processor.console.print(f"[red]Error during redox detection: {str(e)}[/red]")
            logger.error(f"Redox detection failed: {e}", exc_info=True)
            return None, None

    def _prompt_for_json_import(self):
        """
        Prompt user to import redox sites from JSON file if available.
        
        Returns:
            Selected JSON filename or None if user chooses not to import
        """
        import os
        import glob
        
        # Ask if user wants to import
        if not confirm_with_context(
            processor=self.processor,
            prompt="[bold cyan]Import redox sites from previous JSON file?[/bold cyan]",
            default=False,
            module="PDB Filter",
            description="Import redox sites from JSON file"
        ):
            return None
            
        # Find available JSON files - use specific patterns to avoid ProPrep session files
        json_files = []
        patterns = [
            "*_redox_sites.json",  # Primary pattern used by comprehensive redox detector
            "*redox*.json",        # Alternative redox-related files
        ]
        
        for pattern in patterns:
            found_files = glob.glob(pattern)
            json_files.extend(found_files)
        
        # Remove duplicates and sort
        json_files = sorted(list(set(json_files)))
        
        if not json_files:
            self.processor.console.print("[yellow]No JSON files found in current directory.[/yellow]")
            return None
            
        # Show available files in a table
        from rich.table import Table
        
        table = Table(title="Available JSON Files")
        table.add_column("Index", style="cyan")
        table.add_column("Filename", style="green") 
        table.add_column("Size", style="yellow")
        
        for i, json_file in enumerate(json_files, 1):
            try:
                file_size = os.path.getsize(json_file)
                if file_size > 1024:
                    size_str = f"{file_size / 1024:.1f} KB"
                else:
                    size_str = f"{file_size} bytes"
            except OSError:
                size_str = "Unknown"
                
            table.add_row(str(i), json_file, size_str)
            
        self.processor.console.print(table)
        
        # Get user selection
        while True:
            # Build options map
            file_options = {str(i): json_files[i-1] for i in range(1, len(json_files) + 1)}
            file_options["cancel"] = "Cancel import"

            choice = prompt_with_context(
                processor=self.processor,
                prompt="Select JSON file to import",
                choices=[str(i) for i in range(1, len(json_files) + 1)] + ["cancel"],
                default="cancel",
                module="PDB Filter",
                description="Select redox sites JSON file to import",
                options_map=file_options
            )
            
            if choice == "cancel":
                return None
                
            try:
                file_idx = int(choice) - 1
                selected_file = json_files[file_idx]
                self.processor.console.print(f"[green]Selected: {selected_file}[/green]")
                return selected_file
            except (ValueError, IndexError):
                self.processor.console.print(f"[red]Invalid selection: {choice}[/red]")
                
    def _import_redox_sites_from_json(self, json_file):
        """
        Import redox sites from JSON file using the comprehensive detector's import functionality.
        
        Args:
            json_file: Path to JSON file containing redox site data
            
        Returns:
            List of RedoxSite objects or None if import fails
        """
        try:
            # Import the _import_from_json function from comprehensive_redox_detector  
            from .comprehensive_redox_detector import _import_from_json
            
            self.processor.console.print(f"[cyan]Importing redox sites from: {json_file}[/cyan]")
            detected_sites, _ = _import_from_json(json_file)
            
            if detected_sites:
                self.processor.console.print(f"[green]✓ Successfully imported {len(detected_sites)} redox-active site(s)[/green]")
                
                # Show imported sites summary and provide interactive options
                self._show_imported_sites_summary(detected_sites, json_file)
                
            else:
                self.processor.console.print("[yellow]No redox sites found in JSON file[/yellow]")
                
            return detected_sites
            
        except FileNotFoundError:
            self.processor.console.print(f"[red]Error: JSON file not found: {json_file}[/red]")
            return None
        except Exception as e:
            self.processor.console.print(f"[red]Error importing from JSON: {str(e)}[/red]")
            logger.error(f"JSON import failed: {e}", exc_info=True)
            return None

    def _show_imported_sites_summary(self, sites, json_file):
        """
        Show summary of imported redox sites with interactive options.
        
        Args:
            sites: List of imported RedoxSite objects
            json_file: Source JSON filename for context
        """
        from rich.table import Table
        
        # Show overview table of imported sites
        sites_table = Table(title=f"Imported Redox Sites from {json_file}")
        sites_table.add_column("Index", style="cyan")
        sites_table.add_column("Site ID", style="bold green")
        sites_table.add_column("Atoms", style="yellow")
        sites_table.add_column("Bonds", style="magenta")
        sites_table.add_column("Centers", style="red")
        
        for i, site in enumerate(sites, 1):
            sites_table.add_row(
                f"[{i}]",
                site.site_id,
                str(len(site.atoms)),
                str(len(site.bonds)),
                str(len(site.centers))
            )
        
        self.processor.console.print(sites_table)
        
        # Interactive options
        while True:
            self.processor.console.print("\n[bold]Options:[/bold]")
            self.processor.console.print("[grey50]1.[/grey50] View detailed site summaries")
            self.processor.console.print("[grey50]2.[/grey50] Export sites to files")
            self.processor.console.print("[grey50]3.[/grey50] Continue with filtering")
            
            choice = prompt_with_context(
                processor=self.processor,
                prompt="[green]Choose option[/green]",
                choices=["1", "2", "3"],
                default="3",
                module="PDB Filter",
                description="Imported redox sites - select action",
                options_map={
                    "1": "View detailed site summaries",
                    "2": "Export sites to files",
                    "3": "Continue with filtering"
                }
            )
            
            if choice == "1":
                self._show_detailed_site_summaries(sites)
            elif choice == "2":
                self._export_imported_sites(sites, json_file)
            else:  # choice == "3"
                break
                
    def _show_detailed_site_summaries(self, sites):
        """Show detailed summaries for selected sites."""
        selection = prompt_with_context(
            processor=self.processor,
            prompt=f"[green]Select sites to view[/green] (1-{len(sites)}, comma-separated, or 'all')",
            default="all",
            module="PDB Filter",
            description="Select redox sites to view detailed summaries"
        ).strip()
        
        if selection.lower() == "all":
            indices = list(range(len(sites)))
        else:
            try:
                # Parse comma-separated selection
                indices = []
                for part in selection.split(","):
                    part = part.strip()
                    if "-" in part:  # Handle ranges like "1-3"
                        start, end = map(int, part.split("-"))
                        indices.extend(range(start-1, end))
                    else:
                        indices.append(int(part) - 1)
                        
                # Validate indices
                indices = [i for i in indices if 0 <= i < len(sites)]
                
            except ValueError:
                self.processor.console.print("[red]Invalid selection format[/red]")
                return
                
        # Show detailed summaries
        try:
            from .comprehensive_redox_detector import SiteRefinementInterface, ComprehensiveRedoxDetector
            
            # Create a temporary detector instance to access configuration
            detector = ComprehensiveRedoxDetector(console=self.processor.console, processor=self.processor)
            refinement_interface = SiteRefinementInterface(detector.config, console=self.processor.console)
            
            for i in indices:
                site = sites[i]
                self.processor.console.print(f"\n[bold underline]Site Summary: {site.site_id}[/bold underline]")
                refinement_interface._display_site_summary(site)
                
            # Simple continue prompt
            input("\nPress Enter to continue...")
            
        except Exception as e:
            self.processor.console.print(f"[red]Error displaying site summaries: {str(e)}[/red]")
            
    def _export_imported_sites(self, sites, original_json_file):
        """Export imported sites to various formats."""
        from rich.prompt import Confirm

        # Ask if user wants to export
        if not confirm_with_context(
            processor=self.processor,
            prompt="[green]Export redox sites?[/green]",
            default=False,
            module="PDB Filter",
            description="Export redox sites to files"
        ):
            return
            
        # Show export format options
        self.processor.console.print("\n[bold]Export formats:[/bold]")
        self.processor.console.print("[grey50]1.[/grey50] JSON (complete site data)")
        self.processor.console.print("[grey50]2.[/grey50] PDB - all sites in one file")
        self.processor.console.print("[grey50]3.[/grey50] PDB - each site in separate files")
        
        format_choice = prompt_with_context(
            processor=self.processor,
            prompt="[green]Choose export format[/green] (comma-separated list)",
            default="1",
            module="PDB Filter",
            description="Select redox site export format(s)"
        ).strip()
        
        # Parse comma-separated choices
        choices = [choice.strip() for choice in format_choice.split(',') if choice.strip()]
        
        # Validate choices
        valid_choices = {'1', '2', '3'}
        invalid_choices = [c for c in choices if c not in valid_choices]
        if invalid_choices:
            self.processor.console.print(f"[red]Invalid choices: {', '.join(invalid_choices)}. Valid options are 1, 2, 3[/red]")
            return
            
        if not choices:
            self.processor.console.print("[yellow]No export formats selected[/yellow]")
            return
        
        # Get source PDB filename for export
        workspace = self.processor._get_workspace()
        source_pdb = self._get_pdb_file_from_workspace(workspace) or "structure.pdb"
        
        try:
            from .comprehensive_redox_detector import _export_to_json, _export_to_pdb_single, _export_to_pdb_separate
            
            if "1" in choices:
                _export_to_json(sites, source_pdb, self.processor.console)
                self.processor.console.print("[green]✓ JSON export complete[/green]")
                
            if "2" in choices:
                _export_to_pdb_single(sites, source_pdb, self.processor.console)
                self.processor.console.print("[green]✓ Single PDB export complete[/green]")
                
            if "3" in choices:
                _export_to_pdb_separate(sites, source_pdb, self.processor.console)
                self.processor.console.print("[green]✓ Separate PDB export complete[/green]")
                
        except Exception as e:
            self.processor.console.print(f"[red]Export failed: {str(e)}[/red]")


    def _check_chain_redox_involvement(self, chain_id, redox_sites):
        """
        Check if a chain is involved in any redox sites.
        
        Args:
            chain_id: Chain identifier to check
            redox_sites: List of RedoxSite objects from comprehensive detector
            
        Returns:
            Dict with redox site involvement information or None
        """
        if not redox_sites:
            return None
            
        involvement = {
            'centers': [],  # Sites where this chain has the redox center
            'coordinators': [],  # Sites where this chain provides other atoms (ligands)
            'sites': []  # All sites this chain is involved in
        }
        
        for site in redox_sites:
            site_involved = False
            
            # Check if chain contains redox centers
            if hasattr(site, 'centers') and site.centers:
                for center in site.centers:
                    if center.chain == chain_id:
                        involvement['centers'].append(site)
                        site_involved = True
                        break
                        
            # Check if chain provides other atoms (ligands, etc.)
            if hasattr(site, 'atoms') and site.atoms:
                for atom in site.atoms:
                    if atom.chain == chain_id:
                        # Only add to coordinators if not already a center chain
                        if site not in involvement['centers']:
                            involvement['coordinators'].append(site)
                        site_involved = True
                        break
                        
            if site_involved:
                involvement['sites'].append(site)
                
        return involvement if involvement['sites'] else None

    def _check_redox_conflict_for_component(self, chain_id, comp_type, chain_redox_involvement):
        """
        Check if filtering a component type would affect redox sites.
        
        Args:
            chain_id: Chain identifier
            comp_type: Component type being filtered
            chain_redox_involvement: Result from _check_chain_redox_involvement
            
        Returns:
            List of redox conflict information or None
        """
        if not chain_redox_involvement:
            return None
            
        conflicts = []
        
        for site in chain_redox_involvement['sites']:
            center_residues = set()
            other_residues = set()
            
            # Collect unique center residues from this chain that match the component type
            if hasattr(site, 'centers') and site.centers:
                for center in site.centers:
                    if center.chain == chain_id:
                        center_comp_type = self._get_residue_component_type(center.resname)
                        if center_comp_type == comp_type:
                            center_residues.add(f"{center.resname}{center.resid}")
                        
            # Collect unique other residues from this chain that match the component type
            if hasattr(site, 'atoms') and site.atoms:
                for atom in site.atoms:
                    if atom.chain == chain_id:
                        atom_comp_type = self._get_residue_component_type(atom.resname)
                        if atom_comp_type == comp_type:
                            # Don't double-count centers
                            residue_id = f"{atom.resname}{atom.resid}"
                            if residue_id not in center_residues:
                                other_residues.add(residue_id)
                        
            if center_residues or other_residues:
                site_parts = []
                if center_residues:
                    site_parts.append(f"Centers: {', '.join(sorted(center_residues))}")
                if other_residues:
                    site_parts.append(f"Residues: {', '.join(sorted(other_residues))}")
                    
                site_type = getattr(site, 'site_type', 'unknown')
                conflicts.append({
                    'site_id': f"{site.site_id} ({site_type})",
                    'residues': '; '.join(site_parts)
                })
                
        return conflicts if conflicts else None

    def _get_redox_site_residues_for_component(self, chain_id, comp_type, chain_redox_involvement):
        """
        Get residue IDs of redox site residues for a specific component type.
        
        Args:
            chain_id: Chain identifier
            comp_type: Component type
            chain_redox_involvement: Result from _check_chain_redox_involvement
            
        Returns:
            Set of residue IDs
        """
        if not chain_redox_involvement:
            return set()
            
        redox_residues = set()
        
        for site in chain_redox_involvement['sites']:
            # Add center residues that match the component type
            if hasattr(site, 'centers') and site.centers:
                for center in site.centers:
                    if center.chain == chain_id:
                        center_comp_type = self._get_residue_component_type(center.resname)
                        if center_comp_type == comp_type:
                            redox_residues.add(center.resid)
                        
            # Add other atoms from the same chain that match the component type
            if hasattr(site, 'atoms') and site.atoms:
                for atom in site.atoms:
                    if atom.chain == chain_id:
                        atom_comp_type = self._get_residue_component_type(atom.resname)
                        if atom_comp_type == comp_type:
                            redox_residues.add(atom.resid)
                        
        return redox_residues

    def _get_flanking_residues(self, chain, redox_residues, flanking_count):
        """
        Get redox site residues plus flanking residues.
        
        Args:
            chain: BioPython Chain object
            redox_residues: Set of redox site residue IDs
            flanking_count: Number of flanking residues on each side
            
        Returns:
            Set of residue IDs including flanking
        """
        if not redox_residues:
            return set()
            
        # Get all residue IDs in chain order
        all_residues = [residue.id[1] for residue in chain]
        all_residues.sort()
        
        flanking_residues = set(redox_residues)
        
        for coord_resid in redox_residues:
            try:
                coord_idx = all_residues.index(coord_resid)
                
                # Add flanking residues before
                start_idx = max(0, coord_idx - flanking_count)
                for i in range(start_idx, coord_idx):
                    flanking_residues.add(all_residues[i])
                    
                # Add flanking residues after
                end_idx = min(len(all_residues), coord_idx + flanking_count + 1)
                for i in range(coord_idx + 1, end_idx):
                    flanking_residues.add(all_residues[i])
                    
            except ValueError:
                # Residue not found in chain, skip
                continue
                
        return flanking_residues

    def _launch_structure_viewer(self, force: bool = True):
        """
        Launch interactive 3D structure viewer to visualize the structure
        before filtering. Routed through the ``ViewerCoordinator`` singleton
        so the launch shares state with all other coordinator-driven hooks
        across the session (previously instantiated its own
        ``InteractiveStructureViewer``, creating a phantom viewer instance
        that didn't share annotation state with the rest of ProPrep).

        ``force`` defaults to True for the user-initiated terminal path (the
        "View structure?" Y/N gate). The web-shell auto-show path passes
        force=False so the structure is shown passively in the already-docked
        viewer without popping a fresh tab.
        """
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer
        except ImportError:
            self.processor.console.print("[yellow]Structure viewer not available[/yellow]")
            return

        # Get structure file from workspace
        workspace = self.processor._get_workspace()
        structure_file = None

        # Try different workspace keys for structure file
        for key in ['processed_pdb_file', 'local_pdb_file', 'rcsb_pdb_file', 'alphafold_pdb_file', 'alphafill_pdb_file']:
            structure_file = workspace.get(key)
            if structure_file:
                break

        if not structure_file:
            self.processor.console.print("[yellow]No structure file found in workspace[/yellow]")
            return

        # Wire processor for session recording, then route the launch
        # through the coordinator. force=True (terminal path) is a
        # user-initiated view via the "View structure?" Y/N gate above;
        # force=False (web-shell path) shows passively in the docked viewer.
        if force:
            self.processor.console.print(
                f"[cyan]Launching interactive structure viewer for {structure_file}...[/cyan]"
            )
        else:
            self.processor.console.print(
                f"[grey50]Showing {structure_file} in the docked viewer...[/grey50]"
            )
        _viewer.set_processor(self.processor)
        _viewer.show_structure(structure_file, force=force)

    def _get_chain_redox_display(self, chain_id, detected_redox_sites):
        """
        Get display string for chain's redox site involvement.

        Args:
            chain_id: Chain identifier to check
            detected_redox_sites: List of RedoxSite objects

        Returns:
            String with redox involvement information for display
        """
        if not detected_redox_sites:
            return ""

        involvement = self._check_chain_redox_involvement(chain_id, detected_redox_sites)

        if not involvement:
            return ""

        site_count = len(involvement['sites'])
        center_count = len(involvement['centers'])

        display_parts = []

        if center_count > 0:
            display_parts.append(f"[bold red]{center_count} redox center(s)[/bold red]")

        if len(involvement['coordinators']) > 0:
            coord_count = len(involvement['coordinators'])
            if coord_count != center_count:  # Don't double-count if center and coordinator are same
                display_parts.append(f"[bold yellow]{coord_count} site ligand(s)[/bold yellow]")

        if display_parts:
            return f" - [✓ {' + '.join(display_parts)}]"

        return ""

    def _get_residue_component_type(self, resname):
        """
        Get the component type for a residue name.
        
        Args:
            resname: Residue name (e.g., 'HIS', 'HEC', 'ALA')
            
        Returns:
            Component type string
        """
        # Use the same logic as ComponentClassifier but simplified for residue names
        if not self.filter_worker or not self.filter_worker.ccd_parser:
            # Fallback classification
            standard_aa = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
                          'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'}
            if resname in standard_aa:
                return "amino_acid"
            elif resname in ('HOH', 'WAT'):
                return "water"
            else:
                return resname  # Return the residue name as the component type for non-standard
        
        # Use the CCD parser if available
        try:
            if self.filter_worker.ccd_parser.is_standard_amino_acid(resname):
                return "amino_acid"
            elif self.filter_worker.ccd_parser.is_dna_base(resname):
                return "dna_base"
            elif self.filter_worker.ccd_parser.is_rna_base(resname):
                return "rna_base"
            elif resname in ('HOH', 'WAT'):
                return "water"
            else:
                return resname  # Return the residue name as the component type for non-standard
        except:
            # Fallback if CCD parser fails
            return resname

    def _check_chain_exclusion_conflicts(self, excluded_chains, detected_redox_sites):
        """
        Check if excluding chains would break redox sites.
        
        Args:
            excluded_chains: List of chain IDs being excluded
            detected_redox_sites: List of RedoxSite objects
            
        Returns:
            List of conflict descriptions
        """
        if not detected_redox_sites or not excluded_chains:
            return []
            
        conflicts = []
        
        for site in detected_redox_sites:
            affected_chains = []
            
            # Check if any centers are in excluded chains
            if hasattr(site, 'centers') and site.centers:
                for center in site.centers:
                    if center.chain in excluded_chains:
                        affected_chains.append(f"center in chain {center.chain}")
                        
            # Check if any other atoms are in excluded chains
            if hasattr(site, 'atoms') and site.atoms:
                for atom in site.atoms:
                    if atom.chain in excluded_chains:
                        affected_chains.append(f"residues in chain {atom.chain}")
                        break  # Don't repeat for same chain
                        
            if affected_chains:
                site_type = getattr(site, 'site_type', 'unknown')
                conflicts.append(f"Site {site.site_id} ({site_type}) has {', '.join(set(affected_chains))}")
                
        return conflicts

    def _get_required_chains_for_redox_sites(self, selected_chains, detected_redox_sites):
        """
        Get additional chains required to maintain redox site integrity.
        
        Args:
            selected_chains: List of currently selected chain IDs
            detected_redox_sites: List of RedoxSite objects
            
        Returns:
            List of additional chain IDs needed
        """
        if not detected_redox_sites:
            return []
            
        required_chains = set()
        
        for site in detected_redox_sites:
            site_chains = set()
            
            # Get all chains involved in this site
            if hasattr(site, 'centers') and site.centers:
                for center in site.centers:
                    site_chains.add(center.chain)
                    
            if hasattr(site, 'atoms') and site.atoms:
                for atom in site.atoms:
                    site_chains.add(atom.chain)
                    
            # If any chain from this site is selected, all chains must be selected
            if site_chains.intersection(selected_chains):
                required_chains.update(site_chains - set(selected_chains))
                
        return list(required_chains)

    def _check_redox_conflict_for_component_integrated(self, comp_type, chain_composition):
        """
        Check if a component type has redox site involvement using integrated composition data.
        
        Args:
            comp_type: Component type being evaluated
            chain_composition: Pre-computed composition with redox information
            
        Returns:
            List of conflict information or None
        """
        if comp_type not in chain_composition:
            return None
            
        component_data = chain_composition[comp_type]
        
        # Build site-specific residue mapping
        site_residues = {}  # site_id -> set[residue_ids] (using set to avoid duplicates)
        
        for residue_name, residue_data in component_data.items():
            if isinstance(residue_data, dict) and residue_data.get("redox_residues"):
                redox_residues = residue_data["redox_residues"]
                
                for resid, site_list in redox_residues.items():
                    residue_id = f"{residue_name}{resid}"
                    for site in site_list:
                        if site not in site_residues:
                            site_residues[site] = set()
                        site_residues[site].add(residue_id)  # Use set.add() to automatically deduplicate
        
        # Build conflicts with site-specific residue lists
        conflicts = []
        for site, residues in site_residues.items():
            conflicts.append({
                'site_id': site,
                'residues': ', '.join(sorted(residues))  # Convert set to sorted list for display
            })
        
        return conflicts if conflicts else None

    def _get_redox_site_residues_for_component_integrated(self, comp_type, chain_composition):
        """
        Get redox site residue IDs for a component type using integrated composition data.
        
        Args:
            comp_type: Component type
            chain_composition: Pre-computed composition with redox information
            
        Returns:
            Set of residue IDs
        """
        if comp_type not in chain_composition:
            return set()
            
        redox_residues = set()
        component_data = chain_composition[comp_type]
        
        for residue_name, residue_data in component_data.items():
            if isinstance(residue_data, dict) and residue_data.get("redox_residues"):
                redox_residues.update(residue_data["redox_residues"].keys())
        
        return redox_residues

    def _record_local_command(self, command):
        """Record command in local history for module-specific tracking"""
        self._local_command_history.append(command)

    def get_breadcrumb_for_command(self, command):
        """Return breadcrumb string for a command"""
        action_to_run = getattr(command, "_action_to_run", None)
        if action_to_run:
            action_map = {
                "filter_pdb_structure": "Filter",
                "show_filter_status": "Status",
                "export_filter_statistics": "Export",
            }
            action_name = action_map.get(action_to_run, action_to_run)
            return f"{self.NAME} > {action_name}"
        return self.NAME

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

    def _get_structure_from_workspace(self, workspace: Dict[str, Any], pdb_file: str = None, silent: bool = False):
        """Get structure object from workspace using structure selector.

        Args:
            workspace: Current workspace
            pdb_file: Optional PDB file path to match against workspace keys
            silent: If True, suppress console output

        Returns:
            BioPython Structure object, or None if no structure available
        """
        from proprep.utils.structure_selector import StructureSelector, StructureRegistry

        selector = StructureSelector(workspace, self.console, processor=self.processor)

        # If pdb_file provided, try to find matching structure object
        if pdb_file:
            # Map file path to structure object key
            # e.g., rcsb_pdb_file -> rcsb_structure
            for structure_type in StructureRegistry.get_all():
                file_key = structure_type.workspace_key
                file_path = workspace.get(file_key)

                if file_path and str(file_path) == str(pdb_file):
                    # Found matching file, get corresponding structure object key.
                    # Handle both the '<x>_pdb_file'/'<x>_file' convention and the
                    # bare '<x>_pdb' keys (e.g. 'prepared_pdb'), which neither
                    # replace below would otherwise rewrite.
                    structure_key = (file_key
                                     .replace('_pdb_file', '_structure')
                                     .replace('_file', '_structure'))
                    if structure_key == file_key and file_key.endswith('_pdb'):
                        # no suffix matched (e.g. 'prepared_pdb' -> 'prepared_structure')
                        structure_key = file_key[:-len('_pdb')] + '_structure'
                    structure = selector.get_structure_by_key(structure_key, require_exists=False)

                    # get_structure_by_key may return a file-path STRING when the
                    # key holds only a path (no cached Structure object). A string
                    # is NOT a usable structure — treat it as "no cached object"
                    # and fall through to return None so the worker parses the
                    # selected file. (Guards against len(path)-'models' and the
                    # 'str has no attribute id' crash downstream.)
                    if structure is not None and not isinstance(structure, str):
                        if not silent:
                            self.console.print(
                                f"[green]Using {structure_type.display_name} Structure object "
                                f"(workspace key: {structure_key})[/green]"
                            )
                        return structure

                    # The selected file matched a registered key but has no
                    # cached Structure object (e.g. the Biological Assembly
                    # Generator stores only a path, never a parsed structure).
                    # Return None so the worker PARSES THIS FILE, rather than
                    # falling through to the priority fallback below, which
                    # would substitute a cached structure for a DIFFERENT
                    # file — e.g. analyzing the 4-chain asymmetric unit in
                    # place of the selected 12-chain biological assembly.
                    if not silent:
                        self.console.print(
                            f"[grey50]No cached structure for {structure_type.display_name}; "
                            f"parsing the selected file directly.[/grey50]"
                        )
                    return None

        # Fallback: use priority-based structure object selection
        structure = selector.get_structure_object(silent=silent)

        return structure

    def get_workspace_requirements(self) -> List[str]:
        """Get workspace requirements - needs at least one structure loaded"""
        return [
            "rcsb_pdb_file | local_pdb_file | alphafold_pdb_file | alphafold_homolog_pdb_file"
        ]

    def get_workspace_outputs(self) -> List[str]:
        """Get workspace outputs"""
        return [
            "filtered_structure",
            "filtered_pdb_file",
            "filter_selections",
            "is_hplusplus_structure",
        ]

    def can_process(self, workspace: Dict[str, Any]) -> bool:
        """Check if the module can process the current workspace"""
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(workspace, self.console)
        status = selector.get_structure_status()
        return status.get("has_any", False)

    def process(self, workspace):
        """Process the workspace"""
        if self.can_process(workspace):
            pdb_file = self._get_pdb_file_from_workspace(workspace)
            existing_structure = self._get_structure_from_workspace(workspace, pdb_file)

            if pdb_file:
                self.filter_worker = PDBFilterWorker(pdb_file, existing_structure, processor=self.processor)

                auto_filter = workspace.get("auto_filter", False)
                if auto_filter:
                    workspace = self.filter_pdb_structure(workspace, interactive=False)

        return workspace

    def cleanup(self):
        """Clean up module resources"""
        if self.filter_worker:
            self.filter_worker.clear_cache()
        self.filter_worker = None
        self.filtered_structure = None
        self.filter_selections = {}
        self._local_command_history.clear()
