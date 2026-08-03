"""
PDB Processor for MD Simulations

A modular workflow manager for processing PDB structures for MD simulations.
Integrates various structural analysis and modification modules through a menu interface.
"""

import importlib
import json
import os
from pathlib import Path
from typing import Optional

# Rich UI
from rich.console import Console
from rich.panel import Panel
from proprep.utils.prompts import prompt_with_context, confirm_with_context
from rich.status import Status
from rich.table import Table

# Bio Python
from Bio.PDB import PDBParser

# custom modules
from proprep.structure_prep.pdb_filter_worker import PDBFilterTool
from proprep.structure_prep.pdb_loader import (
    PDBMetadataExtractor,
    display_metadata_menu,
    download_pdb,
    load_pdb_interactive,
)

# Note: AmberInputGenerator is auto-registered via @register_module decorator
# from proprep.md_prep.amber_input_generator import AmberInputGenerator

from proprep.utils import MPSAWorkspace, Workspace
from proprep.utils.debug_utils import debug_workspace
from proprep.utils.module_registry import ProcessingModule, register_module, registry


# Import command pattern classes
from .command_history import CommandHistory
from .menu_commands import MainMenuCommand

from .processor import Processor


class PDBProcessor(Processor):
    """
    Main PDB processor class that coordinates the workflow modules.
    """

    def __init__(self):
        """Initialize the PDB processor"""
        self.console = Console()
        self.registry = registry  # Use the external registry
        
        # Initialize command history for breadcrumb tracking
        self.command_history = CommandHistory()
        
        # Initialize workspace - MODIFY THIS SECTION:
        # Original line (COMMENT OUT):
        self.workspace = MPSAWorkspace(debug=False)
        
        # Debug version (ADD THESE LINES):
#       from proprep.utils.mpsa_workspace import MPSAWorkspace
#       from proprep.utils.debug_workspace import DebugWorkspace
        
#       self._original_workspace = MPSAWorkspace(debug=False)
#       debug_log_file = "workspace_debug.log"
#       self.workspace = DebugWorkspace(self._original_workspace, log_file=debug_log_file)
        
#       print(f"\n{'='*60}")
#       print("DEBUG MODE ENABLED - Workspace access is being logged")
#       print(f"Log file: {debug_log_file}")
#       print(f"Monitoring: filtered_structure, repaired_structure, structure")
#       print(f"{'='*60}\n")
        
        # Rest remains the same:
        self.module_instances = {}
        
        # Initialize common workspace variables
        self.pdb_file = None
        self.pdb_id = None
        self.structure = None
        self.metadata = None
        self.filtered_structure = None
        self.filter_selections = None
        
        # Import required modules
        self._import_required_modules()

        # Check for module registration issues (only show errors)
        registered_modules = self.registry.get_module_names()
        if not registered_modules:
            self.console.print("[red]Error: No modules registered[/red]")

        # Discover external plugins via the ``proprep.plugins``
        # entry-point group. Empty when nothing is installed; never
        # raises. Discovered plugins are spliced into the workflow
        # menu under whatever stage they declare.
        self._plugins = {}
        self._discover_and_register_plugins()

        # Sync attributes to workspace
        self._sync_attributes_to_workspace()
        
    def _import_required_modules(self):
        """Import required modules from the package"""
        try:
            # Define modules to import with their new package paths
            modules_to_import = [
                "proprep.structure_prep.pdb_loader",
                "proprep.structure_prep.structure_loader",  # New unified structure loader
                # "proprep.structure_prep.alphafold_predictor",  # Unregistered - replaced by structure_loader
                "proprep.structure_prep.redox_detector_module",
                "proprep.emboss.emboss_module",
                "proprep.structure_prep.homology_searcher",
                "proprep.structure_prep.pdb_filter",
                "proprep.structure_prep.structure_alignment",
                "proprep.structure_prep.vmd_visualizer",
                "proprep.structure_prep.amino_acid_mutator",
                "proprep.structure_prep.structure_completeness",
                "proprep.redoxsite_prep.redoxsite_integration",
                "proprep.structure_prep.md_restraint_manager",
                "proprep.structure_prep.protonation_state_analyzer",
                "proprep.membrane_prep.membrane_builder",
                "proprep.tleap_prep.tleap_input_generator",
                "proprep.md_prep.molecular_dynamics_manager",
                "proprep.qmmm_prep.qmmm_preparator",
                "proprep.forcefield_prep.forcefield_parameterizer",
                "proprep.forcefield_prep.small_molecule_parameterizer",
                # Analysis modules (DISABLED for pre-analysis release - see docs/developer/ANALYSIS_MODULES_REENABLE.md)
                # "proprep.modules.analysis.electron_transfer",
            ]

            for module_name in modules_to_import:
                try:
                    # Import the module to trigger registration
                    importlib.import_module(module_name)
                except ImportError as e:
                    self.console.print(
                        f"[yellow]Warning: Could not import {module_name}: {str(e)}[/yellow]"
                    )

            self.console.print("[green]Successfully imported core modules[/green]")
        except Exception as e:
            self.console.print(f"[red]Error importing modules: {str(e)}[/red]")

    def _discover_and_register_plugins(self):
        """Find installed ProPrep plugins and splice them into the menu.

        Plugins register under the ``proprep.plugins`` entry-point
        group. Each one declares a stage (which may extend
        ``WorkflowStateManager.STAGES``) and a tool that lands under
        that stage in ``WorkflowMenuCommand.STAGE_TOOLS``. Discovery
        is best-effort — any plugin that fails to import, mismatches
        the API version, or returns ``False`` from ``is_available()``
        is silently dropped (with a warning logged via
        ``proprep.plugins``).
        """
        try:
            from proprep.plugins import discover_plugins
            from proprep.utils.workflow_state_manager import WorkflowStateManager
            from proprep.application.menu_commands import WorkflowMenuCommand
        except ImportError as e:
            # The plugins package itself failed to import — surface
            # the error but don't crash startup.
            self.console.print(
                f"[yellow]Plugin infrastructure unavailable: {e}[/yellow]"
            )
            return

        self._plugins = discover_plugins()
        for name, plugin in self._plugins.items():
            try:
                meta = plugin.get_metadata()
                WorkflowStateManager.register_plugin_stage(meta)
                WorkflowMenuCommand.register_plugin_tool(meta, plugin)
            except Exception as e:  # noqa: BLE001 — plugin isolation
                self.console.print(
                    f"[yellow]Plugin {name!r} failed during registration: "
                    f"{e}[/yellow]"
                )

        if self._plugins:
            self.console.print(
                f"[green]Loaded {len(self._plugins)} plugin(s): "
                f"{', '.join(sorted(self._plugins))}[/green]"
            )

    def build_plugin_seed(self, meta) -> dict:
        """Pull the workspace keys a plugin declared into a seed dict.

        ``meta.consumes_workspace_keys`` is the plugin's pull list:
        for each name there, if the host workspace currently holds
        the key, copy the value verbatim into the returned dict.
        Absent keys are omitted (not present as ``None``) so plugins
        can disambiguate "host doesn't have it" from "host has it
        and it's null".

        Generic by design — a new shared concept between ProPrep
        and a plugin is added by the plugin extending its
        ``consumes_workspace_keys`` list; ProPrep's seed builder
        needs no edit.
        """
        seed: dict = {}
        for key in getattr(meta, "consumes_workspace_keys", []):
            if self.workspace.has(key):
                seed[key] = self.workspace.get(key)
        return seed

    def setup_module_dependencies(self):
        """Establish dependencies between modules"""
        completeness_module = self.get_module_instance("Structure Completeness")
        repairer_module = self.get_module_instance("Structure Repairer")
        selector_module = self.get_module_instance("AltLoc Selector")
        metal_detector_module = self.get_module_instance("Metal Site Detector")

        # Define parent-child relationships (for menu navigation)
        if completeness_module and repairer_module:
            completeness_module.repairer_module = repairer_module

        if completeness_module and selector_module:
            completeness_module.selector_module = selector_module

        # Connect Metal Detector to other modules
        if metal_detector_module:
            if repairer_module:
                metal_detector_module.repairer_module = repairer_module

    def _sync_workspace_to_attributes(self):
        """Synchronize workspace values to processor attributes"""
        # List of attributes that should be synchronized
        sync_attrs = [
            "pdb_file",
            "pdb_id",
            "structure",
            "metadata",
            "filtered_structure",
            "filter_selections",
        ]

        # Sync from Workspace object
        for attr in sync_attrs:
            if self.workspace.has(attr):
                if attr == "pdb_file":
                    # Prefer original_pdb_file if available, fallback to pdb_file
                    self.pdb_file = self.workspace.get("original_pdb_file") or self.workspace.get(attr)
                elif attr == "pdb_id":
                    self.pdb_id = self.workspace.get(attr)
                elif attr == "structure":
                    self.structure = self.workspace.get(attr)
                elif attr == "metadata":
                    self.metadata = self.workspace.get(attr)
                elif attr == "filtered_structure":
                    self.filtered_structure = self.workspace.get(attr)
                elif attr == "filter_selections":
                    self.filter_selections = self.workspace.get(attr)

    def _sync_attributes_to_workspace(self):
        """Synchronize processor attributes to workspace"""
        # List of attributes that should be synchronized
        sync_attrs = [
            "pdb_file",
            "pdb_id",
            "structure",
            "metadata",
            "filtered_structure",
            "filter_selections",
        ]

        # Sync to Workspace object
        for attr in sync_attrs:
            value = None
            if attr == "pdb_file":
                value = self.pdb_file
            elif attr == "pdb_id":
                value = self.pdb_id
            elif attr == "structure":
                value = self.structure
            elif attr == "metadata":
                value = self.metadata
            elif attr == "filtered_structure":
                value = self.filtered_structure
            elif attr == "filter_selections":
                value = self.filter_selections

            if value is not None:
                # CRITICAL FIX: Don't overwrite file paths with BioPython objects
                if attr in ["filtered_structure", "repaired_structure"]:
                    # Check what's currently in the workspace
                    current_value = self.workspace.get(attr)

                    # If workspace has a file path and we're trying to set a BioPython object, skip
                    if isinstance(current_value, (str, Path)) and hasattr(value, 'get_atoms'):
                        # Don't overwrite the file path with the BioPython object
                        continue

                # Normal sync for everything else
                self.workspace.set(attr, value)

                # Special handling for pdb_file: also sync to original_pdb_file for consistency
                if attr == "pdb_file" and value is not None:
                    self.workspace.set("original_pdb_file", value)
                
    def _get_workspace(self):
        """Get the current workspace"""
        return self.workspace

    def check_residue_numbering(self):
        """Check if the first residue in each chain starts with number 1."""
        if not self.structure:
            return

        # Dictionary to store results
        numbering_info = {}

        # Iterate through all models, chains, and residues
        for model in self.structure:
            for chain in model:
                # Get the first residue of each chain
                residues = list(chain.get_residues())
                if residues:
                    first_res = residues[0]
                    # Check if the first residue number is not 1
                    if first_res.id[1] != 1:
                        numbering_info[chain.id] = {
                            "first_res_num": first_res.id[1],
                            "first_res_name": first_res.resname,
                            "warning": True,
                        }
                    else:
                        numbering_info[chain.id] = {
                            "first_res_num": first_res.id[1],
                            "first_res_name": first_res.resname,
                            "warning": False,
                        }

        # Store in workspace
        workspace = self._get_workspace()
        workspace.set("residue_numbering_info", numbering_info)

        # Check if there are any warnings
        has_warnings = any(
            info.get("warning", False) for info in numbering_info.values()
        )

        if has_warnings:
            from rich.panel import Panel

            warning_lines = []
            warning_lines.append(
                "[bold yellow]WARNING: Non-standard residue numbering detected[/bold yellow]"
            )
            warning_lines.append("")
            warning_lines.append(
                "The following chains do not start with residue number 1:"
            )

            for chain_id, info in numbering_info.items():
                if info.get("warning", False):
                    warning_lines.append(
                        f"  Chain {chain_id}: First residue is {info['first_res_name']} #{info['first_res_num']}"
                    )

            warning_lines.append("")
            warning_lines.append("[yellow]This may indicate:[/yellow]")
            warning_lines.append("• Missing N-terminal residues")
            warning_lines.append(
                "• Numbering that reflects the presence of a leader sequence in the immature protein"
            )
            warning_lines.append(
                "• PDB numbering that retains the original residue numbers from a reference sequence"
            )

            self.console.print(
                Panel(
                    "\n".join(warning_lines),
                    title="Residue Numbering Warning",
                    border_style="yellow",
                    expand=False,
                )
            )

    def load_pdb_interactive(self):
        """Interactive helper to load a PDB file."""

        pdb_file = load_pdb_interactive(processor=self)
        if pdb_file:
            self.pdb_file = pdb_file

            # Update workspace
            workspace = self._get_workspace()
            workspace.set("original_pdb_file", pdb_file)  # Primary key
            workspace.set("pdb_file", pdb_file)  # Backward compatibility

            # Load metadata
            self.metadata = PDBMetadataExtractor(pdb_file)
            workspace.set("metadata", self.metadata)

            # Try to load structure
            try:
                parser = PDBParser(QUIET=True)
                self.structure = parser.get_structure("protein", pdb_file)

                workspace.set("structure", self.structure)

                # Add residue numbering check right here
                self.check_residue_numbering()

            except Exception as e:
                self.console.print(
                    f"[yellow]Could not load structure: {str(e)}[/yellow]"
                )

            self.console.print(
                f"[green]Successfully loaded PDB file: {pdb_file}[/green]"
            )
            return True

        return False

    def run_metadata_module(self):
        """Interactive helper to view PDB metadata"""

        if not self.metadata:
            self.console.print("[yellow]No PDB metadata loaded yet.[/yellow]")
            return

        display_metadata_menu(self.metadata, self.console, processor=self)

    def run_filter_module(self):
        """Interactive helper to run the PDB filter module"""
        if not self.pdb_file:
            self.console.print("[yellow]No PDB file loaded yet.[/yellow]")
            return

        filter_tool = PDBFilterTool(self.pdb_file, processor=self)
        filtered_structure = filter_tool.interactive_filter()

        if filtered_structure:
            self.filtered_structure = filtered_structure
            self.filter_selections = filter_tool.filter_selections

            # Update workspace
            workspace = self._get_workspace()
            workspace.set("filtered_structure", filtered_structure)
            workspace.set("filter_selections", filter_tool.filter_selections)

            self.console.print("[green]Structure filtered successfully[/green]")

    def get_module_instance(self, name: str) -> Optional[ProcessingModule]:
        """Get or create a module instance and setup processor reference"""
        if name in self.module_instances:
            return self.module_instances[name]

        module = self.registry.get_module(name)
        if module:
            module.set_processor(self)
            module.initialize()
            self.module_instances[name] = module
            return module

        return None

    def run_module_menu(self, module: ProcessingModule):
        """Run a module's menu"""
        module_name = module.NAME
        menu_options = module.get_menu_options()

        workspace = self._get_workspace()
        debug_enabled = workspace.get("debug", False)
        if debug_enabled:
            self.console.print(
                f"[yellow]DEBUG: Module menu options for {module_name}: {list(menu_options.keys())}[/yellow]"
            )

        if not menu_options:
            self.console.print(
                f"[yellow]Module {module_name} has no menu options[/yellow]"
            )
            return

        while True:
            # Display module menu
            self.console.print(f"\n[bold]===== {module_name} Menu =====[/bold]")

            # Process can_process status
            can_process = module.can_process(workspace)
            if not can_process:
                requirements = module.get_workspace_requirements()
                self.console.print(
                    "[yellow]Warning: Missing workspace requirements:[/yellow]"
                )
                for req in requirements:
                    status = "✓" if workspace.get(req) is not None else "✗"
                    self.console.print(f"  {status} {req}")

            # Show menu options
            for i, (option, description) in enumerate(menu_options.items()):
                self.console.print(f"{i+1}. {description}")

            self.console.print(f"{len(menu_options)+1}. Back to main menu")

            # Build options_map for session context
            module_options_map = {
                str(i + 1): description
                for i, (_, description) in enumerate(menu_options.items())
            }
            module_options_map[str(len(menu_options) + 1)] = "Back to main menu"

            # Get user choice
            choice = prompt_with_context(
                self,
                "Enter your choice",
                choices=[str(i + 1) for i in range(len(menu_options) + 1)],
                default="1",
                module=f"{module_name} Menu",
                description=f"Select {module_name} option",
                options_map=module_options_map,
            )

            # Back to main menu
            if int(choice) == len(menu_options) + 1:
                break

            # Execute module option
            option_idx = int(choice) - 1
            option_key = list(menu_options.keys())[option_idx]

            # Check if can process first
            if not can_process and confirm_with_context(
                self,
                "[yellow]Module requirements not met. Try anyway?[/yellow]",
                default=False,
                module=module_name,
                description="Confirm run despite requirements not met",
            ):
                pass
            elif not can_process:
                continue

            # Handle the selected option
            module.handle_menu_option(option_key)

            # Force update from module instance variables to workspace
            if hasattr(module, "processor") and module.processor:
                # Update key processor variables to workspace
                self._sync_attributes_to_workspace()

        # Sync workspace to processor attributes
        self._sync_workspace_to_attributes()

    def run_workspace_options_menu(self):
        """Run the workspace options submenu"""
        workspace = self._get_workspace()

        while True:
            self.console.print("\n[bold]===== Workspace Options =====[/bold]")

            debug_enabled = workspace.get("debug", False)
            options = {
                "1": ("Show workspace status", "workspace_status"),
                "2": ("Show detailed workspace", "workspace_detailed"),
                "3": ("Save workspace", "save_workspace"),
                "4": ("Load workspace", "load_workspace"),
                "5": ("Reset workspace", "reset_workspace"),
                "6": (
                    f"Toggle debug mode (currently {'ON' if debug_enabled else 'OFF'})",
                    "toggle_debug",
                ),
                "7": ("Back to main menu", "back"),
            }

            # Display options
            for key, (description, _) in options.items():
                self.console.print(f"{key}. {description}")

            workspace_options_map = {k: v[0] for k, v in options.items()}

            # Get user choice
            choice = prompt_with_context(
                self,
                "\nEnter your choice",
                choices=list(options.keys()),
                default="1",
                module="Workspace Options",
                description="Select workspace option",
                options_map=workspace_options_map,
            )

            # Process choice
            option_name = options[choice][1]

            if option_name == "back":
                break
            elif option_name == "workspace_status":
                self.display_workspace_status(detailed=False)
            elif option_name == "workspace_detailed":
                self.display_workspace_status(detailed=True)
            elif option_name == "save_workspace":
                filename = prompt_with_context(
                    self,
                    "Enter filename to save workspace",
                    default="workspace.json",
                    module="Workspace Options",
                    description="Filename to save workspace to",
                )
                self.save_workspace(filename)
            elif option_name == "load_workspace":
                filename = prompt_with_context(
                    self,
                    "Enter filename to load workspace",
                    default="workspace.json",
                    module="Workspace Options",
                    description="Filename to load workspace from",
                )
                if os.path.exists(filename):
                    self.load_workspace(filename)
                else:
                    self.console.print(f"[red]File not found: {filename}[/red]")
            elif option_name == "reset_workspace":
                self.reset_workspace()
            elif option_name == "toggle_debug":
                self.workspace.set("debug", not self.workspace.get("debug", False))

                debug_status = "ON" if workspace.get("debug", False) else "OFF"
                self.console.print(f"[green]Debug mode is now {debug_status}[/green]")

    def run_workflow_step(self, step_name: str) -> bool:
        """Run a specific workflow step"""
        module = self.get_module_instance(step_name)
        if not module:
            self.console.print(f"[red]Module {step_name} not found[/red]")
            return False

        workspace = self._get_workspace()

        # Special case for PDB Loader with pdb_id
        if step_name == "PDB Loader" and workspace.get("pdb_id"):
            pdb_id = workspace.get("pdb_id")

            # Download the PDB file
            pdb_file = download_pdb(pdb_id)
            if pdb_file:
                self.console.print(f"[green]Downloaded PDB file: {pdb_file}[/green]")
                self.pdb_file = pdb_file
                self.workspace.set("original_pdb_file", pdb_file)  # Primary key
                self.workspace.set("pdb_file", pdb_file)  # Backward compatibility

                # Load metadata
                self.metadata = PDBMetadataExtractor(pdb_file)
                self.workspace.set("metadata", self.metadata)

                # Try to load structure
                try:
                    from Bio.PDB import PDBParser

                    parser = PDBParser(QUIET=True)
                    self.structure = parser.get_structure("protein", pdb_file)

                    self.workspace.set("structure", self.structure)

                    # Add residue numbering check right here
                    self.check_residue_numbering()

                except Exception as e:
                    self.console.print(
                        f"[yellow]Could not load structure: {str(e)}[/yellow]"
                    )

                self.console.print(f"[green]Successfully loaded PDB {pdb_id}[/green]")
                return True

        # Check if module can process current workspace
        if not module.can_process(workspace):
            requirements = module.get_workspace_requirements()
            self.console.print(
                f"[yellow]Cannot run {step_name} - missing requirements:[/yellow]"
            )
            for req in requirements:
                status = "✓" if workspace.get(req) is not None else "✗"
                self.console.print(f"  {status} {req}")
            return False

        # Process workspace using the module
        with Status(f"Running {step_name}..."):
            updated_workspace = module.process(workspace)

            # No need to update since the Workspace object is passed by reference

        self.console.print(f"[green]Successfully ran {step_name}[/green]")

        # Sync workspace to processor attributes
        self._sync_workspace_to_attributes()

        return True

    def display_workspace_status(self, detailed=False):
        """Display the current workspace status"""
        self.console.print("\n[bold]===== Workspace Status =====[/bold]")

        if detailed:
            self.display_detailed_workspace()
            return

        workspace = self._get_workspace()

        # Use Workspace's dump method
        self.workspace.dump(self.console)

    def display_detailed_workspace(self):
        """Display detailed information about the workspace contents"""
        self.console.print("\n[bold]===== Detailed Workspace Contents =====[/bold]")

        workspace = self._get_workspace()

        # PDB File
        pdb_file = workspace.get("pdb_file")
        self.console.print(
            Panel(
                f"[bright_blue]Path:[/bright_blue] {pdb_file if pdb_file else 'Not loaded'}",
                title="PDB File",
                border_style="green" if pdb_file else "red",
            )
        )

        # PDB ID
        pdb_id = workspace.get("pdb_id")
        if pdb_id:
            self.console.print(
                Panel(
                    f"[bright_blue]ID:[/bright_blue] {pdb_id}", title="PDB ID", border_style="green"
                )
            )

        # Try module-specific display for workspace items
        displayed_keys = set()
        for key, value in workspace.items():
            if key not in ["debug"]:  # Skip debug flag
                # Try to find a module that can display this data
                for module_name, module in self.module_instances.items():
                    try:
                        # Check if module has and handles this key
                        if hasattr(
                            module, "get_workspace_display"
                        ) and module.get_workspace_display(key, value, self.console):
                            displayed_keys.add(key)
                            break
                    except Exception as e:
                        # Log the error but continue processing
                        self.console.print(
                            f"[yellow]Warning: Error in display hook for module {module_name}: {str(e)}[/yellow]"
                        )

        # Display Structure and other standard items using default display if not handled by modules
        structure = workspace.get("structure")
        if structure and "structure" not in displayed_keys:
            # Safely handle structure
            if hasattr(structure, "__iter__") and hasattr(structure, "get_models"):
                # It's a proper structure object
                structure_info = []
                structure_info.append(f"[bright_blue]Model Count:[/bright_blue] {len(structure)}")

                # Count chains and residues per model
                for model_idx, model in enumerate(structure):
                    chain_count = len(list(model))
                    residue_count = sum(len(list(chain)) for chain in model)
                    atom_count = sum(
                        sum(len(list(residue)) for residue in chain) for chain in model
                    )
                    structure_info.append(
                        f"[bright_blue]Model {model_idx}:[/bright_blue] {chain_count} chains, {residue_count} residues, {atom_count} atoms"
                    )

                    # Chain details
                    for chain in model:
                        chain_id = chain.id
                        residue_count = len(list(chain))
                        atom_count = sum(len(list(residue)) for residue in chain)
                        structure_info.append(
                            f"  [bright_blue]Chain {chain_id}:[/bright_blue] {residue_count} residues, {atom_count} atoms"
                        )

                self.console.print(
                    Panel(
                        "\n".join(structure_info),
                        title="Structure",
                        border_style="green",
                    )
                )
            else:
                # Not a proper structure object
                self.console.print(
                    Panel(
                        f"[bright_blue]Type:[/bright_blue] {type(structure).__name__}",
                        title="Structure",
                        border_style="yellow",
                    )
                )

            displayed_keys.add("structure")

        # Filtered Structure
        filtered_structure = workspace.get("filtered_structure")
        if filtered_structure and "filtered_structure" not in displayed_keys:
            filtered_info = []

            # Check if it's a Path object or string (file path)
            if isinstance(filtered_structure, (str, Path)):
                # It's a file path, display info about the file
                filtered_file_path = Path(filtered_structure)
                filtered_info.append(f"[bright_blue]File:[/bright_blue] {filtered_file_path}")

                if filtered_file_path.exists():
                    # Get file size
                    file_size = filtered_file_path.stat().st_size
                    filtered_info.append(
                        f"[bright_blue]File Size:[/bright_blue] {file_size / 1024:.1f} KB"
                    )

                    # Count atoms, chains, etc.
                    try:
                        atom_count = 0
                        hetatm_count = 0
                        chains = set()
                        residues = set()

                        with open(filtered_file_path, "r") as f:
                            for line in f:
                                if line.startswith("ATOM"):
                                    atom_count += 1
                                    chains.add(line[21])
                                    residues.add(
                                        (
                                            line[21],
                                            line[17:20].strip(),
                                            int(line[22:26].strip()),
                                        )
                                    )
                                elif line.startswith("HETATM"):
                                    hetatm_count += 1
                                    chains.add(line[21])
                                    residues.add(
                                        (
                                            line[21],
                                            line[17:20].strip(),
                                            int(line[22:26].strip()),
                                        )
                                    )

                        filtered_info.append(f"[bright_blue]Chains:[/bright_blue] {len(chains)}")
                        filtered_info.append(f"[bright_blue]Residues:[/bright_blue] {len(residues)}")
                        filtered_info.append(
                            f"[bright_blue]Atoms:[/bright_blue] {atom_count} ATOM + {hetatm_count} HETATM = {atom_count + hetatm_count} total"
                        )

                        # List chains
                        if chains:
                            filtered_info.append(
                                f"[bright_blue]Chain IDs:[/bright_blue] {', '.join(sorted(chains))}"
                            )
                    except Exception as e:
                        filtered_info.append(
                            f"[yellow]Error reading file details: {str(e)}[/yellow]"
                        )
                else:
                    filtered_info.append(f"[yellow]File not found[/yellow]")
            # Check if it's a list of PDB lines
            elif isinstance(filtered_structure, list) and all(
                isinstance(line, str) for line in filtered_structure[:10]
            ):
                # It's a list of PDB lines
                filtered_info.append(
                    f"[bright_blue]Format:[/bright_blue] PDB Lines in memory ({len(filtered_structure)} lines)"
                )

                # Count atoms, chains, etc.
                try:
                    atom_count = 0
                    hetatm_count = 0
                    chains = set()
                    residues = set()

                    for line in filtered_structure:
                        if line.startswith("ATOM"):
                            atom_count += 1
                            chains.add(line[21])
                            try:
                                residues.add(
                                    (
                                        line[21],
                                        line[17:20].strip(),
                                        int(line[22:26].strip()),
                                    )
                                )
                            except (IndexError, ValueError):
                                pass
                        elif line.startswith("HETATM"):
                            hetatm_count += 1
                            chains.add(line[21])
                            try:
                                residues.add(
                                    (
                                        line[21],
                                        line[17:20].strip(),
                                        int(line[22:26].strip()),
                                    )
                                )
                            except (IndexError, ValueError):
                                pass

                    filtered_info.append(f"[bright_blue]Chains:[/bright_blue] {len(chains)}")
                    filtered_info.append(f"[bright_blue]Residues:[/bright_blue] {len(residues)}")
                    filtered_info.append(
                        f"[bright_blue]Atoms:[/bright_blue] {atom_count} ATOM + {hetatm_count} HETATM = {atom_count + hetatm_count} total"
                    )

                    # List chains
                    if chains:
                        filtered_info.append(
                            f"[bright_blue]Chain IDs:[/bright_blue] {', '.join(sorted(chains))}"
                        )
                except Exception as e:
                    filtered_info.append(
                        f"[yellow]Error analyzing PDB lines: {str(e)}[/yellow]"
                    )
            else:
                # Assume it's a BioPython structure
                try:
                    # Check if it acts like a BioPython structure
                    if hasattr(filtered_structure, "__iter__") and hasattr(
                        filtered_structure, "get_models"
                    ):
                        filtered_info.append(
                            f"[bright_blue]Model Count:[/bright_blue] {len(filtered_structure)}"
                        )

                        # Count chains and residues per model
                        for model_idx, model in enumerate(filtered_structure):
                            chain_count = len(list(model))
                            residue_count = sum(len(list(chain)) for chain in model)
                            atom_count = sum(
                                sum(len(list(residue)) for residue in chain)
                                for chain in model
                            )
                            filtered_info.append(
                                f"[bright_blue]Model {model_idx}:[/bright_blue] {chain_count} chains, {residue_count} residues, {atom_count} atoms"
                            )

                            # Chain details
                            for chain in model:
                                chain_id = chain.id
                                residue_count = len(list(chain))
                                atom_count = sum(
                                    len(list(residue)) for residue in chain
                                )
                                filtered_info.append(
                                    f"  [bright_blue]Chain {chain_id}:[/bright_blue] {residue_count} residues, {atom_count} atoms"
                                )
                    else:
                        # It's neither a file path nor a proper BioPython structure
                        filtered_info.append(
                            f"[yellow]Type: {type(filtered_structure).__name__} (format not recognized)[/yellow]"
                        )
                except (TypeError, AttributeError) as e:
                    # It's neither a file path nor a proper BioPython structure
                    filtered_info.append(
                        f"[yellow]Type: {type(filtered_structure).__name__} (error: {str(e)})[/yellow]"
                    )

            self.console.print(
                Panel(
                    "\n".join(filtered_info),
                    title="Filtered Structure",
                    border_style="green",
                )
            )
            displayed_keys.add("filtered_structure")

        # Filter Selections
        filter_selections = workspace.get("filter_selections")
        if filter_selections and "filter_selections" not in displayed_keys:
            selection_info = []
            selection_info.append(f"[bright_blue]Chain Count:[/bright_blue] {len(filter_selections)}")

            # Display each selection
            for chain_id, selection in filter_selections.items():
                if isinstance(selection, dict):
                    selection_str = ", ".join(
                        [f"{k}: {v}" for k, v in selection.items()]
                    )
                elif isinstance(selection, list):
                    selection_str = ", ".join(map(str, selection))
                else:
                    selection_str = str(selection)

                if len(selection_str) > 100:
                    selection_str = selection_str[:97] + "..."

                selection_info.append(f"[bright_blue]Chain {chain_id}:[/bright_blue] {selection_str}")

            self.console.print(
                Panel(
                    "\n".join(selection_info),
                    title="Filter Selections",
                    border_style="green",
                )
            )
            displayed_keys.add("filter_selections")

        # Handle other common data types if not already displayed by modules
        # Completeness Results
        completeness_results = workspace.get("completeness_results")
        if completeness_results and "completeness_results" not in displayed_keys:
            comp_info = []

            # Display missing residues
            missing_residues = completeness_results.get("missing_residues", [])
            if missing_residues:
                comp_info.append(
                    f"[bright_blue]Missing Residues Count:[/bright_blue] {len(missing_residues)}"
                )
                # Show a subset if there are many
                if len(missing_residues) > 5:
                    for i in range(min(5, len(missing_residues))):
                        res = missing_residues[i]
                        comp_info.append(
                            f"  {res.get('chain', '?')}: {res.get('resname', '?')} {res.get('resseq', '?')}"
                        )
                    comp_info.append(f"  ... and {len(missing_residues) - 5} more")
                else:
                    for res in missing_residues:
                        comp_info.append(
                            f"  {res.get('chain', '?')}: {res.get('resname', '?')} {res.get('resseq', '?')}"
                        )

            # Display missing atoms
            missing_atoms = completeness_results.get("missing_atoms", [])
            if missing_atoms:
                comp_info.append(
                    f"[bright_blue]Missing Atoms Count:[/bright_blue] {len(missing_atoms)}"
                )
                # Show a subset if there are many
                if len(missing_atoms) > 5:
                    for i in range(min(5, len(missing_atoms))):
                        atom = missing_atoms[i]
                        comp_info.append(
                            f"  {atom.get('chain', '?')}: {atom.get('resname', '?')} {atom.get('resseq', '?')} - {atom.get('name', '?')}"
                        )
                    comp_info.append(f"  ... and {len(missing_atoms) - 5} more")
                else:
                    for atom in missing_atoms:
                        comp_info.append(
                            f"  {atom.get('chain', '?')}: {atom.get('resname', '?')} {atom.get('resseq', '?')} - {atom.get('name', '?')}"
                        )

            # Display alternate locations
            alt_locs = completeness_results.get("alt_locs", [])
            if alt_locs:
                comp_info.append(
                    f"[bright_blue]Alternate Locations Count:[/bright_blue] {len(alt_locs)}"
                )
                # Show a subset if there are many
                if len(alt_locs) > 5:
                    for i in range(min(5, len(alt_locs))):
                        alt = alt_locs[i]
                        comp_info.append(
                            f"  {alt.get('chain', '?')}: {alt.get('resname', '?')} {alt.get('resseq', '?')} - {alt.get('altloc', '?')}"
                        )
                    comp_info.append(f"  ... and {len(alt_locs) - 5} more")
                else:
                    for alt in alt_locs:
                        comp_info.append(
                            f"  {alt.get('chain', '?')}: {alt.get('resname', '?')} {alt.get('resseq', '?')} - {alt.get('altloc', '?')}"
                        )

            if comp_info:
                self.console.print(
                    Panel(
                        "\n".join(comp_info),
                        title="Completeness Results",
                        border_style="green",
                    )
                )
                displayed_keys.add("completeness_results")

        # Metal Sites
        metal_sites = workspace.get("metal_sites")
        if metal_sites and "metal_sites" not in displayed_keys:
            metal_info = []
            metal_info.append(f"[bright_blue]Metal Sites Count:[/bright_blue] {len(metal_sites)}")

            for i, site in enumerate(metal_sites):
                metal_ion = self._get_site_attribute(site, "metal_element", "Unknown")
                chain = self._get_site_attribute(site, "metal_chain", "?")
                resseq = self._get_site_attribute(site, "metal_resid", "?")
                geometry = self._get_site_attribute(site, "geometry", "Unknown")
                coordination = len(self._get_site_attribute(site, "ligands", []))

                metal_info.append(
                    f"[bright_blue]Site {i+1}:[/bright_blue] {metal_ion} (Chain {chain}: {resseq})"
                )
                metal_info.append(
                    f"  Geometry: {geometry}, Coordination: {coordination}"
                )

                # Ligands (if available)
                ligands = self._get_site_attribute(site, "ligands", [])
                if ligands:
                    metal_info.append("  Ligands:")
                    for j, ligand in enumerate(ligands[:5]):  # Show only first 5
                        if hasattr(ligand, "get"):
                            metal_info.append(
                                f"    {ligand.get('chain', '?')}: {ligand.get('resname', '?')} {ligand.get('resid', '?')} - {ligand.get('atom_name', '?')}"
                            )
                        elif hasattr(ligand, "chain"):
                            metal_info.append(
                                f"    {getattr(ligand, 'chain', '?')}: {getattr(ligand, 'resname', '?')} {getattr(ligand, 'resid', '?')} - {getattr(ligand, 'atom_name', '?')}"
                            )
                    if len(ligands) > 5:
                        metal_info.append(f"    ... and {len(ligands) - 5} more")

            self.console.print(
                Panel("\n".join(metal_info), title="Metal Sites", border_style="green")
            )
            displayed_keys.add("metal_sites")

        # Other workspace items not handled by modules or standard handlers
        other_items = {}
        for key, value in workspace.items():
            if key not in displayed_keys and key != "debug":
                other_items[key] = value

        if other_items:
            other_info = []
            for key, value in other_items.items():
                value_type = type(value).__name__

                if value is None:
                    other_info.append(f"[bright_blue]{key}:[/bright_blue] None")
                elif isinstance(value, (str, int, float, bool)):
                    other_info.append(f"[bright_blue]{key}:[/bright_blue] {value} ({value_type})")
                elif hasattr(value, "__len__"):
                    if hasattr(value, "keys"):  # Dictionary-like
                        other_info.append(
                            f"[bright_blue]{key}:[/bright_blue] {len(value)} items ({value_type})"
                        )
                        # Show a few keys as example
                        keys = list(value.keys())[:5]
                        if keys:
                            other_info.append(
                                f"  Keys: {', '.join(str(k) for k in keys)}{' ...' if len(value) > 5 else ''}"
                            )
                    else:  # List-like
                        other_info.append(
                            f"[bright_blue]{key}:[/bright_blue] {len(value)} elements ({value_type})"
                        )
                        # Try to show a few elements
                        if len(value) > 0 and all(
                            isinstance(item, (str, int, float, bool))
                            for item in list(value)[:5]
                        ):
                            items = list(value)[:5]
                            other_info.append(
                                f"  Items: {', '.join(str(item) for item in items)}{' ...' if len(value) > 5 else ''}"
                            )
                else:
                    other_info.append(f"[bright_blue]{key}:[/bright_blue] {value_type}")

            if other_info:
                self.console.print(
                    Panel(
                        "\n".join(other_info),
                        title="Other Workspace Items",
                        border_style="green",
                    )
                )

    def _get_site_attribute(self, site, attribute, default=None):
        """Safely get an attribute from a metal site object, handling both dict and object styles."""
        # Try dictionary-style access first
        if isinstance(site, dict):
            return site.get(attribute, default)

        # Then try attribute-style access
        if hasattr(site, attribute):
            return getattr(site, attribute)

        # Fall back to default
        return default

    def save_workspace(self, filename: str):
        """Save workspace to a file"""
        # Create a serializable version of the workspace
        serializable = {}

        workspace = self._get_workspace()

        # Save basic attributes
        serializable["pdb_file"] = workspace.get("pdb_file")

        # Handle filter_selections
        if workspace.get("filter_selections"):
            serializable["filter_selections"] = workspace.get("filter_selections")

        # Save completeness results if available
        if workspace.get("completeness_results"):
            serializable["completeness_results"] = workspace.get("completeness_results")

        # Write to file
        try:
            with open(filename, "w") as f:
                json.dump(serializable, f, indent=2)
            self.console.print(f"[green]Workspace saved to {filename}[/green]")
            return True
        except Exception as e:
            self.console.print(f"[red]Error saving workspace: {str(e)}[/red]")
            return False

    def load_workspace(self, filename: str):
        """Load workspace from a file"""
        try:
            with open(filename, "r") as f:
                data = json.load(f)

            # Update workspace with loaded data
            workspace = self._get_workspace()

            for key, value in data.items():
                self.workspace.set(key, value)

            # Reload key structures if PDB file is specified
            if "pdb_file" in data and data["pdb_file"]:
                pdb_file = data["pdb_file"]
                if os.path.exists(pdb_file):
                    # Try to reload structures and metadata
                    try:
                        # Get PDB loader module
                        pdb_loader = self.get_module_instance("PDB Loader")
                        if pdb_loader:
                            # Process to reload structure and metadata
                            self.run_workflow_step("PDB Loader")
                    except Exception as e:
                        self.console.print(
                            f"[yellow]Warning: Could not reload structures: {str(e)}[/yellow]"
                        )

            # Update processor attributes
            self._sync_workspace_to_attributes()

            self.console.print(f"[green]Workspace loaded from {filename}[/green]")
            return True
        except Exception as e:
            self.console.print(f"[red]Error loading workspace: {str(e)}[/red]")
            return False

    def reset_workspace(self):
        """Reset the workspace to initial state"""
        if confirm_with_context(
            self,
            "[yellow]Are you sure you want to reset the workspace?[/yellow]",
            default=False,
            module="Workspace Options",
            description="Confirm reset workspace",
        ):
            workspace = self._get_workspace()
            debug_setting = workspace.get("debug", False)

            self.workspace.clear()
            self.workspace.set("debug", debug_setting)

            self.pdb_file = None
            self.pdb_id = None
            self.structure = None
            self.metadata = None
            self.filtered_structure = None
            self.filter_selections = None
            self.console.print("[green]Workspace has been reset[/green]")
            return True
        return False

    def run_main_menu(self):
        """Run the main menu interface using command pattern"""
        from proprep.application.menu_commands import MainMenuCommand, WorkflowMenuCommand

        # Check for workflow shortcuts (e.g., --analysis flag)
        if self.workspace.get("jump_to_analysis", False):
            # Clear the flag so it doesn't keep jumping
            self.workspace.set("jump_to_analysis", False)

            # Jump directly to MD Manager analysis
            md_manager = self.get_module_instance("Molecular Dynamics Manager")
            if md_manager:
                self.console.print("\n[bold bright_blue]Jumping to Simulation Analysis...[/bold bright_blue]\n")
                try:
                    md_manager.handle_menu_option("analyze")
                    self._sync_attributes_to_workspace()
                except Exception as e:
                    self.console.print(f"[red]Error during analysis: {e}[/red]")
                    import traceback
                    traceback.print_exc()
            else:
                self.console.print("[red]Error: Molecular Dynamics Manager not available[/red]")

            # After analysis completes, show the menu (user can exit if they want)
            # Fall through to normal menu below

        # Check for PDB view shortcut (--pdbview flag)
        if self.workspace.get("jump_to_pdbview", False):
            # Clear the flag so it doesn't keep jumping
            self.workspace.set("jump_to_pdbview", False)
            target = self.workspace.get("pdbview_target", "")

            if target:
                loader = self.get_module_instance("Structure Loader")
                workspace = self._get_workspace()
                pdb_file = None

                if len(target) == 4 and target.isalnum():
                    # PDB ID — download via StructureLoader
                    loader._download_and_load_pdb(target.upper(), workspace)
                    pdb_file = workspace.get("rcsb_pdb_file")
                else:
                    # Local file — load via StructureLoader
                    loader._load_local_file_by_path(target, workspace)
                    pdb_file = workspace.get("local_pdb_file") or workspace.get("alphafold_pdb_file")

                # Launch viewer if we have a valid file
                if pdb_file:
                    viewer = self.get_module_instance("Structure Viewer")
                    if viewer:
                        try:
                            viewer.quick_launch([pdb_file])
                        except Exception as e:
                            self.console.print(f"[red]Error launching viewer: {e}[/red]")
                    else:
                        self.console.print("[red]Error: Structure Viewer module not available[/red]")

            # Fall through to normal menu below

        # Loop to allow switching between menu modes
        while True:
            # Check which menu mode to use
            menu_mode = self.workspace.get("menu_mode", "full-menu")

            if menu_mode == "workflow":
                menu_command = WorkflowMenuCommand(self)
            else:
                menu_command = MainMenuCommand(self)

            menu_command.execute_with_error_handling()

            # Check if menu mode changed during execution
            new_menu_mode = self.workspace.get("menu_mode", "full-menu")
            if new_menu_mode == menu_mode:
                # Menu mode didn't change, so user wants to exit
                break
            # Otherwise, loop continues to reload with new menu mode

    def cleanup(self):
        """Cleanup processor resources"""
        self.registry.cleanup()

    def show_breadcrumbs(self):
        """Display current breadcrumb trail."""
        trail = self.command_history.get_breadcrumb_trail()
        self.console.print(f"[grey50]Navigation: {trail}[/grey50]")

    def show_command_history(self, count: int = 10):
        """Display recent command history."""
        recent_commands = self.command_history.get_recent_commands(count)
        if not recent_commands:
            self.console.print("[yellow]No command history available[/yellow]")
            return

        self.console.print(f"\n[bold]===== Recent Commands (last {count}) =====[/bold]")
        for i, cmd in enumerate(recent_commands, 1):
            timestamp = cmd.timestamp.strftime("%H:%M:%S")
            self.console.print(f"{i:2d}. [{timestamp}] {cmd}")

    def clear_breadcrumbs(self):
        """Clear all breadcrumbs."""
        self.command_history.clear_breadcrumbs()
        self.console.print("[green]Breadcrumbs cleared[/green]")
