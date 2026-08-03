"""
Redox Site Detector Module

Standalone module for detecting redox-active sites in protein structures.
This module provides early access to redox site detection in the workflow,
allowing detected sites to be used for structure alignment and other analyses.
"""

import logging
from typing import Dict, List, Any

from rich.console import Console

from proprep.utils.module_registry import ProcessingModule, register_module
from proprep.utils.prompts import prompt_with_context

# Import redox detector
try:
    from .comprehensive_redox_detector import ComprehensiveRedoxDetector
    REDOX_DETECTION_AVAILABLE = True
except ImportError:
    REDOX_DETECTION_AVAILABLE = False
    logger = logging.getLogger(__name__)
    logger.warning("Redox detection module not available.")

logger = logging.getLogger(__name__)


@register_module
class RedoxSiteDetectorModule(ProcessingModule):
    """Module for detecting redox-active sites in protein structures"""

    NAME = "Redox Site Detector"
    DESCRIPTION = "Detect redox-active groups in protein structures"
    VERSION = "1.0.0"
    CATEGORY = "analysis"
    PRIORITY = 1.5  # Between PDB Loader (1) and Compare stage modules (10+)

    def __init__(self):
        """Initialize the Redox Site Detector module"""
        super().__init__()
        self.detector = None
        self.detected_sites = None

    @property
    def console(self) -> Console:
        """Get console from processor or create new one"""
        if self.processor and hasattr(self.processor, 'console'):
            return self.processor.console
        else:
            return Console()

    def get_workspace_requirements(self) -> List[str]:
        """Get workspace requirements.

        A loaded structure is needed only to *detect* redox sites. The other
        menu action — *importing* a redox-sites JSON from a prior session —
        needs nothing in the workspace, and that import is the path back into a
        resumed run where only a prmtop/rst7 pair was reloaded (no PDB). So the
        module is always reachable; the detect action validates its own input
        at run time.
        """
        return []

    def get_workspace_outputs(self) -> List[str]:
        """Get workspace outputs"""
        return [
            "detected_redox_sites",
            "redox_transformer_mappings",
        ]

    def can_process(self, workspace: Dict[str, Any]) -> bool:
        """Check if the module can process the current workspace.

        Available whenever the detection code is importable. We intentionally
        do NOT require a loaded PDB: one menu action detects sites (needs a
        structure, checked at run time), the other imports a redox-sites JSON
        from a prior session (needs nothing). The latter is how a resumed run
        — reloaded from a prmtop/rst7 pair with no PDB — gets its redox sites
        back, so gating on a PDB here would wrongly hide that path.
        """
        return REDOX_DETECTION_AVAILABLE

    def _get_workspace_value(self, workspace, *keys):
        """Helper to get value from workspace trying multiple keys"""
        for key in keys:
            value = workspace.get(key)
            if value is not None:
                return value
        return None

    def process(self, workspace: Dict[str, Any]) -> Dict[str, Any]:
        """Process the workspace"""
        # Check if already detected
        if workspace.get("detected_redox_sites"):
            self.detected_sites = workspace.get("detected_redox_sites")

        return workspace

    def get_menu_options(self) -> Dict[str, str]:
        """Get module menu options"""
        return {
            "import": "Import redox sites from JSON file",
            "detect": "Detect redox-active sites",
            "view": "View detected site summaries",
            "edit": "Edit detected sites",
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

        # Check if redox sites have been detected
        detected_sites = workspace.get("detected_redox_sites")
        sites_detected = detected_sites is not None and len(detected_sites) > 0

        # Option 1: Import redox sites from JSON - always available or completed
        options.append(MenuOption(
            key="1",
            description="Import redox sites from JSON file",
            status=OptionStatus.COMPLETED if sites_detected else OptionStatus.AVAILABLE
        ))

        # Option 2: Detect redox sites - requires a loaded structure. The
        # detection runs against a parsed Structure object, so gate on the
        # same "is a structure available" signal; otherwise the option
        # showed ✓ but failed on selection with "load a structure first".
        from proprep.utils.structure_selector import StructureSelector
        has_structure = StructureSelector(
            workspace, self.console
        ).get_structure_status().get("has_any", False)
        if sites_detected:
            detect_status, detect_dep = OptionStatus.COMPLETED, ""
        elif has_structure:
            detect_status, detect_dep = OptionStatus.AVAILABLE, ""
        else:
            detect_status, detect_dep = OptionStatus.BLOCKED, "Load a structure first"
        options.append(MenuOption(
            key="2",
            description="Detect redox-active sites",
            status=detect_status,
            dependency_text=detect_dep,
        ))

        # Option 3: View detected site summaries - requires detection to be done
        options.append(MenuOption(
            key="3",
            description="View detected site summaries",
            status=OptionStatus.READY if sites_detected else OptionStatus.BLOCKED,
            dependency_text="[Need to detect or import sites first] ○" if not sites_detected else ""
        ))

        # Option 4: Edit site(s) - requires detection to be done
        options.append(MenuOption(
            key="4",
            description="Edit detected sites",
            status=OptionStatus.READY if sites_detected else OptionStatus.BLOCKED,
            dependency_text="[Need to detect or import sites first] ○" if not sites_detected else ""
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
        detected_sites = workspace.get("detected_redox_sites")

        if not detected_sites:
            return "Import from JSON (option 1) or detect sites (option 2) to identify redox-active residues"
        else:
            num_sites = len(detected_sites)
            return f"Found {num_sites} redox site(s). View with option 3, edit with option 4, or press [m] to return to the main menu"

    def handle_menu_option(self, option: str) -> bool:
        """Handle a menu option selection"""
        try:
            if option == "import":
                return self._import_redox_sites()
            elif option == "detect":
                return self._detect_redox_sites_interactive()
            elif option == "view":
                return self._view_detected_sites()
            elif option == "edit":
                return self._edit_redox_sites()
            else:
                logger.warning(f"Unknown menu option: {option}")
                return False

        except Exception as e:
            logger.error(f"Error handling menu option '{option}': {str(e)}")
            return False

    def _show_workflow_overview(self):
        """Display overview of the redox site detection workflow"""
        from rich.panel import Panel
        from rich.text import Text

        # Styling: phase titles are bold blue (the menu header colour); the
        # role column makes the division of labour scannable -- orange "you"
        # (you act) vs purple "ProPrep" (done automatically). Description text
        # uses the terminal default so it stays legible on light or dark
        # backgrounds.
        from rich.table import Table

        TITLE = "bold blue"
        YOU = "bold dark_orange3"      # you act
        PP = "bold dark_violet"        # ProPrep does it automatically

        guide = Table.grid(padding=(0, 2, 0, 0), expand=True)
        guide.add_column(style=TITLE, no_wrap=True)   # phase title
        guide.add_column(no_wrap=True)                # role: you / ProPrep
        guide.add_column(overflow="fold", ratio=1)    # what happens
        phases = [
            ("1. Configuration", "you",
             "set scan parameters, cutoffs, and which metals / residues to include"),
            ("2. Inventory Scan", "ProPrep",
             "finds cofactors, metals, redox amino acids and disulfides; summarizes"),
            ("3. Selection & Grouping", "you",
             "choose centers and group into sites (ProPrep auto-pairs disulfides)"),
            ("4. Templates (optional)", "you",
             "build templates that batch-apply to similar site types"),
            ("5. Site Refinement", "you",
             "set each site's boundaries and bonds (ProPrep applies templates)"),
            ("6. Review & Export", "you",
             "review, then save JSON or export PDBs"),
        ]
        for i, (phase, who, what) in enumerate(phases):
            if i:
                guide.add_row("", "", "")  # blank line between items
            label = "You" if who == "you" else who
            role = Text(label, style=(YOU if who == "you" else PP))
            guide.add_row(phase, role, what)

        panel = Panel(guide, title="[bold blue]How RedoxSite Detection Works[/bold blue]",
                      border_style="blue", padding=(1, 2), expand=True)
        self.console.print(panel)
        self.console.print()

    def _import_redox_sites(self) -> bool:
        """Import redox sites from a JSON file"""
        from proprep.utils.prompts import confirm_with_context
        import os
        import json

        if not REDOX_DETECTION_AVAILABLE:
            self.console.print("[red]Error: Redox detection module not available[/red]")
            return False

        workspace = self.processor.workspace if self.processor else {}
        working_dir = workspace.get("working_directory", os.getcwd())

        self.console.print("\n[bold cyan]Import Redox Sites from JSON[/bold cyan]\n")

        # Use the file browser from comprehensive_redox_detector
        from .comprehensive_redox_detector import display_json_file_menu

        selected_file = display_json_file_menu(
            directory=working_dir,
            console=self.console,
            processor=self.processor
        )

        if not selected_file:
            self.console.print("[grey50]Import cancelled[/grey50]")
            return False

        self.console.print(f"\n[cyan]Loading: {os.path.basename(selected_file)}[/cyan]")

        # Import the JSON file using the function from comprehensive_redox_detector
        try:
            from .comprehensive_redox_detector import _import_from_json

            sites, transformer_mappings = _import_from_json(selected_file)

            if sites:
                self.console.print(f"[green]✓ Successfully imported {len(sites)} redox site(s)[/green]")

                # Store in workspace
                workspace["detected_redox_sites"] = sites
                self.detected_sites = sites

                # Restore the site-type → transformer-name map so the
                # transformation manager's auto-assign branch fires for
                # an imported session, matching same-session behavior.
                if transformer_mappings:
                    workspace["redox_transformer_mappings"] = transformer_mappings
                    self.console.print(
                        f"[grey50]Restored transformer mappings for "
                        f"{len(transformer_mappings)} site type(s)[/grey50]"
                    )

                # Show brief summary
                self.console.print("\n[bold]Imported sites:[/bold]")
                for i, site in enumerate(sites, 1):
                    num_centers = len(site.centers) if hasattr(site, 'centers') else 0
                    num_atoms = len(site.atoms)
                    self.console.print(f"  Site {i}: {num_centers} center(s), {num_atoms} atoms")

                return True
            else:
                self.console.print("[yellow]No sites found in JSON file[/yellow]")
                return False

        except FileNotFoundError as e:
            self.console.print(f"[red]Error: {e}[/red]")
            return False
        except json.JSONDecodeError as e:
            self.console.print(f"[red]Error parsing JSON file: {e}[/red]")
            return False
        except Exception as e:
            self.console.print(f"[red]Error importing redox sites: {e}[/red]")
            logger.error(f"Import error: {e}", exc_info=True)
            return False

    def _detect_redox_sites_interactive(self) -> bool:
        """Run interactive redox site detection with model selection"""
        if not REDOX_DETECTION_AVAILABLE:
            self.console.print("[red]Error: Redox detection module not available[/red]")
            return False

        workspace = self.processor.workspace if self.processor else {}

        # Get structure using the new structure selector
        # This will automatically find and allow user to choose between
        # original_structure, alphafold_structure, etc.
        from proprep.utils.structure_selector import get_interactive_structure_object

        structure, file_path_key = get_interactive_structure_object(
            workspace,
            self.console,
            processor=self.processor  # Pass processor for session recording
        )
        if not structure:
            self.console.print("[yellow]No structure loaded in workspace[/yellow]")
            self.console.print("[grey50]Tip: Use 'Structure Loader' to load a PDB or AlphaFold structure[/grey50]")
            return False

        # Show workflow overview
        self._show_workflow_overview()

        # Get the PDB file path that corresponds to the selected structure
        pdb_file = workspace.get(file_path_key) if file_path_key else None
        if not pdb_file:
            # Fallback to any available file (shouldn't happen, but be safe)
            from proprep.utils.structure_selector import get_priority_pdb_file
            pdb_file = get_priority_pdb_file(workspace, self.console, silent=True) or "structure.pdb"

        # Model selection for multi-model structures (e.g., NMR)
        selected_model_idx = self._get_model_selection(structure)
        selected_model = structure[selected_model_idx]

        # Initialize detector
        self.detector = ComprehensiveRedoxDetector(console=self.console, processor=self.processor)
        self.detector.source_pdb_file = pdb_file

        # Run detection on selected model
        self.console.print("[cyan]Running comprehensive redox site detection...[/cyan]")

        # Create a temporary single-model structure for detection
        from Bio.PDB import Structure as BioStructure
        temp_structure = BioStructure.Structure("temp")
        temp_structure.add(selected_model)

        detected_sites = self.detector.detect_redox_sites(
            structure=temp_structure,
            selected_chains=None,  # Analyze all chains
            interactive=True
        )

        if detected_sites:
            self.detected_sites = detected_sites
            workspace.set("detected_redox_sites", detected_sites)

            # Save transformer mappings to workspace if they were created during detection
            if hasattr(self.detector, 'transformer_mappings') and self.detector.transformer_mappings:
                workspace.set("redox_transformer_mappings", self.detector.transformer_mappings)
                self.console.print(f"[grey50]Saved transformer mappings for {len(self.detector.transformer_mappings)} site type(s)[/grey50]")

            num_sites = len(detected_sites)
            self.console.print(f"[green]✓ Successfully detected {num_sites} redox-active site(s)[/green]")

            # Summary is available via Site Review Options menu (option 1)
            # self._show_sites_summary(detected_sites)

            return True
        else:
            self.console.print("[yellow]No redox sites detected[/yellow]")
            return False

    def _get_model_selection(self, structure) -> int:
        """Prompt user to select a model from multi-model structure"""
        num_models = len(structure)

        if num_models == 1:
            return 0  # Only one model, use it

        self.console.print(f"\n[bold]Structure contains {num_models} models (NMR ensemble)[/bold]")
        self.console.print("Please select which model to use for redox site detection:\n")

        for i, model in enumerate(structure):
            self.console.print(f"  {i+1}. Model {model.id}")

        if self.processor:
            choice = prompt_with_context(
                processor=self.processor,
                prompt="\nSelect model",
                choices=[str(i+1) for i in range(num_models)],
                default="1",
                module="Redox Site Detector",
                description="Select PDB model for detection",
                options_map={str(i+1): f"Model {i}" for i in range(num_models)}
            )
        else:
            choice = prompt_with_context(None,
                "\nSelect model",
                choices=[str(i+1) for i in range(num_models)],
                default="1"
            )

        return int(choice) - 1

    def _view_detected_sites(self) -> bool:
        """View previously detected redox sites"""
        workspace = self.processor.workspace if self.processor else {}
        detected_sites = workspace.get("detected_redox_sites")

        if not detected_sites:
            self.console.print("[yellow]No redox sites have been detected yet[/yellow]")
            self.console.print("Run option 1 first to detect redox sites")
            return False

        self._show_sites_summary(detected_sites)
        return True

    def _show_sites_summary(self, sites):
        """Display summary of detected redox sites"""
        if not sites:
            self.console.print("[yellow]No sites to display[/yellow]")
            return

        # Use the proper detailed display from SiteRefinementInterface
        self._show_detailed_sites(sites)

    def _show_detailed_sites(self, sites):
        """Show detailed information for each site"""
        try:
            from .comprehensive_redox_detector import SiteRefinementInterface, DetectionConfig

            # Create refinement interface for display
            # Even if we don't have a detector, we can create a minimal config
            # since _display_site_summary doesn't actually use the config
            if self.detector:
                config = self.detector.config
            else:
                # Create minimal config with defaults for display purposes
                config = DetectionConfig()

            refinement_interface = SiteRefinementInterface(
                config,
                console=self.console
            )

            for site in sites:
                self.console.print(f"\n[bold underline]Site Details: {site.site_id}[/bold underline]")

                # Display site type if available
                if hasattr(site, 'site_type') and site.site_type:
                    self.console.print(f"[bold]Site Type:[/bold] {site.site_type}\n")

                refinement_interface._display_site_summary(site)

            # Prompt once at the end
            input("\nPress Enter to continue...")

        except Exception as e:
            self.console.print(f"[yellow]Could not display detailed view: {e}[/yellow]")
            logger.error(f"Error displaying sites: {e}", exc_info=True)

    def _remove_redox_site(self) -> bool:
        """Remove one or more redox sites from the detected sites list"""
        from rich.table import Table
        from proprep.utils.prompts import confirm_with_context

        workspace = self.processor.workspace if self.processor else {}
        detected_sites = workspace.get("detected_redox_sites")

        if not detected_sites:
            self.console.print("[yellow]No redox sites have been detected yet[/yellow]")
            self.console.print("Run option 2 to detect redox sites first")
            return False

        if len(detected_sites) == 0:
            self.console.print("[yellow]No sites to remove[/yellow]")
            return False

        self.console.print("\n[bold cyan]Remove Redox Site(s)[/bold cyan]\n")

        # Display current sites in a table
        table = Table(title="Currently Detected Redox Sites")
        table.add_column("Index", style="cyan", width=6)
        table.add_column("Site ID", style="green")
        table.add_column("Centers", style="yellow", justify="right")
        table.add_column("Atoms", style="yellow", justify="right")
        table.add_column("Site Type", style="magenta")

        for i, site in enumerate(detected_sites, 1):
            num_centers = len(site.centers) if hasattr(site, 'centers') else 0
            num_atoms = len(site.atoms) if hasattr(site, 'atoms') else 0
            site_type = site.site_type if hasattr(site, 'site_type') and site.site_type else "N/A"
            table.add_row(str(i), site.site_id, str(num_centers), str(num_atoms), site_type)

        self.console.print(table)

        # Get user selection
        self.console.print("\n[grey50]Enter site index/indices to remove (e.g., '1', '1,3', '1-3', or 'all')[/grey50]")
        self.console.print("[grey50]Enter 'cancel' to go back[/grey50]")

        choice = prompt_with_context(
            self.processor,
            "\nSite(s) to remove",
            module="Redox Site Detector",
            description="Select site(s) to remove from detected sites"
        )

        if choice.lower() == 'cancel':
            self.console.print("[grey50]Remove cancelled[/grey50]")
            return False

        # Parse the selection
        try:
            indices_to_remove = self._parse_selection(choice, len(detected_sites))
        except ValueError as e:
            self.console.print(f"[red]Invalid selection: {e}[/red]")
            return False

        if not indices_to_remove:
            self.console.print("[yellow]No sites selected for removal[/yellow]")
            return False

        # Show what will be removed
        self.console.print("\n[bold]Sites to be removed:[/bold]")
        for idx in sorted(indices_to_remove):
            site = detected_sites[idx - 1]
            num_centers = len(site.centers) if hasattr(site, 'centers') else 0
            self.console.print(f"  {idx}. {site.site_id} ({num_centers} center(s))")

        # Confirm removal
        if self.processor:
            from .pdb_filter import confirm_with_context
            confirm = confirm_with_context(
                processor=self.processor,
                prompt="\nConfirm removal of these site(s)?",
                default=False,
                module="Redox Site Detector",
                description="Confirm site removal"
            )
        else:
            confirm = confirm_with_context(None, "\nConfirm removal of these site(s)?", default=False)

        if not confirm:
            self.console.print("[grey50]Removal cancelled[/grey50]")
            return False

        # Remove the sites (in reverse order to maintain indices)
        for idx in sorted(indices_to_remove, reverse=True):
            removed_site = detected_sites.pop(idx - 1)
            self.console.print(f"[green]✓ Removed {removed_site.site_id}[/green]")

        # Update workspace
        workspace.set("detected_redox_sites", detected_sites)
        self.detected_sites = detected_sites

        remaining = len(detected_sites)
        self.console.print(f"\n[green]✓ Successfully removed {len(indices_to_remove)} site(s)[/green]")
        self.console.print(f"[cyan]Remaining sites: {remaining}[/cyan]")

        return True

    def _edit_redox_sites(self) -> bool:
        """Launch the comprehensive RedoxSite editor"""
        from .redox_site_editor import RedoxSiteEditor

        workspace = self.processor.workspace if self.processor else {}
        detected_sites = workspace.get("detected_redox_sites")

        if not detected_sites:
            self.console.print("[yellow]No redox sites have been detected yet[/yellow]")
            self.console.print("Run option 2 to detect redox sites first")
            return False

        if len(detected_sites) == 0:
            self.console.print("[yellow]No sites to edit[/yellow]")
            return False

        # Launch editor with workspace for structure loading
        editor = RedoxSiteEditor(
            sites=detected_sites,
            console=self.console,
            processor=self.processor,
            workspace=workspace
        )

        modified_sites = editor.run()

        if modified_sites is not None:
            # User saved changes
            workspace.set("detected_redox_sites", modified_sites)
            self.detected_sites = modified_sites
            self.console.print("\n[green]✓ Changes saved to workspace[/green]")
            return True
        else:
            # User cancelled
            self.console.print("\n[grey50]No changes made[/grey50]")
            return False

    def _parse_selection(self, selection: str, max_index: int) -> set:
        """
        Parse user selection string into a set of indices.

        Supports:
        - Single index: '1'
        - Multiple indices: '1,3,5'
        - Ranges: '1-3' (inclusive)
        - Combined: '1,3-5,7'
        - All: 'all'

        Returns:
            Set of 1-based indices

        Raises:
            ValueError: If selection is invalid
        """
        selection = selection.strip().lower()

        if selection == 'all':
            return set(range(1, max_index + 1))

        indices = set()
        parts = selection.split(',')

        for part in parts:
            part = part.strip()

            if '-' in part:
                # Range
                try:
                    start, end = part.split('-')
                    start_idx = int(start.strip())
                    end_idx = int(end.strip())

                    if start_idx < 1 or end_idx > max_index:
                        raise ValueError(f"Range {start_idx}-{end_idx} is out of bounds (1-{max_index})")
                    if start_idx > end_idx:
                        raise ValueError(f"Invalid range: {start_idx}-{end_idx}")

                    indices.update(range(start_idx, end_idx + 1))
                except ValueError as e:
                    if "invalid literal" in str(e):
                        raise ValueError(f"Invalid range format: '{part}'")
                    raise
            else:
                # Single index
                try:
                    idx = int(part)
                    if idx < 1 or idx > max_index:
                        raise ValueError(f"Index {idx} is out of bounds (1-{max_index})")
                    indices.add(idx)
                except ValueError as e:
                    if "invalid literal" in str(e):
                        raise ValueError(f"Invalid index: '{part}'")
                    raise

        return indices
