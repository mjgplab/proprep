"""
Interactive Structure Viewer

Web-based structure viewer with workspace-aware annotations.
Launches a local web server displaying structures from workspace with
optional highlighting of redox sites, titratable residues, and other
workspace-stored annotations.

Features:
- Multi-structure viewing via StructureSelector
- Redox site visualization from detected_redox_sites
- Titratable residue highlighting from protonation_results
- Filter context visualization (kept vs removed portions)
- Bond classification display (disulfide, coordinate, metal-metal, covalent)
- Pre-aligned structure comparison

Author: ProPrep Development Team
Date: 2025-11-14
Version: 3.0.0
"""

import atexit
import logging
import os
from typing import Dict, List, Optional, Any
from pathlib import Path

from rich.console import Console
from proprep.utils.prompts import prompt_with_context, confirm_with_context
from rich.panel import Panel
from rich.table import Table

from proprep.utils.module_registry import ProcessingModule, register_module
from proprep.utils.structure_selector import StructureSelector

logger = logging.getLogger(__name__)


# Process-wide handle to the most-recently-started ViewerServer. The 7
# launch sites all funnel through ``_launch_viewer`` below; this avoids
# leaking ports across repeated launches and lets the web shell rely on a
# stable single viewer per ProPrep run. (Subsequent launches stop the
# previous server before starting a new one — see ``_launch_viewer``.)
_active_viewer_server = None


def _stop_active_viewer_server() -> None:
    """Best-effort stop of the process-wide viewer server. Idempotent."""
    global _active_viewer_server
    server = _active_viewer_server
    if server is None:
        return
    try:
        if hasattr(server, "is_running") and server.is_running():
            server.stop()
    except Exception as exc:  # pragma: no cover - defensive
        logger.debug("Error stopping viewer server during shutdown: %s", exc)
    _active_viewer_server = None


# Run on interpreter exit so we don't leak the bound port if ProPrep
# crashes or is Ctrl-C'd before normal cleanup.
atexit.register(_stop_active_viewer_server)


@register_module
class InteractiveStructureViewer(ProcessingModule):
    """
    Interactive web-based structure viewer with workspace annotations.

    Replaces the legacy VMD Visualizer with a modern web-based approach
    using NGL Viewer. Integrates seamlessly with ProPrep workspace to
    display redox sites, titratable residues, and other annotations.
    """

    NAME = "Structure Viewer"
    DESCRIPTION = "Interactive web-based viewer with workspace annotations"
    CATEGORY = "visualization"
    VERSION = "3.0.0"

    def __init__(self):
        """Initialize the Interactive Structure Viewer."""
        super().__init__()
        self.selected_structures = []
        self.available_annotations = {}
        self.annotation_config = {}
        self.viewer_config = {}
        self.shape_config = {}

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

    @property
    def workspace(self):
        """Get workspace from processor."""
        if hasattr(self, 'processor') and self.processor:
            return self.processor._get_workspace()
        return {}

    def get_workspace_requirements(self) -> List[str]:
        """
        Get workspace requirements.

        Module requires at least one structure loaded by Structure Loader.
        All annotation keys are optional.

        Returns:
            List with structure requirement (OR syntax for any structure source)
        """
        # Require at least one structure from Structure Loader
        return ["rcsb_pdb_file | local_pdb_file | alphafold_pdb_file | alphafill_pdb_file"]

    def availability_note(self, workspace):
        """Menu note when unavailable (○). Mirrors can_process exactly."""
        return None if self.can_process(workspace) else "Needs a loaded structure"

    def can_process(self, workspace) -> bool:
        """
        Check if module can process current workspace.

        Module is always available as a utility if at least one structure exists.

        Args:
            workspace: ProPrep workspace object

        Returns:
            True if at least one structure exists in workspace
        """
        selector = StructureSelector(workspace, self.console)
        available = selector.get_available_structures()

        has_structures = len(available) > 0

        if not has_structures:
            logger.debug("No structures found in workspace - viewer unavailable")

        return has_structures

    def get_menu_options(self) -> Dict[str, str]:
        """
        Get available menu options.

        Returns:
            Dict of option keys and descriptions
        """
        return {
            "launch": "Launch interactive viewer",
            "info": "Show available annotations",
        }

    def handle_menu_option(self, option: str) -> bool:
        """
        Handle menu option selection using command pattern.

        Args:
            option: Selected menu option key

        Returns:
            True if option handled successfully
        """
        if option == "launch":
            from proprep.structure_prep.viewer_commands import LaunchViewerCommand
            command = LaunchViewerCommand(self.processor)
            return command.execute()
        elif option == "info":
            from proprep.structure_prep.viewer_commands import ShowAnnotationInfoCommand
            command = ShowAnnotationInfoCommand(self.processor)
            return command.execute()

        return False

    def quick_launch(self, pdb_files: List[str]) -> bool:
        """
        Launch viewer directly with specified files, bypassing interactive selection.

        This method is used by the --pdbview shortcut to quickly view a structure
        without going through the full workflow prompts.

        Args:
            pdb_files: List of PDB file paths to display

        Returns:
            True if viewer launched successfully
        """
        self.selected_structures = pdb_files
        self.available_annotations = {}  # Skip annotation detection for quick view
        self.annotation_config = {}
        return self._launch_viewer()

    # ========================================================================
    # Structure Selection (Phase 1)
    # ========================================================================

    def _select_structures_for_viewing(self) -> List[str]:
        """
        Let user interactively select one or more structures from workspace.

        Uses StructureSelector._interactive_multi_selection() to allow
        single or multiple structure selection with various input formats:
        - "1" - single selection
        - "1,3,5" - comma-separated
        - "1-3" - range notation
        - "all" - select all structures

        Returns:
            List of PDB file paths selected by user
        """
        selector = StructureSelector(
            workspace=self.workspace,
            console=self.console,
            processor=self.processor  # For session recording
        )

        # Get all available structures
        available = selector.get_available_structures()

        if not available:
            self.console.print("[red]No structures found in workspace[/red]")
            self.console.print("[yellow]Please load a structure first[/yellow]")
            return []

        self.console.print(Panel(
            "[bold cyan]Structure Selection[/bold cyan]\n\n"
            "Select one or more structures to view in the interactive viewer.\n"
            "Structures will be displayed together (superimposed if aligned).",
            title="Interactive Structure Viewer",
            border_style="cyan",
            expand=False
        ))

        # Use the interactive multi-selection method
        selected_files = selector._interactive_multi_selection(available)

        self.selected_structures = selected_files

        return selected_files

    # ========================================================================
    # Annotation Detection (Phase 2)
    # ========================================================================

    def _detect_available_annotations(self) -> Dict[str, Any]:
        """
        Scan workspace for annotation data that can be highlighted.

        Checks for:
        - detected_redox_sites (List[RedoxSite])
        - protonation_results (Dict[str, Dict])
        - filter_selections (Dict[str, Dict[str, List[int]]])
        - aligned_structures (Dict[int, Structure])

        Returns:
            Dict of annotation types and their data with counts and descriptions
        """
        annotations = {}

        # 1. Redox Sites (HIGH PRIORITY)
        detected_redox_sites = self.workspace.get('detected_redox_sites')
        if detected_redox_sites and len(detected_redox_sites) > 0:
            # Option 1: All redox sites combined (existing behavior)
            annotations['redox_sites'] = {
                'count': len(detected_redox_sites),
                'data': detected_redox_sites,
                'description': f"{len(detected_redox_sites)} redox site(s) - combined view",
                'priority': 'HIGH'
            }
            logger.debug(f"Found {len(detected_redox_sites)} redox sites in workspace")

            # Option 2: Individual redox sites (NEW) - only for single structure
            if len(self.selected_structures) == 1:
                annotations['individual_redox_sites'] = {
                    'count': len(detected_redox_sites),
                    'data': detected_redox_sites,
                    'description': f"{len(detected_redox_sites)} redox site(s) - toggle each individually",
                    'priority': 'HIGH',
                    'requires_selection': len(detected_redox_sites) > 10,
                    'single_structure_only': True
                }
                logger.debug(f"Individual redox site viewing available ({len(detected_redox_sites)} sites)")

        # 2. Titratable Residues from Protonation Analysis (HIGH PRIORITY)
        protonation_results = self.workspace.get('protonation_results')
        if protonation_results and len(protonation_results) > 0:
            # Filter out termini if needed for cleaner visualization
            titratable_residues = {
                k: v for k, v in protonation_results.items()
                if v.get('type') not in ['N_TERM', 'C_TERM']
            }

            if titratable_residues:
                annotations['titratable_residues'] = {
                    'count': len(titratable_residues),
                    'data': titratable_residues,
                    'description': f"{len(titratable_residues)} titratable residue(s) (excludes termini)",
                    'priority': 'HIGH',
                    'all_data': protonation_results  # Include termini if user wants them
                }
                logger.info(f"Found {len(titratable_residues)} titratable residues in workspace")

        # 3. Filter Selections - Context Visualization (HIGH PRIORITY)
        filter_selections = self.workspace.get('filter_selections')
        if filter_selections and len(filter_selections) > 0:
            # Count total kept residues across all chains
            total_kept = sum(
                len(residues)
                for chain_data in filter_selections.values()
                for residues in chain_data.values()
                if residues  # Non-empty list
            )

            annotations['filter_selections'] = {
                'count': len(filter_selections),  # Number of chains
                'data': filter_selections,
                'description': f"{len(filter_selections)} chain(s) filtered, {total_kept} residues kept",
                'priority': 'HIGH'
            }
            logger.debug(f"Found filter selections for {len(filter_selections)} chains in workspace")

        # 4. Bond Classifications from Redox Sites (MEDIUM PRIORITY)
        if detected_redox_sites:
            bonds_by_type = self._extract_bonds_by_type(detected_redox_sites)
            if bonds_by_type:
                total_bonds = sum(len(bonds) for bonds in bonds_by_type.values())
                annotations['bond_classifications'] = {
                    'count': total_bonds,
                    'data': bonds_by_type,
                    'description': f"{total_bonds} classified bond(s) (disulfide, coordinate, metal-metal, covalent)",
                    'priority': 'MEDIUM'
                }
                logger.debug(f"Extracted {total_bonds} bonds from redox sites")

        # Note: Aligned structures annotation removed - users already select which structures to view
        # No need to annotate aligned structures separately

        return annotations

    def _extract_bonds_by_type(self, redox_sites) -> Dict[str, List]:
        """
        Extract bonds from RedoxSite objects grouped by chemical type.

        Args:
            redox_sites: List of RedoxSite objects

        Returns:
            Dict mapping chemical_type to list of bond info dicts
        """
        bonds_by_type = {
            'disulfide': [],
            'coordinate': [],
            'metal-metal': [],
            'covalent': []
        }

        for site in redox_sites:
            for bond in site.bonds:
                chemical_type = bond.chemical_type

                if chemical_type in bonds_by_type:
                    bond_info = {
                        'atom1_residue': bond.atom1_residue_info,
                        'atom2_residue': bond.atom2_residue_info,
                        'distance': bond.distance,
                        'bond_type': bond.bond_type,  # intraresidue or interresidue
                        'site_id': site.site_id
                    }
                    bonds_by_type[chemical_type].append(bond_info)

        # Remove empty categories
        bonds_by_type = {k: v for k, v in bonds_by_type.items() if v}

        return bonds_by_type

    # ========================================================================
    # Configuration & Launch (Phase 3-5 - Placeholder)
    # ========================================================================

    def _configure_annotation_display(self) -> Dict[str, Any]:
        """
        Interactive menu for configuring annotation display.

        Allows user to select which annotations to display and configure
        their visual appearance (color, representation style).

        Returns:
            Configuration dict for viewer with annotation settings
        """
        if not self.available_annotations:
            return {}

        config = {}

        self.console.print(Panel(
            "[bold cyan]Annotation Configuration[/bold cyan]\n\n"
            "Configure which annotations to display and their visual style.",
            title="Interactive Structure Viewer",
            border_style="cyan"
        ))

        # Display available annotations
        table = Table(title="Available Annotations")
        table.add_column("#", style="cyan", justify="center")
        table.add_column("Type", style="green")
        table.add_column("Count", style="yellow", justify="right")
        table.add_column("Description", style="white")

        ann_list = list(self.available_annotations.items())
        for idx, (ann_type, ann_data) in enumerate(ann_list, 1):
            table.add_row(
                str(idx),
                ann_type.replace('_', ' ').title(),
                str(ann_data['count']),
                ann_data['description']
            )

        self.console.print(table)

        # Ask which annotations to display
        self.console.print("\n[cyan]Select annotations to display:[/cyan]")
        self.console.print("[grey50]Enter numbers (e.g., 1,2,3), ranges (e.g., 1-3), or 'all'[/grey50]")
        self.console.print("[grey50]Press Enter for 'all'[/grey50]")

        from proprep.application.menu_commands import prompt_with_context
        selection = prompt_with_context(
            processor=self.processor,
            prompt="Select annotations",
            default="all",
            module="Structure Viewer",
            description="Select annotations to display"
        )

        # Parse selection
        selected_indices = self._parse_selection(selection, len(ann_list))

        if not selected_indices:
            self.console.print("[yellow]No annotations selected[/yellow]")
            return {}

        # Configure each selected annotation
        for idx in selected_indices:
            ann_type, ann_data = ann_list[idx - 1]
            ann_config = self._configure_single_annotation(ann_type, ann_data)
            if ann_config:
                # Special handling for individual_redox_sites which returns multiple configs
                if ann_type == 'individual_redox_sites':
                    # Merge individual site configs into main config
                    config.update(ann_config)
                else:
                    config[ann_type] = ann_config

        return config

    def _parse_selection(self, selection: str, max_num: int) -> List[int]:
        """
        Parse user selection string into list of indices.

        Args:
            selection: User input (e.g., "1,2,3", "1-3", "all")
            max_num: Maximum valid number

        Returns:
            List of selected indices (1-based)
        """
        if selection.lower() == "all":
            return list(range(1, max_num + 1))

        indices = set()
        parts = selection.split(",")

        for part in parts:
            part = part.strip()
            if "-" in part:
                # Range notation
                try:
                    start, end = part.split("-", 1)
                    start_idx = int(start.strip())
                    end_idx = int(end.strip())
                    for idx in range(start_idx, end_idx + 1):
                        if 1 <= idx <= max_num:
                            indices.add(idx)
                except ValueError:
                    continue
            else:
                # Single number
                try:
                    idx = int(part)
                    if 1 <= idx <= max_num:
                        indices.add(idx)
                except ValueError:
                    continue

        return sorted(list(indices))

    def _configure_single_annotation(self, ann_type: str, ann_data: Dict) -> Dict[str, Any]:
        """
        Configure display settings for a single annotation type.

        Args:
            ann_type: Annotation type key
            ann_data: Annotation data from workspace

        Returns:
            Configuration dict with color, style, and NGL selection
        """
        from proprep.application.menu_commands import prompt_with_context

        # Special handling for individual_redox_sites
        if ann_type == 'individual_redox_sites':
            return self._configure_individual_redox_sites(ann_data)

        self.console.print(f"\n[bold cyan]Configure: {ann_type.replace('_', ' ').title()}[/bold cyan]")

        # Get default color and style based on annotation type
        default_color, default_style = self._get_default_style(ann_type)

        # Structure selection prompt if multiple structures loaded
        structure_indices = []
        if len(self.selected_structures) > 1:
            structure_options = {}
            for i, struct_file in enumerate(self.selected_structures, 1):
                struct_name = Path(struct_file).stem
                structure_options[str(i)] = struct_name
            structure_options[str(len(self.selected_structures) + 1)] = "All structures"

            # Display structure options
            self.console.print("\n[bold]Target Structure(s):[/bold]")
            for key, value in structure_options.items():
                self.console.print(f"  {key}. {value}")

            structure_choice = prompt_with_context(
                processor=self.processor,
                prompt="Apply annotation to which structure(s)?",
                choices=list(structure_options.keys()),
                default=str(len(self.selected_structures) + 1),  # Default to "All"
                module="Structure Viewer",
                description=f"Target structure(s) for {ann_type}",
                options_map=structure_options
            )

            if structure_choice == str(len(self.selected_structures) + 1):
                # All structures
                structure_indices = list(range(len(self.selected_structures)))
            else:
                # Specific structure (convert to 0-indexed)
                structure_indices = [int(structure_choice) - 1]
        else:
            # Only one structure, apply to it
            structure_indices = [0]

        # Color selection - per structure if multiple structures selected
        color_options = {
            "1": "Red",
            "2": "Blue",
            "3": "Green",
            "4": "Yellow",
            "5": "Orange",
            "6": "Cyan",
            "7": "Magenta",
            "8": "White"
        }

        color_styles = {
            "1": "red",
            "2": "blue",
            "3": "green",
            "4": "yellow",
            "5": "#ff8c00",  # orange
            "6": "cyan",
            "7": "magenta",
            "8": "white"
        }

        color_map = {
            "1": "#ff0000", "2": "#0000ff", "3": "#00ff00", "4": "#ffff00",
            "5": "#ff8c00", "6": "#00ffff", "7": "#ff00ff", "8": "#ffffff"
        }

        # If applying to multiple structures, ask for color per structure
        colors_per_structure = {}
        if len(structure_indices) > 1:
            self.console.print("\n[bold cyan]Configure colors for each structure:[/bold cyan]")

            for idx in structure_indices:
                struct_name = Path(self.selected_structures[idx]).stem

                # Display color options with actual colors
                self.console.print(f"\n[bold]Colors for {struct_name}:[/bold]")
                for key, name in color_options.items():
                    style = color_styles[key]
                    self.console.print(f"  {key}. [{style}]■■■[/{style}] {name}")

                color_choice = prompt_with_context(
                    processor=self.processor,
                    prompt=f"Select color for {struct_name}",
                    choices=list(color_options.keys()),
                    default=str(idx + 1) if str(idx + 1) in color_options else "1",
                    module="Structure Viewer",
                    description=f"Color for {ann_type} on {struct_name}",
                    options_map=color_options
                )

                colors_per_structure[idx] = color_map[color_choice]
        else:
            # Single structure - one color
            # Display color options with actual colors
            self.console.print("\n[bold]Available Colors:[/bold]")
            for key, name in color_options.items():
                style = color_styles[key]
                self.console.print(f"  {key}. [{style}]■■■[/{style}] {name}")

            color_choice = prompt_with_context(
                processor=self.processor,
                prompt="Select color",
                choices=list(color_options.keys()),
                default="1",
                module="Structure Viewer",
                description=f"Color for {ann_type}",
                options_map=color_options
            )

            colors_per_structure[structure_indices[0]] = color_map[color_choice]

        # Representation style selection
        style_options = {
            "1": "Licorice",
            "2": "Ball+Stick",
            "3": "Spacefill",
            "4": "Cartoon",
            "5": "Surface"
        }

        # Display style options
        self.console.print("\n[bold]Representation Styles:[/bold]")
        for key, value in style_options.items():
            self.console.print(f"  {key}. {value}")

        style_choice = prompt_with_context(
            processor=self.processor,
            prompt="Select representation style",
            choices=list(style_options.keys()),
            default="1",
            module="Structure Viewer",
            description=f"Style for {ann_type}",
            options_map=style_options
        )

        style_map = {
            "1": "licorice", "2": "ball+stick", "3": "spacefill",
            "4": "cartoon", "5": "surface"
        }
        style = style_map[style_choice]

        # Generate NGL selection based on annotation type
        ngl_selection = self._generate_ngl_selection(ann_type, ann_data)

        return {
            'colors_per_structure': colors_per_structure,  # Dict mapping structure index to color
            'style': style,
            'ngl_selection': ngl_selection,
            'count': ann_data['count'],
            'priority': ann_data.get('priority', 'MEDIUM'),
            'structure_indices': structure_indices  # Which structures to apply annotation to
        }

    def _configure_individual_redox_sites(self, ann_data: Dict) -> Dict[str, Any]:
        """
        Configure individual redox site visualizations with per-site toggle controls.

        Args:
            ann_data: Annotation data containing detected_redox_sites list

        Returns:
            Dict mapping site IDs to individual annotation configs
        """
        from proprep.application.menu_commands import prompt_with_context

        redox_sites = ann_data['data']
        total_sites = len(redox_sites)

        self.console.print(f"\n[bold cyan]Configure Individual Redox Sites[/bold cyan]")
        self.console.print(f"[grey50]Detected {total_sites} redox site(s)[/grey50]\n")

        # If more than 10 sites, prompt user to select which ones to view
        selected_sites = redox_sites
        if total_sites > 10:
            self.console.print("[yellow]Warning: More than 10 redox sites detected.[/yellow]")
            self.console.print("[grey50]Please specify which sites you want to visualize.[/grey50]\n")

            # Display all available sites
            self.console.print("[bold]Available Sites:[/bold]")
            for site in redox_sites:
                site_id = site.site_id.replace('site_', '')
                self.console.print(f"  Site {site_id}")

            self.console.print("\n[bold]Selection Options:[/bold]")
            self.console.print("  • Enter 'all' to view all sites")
            self.console.print("  • Enter specific site numbers (e.g., '1,2,5')")
            self.console.print("  • Enter ranges (e.g., '1-5,8,10-12')")

            selection = prompt_with_context(
                getattr(self, 'processor', None),
                "\n[cyan]Which sites would you like to view?[/cyan]",
                default="all",
                module="Structure Viewer",
                description="Select redox sites to view (all, list, or ranges)",
            )

            if selection.lower() != 'all':
                # Parse selection
                selected_indices = self._parse_site_selection(selection, total_sites)
                if selected_indices:
                    selected_sites = [redox_sites[i] for i in selected_indices]
                    self.console.print(f"[green]Selected {len(selected_sites)} site(s)[/green]\n")
                else:
                    self.console.print("[yellow]Invalid selection, showing all sites[/yellow]\n")

        # Generate individual annotation configs for each selected site
        site_configs = {}

        for site in selected_sites:
            site_id = site.site_id.replace('site_', '')
            site_key = f"redox_site_{site_id}"

            # Generate NGL selection for this site
            ngl_selection = self._redox_site_to_ngl_selection(site)

            if ngl_selection:
                site_configs[site_key] = {
                    'label': f"Site {site_id}",
                    'color': "element",  # Element coloring for ball+stick
                    'style': "ball+stick",
                    'ngl_selection': ngl_selection,
                    'site_id': site_id,
                    'structure_indices': [0]  # Single structure only
                }

        self.console.print(f"[green]Configured {len(site_configs)} redox site(s) for visualization[/green]")

        return site_configs

    def _parse_site_selection(self, selection: str, max_sites: int) -> List[int]:
        """
        Parse user's site selection string into list of 0-indexed site indices.

        Args:
            selection: User input string (e.g., "1,2,5" or "1-5,8")
            max_sites: Total number of sites available

        Returns:
            List of 0-indexed site indices
        """
        indices = set()

        # Split by comma
        parts = selection.split(',')

        for part in parts:
            part = part.strip()

            # Check for range (e.g., "1-5")
            if '-' in part:
                try:
                    start, end = part.split('-')
                    start_idx = int(start.strip()) - 1  # Convert to 0-indexed
                    end_idx = int(end.strip()) - 1

                    # Validate range
                    if 0 <= start_idx < max_sites and 0 <= end_idx < max_sites:
                        indices.update(range(start_idx, end_idx + 1))
                    else:
                        self.console.print(f"[yellow]Warning: Range '{part}' out of bounds, skipping[/yellow]")
                except ValueError:
                    self.console.print(f"[yellow]Warning: Invalid range '{part}', skipping[/yellow]")
            else:
                # Single number
                try:
                    idx = int(part) - 1  # Convert to 0-indexed
                    if 0 <= idx < max_sites:
                        indices.add(idx)
                    else:
                        self.console.print(f"[yellow]Warning: Site {part} out of bounds, skipping[/yellow]")
                except ValueError:
                    self.console.print(f"[yellow]Warning: Invalid site number '{part}', skipping[/yellow]")

        return sorted(list(indices))

    def _get_default_style(self, ann_type: str) -> tuple:
        """
        Get default color and style for annotation type.

        Args:
            ann_type: Annotation type key

        Returns:
            Tuple of (color_hex, style_name)
        """
        defaults = {
            'redox_sites': ("#ff0000", "ball+stick"),
            'titratable_residues': ("#ffff00", "licorice"),
            'filter_selections': ("#808080", "cartoon"),
            'bond_classifications': ("#00ffff", "licorice")
        }
        return defaults.get(ann_type, ("#ffffff", "licorice"))

    def _generate_ngl_selection(self, ann_type: str, ann_data: Dict) -> str:
        """
        Generate NGL selection string for annotation type.

        Args:
            ann_type: Annotation type key
            ann_data: Annotation data from workspace

        Returns:
            NGL selection string
        """
        if ann_type == 'redox_sites':
            # Combine all redox sites
            selections = []
            for site in ann_data['data']:
                site_sel = self._redox_site_to_ngl_selection(site)
                if site_sel:
                    selections.append(f"({site_sel})")
            return " or ".join(selections) if selections else ""

        elif ann_type == 'titratable_residues':
            return self._protonation_results_to_ngl_selection(ann_data['data'])

        elif ann_type == 'filter_selections':
            ngl_sels = self._filter_selections_to_ngl(ann_data['data'])
            return ngl_sels.get('kept', '')

        elif ann_type == 'bond_classifications':
            # Combine all bond types
            selections = []
            for bond_type, bonds in ann_data['data'].items():
                bond_sel = self._bonds_to_ngl_selection(bonds)
                if bond_sel:
                    selections.append(f"({bond_sel})")
            return " or ".join(selections) if selections else ""

        return ""

    def _launch_viewer_workflow(self) -> bool:
        """
        Main workflow: Select structures, configure annotations, launch viewer.

        This is the entry point called by LaunchViewerCommand.

        Returns:
            True if workflow completed successfully
        """
        # Step 1: Select structures
        selected_files = self._select_structures_for_viewing()

        if not selected_files:
            self.console.print("[yellow]No structures selected - viewer not launched[/yellow]")
            return False

        # Step 2: Detect available annotations (Phase 2)
        self.available_annotations = self._detect_available_annotations()

        # Step 3: Configure annotations if any found (Phase 4)
        if self.available_annotations:
            self.annotation_config = self._configure_annotation_display()

        # Step 4: Launch viewer (Phase 5)
        return self._launch_viewer()

    def _launch_viewer(self, open_browser: bool = True) -> bool:
        """
        Launch the web-based viewer with HTTP server.

        Creates configuration JSON, starts HTTP server, and (in CLI mode)
        opens a browser tab. Stops any previously-launched server first so
        the port (8765) is reused — the 7 viewer launch sites in ProPrep
        accumulate over a run, and each one creating a fresh server would
        leak ports.

        ``open_browser=False`` starts/relaunches the server without popping
        a browser tab. The coordinator passes this when it is *relaunching*
        an already-running viewer purely to swap the bound structure file
        (the server binds its file list at construction): in CLI the
        existing tab re-sources itself from the reused port via app.js's
        cache-buster, so a fresh tab would be a redundant, unbidden window.
        A genuinely fresh launch — or an explicit user "open the viewer"
        request — leaves this True. (In web-shell mode ``ViewerServer.start``
        suppresses the tab regardless, sourcing the iframe instead.)

        Returns:
            True if viewer launched successfully
        """
        # Headless batch replay (PROPREP_BATCH): never spin up a real viewer
        # server or browser tab. The batch runner may set PROPREP_WEB_SHELL to
        # match a web-recorded template so mode-gated PROMPTS stay in sync, but
        # a background batch has no browser to serve — launching one server per
        # protein would leak threads/ports. Prompt suppression still happens
        # upstream; this only stubs out the actual launch.
        if os.environ.get("PROPREP_BATCH"):
            return False

        from proprep.structure_prep import interactive_structure_viewer as _module
        from proprep.structure_prep.viewer_server import ViewerServer

        # Tear down any prior viewer server before claiming the port.
        _stop_active_viewer_server()

        # Build complete configuration for viewer
        viewer_config = self._build_viewer_config()

        # Create server instance
        try:
            server = ViewerServer(
                config=viewer_config,
                structure_files=self.selected_structures,
                port=8765
            )

            # Start server. ViewerServer.start() additionally consults
            # PROPREP_WEB_SHELL to suppress the tab in web-shell mode; in
            # CLI it honours open_browser as passed.
            success = server.start(open_browser=open_browser)

            if not success:
                self.console.print("[red]Failed to start viewer server[/red]")
                return False

            # Display success message
            url = server.get_url()
            n = len(self.selected_structures)
            self.console.print(
                f"[green]✓ Structure viewer launched at {url} ({n} structure{'s' if n != 1 else ''})[/green]"
            )

            # Keep reference to server so it doesn't get garbage collected
            self.server = server
            _module._active_viewer_server = server

            return True

        except Exception as e:
            logger.error(f"Error launching viewer: {e}")
            self.console.print(f"[red]Error launching viewer: {e}[/red]")
            return False

    def update_annotations(
        self,
        annotation_config: Optional[Dict[str, Any]] = None,
        viewer_config: Optional[Dict[str, Any]] = None,
    ) -> bool:
        """Push new annotations / display flags to a running viewer.

        Mutates this viewer's ``annotation_config`` (and optionally
        ``viewer_config``), rebuilds the per-structure representation list via
        ``_build_viewer_config()``, and pushes it to the live HTTP server. The
        browser's poll loop will pick up the version bump and rebuild
        representations on the existing stage without a page reload.

        Returns True if the live config was updated, False if no viewer is
        running (in which case the caller can decide whether to launch one).
        """
        if annotation_config is not None:
            self.annotation_config = annotation_config
        if viewer_config is not None:
            self.viewer_config = viewer_config

        server = getattr(self, "server", None)
        if server is None or not server.is_running():
            return False

        new_config = self._build_viewer_config()
        server.update_config(new_config)
        return True

    def _build_viewer_config(self) -> Dict:
        """
        Build complete configuration dict for the viewer.

        Produces a per-structure representation list that unifies default
        representations (protein, ligands, ions, waters) with annotation
        representations. The HTML viewer uses this to populate an interactive
        representation manager.

        Returns:
            Configuration dict with structures (each containing representations)
            and focused_mode flag for auto-centering behavior.
        """
        vc = self.viewer_config if self.viewer_config else {}
        show_protein = vc.get('show_protein_ribbon', True)
        show_ligands = vc.get('show_ligands', True)
        show_ions = vc.get('show_ions', False)
        show_waters = vc.get('show_waters', False)
        rep_opacity = vc.get('rep_opacity') or {}

        # Focused mode: when protein and ligands are both hidden, the viewer
        # auto-centers on visible annotation representations instead of the
        # whole structure. Callers (e.g. the ViewerCoordinator) can also
        # request centering with the protein still visible by setting
        # ``focused_mode=True`` directly in viewer_config — useful for the
        # redox detector's "narrow to current site" hook.
        focused_mode = vc.get('focused_mode', not show_protein and not show_ligands)

        structure_colors = [
            '#0078d4', '#ff8c00', '#00a86b', '#9b59b6', '#e74c3c', '#f39c12'
        ]
        use_fixed_colors = len(self.selected_structures) > 1

        structures = []
        for idx, structure_file in enumerate(self.selected_structures):
            structure_name = os.path.basename(structure_file).replace('.pdb', '')
            protein_color = (
                structure_colors[idx % len(structure_colors)]
                if use_fixed_colors else 'residueindex'
            )

            # Default representations
            reps = [
                {
                    'id': 'default_protein',
                    'label': 'Protein',
                    'selection': 'protein',
                    'style': 'cartoon',
                    'color': protein_color,
                    'visible': show_protein,
                },
                {
                    'id': 'default_ligands',
                    'label': 'Ligands',
                    'selection': 'hetero and not (water or ion)',
                    'style': 'ball+stick',
                    'color': 'element',
                    'visible': show_ligands,
                },
                {
                    'id': 'default_ions',
                    'label': 'Ions',
                    'selection': 'ion',
                    'style': 'spacefill',
                    'color': 'element',
                    'visible': show_ions,
                },
                {
                    'id': 'default_waters',
                    'label': 'Waters',
                    'selection': 'water',
                    'style': 'ball+stick',
                    'color': 'element',
                    'visible': show_waters,
                },
                {
                    'id': 'default_nonstandard',
                    'label': 'Non-standard residues',
                    'selection': 'not protein and not hetero and not water and not ion',
                    'style': 'ball+stick',
                    'color': 'element',
                    'visible': True,
                },
            ]

            # Optional per-representation translucency. viewer_config may
            # carry {'rep_opacity': {'default_protein': 0.4}} to render a
            # default rep see-through — useful when the point of the view is
            # something buried inside the protein (a metal site, a channel)
            # that a solid cartoon would hide. Values are 0.0-1.0; the
            # browser exposes the same knob as a per-rep slider, so this only
            # sets the starting state.
            for rep in reps:
                opacity = rep_opacity.get(rep['id'])
                if opacity is not None:
                    rep['opacity'] = max(0.0, min(1.0, float(opacity)))

            # Annotation representations for this structure
            for ann_type, ann_config in self.annotation_config.items():
                structure_indices = ann_config.get('structure_indices', [0])
                if idx not in structure_indices:
                    continue

                # Resolve color for this structure
                if 'colors_per_structure' in ann_config:
                    color = ann_config['colors_per_structure'].get(
                        idx, '#ffffff'
                    )
                elif 'color' in ann_config:
                    color = ann_config['color']
                else:
                    color = '#ffffff'

                label = ann_config.get(
                    'label',
                    ann_type.replace('_', ' ').title()
                )

                rep = {
                    'id': f'ann_{ann_type}',
                    'label': label,
                    'selection': ann_config.get('ngl_selection', ''),
                    'style': ann_config.get('style', 'ball+stick'),
                    'color': color,
                    'visible': True,
                }
                # Optional NGL params (halo overlay uses opacity/scale; the
                # bond style uses atom_pairs/show_labels for NGL's distance
                # representation).
                if 'opacity' in ann_config:
                    rep['opacity'] = ann_config['opacity']
                if 'scale' in ann_config:
                    rep['scale'] = ann_config['scale']
                if 'atom_pairs' in ann_config:
                    rep['atom_pairs'] = ann_config['atom_pairs']
                if 'show_labels' in ann_config:
                    rep['show_labels'] = ann_config['show_labels']
                reps.append(rep)

            structures.append({
                'index': idx,
                'name': structure_name,
                'representations': reps,
            })

        # Coordinator-managed shapes (e.g. centroid spheres) live alongside
        # the per-structure reps — they're attached to the stage via NGL's
        # Shape API, not to a structure component, so they have no
        # structure index. Pass them through as a top-level list keyed by
        # label.
        shape_cfg = getattr(self, 'shape_config', None) or {}
        shapes = []
        for label, sc in shape_cfg.items():
            shapes.append({
                'label': label,
                'type': sc.get('type', 'sphere'),
                'coords': sc.get('coords', [0.0, 0.0, 0.0]),
                'radius': sc.get('radius', 1.0),
                'color': sc.get('color', '#ffaa00'),
                'opacity': sc.get('opacity', 0.5),
            })

        return {
            'structures': structures,
            'shapes': shapes,
            'focused_mode': focused_mode,
        }

    # ========================================================================
    # NGL Selection Conversion Helpers
    # ========================================================================

    def _redox_site_to_ngl_selection(self, redox_site) -> str:
        """
        Convert RedoxSite.residue_groups to NGL selection string.

        Args:
            redox_site: RedoxSite object

        Returns:
            NGL selection string like ":A and 150 or :A and 151 or :B and 45"
        """
        selections = []
        for (chain_id, resid, icode), coords_list in redox_site.residue_groups.items():
            selections.append(f":{chain_id} and {resid}")

        return " or ".join(selections) if selections else ""

    def _protonation_results_to_ngl_selection(self, protonation_results: Dict) -> str:
        """
        Convert protonation_results to NGL selection for titratable residues.

        Args:
            protonation_results: Dict of titratable residues

        Returns:
            NGL selection string
        """
        selections = []
        for key, res_data in protonation_results.items():
            chain = res_data["chain"]
            resnum = res_data["number"]
            # Exclude termini (handled separately if needed)
            if res_data["type"] not in ["N_TERM", "C_TERM"]:
                selections.append(f":{chain} and {resnum}")

        return " or ".join(selections) if selections else ""

    def _filter_selections_to_ngl(self, filter_selections: Dict) -> Dict[str, str]:
        """
        Convert filter_selections to NGL selections for kept vs removed.

        Args:
            filter_selections: Dict of filter selections from workspace

        Returns:
            Dict with 'kept' and 'removed' NGL selections
        """
        kept_selections = []

        for chain_id, components in filter_selections.items():
            for comp_type, residue_list in components.items():
                if residue_list:  # Non-empty = kept
                    for resid in residue_list:
                        kept_selections.append(f":{chain_id} and {resid}")

        kept_selection = " or ".join(kept_selections) if kept_selections else ""

        # Removed selection is the inverse
        removed_selection = f"not ({kept_selection})" if kept_selection else "all"

        return {
            'kept': kept_selection,
            'removed': removed_selection
        }

    def _bonds_to_ngl_selection(self, bonds: List[Dict]) -> str:
        """
        Convert bond list to NGL selection for residues involved.

        Args:
            bonds: List of bond info dicts

        Returns:
            NGL selection string
        """
        residues = set()

        for bond in bonds:
            # Add both residues involved in bond
            res1 = (bond['atom1_residue']['chain'],
                    bond['atom1_residue']['resid'])
            res2 = (bond['atom2_residue']['chain'],
                    bond['atom2_residue']['resid'])
            residues.add(res1)
            residues.add(res2)

        # Convert to NGL selection
        selections = [f":{chain} and {resid}" for chain, resid in residues]
        return " or ".join(selections) if selections else ""

    # ========================================================================
    # Utility Methods
    # ========================================================================

    def _show_annotation_info(self) -> bool:
        """
        Display information about available annotations in workspace.

        Returns:
            True if info displayed successfully
        """
        annotations = self._detect_available_annotations()

        if not annotations:
            self.console.print(Panel(
                "[yellow]No annotations found in workspace[/yellow]\n\n"
                "Available annotation types:\n"
                "- Redox Sites (from Redox Detector)\n"
                "- Titratable Residues (from Protonation Analyzer)\n"
                "- Filter Selections (from PDB Filter)\n"
                "- Bond Classifications (from Redox Sites)\n"
                "- Aligned Structures (from Structure Aligner)",
                title="Available Annotations",
                border_style="yellow"
            ))
            return True

        # Display found annotations grouped by priority
        high_priority = []
        medium_priority = []

        for ann_type, ann_data in annotations.items():
            entry = (ann_type, ann_data)
            if ann_data.get('priority') == 'HIGH':
                high_priority.append(entry)
            else:
                medium_priority.append(entry)

        # Display high priority annotations
        if high_priority:
            table = Table(title="[bold cyan]High Priority Annotations[/bold cyan]")
            table.add_column("Type", style="cyan")
            table.add_column("Count", style="green", justify="right")
            table.add_column("Description", style="white")

            for ann_type, ann_data in high_priority:
                table.add_row(
                    ann_type.replace('_', ' ').title(),
                    str(ann_data.get('count', 'N/A')),
                    ann_data.get('description', '')
                )

            self.console.print(table)

        # Display medium priority annotations
        if medium_priority:
            table = Table(title="[bold yellow]Optional Annotations[/bold yellow]")
            table.add_column("Type", style="yellow")
            table.add_column("Count", style="green", justify="right")
            table.add_column("Description", style="white")

            for ann_type, ann_data in medium_priority:
                table.add_row(
                    ann_type.replace('_', ' ').title(),
                    str(ann_data.get('count', 'N/A')),
                    ann_data.get('description', '')
                )

            self.console.print(table)

        return True
