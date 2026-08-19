#!/usr/bin/env python3
"""
RedoxSite Transformation Manager

Main integration module that replaces the MetalSite transformation system.
Provides the same interface and functionality using RedoxSite objects.

Author: Claude Code Implementation
"""

import copy
import json
import logging
import os
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from rich.console import Console
from rich.panel import Panel
from rich.prompt import Confirm, Prompt
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn

from proprep.redoxsite_prep.transformation.redox_transformer_framework import (
    RedoxSiteTransformerSelector,
    RedoxSiteSpaceAnalyzer,
    TransformationExecutor
)
from proprep.redoxsite_prep.transformation.redox_id_mapper import RedoxSiteIDMapper
from proprep.utils.prompts import prompt_with_context, confirm_with_context

logger = logging.getLogger(__name__)

# Verbosity level constants
class VerbosityLevel:
    """Output verbosity levels for transformation process"""
    FULL = "full"        # All panels, plans, and details
    COMPACT = "compact"  # Summaries and metrics only
    MINIMAL = "minimal"  # One-line updates per site

def parse_bulk_assignment(assignment_str: str, num_sites: int) -> Dict[int, str]:
    """
    Parse bulk transformer assignment string.

    Supports formats:
    - "1-61:heme_bis_his_c_type" - range assignment
    - "1,2,3:no_transformation" - list assignment
    - "1-3:heme, 4-6:no_transformation" - multiple groups
    - "1,3-5,7:disulfide" - mixed list and range

    Args:
        assignment_str: User input string
        num_sites: Total number of sites for validation

    Returns:
        Dict mapping site number (1-indexed) to transformer name

    Raises:
        ValueError: If parsing fails or validation fails
    """
    assignments = {}

    # Split by comma to get individual groups
    groups = [g.strip() for g in assignment_str.split(',')]

    current_transformer = None
    current_sites = []

    for group in groups:
        # Check if this group contains a colon (new transformer assignment)
        if ':' in group:
            # Save previous transformer assignment if exists
            if current_transformer and current_sites:
                for site_num in current_sites:
                    assignments[site_num] = current_transformer
                current_sites = []

            # Parse new assignment
            sites_part, transformer = group.split(':', 1)
            current_transformer = transformer.strip().strip('"\'')
            sites_part = sites_part.strip()

            # Parse site numbers/ranges
            for site_spec in sites_part.split(','):
                site_spec = site_spec.strip()
                if '-' in site_spec:
                    # Range: "1-5"
                    start, end = site_spec.split('-', 1)
                    start_num = int(start.strip())
                    end_num = int(end.strip())
                    if start_num > end_num:
                        raise ValueError(f"Invalid range: {site_spec} (start > end)")
                    if start_num < 1 or end_num > num_sites:
                        raise ValueError(f"Range {site_spec} out of bounds (1-{num_sites})")
                    current_sites.extend(range(start_num, end_num + 1))
                else:
                    # Single site
                    site_num = int(site_spec.strip())
                    if site_num < 1 or site_num > num_sites:
                        raise ValueError(f"Site {site_num} out of bounds (1-{num_sites})")
                    current_sites.append(site_num)
        else:
            # This is a continuation of sites for current transformer
            if not current_transformer:
                raise ValueError(f"No transformer specified for sites: {group}")

            site_spec = group.strip()
            if '-' in site_spec:
                # Range
                start, end = site_spec.split('-', 1)
                start_num = int(start.strip())
                end_num = int(end.strip())
                if start_num > end_num:
                    raise ValueError(f"Invalid range: {site_spec} (start > end)")
                if start_num < 1 or end_num > num_sites:
                    raise ValueError(f"Range {site_spec} out of bounds (1-{num_sites})")
                current_sites.extend(range(start_num, end_num + 1))
            else:
                # Single site
                site_num = int(site_spec.strip())
                if site_num < 1 or site_num > num_sites:
                    raise ValueError(f"Site {site_num} out of bounds (1-{num_sites})")
                current_sites.append(site_num)

    # Save final transformer assignment
    if current_transformer and current_sites:
        for site_num in current_sites:
            if site_num in assignments:
                raise ValueError(f"Site {site_num} assigned multiple times")
            assignments[site_num] = current_transformer

    return assignments


def export_redox_sites_to_json(redox_sites: List, pdb_file,
                               transformer_info=None) -> Path:
    """
    Export RedoxSite objects to a JSON file keyed off the given PDB filename.

    The output file is written to `{pdb_stem}_redox_sites_updated.json` in the
    same directory as `pdb_file`. Residue IDs in the export reflect the current
    state of the RedoxSite objects, so callers are responsible for synchronizing
    sites before calling this helper if the structure has been reordered.

    `transformer_info`, when provided, is the per-site transformer-assignment
    list (the same structure stored under the workspace key ``transformer_info``
    by RedoxTransformationManager._build_transformer_info). It is embedded as a
    top-level ``transformer_info`` key so the Topology Generator can recover the
    site→transformer→forcefield mapping from disk after a ProPrep restart
    (when the in-memory workspace copy is gone). Omitted/None → not written
    (backward compatible with older JSON consumers).

    Returns the output Path.
    """
    export_data = {
        "source_pdb": Path(pdb_file).name,
        "sites": []
    }
    if transformer_info:
        export_data["transformer_info"] = transformer_info

    for site in redox_sites:
        site_data = {
            "site_id": site.site_id,
            "structure_id": site.structure_id,
            "site_type": site.site_type,
            "centers": [],
            "atoms": [],
            "bonds": [],
            "coord_to_pdb": {}
        }

        for center in site.centers:
            center_data = {
                "chain": center.chain,
                "resname": center.resname,
                "resid": center.resid,
                "atom_name": center.atom_name,
                "element": center.element,
                "center_type": center.center_type.value if hasattr(center.center_type, 'value') else str(center.center_type),
                "coordinates": [round(float(x), 3) for x in center.coords]
            }
            site_data["centers"].append(center_data)

        for atom in site.atoms:
            atom_data = {
                "chain": atom.chain,
                "resname": atom.resname,
                "resid": atom.resid,
                "atom_name": atom.atom_name,
                "element": atom.element,
                "coordinates": [round(float(x), 3) for x in atom.coords]
            }
            site_data["atoms"].append(atom_data)

        for bond in site.bonds:
            bond_data = {
                "atom1": bond.atom1_residue_info,
                "atom2": bond.atom2_residue_info,
                "bond_type": bond.bond_type,
                "chemical_type": bond.chemical_type,
                "distance": float(bond.distance),
                "atom1_element": bond.atom1_element,
                "atom2_element": bond.atom2_element,
                "atom1_coordinates": [round(float(x), 3) for x in bond.atom1_coords],
                "atom2_coordinates": [round(float(x), 3) for x in bond.atom2_coords]
            }
            site_data["bonds"].append(bond_data)

        for coord_tuple, pdb_info in site.coord_to_pdb.items():
            coord_str = f"{coord_tuple[0]:.3f},{coord_tuple[1]:.3f},{coord_tuple[2]:.3f}"
            site_data["coord_to_pdb"][coord_str] = pdb_info

        export_data["sites"].append(site_data)

    base_name = Path(pdb_file).stem
    output_file = Path(pdb_file).parent / f"{base_name}_redox_sites_updated.json"

    with open(output_file, 'w') as f:
        json.dump(export_data, f, indent=2)

    return output_file


class RedoxTransformationManager:
    """
    Main manager for RedoxSite transformations
    
    Replaces the MetalSite transformation workflow with RedoxSite coordinate-based tracking.
    Provides the same user interface and functionality as the original system.
    """
    
    def __init__(self, processor=None, console: Console = None):
        self.processor = processor
        self.console = console or Console()

        # Initialize components
        self.transformer_selector = RedoxSiteTransformerSelector(self.console)

        # Register any auto-emitted rename transformers from ~/.proprep/transformers.
        # Built-ins register via static imports in transformers/__init__.py;
        # nothing loads the user dir, so do it here (idempotent, defensive).
        try:
            from proprep.redoxsite_prep.transformation.auto_rename import (
                load_user_transformers,
            )
            newly = load_user_transformers()
            if newly:
                logger.debug("Loaded %d user transformer(s): %s", len(newly), newly)
        except Exception as exc:  # pragma: no cover - never block the preparer
            logger.warning("Could not load user transformers: %s", exc)

        self.space_analyzer = RedoxSiteSpaceAnalyzer(self.console)
        self.id_mapper = RedoxSiteIDMapper(self.console)
        self.transformation_executor = TransformationExecutor(self.console, temp_dir=None)  # Will be set when temp_dir is created

        # State
        self.redox_sites: List = []
        self.selected_transformers: Dict[str, str] = {}
        self.transformation_parameters: Dict[str, Dict[str, Any]] = {}
        # Batch protonation config: site_id -> state_key ("redox_spin") ->
        # {ph_treatment, protonation_<role>...}. Bound per redox state so a
        # site's reduced/oxidized forms may carry different propionate
        # protonation; each microstate inherits the treatment of its states.
        self.site_state_protonation: Dict[str, Dict[str, Dict[str, Any]]] = {}
        self.structure_lines: List[str] = []
        self.id_mapping: Dict[Tuple[str, int], int] = {}

        # Original state tracking (before any transformations)
        self.original_redox_sites: List = []  # Deep copy of sites before any changes
        self.original_atom_states: Dict[str, Dict] = {}  # site_id -> {coord -> atom_info}
        self.transformation_tracking: Dict[str, Dict] = {}  # site_id -> transformation tracking data
        self.executed_transformations: Dict[str, List] = {}  # site_id -> executed transformation list with real atom counts

        # Output tracking
        self.temp_dir: Optional[Path] = None
        self.transformation_reports: Dict[str, Any] = {}
        self.transformation_summary: Dict[str, Any] = {}

        # Batch mode flag (suppresses verbose output) - kept for backwards compatibility
        self._batch_mode: bool = False

        # Verbosity level for output control
        self.verbosity: str = VerbosityLevel.FULL
    
    def transform_redox_sites(self) -> bool:
        """
        Main entry point for RedoxSite transformation process
        
        Replaces the existing transform_metal_sites() method with RedoxSite workflow.
        
        Returns:
            Success status
        """
        try:
            self.console.print("\n[bold blue]===== RedoxSite Transformation Process =====[/bold blue]")

            # Display educational overview
            self._display_transformation_overview()

            # Step 1: Load RedoxSites and structure
            if not self._load_redox_sites_and_structure():
                return False

            # Establish the per-site visual baseline. Show the structure
            # we're about to transform and overlay one ball+stick rep
            # per detected site, palette-coloured by site index. The
            # transformer-assignment prompt later asks the user to type
            # "site N: transformer X" — those numbers map straight back
            # to these colours. Labels keyed by site_id so subsequent
            # hooks (ambiguity, component matching) can refer to the
            # same site without re-establishing identity.
            self._show_redox_site_baseline()

            # Step 1.5: Check site count and offer verbosity options if needed
            self._configure_verbosity()

            # Step 2: Site Selection and Transformer Assignment
            self.console.print("\n[bold]Site Selection[/bold]")
            if not self._configure_site_transformers():
                return False
            
            # Step 3: Space Analysis and ID Mapping
            self.console.print("\n[bold]Residue ID Mapping[/bold]")
            if not self._handle_id_mapping():
                return False
            
            # Step 4: Component Matching and Transformation
            if not self._perform_transformations():
                return False
            
            # Step 5: Site Updates and Output Generation
            if not self._finalize_transformations():
                return False
            
            self.console.print("\n[bold green]✓ Transformation complete[/bold green]")
            return True
            
        except Exception as e:
            import traceback
            logger.error(f"Transformation failed: {e}")
            logger.error(f"Traceback: {traceback.format_exc()}")
            self.console.print(f"[red]Transformation failed: {e}[/red]")
            self.console.print(f"[red]Traceback:[/red]")
            self.console.print(traceback.format_exc())
            return False

    def _show_redox_site_baseline(self):
        """Show the structure and overlay one ball+stick rep per detected
        redox site, palette-coloured by site index.

        All viewer interactions are best-effort — a failure here must
        never block the transformation flow.
        """
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer

            structure_file = self._get_current_structure_file()
            if not structure_file:
                return

            # Auto-fired baseline at start of the transformation flow. No
            # force= so it stays silent in CLI when no viewer is open; pushes
            # live to one that is. Web shell auto-launches via env flag.
            _viewer.show_structure(structure_file)

            for idx, site in enumerate(self.redox_sites, 1):
                pairs = set()
                for c in getattr(site, 'centers', []) or []:
                    chain = getattr(c, 'chain', None)
                    resid = getattr(c, 'resid', None)
                    if chain and resid is not None:
                        pairs.add((chain, resid))
                for a in getattr(site, 'atoms', []) or []:
                    chain = getattr(a, 'chain', None)
                    resid = getattr(a, 'resid', None)
                    if chain and resid is not None:
                        pairs.add((chain, resid))
                if not pairs:
                    continue
                clauses = [f"(:{c} and {r})" for c, r in sorted(pairs)]
                site_id = getattr(site, 'site_id', f'site_{idx}')
                _viewer.highlight(
                    " or ".join(clauses),
                    style="ball+stick",
                    color=f"palette:{idx}",
                    label=f"redox_site_{site_id}",
                )
        except Exception as exc:
            logger.debug("redox baseline viewer hook silenced: %s", exc)

    def _display_transformation_overview(self):
        """Display educational overview of the transformation process.

        Phase titles are bold blue; the role column makes the division of
        labour scannable -- orange "you" (you act) vs purple "ProPrep" (done
        automatically). The takeaway is the asymmetry: you assign transformers
        (Step 1), ProPrep runs everything after.
        """
        from rich.panel import Panel
        from rich.console import Group
        from rich.table import Table
        from rich.text import Text

        TITLE = "bold blue"
        YOU = "bold dark_orange3"      # you act
        PP = "bold dark_violet"        # ProPrep does it automatically

        intro = Text(
            "You assign a transformer to each redox site; ProPrep does the rest.",
            style="grey50")

        guide = Table.grid(padding=(0, 2, 0, 0), expand=True)
        guide.add_column(style=TITLE, no_wrap=True)   # step title
        guide.add_column(no_wrap=True)                # role: you / ProPrep
        guide.add_column(overflow="fold", ratio=1)    # what happens
        steps = [
            ("1. Transformer Assignment", "you",
             "pick which transformer applies to each redox site; choose "
             "'no_transformation' for sites already forcefield-compliant"),
            ("2. Residue ID Mapping", "ProPrep",
             "allocates fresh residue IDs for every site (all sites renumbered) "
             "so transformations never collide"),
            ("3. Component Matching", "ProPrep",
             "identifies the atoms and residues each transformer needs and "
             "validates that the site matches its requirements"),
            ("4. Apply Transformations", "ProPrep",
             "executes each transformation, tracks atom movements through "
             "coordinate mapping, and applies the residue-ID remapping"),
            ("5. Update Workspace", "ProPrep",
             "updates RedoxSite definitions to the transformed state, writes "
             "mapping and validation reports, and saves the structure"),
        ]
        for i, (step, who, what) in enumerate(steps):
            if i:
                guide.add_row("", "", "")  # blank line between items
            label = "You" if who == "you" else who
            role = Text(label, style=(YOU if who == "you" else PP))
            guide.add_row(step, role, what)

        self.console.print()
        self.console.print(Panel(
            Group(intro, Text(""), guide),
            title="[bold blue]How RedoxSite Transformation Works[/bold blue]",
            border_style="blue", padding=(1, 2), expand=True))
        self.console.print()

    def _configure_verbosity(self):
        """Configure output verbosity based on number of sites"""
        # If batch mode is set, use minimal verbosity
        if self._batch_mode:
            self.verbosity = VerbosityLevel.MINIMAL
            self.transformation_executor.verbose = False
            return

        num_sites = len(self.redox_sites)

        # Threshold: offer verbosity options for 20+ sites
        if num_sites < 20:
            self.verbosity = VerbosityLevel.FULL
            self.transformation_executor.verbose = True
            return

        # Show verbosity options
        self.console.print(f"\n[bold]Output Verbosity Configuration[/bold]")
        self.console.print(f"You have {num_sites} redox sites to transform.")
        self.console.print()
        self.console.print("Choose output verbosity:")
        self.console.print("  [cyan](f)ull[/cyan]    - Show all transformation details (educational, more output)")
        self.console.print("  [cyan](c)ompact[/cyan] - Show summaries and metrics only (recommended for >20 sites)")
        self.console.print("  [cyan](m)inimal[/cyan] - One line per site (fastest, least output)")
        self.console.print()
        self.console.print("[grey50]Note: Full transformation details always saved to log files[/grey50]")
        self.console.print()

        choice = prompt_with_context(
            processor=self.processor,
            prompt="Verbosity",
            choices=["f", "c", "m"],
            default="c",
            module="Redox Site Transformation",
            description="Select output verbosity level"
        )

        # Map choice to verbosity level
        verbosity_map = {
            "f": VerbosityLevel.FULL,
            "c": VerbosityLevel.COMPACT,
            "m": VerbosityLevel.MINIMAL
        }

        self.verbosity = verbosity_map.get(choice, VerbosityLevel.COMPACT)
        self.console.print(f"[green]✓ Using {self.verbosity} verbosity mode[/green]")

        # Update transformation executor verbose flag
        self.transformation_executor.verbose = (self.verbosity == VerbosityLevel.FULL)

        # Show log location if not in full mode
        if self.verbosity != VerbosityLevel.FULL and self.temp_dir:
            self.console.print(f"[grey50]Full details saved to: {self.temp_dir / 'transformation_reports'}[/grey50]")
        self.console.print()

    def _load_redox_sites_and_structure(self) -> bool:
        """Load RedoxSites from workspace and current structure"""
        if not self.processor:
            self.console.print("[red]No processor available[/red]")
            return False
        
        # Get workspace
        workspace = self.processor.workspace
        
        # Load RedoxSites
        self.redox_sites = workspace.get("detected_redox_sites", [])
        if not self.redox_sites:
            self.console.print("[yellow]No redox sites detected yet. Please run detection first.[/yellow]")
            return False
        
        # CRITICAL: Capture original states before ANY modifications
        self._capture_original_states()
        
        # Store the original structure before ANY transformations for coordinate mapping
        structure_file = None
        if not workspace.get("untransformed_structure"):
            structure_to_transform = self._get_structure_to_transform()
            if structure_to_transform:
                workspace.set("untransformed_structure", structure_to_transform)
                structure_file = structure_to_transform  # Cache to avoid duplicate message

        # Load current structure (use cached if available)
        if not structure_file:
            structure_file = self._get_current_structure_file()
        if not structure_file or not Path(structure_file).exists():
            self.console.print("[red]No structure file available for transformation[/red]")
            return False
        
        # Read structure lines
        with open(structure_file, 'r') as f:
            self.structure_lines = f.readlines()

        self.console.print(f"[green]Loaded {len(self.redox_sites)} redox sites and structure file[/green]")

        # Create temporary directory for reports
        self.temp_dir = Path(tempfile.mkdtemp(prefix="redox_transformation_"))
        self.temp_dir.mkdir(exist_ok=True)
        (self.temp_dir / "transformation_reports").mkdir(exist_ok=True)
        
        # Update transformation executor with temp directory
        self.transformation_executor.temp_dir = self.temp_dir
        
        return True
    
    def _capture_original_states(self):
        """Capture complete original state of RedoxSites before any transformations or ID mapping"""
        # Capture original states (silently)
        
        # Deep copy the original sites
        import copy
        self.original_redox_sites = copy.deepcopy(self.redox_sites)
        
        # Extract original atom states for each site
        for site in self.redox_sites:
            site_id = site.site_id
            
            # Store original atom information by coordinates
            original_atoms = {}
            for coords, pdb_info in site.coord_to_pdb.items():
                original_atoms[coords] = {
                    'chain': pdb_info.get('chain'),
                    'resname': pdb_info.get('resname'), 
                    'resid': pdb_info.get('resid'),
                    'insertion_code': pdb_info.get('insertion_code', ''),
                    'atom_name': pdb_info.get('atom_name')
                }
            
            self.original_atom_states[site_id] = original_atoms
            
            # Initialize transformation tracking for this site
            self.transformation_tracking[site_id] = {
                'atom_movements': {},  # coord -> list of transformation steps
                'residue_changes': {},  # (chain, resid) -> list of changes
                'splits': [],  # List of residue splits
                'renames': [],  # List of residue renames
                'total_atoms_tracked': len(original_atoms)
            }
        
        total_sites = len(self.redox_sites)
        total_atoms = sum(len(atoms) for atoms in self.original_atom_states.values())
        self.console.print(f"[green]✅ Captured original states: {total_sites} sites, {total_atoms} atoms[/green]")
    
    # Mapping from workspace key to source name for display
    _KEY_TO_SOURCE_NAME = {
        "protonation_pdb_file": "protonation-updated",
        "structure_with_prot_resnames": "protonation-updated",
        "repaired_pdb_file": "repaired",
        "filtered_pdb_file": "filtered",
        "rcsb_pdb_file": "RCSB PDB",
        "local_pdb_file": "local",
        "alphafold_pdb_file": "AlphaFold",
        "alphafold_homolog_pdb_file": "AlphaFold homolog",
        "aligned_target_pdb_file": "aligned target",
        "aligned_ref_pdb_file": "aligned reference",
    }

    def _get_structure_to_transform(self) -> Optional[str]:
        """Get the structure that will be transformed using correct priority order.

        Uses StructureSelector with custom priority and exclude_keys to:
        - Prioritize protonation-updated -> repaired -> filtered -> others
        - Exclude transformed_pdb_file (can't transform already-transformed structure)
        """
        if not self.processor:
            return None

        workspace = self.processor.workspace
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(workspace, self.console, processor=self.processor)

        result = selector.get_structure(
            priority_override=[
                "protonation_pdb_file",
                "structure_with_prot_resnames",  # Legacy key
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
            exclude_keys=["transformed_pdb_file"],  # Can't transform already-transformed
            interactive=True,
            return_key=True,
            silent=True,  # We'll display our own messages
        )

        if result is None:
            self.console.print("[red]No valid PDB file found in workspace[/red]")
            self.console.print("[yellow]Please load a PDB file first[/yellow]")
            return None

        structure_file, workspace_key = result

        # Display which structure we're using
        source_name = self._KEY_TO_SOURCE_NAME.get(workspace_key, "selected")
        self.console.print(
            f"[green]Using {source_name} PDB file for transformation: {structure_file}[/green]"
        )

        return structure_file
    
    def _get_current_structure_file(self) -> Optional[str]:
        """Get the current structure file from workspace"""
        return self._get_structure_to_transform()
    
    def _configure_site_transformers(self) -> bool:
        """Configure transformers for each site (matching original UI)"""
        self.console.print("\n[bold]═══ Configure Site Transformers ═══[/bold]")

        self.console.print(f"\nDetected {len(self.redox_sites)} redox sites requiring configuration.")

        # Display sites table (only once at the beginning)
        self._display_sites_table()

        # Track which step to start from (1 or 2)
        start_step = 1

        while True:
            # Step 1: Assign Transformers
            if start_step <= 1:
                self.console.print("\n[bold]Step 1: Assign Transformers[/bold]")
                self.console.print("──────────────────────────────")

                if not self._assign_transformers():
                    return False

            # Step 2: Configure Parameters
            if start_step <= 2:
                self.console.print("\n[bold]Step 2: Configure Parameters[/bold]")
                self.console.print("──────────────────────────────")

                if not self._configure_parameters():
                    return False

            # Reset start_step for subsequent iterations
            start_step = 1

            # Display configuration summary
            self._display_configuration_summary()

            # Confirm configuration
            choice = prompt_with_context(
                processor=self.processor,
                prompt="\nOptions:\n(a)pply | (e)dit step | (r)estart | (c)ancel\nChoice [a/e/r/c]",
                choices=["a", "e", "r", "c"],
                default="a",
                module="Redox Site Transformation",
                description="Select transformation action",
                options_map={
                    "a": "Apply transformations",
                    "e": "Edit step",
                    "r": "Restart",
                    "c": "Cancel"
                }
            )

            if choice.lower() == "a":
                # Apply transformations
                break
            elif choice.lower() == "e":
                # Edit step - ask which step to edit
                step_choice = prompt_with_context(
                    processor=self.processor,
                    prompt="Which step to edit?\n(1) Assign Transformers\n(2) Configure Parameters\nChoice [1/2]",
                    choices=["1", "2"],
                    default="2",
                    module="Redox Site Transformation",
                    description="Select step to edit"
                )

                if step_choice == "1":
                    # Clear selections and restart from Step 1
                    self.selected_transformers.clear()
                    self.transformation_parameters.clear()
                    self.console.print("[yellow]Restarting from Step 1...[/yellow]")
                    start_step = 1
                else:
                    # Keep transformer selections, only reconfigure parameters
                    self.transformation_parameters.clear()
                    self.console.print("[yellow]Restarting from Step 2...[/yellow]")
                    start_step = 2
                continue
            elif choice.lower() == "r":
                # Restart - clear all selections
                self.selected_transformers.clear()
                self.transformation_parameters.clear()
                self.console.print("[yellow]Restarting configuration...[/yellow]")
                start_step = 1
                continue
            else:  # choice == "c"
                # Cancel
                self.console.print("[yellow]Transformation cancelled[/yellow]")
                return False

        # Store transformer information in workspace for use by tLEaP template
        self._store_transformer_info_in_workspace()

        return True
    
    def _display_sites_table(self):
        """Display sites table with detailed residue composition"""
        table = Table()
        table.add_column("Site #", style="cyan")
        table.add_column("Center", style="yellow")
        table.add_column("Location", style="green")
        table.add_column("Other Residues", style="blue", width=25)
        table.add_column("Compatible Transformers", style="magenta", width=30)
        
        no_match_sites: List[Tuple[int, Any, Dict[str, Any]]] = []

        for i, site in enumerate(self.redox_sites, 1):
            # Get main center info
            if site.centers:
                center = site.centers[0]
                # Center info: use atom_name (e.g. "FE"), fallback to element, then resname
                center_name = center.atom_name or center.element or center.resname or "?"
                location = f"Chain {center.chain}, {center.resname} {center.resid}"
            else:
                center_name = "Unknown"
                location = f"Site {site.site_id}"

            # Get residue composition (excluding the center residue)
            residue_composition = self._get_residue_composition(site)

            # Evaluate every transformer once and reuse the result for the table
            # cell and (if it's a no-match site) the diagnostic block below.
            evaluations = self._evaluate_site_against_all(site)
            compatible = [name for name, ev in evaluations.items() if ev.is_valid]
            compatible_transformers = ", ".join(compatible) if compatible else "None"
            if not compatible:
                no_match_sites.append((i, site, evaluations))

            table.add_row(
                str(i),
                center_name,
                location,
                residue_composition,
                compatible_transformers
            )

        self.console.print(table)
        self.console.print("[grey50]Note: All sites are compatible with 'no_transformation' (pass-through, no modifications).[/grey50]")
        self.console.print("[grey50]The table above shows only specific transformers capable of modifying each site.[/grey50]")

        if no_match_sites:
            # The full per-transformer rejection detail grows with the number of
            # registered transformers, so write it to a report file in the working
            # directory rather than flooding the screen. The on-demand `why <site#>`
            # command still prints a single site's diagnosis inline.
            report_lines: List[str] = [
                "Transformer compatibility report",
                "================================",
                "",
                "The sites below had no specific compatible transformer (only the",
                "no_transformation pass-through applies). Each block lists every",
                "rejected transformer and the requirements it failed.",
                "",
            ]
            for site_idx, site, evals in no_match_sites:
                report_lines.extend(
                    self._render_site_diagnosis_lines(site_idx, site, evaluations=evals)
                )
                report_lines.append("")

            report_path = Path(os.getcwd()) / "transformer_compatibility_report.txt"
            try:
                report_path.write_text("\n".join(report_lines))
                self.console.print(
                    f"\n[grey50]{len(no_match_sites)} site(s) have no specific compatible "
                    f"transformer; detailed diagnostics written to "
                    f"[cyan]{report_path.name}[/cyan].[/grey50]"
                )
            except OSError as e:
                # If the report can't be written, fall back to inline output so the
                # information isn't lost.
                self.console.print(
                    f"[yellow]Could not write compatibility report ({e}); "
                    f"showing diagnostics inline instead.[/yellow]"
                )
                for site_idx, site, evals in no_match_sites:
                    self._format_site_diagnosis(site_idx, site, evaluations=evals)

            self.console.print(
                "[grey50]Tip: at the assignment prompt enter "
                "[cyan]why <site#>[/cyan] to see compatibility for any site.[/grey50]"
            )
    
    def _get_residue_composition(self, site) -> str:
        """Get a summary of residue types in the site (excluding center residue)"""
        residue_counts = {}
        center_residue_keys = set()
        
        # Identify center residue(s) to exclude
        for center in site.centers:
            center_key = (center.chain, center.resid, getattr(center, 'insertion_code', ''))
            center_residue_keys.add(center_key)
        
        # Count all other residues
        for (chain, resid, insertion_code), coords_list in site.residue_groups.items():
            residue_key = (chain, resid, insertion_code)
            if residue_key not in center_residue_keys:
                # Get residue name from first atom
                for coords in coords_list:
                    atom_info = site.coord_to_pdb.get(coords, {})
                    resname = atom_info.get('resname', 'UNK')
                    if resname != 'UNK':
                        residue_counts[resname] = residue_counts.get(resname, 0) + 1
                        break
        
        # Format as "2×HIS, 2×CYS"
        if residue_counts:
            parts = [f"{count}×{resname}" for resname, count in sorted(residue_counts.items())]
            return ", ".join(parts)
        else:
            return "None"
    
    def _get_compatible_transformers(self, site) -> str:
        """Comma-separated names of compatible transformers (excludes no_transformation)."""
        evaluations = self._evaluate_site_against_all(site)
        compatible = [name for name, ev in evaluations.items() if ev.is_valid]
        return ", ".join(compatible) if compatible else "None"

    def _evaluate_site_against_all(self, site) -> Dict[str, Any]:
        """Run every registered transformer's evaluator against ``site``.

        Excludes ``no_transformation`` (the implicit pass-through). When a
        transformer's ``evaluate_redox_site`` raises, a synthetic
        TransformerEvaluation with ``is_valid=False`` and the exception
        message is recorded so the diagnostic UI can still report it.
        """
        from proprep.redoxsite_prep.transformation.redox_transformer_framework import (
            TransformerEvaluation,
        )
        results: Dict[str, Any] = {}
        for name, cls in self.transformer_selector.registry.get_all_transformers().items():
            if name == "no_transformation":
                continue
            try:
                results[name] = cls.evaluate_redox_site(site)
            except Exception as e:
                results[name] = TransformerEvaluation(
                    is_valid=False,
                    confidence=0.0,
                    description="(evaluation raised)",
                    requirements_met=0,
                    total_requirements=0,
                    details=[],
                    error_message=f"{type(e).__name__}: {e}",
                )
        return results

    def _format_site_diagnosis(
        self,
        site_index_1based: int,
        site,
        evaluations: Optional[Dict[str, Any]] = None,
        include_compatible: bool = False,
    ) -> None:
        """Print a per-site transformer-compatibility diagnosis block.

        Always lists rejected transformers with the requirements that failed
        (description plus expected/found values). When ``include_compatible``
        is True, also prints a green checkmark line per matching transformer.
        """
        if evaluations is None:
            evaluations = self._evaluate_site_against_all(site)

        if site.centers:
            c = site.centers[0]
            center_name = c.atom_name or c.element or c.resname or "?"
            loc = f"Chain {c.chain}, {c.resname} {c.resid}"
        else:
            center_name = "Unknown"
            loc = f"site_id={site.site_id}"

        self.console.print(
            f"\n[bold]Site {site_index_1based} ({center_name} — {loc}) "
            f"transformer compatibility:[/bold]"
        )

        if not evaluations:
            self.console.print(
                "  [grey50]No specific transformers registered (only no_transformation).[/grey50]"
            )
            return

        any_compatible = any(ev.is_valid for ev in evaluations.values())
        for name, ev in evaluations.items():
            if ev.is_valid:
                if include_compatible:
                    self.console.print(f"  [green]✓ {name}[/green]")
                continue
            summary = ev.error_message or "requirements not met"
            if ev.total_requirements:
                summary = (
                    f"{ev.requirements_met}/{ev.total_requirements} requirements met — {summary}"
                )
            self.console.print(f"  [red]✗ {name}[/red]: {summary}")
            for d in ev.details:
                if d.passed:
                    continue
                line = f"      - {d.description}"
                if d.value_expected is not None or d.value_found is not None:
                    line += f" (expected {d.value_expected}, got {d.value_found})"
                if d.error_message:
                    line += f" — {d.error_message}"
                self.console.print(line)
        if not any_compatible and not include_compatible:
            self.console.print(
                "  [grey50]Only [italic]no_transformation[/italic] (pass-through) is compatible.[/grey50]"
            )
    
    def _render_site_diagnosis_lines(
        self,
        site_index_1based: int,
        site,
        evaluations: Optional[Dict[str, Any]] = None,
        include_compatible: bool = False,
    ) -> List[str]:
        """Build the per-site transformer-compatibility diagnosis as plain-text
        lines (no rich markup), for writing to the compatibility report file.

        Mirrors the on-screen layout of ``_format_site_diagnosis`` but returns
        strings instead of printing, so the verbose dump can be written to a file
        while the interactive ``why <site#>`` view stays on the console.
        """
        if evaluations is None:
            evaluations = self._evaluate_site_against_all(site)

        if site.centers:
            c = site.centers[0]
            center_name = c.atom_name or c.element or c.resname or "?"
            loc = f"Chain {c.chain}, {c.resname} {c.resid}"
        else:
            center_name = "Unknown"
            loc = f"site_id={site.site_id}"

        lines: List[str] = [
            f"Site {site_index_1based} ({center_name} — {loc}) transformer compatibility:"
        ]

        if not evaluations:
            lines.append("  No specific transformers registered (only no_transformation).")
            return lines

        any_compatible = any(ev.is_valid for ev in evaluations.values())
        for name, ev in evaluations.items():
            if ev.is_valid:
                if include_compatible:
                    lines.append(f"  ✓ {name}")
                continue
            summary = ev.error_message or "requirements not met"
            if ev.total_requirements:
                summary = (
                    f"{ev.requirements_met}/{ev.total_requirements} requirements met — {summary}"
                )
            lines.append(f"  ✗ {name}: {summary}")
            for d in ev.details:
                if d.passed:
                    continue
                line = f"      - {d.description}"
                if d.value_expected is not None or d.value_found is not None:
                    line += f" (expected {d.value_expected}, got {d.value_found})"
                if d.error_message:
                    line += f" — {d.error_message}"
                lines.append(line)
        if not any_compatible and not include_compatible:
            lines.append("  Only no_transformation (pass-through) is compatible.")
        return lines

    def _suggest_transformer_type(self, site) -> str:
        """Suggest the best transformer type for a site"""
        candidates = []
        try:
            for transformer_name, transformer_class in self.transformer_selector.registry.get_all_transformers().items():
                evaluation = transformer_class.evaluate_redox_site(site)
                if evaluation.is_valid:
                    candidates.append((transformer_name, evaluation.confidence))
        except Exception:
            pass
        
        if candidates:
            # Sort by confidence and return the best
            candidates.sort(key=lambda x: x[1], reverse=True)
            return candidates[0][0]
        
        return "unknown"
    
    def _assign_transformers(self) -> bool:
        """Assign transformers to sites with pre-validated compatibility"""
        # Check for pre-existing transformer mappings from detection phase
        workspace = self.processor._get_workspace()
        transformer_mappings = workspace.get('redox_transformer_mappings', {})
        
        # Build compatibility matrix for all sites
        site_compatibility = {}
        all_compatible_transformers = set()
        
        for i, site in enumerate(self.redox_sites, 1):
            compatible_transformers = []
            for transformer_name, transformer_class in self.transformer_selector.registry.get_all_transformers().items():
                try:
                    evaluation = transformer_class.evaluate_redox_site(site)
                    if evaluation.is_valid:
                        compatible_transformers.append(transformer_name)
                        all_compatible_transformers.add(transformer_name)
                except Exception:
                    continue
            site_compatibility[i] = compatible_transformers
        
        if not all_compatible_transformers:
            self.console.print("[red]No transformers are compatible with any of the detected sites[/red]")
            self._show_compatibility_diagnostics(site_compatibility)
            return False
        
        # Try automatic assignment using transformer mappings
        if transformer_mappings:
            auto_assigned = self._try_automatic_assignment(transformer_mappings, site_compatibility)
            if auto_assigned:
                return True
        else:
            logger.debug("No transformer mappings in workspace - falling back to manual assignment")

        # Display only compatible transformers
        self.console.print("\nAvailable compatible transformers:")
        transformer_list = sorted(list(all_compatible_transformers))
        for i, name in enumerate(transformer_list, 1):
            # Show which sites are compatible with this transformer
            compatible_sites = [site_num for site_num, compat_list in site_compatibility.items() if name in compat_list]
            sites_str = ",".join(map(str, compatible_sites))
            self.console.print(f"{i}. {name} (compatible with sites: {sites_str})")

        # Show example format tailored to the actual number of sites and
        # available transformers.
        n_sites = len(self.redox_sites)
        self.console.print("\nExample format:")
        if n_sites == 1:
            self.console.print(f'  1:"{transformer_list[0]}"', highlight=False)
            if len(transformer_list) > 1:
                self.console.print(
                    f'(With one site you may also enter just the transformer index, '
                    f'e.g. {transformer_list.index(transformer_list[0]) + 1}.)',
                    highlight=False,
                )
            else:
                self.console.print('(With one site and one transformer you may also enter just "1".)', highlight=False)
        elif len(transformer_list) == 1:
            self.console.print(f'  1-{n_sites}:"{transformer_list[0]}"', highlight=False)
        else:
            mid = max(1, n_sites // 2)
            self.console.print(
                f'  1-{mid}:"{transformer_list[0]}", '
                f'{mid + 1}-{n_sites}:"{transformer_list[1]}"',
                highlight=False,
            )

        # Interactive assignment with validation
        n_sites_total = len(self.redox_sites)
        self.console.print(
            "Enter [cyan]why <site#>[/cyan] to see why a transformer was rejected "
            "for a given site."
        )
        while True:
            assignment_input = prompt_with_context(
                processor=self.processor,
                prompt="\nEnter transformer assignments",
                module="Redox Site Transformation",
                description="Assign transformers to sites"
            ).strip()

            if not assignment_input:
                self.console.print("[red]No assignments provided[/red]")
                return False

            # Diagnostic command: 'why N' prints full per-transformer evaluation
            # for site N (compatible AND rejected) and re-prompts.
            lowered = assignment_input.lower()
            if lowered == "why" or lowered.startswith("why "):
                arg = assignment_input[3:].strip() if lowered != "why" else ""
                if not arg:
                    self.console.print(
                        f"[yellow]Usage: why <site_number> "
                        f"(1-{n_sites_total})[/yellow]"
                    )
                    continue
                try:
                    site_num = int(arg)
                except ValueError:
                    self.console.print(
                        f"[yellow]Site number must be an integer "
                        f"(1-{n_sites_total}); got {arg!r}[/yellow]"
                    )
                    continue
                if not (1 <= site_num <= n_sites_total):
                    self.console.print(
                        f"[yellow]Site number out of range: {site_num} "
                        f"(valid: 1-{n_sites_total})[/yellow]"
                    )
                    continue
                self._format_site_diagnosis(
                    site_num, self.redox_sites[site_num - 1], include_compatible=True
                )
                continue
            break
        
        # Parse assignments with pre-validated compatibility
        assignments = []
        n_sites = len(self.redox_sites)
        for assignment in assignment_input.split(','):
            assignment = assignment.strip()
            if not assignment:
                continue
            if ':' not in assignment:
                # Shorthand: with a single site, allow either bare transformer
                # index (e.g. "1") or bare transformer name. Maps the entry to
                # site 1 with that transformer.
                if n_sites == 1:
                    if assignment.isdigit():
                        idx = int(assignment) - 1
                        if 0 <= idx < len(transformer_list):
                            assignment = f'1:"{transformer_list[idx]}"'
                        else:
                            self.console.print(
                                f"[yellow]Invalid transformer index: {assignment}. "
                                f"Expected 1-{len(transformer_list)}.[/yellow]"
                            )
                            continue
                    elif assignment in all_compatible_transformers:
                        assignment = f'1:"{assignment}"'
                    else:
                        self.console.print(
                            f"[yellow]Invalid entry '{assignment}': not a transformer index "
                            f"(1-{len(transformer_list)}) or compatible transformer name.[/yellow]"
                        )
                        continue
                else:
                    self.console.print(
                        f"[yellow]Skipping entry '{assignment}': missing ':' separator. "
                        f"Expected format <site_range>:<transformer>, e.g. 1-3:\"{transformer_list[0]}\".[/yellow]"
                    )
                    continue

            site_range, transformer_spec = assignment.split(':', 1)
            site_range = site_range.strip()
            transformer_spec = transformer_spec.strip()

            try:
                # Handle both quoted names and numeric indices for backwards compatibility
                if transformer_spec.startswith('"') and transformer_spec.endswith('"'):
                    # Quoted transformer name (new format)
                    transformer_name = transformer_spec.strip('"')
                elif transformer_spec.startswith("'") and transformer_spec.endswith("'"):
                    # Single-quoted transformer name
                    transformer_name = transformer_spec.strip("'")
                else:
                    # Numeric index (old format - backwards compatibility)
                    try:
                        transformer_idx = int(transformer_spec) - 1
                        if not (0 <= transformer_idx < len(transformer_list)):
                            self.console.print(f"[yellow]Invalid transformer index: {transformer_idx + 1}[/yellow]")
                            continue
                        transformer_name = transformer_list[transformer_idx]
                    except ValueError:
                        # Not a number, treat as unquoted name
                        transformer_name = transformer_spec

                # Validate transformer exists
                if transformer_name not in all_compatible_transformers:
                    self.console.print(f"[yellow]Unknown transformer: {transformer_name}[/yellow]")
                    continue

                # Parse site range (supports ranges like 1-10 or individual sites like 5,7,9)
                site_numbers = []
                if '-' in site_range:
                    start, end = site_range.split('-', 1)
                    start_idx = int(start.strip())
                    end_idx = int(end.strip())
                    site_numbers = list(range(start_idx, end_idx + 1))
                else:
                    site_numbers = [int(site_range)]

                # Validate site numbers and compatibility
                for site_num in site_numbers:
                    if 1 <= site_num <= len(self.redox_sites):
                        # Check pre-calculated compatibility
                        if transformer_name in site_compatibility[site_num]:
                            site = self.redox_sites[site_num - 1]
                            self.selected_transformers[site.site_id] = transformer_name
                            assignments.append((site_num, transformer_name))
                        else:
                            self.console.print(f"[yellow]Warning: {transformer_name} not compatible with site {site_num}[/yellow]")
                    else:
                        self.console.print(f"[yellow]Invalid site number: {site_num}[/yellow]")

            except (ValueError, IndexError) as e:
                self.console.print(f"[yellow]Invalid assignment: {assignment} - {e}[/yellow]")
                continue
        
        if not self.selected_transformers:
            self.console.print("[red]No valid assignments made[/red]")
            return False
        
        # Display assignments
        self.console.print("\n✓ Assignments:")
        assignment_groups = {}
        for site_id, transformer in self.selected_transformers.items():
            if transformer not in assignment_groups:
                assignment_groups[transformer] = []
            # Find site number
            site_num = None
            for i, site in enumerate(self.redox_sites):
                if site.site_id == site_id:
                    site_num = i + 1
                    break
            if site_num:
                assignment_groups[transformer].append(site_num)
        
        for transformer, site_nums in assignment_groups.items():
            site_nums.sort()
            if len(site_nums) > 1 and all(site_nums[i] + 1 == site_nums[i+1] for i in range(len(site_nums)-1)):
                # Consecutive range
                sites_str = f"{site_nums[0]}-{site_nums[-1]}"
            else:
                sites_str = ",".join(map(str, site_nums))
            self.console.print(f"  Sites {sites_str}: {transformer}")
        
        return True
    
    def _try_automatic_assignment(self, transformer_mappings: Dict[str, str], site_compatibility: Dict[int, List[str]]) -> bool:
        """
        Try to automatically assign transformers based on detection phase mappings.
        
        Args:
            transformer_mappings: Dict mapping site types to transformer names from detection
            site_compatibility: Dict mapping site numbers to compatible transformer lists
            
        Returns:
            True if all sites were successfully auto-assigned, False otherwise
        """
        self.console.print("[cyan]Found transformer mappings from detection phase. Attempting automatic assignment...[/cyan]")
        self.console.print()
        self.console.print("[grey50]Note: 'No mapping found' occurs when the template name from detection doesn't match a transformer.[/grey50]")
        self.console.print("[grey50]This can happen, for example, when you assign distinct template names that don't match a transformer[/grey50]")
        self.console.print("[grey50]to group sites with different detection criteria. If desired, all those sites can be treated[/grey50]")
        self.console.print("[grey50]with the same or different transformers.[/grey50]")
        self.console.print("[grey50]Example: 'ca_binding_site_type1' and 'ca_binding_site_type2' both use 'no_transformation'[/grey50]", highlight=False)
        self.console.print()

        auto_assignments = {}
        failed_assignments = []
        
        for i, site in enumerate(self.redox_sites, 1):
            site_type = getattr(site, 'site_type', None)
            
            if site_type and site_type in transformer_mappings:
                suggested_transformer = transformer_mappings[site_type]
                
                # Check if the suggested transformer is compatible with this site
                if suggested_transformer in site_compatibility[i]:
                    auto_assignments[site.site_id] = suggested_transformer
                    self.console.print(f"  Site {i} ({site_type}): {suggested_transformer} ✓")
                else:
                    failed_assignments.append((i, site_type, suggested_transformer))
                    self.console.print(f"  Site {i} ({site_type}): {suggested_transformer} [red]✗ Not compatible[/red]")
            else:
                failed_assignments.append((i, site_type or "unknown", "no mapping"))
                self.console.print(f"  Site {i} ({site_type or 'unknown'}): [yellow]No mapping found[/yellow]")
        
        if not failed_assignments:
            # All sites successfully auto-assigned
            self.selected_transformers = auto_assignments
            self.console.print(f"\n[green]✓ Successfully auto-assigned transformers for all {len(auto_assignments)} sites[/green]")
            
            # Display assignment summary
            self.console.print("\n✓ Automatic Assignments:")
            assignment_groups = {}
            for site_id, transformer in self.selected_transformers.items():
                if transformer not in assignment_groups:
                    assignment_groups[transformer] = []
                # Find site number
                site_num = None
                for j, site in enumerate(self.redox_sites):
                    if site.site_id == site_id:
                        site_num = j + 1
                        break
                if site_num:
                    assignment_groups[transformer].append(site_num)
            
            for transformer, site_nums in assignment_groups.items():
                site_nums.sort()
                if len(site_nums) > 1 and all(site_nums[k] + 1 == site_nums[k+1] for k in range(len(site_nums)-1)):
                    # Consecutive range
                    sites_str = f"{site_nums[0]}-{site_nums[-1]}"
                else:
                    sites_str = ",".join(map(str, site_nums))
                self.console.print(f"  Sites {sites_str}: {transformer}")
            
            return True
        else:
            # Some assignments failed - fall back to manual assignment
            self.console.print(f"\n[yellow]Automatic assignment failed for {len(failed_assignments)} sites. Falling back to manual assignment.[/yellow]")
            return False
    
    def _show_compatibility_diagnostics(self, site_compatibility: Dict[int, List[str]]):
        """Show detailed compatibility diagnostics when no transformers are compatible"""
        self.console.print("\n[bold yellow]Compatibility Diagnostics:[/bold yellow]")
        
        # Show what each transformer requires vs what each site has
        all_transformers = self.transformer_selector.registry.get_all_transformers()
        
        for i, site in enumerate(self.redox_sites, 1):
            self.console.print(f"\n[bold]Site {i} Analysis:[/bold]")
            
            # Show site composition
            residue_composition = self._get_residue_composition(site)
            center_info = "Unknown"
            if site.centers:
                center = site.centers[0]
                center_info = f"{center.element if hasattr(center, 'element') else center.resname} at {center.chain}:{center.resid}"
            
            self.console.print(f"  Center: {center_info}")
            self.console.print(f"  Other Residues: {residue_composition}")
            
            # Test against each transformer and show why it failed
            for transformer_name, transformer_class in all_transformers.items():
                try:
                    evaluation = transformer_class.evaluate_redox_site(site)
                    if not evaluation.is_valid:
                        self.console.print(f"  [red]✗ {transformer_name}:[/red] {evaluation.error_message or 'Requirements not met'}")
                        
                        # Show specific requirement failures
                        if evaluation.details:
                            for detail in evaluation.details:
                                if not detail.passed:
                                    self.console.print(f"    - {detail.description}: Expected {detail.value_expected}, got {detail.value_found}")
                except Exception as e:
                    self.console.print(f"  [red]✗ {transformer_name}:[/red] Evaluation failed: {e}")
        
        self.console.print(f"\n[yellow]Consider creating a custom transformer for these site types or check if the site detection parameters need adjustment.[/yellow]")
    
    def _configure_parameters(self) -> bool:
        """Configure parameters for each transformer with full interactive interface"""
        # Brief preamble explaining the flow.
        self.console.print(
            "[grey50]For each parameter we'll show the valid options (with descriptions where "
            "available) and let you assign sites either in one shot using batch syntax "
            "(<site_range>:<value>, ...) or one option at a time. Sites left unassigned at "
            "the end take the parameter's default.[/grey50]"
        )

        # Group sites by transformer
        by_transformer = {}
        for site_id, transformer_name in self.selected_transformers.items():
            if transformer_name not in by_transformer:
                by_transformer[transformer_name] = []
            by_transformer[transformer_name].append(site_id)
        
        # Configure parameters for each transformer group
        for transformer_name, site_ids in by_transformer.items():
            transformer_class = self.transformer_selector.registry.get_transformer(transformer_name)
            if not transformer_class:
                continue
            
            # Get site numbers for display
            site_numbers = []
            for site_id in site_ids:
                for i, site in enumerate(self.redox_sites):
                    if site.site_id == site_id:
                        site_numbers.append(i + 1)
                        break
            
            self.console.print(f"\n[Transformer: {transformer_name}] Sites: {', '.join(map(str, site_numbers))}")

            param_defs = transformer_class.get_parameter_definitions()

            # Check if transformer has no parameters
            if not param_defs:
                self.console.print("[grey50]No parameters need to be set.[/grey50]")
                # Initialize empty parameters for these sites
                for site_id in site_ids:
                    self.transformation_parameters[site_id] = {}
                continue

            # Initialize site parameters
            for site_id in site_ids:
                self.transformation_parameters[site_id] = {}

            for param_name, param_def in param_defs.items():
                if param_def["type"] == "choice":
                    self._configure_choice_parameter_interactive(
                        param_name, param_def, site_ids, site_numbers, transformer_name
                    )
                elif param_def["type"] == "fixed":
                    self._configure_fixed_parameter(
                        param_name, param_def, site_ids, site_numbers, transformer_name
                    )
        
        return True
    
    def _configure_choice_parameter_interactive(self, param_name: str, param_def: Dict[str, Any],
                                              site_ids: List[str], site_numbers: List[int],
                                              transformer_name: str):
        """Configure a choice parameter, with per-site valid-option filtering,
        an info panel showing per-option descriptions, and a batch-syntax
        entry point.

        Flow:
        1. Compute valid options per site (transformer's get_valid_options
           uses already-chosen params for that site, so e.g. a site already
           set to redox_state='oxidized' only sees FO1-FO6 for spin_state).
        2. Render an info panel listing each option in the union with its
           description (when get_option_description returns text).
        3. Offer a batch-syntax prompt accepting <site_range>:<value> pairs.
           Empty input falls through to the existing per-option loop, so
           recorded sessions keep replaying through the per-option prompts.
        4. Per-option fallback: only iterate options valid for at least one
           still-unassigned site.
        5. Any site still unassigned at the end gets the parameter's default
           (when that default is valid for the site).
        """
        from proprep.redoxsite_prep.transformation.redox_transformer_framework import redox_transformer_registry

        from rich.table import Table

        description = param_def["description"]
        static_options = list(param_def.get("options") or [])
        default = param_def.get("default")

        transformer_class = redox_transformer_registry.get_transformer(transformer_name)

        # Per-site valid options. With the base default impl this just returns
        # the static options list for every site; transformers that override
        # get_valid_options can return narrower lists based on each site's
        # already-chosen parameters.
        valid_options_by_site: Dict[str, List[Any]] = {}
        for site_id in site_ids:
            current_params = self.transformation_parameters.get(site_id, {}).copy()
            current_params.pop(param_name, None)
            try:
                valid = transformer_class.get_valid_options(param_name, current_params)
            except Exception as e:
                logger.warning(
                    "get_valid_options failed for %s/%s on site %s: %s. "
                    "Falling back to static options.",
                    transformer_name, param_name, site_id, e,
                )
                valid = list(static_options)
            valid_options_by_site[site_id] = list(valid) or list(static_options)

        # Union of valid options across sites, preserving the order they
        # appear in static_options (or get_valid_options output).
        seen = set()
        union_options: List[Any] = []
        for site_id in site_ids:
            for opt in valid_options_by_site[site_id]:
                if opt not in seen:
                    seen.add(opt)
                    union_options.append(opt)

        if not union_options:
            # Dim, not red: an empty option set here is usually an intentional
            # gate (e.g. a per-ring protonation_<role> choice is only applicable
            # under fixed_pH), not an error — red made users read it as a failure.
            self.console.print(
                f"[grey50]No options to set for parameter '{param_name}' in this "
                f"configuration. Skipping.[/grey50]"
            )
            return

        self.console.print(f"\n{description} ({'/'.join(map(str, union_options))}):")

        # Info panel: only render if at least one option has a description.
        descriptions: Dict[Any, str] = {}
        for opt in union_options:
            try:
                descriptions[opt] = transformer_class.get_option_description(param_name, opt) or ""
            except Exception as e:
                logger.warning(
                    "get_option_description failed for %s/%s/%s: %s",
                    transformer_name, param_name, opt, e,
                )
                descriptions[opt] = ""
        if any(descriptions.values()):
            info = Table(show_header=True, header_style="grey50", box=None, padding=(0, 1))
            info.add_column("Option", style="cyan", no_wrap=True)
            info.add_column("Description", style="grey50")
            for opt in union_options:
                info.add_row(str(opt), descriptions.get(opt, ""))
            self.console.print(info)

        # If sites have different valid-option subsets, surface that so the
        # user knows their choices may not apply universally.
        site_subsets = {site_id: tuple(valid_options_by_site[site_id]) for site_id in site_ids}
        unique_subsets = set(site_subsets.values())
        if len(unique_subsets) > 1:
            self.console.print(
                "[grey50]Note: valid options differ across sites based on their "
                "already-chosen parameters. Each site will only accept options "
                "in its own valid set.[/grey50]"
            )

        assigned_sites: set = set()

        # Batch-syntax entry point. Accepts the same <site_range>:<value>
        # format as Step 1's transformer assignment. Quoted/unquoted values
        # are both accepted. Empty input drops into the per-option fallback,
        # preserving session-replay behavior for sessions recorded under the
        # old flow.
        site_numbers_str = ",".join(map(str, site_numbers))
        n_sites = len(site_numbers)
        if n_sites == 1:
            batch_example = f'1:{union_options[0]}'
        elif len(union_options) >= 2:
            mid = max(1, n_sites // 2)
            batch_example = (
                f'1-{mid}:{union_options[0]}, '
                f'{mid + 1}-{site_numbers[-1]}:{union_options[1]}'
            )
        else:
            batch_example = f'{site_numbers[0]}-{site_numbers[-1]}:{union_options[0]}'
        # Lay the help out as short, label-aligned lines instead of one long
        # wrapping sentence. Reads faster and survives narrow terminals
        # without breaking mid-example.
        if n_sites == 1:
            # highlight=False: the example contains numbers and "<sites>:<value>"
            # tokens that Rich's ReprHighlighter would otherwise colour
            # (numbers cyan, value tokens green), making the help look styled.
            self.console.print(
                f"  Enter a value ({' or '.join(union_options[:3])}"
                f"{' …' if len(union_options) > 3 else ''}), "
                f"or batch syntax {batch_example}",
                highlight=False,
            )
            self.console.print("  Press Enter to fall back to per-option prompts.")
        else:
            second_value = union_options[1] if len(union_options) >= 2 else union_options[0]
            self.console.print(
                "  Syntax:  <sites>:<value>, ...    sites = N, N-M, or N,M",
                highlight=False,
            )
            self.console.print(f"  Example: {batch_example}", highlight=False)
            self.console.print(
                f"           1,3:{union_options[0]}, 2,4:{second_value}",
                highlight=False,
            )
            self.console.print("  Press Enter to fall back to per-option prompts.")
        self.console.print()

        batch_input = prompt_with_context(
            processor=self.processor,
            prompt=f"  Batch assignment for {param_name} (sites {site_numbers_str})",
            default="",
            module="Redox Site Transformation",
            description=f"Batch-assign '{param_name}' across sites; empty drops to per-option prompts",
        ).strip()

        if batch_input:
            # Single-site convenience: if there's exactly one site and the user
            # typed a bare value matching one of the options, treat it as
            # "<site>:<value>" so they don't have to prefix the site number.
            if (n_sites == 1
                    and ':' not in batch_input
                    and batch_input.strip().strip("'").strip('"') in union_options):
                value = batch_input.strip().strip("'").strip('"')
                batch_input = f"{site_numbers[0]}:{value}"
            self._apply_batch_parameter_assignment(
                batch_input=batch_input,
                param_name=param_name,
                site_ids=site_ids,
                site_numbers=site_numbers,
                valid_options_by_site=valid_options_by_site,
                union_options=union_options,
                assigned_sites=assigned_sites,
            )

        # Per-option fallback for any sites still unassigned. Only iterate
        # options that remain valid for at least one unassigned site.
        for option in union_options:
            unassigned_numbers = []
            for i, site_id in enumerate(site_ids):
                if site_id in assigned_sites:
                    continue
                if option in valid_options_by_site[site_id]:
                    unassigned_numbers.append(site_numbers[i])
            if not unassigned_numbers:
                continue

            unassigned_str = ",".join(map(str, unassigned_numbers))
            site_input = prompt_with_context(
                processor=self.processor,
                prompt=f"  Sites for '{option}' (numbers/ranges, empty=skip)",
                default="",
                module="Redox Site Transformation",
                description=f"Which of sites ({unassigned_str}) should get '{option}'?",
            ).strip()

            if not site_input:
                continue

            try:
                assigned_numbers = []
                for part in site_input.split(','):
                    part = part.strip()
                    if not part:
                        continue
                    if '-' in part:
                        start, end = part.split('-', 1)
                        assigned_numbers.extend(range(int(start), int(end) + 1))
                    else:
                        assigned_numbers.append(int(part))

                valid_assignments = []
                for site_num in assigned_numbers:
                    if site_num not in site_numbers:
                        self.console.print(
                            f"    [yellow]Site {site_num} is not in this transformer group; "
                            f"skipping.[/yellow]"
                        )
                        continue
                    site_idx = site_numbers.index(site_num)
                    site_id = site_ids[site_idx]
                    if site_id in assigned_sites:
                        continue
                    if option not in valid_options_by_site[site_id]:
                        self.console.print(
                            f"    [yellow]Option '{option}' is not valid for site {site_num} "
                            f"(valid: {valid_options_by_site[site_id]}); skipping.[/yellow]"
                        )
                        continue
                    self.transformation_parameters[site_id][param_name] = option
                    assigned_sites.add(site_id)
                    valid_assignments.append(site_num)

                if valid_assignments:
                    self.console.print(
                        f"    Assigned {option} to sites: {','.join(map(str, valid_assignments))}"
                    )

            except ValueError:
                self.console.print(f"    [yellow]Invalid input: {site_input}[/yellow]")

        # Auto-assign default to any sites still unassigned, when the default
        # is valid for that site.
        for i, site_id in enumerate(site_ids):
            if site_id in assigned_sites:
                continue
            site_default = default
            if site_default not in valid_options_by_site[site_id]:
                if valid_options_by_site[site_id]:
                    site_default = valid_options_by_site[site_id][0]
                else:
                    self.console.print(
                        f"  [yellow]Site {site_numbers[i]}: no valid value to auto-assign for "
                        f"'{param_name}'.[/yellow]"
                    )
                    continue
            self.transformation_parameters[site_id][param_name] = site_default

        self.console.print()
        self.console.print("  ✓ All sites assigned")

    def _apply_batch_parameter_assignment(self, batch_input: str, param_name: str,
                                          site_ids: List[str], site_numbers: List[int],
                                          valid_options_by_site: Dict[str, List[Any]],
                                          union_options: List[Any],
                                          assigned_sites: set) -> None:
        """Parse and apply a batch-syntax parameter assignment.

        Splits ``batch_input`` on commas. Each token is either a *site fragment*
        (no ``:``) or a *site:value* pair. Site fragments accumulate into a
        pending-sites buffer that empties whenever a ``site:value`` token is
        seen — the value applies to (pending sites + the sites in the value's
        own site_range).

        This lets the user write any of:

            1-2:reduced, 3-4:oxidized              (canonical)
            1:reduced, 2:reduced, 3:oxidized        (each entry standalone)
            1,3:oxidized, 2,4:reduced               (comma-list of sites)
            1, 3:oxidized, 2, 4:reduced             (same, with spaces)

        Site fragments support ``N`` and ``N-M`` range forms. Value can be
        quoted or unquoted; must be one of ``union_options``. Mutates
        ``self.transformation_parameters`` and ``assigned_sites`` in place.
        """

        def parse_site_fragment(frag: str) -> Optional[List[int]]:
            """Parse a single site fragment ('1' or '1-3') into a list of ints.
            Returns None on parse failure (caller logs)."""
            frag = frag.strip()
            if not frag:
                return []
            try:
                if '-' in frag:
                    start, end = frag.split('-', 1)
                    return list(range(int(start.strip()), int(end.strip()) + 1))
                return [int(frag)]
            except ValueError:
                return None

        pending_sites: List[int] = []  # Site numbers accumulated since last :value
        for raw in batch_input.split(','):
            entry = raw.strip()
            if not entry:
                continue

            if ':' not in entry:
                # Pure site fragment — append to pending buffer
                nums = parse_site_fragment(entry)
                if nums is None:
                    self.console.print(
                        f"    [yellow]Skipping batch entry '{entry}': invalid site range "
                        f"(expected N or N-M).[/yellow]"
                    )
                    continue
                pending_sites.extend(nums)
                continue

            # entry contains ':' — site:value pair. The site portion is its own
            # fragment; combine with anything in the pending buffer.
            site_range, value_spec = entry.split(':', 1)
            value_spec = value_spec.strip()
            if (value_spec.startswith('"') and value_spec.endswith('"')) or \
               (value_spec.startswith("'") and value_spec.endswith("'")):
                value_spec = value_spec[1:-1]

            if value_spec not in union_options:
                self.console.print(
                    f"    [yellow]Skipping batch entry '{entry}': value '{value_spec}' is not "
                    f"a valid option (valid: {union_options}).[/yellow]"
                )
                pending_sites = []
                continue

            nums = parse_site_fragment(site_range)
            if nums is None:
                self.console.print(
                    f"    [yellow]Skipping batch entry '{entry}': invalid site range "
                    f"'{site_range.strip()}'.[/yellow]"
                )
                pending_sites = []
                continue

            combined = pending_sites + nums
            pending_sites = []

            applied = []
            for site_num in combined:
                if site_num not in site_numbers:
                    self.console.print(
                        f"    [yellow]Site {site_num} is not in this transformer group; "
                        f"skipping.[/yellow]"
                    )
                    continue
                idx = site_numbers.index(site_num)
                site_id = site_ids[idx]
                if site_id in assigned_sites:
                    self.console.print(
                        f"    [yellow]Site {site_num} already assigned; skipping in this "
                        f"batch entry.[/yellow]"
                    )
                    continue
                if value_spec not in valid_options_by_site[site_id]:
                    self.console.print(
                        f"    [yellow]'{value_spec}' is not a valid option for site {site_num} "
                        f"(valid: {valid_options_by_site[site_id]}); skipping.[/yellow]"
                    )
                    continue
                self.transformation_parameters[site_id][param_name] = value_spec
                assigned_sites.add(site_id)
                applied.append(site_num)
            if applied:
                self.console.print(
                    f"    Assigned {value_spec} to sites: {','.join(map(str, applied))}"
                )

        # Trailing pending sites with no :value attached — warn the user.
        if pending_sites:
            self.console.print(
                f"    [yellow]Trailing site numbers {pending_sites} had no ':<value>' "
                f"attached; skipping. (Sites need a :value pair to be assigned.)[/yellow]"
            )
    
    def _configure_fixed_parameter(self, param_name: str, param_def: Dict[str, Any], 
                                  site_ids: List[str], site_numbers: List[int], 
                                  transformer_name: str):
        """Configure a fixed parameter that has only one valid value"""
        fixed_value = param_def["value"]
        description = param_def["description"]
        note = param_def.get("note", "")
        
        self.console.print(f"\n{description}:")
        if note:
            self.console.print(f"  [grey50]{note}[/grey50]")
        
        sites_str = ",".join(map(str, site_numbers))
        self.console.print(f"  [green]Setting '{fixed_value}' for all sites ({sites_str})[/green]")
        
        # Assign fixed value to all sites
        for site_id in site_ids:
            self.transformation_parameters[site_id][param_name] = fixed_value
        
        self.console.print("  ✓ All sites configured")
    
    def _display_configuration_summary(self):
        """Display configuration summary"""
        self.console.print("\n[bold]═══ Configuration Summary ═══[/bold]")
        
        # Group by parameters for display
        param_groups = {}
        for site_id, params in self.transformation_parameters.items():
            transformer = self.selected_transformers[site_id]
            param_key = (transformer, tuple(sorted(params.items())))
            if param_key not in param_groups:
                param_groups[param_key] = []
            param_groups[param_key].append(site_id)
        
        # Convert site IDs to site numbers for display
        for (transformer, param_items), site_ids in param_groups.items():
            site_numbers = []
            for site_id in site_ids:
                for i, site in enumerate(self.redox_sites):
                    if site.site_id == site_id:
                        site_numbers.append(i + 1)
                        break
            
            param_str = ", ".join(f"{k}={v}" for k, v in param_items)
            sites_str = ",".join(map(str, sorted(site_numbers)))
            self.console.print(f"Sites {sites_str}: {transformer} ({param_str})")
    
    def _handle_id_mapping(self) -> bool:
        """Handle residue ID space analysis and mapping"""
        # Analyze space requirements
        space_analysis = self.space_analyzer.analyze_transformation_space(
            self.redox_sites, self.selected_transformers
        )
        
        if not space_analysis['needs_id_mapping']:
            self.console.print("[green]No ID conflicts detected[/green]")
            return True
        
        # Generate ID mapping
        self.id_mapping = self.id_mapper.generate_id_mapping(
            self.redox_sites, 
            space_analysis['space_requirements'], 
            space_analysis['conflicts']
        )
        
        if self.id_mapping:
            # Save ID mapping to file
            self._save_id_mapping()
            
            # Apply ID mapping to structure
            self.structure_lines = self.id_mapper.apply_id_mapping(self.structure_lines, self.id_mapping)
            
            # CRITICAL: Update RedoxSite objects with new ID mapping
            self._update_redox_sites_with_id_mapping()
            
            # Update components in transformation parameters
            for site_id, transformer_name in self.selected_transformers.items():
                transformer_class = self.transformer_selector.registry.get_transformer(transformer_name)
                if transformer_class and site_id in space_analysis['space_requirements']:
                    req = space_analysis['space_requirements'][site_id]
                    components = req['components']
                    updated_components = transformer_class.update_components_with_id_mapping(
                        components, self.id_mapping
                    )
                    # Store updated components for transformation
                    if site_id not in self.transformation_parameters:
                        self.transformation_parameters[site_id] = {}
                    self.transformation_parameters[site_id]['_components'] = updated_components
        
        return True
    
    def _update_redox_sites_with_id_mapping(self):
        """Update RedoxSite objects to use new ID mapping"""
        if not self.id_mapping:
            return
        
        # Update RedoxSite objects with new ID mapping (silently)
        
        for site in self.redox_sites:
            # Update all atoms in the site
            for atom in site.atoms:
                mapping_key = (atom.chain, atom.resid)
                if mapping_key in self.id_mapping:
                    new_resid = self.id_mapping[mapping_key]
                    # Update atom silently
                    atom.resid = new_resid
            
            # Update all centers in the site  
            for center in site.centers:
                mapping_key = (center.chain, center.resid)
                if mapping_key in self.id_mapping:
                    new_resid = self.id_mapping[mapping_key]
                    # Update center silently
                    center.resid = new_resid
            
            # Update residue_groups dictionary
            updated_residue_groups = {}
            for (chain, resid, insertion_code), coords_list in site.residue_groups.items():
                mapping_key = (chain, resid)
                if mapping_key in self.id_mapping:
                    new_resid = self.id_mapping[mapping_key]
                    new_key = (chain, new_resid, insertion_code)
                    updated_residue_groups[new_key] = coords_list
                    # Update residue group silently
                else:
                    updated_residue_groups[(chain, resid, insertion_code)] = coords_list
            
            site.residue_groups = updated_residue_groups
            
            # Update coord_to_pdb mapping with new residue IDs
            for coords, pdb_info in site.coord_to_pdb.items():
                if 'chain' in pdb_info and 'resid' in pdb_info:
                    mapping_key = (pdb_info['chain'], pdb_info['resid'])
                    if mapping_key in self.id_mapping:
                        new_resid = self.id_mapping[mapping_key]
                        pdb_info['resid'] = new_resid
        
        if self.verbosity == VerbosityLevel.FULL:
            self.console.print(f"[green]Updated {len(self.redox_sites)} RedoxSite objects with new ID mapping[/green]")

    def _save_id_mapping(self):
        """Save ID mapping to file"""
        if not self.temp_dir or not self.id_mapping:
            return

        mapping_file = self.temp_dir / "transformation_reports" / "residue_id_mapping.txt"

        with open(mapping_file, 'w') as f:
            f.write("Redox Site Residue ID Mapping\n")
            f.write("=" * 40 + "\n\n")

            for (chain, original_id), new_id in sorted(self.id_mapping.items()):
                f.write(f"Chain {chain}: {original_id} -> {new_id}\n")

        if self.verbosity == VerbosityLevel.FULL:
            self.console.print(f"ID mapping saved to {mapping_file}")
    
    def _perform_transformations(self) -> bool:
        """Perform component matching and transformations for each site"""
        # Phase 1: Collect component matchings for all sites
        site_data = []  # List of (site, site_id, transformer_name, transformer_class, components)

        for site in self.redox_sites:
            site_id = site.site_id
            transformer_name = self.selected_transformers.get(site_id)

            if not transformer_name:
                site_data.append((site, site_id, None, None, None))
                continue

            transformer_class = self.transformer_selector.registry.get_transformer(transformer_name)
            if not transformer_class:
                site_data.append((site, site_id, transformer_name, None, None))
                continue

            # Use pre-calculated components if available, otherwise calculate
            if '_components' in self.transformation_parameters.get(site_id, {}):
                components = self.transformation_parameters[site_id]['_components']
                missing = []
            else:
                components, missing = transformer_class.match_components(site)

            if missing:
                self.console.print(f"[yellow]Missing components for {site_id}: {missing}[/yellow]")
                site_data.append((site, site_id, transformer_name, transformer_class, None))
                continue

            # If the transformer declared any ambiguous role assignments,
            # prompt the user to resolve them before the summary table is
            # built. Cache resolved components so subsequent microstate
            # passes don't re-prompt.
            if isinstance(components, dict) and components.get("_ambiguous"):
                components = self._resolve_ambiguities(site_id, components)
                self.transformation_parameters.setdefault(site_id, {})['_components'] = components

            site_data.append((site, site_id, transformer_name, transformer_class, components))

        # Phase 2: Show summary table and get sites to review
        sites_to_review = self._display_component_matching_summary_and_get_review(site_data)

        # Phase 3: Review selected sites
        sites_to_skip = set()
        if sites_to_review:
            sites_to_skip = self._review_selected_sites(site_data, sites_to_review)

        # Phase 4: Apply transformations to all non-skipped sites
        transformed_sites = []

        # Add visual separation and title before transformation output
        if self.verbosity in (VerbosityLevel.COMPACT, VerbosityLevel.MINIMAL):
            self.console.print()
            self.console.print("[bold]Applying transformations:[/bold]")

        for site, site_id, transformer_name, transformer_class, components in site_data:
            # Skip sites without transformer or components, or explicitly skipped
            if not transformer_class or components is None or site_id in sites_to_skip:
                transformed_sites.append(site)
                continue

            # Get transformation parameters
            site_params = self.transformation_parameters.get(site_id, {})
            # Remove internal components key
            site_params = {k: v for k, v in site_params.items() if not k.startswith('_')}

            # Generate and display transformation sequence
            transformations = transformer_class.get_transformation_sequence(components, site_params)
            self._display_transformation_steps(site_id, transformations)

            # Save transformation details
            self._save_transformation_details(site_id, transformations)

            # Apply transformations
            if self.verbosity == VerbosityLevel.FULL:
                self.console.print(f"\n[bold]Applying Transformations for {site_id}[/bold]")

                # Show temp directory once for debugging
                if self.temp_dir and site_id == list(self.selected_transformers.keys())[0]:  # Only for first site
                    self.console.print(f"[grey50]🔧 For debugging purposes, intermediate PDBs saved to: {self.temp_dir}[/grey50]")

            # Track initial structure
            original_site_id = self._generate_original_site_id(site)
            if self.verbosity == VerbosityLevel.FULL:
                self.console.print(f"\n🔄 Starting transformation for: {original_site_id}")

            updated_structure_lines, updated_site = self.transformation_executor.apply_transformation_sequence(
                site, transformations, self.structure_lines
            )

            self.structure_lines = updated_structure_lines
            transformed_sites.append(updated_site)

            # CRITICAL: Store the executed transformations with real atom counts for enhanced report
            # Also calculate unique atoms affected and per-atom tracking
            unique_atoms_affected = self._calculate_unique_atoms_affected(site_id, transformations)
            self.executed_transformations[site_id] = {
                'transformations': transformations,
                'unique_atoms_affected': unique_atoms_affected,
                'per_atom_tracking': self._extract_per_atom_tracking(site_id, transformations)
            }

            # Display transformation completion
            self._display_transformation_completion(site_id, transformations, original_site_id)

        # Update redox sites
        self.redox_sites = transformed_sites
        
        return True
    
    def _resolve_ambiguities(self, site_id: str, components: Dict[str, Any]) -> Dict[str, Any]:
        """Resolve transformer-declared ambiguous role assignments by prompting.

        Mutates and returns ``components``. See
        :py:meth:`RedoxSiteTransformerBase.match_components` for the contract
        of the ``"_ambiguous"`` block. After resolution the ``"_ambiguous"``
        key is removed; resolved values are written as ``<role>_id`` /
        ``<role>_chain`` for single-slot roles and ``<role>_id_<n>`` /
        ``<role>_chain_<n>`` (1-indexed) for multi-slot roles.
        """
        ambiguous = components.pop("_ambiguous", None)
        if not ambiguous:
            return components

        for block in ambiguous:
            label = block.get("label", "Role assignment")
            description = block.get("description", "")
            candidates = list(block.get("candidates") or [])
            roles = block.get("roles") or {}

            if not candidates or not roles:
                continue

            total_needed = sum(roles.values())
            if total_needed != len(candidates):
                logger.warning(
                    "Ambiguous block '%s' for site %s: roles request %d slots "
                    "but %d candidates supplied. Skipping resolution.",
                    label, site_id, total_needed, len(candidates),
                )
                continue

            self.console.print()
            self.console.print(
                f"[bold]Resolve role assignment for {site_id}: {label}[/bold]"
            )
            if description:
                self.console.print(f"[grey50]{description}[/grey50]")
            diagnostic = block.get("diagnostic")
            if diagnostic:
                # Escape Rich markup: fingerprints contain '[' from '[grey50]'-safe
                # role text is fine, but signatures/atom labels may include braces.
                from rich.markup import escape as _rich_escape
                self.console.print(f"[grey50]{_rich_escape(str(diagnostic))}[/grey50]")

            remaining = list(candidates)
            # Track per-candidate viewer labels so we can clear them
            # before each next slot's prompt (and after the whole
            # ambiguous block ends).
            _ambig_labels: list = []

            def _clear_ambig_labels():
                if not _ambig_labels:
                    return
                try:
                    from proprep.structure_prep.viewer_coordinator import (
                        viewer as _viewer,
                    )
                    for lbl in _ambig_labels:
                        _viewer.unhighlight(lbl)
                except Exception as exc:
                    logger.debug("ambig cleanup silenced: %s", exc)
                _ambig_labels.clear()

            for role_name, count in roles.items():
                for slot in range(1, count + 1):
                    if len(remaining) == 1:
                        pick = remaining.pop(0)
                    else:
                        self.console.print()
                        options_map = {}
                        for i, cand in enumerate(remaining, 1):
                            disp = cand.get("display") or (
                                f"{cand.get('chain', '?')}:{cand.get('resid', '?')}"
                            )
                            self.console.print(f"  {i}. {disp}")
                            options_map[str(i)] = disp
                        if count > 1:
                            prompt_text = (
                                f"Pick candidate for {role_name} (slot {slot}/{count})"
                            )
                        else:
                            prompt_text = f"Pick candidate for {role_name}"
                        choices = [str(i) for i in range(1, len(remaining) + 1)]
                        # Highlight each remaining candidate with the
                        # palette colour matching its 1./2./3. index in
                        # the prompt, so the user can pick by visual
                        # inspection rather than reading chain:resid
                        # strings. Old labels from the prior slot are
                        # cleared first so colours stay aligned with the
                        # current numbering.
                        _clear_ambig_labels()
                        try:
                            from proprep.structure_prep.viewer_coordinator import (
                                viewer as _viewer,
                            )
                            for i, cand in enumerate(remaining, 1):
                                cand_chain = cand.get("chain")
                                cand_resid = cand.get("resid")
                                if not cand_chain or cand_resid is None:
                                    continue
                                lbl = (
                                    f"redox_ambig_{site_id}_"
                                    f"{role_name}_slot{slot}_{i}"
                                )
                                _viewer.highlight(
                                    f":{cand_chain} and {cand_resid}",
                                    style="ball+stick",
                                    color=f"palette:{i}",
                                    label=lbl,
                                )
                                _ambig_labels.append(lbl)
                        except Exception as exc:
                            logger.debug(
                                "ambig candidate viewer hook silenced: %s", exc
                            )
                        choice = prompt_with_context(
                            processor=self.processor,
                            prompt=prompt_text,
                            choices=choices,
                            default=choices[0],
                            module="Redox Site Transformation",
                            description=(
                                f"Resolve ambiguous role '{role_name}' for site "
                                f"{site_id} ({label})"
                            ),
                            options_map=options_map,
                        )
                        pick = remaining.pop(int(choice) - 1)

                    if count == 1:
                        components[f"{role_name}_id"] = pick.get("resid")
                        components[f"{role_name}_chain"] = pick.get("chain")
                    else:
                        components[f"{role_name}_id_{slot}"] = pick.get("resid")
                        components[f"{role_name}_chain_{slot}"] = pick.get("chain")

            # Done with this ambiguous block — clear any per-candidate
            # highlights so they don't leak into the next block or back
            # to the matching-summary view.
            _clear_ambig_labels()

        return components

    def _generate_original_site_id(self, site) -> str:
        """Generate original site ID for tracking"""
        if site.centers:
            center = site.centers[0]
            return f"{center.chain}_{center.resname}_{center.resid}"
        return site.site_id

    def _display_component_matching_summary_and_get_review(self, site_data: List) -> List[str]:
        """Display summary table of all component matchings and get sites to review"""
        from rich.table import Table
        from rich.prompt import Prompt

        # Filter sites with valid components
        valid_sites = [(site, site_id, transformer_name, components)
                       for site, site_id, transformer_name, transformer_class, components in site_data
                       if components is not None]

        if not valid_sites:
            return []

        # In minimal/compact mode, skip the detailed table and prompt
        if self.verbosity in (VerbosityLevel.MINIMAL, VerbosityLevel.COMPACT):
            return []

        self.console.print("\n[bold]═══ Component Matching Summary ═══[/bold]\n")
        self.console.print(f"All {len(valid_sites)} sites detected with component matching:\n")

        # Determine all unique component keys across all sites
        all_keys = set()
        for _, _, _, components in valid_sites:
            all_keys.update(components.keys())

        # Sort keys for consistent ordering
        sorted_keys = sorted(all_keys)

        # Create summary table with dynamic columns
        # Merge _chain/_id pairs into single "Chain:Id" columns
        # e.g., center_chain + center_id -> "Center" column showing "A 81"
        merged_columns = []  # list of (display_name, key_or_pair)
        seen_keys = set()

        for key in sorted_keys:
            if key in seen_keys:
                continue
            if key.endswith('_chain'):
                base = key[:-6]  # strip '_chain'
                id_key = f"{base}_id"
                if id_key in all_keys:
                    col_name = base.replace('_', ' ').title()
                    merged_columns.append((col_name, (key, id_key)))
                    seen_keys.add(key)
                    seen_keys.add(id_key)
                else:
                    col_name = key.replace('_', ' ').title()
                    merged_columns.append((col_name, key))
                    seen_keys.add(key)
            elif key.endswith('_id'):
                base = key[:-3]  # strip '_id'
                chain_key = f"{base}_chain"
                if chain_key in all_keys:
                    # Will be handled when we encounter the _chain key
                    continue
                else:
                    col_name = key.replace('_', ' ').title()
                    merged_columns.append((col_name, key))
                    seen_keys.add(key)
            else:
                col_name = key.replace('_', ' ').title()
                merged_columns.append((col_name, key))
                seen_keys.add(key)

        table = Table(title="Redox Site Component Matching Summary")
        table.add_column("Site ID", style="cyan")
        table.add_column("Type", style="yellow")

        for col_name, _ in merged_columns:
            table.add_column(col_name, style="green")

        # Add rows
        for site, site_id, transformer_name, components in valid_sites:
            row_data = [site_id, transformer_name]

            # Special handling for transformers with no components
            if not components:
                for _ in merged_columns:
                    row_data.append("N/A")
            else:
                for _, key_or_pair in merged_columns:
                    if isinstance(key_or_pair, tuple):
                        chain_key, id_key = key_or_pair
                        chain = components.get(chain_key, '—')
                        res_id = components.get(id_key, '—')
                        row_data.append(f"{chain} {res_id}")
                    else:
                        value = components.get(key_or_pair, '—')
                        row_data.append(str(value))

            table.add_row(*row_data)

        self.console.print(table)
        self.console.print()

        # Highlight matched components across all sites coloured by
        # role — every "Cys_SG" the same colour, every "center" the
        # same colour, etc. This makes role-mismatches across a
        # multi-site system (e.g. one heme's "histidine_proximal"
        # accidentally pointing at the wrong residue) visually
        # obvious. Up to 12 distinct roles via the Paired-12 palette.
        # Labels keyed on the role base name so the per-site review
        # hook can re-use the same colour scheme.
        self._matching_role_labels = self._highlight_matched_components_by_role(valid_sites)

        # Prompt for sites to review
        response = prompt_with_context(
            processor=self.processor,
            prompt="Review specific sites?",
            default="no",
            module="Redox Site Transformation",
            description="Review specific sites"
        )
        response = response.strip().lower()

        if response in ("no", "n", ""):
            self.console.print("[green]✓ All component matchings accepted[/green]\n")
            self._clear_matching_role_highlights()
            return []
        elif response == "all":
            return [site_id for _, site_id, _, _ in valid_sites]
        else:
            # Parse comma-separated site numbers or IDs
            try:
                # Support both "1,2,3" and "site_1,site_2,site_3"
                site_refs = [s.strip() for s in response.split(',')]
                sites_to_review = []

                for ref in site_refs:
                    if ref.startswith('site_'):
                        sites_to_review.append(ref)
                    else:
                        # Assume it's a number, convert to site_N
                        sites_to_review.append(f"site_{ref}")

                return sites_to_review
            except Exception as e:
                self.console.print(f"[yellow]Could not parse site selection: {e}[/yellow]")
                return []

    def _highlight_matched_components_by_role(
        self, valid_sites: List
    ) -> List[str]:
        """Highlight every matched component across all sites with one
        palette colour per role. Returns the list of viewer labels
        applied so the caller can later clear them.
        """
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer
        except Exception as exc:
            logger.debug("matching summary viewer hook silenced: %s", exc)
            return []

        # Collect every (chain, resid) per role base across all sites.
        # Walk only the chain-bearing keys; pair each with its
        # corresponding id key. This avoids double-counting.
        role_residues: Dict[str, set] = {}
        for _, _site_id, _xform, components in valid_sites:
            if not components:
                continue
            for key, chain_val in list(components.items()):
                if key.endswith("_chain"):
                    base = key[: -len("_chain")]
                    id_key = f"{base}_id"
                elif "_chain_" in key and key.rsplit("_", 1)[-1].isdigit():
                    slot = key.rsplit("_", 1)[-1]
                    base = key[: -(len("_chain_") + len(slot))]
                    id_key = f"{base}_id_{slot}"
                else:
                    continue
                resid_val = components.get(id_key)
                if not chain_val or resid_val is None:
                    continue
                role_residues.setdefault(base, set()).add((str(chain_val), resid_val))

        if not role_residues:
            return []

        applied_labels: List[str] = []
        # Sort roles for stable palette assignment across reruns.
        for idx, role in enumerate(sorted(role_residues.keys()), 1):
            pairs = role_residues[role]
            if not pairs:
                continue
            clauses = [f"(:{c} and {r})" for c, r in sorted(pairs)]
            label = f"redox_role_{role}"
            try:
                _viewer.highlight(
                    " or ".join(clauses),
                    style="ball+stick",
                    color=f"palette:{idx}",
                    label=label,
                )
                applied_labels.append(label)
            except Exception as exc:
                logger.debug("role highlight silenced (%s): %s", role, exc)
        return applied_labels

    def _clear_matching_role_highlights(self) -> None:
        """Remove all per-role highlights applied by the matching summary
        view. Safe to call when no highlights are present.
        """
        labels = getattr(self, "_matching_role_labels", None) or []
        if not labels:
            return
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer
            for lbl in labels:
                _viewer.unhighlight(lbl)
        except Exception as exc:
            logger.debug("matching role cleanup silenced: %s", exc)
        self._matching_role_labels = []

    def _review_selected_sites(self, site_data: List, sites_to_review: List[str]) -> set:
        """Review selected sites and return set of sites to skip"""
        from rich.prompt import Confirm, Prompt

        sites_to_skip = set()
        site_dict = {site_id: (site, site_id, transformer_name, transformer_class, components)
                     for site, site_id, transformer_name, transformer_class, components in site_data}

        self.console.print("\n[bold]Reviewing selected sites...[/bold]\n")

        for site_id in sites_to_review:
            if site_id not in site_dict:
                self.console.print(f"[yellow]Site {site_id} not found, skipping[/yellow]")
                continue

            site, _, transformer_name, transformer_class, components = site_dict[site_id]

            if components is None:
                self.console.print(f"[yellow]Site {site_id} has no valid components, skipping[/yellow]")
                sites_to_skip.add(site_id)
                continue

            # Display detailed component matching
            self._display_component_matching(site_id, transformer_name, components)

            # Focus the viewer on this site's matched components. The
            # per-role colours from the matching summary stay visible
            # in the background, so the user can see this site
            # zoomed-in while still spotting the matched components on
            # other sites for context. Replaces the prior site's
            # focus on each iteration.
            self._focus_review_on_site(site_id, components)

            # Confirm component matching
            if not confirm_with_context(
                processor=self.processor,
                prompt="Is this component matching correct?",
                default=True,
                module="Redox Site Transformation",
                description="Confirm component matching"
            ):
                self.console.print(f"[yellow]⚠ Component matching rejected for {site_id}[/yellow]\n")

                # Ask what to do
                self.console.print("Options:")
                self.console.print("  1. Skip this site (continue with other sites)")
                self.console.print("  2. Cancel all transformations")

                choice = prompt_with_context(
                    processor=self.processor,
                    prompt="Choose option",
                    choices=["1", "2"],
                    default="1",
                    module="Redox Site Transformation",
                    description="Choose action after rejection",
                    options_map={
                        "1": "Skip this site",
                        "2": "Cancel all transformations"
                    }
                )

                if choice == "1":
                    sites_to_skip.add(site_id)
                    self.console.print(f"[yellow]Skipping {site_id}[/yellow]\n")
                elif choice == "2":
                    self.console.print("[red]Cancelling all transformations[/red]")
                    # Return all sites as skipped
                    return set(s[1] for s in site_data if s[4] is not None)

        if sites_to_skip:
            num_proceeding = len([s for s in site_data if s[4] is not None]) - len(sites_to_skip)
            self.console.print(f"[green]✓ Proceeding with {num_proceeding} sites (skipping {len(sites_to_skip)})[/green]\n")
        else:
            self.console.print("[green]✓ All reviewed component matchings confirmed[/green]\n")

        # Done with the review loop — clear both the per-review focus
        # and the per-role summary highlights so the next stage starts
        # with a clean viewer.
        self._clear_review_focus()
        self._clear_matching_role_highlights()

        return sites_to_skip

    def _focus_review_on_site(self, site_id: str, components: Dict[str, Any]) -> None:
        """Camera-focus the viewer on a single site's matched components.

        The per-role highlights from the matching summary stay visible;
        this hook only adds a focused-mode rep that auto-centres the
        camera on this site's atoms. Single shared label so re-firing
        per site replaces the prior focus.
        """
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer
        except Exception as exc:
            logger.debug("review focus viewer hook silenced: %s", exc)
            return

        pairs = set()
        for key, chain_val in components.items():
            if key.endswith("_chain"):
                base = key[: -len("_chain")]
                resid_val = components.get(f"{base}_id")
            elif "_chain_" in key and key.rsplit("_", 1)[-1].isdigit():
                slot = key.rsplit("_", 1)[-1]
                base = key[: -(len("_chain_") + len(slot))]
                resid_val = components.get(f"{base}_id_{slot}")
            else:
                continue
            if not chain_val or resid_val is None:
                continue
            pairs.add((str(chain_val), resid_val))

        self._clear_review_focus()
        if not pairs:
            return
        clauses = [f"(:{c} and {r})" for c, r in sorted(pairs)]
        try:
            _viewer.highlight(
                " or ".join(clauses),
                style="ball+stick",
                color="#ffff00",
                label=self._review_focus_label(site_id),
                focused=True,
            )
            self._active_review_label = self._review_focus_label(site_id)
        except Exception as exc:
            logger.debug("review focus highlight silenced: %s", exc)

    @staticmethod
    def _review_focus_label(site_id: str) -> str:
        return f"redox_review_{site_id}"

    def _clear_review_focus(self) -> None:
        """Remove the active per-site review focus highlight."""
        label = getattr(self, "_active_review_label", None)
        if not label:
            return
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer
            _viewer.unhighlight(label)
        except Exception as exc:
            logger.debug("review focus cleanup silenced: %s", exc)
        self._active_review_label = None

    def _display_component_matching(self, site_id: str, transformer_name: str, components: Dict[str, Any]):
        """Display component matching results"""
        panel_content = f"Matched Components:\n"
        for key, value in components.items():
            panel_content += f"  ✓ {key}: {value}\n"
        
        panel = Panel(
            panel_content.strip(),
            title=f"Component Matching for {site_id} ({transformer_name})",
            border_style="blue"
        )
        self.console.print(panel)
    
    def _display_transformation_steps(self, site_id: str, transformations: List[Dict[str, Any]]):
        """Display transformation steps in narrative format"""
        if self.verbosity != VerbosityLevel.FULL:
            return  # Only show in full verbosity mode

        from rich.panel import Panel

        # Build narrative content
        content = []
        content.append(f"[bold cyan]Transformation Plan for {site_id}[/bold cyan]\n")

        for i, transform in enumerate(transformations, 1):
            selector = transform["selector"]
            action = transform["action"]
            description = transform["description"]

            # Build a readable one-line summary
            # Extract key info from selector
            chain = selector.get("chain_id", "")
            resname = selector.get("residue_name", "")
            resid = selector.get("residue_id", "")
            atom_names = selector.get("atom_names", [])

            # Build selector string
            if resname and resid:
                selector_str = f"{chain}:{resname} {resid}"
                if atom_names:
                    if len(atom_names) <= 3:
                        selector_str += f" [{', '.join(atom_names)}]"
                    else:
                        selector_str += f" [{len(atom_names)} atoms]"
            elif resid:
                selector_str = f"{chain}:{resid}"
            else:
                # Fallback to showing all selector items
                selector_str = ", ".join(f"{k}={v}" for k, v in selector.items() if v)

            # Build action string
            action_parts = []
            if "change_residue_name" in action:
                action_parts.append(f"→ {action['change_residue_name']}")
            if "change_residue_id" in action and action["change_residue_id"] != resid:
                action_parts.append(f"→ {action['change_residue_id']}")
            if "rename_atoms" in action:
                renames = action["rename_atoms"]
                if len(renames) <= 2:
                    rename_str = ", ".join(f"{old}→{new}" for old, new in renames.items())
                    action_parts.append(f"[{rename_str}]")
                else:
                    action_parts.append(f"[rename {len(renames)} atoms]")
            if "convert_to_hetatm" in action and action["convert_to_hetatm"]:
                action_parts.append("(HETATM)")

            action_str = " ".join(action_parts) if action_parts else "modify"

            # Format the line
            content.append(f"[cyan]{i:2d}.[/cyan] {description}")
            content.append(f"    [yellow]{selector_str}[/yellow] [green]{action_str}[/green]\n")

        content.append(f"[grey50]Total: {len(transformations)} transformation operations[/grey50]")

        # Display in a panel
        panel = Panel(
            "\n".join(content),
            border_style="cyan",
            padding=(1, 2),
            expand=False
        )
        self.console.print(panel)
    
    def _save_transformation_details(self, site_id: str, transformations: List[Dict[str, Any]]):
        """Save transformation details to file"""
        if not self.temp_dir:
            return

        details_file = self.temp_dir / "transformation_reports" / f"transformations_{site_id}.json"

        with open(details_file, 'w') as f:
            json.dump(transformations, f, indent=2)

        if self.verbosity == VerbosityLevel.FULL:
            self.console.print(f"Transformations saved to {details_file}")
    
    def _display_transformation_completion(self, site_id: str, transformations: List[Dict[str, Any]], original_site_id: str):
        """Display atom-level transformation tracking for the site"""
        # Get per-atom tracking data if available
        executed_data = self.executed_transformations.get(site_id, {})
        if isinstance(executed_data, dict):
            per_atom_tracking = executed_data.get('per_atom_tracking', {})
        else:
            per_atom_tracking = {}

        # Analyze the transformations to build atom movement tracking
        atom_transformations = self._analyze_atom_transformations_with_per_atom_data(site_id, transformations, per_atom_tracking)

        # Calculate summary stats
        total_atoms_tracked = sum(t.get('total_atoms', t.get('atom_count', 0)) for t in atom_transformations)
        transformation_count = len(atom_transformations)
        total_lines_modified = sum(t.get('lines_modified', 0) for t in transformations)

        # Handle different verbosity levels
        if self.verbosity == VerbosityLevel.MINIMAL:
            # Minimal: One-line summary with metrics
            self.console.print(f"  ✓ {original_site_id}: {total_atoms_tracked} atoms transformed across {transformation_count} residues | {total_lines_modified} PDB lines modified")
            return
        elif self.verbosity == VerbosityLevel.COMPACT:
            # Compact: Just show metrics line, skip the panel
            self.console.print(f"\n[bold]✓ {original_site_id}[/bold]")
            self.console.print(f"[bold]Metrics:[/bold] {total_atoms_tracked} atoms transformed across {transformation_count} residues | {total_lines_modified} PDB lines modified in {len(transformations)} operations\n")
            return

        # Full verbosity: Build the complete display format
        panel_content = f"\nAtom Transformations:\n\n"

        for transformation in atom_transformations:
            if transformation['type'] == 'split':
                original_res = transformation['original']
                splits = transformation['splits']
                total_atoms = transformation['total_atoms']

                panel_content += f"{original_res['chain']}:{original_res['resname']} {original_res['resid']} ({total_atoms} atoms) -> SPLIT:\n"

                for i, split in enumerate(splits):
                    prefix = "├─" if i < len(splits) - 1 else "└─"
                    atom_detail = split.get('atom_detail', f"{split['atom_count']} atoms")
                    line_prefix = f"  {prefix} {split['chain']}:{split['resname']} {split['resid']}: "
                    # Wrap long atom lists with indentation aligned to content
                    continuation_indent = " " * len(line_prefix)
                    # Use terminal width minus panel border overhead (4 chars for borders + padding)
                    try:
                        max_width = self.console.width - 4
                    except Exception:
                        max_width = 120
                    content_width = max_width - len(line_prefix)
                    if content_width > 20 and len(atom_detail) > content_width:
                        # Split atom detail at commas to wrap neatly
                        atoms = atom_detail.split(", ")
                        lines = []
                        current_line = ""
                        for atom in atoms:
                            candidate = f"{current_line}, {atom}" if current_line else atom
                            if current_line and len(candidate) > content_width:
                                lines.append(current_line + ",")
                                current_line = atom
                            else:
                                current_line = candidate
                        if current_line:
                            lines.append(current_line)
                        atom_detail = ("\n" + continuation_indent).join(lines)
                    panel_content += f"{line_prefix}{atom_detail}\n"

                panel_content += "\n"

            elif transformation['type'] == 'rename':
                original_res = transformation['original']
                final_res = transformation['final']
                atom_count = transformation['atom_count']

                panel_content += f"{original_res['chain']}:{original_res['resname']} {original_res['resid']} ({atom_count} atoms) -> "
                panel_content += f"{final_res['chain']}:{final_res['resname']} {final_res['resid']} ({atom_count} atoms renamed)\n\n"

            elif transformation['type'] == 'unchanged':
                res = transformation.get('residue') or transformation.get('original')
                atom_count = transformation['atom_count']

                panel_content += f"{res['chain']}:{res['resname']} {res['resid']} ({atom_count} atoms) -> UNCHANGED\n\n"

        panel_content += f"Summary: {total_atoms_tracked} atoms tracked through {transformation_count} residue transformations"

        panel = Panel(
            panel_content.strip(),
            title=f"✅ Transformation Complete: {original_site_id}",
            border_style="green",
            expand=False
        )
        self.console.print(panel)

        # Show a simple metrics summary with clear labels
        self.console.print(f"\n[bold]Metrics:[/bold] {total_atoms_tracked} atoms transformed across {transformation_count} residues | {total_lines_modified} PDB lines modified in {len(transformations)} operations")

        # Save detailed log
        self._save_transformation_log(site_id, transformations, original_site_id)
    
    def _save_transformation_log(self, site_id: str, transformations: List[Dict[str, Any]], original_site_id: str):
        """Save detailed transformation log"""
        if not self.temp_dir:
            return
        
        log_dir = self.temp_dir / "transformation_reports" / "transformation_details" / f"site_{site_id}"
        log_dir.mkdir(parents=True, exist_ok=True)
        
        log_file = log_dir / "transformation_log.txt"
        
        with open(log_file, 'w') as f:
            f.write(f"Transformation Log for Site {site_id}\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Original Site ID: {original_site_id}\n\n")
            
            for i, transform in enumerate(transformations, 1):
                f.write(f"Step {i}: {transform['description']}\n")
                f.write(f"  Lines modified: {transform.get('lines_modified', 0)}\n")
                f.write(f"  Selector: {transform['selector']}\n")
                f.write(f"  Action: {transform['action']}\n\n")
        
        self.console.print(f"  Detailed log saved to\n{log_file}")
    
    def _finalize_transformations(self) -> bool:
        """Finalize transformations and generate outputs"""
        # Update redox site definitions
        self._update_redox_site_definitions()

        # Process PDB for final output (reorder, renumber, update sites)
        self._process_final_pdb()
        
        # Save transformed structure to working directory
        import os
        output_file = Path(os.getcwd()) / "final_transformed_structure.pdb"
        
        with open(output_file, 'w') as f:
            f.writelines(self.structure_lines)
        
        # Structure saving message will be in the completion summary
        
        # Create mapping report (also saved to working directory)
        self._create_residue_mapping_report()
        
        # Display completion message
        self._display_completion_message(output_file)
        
        # Update workspace with transformed structure
        if self.processor:
            workspace = self.processor.workspace
            workspace.set("detected_redox_sites", self.redox_sites)

            # Store the transformed structure with the correct workspace key
            # This ensures it will be picked up by tLEaP and other downstream processes
            workspace.set("transformed_pdb_file", str(output_file))

            # Also update the generic transformed_structure key for backward
            # compatibility. NOTE: this stores a PATH string (not a Structure
            # object) because the modAA parameterizer reads it as a path. That
            # breaks the usual *_structure = Structure contract, so
            # StructureSelector.get_structure_object skips string-valued
            # *_structure keys rather than iterating them and crashing on .id.
            workspace.set("transformed_structure", str(output_file))

            self.console.print(f"[green]✅ Workspace updated with transformed_pdb_file: {output_file.name}[/green]")

            # Export updated RedoxSites to JSON with renumbered residues
            self._export_updated_redox_sites(output_file)

        # Validate workspace consistency
        self._validate_workspace_consistency()

        return True
    
    def batch_transformation(self) -> bool:
        """
        Batch transformation mode - generate redox microstates

        Multi-step workflow:
        1. Load structure and sites
        2. Setup transformer for each site
        3. Select redox/spin state combinations for each site
        4. View available microstates
        5. Select specific microstates or generate all
        6. Configure propionate/pH treatment for the states actually used
        7. Generate selected microstates
        """
        self.console.print("\n[bold cyan]🔄 Batch Redox Microstate Generation[/bold cyan]")

        # Step 1: Load structure and sites
        if not self._load_redox_sites_and_structure():
            return False

        self.console.print(f"\n[bold]Found {len(self.redox_sites)} redox sites[/bold]")

        # Step 2: Setup transformer for each site (interactive)
        if not self._setup_site_transformers():
            return False

        # Step 3: Select redox/spin state combinations for each site
        if not self._select_state_combinations():
            return False

        # Step 4: Generate all microstate combinations from selected states
        microstates = self._generate_microstate_combinations()

        if not microstates:
            self.console.print("[red]No microstates could be generated[/red]")
            return False

        # Step 5: Display available microstates
        self._display_microstates(microstates)

        # Step 6: Let user select which microstates to generate
        selected_microstates = self._select_microstates_to_generate(microstates)

        if not selected_microstates:
            self.console.print("[yellow]No microstates selected for generation[/yellow]")
            return False

        # Step 7: Propionate / pH treatment, bound per redox state. Runs AFTER
        # selection so it only prompts for the (site, state) pairs the chosen
        # microstates actually use. Mirrors the single-microstate path's prompt
        # so batch no longer silently defaults the propionate protonation.
        if not self._configure_state_protonation(selected_microstates):
            return False

        # Record per-(site,state) transformer info (incl. the Step-6 ph_treatment
        # and its resolved forcefield set) in the workspace so the Topology
        # Generator pre-selects / filters the FF-set picker to the propionate
        # treatment chosen here instead of offering every set.
        self._store_batch_transformer_info_in_workspace(selected_microstates)

        # Step 8: Generate the selected microstates
        return self._generate_selected_microstates(selected_microstates)
    
    def _setup_site_transformers(self) -> bool:
        """Setup transformer for each redox site using bulk assignment (Step 2)"""
        self.console.print("\n[bold]Step 2: Setup Transformers for Each Site[/bold]")

        self.site_transformer_assignments = {}

        # Build compatibility map for all sites
        site_compatibility = {}  # site_num -> {site_id, compatible_transformers}
        all_transformer_names = set()

        for i, site in enumerate(self.redox_sites):
            site_num = i + 1
            site_id = site.site_id if hasattr(site, 'site_id') else f'site_{i}'

            # Get compatible transformers for this site
            compatible = self.transformer_selector.get_compatible_transformers(site)

            if not compatible:
                self.console.print(f"[red]No compatible transformers found for {site_id}[/red]")
                return False

            site_compatibility[site_num] = {
                'site_id': site_id,
                'compatible': compatible
            }

            for t in compatible:
                all_transformer_names.add(t['name'])

        # Check for pre-existing transformer mappings from detection phase
        workspace = self.processor._get_workspace()
        transformer_mappings = workspace.get('redox_transformer_mappings', {})

        if transformer_mappings:
            self.console.print("[cyan]Found transformer mappings from detection phase. Attempting automatic assignment...[/cyan]")

            auto_assignments = {}
            failed = []

            for site_num in sorted(site_compatibility.keys()):
                info = site_compatibility[site_num]
                site_id = info['site_id']
                site = self.redox_sites[site_num - 1]
                site_type = getattr(site, 'site_type', None)
                compatible_names = [t['name'] for t in info['compatible']]

                if site_type and site_type in transformer_mappings:
                    suggested = transformer_mappings[site_type]
                    if suggested in compatible_names:
                        auto_assignments[site_id] = suggested
                        self.console.print(f"  Site {site_num} ({site_type}): {suggested} ✓")
                    else:
                        failed.append((site_num, site_type, suggested))
                        self.console.print(f"  Site {site_num} ({site_type}): {suggested} [red]✗ Not compatible[/red]")
                else:
                    failed.append((site_num, site_type or "unknown", "no mapping"))
                    self.console.print(f"  Site {site_num} ({site_type or 'unknown'}): [yellow]No mapping found[/yellow]")

            if not failed:
                self.site_transformer_assignments = auto_assignments
                self.console.print(f"\n[green]✓ Auto-assigned transformers for all {len(auto_assignments)} sites[/green]")
                return True
            else:
                self.console.print(f"\n[yellow]{len(failed)} site(s) could not be auto-assigned. Falling back to manual assignment.[/yellow]")

        # Display summary table of all sites and their compatible transformers
        self.console.print(f"\n[cyan]Found {len(self.redox_sites)} sites requiring transformer assignment[/cyan]\n")

        table = Table(title="Site Compatibility Overview")
        table.add_column("Site #", style="cyan", justify="right")
        table.add_column("Site ID", style="yellow")
        table.add_column("Compatible Transformers", style="green")

        for site_num in sorted(site_compatibility.keys()):
            info = site_compatibility[site_num]
            compatible_str = ", ".join([f"{t['name']} ({t['confidence']:.2f})" for t in info['compatible']])
            table.add_row(str(site_num), info['site_id'], compatible_str)

        self.console.print(table)

        # Show available transformers
        self.console.print(f"\n[bold]Available transformers:[/bold] {', '.join(sorted(all_transformer_names))}")

        # Generate smart default suggestion
        # For each site, prefer transformer with 1.0 confidence that is NOT no_transformation
        best_transformers = []
        for site_num in sorted(site_compatibility.keys()):
            compatible = site_compatibility[site_num]['compatible']

            # Find all transformers with 1.0 confidence
            perfect_matches = [t for t in compatible if t['confidence'] == 1.0]

            if perfect_matches:
                # Prefer non-no_transformation if available
                non_no_transform = [t for t in perfect_matches if t['name'] != 'no_transformation']
                if non_no_transform:
                    best_transformers.append(non_no_transform[0]['name'])
                else:
                    best_transformers.append(perfect_matches[0]['name'])
            else:
                # No perfect match, use highest confidence
                best_transformers.append(compatible[0]['name'])

        # Build default suggestion
        if len(set(best_transformers)) == 1:
            # All sites have same best transformer
            default_suggestion = f'1-{len(self.redox_sites)}:"{best_transformers[0]}"'
        else:
            # Mixed - group by transformer
            from collections import Counter
            transformer_counts = Counter(best_transformers)

            # Build groups
            groups = {}
            for site_num, transformer in enumerate(best_transformers, 1):
                if transformer not in groups:
                    groups[transformer] = []
                groups[transformer].append(site_num)

            # Format as ranges
            group_strs = []
            for transformer in sorted(groups.keys()):
                sites = groups[transformer]
                # Convert to ranges
                ranges = []
                start = sites[0]
                end = sites[0]
                for site_num in sites[1:]:
                    if site_num == end + 1:
                        end = site_num
                    else:
                        if start == end:
                            ranges.append(str(start))
                        else:
                            ranges.append(f"{start}-{end}")
                        start = end = site_num
                # Add final range
                if start == end:
                    ranges.append(str(start))
                else:
                    ranges.append(f"{start}-{end}")

                group_strs.append(f'{",".join(ranges)}:"{transformer}"')

            default_suggestion = ", ".join(group_strs)

        # Prompt for bulk assignment
        self.console.print("\n[bold]Bulk Transformer Assignment[/bold]")
        self.console.print("Enter assignments using the format:")
        self.console.print('  [cyan]<site-range>:"<transformer>"[/cyan]')
        self.console.print("\nExamples:")
        self.console.print('  1-61:"heme_bis_his_c_type"', highlight=False)
        self.console.print('  1,2,3:"no_transformation"', highlight=False)
        self.console.print('  1-3:"heme_bis_his_c_type", 4-6:"no_transformation"', highlight=False)
        self.console.print(f"\n[yellow]Suggested:[/yellow] {default_suggestion}", highlight=False)

        while True:
            assignment_str = prompt_with_context(
                processor=self.processor,
                prompt="\nEnter assignment",
                default=default_suggestion,
                module="Redox Site Transformation",
                description="Bulk transformer assignment"
            )

            try:
                # Parse the assignment string
                assignments = parse_bulk_assignment(assignment_str, len(self.redox_sites))

                # Validate that all sites are assigned
                if len(assignments) != len(self.redox_sites):
                    missing = set(range(1, len(self.redox_sites) + 1)) - set(assignments.keys())
                    self.console.print(f"[red]Error: Not all sites assigned. Missing sites: {sorted(missing)}[/red]")
                    continue

                # Validate transformer names and compatibility
                validation_errors = []
                for site_num, transformer_name in assignments.items():
                    site_info = site_compatibility[site_num]
                    compatible_names = [t['name'] for t in site_info['compatible']]

                    if transformer_name not in compatible_names:
                        validation_errors.append(
                            f"Site {site_num} ({site_info['site_id']}): '{transformer_name}' not compatible. "
                            f"Options: {', '.join(compatible_names)}"
                        )

                if validation_errors:
                    self.console.print("[red]Validation errors:[/red]")
                    for error in validation_errors:
                        self.console.print(f"  [red]• {error}[/red]")
                    continue

                # Show confirmation table
                self.console.print("\n[bold green]Assignment Preview:[/bold green]")
                preview_table = Table()
                preview_table.add_column("Site #", style="cyan")
                preview_table.add_column("Site ID", style="yellow")
                preview_table.add_column("Assigned Transformer", style="green")

                for site_num in sorted(assignments.keys()):
                    site_id = site_compatibility[site_num]['site_id']
                    transformer = assignments[site_num]
                    preview_table.add_row(str(site_num), site_id, transformer)

                self.console.print(preview_table)

                # Confirm
                if confirm_with_context(
                    processor=self.processor,
                    prompt="Accept these assignments?",
                    default=True,
                    module="Redox Site Transformation",
                    description="Confirm bulk transformer assignment"
                ):
                    # Apply assignments
                    for site_num, transformer_name in assignments.items():
                        site_id = site_compatibility[site_num]['site_id']
                        self.site_transformer_assignments[site_id] = transformer_name

                    self.console.print(f"[bold green]✓ Assigned transformers to {len(self.redox_sites)} sites[/bold green]")
                    return True
                else:
                    self.console.print("[yellow]Retrying assignment...[/yellow]")

            except ValueError as e:
                self.console.print(f"[red]Parse error: {e}[/red]")
                self.console.print("[yellow]Please try again with correct format[/yellow]")
                continue
            except Exception as e:
                self.console.print(f"[red]Unexpected error: {e}[/red]")
                logger.error(f"Bulk assignment error: {e}")
                continue

        return False  # Should never reach here

    def _available_states_for_transformer(self, transformer_name: str) -> List[Dict[str, Any]]:
        """Enumerate the (redox_state, spin_state) options a transformer offers.

        Single source of truth shared by the protonation step (2b) and the
        state-selection step (3), so the two can never disagree on which states
        a site has. Reads the transformer's forcefield metadata; falls back to a
        reduced/oxidized low-spin pair when metadata is unavailable, and returns
        a single 'unchanged' state for the no-op transformer.
        """
        from proprep.forcefield_params import load_forcefield_metadata
        from proprep.redoxsite_prep.transformation.redox_transformer_framework import redox_transformer_registry

        transformer_class = redox_transformer_registry.get_transformer(transformer_name)

        if transformer_name == 'no_transformation':
            return [
                {'redox_state': 'unchanged', 'spin_state': 'unchanged', 'is_valid': True, 'is_default': True}
            ]
        if not transformer_class or not getattr(transformer_class, 'FORCEFIELD_PATH', None):
            return [
                {'redox_state': 'reduced', 'spin_state': 'low_spin', 'is_valid': True},
                {'redox_state': 'oxidized', 'spin_state': 'low_spin', 'is_valid': True}
            ]
        try:
            metadata = load_forcefield_metadata(transformer_class.FORCEFIELD_PATH)
            available_states = []
            for redox_state, redox_data in metadata['redox_states'].items():
                for spin_state, spin_data in redox_data['spin_states'].items():
                    available_states.append({
                        'redox_state': redox_state,
                        'spin_state': spin_state,
                        'is_valid': spin_data.get('is_valid', True),
                        'is_default': spin_data.get('is_default', False)
                    })
            return available_states
        except Exception as e:
            logger.warning(f"Failed to load forcefield metadata for {transformer_name}: {e}")
            return [
                {'redox_state': 'reduced', 'spin_state': 'low_spin', 'is_valid': True},
                {'redox_state': 'oxidized', 'spin_state': 'low_spin', 'is_valid': True}
            ]

    def _configure_state_protonation(self, microstates: List[Dict]) -> bool:
        """Step 6: configure propionate / pH treatment per (site, redox state).

        Model A ("bind to each redox state"): the propionate protonation choice
        is attached to each (site, redox/spin state), so a site's reduced and
        oxidized forms may carry different protonation, and every microstate
        inherits the treatment of the states it is built from.

        Runs after microstate selection and is scoped to the (site, state) pairs
        the ``microstates`` actually use, so it never prompts for a state the
        chosen microstates don't contain. Only transformers whose metadata
        declares a constant_pH/fixed_pH fork (i.e. expose a ``ph_treatment``
        parameter) are prompted; cofactors with a single declared treatment are
        left untouched. Reuses the single-path interactive configurator so the
        prompt, per-option help, per-ring gating (protonation_<role> valid only
        under fixed_pH) and batch site-syntax are identical to option 2.
        """
        from proprep.redoxsite_prep.transformation.redox_transformer_framework import redox_transformer_registry

        self.site_state_protonation = {}

        # (site_id -> display number), stable across the flow.
        site_number: Dict[str, int] = {}
        for i, site in enumerate(self.redox_sites):
            sid = site.site_id if hasattr(site, 'site_id') else f'site_{i}'
            site_number[sid] = i + 1

        # Exactly which (site, redox/spin) pairs the chosen microstates use.
        used: Dict[str, set] = {}
        for ms in microstates:
            for sid, st in ms.items():
                used.setdefault(sid, set()).add((st['redox_state'], st['spin_state']))

        # Resolve each transformer once: (tclass, param_defs, ph_param, prefix)
        # or None when it exposes no pH fork.
        fork_info: Dict[str, Optional[Tuple]] = {}

        def _fork(tname):
            if tname not in fork_info:
                tclass = redox_transformer_registry.get_transformer(tname)
                pdefs = tclass.get_parameter_definitions() if tclass else {}
                ph_param = getattr(tclass, 'PH_TREATMENT_PARAM', 'ph_treatment') if tclass else None
                if tclass and ph_param in pdefs:
                    fork_info[tname] = (
                        tclass, pdefs, ph_param,
                        getattr(tclass, 'PROTONATION_PARAM_PREFIX', 'protonation_'),
                    )
                else:
                    fork_info[tname] = None
            return fork_info[tname]

        # transformer -> {(redox, spin) -> [site_ids in site order]}, restricted
        # to fork transformers and the states actually used.
        group_state_sites: Dict[str, Dict[Tuple[str, str], List[str]]] = {}
        for i, site in enumerate(self.redox_sites):
            sid = site.site_id if hasattr(site, 'site_id') else f'site_{i}'
            if sid not in used:
                continue
            tname = self.site_transformer_assignments.get(sid)
            if not tname or _fork(tname) is None:
                continue
            for state in used[sid]:
                group_state_sites.setdefault(tname, {}).setdefault(state, [])
                if sid not in group_state_sites[tname][state]:
                    group_state_sites[tname][state].append(sid)

        if not group_state_sites:
            return True  # nothing titratable in the selected microstates; silent

        self.console.print("\n[bold]Step 6: Propionate / pH Treatment (per redox state)[/bold]")
        self.console.print(
            "[grey50]Protonation is bound to each redox state, so a site's reduced and "
            "oxidized forms may differ. Each selected microstate inherits the treatment "
            "of the states it is built from. Only the states your selection uses are "
            "shown; leave a prompt empty to take the default (constant_pH / titratable "
            "PRN).[/grey50]"
        )

        # Scratch store used by the reused configurator; cleared at the end so it
        # never leaks into microstate generation (which builds its own params).
        self.transformation_parameters = {}

        for tname in sorted(group_state_sites):
            tclass, param_defs, ph_param, proto_prefix = _fork(tname)
            # ph_treatment first (it gates the per-ring choices), then the roles.
            ordered_params = [ph_param] + [n for n in param_defs if n.startswith(proto_prefix)]
            state_map = group_state_sites[tname]

            all_nums = sorted(site_number[s] for sids in state_map.values() for s in sids)
            self.console.print(
                f"\n[Transformer: {tname}] Sites: "
                f"{', '.join(map(str, sorted(set(all_nums))))}"
            )

            for redox, spin in sorted(state_map):
                site_ids = state_map[(redox, spin)]
                site_numbers = [site_number[sid] for sid in site_ids]
                self.console.print(
                    f"\n[bold cyan]State: {redox}/{spin}[/bold cyan] "
                    f"[grey50](sites {', '.join(map(str, site_numbers))})[/grey50]"
                )

                # Seed the redox/spin context so get_valid_options resolves the
                # right treatments (and gates protonation_<role> once ph_treatment
                # is set) for each site.
                for sid in site_ids:
                    self.transformation_parameters[sid] = {
                        'redox_state': redox, 'spin_state': spin
                    }

                for pname in ordered_params:
                    pdef = param_defs.get(pname)
                    if not pdef or pdef.get('type') != 'choice':
                        continue
                    self._configure_choice_parameter_interactive(
                        pname, pdef, site_ids, site_numbers, tname
                    )

                # Snapshot protonation-only keys for this (site, state).
                state_key = f"{redox}_{spin}"
                for sid in site_ids:
                    params = self.transformation_parameters.get(sid, {})
                    self.site_state_protonation.setdefault(sid, {})[state_key] = {
                        k: v for k, v in params.items()
                        if k == ph_param or k.startswith(proto_prefix)
                    }

        self.transformation_parameters = {}
        return True

    def _select_state_combinations(self) -> bool:
        """Select redox/spin state combinations for each site (Step 3)"""
        self.console.print("\n[bold]Step 3: Select Redox/Spin State Combinations[/bold]")

        # Build state information for all sites
        site_state_info = {}  # site_number -> {site_id, transformer_name, available_states}

        for i, site in enumerate(self.redox_sites):
            site_num = i + 1
            site_id = site.site_id if hasattr(site, 'site_id') else f'site_{i}'
            transformer_name = self.site_transformer_assignments[site_id]

            available_states = self._available_states_for_transformer(transformer_name)

            site_state_info[site_num] = {
                'site_id': site_id,
                'transformer_name': transformer_name,
                'available_states': available_states
            }

        # Display combined table of all sites and their states
        table = Table(title="Available Redox/Spin States for All Sites")
        table.add_column("Site #", style="cyan", width=7)
        table.add_column("Center(s)", style="yellow", width=22)
        table.add_column("Transformer", style="magenta", width=20)
        table.add_column("State Options", style="green")

        default_selections = []

        for site_num in sorted(site_state_info.keys()):
            info = site_state_info[site_num]

            # Build center description from RedoxSite centers
            site_idx = site_num - 1
            center_desc = ""
            if site_idx < len(self.redox_sites):
                site = self.redox_sites[site_idx]
                if hasattr(site, 'centers') and site.centers:
                    center_parts = []
                    for c in site.centers:
                        center_parts.append(f"{c.chain}:{c.resname}{c.resid}")
                    center_desc = ", ".join(center_parts)
                elif hasattr(site, 'atoms') and site.atoms:
                    # Fallback: summarize from atoms
                    residues = set()
                    for a in site.atoms:
                        residues.add(f"{a.chain}:{a.resname}{a.resid}")
                    center_desc = ", ".join(sorted(residues)[:3])
                    if len(residues) > 3:
                        center_desc += f" +{len(residues)-3}"

            # Format state options
            state_lines = []
            valid_state_nums = []
            for state_idx, state in enumerate(info['available_states'], 1):
                valid_marker = "✓" if state['is_valid'] else "✗"
                default_marker = "★" if state.get('is_default', False) else " "
                state_lines.append(
                    f"{state_idx}: {state['redox_state']}/{state['spin_state']} {valid_marker}{default_marker}"
                )
                if state['is_valid']:
                    valid_state_nums.append(str(state_idx))

            state_options = "\n".join(state_lines)

            table.add_row(
                str(site_num),
                center_desc,
                info['transformer_name'],
                state_options
            )

            # Build default selection (all valid states for each site)
            default_selections.append(f"{site_num}:{','.join(valid_state_nums)}")

        self.console.print(table)

        # Show instructions and default
        self.console.print("\n[bold]Selection Syntax:[/bold]")
        self.console.print("  [site#]:[state1],[state2],... [site#]:[state1],[state2],...", highlight=False)
        self.console.print("\n[bold]Standard Examples:[/bold]")
        self.console.print("  1:1,2 2:1,2 3:1      - Sites 1&2 use states 1&2, site 3 uses only state 1", highlight=False)
        self.console.print("  1-5:1,2 6:1          - Sites 1-5 use states 1&2, site 6 uses only state 1", highlight=False)
        self.console.print("  all:1,2              - All sites use states 1&2", highlight=False)
        self.console.print("\n[bold cyan]Systematic Selections (Recommended for many sites):[/bold cyan]")
        self.console.print("  single:[ref_state]          - Single-site variants (vary one site at a time from reference)", highlight=False)
        self.console.print("  single:[ref_state][X-Y]     - Single-site variants only for sites X through Y", highlight=False)
        self.console.print("  pairs:[ref_state]           - Pairwise variants (vary all pairs from reference)", highlight=False)
        self.console.print("  pairs:[ref_state][X-Y]      - Pairwise variants only for sites X through Y", highlight=False)
        self.console.print("\n[bold]Systematic Examples:[/bold]")
        self.console.print("  single:1                    - All sites in state 1, except vary each site through all states", highlight=False)
        self.console.print("  single:2[1-23]              - All sites in state 2, except vary sites 1-23 one at a time", highlight=False)
        self.console.print("  pairs:1                     - All sites in state 1, except vary all pairs through all states", highlight=False)
        self.console.print("  pairs:1[5-10]               - All sites in state 1, except vary sites 5-10 in pairs", highlight=False)
        self.console.print("  single:1 pairs:1            - Both single-site and pairwise variants", highlight=False)

        # Calculate what the "all sites all states" option would generate
        num_sites = len(site_state_info)
        if num_sites > 10:
            self.console.print(f"\n[yellow]⚠ Warning: You have {num_sites} sites. Using 'all:1,2' would generate 2^{num_sites} = {2**num_sites:,} combinations![/yellow]", highlight=False)
            self.console.print(f"[yellow]  Consider using 'single:1' ({num_sites + 1} microstates) or 'pairs:1' (~{num_sites*(num_sites-1)*3//2:,} microstates) instead.[/yellow]", highlight=False)

        # Loop until user provides valid input
        selection = None
        while not selection or not selection.strip():
            selection = prompt_with_context(
                processor=self.processor,
                prompt="\n[bold]Enter selection (required):[/bold]",
                default="",
                module="Redox Site Transformation",
                description="Select state combinations"
            )

            if not selection or not selection.strip():
                self.console.print("[red]Error: Selection is required. Please specify states for microstate generation.[/red]")
                self.console.print("[yellow]Try 'single:1' for single-site variants or 'pairs:1' for pairwise variants.[/yellow]", highlight=False)

        # Parse the selection
        self.site_state_selections = {}

        try:
            self.site_state_selections = self._parse_state_selection(
                selection, site_state_info
            )

            # Show summary
            self.console.print("\n[bold green]Selection Summary:[/bold green]")
            for site_num in sorted(site_state_info.keys()):
                info = site_state_info[site_num]
                site_id = info['site_id']
                selected = self.site_state_selections.get(site_id, [])
                state_desc = ", ".join([f"{s['redox_state']}/{s['spin_state']}" for s in selected])
                self.console.print(f"  Site {site_num} ({site_id}): {len(selected)} states - {state_desc}")

            return True

        except Exception as e:
            self.console.print(f"[red]Error parsing selection: {e}[/red]")
            return False

    def _parse_state_selection(self, selection: str, site_state_info: Dict) -> Dict:
        """Parse state selection string into site_state_selections dict"""
        result = {}

        # Check for special systematic selections first
        if 'single:' in selection or 'pairs:' in selection:
            return self._parse_systematic_selection(selection, site_state_info)

        # Split by spaces to get individual site selections
        parts = selection.strip().split()

        for part in parts:
            if ':' not in part:
                continue

            site_spec, state_spec = part.split(':', 1)

            # Parse site specification (can be single number, range, or "all")
            if site_spec.lower() == 'all':
                site_nums = list(site_state_info.keys())
            elif '-' in site_spec:
                start, end = site_spec.split('-', 1)
                site_nums = list(range(int(start), int(end) + 1))
            else:
                site_nums = [int(site_spec)]

            # Parse state specification (comma-separated state numbers)
            state_indices = [int(x.strip()) - 1 for x in state_spec.split(',') if x.strip()]

            # Apply to all specified sites
            for site_num in site_nums:
                if site_num not in site_state_info:
                    continue

                info = site_state_info[site_num]
                site_id = info['site_id']
                available_states = info['available_states']

                selected_states = []
                for idx in state_indices:
                    if 0 <= idx < len(available_states):
                        state = available_states[idx]
                        selected_states.append({
                            'redox_state': state['redox_state'],
                            'spin_state': state['spin_state']
                        })

                result[site_id] = selected_states

        return result

    def _parse_systematic_selection(self, selection: str, site_state_info: Dict) -> Dict:
        """
        Parse systematic selection patterns like 'single:1' or 'pairs:1'

        Supports site ranges:
        - single:2              - All sites vary, reference is state 2
        - single:2[1-23]        - Only sites 1-23 vary, reference is state 2
        - pairs:1[5-10]         - Only sites 5-10 vary in pairs, reference is state 1

        Returns a dict where each site maps to a list of states it should explore,
        which will later be expanded into specific microstate combinations.
        """
        import re

        # Extract single: and pairs: patterns with optional site range
        single_match = re.search(r'single:(\d+)(?:\[(\d+)-(\d+)\])?', selection)
        pairs_match = re.search(r'pairs:(\d+)(?:\[(\d+)-(\d+)\])?', selection)

        if not single_match and not pairs_match:
            raise ValueError("No valid systematic selection found")

        # Store systematic selection metadata that will be used during microstate generation
        # We'll store this in a special format that _generate_microstate_combinations can detect
        result = {}

        # Mark that we're using systematic selection
        result['__systematic__'] = {
            'type': [],
            'site_state_info': site_state_info
        }

        if single_match:
            ref_state_idx = int(single_match.group(1)) - 1
            # Check if site range is specified
            if single_match.group(2) and single_match.group(3):
                site_range = (int(single_match.group(2)), int(single_match.group(3)))
            else:
                site_range = None
            result['__systematic__']['type'].append(('single', ref_state_idx, site_range))

        if pairs_match:
            ref_state_idx = int(pairs_match.group(1)) - 1
            # Check if site range is specified
            if pairs_match.group(2) and pairs_match.group(3):
                site_range = (int(pairs_match.group(2)), int(pairs_match.group(3)))
            else:
                site_range = None
            result['__systematic__']['type'].append(('pairs', ref_state_idx, site_range))

        # Store all available states for each site (needed for systematic generation)
        for site_num, info in site_state_info.items():
            site_id = info['site_id']
            result[site_id] = info['available_states']

        return result

    def _generate_microstate_combinations(self) -> List[Dict]:
        """Generate all microstate combinations from selected states (Step 4)"""
        from itertools import product, combinations
        import math

        # Check if systematic selection was used
        if '__systematic__' in self.site_state_selections:
            return self._generate_systematic_microstates()

        # Get site IDs in order
        site_ids = [site.site_id if hasattr(site, 'site_id') else f'site_{i}'
                   for i, site in enumerate(self.redox_sites)]

        # Get state options for each site
        site_options = [self.site_state_selections[site_id] for site_id in site_ids]

        # Calculate total number of combinations
        total_combinations = math.prod(len(options) for options in site_options)

        self.console.print(f"\n[bold]Calculating microstate combinations...[/bold]")
        self.console.print(f"  Sites: {len(site_ids)}")
        self.console.print(f"  Total combinations: {total_combinations:,}")

        # Sanity check - warn if too many combinations
        if total_combinations > 1_000_000:
            self.console.print(f"\n[bold red]WARNING: {total_combinations:,} combinations is extremely large![/bold red]")
            self.console.print("[yellow]This will consume significant memory and may crash the system.[/yellow]")
            self.console.print("\n[bold]Suggestions:[/bold]")
            self.console.print("  1. Reduce the number of states per site")
            self.console.print("  2. Select fewer sites to vary")
            self.console.print("  3. Generate microstates in smaller batches")

            if not confirm_with_context(
                processor=self.processor,
                prompt="\nContinue anyway?",
                default=False,
                module="Redox Site Transformation",
                description="Continue with large combination count"
            ):
                return []

        elif total_combinations > 10_000:
            self.console.print(f"\n[yellow]Warning: {total_combinations:,} combinations is quite large[/yellow]")
            if not confirm_with_context(
                processor=self.processor,
                prompt="Continue?",
                default=True,
                module="Redox Site Transformation",
                description="Continue with large combination count"
            ):
                return []

        # Generate Cartesian product of all combinations
        self.console.print("[grey50]Generating combinations...[/grey50]")
        combinations_list = list(product(*site_options))

        # Convert to list of microstate dictionaries
        microstates = []
        for combo in combinations_list:
            microstate = {}
            for site_id, state in zip(site_ids, combo):
                microstate[site_id] = state
            microstates.append(microstate)

        return microstates

    def _generate_systematic_microstates(self) -> List[Dict]:
        """Generate systematic single-site and/or pairwise microstates"""
        from itertools import combinations, product

        systematic_info = self.site_state_selections['__systematic__']
        site_state_info = systematic_info['site_state_info']
        selection_types = systematic_info['type']

        # Get site IDs in order
        site_ids = [site.site_id if hasattr(site, 'site_id') else f'site_{i}'
                   for i, site in enumerate(self.redox_sites)]

        all_microstates = []

        for selection_tuple in selection_types:
            sel_type = selection_tuple[0]
            ref_state_idx = selection_tuple[1]
            site_range = selection_tuple[2] if len(selection_tuple) > 2 else None

            if sel_type == 'single':
                microstates = self._generate_single_site_variants(site_ids, site_state_info, ref_state_idx, site_range)
                all_microstates.extend(microstates)
                if site_range:
                    self.console.print(f"[green]Generated {len(microstates)} single-site variants (sites {site_range[0]}-{site_range[1]})[/green]")
                else:
                    self.console.print(f"[green]Generated {len(microstates)} single-site variants[/green]")

            elif sel_type == 'pairs':
                microstates = self._generate_pairwise_variants(site_ids, site_state_info, ref_state_idx, site_range)
                all_microstates.extend(microstates)
                if site_range:
                    self.console.print(f"[green]Generated {len(microstates)} pairwise variants (sites {site_range[0]}-{site_range[1]})[/green]")
                else:
                    self.console.print(f"[green]Generated {len(microstates)} pairwise variants[/green]")

        # Remove duplicates (in case single and pairs overlap)
        unique_microstates = []
        seen = set()
        for ms in all_microstates:
            # Create a hashable representation
            key = tuple(sorted((sid, s['redox_state'], s['spin_state']) for sid, s in ms.items()))
            if key not in seen:
                seen.add(key)
                unique_microstates.append(ms)

        return unique_microstates

    def _generate_single_site_variants(self, site_ids: List[str], site_state_info: Dict, ref_state_idx: int, site_range: Optional[Tuple[int, int]] = None) -> List[Dict]:
        """
        Generate single-site variants: all sites at reference state except vary one site at a time

        For N sites with M states each:
        - 1 reference microstate (all sites at ref_state)
        - N × (M-1) variants where one site differs from reference

        Args:
            site_ids: List of site IDs
            site_state_info: Dict with site state information
            ref_state_idx: Index of reference state (0-indexed)
            site_range: Optional tuple (start, end) to limit which sites vary (1-indexed, inclusive)
        """
        microstates = []

        # Create reference state (all sites at reference)
        ref_microstate = {}
        site_actual_ref = {}  # Track actual ref index per site (may differ if site has fewer states)
        for i, site_id in enumerate(site_ids):
            site_num = i + 1
            info = site_state_info[site_num]
            available_states = info['available_states']

            # Use requested ref state, or fall back to first state if site has fewer states
            actual_ref = ref_state_idx if ref_state_idx < len(available_states) else 0
            site_actual_ref[site_id] = actual_ref

            ref_state = available_states[actual_ref]
            ref_microstate[site_id] = {
                'redox_state': ref_state['redox_state'],
                'spin_state': ref_state['spin_state']
            }

        # Add reference microstate
        microstates.append(ref_microstate)

        # Generate single-site variants
        for i, site_id in enumerate(site_ids):
            site_num = i + 1

            # Check if this site is in the allowed range
            if site_range:
                if site_num < site_range[0] or site_num > site_range[1]:
                    continue  # Skip sites outside the range

            info = site_state_info[site_num]
            available_states = info['available_states']
            actual_ref = site_actual_ref[site_id]

            # Try each state for this site (except its reference)
            for state_idx, state in enumerate(available_states):
                if state_idx == actual_ref:
                    continue  # Skip reference state

                # Create variant with this site in different state
                variant = ref_microstate.copy()
                variant[site_id] = {
                    'redox_state': state['redox_state'],
                    'spin_state': state['spin_state']
                }
                microstates.append(variant)

        return microstates

    def _generate_pairwise_variants(self, site_ids: List[str], site_state_info: Dict, ref_state_idx: int, site_range: Optional[Tuple[int, int]] = None) -> List[Dict]:
        """
        Generate pairwise variants: reference + all sites at reference except vary two at a time

        Args:
            site_ids: List of site IDs
            site_state_info: Dict with site state information
            ref_state_idx: Index of reference state (0-indexed)
            site_range: Optional tuple (start, end) to limit which sites vary (1-indexed, inclusive)

        For N sites with 2 states each:
        - 1 reference microstate (all sites at ref_state)
        - C(N, 2) pairs = N*(N-1)/2 pairs
        - 3 non-ref combinations per pair (ref-var, var-ref, var-var)
        - Total: 1 + C(N, 2) × 3 microstates
        """
        from itertools import combinations

        microstates = []

        # Create reference state
        ref_microstate = {}
        site_actual_ref = {}  # Track actual ref index per site
        for i, site_id in enumerate(site_ids):
            site_num = i + 1
            info = site_state_info[site_num]
            available_states = info['available_states']

            actual_ref = ref_state_idx if ref_state_idx < len(available_states) else 0
            site_actual_ref[site_id] = actual_ref

            ref_state = available_states[actual_ref]
            ref_microstate[site_id] = {
                'redox_state': ref_state['redox_state'],
                'spin_state': ref_state['spin_state']
            }

        # Add reference microstate so `pairs:` alone includes the all-reference state.
        # When combined with `single:`, the dedup pass in _generate_microstate_combinations
        # collapses the overlap.
        microstates.append(ref_microstate)

        # Filter sites based on range if specified
        if site_range:
            filtered_enum_sites = [(i, sid) for i, sid in enumerate(site_ids)
                                   if site_range[0] <= (i + 1) <= site_range[1]]
        else:
            filtered_enum_sites = list(enumerate(site_ids))

        # Generate all pairs of sites (within the filtered range)
        site_pairs = list(combinations(filtered_enum_sites, 2))

        for (idx1, site_id1), (idx2, site_id2) in site_pairs:
            site_num1 = idx1 + 1
            site_num2 = idx2 + 1

            info1 = site_state_info[site_num1]
            info2 = site_state_info[site_num2]

            states1 = info1['available_states']
            states2 = info2['available_states']

            actual_ref1 = site_actual_ref[site_id1]
            actual_ref2 = site_actual_ref[site_id2]

            # Generate all combinations of states for this pair
            for state1 in states1:
                for state2 in states2:
                    # Skip if both are reference (already covered in single-site)
                    state1_is_ref = (states1.index(state1) == actual_ref1)
                    state2_is_ref = (states2.index(state2) == actual_ref2)

                    if state1_is_ref and state2_is_ref:
                        continue

                    # Create variant
                    variant = ref_microstate.copy()
                    variant[site_id1] = {
                        'redox_state': state1['redox_state'],
                        'spin_state': state1['spin_state']
                    }
                    variant[site_id2] = {
                        'redox_state': state2['redox_state'],
                        'spin_state': state2['spin_state']
                    }
                    microstates.append(variant)

        return microstates

    def _display_microstates(self, microstates: List[Dict]):
        """Display available microstates (Step 5)"""
        self.console.print(f"\n[bold]Step 4: Available Microstates[/bold]")
        self.console.print(f"Total microstates: {len(microstates):,}")

        # Explain the encoding
        self.console.print("\n[bold]Microstate Encoding:[/bold]")
        self.console.print("  Each site is represented by a 2-letter code:")
        self.console.print("    First letter: R=Reduced, O=Oxidized")
        self.console.print("    Second letter: L=Low spin, H=High spin, I=Intermediate spin")
        self.console.print(f"  Example: RLOL = Site1:Reduced/Low, Site2:Oxidized/Low", highlight=False)

        # Show all microstates
        self.console.print(f"\n[bold]Available Microstates:[/bold]")
        for i, microstate in enumerate(microstates):
            code, _ = self._format_microstate(microstate)
            self.console.print(f"  {i+1}. {code}")

    def _format_microstate(self, microstate: Dict) -> Tuple[str, str]:
        """Format microstate for display"""
        codes = []
        descriptions = []

        for site_id, state in microstate.items():
            # Create short code: R/O for redox, L/H/I for spin
            redox_code = 'R' if state['redox_state'] == 'reduced' else 'O'
            spin_state = state['spin_state']
            if 'low' in spin_state:
                spin_code = 'L'
            elif 'high' in spin_state:
                spin_code = 'H'
            elif 'intermediate' in spin_state:
                spin_code = 'I'
            else:
                spin_code = 'L'

            codes.append(f"{redox_code}{spin_code}")
            descriptions.append(f"{site_id}:{state['redox_state']}_{state['spin_state']}")

        return ''.join(codes), ', '.join(descriptions)

    def _select_microstates_to_generate(self, microstates: List[Dict]) -> List[Dict]:
        """Let user select which microstates to generate (Step 5)"""
        self.console.print(f"\n[bold]Step 5: Confirm Microstates to Generate[/bold]")
        self.console.print(f"Your selection in Step 3 produced {len(microstates)} microstates (listed above).")
        self.console.print("1. Generate all of them", highlight=False)
        self.console.print("2. Pick a subset by number (e.g. 1,3,5 or 1-4,7)", highlight=False)

        choice = prompt_with_context(
            processor=self.processor,
            prompt="Choose",
            choices=["1", "2"],
            default="1",
            module="Redox Site Transformation",
            description="Confirm microstates to generate",
            options_map={
                "1": "Generate all of them",
                "2": "Pick a subset by number"
            }
        )

        if choice == "1":
            # Warn if too many
            if len(microstates) > 100:
                self.console.print(f"[yellow]Warning: This will generate {len(microstates)} microstate files[/yellow]")
                if not confirm_with_context(
                    processor=self.processor,
                    prompt="Continue?",
                    module="Redox Site Transformation",
                    description="Continue with many microstates"
                ):
                    return []
            return microstates

        elif choice == "2":
            return self._select_specific_microstates(microstates)

        return []

    def _select_specific_microstates(self, microstates: List[Dict]) -> List[Dict]:
        """Interactive selection of specific microstates"""
        self.console.print("\nEnter microstate numbers (comma-separated or ranges like 1-5):", highlight=False)
        self.console.print("Example: 1,3,5 or 1-4,7,8", highlight=False)

        selection = prompt_with_context(
            processor=self.processor,
            prompt="Microstate numbers",
            module="Redox Site Transformation",
            description="Enter microstate numbers"
        )
        indices = self._parse_selection_string(selection, len(microstates))

        selected = [microstates[i] for i in indices]
        self.console.print(f"[green]Selected {len(selected)} microstates[/green]")

        return selected

    def _select_pattern_microstates(self, microstates: List[Dict]) -> List[Dict]:
        """Select microstates by pattern"""
        self.console.print("\n[bold]Available Patterns:[/bold]")
        self.console.print("1. All reduced", highlight=False)
        self.console.print("2. All oxidized", highlight=False)
        self.console.print("3. Single oxidized (all others reduced)", highlight=False)
        self.console.print("4. Half-half (approx. equal reduced/oxidized)", highlight=False)

        choice = prompt_with_context(
            processor=self.processor,
            prompt="Select pattern",
            choices=["1", "2", "3", "4"],
            module="Redox Site Transformation",
            description="Select microstate pattern",
            options_map={
                "1": "All reduced",
                "2": "All oxidized",
                "3": "Single oxidized",
                "4": "Half-half"
            }
        )

        filtered = []
        for microstate in microstates:
            states = list(microstate.values())
            reduced_count = sum(1 for s in states if s['redox_state'] == 'reduced')
            oxidized_count = len(states) - reduced_count

            if choice == "1" and reduced_count == len(states):
                filtered.append(microstate)
            elif choice == "2" and oxidized_count == len(states):
                filtered.append(microstate)
            elif choice == "3" and oxidized_count == 1:
                filtered.append(microstate)
            elif choice == "4" and abs(reduced_count - oxidized_count) <= 1:
                filtered.append(microstate)

        self.console.print(f"[green]Found {len(filtered)} matching microstates[/green]")
        return filtered

    def _parse_selection_string(self, selection: str, max_num: int) -> List[int]:
        """Parse selection string like '1,3,5' or '1-4,7,8' into list of indices"""
        indices = []

        for part in selection.split(','):
            part = part.strip()
            if '-' in part:
                start, end = part.split('-', 1)
                start_idx = int(start) - 1
                end_idx = int(end) - 1
                indices.extend(range(start_idx, end_idx + 1))
            else:
                idx = int(part) - 1
                indices.append(idx)

        valid_indices = [i for i in indices if 0 <= i < max_num]
        return sorted(list(set(valid_indices)))

    def _generate_selected_microstates(self, microstates: List[Dict]) -> bool:
        """Generate the selected microstates (Step 7)"""
        self.console.print(f"\n[bold]Step 7: Generating {len(microstates)} Microstates[/bold]")

        # Inform user about memory-efficient processing
        self.console.print("[grey50]ℹ️  Using memory-efficient mode: intermediate debug files disabled, automatic cleanup enabled[/grey50]")

        # Get structure file (needed for transformation input)
        structure_file = self._get_current_structure_file()
        if not structure_file:
            self.console.print("[red]No structure file available[/red]")
            return False
        # Use a consistent base name for transformed outputs (not the input filename)
        base_name = "transformed_microstate"

        successful = 0
        failed = 0
        skipped = 0  # Track files that already exist and were skipped
        generated_files = []  # Track all generated microstate files

        # Determine if we should show compact output (many sites or many microstates)
        num_sites = len(next(iter(microstates), {}))
        num_microstates = len(microstates)
        show_compact = num_sites > 10 or num_microstates > 20

        # Use progress bar for many microstates
        if show_compact and num_microstates > 10:
            from rich.progress import Progress, SpinnerColumn, BarColumn, TaskProgressColumn, TimeElapsedColumn, TimeRemainingColumn

            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TaskProgressColumn(),
                TimeElapsedColumn(),
                TimeRemainingColumn(),
                console=self.console,
                speed_estimate_period=600,  # Use longer window for slow tasks
            ) as progress:
                task = progress.add_task(
                    f"Generating {num_microstates} microstates...",
                    total=num_microstates
                )

                for i, microstate in enumerate(microstates, 1):
                    code, desc = self._format_microstate(microstate)
                    output_file = f"{base_name}_{i:03d}.pdb"
                    output_path = Path(output_file)

                    # Update progress description
                    progress.update(task, description=f"[cyan]{i}/{num_microstates}:[/cyan] {desc[:50]}...")

                    # Check if file already exists
                    is_final = (i == num_microstates)
                    if output_path.exists() and output_path.stat().st_size > 0:
                        skipped += 1
                        generated_files.append({
                            'file': str(output_path.resolve()),
                            'code': code,
                            'description': desc,
                            'microstate': microstate,
                            'skipped': True
                        })
                        # For the final microstate, update RedoxSites from existing PDB
                        # to ensure workspace has transformed residue names
                        if is_final:
                            self._update_redox_sites_from_existing_pdb(
                                str(output_path), self.redox_sites
                            )
                            workspace = self.processor.workspace
                            workspace.set("detected_redox_sites", self.redox_sites)
                        progress.advance(task)
                        continue

                    # Generate this microstate
                    success = self._generate_single_microstate(microstate, output_file, is_final)

                    if success:
                        successful += 1
                        generated_files.append({
                            'file': str(output_path.resolve()),
                            'code': code,
                            'description': desc,
                            'microstate': microstate
                        })
                    else:
                        failed += 1

                    progress.advance(task)

            # Show skip summary after progress bar completes
            if skipped > 0:
                self.console.print(f"[grey50]  Skipped {skipped} existing file(s)[/grey50]")
        else:
            # Non-compact: detailed output for each microstate
            for i, microstate in enumerate(microstates, 1):
                code, desc = self._format_microstate(microstate)

                if show_compact:
                    # Compact: One-line format with just the code/description
                    self.console.print(f"\n[cyan]{i}/{num_microstates}:[/cyan] {desc}")
                else:
                    # Verbose: Show full state list
                    self.console.print(f"\n[cyan]Generating {i}/{num_microstates}: {code}[/cyan]")
                    self.console.print("  States:")
                    for site_id, state in microstate.items():
                        state_desc = f"{state['redox_state']}_{state['spin_state']}"
                        self.console.print(f"    {site_id:>10}: {state_desc}")

                # Generate filename
                output_file = f"{base_name}_{i:03d}.pdb"

                # Check if file already exists and is valid
                output_path = Path(output_file)
                is_final = (i == num_microstates)
                if output_path.exists() and output_path.stat().st_size > 0:
                    skipped += 1
                    generated_files.append({
                        'file': str(output_path.resolve()),
                        'code': code,
                        'description': desc,
                        'microstate': microstate,
                        'skipped': True
                    })
                    # For the final microstate, update RedoxSites from existing PDB
                    # to ensure workspace has transformed residue names
                    if is_final:
                        self._update_redox_sites_from_existing_pdb(
                            str(output_path), self.redox_sites
                        )
                        workspace = self.processor.workspace
                        workspace.set("detected_redox_sites", self.redox_sites)
                    self.console.print(f"[grey50]✓ Already exists (skipped): {output_file}[/grey50]")
                    continue

                # Generate this microstate
                success = self._generate_single_microstate(microstate, output_file, is_final)

                if success:
                    successful += 1
                    generated_files.append({
                        'file': str(output_path.resolve()),
                        'code': code,
                        'description': desc,
                        'microstate': microstate
                    })
                    self.console.print(f"[green]✓ Generated: {output_file}[/green]")
                else:
                    failed += 1
                    self.console.print(f"[red]✗ Failed: {output_file}[/red]")

        self.console.print(f"\n[bold]Generation Summary:[/bold]")
        self.console.print(f"  Newly generated: {successful}")
        if skipped > 0:
            self.console.print(f"  Already existed (skipped): {skipped}")
        if failed > 0:
            self.console.print(f"  Failed: {failed}")
        self.console.print(f"  Total available: {successful + skipped}")

        # Save generated microstate files to workspace
        if generated_files:
            workspace = self.processor.workspace
            workspace.set("generated_microstate_pdbs", generated_files)
            self.console.print(f"\n[grey50]Saved {len(generated_files)} microstate file references to workspace[/grey50]")

        # Generate metadata file for tLEaP input generation (includes both generated and skipped)
        total_available = successful + skipped
        if total_available > 0:
            # Use generated_files to get the microstates that have PDB files available
            available_microstates = [gf['microstate'] for gf in generated_files]
            self._generate_microstate_metadata(available_microstates, base_name)

        return total_available > 0

    def _generate_microstate_metadata(self, microstates: List[Dict], base_name: str):
        """
        Generate metadata JSON file for tLEaP input generation

        Args:
            microstates: List of successfully generated microstates
            base_name: Base filename for output files
        """
        import json
        from datetime import datetime
        from proprep.forcefield_params import load_forcefield_metadata

        try:
            metadata = {
                "generation_date": datetime.now().isoformat(),
                "original_structure": base_name,
                "total_microstates": len(microstates),
                "successful_generations": len(microstates),
                "microstates": []
            }

            for ms_idx, microstate in enumerate(microstates, 1):
                code, _ = self._format_microstate(microstate)
                filename = f"{base_name}_{ms_idx:03d}.pdb"
                redox_sites_json = f"{base_name}_{ms_idx:03d}_redox_sites_updated.json"

                microstate_info = {
                    "filename": filename,
                    "redox_sites_json": redox_sites_json,
                    "code": code,
                    "sites": []
                }

                # Process each site in this microstate
                for site_id, state in microstate.items():
                    transformer_type = self.site_transformer_assignments[site_id]
                    redox_state = state['redox_state']
                    spin_state = state['spin_state']

                    # Get forcefield information
                    atom_types = []
                    residue_name = "UNK"
                    spin_multiplicity = None

                    try:
                        from proprep.redoxsite_prep.transformation.redox_transformer_framework import redox_transformer_registry
                        transformer_class = redox_transformer_registry.get_transformer(transformer_type)

                        if transformer_class and getattr(transformer_class, 'FORCEFIELD_PATH', None):
                            ff_metadata = load_forcefield_metadata(transformer_class.FORCEFIELD_PATH)
                            spin_data = ff_metadata['redox_states'][redox_state]['spin_states'][spin_state]

                            atom_types = spin_data.get('atom_types', [])
                            residue_name = spin_data.get('residue_name', 'UNK')
                            spin_multiplicity = spin_data.get('spin_multiplicity')

                    except Exception as e:
                        logger.warning(f"Could not load forcefield info for {transformer_type}: {e}")

                    # Propionate/pH treatment chosen per redox state (Step 6) and
                    # the forcefield set it implies, so the metadata is
                    # self-describing per microstate (the heme lib differs between
                    # the fixed-pH and constant-pH sets).
                    state_params = self._resolve_state_params(site_id, redox_state, spin_state)
                    ff_set = None
                    try:
                        from proprep.redoxsite_prep.transformation.redox_transformer_framework import redox_transformer_registry
                        tclass = redox_transformer_registry.get_transformer(transformer_type)
                        if tclass and hasattr(tclass, 'select_forcefield_set_name'):
                            ff_set = tclass.select_forcefield_set_name(state_params)
                    except Exception:
                        pass

                    site_info = {
                        "site_id": site_id,
                        "transformer_type": transformer_type,
                        "redox_state": redox_state,
                        "spin_state": spin_state,
                        "residue_name": residue_name,
                        "ph_treatment": state_params.get('ph_treatment'),
                        "forcefield_set": ff_set,
                        "forcefield_info": {
                            "atom_types": atom_types,
                            "spin_multiplicity": spin_multiplicity
                        }
                    }

                    microstate_info["sites"].append(site_info)

                metadata["microstates"].append(microstate_info)

            # Write metadata file
            metadata_file = Path("microstate_metadata.json")
            with open(metadata_file, 'w') as f:
                json.dump(metadata, indent=2, fp=f)

            # Expose the resolved path in workspace so downstream
            # consumers (including the etanalyze plugin's auto-load)
            # don't need to guess where the file landed. ``resolve()``
            # captures the absolute path at write time — if the user
            # later moves the file, the workspace value goes stale,
            # which is the user's responsibility to refresh.
            self.processor.workspace.set(
                "microstate_metadata_path", str(metadata_file.resolve()),
            )

            self.console.print(f"\n[green]✓ Microstate metadata saved: {metadata_file}[/green]")
            self.console.print(f"[grey50]  Use this file for tLEaP input generation (option 3 in Topology Generator menu)[/grey50]")

        except Exception as e:
            logger.error(f"Failed to generate microstate metadata: {e}")
            self.console.print(f"[yellow]Warning: Could not generate metadata file: {e}[/yellow]")

    def _generate_single_microstate(self, microstate: Dict, output_filename: str, is_final: bool) -> bool:
        """
        Generate a single microstate PDB by running transformation with specific redox/spin states

        Args:
            microstate: Dict mapping site_id -> {'redox_state': str, 'spin_state': str}
            output_filename: Output PDB filename
            is_final: If True, keep workspace in transformed state; otherwise restore

        Returns:
            True if successful
        """
        try:
            import copy

            # Save original workspace state
            workspace = self.processor.workspace
            original_structure = workspace.get("untransformed_structure") or self._get_current_structure_file()

            # Create deep copy of redox sites for this transformation
            sites_for_transformation = copy.deepcopy(self.original_redox_sites)

            # Reset workspace to untransformed state
            workspace.set("detected_redox_sites", sites_for_transformation)

            # Clear transformation state
            if workspace.has("transformed_pdb_file"):
                workspace.set("transformed_pdb_file", None)

            # Create new transformation manager instance for this microstate
            microstate_manager = RedoxTransformationManager(self.processor, self.console)

            # Load sites and structure
            microstate_manager.redox_sites = sites_for_transformation
            microstate_manager.original_redox_sites = copy.deepcopy(sites_for_transformation)

            # Read structure
            with open(original_structure, 'r') as f:
                microstate_manager.structure_lines = f.readlines()

            # Set up temp directory for intermediate files
            import tempfile
            microstate_manager.temp_dir = Path(tempfile.mkdtemp(prefix=f"microstate_{Path(output_filename).stem}_"))
            microstate_manager.temp_dir.mkdir(exist_ok=True)
            (microstate_manager.temp_dir / "transformation_reports").mkdir(exist_ok=True)

            # Update transformation executor with temp directory
            microstate_manager.transformation_executor.temp_dir = microstate_manager.temp_dir

            # Set components to non-verbose mode for batch processing
            microstate_manager.id_mapper.verbose = False
            microstate_manager.transformation_executor.verbose = False

            # Capture original states (suppress output)
            import io
            import sys
            old_stdout = sys.stdout
            sys.stdout = io.StringIO()
            try:
                microstate_manager._capture_original_states()
            finally:
                sys.stdout = old_stdout

            # Set up transformer assignments from batch configuration
            microstate_manager.selected_transformers = {}
            microstate_manager.transformation_parameters = {}

            for site_id, states in microstate.items():
                # Get transformer name for this site
                transformer_name = self.site_transformer_assignments[site_id]
                microstate_manager.selected_transformers[site_id] = transformer_name

                # Set parameters including redox/spin state for this microstate,
                # then merge the propionate/pH treatment configured for THIS
                # site's chosen state (Step 6). Bound per redox state, so
                # reduced vs oxidized may carry different PRN/PRD/PRP. Absent =>
                # the transformer's own default (matches the single path).
                params = {
                    'redox_state': states['redox_state'],
                    'spin_state': states['spin_state']
                }
                state_key = f"{states['redox_state']}_{states['spin_state']}"
                params.update(
                    self.site_state_protonation.get(site_id, {}).get(state_key, {})
                )
                microstate_manager.transformation_parameters[site_id] = params

            # Enable batch mode to suppress verbose output
            microstate_manager._batch_mode = True
            # Configure verbosity based on batch mode (sets MINIMAL verbosity)
            microstate_manager._configure_verbosity()

            # Run ID mapping analysis (important to avoid residue ID conflicts!)
            # Suppress output in batch mode
            if not microstate_manager._handle_id_mapping():
                microstate_manager.id_mapping = {}

            # Perform transformations (verbose output suppressed by batch mode)
            if not microstate_manager._perform_transformations():
                return False

            # Process final PDB: reorder, renumber, update sites (silent in batch mode)
            microstate_manager._process_final_pdb()

            # Save output with REMARK header documenting the microstate
            output_path = Path(output_filename)
            code, _ = self._format_microstate(microstate)
            remark_lines = [
                f"REMARK   MICROSTATE {output_path.stem}\n",
                f"REMARK   CODE {code}\n",
            ]
            for site_id, states in microstate.items():
                transformer = self.site_transformer_assignments.get(site_id, 'unknown')
                state_key = f"{states['redox_state']}_{states['spin_state']}"
                proton = self.site_state_protonation.get(site_id, {}).get(state_key, {})
                proton_note = ""
                if proton:
                    proton_note = " " + " ".join(f"{k}={v}" for k, v in sorted(proton.items()))
                remark_lines.append(
                    f"REMARK   SITE {site_id} {states['redox_state']}/{states['spin_state']} "
                    f"({transformer}){proton_note}\n"
                )
            with open(output_path, 'w') as f:
                f.writelines(remark_lines)
                f.writelines(microstate_manager.structure_lines)

            # Export updated RedoxSites JSON for this microstate
            microstate_manager._export_updated_redox_sites(output_path)

            # Clean up intermediate structure files to free disk space
            microstate_manager.transformation_executor.cleanup_intermediate_structures()

            # Clean up temporary directory for this microstate
            try:
                import shutil
                if microstate_manager.temp_dir and microstate_manager.temp_dir.exists():
                    shutil.rmtree(microstate_manager.temp_dir)
            except Exception as e:
                logger.warning(f"Failed to cleanup temporary directory: {e}")

            # Restore workspace if not final; if final, keep transformed state
            if not is_final:
                workspace.set("detected_redox_sites", copy.deepcopy(self.original_redox_sites))
            else:
                # Final microstate: save the transformed redox sites to workspace
                # so downstream steps (tLEaP, etc.) see updated residue names/bonds
                workspace.set("detected_redox_sites", microstate_manager.redox_sites)
                self.redox_sites = microstate_manager.redox_sites

            return True

        except Exception as e:
            logger.error(f"Error generating microstate {output_filename}: {e}")
            self.console.print(f"[red]  Error: {e}[/red]")
            import traceback
            logger.error(traceback.format_exc())
            return False

    def _update_redox_site_definitions(self):
        """Update RedoxSite objects after transformations"""
        self.console.print("\n🔗 Updating Redox Site Definitions")
        self.console.print("[grey50]Matching atom coordinates to update site definitions (may take time for large structures)...[/grey50]")

        updated_count = 0

        # Use progress bar for large site counts
        if len(self.redox_sites) >= 20:
            from rich.progress import BarColumn, TaskProgressColumn, TimeRemainingColumn
            with Progress(
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TaskProgressColumn(),
                TimeRemainingColumn(),
                console=self.console,
                transient=False
            ) as progress:
                task = progress.add_task(f"Updating {len(self.redox_sites)} sites...", total=len(self.redox_sites))

                for site in self.redox_sites:
                    self._update_single_site_definition(site)
                    updated_count += 1
                    progress.update(task, advance=1)
        else:
            for site in self.redox_sites:
                self._update_single_site_definition(site)
                updated_count += 1

        self.console.print(f"[green]✅ Updated {updated_count} redox site definitions[/green]")

    def _update_single_site_definition(self, site):
        """Update a single site's definition with transformed structure info"""
        # Update center information
        if site.centers:
            center = site.centers[0]

            # Apply ID mapping if exists
            mapping_key = (center.chain, center.resid)
            if mapping_key in self.id_mapping:
                center.resid = self.id_mapping[mapping_key]

        # Update all atoms with current structure info using coordinate matching
        if hasattr(site, 'atoms') and site.atoms:
            for atom in site.atoms:
                if hasattr(atom, 'coords'):
                    current_info = self._get_current_atom_info_from_coords(atom.coords)
                    if current_info:
                        # Update the atom's PDB identification to match transformed structure
                        atom.chain = current_info['chain']
                        atom.resname = current_info['resname']
                        atom.resid = current_info['resid']
                        atom.atom_name = current_info['atom_name']

        # Update bonds with current structure info
        self._update_bonds_with_current_structure(site)

    def _update_bonds_with_current_structure(self, site):
        """Update all bonds in a site with current structure residue info"""
        if not hasattr(site, 'bonds') or not site.bonds:
            return
            
        for bond in site.bonds:
            # Update atom1_residue_info with current structure state
            if hasattr(bond, 'atom1_coords') and hasattr(bond, 'atom1_residue_info'):
                current_info = self._get_current_atom_info_from_coords(bond.atom1_coords)
                if current_info and bond.atom1_residue_info:
                    bond.atom1_residue_info.update(current_info)
            
            # Update atom2_residue_info with current structure state  
            if hasattr(bond, 'atom2_coords') and hasattr(bond, 'atom2_residue_info'):
                current_info = self._get_current_atom_info_from_coords(bond.atom2_coords)
                if current_info and bond.atom2_residue_info:
                    bond.atom2_residue_info.update(current_info)
    
    def _create_residue_mapping_report(self):
        """Create enhanced residue mapping report in deprecated format"""
        self.console.print("\n[bold]Creating Enhanced Residue Mapping Report[/bold]")
        
        # Generate the enhanced mapping report
        report_content = self._generate_enhanced_mapping_report()
        
        # Save to working directory
        import os
        report_file = Path(os.getcwd()) / "enhanced_residue_mapping.txt"
        with open(report_file, 'w') as f:
            f.write(report_content)
        
        self.console.print(f"[green]✅ Enhanced mapping report saved: {report_file.name}[/green]")
    
    def _display_completion_message(self, output_file: Path):
        """Display final completion message"""
        separator = "═" * 60
        self.console.print(f"\n{separator}")
        self.console.print("🎉 TRANSFORMATION COMPLETED SUCCESSFULLY")
        self.console.print(separator)
        self.console.print(f"✅ Successfully transformed {len(self.selected_transformers)} redox sites")
        self.console.print("🔄 Updated structure in workspace")
        self.console.print("\n📁 Output Files saved to working directory:")
        self.console.print(f"  📄 {output_file.name}")
        self.console.print(f"  📊 enhanced_residue_mapping.txt")
        self.console.print(separator)
    
    def _analyze_atom_transformations(self, site_id: str, transformations: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Analyze transformations to extract atom-level movement information"""
        atom_transformations = []
        
        # Get original atom states for this site
        original_atoms = self.original_atom_states.get(site_id, {})
        if not original_atoms:
            return atom_transformations
        
        # Group transformations by original residue
        residue_transformations = {}
        
        for transform in transformations:
            selector = transform['selector']
            action = transform['action']
            
            chain_id = selector.get('chain_id')
            residue_id = selector.get('residue_id')
            
            if not chain_id or not residue_id:
                continue
            
            # Find original residue name from our captured states
            original_resname = None
            for coords, atom_info in original_atoms.items():
                if (atom_info['chain'] == chain_id and 
                    atom_info['resid'] == residue_id):
                    original_resname = atom_info['resname']
                    break
            
            if not original_resname:
                original_resname = selector.get('residue_name', 'UNK')
            
            res_key = (chain_id, residue_id, original_resname)
            
            if res_key not in residue_transformations:
                residue_transformations[res_key] = []
            
            residue_transformations[res_key].append({
                'transform': transform,
                'final_resname': action.get('change_residue_name', original_resname),
                'atoms_modified': transform.get('lines_modified', 0)
            })
        
        # Analyze each residue's transformation pattern
        for (chain_id, residue_id, original_resname), transforms in residue_transformations.items():
            
            # Count original atoms in this residue
            original_atom_count = sum(1 for atom_info in original_atoms.values() 
                                    if (atom_info['chain'] == chain_id and 
                                        atom_info['resid'] == residue_id and 
                                        atom_info['resname'] == original_resname))
            
            # Determine transformation type
            final_residue_names = set(t['final_resname'] for t in transforms)
            
            if len(final_residue_names) == 1 and list(final_residue_names)[0] == original_resname:
                # No change
                atom_transformations.append({
                    'type': 'unchanged',
                    'residue': {
                        'chain': chain_id,
                        'resname': original_resname,
                        'resid': residue_id
                    },
                    'atom_count': original_atom_count
                })
            
            elif len(final_residue_names) == 1:
                # Simple rename
                final_resname = list(final_residue_names)[0]
                total_atoms_modified = sum(t['atoms_modified'] for t in transforms)
                
                atom_transformations.append({
                    'type': 'rename',
                    'original': {
                        'chain': chain_id,
                        'resname': original_resname,
                        'resid': residue_id
                    },
                    'final': {
                        'chain': chain_id,
                        'resname': final_resname,
                        'resid': residue_id
                    },
                    'atom_count': total_atoms_modified
                })
            
            else:
                # Split into multiple residues
                splits = []
                
                for final_resname in final_residue_names:
                    # Find transforms for this final residue
                    final_transforms = [t for t in transforms if t['final_resname'] == final_resname]
                    total_atoms = sum(t['atoms_modified'] for t in final_transforms)
                    
                    # Try to determine the final residue ID (may be different due to ID mapping)
                    final_resid = residue_id  # Default to original
                    if final_transforms:
                        # Look for any ID changes in the transformations
                        for t in final_transforms:
                            action = t['transform']['action']
                            if 'change_residue_id' in action:
                                final_resid = action['change_residue_id']
                                break
                    
                    # Generic atom detail string - no hardcoded chemistry assumptions
                    atom_detail = f"{total_atoms} atoms"
                    
                    splits.append({
                        'chain': chain_id,
                        'resname': final_resname,
                        'resid': final_resid,
                        'atom_count': total_atoms,
                        'atom_detail': atom_detail
                    })
                
                atom_transformations.append({
                    'type': 'split',
                    'original': {
                        'chain': chain_id,
                        'resname': original_resname,
                        'resid': residue_id
                    },
                    'splits': splits,
                    'total_atoms': original_atom_count
                })
        
        return atom_transformations
    
    def _generate_enhanced_mapping_report(self) -> str:
        """Generate enhanced mapping report in the deprecated format"""
        report_lines = []
        
        # Header
        report_lines.append("# Enhanced Residue Mapping: Original Structure -> Transformed Structure")
        report_lines.append("# Generated by RedoxSite Transformer")
        
        # Calculate statistics
        total_original_atoms = sum(len(atoms) for atoms in self.original_atom_states.values())
        total_original_residues = 0
        total_final_residues = 0
        atoms_tracked = 0
        
        # Count original residues
        original_residues = set()
        for site_id, atoms in self.original_atom_states.items():
            for atom_info in atoms.values():
                res_key = (atom_info['chain'], atom_info['resid'], atom_info['resname'])
                original_residues.add(res_key)
        total_original_residues = len(original_residues)
        
        # Use the stored executed transformations with unique atom counts
        all_transformations = []
        for site_id in self.selected_transformers.keys():
            executed_data = self.executed_transformations.get(site_id, {})
            if isinstance(executed_data, dict):
                transformations = executed_data.get('transformations', [])
                unique_atoms = executed_data.get('unique_atoms_affected', 0)
                atoms_tracked += unique_atoms
            else:
                # Fallback for old format
                transformations = executed_data
                for transform in transformations:
                    atoms_tracked += transform.get('lines_modified', 0)
            
            for transform in transformations:
                all_transformations.append(transform)
        
        # Count final residues (approximate)
        final_residues = set()
        for transform in all_transformations:
            action = transform.get('action', {})
            selector = transform['selector']
            if 'change_residue_name' in action:
                chain = selector.get('chain_id')
                resid = action.get('change_residue_id', selector.get('residue_id'))
                resname = action['change_residue_name']
                final_residues.add((chain, resid, resname))
        total_final_residues = len(final_residues)
        
        report_lines.append(f"# Original: {total_original_atoms} atoms in {total_original_residues} residues")
        report_lines.append(f"# Transformed: {total_original_atoms} atoms in {total_final_residues} residues")
        report_lines.append(f"# Atoms tracked: {atoms_tracked}/{total_original_atoms}")
        report_lines.append("#")
        report_lines.append("# Format:")
        report_lines.append("# Simple mapping: CHAIN:RESNAME RESID -> CHAIN:RESNAME RESID")
        report_lines.append("# Split residue: CHAIN:RESNAME RESID -> SPLIT:")
        report_lines.append("#   - CHAIN:RESNAME RESID: atoms...")
        report_lines.append("")
        
        # Group by chain using executed transformations
        chains_data = {}
        for site_id in self.selected_transformers.keys():
            executed_data = self.executed_transformations.get(site_id, {})
            if isinstance(executed_data, dict):
                executed_transforms = executed_data.get('transformations', [])
                per_atom_tracking = executed_data.get('per_atom_tracking', {})
            else:
                # Fallback for old format
                executed_transforms = executed_data
                per_atom_tracking = {}
            
            atom_transformations = self._analyze_atom_transformations_with_per_atom_data(site_id, executed_transforms, per_atom_tracking)
            
            for transformation in atom_transformations:
                if transformation['type'] in ['rename', 'unchanged']:
                    res_info = transformation.get('residue') or transformation.get('original')
                    chain = res_info['chain']
                elif transformation['type'] == 'split':
                    chain = transformation['original']['chain']
                else:
                    continue
                
                if chain not in chains_data:
                    chains_data[chain] = []
                chains_data[chain].append(transformation)
        
        # Generate report by chain
        for chain in sorted(chains_data.keys()):
            report_lines.append(f"## Chain {chain}")
            report_lines.append("")
            
            for transformation in chains_data[chain]:
                if transformation['type'] == 'unchanged':
                    res = transformation.get('residue') or transformation.get('original')
                    if not res:
                        continue
                    report_lines.append(f"{res['chain']}:{res['resname']} {res['resid']} -> UNCHANGED")
                
                elif transformation['type'] == 'rename':
                    orig = transformation['original']
                    final = transformation['final']
                    report_lines.append(f"{orig['chain']}:{orig['resname']} {orig['resid']} -> {final['chain']}:{final['resname']} {final['resid']} (renamed from {orig['resname']})")
                
                elif transformation['type'] == 'split':
                    orig = transformation['original']
                    splits = transformation['splits']
                    report_lines.append(f"{orig['chain']}:{orig['resname']} {orig['resid']} -> SPLIT:")
                    
                    for split in splits:
                        report_lines.append(f"  - {split['chain']}:{split['resname']} {split['resid']}: {split['atom_detail']}")
                
                report_lines.append("")
        
        return "\n".join(report_lines)
    
    def _get_transformations_for_site(self, site_id: str) -> List[Dict[str, Any]]:
        """Get transformation sequence for a specific site"""
        transformer_name = self.selected_transformers.get(site_id)
        if not transformer_name:
            return []
        
        transformer_class = self.transformer_selector.registry.get_transformer(transformer_name)
        if not transformer_class:
            return []
        
        # Find the site
        site = None
        for s in self.redox_sites:
            if s.site_id == site_id:
                site = s
                break
        
        if not site:
            return []
        
        # Get components and parameters
        site_params = self.transformation_parameters.get(site_id, {})
        
        if '_components' in site_params:
            components = site_params['_components']
        else:
            try:
                components, _ = transformer_class.match_components(site)
            except:
                return []
        
        try:
            clean_params = {k: v for k, v in site_params.items() if not k.startswith('_')}
            return transformer_class.get_transformation_sequence(components, clean_params)
        except:
            return []
    
    def _calculate_unique_atoms_affected(self, site_id: str, transformations: List[Dict[str, Any]]) -> int:
        """Calculate unique atoms affected using coordinates as permanent identifiers"""
        original_atoms = self.original_atom_states.get(site_id, {})
        if not original_atoms:
            return 0
        
        # All atoms in the original redox site are considered "affected" since they were part of transformations
        # Even if some steps show 0 lines_modified, the atoms were still processed through the transformation system
        return len(original_atoms)
    
    def _extract_per_atom_tracking(self, site_id: str, transformations: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Extract per-atom tracking by comparing original vs final RedoxSite states"""
        original_atoms = self.original_atom_states.get(site_id, {})
        if not original_atoms:
            return {}
        
        # Get the current RedoxSite after all transformations
        current_site = None
        for site in self.redox_sites:
            if site.site_id == site_id:
                current_site = site
                break
        
        if not current_site:
            return {}
        
        # Compare original vs final state using coordinates as the permanent identifier
        atom_tracking = {}
        missing_atoms = []
        
        for coords, original_info in original_atoms.items():
            # Find this coordinate in the current RedoxSite
            current_info = current_site.coord_to_pdb.get(coords, {})
            
            if current_info:
                atom_key = f"{coords[0]:.3f},{coords[1]:.3f},{coords[2]:.3f}"
                atom_tracking[atom_key] = {
                    'coords': coords,
                    'original': {
                        'chain': original_info['chain'],
                        'resname': original_info['resname'],
                        'resid': original_info['resid'],
                        'atom_name': original_info['atom_name']
                    },
                    'final': {
                        'chain': current_info.get('chain', original_info['chain']),
                        'resname': current_info.get('resname', original_info['resname']),
                        'resid': current_info.get('resid', original_info['resid']),
                        'atom_name': current_info.get('atom_name', original_info['atom_name'])
                    }
                }
            else:
                # MISSING ATOM - this should never happen!
                missing_atoms.append({
                    'coords': coords,
                    'original': original_info
                })
        
        # DEBUG: Report missing atoms
        if missing_atoms:
            self.console.print(f"[red]🚨 MISSING ATOMS DETECTED! {len(missing_atoms)} atoms missing from final RedoxSite:[/red]")
            for missing in missing_atoms:
                orig = missing['original']
                coords = missing['coords']
                self.console.print(f"[red]  Missing: {orig['chain']}:{orig['resname']}{orig['resid']}:{orig['atom_name']} at {coords}[/red]")
        
        # Debug output suppressed for cleaner interface
        
        return atom_tracking
    
    def _analyze_atom_transformations_with_per_atom_data(self, site_id: str, transformations: List[Dict[str, Any]], per_atom_tracking: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Enhanced atom transformation analysis using coordinate-based tracking"""
        original_atoms = self.original_atom_states.get(site_id, {})
        if not original_atoms or not per_atom_tracking:
            return self._analyze_atom_transformations(site_id, transformations)

        # Group coordinates by their original residue and final destinations
        residue_transformations = {}

        for atom_key, atom_data in per_atom_tracking.items():
            original = atom_data['original']
            current = atom_data['final']  # Final state after all transformations

            # Skip virtual centroid entries (atom_name None on both sides).
            # Centroids are whole-residue position markers, not real atoms; their
            # final-state resname can drift from the real atoms' final resname,
            # which would falsely classify a plain residue rename as a SPLIT
            # and inflate the atom count by 1.
            if original.get('atom_name') is None and current.get('atom_name') is None:
                continue
            
            # Group by original residue
            orig_key = (original['chain'], original['resid'], original['resname'])
            if orig_key not in residue_transformations:
                residue_transformations[orig_key] = {
                    'atoms': [],
                    'final_destinations': {}
                }
            
            # Track final destination for this atom
            final_key = (current['chain'], current['resid'], current['resname'])
            if final_key not in residue_transformations[orig_key]['final_destinations']:
                residue_transformations[orig_key]['final_destinations'][final_key] = []
            
            # Add atom info with any renames
            atom_info = {
                'original_name': original['atom_name'],
                'final_name': current['atom_name'],
                'renamed': original['atom_name'] != current['atom_name']
            }
            residue_transformations[orig_key]['final_destinations'][final_key].append(atom_info)
            residue_transformations[orig_key]['atoms'].append(atom_info)
        
        # Convert to transformation analysis format
        atom_transformations = []
        
        for (orig_chain, orig_resid, orig_resname), res_data in residue_transformations.items():
            total_atoms = len(res_data['atoms'])
            final_destinations = res_data['final_destinations']
            
            if len(final_destinations) == 1:
                # Single destination - either unchanged or renamed
                final_key, atoms = list(final_destinations.items())[0]
                final_chain, final_resid, final_resname = final_key
                
                if orig_resname == final_resname and orig_resid == final_resid:
                    # Unchanged or just atom renames
                    transformation_type = 'unchanged'
                    if any(atom['renamed'] for atom in atoms):
                        transformation_type = 'rename'
                        
                    atom_transformations.append({
                        'type': transformation_type,
                        'original': {
                            'chain': orig_chain,
                            'resname': orig_resname,
                            'resid': orig_resid
                        },
                        'final': {
                            'chain': final_chain,
                            'resname': final_resname,
                            'resid': final_resid
                        },
                        'atom_count': total_atoms,
                        'atom_renames': [f"{a['original_name']}->{a['final_name']}" 
                                       for a in atoms if a['renamed']]
                    })
                else:
                    # Simple rename/move
                    atom_transformations.append({
                        'type': 'rename',
                        'original': {
                            'chain': orig_chain,
                            'resname': orig_resname,
                            'resid': orig_resid
                        },
                        'final': {
                            'chain': final_chain,
                            'resname': final_resname,
                            'resid': final_resid
                        },
                        'atom_count': total_atoms
                    })
            else:
                # Multiple destinations - split
                splits = []
                for final_key, atoms in final_destinations.items():
                    final_chain, final_resid, final_resname = final_key

                    # Filter out centroid entries (atom_name=None) - these are
                    # virtual positions, not real atom transformations
                    real_atoms = [a for a in atoms
                                  if a['original_name'] is not None or a['final_name'] is not None]
                    if not real_atoms:
                        continue

                    # Build atom detail string
                    renamed_atoms = [f"{a['original_name']}->{a['final_name']}"
                                   for a in real_atoms if a['renamed']]
                    unchanged_atoms = [a['original_name'] for a in real_atoms if not a['renamed']]

                    if renamed_atoms:
                        atom_detail = ", ".join(renamed_atoms)
                    elif unchanged_atoms:
                        atom_detail = ", ".join(unchanged_atoms)
                    else:
                        atom_detail = f"{len(real_atoms)} atoms"

                    splits.append({
                        'chain': final_chain,
                        'resname': final_resname,
                        'resid': final_resid,
                        'atom_count': len(real_atoms),
                        'atom_detail': atom_detail
                    })
                
                atom_transformations.append({
                    'type': 'split',
                    'original': {
                        'chain': orig_chain,
                        'resname': orig_resname,
                        'resid': orig_resid
                    },
                    'splits': splits,
                    'total_atoms': total_atoms
                })
        
        return atom_transformations
    
    def _process_final_pdb(self):
        """Process final PDB structure: reorder for TLEaP, renumber residues, update sites"""
        # 1. Reorder for TLEaP compatibility (ATOM before HETATM)
        self.structure_lines = self._reorder_for_tleap(self.structure_lines)
        
        # 2. Renumber residues sequentially starting at 1
        residue_mapping = self._renumber_residues_sequential()
        
        # 3. Update RedoxSite objects with new residue IDs via coordinate matching
        self._update_redox_sites_with_renumbering(residue_mapping)
    
    # Caps are protein residues recorded as HETATM. That is the only place in
    # a structure where the record type disagrees with what the residue IS, and
    # it is what made the old grouping wrong.
    _CAP_RESNAMES = ("ACE", "NME")

    def _reorder_for_tleap(self, pdb_lines):
        """Reorder PDB lines for tLEaP: real HETATM after the protein chains.

        Caps stay where they are. The previous version binned every ACE into
        one list and every NME into another, then wrote
        ``all ACE / all protein / all NME`` per chain -- which assumes a cap can
        only be terminal.

        An INTERNAL cap breaks that assumption. Capping an unfilled gap puts an
        NME after the last residue before it and an ACE before the first
        residue after it, mid-chain. Binning them moved the pair to opposite
        ends of the chain: the ACE landed before residue 1, where tLEaP bonded
        it to the N-terminus 70 A away, and the gap it had been guarding was
        left open for a 20 A peptide bond.

        A TER is emitted after each NME, which is what separates the two
        fragments a capped gap creates. Without it tLEaP bonds straight across.
        """
        header_lines = []
        chain_streams = {}          # chain -> [ATOM and cap lines, in file order]
        hetatm_lines_by_chain = {}  # chain -> [real HETATM lines]
        end_lines = []

        for line in pdb_lines:
            if line.startswith("ATOM"):
                chain_streams.setdefault(line[21], []).append(line)
            elif line.startswith("HETATM"):
                chain_id = line[21]
                if line[17:20].strip() in self._CAP_RESNAMES:
                    # Part of the peptide chain; keep it in sequence.
                    chain_streams.setdefault(chain_id, []).append(line)
                else:
                    hetatm_lines_by_chain.setdefault(chain_id, []).append(line)
            elif line.startswith("TER"):
                # Rebuilt below from the cap layout.
                continue
            elif line.startswith("END"):
                end_lines.append(line)
            else:
                header_lines.append(line)

        reordered_lines = list(header_lines)

        for chain_id in sorted(chain_streams.keys()):
            stream = chain_streams[chain_id]
            for i, line in enumerate(stream):
                reordered_lines.append(line)
                # A TER goes after the LAST atom of an NME: that cap ends a
                # fragment, and whatever follows starts a new one.
                if line[17:20].strip() == "NME":
                    nxt = stream[i + 1] if i + 1 < len(stream) else None
                    same_residue = (
                        nxt is not None
                        and nxt[17:20].strip() == "NME"
                        and nxt[22:27] == line[22:27]
                    )
                    if not same_residue:
                        reordered_lines.append("TER\n")
            reordered_lines.append("TER\n")

        # Real HETATM (ligands, cofactors, waters) after every protein chain,
        # grouped by residue so tLEaP sees each as one unit.
        for chain_id in sorted(hetatm_lines_by_chain.keys()):
            hetatm_by_residue = {}
            for line in hetatm_lines_by_chain[chain_id]:
                residue_name = line[17:20].strip()
                try:
                    res_key = (int(line[22:26]), residue_name)
                except ValueError:
                    res_key = (line[22:26].strip(), residue_name)
                hetatm_by_residue.setdefault(res_key, []).append(line)

            for res_key in sorted(
                hetatm_by_residue.keys(),
                key=lambda x: (x[0] if isinstance(x[0], int) else 9999, x[1]),
            ):
                reordered_lines.extend(hetatm_by_residue[res_key])

        reordered_lines.extend(end_lines)
        return reordered_lines

    def _renumber_residues_sequential(self):
        """Renumber all residues sequentially starting at 1, returns mapping"""
        residue_mapping = {}  # (chain, old_resid, resname) -> new_resid
        seen_residues = set()
        current_id = 1
        
        # First pass: collect unique residues maintaining chain order
        ordered_residues = []
        for line in self.structure_lines:
            if line.startswith(("ATOM", "HETATM")):
                chain_id = line[21]
                try:
                    residue_id = int(line[22:26])
                    residue_name = line[17:20].strip()
                    res_key = (chain_id, residue_id, residue_name)
                    
                    if res_key not in seen_residues:
                        ordered_residues.append(res_key)
                        seen_residues.add(res_key)
                except ValueError:
                    pass
        
        # Create mapping
        for res_key in ordered_residues:
            residue_mapping[res_key] = current_id
            current_id += 1
        
        # Apply renumbering to structure lines
        new_lines = []
        for line in self.structure_lines:
            if line.startswith(("ATOM", "HETATM")):
                chain_id = line[21]
                try:
                    residue_id = int(line[22:26])
                    residue_name = line[17:20].strip()
                    res_key = (chain_id, residue_id, residue_name)
                    
                    if res_key in residue_mapping:
                        new_resid = residue_mapping[res_key]
                        # Replace residue number in PDB line
                        new_line = line[:22] + f"{new_resid:4d}" + line[26:]
                        new_lines.append(new_line)
                    else:
                        new_lines.append(line)
                except ValueError:
                    new_lines.append(line)
            else:
                new_lines.append(line)
        
        self.structure_lines = new_lines
        return residue_mapping
    
    def _update_redox_sites_with_renumbering(self, residue_mapping):
        """Update RedoxSite objects with new residue IDs via coordinate matching"""
        for site in self.redox_sites:
            # Update centers
            if site.centers:
                for center in site.centers:
                    self._update_residue_with_mapping(center, residue_mapping)
            
            # Update all atoms in the site
            for atom in site.atoms:
                self._update_residue_with_mapping(atom, residue_mapping)
            
            # Update bonds residue info
            for bond in site.bonds:
                self._update_bond_residue_info(bond, residue_mapping)
            
            # CRITICAL: Update coord_to_pdb mapping with new residue IDs
            for coord, pdb_info in site.coord_to_pdb.items():
                old_key = (pdb_info.get('chain'), pdb_info.get('resid'), pdb_info.get('resname'))
                if old_key in residue_mapping:
                    pdb_info['resid'] = residue_mapping[old_key]
                
    def _update_residue_with_mapping(self, residue, residue_mapping):
        """Update a single residue object with new ID from mapping"""
        if hasattr(residue, 'chain') and hasattr(residue, 'resid') and hasattr(residue, 'resname'):
            old_key = (residue.chain, residue.resid, residue.resname)
            if old_key in residue_mapping:
                residue.resid = residue_mapping[old_key]
    
    def _update_bond_residue_info(self, bond, residue_mapping):
        """Update bond residue info dictionaries with new residue IDs and names from current structure"""
        # Update atom1_residue_info
        if hasattr(bond, 'atom1_residue_info') and bond.atom1_residue_info:
            info = bond.atom1_residue_info
            if hasattr(bond, 'atom1_coords'):
                # Find current state of this atom in the transformed structure
                current_info = self._get_current_atom_info_from_coords(bond.atom1_coords)
                if current_info:
                    info.update(current_info)
        
        # Update atom2_residue_info  
        if hasattr(bond, 'atom2_residue_info') and bond.atom2_residue_info:
            info = bond.atom2_residue_info
            if hasattr(bond, 'atom2_coords'):
                # Find current state of this atom in the transformed structure
                current_info = self._get_current_atom_info_from_coords(bond.atom2_coords)
                if current_info:
                    info.update(current_info)
    
    def _get_current_atom_info_from_coords(self, coords):
        """Get current residue info for atom at given coordinates from transformed structure"""
        if not hasattr(self, 'structure_lines') or not self.structure_lines:
            return None
            
        target_x, target_y, target_z = coords
        tolerance = 0.001
        
        for line in self.structure_lines:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    
                    # Check if coordinates match within tolerance
                    distance = ((target_x - x)**2 + (target_y - y)**2 + (target_z - z)**2)**0.5
                    if distance <= tolerance:
                        # Extract current atom info including atom name
                        return {
                            'chain': line[21],
                            'resname': line[17:20].strip(),
                            'resid': int(line[22:26]),
                            'atom_name': line[12:16].strip(),
                        }
                except (ValueError, IndexError):
                    continue
        return None

    def _update_redox_sites_from_existing_pdb(self, pdb_file: str, sites_to_update: list = None) -> bool:
        """
        Update RedoxSite definitions from an existing PDB file using coordinate-based matching.

        This is used when microstate PDB files already exist and transformation was skipped.
        The method reads the PDB file and updates residue names, IDs, and atom names in the
        RedoxSite objects based on coordinate matching.

        Args:
            pdb_file: Path to the existing PDB file
            sites_to_update: List of RedoxSite objects to update. If None, uses self.redox_sites

        Returns:
            True if successful, False otherwise
        """
        # Delegate to the module-level, reusable coordinate-based sync so the
        # modified-AA parameterizer and any other caller share one
        # implementation. Behavior is identical (coords are the join key).
        from proprep.structure_prep.comprehensive_redox_detector import (
            sync_redox_sites_from_pdb,
        )
        sites = sites_to_update if sites_to_update is not None else self.redox_sites
        return sync_redox_sites_from_pdb(pdb_file, sites)

    def _export_updated_redox_sites(self, pdb_file):
        """
        Export updated RedoxSites to JSON file with renumbered residues.

        This ensures that if the user loads the transformed structure in a new session,
        they can also load the updated RedoxSites JSON that matches the new residue numbering.
        """
        # Embed the transformer assignments so they survive a ProPrep restart.
        # (Transformers are parameterization recipes — not all transformed sites
        # are redox-active; e.g. NADPH/FMN/FAD/Ca²⁺ carry a single fixed state.
        # The Topology Generator reads this back to recover the site→forcefield
        # mapping when the session-only workspace copy is gone.)
        output_file = export_redox_sites_to_json(
            self.redox_sites, pdb_file,
            transformer_info=self._build_transformer_info())

        if not getattr(self, '_batch_mode', False):
            self.console.print(f"[green]✅ Exported updated RedoxSites to {output_file.name}[/green]")
            self.console.print(f"[grey50]This JSON file has residue IDs matching the transformed structure[/grey50]")

    def _validate_workspace_consistency(self):
        """Comprehensive validation that workspace RedoxSites match the final PDB structure"""
        if not self.processor or not self.processor.workspace:
            return
            
        workspace = self.processor.workspace
        
        # Get workspace data
        workspace_sites = workspace.get("detected_redox_sites", [])
        pdb_file_path = workspace.get("transformed_pdb_file")
        
        if not workspace_sites or not pdb_file_path:
            self.console.print("[yellow]⚠️ No workspace data to validate[/yellow]")
            return
            
        self.console.print("\n🔍 Validating workspace consistency...")
        self.console.print("[grey50]Parsing PDB structure and validating all atoms (may take time for large structures)...[/grey50]")

        # Load and parse the final PDB file
        try:
            with open(pdb_file_path, 'r') as f:
                pdb_lines = f.readlines()
        except Exception as e:
            self.console.print(f"[red]❌ Failed to load PDB file: {e}[/red]")
            return
            
        # Parse PDB structure into coordinate lookup
        pdb_atoms = {}  # (x,y,z) -> {chain, resname, resid, atom_name, element}
        
        for line in pdb_lines:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    chain = line[21]
                    resname = line[17:20].strip()
                    resid = int(line[22:26])
                    atom_name = line[12:16].strip()
                    element = line[76:78].strip() or atom_name[0]
                    
                    x = float(line[30:38])
                    y = float(line[38:46]) 
                    z = float(line[46:54])
                    coords = (x, y, z)  # Store exact coordinates for tolerance-based matching
                    
                    pdb_atoms[coords] = {
                        'chain': chain,
                        'resname': resname,
                        'resid': resid,
                        'atom_name': atom_name,
                        'element': element
                    }
                except (ValueError, IndexError):
                    continue
        
        # Validation results
        validation_errors = []
        atoms_validated = 0
        centers_validated = 0
        bonds_validated = 0

        # Validate each RedoxSite with progress bar for large structures
        if len(workspace_sites) >= 20:
            from rich.progress import BarColumn, TaskProgressColumn, TimeRemainingColumn
            with Progress(
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TaskProgressColumn(),
                TimeRemainingColumn(),
                console=self.console,
                transient=False
            ) as progress:
                task = progress.add_task(f"Validating {len(workspace_sites)} sites...", total=len(workspace_sites))

                for site_idx, site in enumerate(workspace_sites, 1):
                    site_errors, site_atoms, site_centers, site_bonds = self._validate_single_site(
                        site, site_idx, pdb_atoms
                    )
                    validation_errors.extend(site_errors)
                    atoms_validated += site_atoms
                    centers_validated += site_centers
                    bonds_validated += site_bonds
                    progress.update(task, advance=1)
        else:
            for site_idx, site in enumerate(workspace_sites, 1):
                site_errors, site_atoms, site_centers, site_bonds = self._validate_single_site(
                    site, site_idx, pdb_atoms
                )
                validation_errors.extend(site_errors)
                atoms_validated += site_atoms
                centers_validated += site_centers
                bonds_validated += site_bonds
        
        # Report results
        if validation_errors:
            self.console.print(f"[red]❌ Validation failed with {len(validation_errors)} errors:[/red]")
            for error in validation_errors[:10]:  # Show first 10 errors
                self.console.print(f"  • {error}")
            if len(validation_errors) > 10:
                self.console.print(f"  • ... and {len(validation_errors) - 10} more errors")
        else:
            self.console.print(f"[green]✅ Workspace validation passed[/green]")
            self.console.print(f"  • {centers_validated} centers validated")
            self.console.print(f"  • {atoms_validated} atoms validated")
            self.console.print(f"  • {bonds_validated} bonds validated")
            self.console.print(f"  • Total components: {centers_validated + atoms_validated}")
            self.console.print(f"  • All RedoxSite definitions match final PDB structure")

    def _validate_single_site(self, site, site_idx, pdb_atoms):
        """Validate a single RedoxSite and return errors and counts"""
        site_errors = []
        atoms_count = 0
        centers_count = 0
        bonds_count = 0

        # Validate centers
        if hasattr(site, 'centers') and site.centers:
            for center_idx, center in enumerate(site.centers):
                center_errors = self._validate_site_component(center, pdb_atoms, f"Site {site_idx} Center {center_idx}")
                site_errors.extend(center_errors)
                if not center_errors:
                    centers_count += 1

        # Validate all atoms
        if hasattr(site, 'atoms') and site.atoms:
            for atom_idx, atom in enumerate(site.atoms):
                atom_errors = self._validate_site_component(atom, pdb_atoms, f"Site {site_idx} Atom {atom_idx}")
                site_errors.extend(atom_errors)
                if not atom_errors:
                    atoms_count += 1

        # Validate bonds
        if hasattr(site, 'bonds') and site.bonds:
            for bond_idx, bond in enumerate(site.bonds):
                bond_errors = self._validate_bond_consistency(bond, pdb_atoms, f"Site {site_idx} Bond {bond_idx}")
                site_errors.extend(bond_errors)
                if not bond_errors:
                    bonds_count += 1

        return site_errors, atoms_count, centers_count, bonds_count

    def _validate_site_component(self, component, pdb_atoms, context):
        """Validate a single site component (center, atom) against PDB structure"""
        errors = []
        
        if not hasattr(component, 'coords'):
            errors.append(f"{context}: Missing coordinates")
            return errors
            
        coords = component.coords
        if isinstance(coords, (list, tuple)) and len(coords) == 3:
            coords = (round(coords[0], 3), round(coords[1], 3), round(coords[2], 3))
        else:
            errors.append(f"{context}: Invalid coordinates format")
            return errors

        # Skip coordinate validation for whole-residue centers (centroids)
        # Centroids are calculated positions, not actual atom coordinates
        if hasattr(component, 'atom_name') and component.atom_name is None:
            # This is a centroid - skip PDB coordinate lookup
            # Just validate that the residue exists in the structure
            return errors

        # Find matching coordinates in PDB using tolerance-based matching
        pdb_data = self._find_matching_pdb_atom(coords, pdb_atoms)
        if not pdb_data:
            errors.append(f"{context}: Coordinates {coords} not found in PDB structure")
            return errors
        
        # Validate chain
        if hasattr(component, 'chain') and component.chain != pdb_data['chain']:
            errors.append(f"{context}: Chain mismatch - RedoxSite: {component.chain}, PDB: {pdb_data['chain']}")
            
        # Validate resname
        if hasattr(component, 'resname') and component.resname != pdb_data['resname']:
            errors.append(f"{context}: Resname mismatch - RedoxSite: {component.resname}, PDB: {pdb_data['resname']}")
            
        # Validate resid
        if hasattr(component, 'resid') and component.resid != pdb_data['resid']:
            errors.append(f"{context}: Resid mismatch - RedoxSite: {component.resid}, PDB: {pdb_data['resid']}")
            
        # Validate atom name
        if hasattr(component, 'atom_name') and component.atom_name != pdb_data['atom_name']:
            errors.append(f"{context}: Atom name mismatch - RedoxSite: {component.atom_name}, PDB: {pdb_data['atom_name']}")
            
        # Validate element (skip if RedoxSite element is None as it may not be populated).
        # Compare case-insensitively: the PDB writer emits the standard element symbol
        # (e.g. "Ca", "Fe"), while RedoxSite stores it uppercased ("CA", "FE"). Without
        # normalization, every two-letter element produces a spurious mismatch.
        if (hasattr(component, 'element') and component.element is not None
                and component.element.strip().upper() != pdb_data['element'].strip().upper()):
            errors.append(f"{context}: Element mismatch - RedoxSite: {component.element}, PDB: {pdb_data['element']}")
            
        return errors
    
    def _validate_bond_consistency(self, bond, pdb_atoms, context):
        """Validate bond residue info against PDB structure"""
        errors = []
        
        # Validate atom1 coordinates exist
        if hasattr(bond, 'atom1_coords'):
            coords1 = bond.atom1_coords
            if isinstance(coords1, (list, tuple)) and len(coords1) == 3:
                pdb_data1 = self._find_matching_pdb_atom(coords1, pdb_atoms)
                if not pdb_data1:
                    errors.append(f"{context}: Atom1 coordinates {coords1} not found in PDB")
                elif hasattr(bond, 'atom1_residue_info') and bond.atom1_residue_info:
                    # Validate atom1 residue info
                    pdb_data = pdb_data1
                    info = bond.atom1_residue_info
                    
                    if 'resid' in info and info['resid'] != pdb_data['resid']:
                        errors.append(f"{context}: Atom1 resid mismatch - Bond: {info['resid']}, PDB: {pdb_data['resid']}")
                    if 'chain' in info and info['chain'] != pdb_data['chain']:
                        errors.append(f"{context}: Atom1 chain mismatch - Bond: {info['chain']}, PDB: {pdb_data['chain']}")
                    if 'resname' in info and info['resname'] != pdb_data['resname']:
                        errors.append(f"{context}: Atom1 resname mismatch - Bond: {info['resname']}, PDB: {pdb_data['resname']}")
        
        # Validate atom2 coordinates exist
        if hasattr(bond, 'atom2_coords'):
            coords2 = bond.atom2_coords
            if isinstance(coords2, (list, tuple)) and len(coords2) == 3:
                pdb_data2 = self._find_matching_pdb_atom(coords2, pdb_atoms)
                if not pdb_data2:
                    errors.append(f"{context}: Atom2 coordinates {coords2} not found in PDB")
                elif hasattr(bond, 'atom2_residue_info') and bond.atom2_residue_info:
                    # Validate atom2 residue info
                    pdb_data = pdb_data2
                    info = bond.atom2_residue_info
                    
                    if 'resid' in info and info['resid'] != pdb_data['resid']:
                        errors.append(f"{context}: Atom2 resid mismatch - Bond: {info['resid']}, PDB: {pdb_data['resid']}")
                    if 'chain' in info and info['chain'] != pdb_data['chain']:
                        errors.append(f"{context}: Atom2 chain mismatch - Bond: {info['chain']}, PDB: {pdb_data['chain']}")
                    if 'resname' in info and info['resname'] != pdb_data['resname']:
                        errors.append(f"{context}: Atom2 resname mismatch - Bond: {info['resname']}, PDB: {pdb_data['resname']}")
        
        return errors
    
    def _find_matching_pdb_atom(self, target_coords, pdb_atoms, tolerance=0.001):
        """Find PDB atom with coordinates matching target within tolerance"""
        target_x, target_y, target_z = target_coords
        
        for pdb_coords, pdb_data in pdb_atoms.items():
            pdb_x, pdb_y, pdb_z = pdb_coords
            distance = ((target_x - pdb_x)**2 + (target_y - pdb_y)**2 + (target_z - pdb_z)**2)**0.5
            if distance <= tolerance:
                return pdb_data
        return None

    def _store_transformer_info_in_workspace(self):
        """Store transformer information in workspace for single state template to use"""
        if not self.processor:
            return

        workspace = self.processor._get_workspace()
        transformer_info = self._build_transformer_info()
        workspace.update({'transformer_info': transformer_info})
        self.console.print(f"[green]Stored transformer info for {len(transformer_info)} sites in workspace[/green]")

    def _build_transformer_info(self):
        """Build the per-site transformer_info list (transformer type, redox/spin
        state, forcefield path, atom types) from the current site→transformer
        assignments. This is the same structure stored under the workspace key
        ``transformer_info`` and — so the assignments survive a ProPrep restart —
        embedded in the exported redox-sites JSON (see export_redox_sites_to_json).
        """
        # Build transformer info for each redox site
        transformer_info = []

        for site in self.redox_sites:
            # Get primary center (first one if multiple)
            primary_center = site.centers[0] if site.centers else None
            
            site_info = {
                'site_id': site.site_id,
                'center_type': primary_center.center_type.value if primary_center else 'unknown',
                'metal_element': primary_center.element if primary_center else '',
                'residue_name': primary_center.resname if primary_center else '',
                'residue_id': primary_center.resid if primary_center else '',
                'chain': primary_center.chain if primary_center else '',
                'atom_name': primary_center.atom_name if primary_center else ''
            }
            
            # Add transformer information if this site has a selected transformer
            if site.site_id in self.selected_transformers:
                transformer_class_name = self.selected_transformers[site.site_id]
                parameters = self.transformation_parameters.get(site.site_id, {})

                # Extract common parameters
                redox_state = parameters.get('redox_state', 'unknown')
                spin_state = parameters.get('spin_state', 'unknown')

                # Get atom types from metadata (similar to metallo_prep)
                atom_types = self._get_atom_types_from_metadata(transformer_class_name, redox_state, spin_state)

                # Thread the transformer's FORCEFIELD_PATH so downstream code
                # (topology generator's prereq-leaprc emission, Layer 1/2/3)
                # can look up metadata without re-resolving by class name.
                cofactor_path = None
                forcefield_set = None
                try:
                    from proprep.redoxsite_prep.transformation.redox_transformer_framework import redox_transformer_registry
                    transformer_class = redox_transformer_registry.get_transformer(transformer_class_name)
                    if transformer_class is not None:
                        cofactor_path = getattr(transformer_class, "FORCEFIELD_PATH", None)
                        # Resolve the forcefield set the parameters imply (e.g. the
                        # pH-treatment fork's fixed-pH vs constant-pH set) so the
                        # Topology Generator can honor the Stage-1 choice: pre-select
                        # the matching set and apply its per-set leaprc prerequisites.
                        if hasattr(transformer_class, "select_forcefield_set_name"):
                            forcefield_set = transformer_class.select_forcefield_set_name(parameters)
                except Exception:
                    pass

                site_info.update({
                    'has_transformer': True,
                    'transformer_type': transformer_class_name,
                    'cofactor_path': cofactor_path,
                    'redox_state': redox_state,
                    'spin_state': spin_state,
                    'parameters': parameters,
                    'atom_types': atom_types,
                    # pH-treatment dimension + resolved set (None for cofactors
                    # without a fork) — consumed by the topology generator.
                    'ph_treatment': parameters.get('ph_treatment'),
                    'forcefield_set': forcefield_set,
                })
            else:
                site_info.update({
                    'has_transformer': False,
                    'transformer_type': None,
                    'cofactor_path': None,
                    'redox_state': None,
                    'spin_state': None,
                    'parameters': {},
                    'atom_types': []
                })
            
            transformer_info.append(site_info)

        return transformer_info

    def _resolve_state_params(self, site_id: str, redox_state: str, spin_state: str) -> Dict[str, Any]:
        """The full transformer parameter dict for a (site, redox/spin) state in
        batch mode: redox/spin plus the propionate/pH treatment configured in
        Step 6 (empty when the site's transformer has no pH fork)."""
        params = {'redox_state': redox_state, 'spin_state': spin_state}
        params.update(
            self.site_state_protonation.get(site_id, {}).get(f"{redox_state}_{spin_state}", {})
        )
        return params

    def _build_batch_transformer_info(self, microstates: List[Dict]) -> List[Dict[str, Any]]:
        """Per-(site, state) transformer_info for the batch flow.

        Unlike the single-path :py:meth:`_build_transformer_info` (one entry per
        site, one state), a batch microstate set exercises several states per
        site, and — since Step 6 binds the propionate treatment per redox state —
        each (site, state) can imply a different forcefield set. The Topology
        Generator's FF-set picker matches on (transformer_type, redox, spin) and
        reads ``ph_treatment`` / ``forcefield_set`` off these entries
        (``_preferred_ph_treatment_for_combo`` / ``_preferred_ff_set_for_combo``),
        so emitting one entry per used (site, state) lets it pre-select and filter
        to the treatment chosen here instead of offering every set.
        """
        from proprep.redoxsite_prep.transformation.redox_transformer_framework import redox_transformer_registry

        # site_id -> set of (redox, spin) actually used by the selected microstates
        used: Dict[str, set] = {}
        for ms in microstates:
            for sid, st in ms.items():
                used.setdefault(sid, set()).add((st['redox_state'], st['spin_state']))

        info: List[Dict[str, Any]] = []
        # The Topology Generator's combo readers key on
        # (transformer_type, redox, spin, ph_treatment, forcefield_set); several
        # identical sites (e.g. 4 equivalent hemes) collapse to the same combo,
        # so dedup to keep the stored info tidy without losing any distinction
        # the readers or prereq collectors care about.
        seen: set = set()
        for site in self.redox_sites:
            sid = getattr(site, 'site_id', None)
            if sid not in used:
                continue
            tname = self.site_transformer_assignments.get(sid)
            if not tname:
                continue
            tclass = redox_transformer_registry.get_transformer(tname)
            primary = site.centers[0] if getattr(site, 'centers', None) else None
            for redox, spin in sorted(used[sid]):
                params = self._resolve_state_params(sid, redox, spin)
                ff_set = None
                if tclass and hasattr(tclass, 'select_forcefield_set_name'):
                    try:
                        ff_set = tclass.select_forcefield_set_name(params)
                    except Exception:
                        pass
                dedup_key = (tname, redox, spin, params.get('ph_treatment'), ff_set)
                if dedup_key in seen:
                    continue
                seen.add(dedup_key)
                info.append({
                    'site_id': sid,
                    'has_transformer': True,
                    'transformer_type': tname,
                    'cofactor_path': getattr(tclass, 'FORCEFIELD_PATH', None) if tclass else None,
                    'redox_state': redox,
                    'spin_state': spin,
                    'parameters': params,
                    'ph_treatment': params.get('ph_treatment'),
                    'forcefield_set': ff_set,
                    'residue_name': primary.resname if primary else '',
                    'chain': primary.chain if primary else '',
                    'residue_id': primary.resid if primary else '',
                })
        return info

    def _store_batch_transformer_info_in_workspace(self, microstates: List[Dict]):
        """Persist batch per-(site, state) transformer_info so the Topology
        Generator can honor the Step-6 propionate/pH choice (pre-select + filter
        the FF-set picker). Best-effort: a failure here must not block generation.
        """
        try:
            info = self._build_batch_transformer_info(microstates)
            if info:
                self.processor.workspace.update({'transformer_info': info})
        except Exception as e:
            logger.warning(f"Could not store batch transformer_info: {e}")

    def _get_atom_types_from_metadata(self, transformer_class_name, redox_state, spin_state):
        """Get atom types from transformer metadata file"""
        atom_types = []
        try:
            # Import the transformer registry from the correct location
            from proprep.redoxsite_prep.transformation.redox_transformer_framework import redox_transformer_registry
            transformer_class = redox_transformer_registry.get_transformer(transformer_class_name)
            
            if transformer_class and hasattr(transformer_class, 'FORCEFIELD_PATH') and transformer_class.FORCEFIELD_PATH is not None:
                from proprep.forcefield_params.loader import (
                    get_forcefield_base_path, get_user_forcefield_base_path,
                )

                # Check the bundled library first, then the user library
                # (~/.proprep) — auto-emitted reuse transformers deposit there.
                metadata_path = get_forcefield_base_path() / transformer_class.FORCEFIELD_PATH / 'metadata.json'
                if not metadata_path.exists():
                    metadata_path = get_user_forcefield_base_path() / transformer_class.FORCEFIELD_PATH / 'metadata.json'

                if metadata_path.exists():
                    import json
                    with open(metadata_path, 'r') as f:
                        metadata = json.load(f)
                    
                    # Navigate: redox_states -> redox_state -> spin_states -> spin_state -> atom_types
                    redox_data = metadata.get('redox_states', {}).get(redox_state, {})
                    if redox_data:
                        spin_data = redox_data.get('spin_states', {}).get(spin_state, {})
                        if spin_data:
                            atom_types = spin_data.get('atom_types', [])
                            
                            if atom_types:
                                self.console.print(f"[grey50]Found {len(atom_types)} atom types for {transformer_class_name} from metadata[/grey50]")
                else:
                    self.console.print(f"[yellow]Metadata file not found: {metadata_path}[/yellow]")
            else:
                # Check if transformer class exists but has no FORCEFIELD_PATH
                if transformer_class and transformer_class.FORCEFIELD_PATH is None:
                    self.console.print(f"[grey50]No forcefield metadata needed for {transformer_class_name}[/grey50]")
                else:
                    self.console.print(f"[yellow]No transformer class or FORCEFIELD_PATH for {transformer_class_name}[/yellow]")
                        
        except Exception as e:
            self.console.print(f"[yellow]Could not get atom types for {transformer_class_name}: {e}[/yellow]")
        
        return atom_types