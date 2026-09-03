"""
ONIOM QM/MM Preparator Module

Generates Gaussian ONIOM QM/MM input files from RedoxSite objects.
Provides a user-friendly interface for configuring QM/MM calculations with
proper layer assignment, link atom placement, and force field parameterization.
"""

import logging
import os
from typing import Any, Dict, Iterable, List, Optional, Tuple
from pathlib import Path

from rich.panel import Panel
from proprep.utils.prompts import prompt_with_context, confirm_with_context, int_prompt_with_context
from rich.table import Table

from proprep.utils.module_registry import ProcessingModule
from proprep.utils.workflow_checklist import WorkflowChecklist, WorkflowStep
from proprep.oniom_prep.oniom_preparator import ONIOMPreparator
from proprep.oniom_prep.data_structures import ONIOMLayer
from proprep.structure_prep.comprehensive_redox_detector import RedoxSite

# Setup logging
logger = logging.getLogger(__name__)


# =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
# Workflow Steps
# =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#

ONIOM_WORKFLOW_STEPS = [
    WorkflowStep(
        id="oniom-1", name="Select RedoxSite",
        description="Choose which redox site(s) to include in the QM region",
        handler="_checklist_select_redox_site",
        section="Setup",
    ),
    WorkflowStep(
        id="oniom-2", name="Configure ONIOM Layers",
        description="Assign residues to HIGH/MEDIUM/LOW layers and set freezing options",
        handler="_checklist_configure_layers",
        section="Setup",
        dependencies=["oniom-1"],
    ),
    WorkflowStep(
        id="oniom-3", name="Configure MM Settings & Atom Typing",
        description="Select force field, run atom typing, analyze charges and boundaries",
        handler="_checklist_configure_mm",
        section="Configuration",
        dependencies=["oniom-1", "oniom-2"],
    ),
    WorkflowStep(
        id="oniom-4", name="Configure QM Settings",
        description="Set QM functional, basis set, charge, and multiplicity",
        handler="_checklist_configure_qm",
        section="Configuration",
        dependencies=["oniom-3"],
    ),
    WorkflowStep(
        id="oniom-5", name="Configure Job Settings",
        description="Set job type, processors, memory, and additional keywords",
        handler="_checklist_configure_job",
        section="Configuration",
        optional=True,
    ),
    WorkflowStep(
        id="oniom-6", name="Prepare & Validate ONIOM Setup",
        description="Build ONIOM configuration, assign layers, place link atoms, validate",
        handler="_checklist_prepare_and_validate",
        section="Generation",
        dependencies=["oniom-2", "oniom-3", "oniom-4"],
    ),
    WorkflowStep(
        id="oniom-7", name="Write Gaussian Input File",
        description="Write production or diagnostic ONIOM input (.com) file",
        handler="_checklist_write_input",
        section="Generation",
        dependencies=["oniom-6"],
    ),
]


def suggested_model_charge(charges: Dict[int, float],
                           whole_residues: Iterable[Iterable[int]],
                           trimmed_fragments: Iterable[Iterable[int]]) -> int:
    """Formal charge of an ONIOM model system from MM partial charges.

    ``whole_residues`` and ``trimmed_fragments`` are groups of atom indices;
    each group is one fragment whose partial charges are summed and rounded
    to its own integer, then the integers are added. An AMBER residue sums to
    an exact integer, a side chain cut at CA-CB does not (ASP -0.86, GLU
    -0.88, HIP +0.94, LYS +1.03), so rounding must happen per fragment:
    four trimmed carboxylates are -4, not round(-3.43) = -3.
    """
    total = 0
    for group in list(whole_residues) + list(trimmed_fragments):
        total += round(sum(charges.get(idx, 0.0) for idx in group))
    return total


class ONIOMQMMMPreparator(ProcessingModule):
    """Module for generating Gaussian ONIOM QM/MM input files.

    Typically instantiated and driven by the unified QM/MM Preparator,
    which handles topology loading and frame extraction upstream.
    """

    NAME = "ONIOM QM/MM Preparator"
    DESCRIPTION = "Generate Gaussian ONIOM QM/MM input files from RedoxSite"
    CATEGORY = "preparation"
    VERSION = "1.0.0"

    def __init__(self):
        """Initialize the ONIOM QM/MM Preparator module"""
        self.processor = None

        # Configuration storage
        self.selected_redox_site: Optional[RedoxSite] = None
        self.selected_site_index: Optional[int] = None
        self.selected_site_indices: List[int] = []  # 0-based indices of selected sites

        # Layer configuration (parmed residue indices)
        self.high_residue_indices: List[int] = []
        self.medium_residue_indices: List[int] = []
        self.n_layers: int = 2
        self.freeze_distance_cutoff: Optional[float] = None
        self.frozen_residue_indices: List[int] = []

        # QM settings
        self.qm_functional: str = "B3LYP"
        self.qm_basis_set: str = "6-31G*"
        self.qm_charge: int = 0
        self.qm_multiplicity: int = 1

        # MM settings
        self.mm_forcefield: str = "AMBER"
        self.force_field_name: str = "prmtop"

        # Individual atoms forced into HIGH (for CA-CB sidechain-only cuts)
        self.high_atom_indices: set = set()

        # Atom types and charges keyed by parmed atom index (0-based)
        self.atom_types: Dict[int, str] = {}
        self.charges: Dict[int, float] = {}
        self.layer_membership: Dict[int, ONIOMLayer] = {}
        self.boundary_atom_indices: List[int] = []
        self.high_layer_charge: float = 0.0
        self.low_layer_charge: float = 0.0
        self.medium_layer_charge: float = 0.0
        self.suggested_qm_charge: int = 0
        self.atom_typing_complete: bool = False

        # Medium layer settings (3-layer only)
        self.medium_method: str = "HF/3-21G"
        self.medium_charge: int = 0
        self.medium_multiplicity: int = 1

        # Job settings
        self.job_type: str = "Opt"
        self.n_processors: int = 4
        self.memory_gb: int = 8
        self.additional_keywords: str = ""

        # ONIOMPreparator instance (created during preparation)
        self.preparator: Optional[ONIOMPreparator] = None

        # parmed structure (set by set_frame_data or _extract_types_from_prmtop)
        self._parm = None

        # Upstream frame data (set by the unified QM/MM Preparator)
        self._using_upstream_frames: bool = False
        self._upstream_prmtop: str = ""
        self._upstream_pdb_paths: List[str] = []
        self._upstream_labels: List[str] = []

    def set_frame_data(self, prmtop_path: str, pdb_paths: List[str],
                       labels: List[str], parm=None):
        """Accept frame data from the unified QM/MM Preparator."""
        self._upstream_prmtop = prmtop_path
        self._upstream_pdb_paths = list(pdb_paths)
        self._upstream_labels = list(labels)
        self._using_upstream_frames = True
        if parm is not None:
            self._parm = parm

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Helper methods to centralize interactions with the workspace

    def get_workspace(self):
        """Get the current workspace object"""
        return self.processor.workspace

    def get_from_workspace(self, key, default=None):
        """Get values from the processor's workspace"""
        return self.processor.workspace.get(key, default)

    def update_workspace(self, key, value):
        """Update the processor's workspace"""
        self.processor.workspace.set(key, value)

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Module Menu

    def get_workspace_outputs(self) -> List[str]:
        """Get workspace outputs"""
        return ["oniom_setup", "oniom_input_file"]

    def get_menu_options(self) -> Dict[str, str]:
        """Get module menu options"""
        return {
            "run_workflow": "Run ONIOM QM/MM workflow (guided checklist)",
            "view_config": "View current ONIOM configuration",
            "reset": "Reset configuration to defaults",
        }

    def get_menu_suggestion(self, workspace):
        """Get a suggestion for the next recommended action."""
        redox_sites = workspace.get("detected_redox_sites", [])
        if not redox_sites:
            return "No redox sites detected. Run RedoxSite Detector first."

        oniom_input_written = workspace.get("oniom_input_file") is not None
        if oniom_input_written:
            return "✓ ONIOM input file written. Setup complete"

        return "Run the ONIOM workflow (option 1) to prepare Gaussian input files"

    def handle_menu_option(self, option: str) -> bool:
        """Handle a menu option selection"""
        if option == "run_workflow":
            return self._run_oniom_workflow()
        elif option == "view_config":
            return self.view_current_setup()
        elif option == "reset":
            return self.reset_configuration()
        return False

    def _run_oniom_workflow(self) -> bool:
        """Run the ONIOM preparation workflow via interactive checklist."""
        workspace = self.get_workspace()
        base_dir = workspace.get('working_directory', os.getcwd())
        state_dir = Path(base_dir) / "oniom_inputs"
        state_dir.mkdir(exist_ok=True)

        checklist = WorkflowChecklist(
            steps=ONIOM_WORKFLOW_STEPS,
            executor=self,
            processor=self.processor,
            workflow_name="ONIOM QM/MM Preparation",
            console=self.processor.console,
            state_dir=state_dir,
        )
        return checklist.run()

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # WorkflowChecklist Handler Methods
    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#

    def _checklist_select_redox_site(self):
        """Handler: select RedoxSite(s) for ONIOM calculation."""
        success = self.select_redox_site()
        if not success:
            raise RuntimeError("RedoxSite selection cancelled or failed")

        n_atoms = len(self.selected_redox_site.atoms)
        n_sites = len(self.selected_site_indices)
        if n_sites == 1:
            return {'summary': f"Selected RedoxSite {self.selected_site_indices[0] + 1} ({n_atoms} atoms)"}
        else:
            labels = ", ".join(str(i + 1) for i in self.selected_site_indices)
            return {'summary': f"Selected RedoxSites {labels} (merged: {n_atoms} atoms)"}

    def _checklist_configure_layers(self):
        """Handler: configure ONIOM layer assignments."""
        success = self.configure_layers()
        if not success:
            raise RuntimeError("Layer configuration cancelled or failed")

        n_high = len(self.high_residue_indices)
        n_medium = len(self.medium_residue_indices) if self.n_layers == 3 else 0
        summary = f"{self.n_layers}-layer ONIOM, {n_high} HIGH residues"
        if n_medium:
            summary += f", {n_medium} MEDIUM residues"
        if self.freeze_distance_cutoff:
            summary += f", freeze beyond {self.freeze_distance_cutoff}Å"
        return {'summary': summary}

    def _checklist_configure_mm(self):
        """Handler: configure MM settings and run atom typing."""
        success = self.configure_mm_settings()
        if not success:
            raise RuntimeError("MM configuration cancelled or failed")

        # Offer to visualize layers
        if confirm_with_context(
            self.processor, "\nView ONIOM layers in 3D structure viewer?",
            default=False,
            module="ONIOM QM/MM Preparator",
            description="View ONIOM layer assignments in structure viewer",
        ):
            self.launch_oniom_viewer()

        n_types = len(set(self.atom_types.values()))
        return {'summary': f"AMBER atom typing complete ({n_types} unique types)"}

    def _checklist_configure_qm(self):
        """Handler: configure QM calculation settings."""
        success = self.configure_qm_settings()
        if not success:
            raise RuntimeError("QM configuration cancelled or failed")
        return {'summary': f"{self.qm_functional}/{self.qm_basis_set}, charge={self.qm_charge}, mult={self.qm_multiplicity}"}

    def _checklist_configure_job(self):
        """Handler: configure job execution settings."""
        success = self.configure_job_settings()
        if not success:
            raise RuntimeError("Job configuration cancelled or failed")
        summary = f"{self.job_type}, {self.n_processors} procs, {self.memory_gb}GB"
        if self.additional_keywords:
            summary += f", keywords: {self.additional_keywords}"
        return {'summary': summary}

    def _checklist_prepare_and_validate(self):
        """Handler: prepare ONIOM setup and show validation."""
        success = self.prepare_oniom_setup()
        if not success:
            raise RuntimeError("ONIOM preparation failed")

        # Automatically show validation report
        self.processor.console.print("\n[grey50]Showing validation report...[/grey50]")
        self.view_validation_report()

        return {'summary': "ONIOM setup prepared and validated"}

    def _checklist_write_input(self):
        """Handler: write Gaussian ONIOM input file (production or diagnostic)."""
        console = self.processor.console

        console.print("\n[bold cyan]Write Gaussian Input[/bold cyan]")
        console.print("1. Write production input file", highlight=False)
        console.print(
            "2. Write diagnostic file (ONIOM=OnlyInputFiles)",
            highlight=False
        )
        console.print(
            "   [grey50]Gaussian will write separate .com files for each sub-calculation[/grey50]",
            highlight=False
        )
        console.print(
            "   [grey50]so you can verify layer assignments. No calculation is run.[/grey50]",
            highlight=False
        )

        write_choice = prompt_with_context(
            self.processor, "Select option",
            choices=["1", "2"], default="1",
            module="ONIOM QM/MM Preparator",
            description="Write production or diagnostic input",
            options_map={"1": "Production input", "2": "Diagnostic file"},
        )

        if write_choice == "1":
            success = self.write_oniom_input()
            if not success:
                raise RuntimeError("Failed to write ONIOM input file")
            summary_msg = "Production .com file written"
        else:
            success = self.write_diagnostic_input()
            if not success:
                raise RuntimeError("Failed to write diagnostic file")
            summary_msg = "Diagnostic .com file written (OnlyInputFiles)"

        # Automatically write summary report alongside the input file
        self.write_summary_report()

        return {'summary': summary_msg}

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # RedoxSite Selection

    def select_redox_site(self) -> bool:
        """
        Select a RedoxSite from detected sites in workspace.

        Returns:
            True if site selected successfully, False otherwise
        """
        redox_sites = self.get_from_workspace("detected_redox_sites")

        if not redox_sites:
            self.processor.console.print("\n[red]No RedoxSite objects found in workspace.[/red]")
            self.processor.console.print("Please run 'Comprehensive RedoxSite Detector' first.")
            return False

        self.processor.console.print(f"\n[bold cyan]Found {len(redox_sites)} RedoxSite(s)[/bold cyan]\n")

        # Display available sites
        table = Table(title="Available RedoxSites")
        table.add_column("Index", style="cyan", width=6)
        table.add_column("Redox Centers", style="green")
        table.add_column("Residues", style="yellow")
        table.add_column("Atoms", style="magenta")
        table.add_column("Bonds", style="blue")

        from proprep.structure_prep.comprehensive_redox_detector import CenterType

        for i, site in enumerate(redox_sites, 1):
            # Describe the redox centers in this site
            center_description = []

            # Check if centers attribute exists and what it contains
            if hasattr(site, 'centers') and site.centers:
                # Group centers by type - compare by enum value string to avoid enum identity issues
                metal_centers = [c for c in site.centers if c.center_type.value == "metal_ion"]
                organometallic_centers = [c for c in site.centers if c.center_type.value == "organometallic_cofactor"]
                organic_centers = [c for c in site.centers if c.center_type.value == "organic_cofactor"]
                amino_centers = [c for c in site.centers if c.center_type.value == "redox_amino_acid"]

                if metal_centers:
                    elements = ", ".join(sorted(set(c.element for c in metal_centers if c.element)))
                    center_description.append(f"{len(metal_centers)} metal ({elements})")
                if organometallic_centers:
                    resnames = ", ".join(sorted(set(c.resname for c in organometallic_centers)))
                    center_description.append(f"{len(organometallic_centers)} organometallic ({resnames})")
                if organic_centers:
                    resnames = ", ".join(sorted(set(c.resname for c in organic_centers)))
                    center_description.append(f"{len(organic_centers)} cofactor ({resnames})")
                if amino_centers:
                    resnames = ", ".join(sorted(set(c.resname for c in amino_centers)))
                    center_description.append(f"{len(amino_centers)} amino ({resnames})")

            center_str = "; ".join(center_description) if center_description else "No centers"

            # Count unique residues
            residues = set((atom.chain, atom.resid, atom.insertion_code) for atom in site.atoms)

            table.add_row(
                str(i),
                center_str,
                str(len(residues)),
                str(len(site.atoms)),
                str(len(site.bonds))
            )

        self.processor.console.print(table)

        # Hook A: highlight every redox site listed in the table
        # palette-coloured by row number (the same number the user
        # types in below). Lets the user confirm site identity
        # spatially before picking. Best-effort — the user is
        # configuring a QM region either way.
        try:
            from proprep.structure_prep.viewer_coordinator import (
                viewer as _viewer,
            )
            pdb_file = self._upstream_pdb_paths[0] if self._upstream_pdb_paths else None
            if pdb_file:
                # Auto-fired workflow hook (palette-colour every candidate
                # site before the user picks). No force= so it stays silent
                # in CLI when no viewer is open; pushes live to one that is.
                _viewer.show_structure(pdb_file)
            applied = []
            for i, site in enumerate(redox_sites, 1):
                pairs = set()
                for atom in getattr(site, 'atoms', []) or []:
                    chain = getattr(atom, 'chain', None)
                    resid = getattr(atom, 'resid', None)
                    if chain and resid is not None:
                        pairs.add((chain, resid))
                if not pairs:
                    continue
                clauses = [f"(:{c} and {r})" for c, r in sorted(pairs)]
                label = f"oniom_pick_site_{i}"
                _viewer.highlight(
                    " or ".join(clauses),
                    style="ball+stick",
                    color=f"palette:{i}",
                    label=label,
                )
                applied.append(label)
            self._oniom_pick_labels = applied
        except Exception:
            self._oniom_pick_labels = []

        # Prompt user to select site(s)
        while True:
            try:
                selection_str = prompt_with_context(
                    self.processor,
                    "\nSelect RedoxSite(s) for ONIOM calculation (e.g. 1, 1-3, 1,3,5, or all)",
                    default="1",
                    module="ONIOM QM/MM Preparator",
                    description="Select RedoxSite(s) for ONIOM calculation",
                )

                # Parse selection string into list of 1-based indices
                indices = self._parse_site_selection(selection_str, len(redox_sites))

                if indices is None:
                    self.processor.console.print(
                        f"[red]Invalid selection. Please choose from 1-{len(redox_sites)}[/red]"
                    )
                    continue

                # Store 0-based indices
                self.selected_site_indices = [i - 1 for i in indices]
                self.selected_site_index = self.selected_site_indices[0]

                if len(indices) == 1:
                    # Single site selection
                    self.selected_redox_site = redox_sites[self.selected_site_index]
                    self.processor.console.print(
                        f"\n[green]✓ Selected RedoxSite {indices[0]} with "
                        f"{len(self.selected_redox_site.atoms)} atoms[/green]"
                    )
                else:
                    # Multiple sites - merge into a combined RedoxSite
                    self.selected_redox_site = self._merge_redox_sites(
                        [redox_sites[i] for i in self.selected_site_indices]
                    )
                    site_labels = ", ".join(str(i) for i in indices)
                    self.processor.console.print(
                        f"\n[green]✓ Selected RedoxSites {site_labels} "
                        f"(merged: {len(self.selected_redox_site.atoms)} atoms)[/green]"
                    )
                # Drop the per-row pick labels — Hook B will re-show the
                # selected site as the HIGH layer below.
                self._clear_oniom_pick_labels()
                return True

            except KeyboardInterrupt:
                self._clear_oniom_pick_labels()
                self.processor.console.print("[yellow]Selection cancelled[/yellow]")
                return False

    def _clear_oniom_pick_labels(self) -> None:
        """Remove the per-row redox-site pick highlights (Hook A)."""
        labels = getattr(self, "_oniom_pick_labels", None) or []
        if not labels:
            return
        try:
            from proprep.structure_prep.viewer_coordinator import (
                viewer as _viewer,
            )
            for lbl in labels:
                _viewer.unhighlight(lbl)
        except Exception:
            pass
        self._oniom_pick_labels = []

    @staticmethod
    def _parse_site_selection(selection_str: str, n_sites: int) -> Optional[List[int]]:
        """
        Parse a site selection string into a sorted list of unique 1-based indices.

        Supports: single int, comma-separated, ranges (e.g. 1-3), or 'all'.
        Returns None if any index is out of range or the string is unparseable.
        """
        selection_str = selection_str.strip().lower()
        if selection_str == "all":
            return list(range(1, n_sites + 1))

        indices = set()
        for part in selection_str.split(","):
            part = part.strip()
            if "-" in part:
                bounds = part.split("-", 1)
                try:
                    start, end = int(bounds[0].strip()), int(bounds[1].strip())
                except ValueError:
                    return None
                if start > end:
                    return None
                indices.update(range(start, end + 1))
            else:
                try:
                    indices.add(int(part))
                except ValueError:
                    return None

        if not indices or any(i < 1 or i > n_sites for i in indices):
            return None

        return sorted(indices)

    @staticmethod
    def _merge_redox_sites(sites: List[RedoxSite]) -> RedoxSite:
        """Merge multiple RedoxSite objects into a single combined site."""
        merged = RedoxSite(
            site_id="_".join(s.site_id for s in sites),
            structure_id=sites[0].structure_id
        )

        seen_coords = set()
        for site in sites:
            for center in site.centers:
                if center.coords not in seen_coords:
                    merged.add_center(center)
                    seen_coords.add(center.coords)
            for atom in site.atoms:
                if atom.coords not in seen_coords:
                    merged.add_atom(atom)
                    seen_coords.add(atom.coords)
            for bond in site.bonds:
                merged.bonds.append(bond)

        merged.detection_method = "merged"
        merged.site_type = sites[0].site_type
        return merged

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Layer Configuration

    def configure_layers(self) -> bool:
        """
        Configure ONIOM layer assignments interactively.

        Returns:
            True if configuration successful, False otherwise
        """
        if self.selected_redox_site is None:
            self.processor.console.print("\n[red]No RedoxSite selected.[/red]")
            self.processor.console.print("Please select a RedoxSite first.")
            return False

        self.processor.console.print("\n[bold cyan]ONIOM Layer Configuration[/bold cyan]\n")

        # Educational panel
        self.processor.console.print(Panel(
            "[bold]Understanding ONIOM Layers:[/bold]\n\n"
            "• [green]HIGH layer[/green]: QM region - treated with quantum mechanics\n"
            "  (e.g., B3LYP/6-31G*). Include the active site and chemically important residues.\n\n"
            "• [yellow]MEDIUM layer[/yellow]: Intermediate - treated with semi-empirical or lower-level QM\n"
            "  (optional, 3-layer ONIOM only). Use for large QM regions.\n\n"
            "• [blue]LOW layer[/blue]: MM region - treated with molecular mechanics\n"
            "  (AMBER forcefield). Includes bulk protein and solvent.\n\n"
            "[grey50]ProPrep automatically adds alpha carbons from neighboring residues to prevent\n"
            "cutting peptide bonds, and merges nearby fragments for chemical sensibility.[/grey50]",
            title="Layer Selection Guide",
            border_style="cyan",
            expand=False,
        ))

        # Ask number of layers
        while True:
            self.n_layers = int_prompt_with_context(
                self.processor,
                "\nNumber of ONIOM layers (2 or 3)",
                default=2,
                module="ONIOM QM/MM Preparator",
                description="Number of ONIOM layers",
            )
            if self.n_layers in (2, 3):
                break
            self.processor.console.print("[red]Please enter 2 or 3[/red]")

        # Configure HIGH layer - auto-populate from RedoxSite
        self.processor.console.print("\n[bold green]HIGH Layer (QM Region)[/bold green]")

        # Map RedoxSite resid (1-based PDB numbering) to parmed residue index
        # (0-based). parmed internally uses 0-based indexing, so parm_idx = resid - 1.
        n_parm_residues = len(self._parm.residues)
        site_resids = set(atom.resid for atom in self.selected_redox_site.atoms)
        self.high_residue_indices = [
            rid - 1 for rid in site_resids
            if 0 < rid <= n_parm_residues
        ]

        self.processor.console.print(
            f"\n[green]✓ Automatically added {len(self.high_residue_indices)} residues "
            f"from RedoxSite to HIGH layer:[/green]"
        )
        for res_idx in sorted(self.high_residue_indices):
            res = self._parm.residues[res_idx]
            n_atoms = len(res.atoms)
            self.processor.console.print(
                f"    [green]{res.name}{res.number}[/green] [grey50]({n_atoms} atoms)[/grey50]"
            )

        # Offer to modify residue list
        self.processor.console.print(
            "\n[grey50]You can add or remove entire residues from the HIGH layer,[/grey50]"
        )
        self.processor.console.print(
            "[grey50]or trim specific residues to include only their sidechain[/grey50]"
        )
        self.processor.console.print(
            "[grey50](backbone stays in LOW, QM/MM cut at CA-CB bond).[/grey50]"
        )

        if confirm_with_context(
            self.processor, "\nAdd or remove whole residues?", default=False,
            module="ONIOM QM/MM Preparator",
            description="Add or remove whole residues from HIGH layer",
        ):
            self._configure_residue_index_list("HIGH", self.high_residue_indices)

        if confirm_with_context(
            self.processor,
            "Trim any residues to sidechain only (CA-CB cut)?",
            default=False,
            module="ONIOM QM/MM Preparator",
            description="Trim residues to sidechain only for CA-CB boundary cut",
        ):
            self._trim_residues_to_sidechain()

        # Configure MEDIUM layer if 3-layer
        if self.n_layers == 3:
            self.processor.console.print("\n[bold yellow]MEDIUM Layer (Intermediate Region)[/bold yellow]")
            self.processor.console.print("Enter residues for MEDIUM layer treatment.")
            self._configure_residue_index_list("MEDIUM", self.medium_residue_indices)

        # Configure freezing
        freeze_str = prompt_with_context(
            self.processor,
            "\nFreeze atoms beyond distance from QM region (Å, 0 for none)",
            default="0",
            module="ONIOM QM/MM Preparator",
            description="Freeze distance cutoff in Angstroms (0 to disable)",
        )
        try:
            freeze_val = float(freeze_str)
        except ValueError:
            freeze_val = 0.0

        if freeze_val > 0:
            self.freeze_distance_cutoff = freeze_val
            self.processor.console.print(
                f"[green]✓ Will freeze atoms beyond {self.freeze_distance_cutoff}Å[/green]"
            )

            # Hook C: visualise the freeze-cutoff result before the
            # user commits — show active LOW (within cutoff of HIGH
            # geometric centre) in element colours so the frozen
            # residues are visible by absence (they fall back to the
            # default ribbon rep). Computed from high_residue_indices
            # because layer_membership isn't populated until MM
            # atom typing.
            try:
                from proprep.structure_prep.viewer_coordinator import (
                    viewer as _viewer,
                )
                from .structure_utils import (
                    calculate_geometric_center, calculate_distance,
                )

                high_coords = []
                high_resnum_set = set()
                for res_idx in self.high_residue_indices:
                    res = self._parm.residues[res_idx]
                    high_resnum_set.add(res.idx + 1)
                    for atom in res.atoms:
                        high_coords.append((atom.xx, atom.xy, atom.xz))

                if high_coords:
                    centre = calculate_geometric_center(high_coords)
                    active_low_resnums = set()
                    for res in self._parm.residues:
                        if res.idx in self.high_residue_indices:
                            continue
                        for atom in res.atoms:
                            if calculate_distance(
                                (atom.xx, atom.xy, atom.xz), centre,
                            ) <= self.freeze_distance_cutoff:
                                active_low_resnums.add(res.idx + 1)
                                break

                    _viewer.unhighlight("oniom_active_low")
                    if active_low_resnums:
                        sel = "(" + " or ".join(
                            str(r) for r in sorted(active_low_resnums)
                        ) + ")"
                        _viewer.highlight(
                            sel, style="ball+stick", color="element",
                            label="oniom_active_low",
                        )
            except Exception:
                pass
        else:
            self.freeze_distance_cutoff = None

        self.processor.console.print("\n[green]✓ Layer configuration complete[/green]")

        # Hook B: preview the HIGH layer assignment in magenta before
        # the user moves on to MM configuration. This is the
        # standalone "what's in the QM region" view; the full
        # multi-layer view (HIGH + boundary + active LOW) happens
        # after MM atom typing via launch_oniom_viewer.
        try:
            from proprep.structure_prep.viewer_coordinator import (
                viewer as _viewer,
            )
            pdb_file = self._upstream_pdb_paths[0] if self._upstream_pdb_paths else None
            if pdb_file:
                _viewer.show_structure(pdb_file)

            high_resnums = sorted(idx + 1 for idx in self.high_residue_indices)
            _viewer.unhighlight("oniom_high")
            if high_resnums:
                sel = "(" + " or ".join(str(r) for r in high_resnums) + ")"
                _viewer.highlight(
                    sel, style="ball+stick", color="#ff00ff",
                    label="oniom_high",
                )

            if self.n_layers == 3 and self.medium_residue_indices:
                med_resnums = sorted(
                    idx + 1 for idx in self.medium_residue_indices
                )
                _viewer.unhighlight("oniom_medium")
                if med_resnums:
                    med_sel = "(" + " or ".join(str(r) for r in med_resnums) + ")"
                    _viewer.highlight(
                        med_sel, style="ball+stick", color="#ffff00",
                        label="oniom_medium",
                    )
        except Exception:
            pass

        return True

    def _configure_residue_index_list(self, layer_name: str, index_list: List[int]):
        """
        Interactively add or remove residues from a layer by residue number.

        Args:
            layer_name: Name of the layer (for display)
            index_list: List of parmed residue indices to modify in place
        """
        # Build resnum → parm residue index lookup
        resnum_to_idx = {res.number: res.idx for res in self._parm.residues}

        self.processor.console.print(
            f"\nEnter residue numbers to add (e.g., '100')"
        )
        self.processor.console.print(
            "Prefix with [bold]-[/bold] to remove (e.g., '-100')"
        )
        self.processor.console.print("Press Enter with empty input when done.\n")

        while True:
            residue_input = prompt_with_context(
                self.processor,
                f"  {layer_name} residue number (or Enter to finish)", default="",
                module="ONIOM QM/MM Preparator",
                description=f"Add/remove residue in {layer_name} layer",
            )

            if not residue_input:
                break

            try:
                raw = residue_input.strip()
                removing = raw.startswith("-")
                if removing:
                    raw = raw[1:]

                resnum = int(raw.strip())
                res_idx = resnum_to_idx.get(resnum)

                if res_idx is None:
                    self.processor.console.print(
                        f"[red]Residue {resnum} not found in topology[/red]"
                    )
                    continue

                res = self._parm.residues[res_idx]

                if removing:
                    if res_idx in index_list:
                        index_list.remove(res_idx)
                        self.processor.console.print(
                            f"[yellow]✓ Removed {res.name}{resnum}[/yellow]"
                        )
                    else:
                        self.processor.console.print(
                            f"[yellow]{res.name}{resnum} not in {layer_name} list[/yellow]"
                        )
                    continue

                if res_idx in index_list:
                    self.processor.console.print(
                        f"[yellow]Residue {res.name}{resnum} already in list[/yellow]"
                    )
                    continue

                index_list.append(res_idx)
                self.processor.console.print(
                    f"[green]✓ Added {res.name}{resnum} ({len(res.atoms)} atoms)[/green]"
                )

            except ValueError:
                self.processor.console.print("[red]Invalid input. Enter a residue number.[/red]")

        self.processor.console.print(
            f"\n[green]{layer_name} layer: {len(index_list)} residues[/green]"
        )

    # Backbone atom names that stay in LOW during a sidechain-only trim
    BACKBONE_ATOMS = {'N', 'H', 'CA', 'HA', 'HA2', 'HA3', 'C', 'O'}

    def _trim_residues_to_sidechain(self):
        """
        Trim selected HIGH residues to sidechain-only (CA-CB cut).

        For each trimmed residue:
        - Remove the residue from high_residue_indices
        - Add only sidechain atoms (CB onward) to high_atom_indices
        - Backbone atoms (N, H, CA, HA, C, O) remain LOW
        - The CA-CB bond becomes the QM/MM boundary (handled by boundary detection)
        """
        console = self.processor.console

        if not self.high_residue_indices:
            console.print("[yellow]No HIGH residues to trim[/yellow]")
            return

        console.print("\n[bold]Select residues to trim to sidechain only:[/bold]\n")

        # Show eligible residues (protein residues with sidechains)
        eligible = []
        for res_idx in sorted(self.high_residue_indices):
            res = self._parm.residues[res_idx]
            atom_names = {a.name for a in res.atoms}
            has_sidechain = bool(atom_names - self.BACKBONE_ATOMS)
            has_backbone = 'CA' in atom_names and 'CB' in atom_names
            if has_sidechain and has_backbone:
                eligible.append(res_idx)
                console.print(
                    f"  [{len(eligible)}] {res.name}{res.idx + 1} "
                    f"({len(atom_names - self.BACKBONE_ATOMS)} sidechain atoms)"
                )

        if not eligible:
            console.print("[yellow]No eligible residues (need both backbone and sidechain)[/yellow]")
            return

        console.print(
            "\nEnter residue numbers to trim (e.g., '1' or '1,3'), "
            "or 'all', or Enter to skip."
        )

        choice = prompt_with_context(
            self.processor,
            "Residues to trim",
            default="",
            module="ONIOM QM/MM Preparator",
            description="Select residues for sidechain-only CA-CB cut",
        ).strip()

        if not choice:
            return

        # Parse selection
        if choice.lower() == 'all':
            selected = list(range(1, len(eligible) + 1))
        else:
            selected = []
            for part in choice.split(','):
                try:
                    selected.append(int(part.strip()))
                except ValueError:
                    pass

        trimmed_count = 0
        for sel_num in selected:
            if sel_num < 1 or sel_num > len(eligible):
                console.print(f"[yellow]Invalid selection: {sel_num}[/yellow]")
                continue

            res_idx = eligible[sel_num - 1]
            res = self._parm.residues[res_idx]

            # Move sidechain atoms to high_atom_indices
            sidechain_atoms = []
            for atom in res.atoms:
                if atom.name not in self.BACKBONE_ATOMS:
                    self.high_atom_indices.add(atom.idx)
                    sidechain_atoms.append(atom.name)

            # Remove from whole-residue HIGH list
            if res_idx in self.high_residue_indices:
                self.high_residue_indices.remove(res_idx)

            trimmed_count += 1
            console.print(
                f"  [green]✓ {res.name}{res.idx + 1}: sidechain only "
                f"({', '.join(sidechain_atoms)})[/green]"
            )

        if trimmed_count:
            console.print(
                f"\n[green]✓ Trimmed {trimmed_count} residue(s) to sidechain only. "
                f"CA-CB bonds will be the QM/MM boundary.[/green]"
            )

    def launch_oniom_viewer(self):
        """
        Show ONIOM layers in the coordinator-managed 3D viewer.

        Shows:
        - HIGH layer (active): magenta ball+stick
        - Boundary atoms: yellow ball+stick
        - LOW layer (active, near QM): element-colored ball+stick
        - LOW layer (frozen): default protein rep (cartoon)

        Routed through the coordinator with force=True so this view
        shares the rep manager with all other module hooks. The
        per-layer labels stay live across subsequent hooks (HIGH-layer
        hook in configure_layers, boundary hook in
        _detect_boundary_atoms, etc.) so calling this method gives
        the user a "see everything together" view at any point.
        """
        pdb_file = self._upstream_pdb_paths[0] if self._upstream_pdb_paths else None
        if not pdb_file:
            self.processor.console.print("[yellow]No PDB file available for viewing[/yellow]")
            return

        try:
            from proprep.structure_prep.viewer_coordinator import (
                viewer as _viewer,
            )
        except Exception:
            self.processor.console.print("[yellow]Viewer coordinator unavailable[/yellow]")
            return

        # HIGH layer residues (original + expanded), boundary, and
        # active LOW residues — same logic as the standalone version.
        high_resnums: set = set()
        boundary_resnums: set = set()
        active_low_resnums: set = set()

        if self.layer_membership and self._parm:
            boundary_set = set(self.boundary_atom_indices)

            for atom_idx, layer in self.layer_membership.items():
                atom = self._parm.atoms[atom_idx]
                resnum = atom.residue.idx + 1  # 1-based for NGL

                if layer == ONIOMLayer.HIGH:
                    high_resnums.add(resnum)
                elif atom_idx in boundary_set:
                    boundary_resnums.add(resnum)

            # Compute frozen vs active LOW residues using freeze cutoff
            if self.freeze_distance_cutoff:
                from .structure_utils import calculate_geometric_center, calculate_distance

                high_coords = [
                    (self._parm.atoms[idx].xx, self._parm.atoms[idx].xy, self._parm.atoms[idx].xz)
                    for idx, layer in self.layer_membership.items()
                    if layer == ONIOMLayer.HIGH
                ]
                if high_coords:
                    center = calculate_geometric_center(high_coords)
                    for atom_idx, layer in self.layer_membership.items():
                        if layer == ONIOMLayer.LOW:
                            atom = self._parm.atoms[atom_idx]
                            resnum = atom.residue.idx + 1
                            dist = calculate_distance(
                                (atom.xx, atom.xy, atom.xz), center)
                            if dist <= self.freeze_distance_cutoff:
                                active_low_resnums.add(resnum)
        elif self.high_residue_indices and self._parm:
            # Before MM config: just show selected residues
            for res_idx in self.high_residue_indices:
                high_resnums.add(res_idx + 1)

        # Helper: build an NGL residue-list selection like
        # ``(123 or 124 or 125)`` from a set of residue numbers.
        def _ngl_resnum_clause(resnums):
            if not resnums:
                return None
            return "(" + " or ".join(str(r) for r in sorted(resnums)) + ")"

        # Refresh: drop prior layer reps before redrawing so a
        # second view call replaces rather than accumulates.
        for label in (
            "oniom_high", "oniom_boundary",
            "oniom_active_low", "oniom_neighborhood",
        ):
            _viewer.unhighlight(label)

        _viewer.show_structure(pdb_file, force=True)

        if high_resnums:
            sel = _ngl_resnum_clause(high_resnums)
            if sel:
                _viewer.highlight(
                    sel, style="ball+stick", color="#ff00ff",
                    label="oniom_high",
                )

        if boundary_resnums:
            sel = _ngl_resnum_clause(boundary_resnums)
            if sel:
                _viewer.highlight(
                    sel, style="ball+stick", color="#ffff00",
                    label="oniom_boundary",
                )

        if active_low_resnums:
            # Subtract HIGH and boundary so each rep covers a distinct
            # residue set (otherwise the lower-priority element colour
            # would clobber the magenta/yellow visually).
            display_set = active_low_resnums - high_resnums - boundary_resnums
            sel = _ngl_resnum_clause(display_set)
            if sel:
                _viewer.highlight(
                    sel, style="ball+stick", color="element",
                    label="oniom_active_low",
                )

        self.processor.console.print(
            "[grey50]ONIOM layers shown in viewer (HIGH magenta, boundary yellow, "
            "active LOW element-coloured). Toggle reps in the rep manager.[/grey50]"
        )

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # QM Settings Configuration

    def configure_qm_settings(self) -> bool:
        """
        Configure QM calculation settings interactively.

        Now informed by MM atom typing - shows charge context from HIGH layer.

        Returns:
            True if configuration successful
        """
        self.processor.console.print("\n[bold cyan]QM Calculation Settings[/bold cyan]\n")

        # Show context from MM configuration if available
        if self.atom_typing_complete:
            self.processor.console.print(Panel(
                f"[bold]HIGH Layer Summary[/bold]\n"
                f"  Atoms:         {sum(1 for layer in self.layer_membership.values() if layer == ONIOMLayer.HIGH)}\n"
                f"  Total Charge:  {self.high_layer_charge:.4f}\n"
                f"  Suggested QM Charge: {self.suggested_qm_charge}",
                title="Context from MM Configuration",
                border_style="grey50",
                expand=False,
            ))
            self.processor.console.print()

        # Functional
        self.qm_functional = prompt_with_context(
            self.processor,
            "QM functional (e.g., B3LYP, M06-2X, wB97X-D)",
            default=self.qm_functional,
            module="ONIOM QM/MM Preparator",
            description="QM functional",
        )

        # Basis set
        self.qm_basis_set = prompt_with_context(
            self.processor,
            "QM basis set (e.g., 6-31G*, 6-311G(d,p), def2-SVP)",
            default=self.qm_basis_set,
            module="ONIOM QM/MM Preparator",
            description="QM basis set",
        )

        # Charge - use suggested charge if available
        charge_prompt = "QM region charge"
        if self.atom_typing_complete:
            charge_prompt = f"QM region charge (HIGH layer total: {self.high_layer_charge:.2f})"
            charge_default = self.suggested_qm_charge
        else:
            charge_default = self.qm_charge

        self.qm_charge = int_prompt_with_context(
            self.processor,
            charge_prompt,
            default=charge_default,
            module="ONIOM QM/MM Preparator",
            description="QM region charge",
        )

        # Multiplicity - show educational prompt
        self.processor.console.print("\n[grey50]Multiplicity Guide:[/grey50]")
        self.processor.console.print("[grey50]  1 = singlet (paired electrons)[/grey50]")
        self.processor.console.print("[grey50]  2 = doublet (1 unpaired electron)[/grey50]")
        self.processor.console.print("[grey50]  3 = triplet (2 unpaired electrons)[/grey50]")

        self.qm_multiplicity = int_prompt_with_context(
            self.processor,
            "QM region multiplicity",
            default=self.qm_multiplicity,
            module="ONIOM QM/MM Preparator",
            description="QM region multiplicity",
        )

        # Open-shell / broken-symmetry advisory.
        # ProPrep doesn't autodetect open-shell HIGH regions, so flag the
        # caveat unconditionally and let the user disregard it.
        self.processor.console.print(
            Panel.fit(
                "If the HIGH layer contains transition metals or other open-\n"
                "shell centers (Fe, Mn, Co, Ni, Cu, Cr, ...), a default SCF\n"
                "with the chosen multiplicity may converge to the wrong\n"
                "electronic state — broken-symmetry / antiferromagnetically\n"
                "coupled clusters (e.g. [4Fe-4S]) in particular.\n\n"
                "Gaussian's recommended antiferromagnetic-coupling (AFC)\n"
                "procedure is a multi-step job sequence:\n"
                "  1. High-spin Opt of the cluster fragment\n"
                "  2. Stability check (stable=opt) + re-optimize\n"
                "  3. Fragment guess for the BS state:\n"
                "     guess=(fragment=N,only)  with per-fragment\n"
                "     charges and spins\n"
                "  4. Final stable=opt on the BS state\n\n"
                "Full procedure: https://gaussian.com/afc/\n\n"
                "ProPrep does NOT add any of this scaffolding. For an\n"
                "open-shell HIGH region you'll typically converge the QM\n"
                "fragment to a BS reference separately, then read the\n"
                "checkpoint into the ONIOM input (Guess=Read).",
                title="[yellow]Open-shell / BS-DFT advisory[/yellow]",
                border_style="yellow",
            )
        )

        self.processor.console.print(
            f"\n[green]✓ QM: {self.qm_functional}/{self.qm_basis_set}, "
            f"charge={self.qm_charge}, mult={self.qm_multiplicity}[/green]"
        )

        return True

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # MM Settings Configuration

    def configure_mm_settings(self) -> bool:
        """
        Configure MM force field settings by reading atom types, charges, and
        parameters directly from the AMBER prmtop file.

        When driven by the unified QM/MM Preparator, the prmtop is available
        from upstream frame extraction. All atom types, partial charges, and
        bonded/nonbonded parameters are already assigned by tLEaP and stored
        in the topology, so no manual force field selection or atom typing
        from AMBERHOME files is required.

        Steps:
        1. Load prmtop with parmed → extract atom types and charges
        2. Analyze layer membership using parmed residue indices
        3. Analyze layer membership (HIGH/MEDIUM/LOW for each atom)
        4. Detect boundary atoms (HIGH layer atoms bonded to LOW layer atoms)
        5. Display charge breakdown tables
        6. Boundary charge handling
        7. Suggest the model-system formal charge (per-fragment rounding)

        Returns:
            True if configuration successful
        """
        # Validate prerequisites
        if self.selected_redox_site is None:
            self.processor.console.print("\n[red]No RedoxSite selected.[/red]")
            self.processor.console.print("Please select a RedoxSite first.")
            return False

        if not self.high_residue_indices:
            self.processor.console.print("\n[red]No HIGH layer residues configured.[/red]")
            self.processor.console.print("Please configure layers first.")
            return False

        self.processor.console.print("\n[bold cyan]═══ MM Force Field Configuration & Atom Typing ═══[/bold cyan]\n")

        # Step 1: Extract atom types and charges from prmtop
        if not self._extract_types_from_prmtop():
            return False

        # Step 3: Analyze layer membership
        self._analyze_layer_membership()

        # Step 4: Detect boundary atoms
        self._detect_boundary_atoms()

        # Step 5: Display charge breakdown tables
        self._display_charge_breakdown()

        # DEBUG: Write per-atom charge assignments to file
        self._write_debug_charge_file()

        # Step 6: Boundary charge handling
        if self.boundary_atom_indices:
            self._handle_boundary_charges()

        # Step 7: Calculate suggested QM charge
        self._calculate_suggested_qm_charge()

        # Step 8: Display summary
        self._display_mm_summary()

        # Mark atom typing as complete
        self.atom_typing_complete = True

        self.processor.console.print("\n[green]✓ MM configuration and atom typing complete[/green]")
        return True

    def _extract_types_from_prmtop(self) -> bool:
        """Extract atom types and charges from the AMBER prmtop using parmed.

        Populates self.atom_types and self.charges keyed by parmed atom
        index (int, 0-based).
        """
        import parmed

        if self._parm is None:
            prmtop_path = getattr(self, '_upstream_prmtop', '')
            pdb_path = (self._upstream_pdb_paths[0]
                        if self._upstream_pdb_paths else '')
            if not prmtop_path or not pdb_path:
                self.processor.console.print(
                    "[red]No prmtop/coordinates available. Run via the QM/MM Preparator.[/red]")
                return False
            self._parm = parmed.load_file(prmtop_path, pdb_path)

        self.processor.console.print("[bold]Step 1: Reading atom types and charges from prmtop[/bold]\n")

        try:
            parm = self._parm
            self.atom_types.clear()
            self.charges.clear()

            for atom in parm.atoms:
                self.atom_types[atom.idx] = atom.type
                self.charges[atom.idx] = atom.charge

            n_atoms = len(self.atom_types)
            n_types = len(set(self.atom_types.values()))
            self.processor.console.print(
                f"[green]✓ Extracted {n_atoms} atoms, {n_types} unique atom types from prmtop[/green]"
            )

            from collections import Counter
            type_counts = Counter(self.atom_types.values())
            table = Table(title="Atom Type Summary (top 15)")
            table.add_column("Type", style="cyan")
            table.add_column("Count", style="green")
            for atype, count in type_counts.most_common(15):
                table.add_row(atype, str(count))
            self.processor.console.print(table)

            return True

        except Exception as e:
            self.processor.console.print(f"[red]Error reading prmtop: {e}[/red]")
            logger.error(f"Prmtop extraction failed: {e}", exc_info=True)
            return False

    def _get_full_structure_pdb(self) -> Optional[str]:
        """
        Step 2: Get PDB file covering the FULL structure.

        When driven by the unified QM/MM Preparator, uses the PDB from
        upstream frame extraction. Otherwise falls back to workspace priority.

        Returns:
            Path to PDB file or None if not found
        """
        self.processor.console.print("\n[bold]Step 2: Locating Full Structure PDB[/bold]\n")

        if self._upstream_pdb_paths:
            pdb_file = self._upstream_pdb_paths[0]
            self.processor.console.print(
                f"[green]Using PDB from frame extraction: "
                f"{os.path.basename(pdb_file)}[/green]"
            )
            return pdb_file

        self.processor.console.print("[red]No PDB available. Run via the QM/MM Preparator.[/red]")
        return None

    def _analyze_layer_membership(self):
        """
        Analyze layer membership for each atom using parmed residue indices
        and individual atom overrides (for sidechain-only CA-CB cuts).
        """
        self.processor.console.print("\n[bold]Analyzing Layer Membership[/bold]\n")

        self.layer_membership.clear()

        high_set = set(self.high_residue_indices)
        medium_set = set(self.medium_residue_indices) if self.n_layers == 3 else set()

        high_count = 0
        medium_count = 0
        low_count = 0

        for residue in self._parm.residues:
            if residue.idx in high_set:
                # Whole residue is HIGH
                res_layer = ONIOMLayer.HIGH
            elif residue.idx in medium_set:
                res_layer = ONIOMLayer.MEDIUM
            else:
                res_layer = ONIOMLayer.LOW

            for atom in residue.atoms:
                # Individual atom override (sidechain-only trim)
                if atom.idx in self.high_atom_indices:
                    layer = ONIOMLayer.HIGH
                else:
                    layer = res_layer

                self.layer_membership[atom.idx] = layer
                if layer == ONIOMLayer.HIGH:
                    high_count += 1
                elif layer == ONIOMLayer.MEDIUM:
                    medium_count += 1
                else:
                    low_count += 1

        self.processor.console.print(f"[green]✓ Layer Analysis Complete:[/green]")
        self.processor.console.print(f"    HIGH layer: {high_count} atoms")
        if self.n_layers == 3:
            self.processor.console.print(f"    MEDIUM layer: {medium_count} atoms")
        self.processor.console.print(f"    LOW layer: {low_count} atoms")

    def _detect_boundary_atoms(self):
        """
        Detect QM/MM boundary atoms and expand the HIGH region using parmed.

        Full-residue backbone boundaries use a SINGLE-bond cut per side,
        leaving the flanking Cα in the LOW layer (one link H per cut):

        - N-terminal flank (LOW residue Y *before* a HIGH residue): C-cut —
          promote Y's {C, O} to HIGH and sever Cα(Y)–C(Y). The link H caps
          the carbonyl C → formamide terminus H–C(=O)–N–…
        - C-terminal flank (LOW residue Y *after* a HIGH residue): N-cut —
          promote Y's {N, H} to HIGH and sever N(Y)–Cα(Y). The link H caps
          the amide N → primary-amide terminus …–C(=O)–N(H)–H.

        This is the standard scaled hydrogen link atom (Gaussian default
        scale 0.7239 → ~C–H distance; cap charge zeroed). The cap's MM
        parameters are emitted explicitly as the HX block in oniom_writer,
        so the formamide cap is fully defined. It replaces the earlier
        acetamide / N-methylamide scheme, which promoted Cα into HIGH and
        put two link H's on one Cα — the minimal single-bond cut is
        preferred now that HX supplies the cap parameters.

        Pro is auto-promoted whole (its ring N can't be cleanly capped).
        A lone LOW residue flanked by HIGH on both sides would put two link
        H's on its Cα, so it is promoted whole as well. The loop re-runs
        after any whole-residue promotion to expose new boundaries.

        Sidechain-only mode (Cα stays LOW, Cβ+ in HIGH) is unchanged — a
        single Cα–Cβ cut with the link H at Cβ.
        """
        self.processor.console.print(
            "\n[bold]Detecting QM/MM Boundary Atoms & Expanding HIGH Region[/bold]\n"
        )

        self.boundary_atom_indices.clear()

        # --- 1. Build peptide-bond map (C of res i bonded to N of res i+1)
        peptide_next: Dict[int, int] = {}
        peptide_prev: Dict[int, int] = {}
        for bond in self._parm.bonds:
            a1, a2 = bond.atom1, bond.atom2
            if a1.residue.idx == a2.residue.idx:
                continue
            if a1.name == 'C' and a2.name == 'N':
                peptide_next[a1.residue.idx] = a2.residue.idx
                peptide_prev[a2.residue.idx] = a1.residue.idx
            elif a2.name == 'C' and a1.name == 'N':
                peptide_next[a2.residue.idx] = a1.residue.idx
                peptide_prev[a1.residue.idx] = a2.residue.idx

        PRO_NAMES = {"PRO", "CPRO", "NPRO"}

        # --- 2. Plan single-bond caps + whole-residue promotions.
        #
        # caps[neighbor_idx] = set of directions ("prev"/"next") the LOW
        # residue is cut from. A lone LOW residue cut from BOTH sides (a gap
        # between two HIGH residues) would put two link H's on one Cα, so it
        # is promoted whole instead; Pro is always promoted whole. Each
        # promotion can expose new boundaries, so iterate to a fixed point.
        high_runtime = set(self.high_residue_indices)
        pro_promoted: List[int] = []
        gap_promoted: List[int] = []

        caps: Dict[int, set] = {}
        while True:
            changed = False
            caps = {}
            for res_idx in list(high_runtime):
                for direction, neighbor_map in (
                    ("prev", peptide_prev),
                    ("next", peptide_next),
                ):
                    neighbor_idx = neighbor_map.get(res_idx)
                    if neighbor_idx is None or neighbor_idx in high_runtime:
                        continue
                    # Pro can't be cleanly capped; promote whole residue.
                    if self._parm.residues[neighbor_idx].name in PRO_NAMES:
                        high_runtime.add(neighbor_idx)
                        if neighbor_idx not in pro_promoted:
                            pro_promoted.append(neighbor_idx)
                        changed = True
                        continue
                    caps.setdefault(neighbor_idx, set()).add(direction)

            # Lone LOW residue cut from both sides -> promote whole.
            for neighbor_idx, dirs in caps.items():
                if len(dirs) == 2 and neighbor_idx not in high_runtime:
                    high_runtime.add(neighbor_idx)
                    gap_promoted.append(neighbor_idx)
                    changed = True

            if not changed:
                break

        # --- 3. Apply layer changes.
        for res_idx in pro_promoted + gap_promoted:
            for atom in self._parm.residues[res_idx].atoms:
                self.layer_membership[atom.idx] = ONIOMLayer.HIGH

        # Promote the cap atoms: {C,O} for a C-cut (prev flank -> formamide),
        # {N,H} for an N-cut (next flank -> primary amide). The flanking Cα
        # stays LOW and becomes the single boundary atom for that cut.
        n_formamide = 0
        n_primary_amide = 0
        cap_atoms_promoted = 0
        for neighbor_idx, dirs in caps.items():
            names = {a.name: a for a in self._parm.residues[neighbor_idx].atoms}
            if next(iter(dirs)) == "prev":
                promote = ('C', 'O')
                n_formamide += 1
            else:
                promote = ('N', 'H')
                n_primary_amide += 1
            for nm in promote:
                a = names.get(nm)
                if a is not None:
                    self.layer_membership[a.idx] = ONIOMLayer.HIGH
                    cap_atoms_promoted += 1

        # --- 4. Identify LOW boundary atoms: one flanking Cα per cut, which
        #        carries the single link H bonded to the promoted carbonyl C
        #        (C-cut) or amide N (N-cut).
        boundary_set: set = set()
        for neighbor_idx in caps:
            ca = next((a for a in self._parm.residues[neighbor_idx].atoms
                       if a.name == 'CA'), None)
            if ca is not None and self.layer_membership.get(ca.idx) == ONIOMLayer.LOW:
                boundary_set.add(ca.idx)

        # Sidechain-only cuts: CA(LOW) bonded to CB(HIGH via high_atom_indices).
        sidechain_boundaries = 0
        if self.high_atom_indices:
            for bond in self._parm.bonds:
                a1, a2 = bond.atom1, bond.atom2
                for ca, other in ((a1, a2), (a2, a1)):
                    if (ca.name == 'CA'
                            and self.layer_membership.get(ca.idx) == ONIOMLayer.LOW
                            and other.idx in self.high_atom_indices):
                        if ca.idx not in boundary_set:
                            boundary_set.add(ca.idx)
                            sidechain_boundaries += 1

        self.boundary_atom_indices = sorted(boundary_set)

        # --- 5. Report
        if self.boundary_atom_indices:
            self.processor.console.print(
                f"[green]✓ Found {len(self.boundary_atom_indices)} boundary atom(s)[/green]"
            )
            n_promoted = cap_atoms_promoted + sum(
                len(self._parm.residues[r].atoms)
                for r in pro_promoted + gap_promoted
            )
            if n_promoted:
                self.processor.console.print(
                    f"[green]✓ Expanded HIGH by {n_promoted} atom(s) "
                    f"(formyl/amide caps + whole-residue promotions)[/green]",
                    highlight=False,
                )
            if n_formamide:
                self.processor.console.print(
                    f"  {n_formamide} formamide cap(s) (Cα–C cut, Cα stays LOW)"
                )
            if n_primary_amide:
                self.processor.console.print(
                    f"  {n_primary_amide} primary-amide cap(s) (N–Cα cut, Cα stays LOW)"
                )
            if sidechain_boundaries:
                self.processor.console.print(
                    f"  {sidechain_boundaries} sidechain cut(s) (Cα–Cβ, Cα stays LOW)"
                )
            if pro_promoted:
                pro_names = ', '.join(
                    f"{self._parm.residues[i].name}{self._parm.residues[i].number}"
                    for i in pro_promoted
                )
                self.processor.console.print(
                    f"[yellow]⚠ Proline auto-promoted whole "
                    f"(ring N can't be cleanly capped): {pro_names}[/yellow]"
                )
            if gap_promoted:
                gap_names = ', '.join(
                    f"{self._parm.residues[i].name}{self._parm.residues[i].number}"
                    for i in gap_promoted
                )
                self.processor.console.print(
                    f"[green]✓ Promoted {len(gap_promoted)} lone gap residue(s) "
                    f"whole (HIGH on both sides): {gap_names}[/green]",
                    highlight=False,
                )
        else:
            self.processor.console.print("[green]No boundary atoms detected[/green]")

        # Hook D: highlight the detected boundary atoms in yellow so
        # the user sees exactly where the QM/MM cuts will fall. Single
        # label that the post-MM launch_oniom_viewer call refreshes
        # later with the unified layer view.
        try:
            from proprep.structure_prep.viewer_coordinator import (
                viewer as _viewer,
            )
            _viewer.unhighlight("oniom_boundary")
            if self.boundary_atom_indices:
                # Group by residue + atom name so the selection is
                # readable: "(:A and 123 and .CA) or (:A and 124 and .CB) or ..."
                clauses = []
                seen_pairs = set()
                for atom_idx in self.boundary_atom_indices:
                    atom = self._parm.atoms[atom_idx]
                    resnum = atom.residue.idx + 1
                    key = (resnum, atom.name)
                    if key in seen_pairs:
                        continue
                    seen_pairs.add(key)
                    clauses.append(f"({resnum} and .{atom.name})")
                if clauses:
                    _viewer.highlight(
                        " or ".join(clauses),
                        style="ball+stick",
                        color="#ffff00",
                        label="oniom_boundary",
                    )
        except Exception:
            pass

    def _calculate_layer_charges(self):
        """Calculate total charges for each layer."""
        self.high_layer_charge = 0.0
        self.medium_layer_charge = 0.0
        self.low_layer_charge = 0.0

        for atom_idx, charge in self.charges.items():
            layer = self.layer_membership.get(atom_idx, ONIOMLayer.LOW)
            if layer == ONIOMLayer.HIGH:
                self.high_layer_charge += charge
            elif layer == ONIOMLayer.MEDIUM:
                self.medium_layer_charge += charge
            else:
                self.low_layer_charge += charge

    def _display_charge_breakdown(self):
        """
        Step 6: Display charge breakdown tables for each layer.
        Shows residue-level summaries for each layer, plus detailed boundary atom table.
        """
        self.processor.console.print("\n[bold]Step 6: Charge Breakdown Analysis[/bold]\n")

        # Calculate layer charges first
        self._calculate_layer_charges()

        # Display HIGH layer table
        self._display_layer_charge_table(ONIOMLayer.HIGH, "HIGH Layer (QM Region)")

        # Display MEDIUM layer table if 3-layer ONIOM
        if self.n_layers == 3:
            self._display_layer_charge_table(ONIOMLayer.MEDIUM, "MEDIUM Layer")

        # Display LOW layer summary (per-chain, not per-residue — too many residues)
        self._display_low_layer_summary()

        # Display boundary atom details if any exist
        if self.boundary_atom_indices:
            self._display_boundary_atom_table()

    def _display_layer_charge_table(self, layer: ONIOMLayer, title: str):
        """Display charge breakdown table for a specific layer at residue level."""
        from rich.table import Table
        from collections import defaultdict

        residue_charges = defaultdict(float)
        residue_atoms = defaultdict(int)

        for res in self._parm.residues:
            for atom in res.atoms:
                if self.layer_membership.get(atom.idx) == layer:
                    key = (res.chain, res.number, res.name)
                    residue_charges[key] += self.charges.get(atom.idx, 0.0)
                    residue_atoms[key] += 1

        if not residue_charges:
            self.processor.console.print(f"\n[grey50]{title}: No residues[/grey50]")
            return

        sorted_residues = sorted(residue_charges.keys(), key=lambda x: (x[0], x[1]))

        table = Table(title=title, show_header=True, header_style="bold cyan")
        table.add_column("Chain", style="grey50", width=6)
        table.add_column("ResID", style="magenta", width=7)
        table.add_column("ResName", style="yellow", width=8)
        table.add_column("Atoms", style="cyan", justify="right", width=6)
        table.add_column("Net Charge", style="green", justify="right", width=11)
        table.add_column("Running Sum", style="green", justify="right", width=13)

        running_sum = 0.0
        max_display = 100

        for i, residue_key in enumerate(sorted_residues):
            chain_id, resid, resname = residue_key
            net_charge = residue_charges[residue_key]
            running_sum += net_charge

            if i < max_display:
                table.add_row(
                    chain_id, str(resid), resname,
                    str(residue_atoms[residue_key]),
                    f"{net_charge:9.4f}", f"{running_sum:11.4f}",
                )

        if len(sorted_residues) > max_display:
            table.add_row("[grey50]...[/grey50]",
                          f"[grey50]+{len(sorted_residues) - max_display} more[/grey50]",
                          "", "", "", "")

        self.processor.console.print(table)

        if layer == ONIOMLayer.HIGH:
            self.processor.console.print(f"\n[bold green]Total HIGH Layer Charge: {self.high_layer_charge:.4f}[/bold green]")
        elif layer == ONIOMLayer.MEDIUM:
            self.processor.console.print(f"\n[bold yellow]Total MEDIUM Layer Charge: {self.medium_layer_charge:.4f}[/bold yellow]")

    def _display_low_layer_summary(self):
        """Display a per-residue-name summary of the LOW layer charges."""
        from rich.table import Table
        from collections import defaultdict

        resname_data = defaultdict(lambda: {'count': 0, 'atoms': 0, 'charge': 0.0})

        for res in self._parm.residues:
            res_has_low = False
            for atom in res.atoms:
                if self.layer_membership.get(atom.idx) == ONIOMLayer.LOW:
                    resname_data[res.name]['atoms'] += 1
                    resname_data[res.name]['charge'] += self.charges.get(atom.idx, 0.0)
                    res_has_low = True
            if res_has_low:
                resname_data[res.name]['count'] += 1

        if not resname_data:
            self.processor.console.print("\n[grey50]LOW Layer (MM Region): No residues[/grey50]")
            return

        table = Table(title="LOW Layer (MM Region) — Summary", show_header=True, header_style="bold cyan")
        table.add_column("ResName", style="yellow", width=10)
        table.add_column("Residues", style="magenta", justify="right", width=10)
        table.add_column("Atoms", style="cyan", justify="right", width=8)
        table.add_column("Net Charge", style="green", justify="right", width=12)

        total_residues = 0
        total_atoms = 0
        total_charge = 0.0

        for resname in sorted(resname_data.keys()):
            data = resname_data[resname]
            total_residues += data['count']
            total_atoms += data['atoms']
            total_charge += data['charge']
            table.add_row(resname, str(data['count']), str(data['atoms']), f"{data['charge']:10.4f}")

        table.add_section()
        table.add_row("[bold]Total[/bold]", f"[bold]{total_residues}[/bold]",
                       f"[bold]{total_atoms}[/bold]", f"[bold]{total_charge:10.4f}[/bold]")

        self.processor.console.print(table)
        self.processor.console.print(f"\n[bold blue]Total LOW Layer Charge: {total_charge:.4f}[/bold blue]")

    def _display_boundary_atom_table(self):
        """Display detailed atom-level breakdown of boundary atoms."""
        from rich.table import Table

        self.processor.console.print("\n")

        high_boundary = []
        low_boundary = []

        for bc_idx in self.boundary_atom_indices:
            atom = self._parm.atoms[bc_idx]
            res = atom.residue
            info = {
                'idx': bc_idx,
                'name': atom.name,
                'resname': res.name,
                'resid': res.number,
                'charge': self.charges.get(bc_idx, 0.0),
                'atom_type': self.atom_types.get(bc_idx, "???"),
                'layer': self.layer_membership.get(bc_idx, ONIOMLayer.LOW),
            }
            if info['layer'] == ONIOMLayer.HIGH:
                high_boundary.append(info)
            else:
                low_boundary.append(info)

        table = Table(title="QM/MM Boundary Atoms", show_header=True, header_style="bold yellow")
        table.add_column("Layer", style="cyan", width=6)
        table.add_column("Index", style="grey50", width=7)
        table.add_column("ResID", style="magenta", width=7)
        table.add_column("ResName", style="yellow", width=8)
        table.add_column("Atom", style="cyan", width=6)
        table.add_column("Type", style="blue", width=6)
        table.add_column("Charge", style="green", justify="right", width=9)
        table.add_column("Role", style="yellow", width=20)

        for info in sorted(high_boundary, key=lambda x: x['idx']):
            table.add_row("HIGH", str(info['idx']), str(info['resid']),
                          info['resname'], info['name'], info['atom_type'],
                          f"{info['charge']:8.4f}", "Expanded QM region")

        for info in sorted(low_boundary, key=lambda x: x['idx']):
            table.add_row("LOW", str(info['idx']), str(info['resid']),
                          info['resname'], info['name'], info['atom_type'],
                          f"{info['charge']:8.4f}", "Link atom placement")

        self.processor.console.print(table)
        self.processor.console.print(
            f"\n[grey50]Found {len(high_boundary)} boundary atoms in HIGH layer "
            f"and {len(low_boundary)} boundary atoms in LOW layer[/grey50]"
        )
        self.processor.console.print(
            "[grey50]Boundary atoms at QM/MM interface may require charge adjustment to avoid overpolarization[/grey50]\n"
        )

    def _write_debug_charge_file(self):
        """Write per-atom charge assignments to file for debugging."""
        import os
        from collections import defaultdict

        workspace = self.get_workspace()
        working_dir = workspace.get('working_directory', '.')
        output_path = os.path.join(working_dir, "debug_atom_charges.txt")

        self.processor.console.print(f"\n[grey50]Writing debug charge file to: {output_path}[/grey50]")

        boundary_set = set(self.boundary_atom_indices)

        with open(output_path, 'w') as f:
            f.write("=" * 120 + "\n")
            f.write("PER-ATOM CHARGE ASSIGNMENTS - DEBUG OUTPUT\n")
            f.write("=" * 120 + "\n\n")
            f.write(f"Total atoms: {len(self._parm.atoms)}\n")
            f.write(f"Source: prmtop\n")
            f.write(f"HIGH layer charge: {self.high_layer_charge:.4f}\n")
            f.write(f"LOW layer charge: {self.low_layer_charge:.4f}\n")
            f.write(f"Total system charge: {self.high_layer_charge + self.low_layer_charge:.4f}\n\n")

            f.write(f"{'Index':<8} {'ResID':<6} {'ResName':<8} {'Atom':<6} {'Type':<6} {'Charge':>10} {'Layer':<6} {'Boundary':<10}\n")
            f.write("-" * 120 + "\n")

            for atom in self._parm.atoms:
                layer = self.layer_membership.get(atom.idx, ONIOMLayer.LOW)
                f.write(
                    f"{atom.idx:<8} "
                    f"{atom.residue.number:<6} "
                    f"{atom.residue.name:<8} "
                    f"{atom.name:<6} "
                    f"{atom.type:<6} "
                    f"{atom.charge:>10.6f} "
                    f"{layer.name:<6} "
                    f"{'YES' if atom.idx in boundary_set else 'NO':<10}\n"
                )

            # Residue-level summary
            f.write("\n" + "=" * 120 + "\n")
            f.write("RESIDUE-LEVEL CHARGE SUMMARY\n")
            f.write("=" * 120 + "\n\n")

            f.write(f"{'ResID':<6} {'ResName':<8} {'Atoms':<6} {'Net Charge':>12} {'Layer':<6}\n")
            f.write("-" * 120 + "\n")

            for res in self._parm.residues:
                res_charge = sum(a.charge for a in res.atoms)
                res_layer = self.layer_membership.get(res.atoms[0].idx, ONIOMLayer.LOW) if res.atoms else ONIOMLayer.LOW
                f.write(
                    f"{res.number:<6} "
                    f"{res.name:<8} "
                    f"{len(res.atoms):<6} "
                    f"{res_charge:>12.6f} "
                    f"{res_layer.name:<6}\n"
                )

        self.processor.console.print(f"[green]✓ Debug charge file written: {output_path}[/green]\n")

    def _handle_boundary_charges(self):
        """
        Handle QM/MM boundary charge treatment.

        Offers four approaches:
        1. No adjustment (mechanical embedding or user handles separately)
        2. Gaussian ScaleCharge keyword (recommended for electronic embedding)
        3. Charge-shift redistribution (zero link hosts, redistribute to MM neighbors)
        4. Manual per-atom override (advanced)
        """
        console = self.processor.console

        console.print("\n[bold cyan]QM/MM Boundary Charge Treatment[/bold cyan]\n")

        console.print(Panel(
            "In electronic embedding (EE), MM point charges polarize the QM wavefunction.\n"
            "Charges on atoms near the QM/MM boundary can cause [bold]overpolarization[/bold]\n"
            "because they are unphysically close to the QM density.\n\n"
            "This step lets you choose how to handle these boundary charges.\n"
            "For mechanical embedding, no adjustment is needed.",
            title="Why Boundary Charges Matter",
            border_style="grey50",
            expand=False,
        ))

        # Show current boundary atoms
        table = Table(title="Current Boundary Atoms (Link Atom Hosts)", show_header=True)
        table.add_column("Atom", style="yellow")
        table.add_column("Type", style="cyan")
        table.add_column("Charge", style="magenta", justify="right")
        for bc_idx in self.boundary_atom_indices:
            atom = self._parm.atoms[bc_idx]
            res = atom.residue
            charge = self.charges.get(bc_idx, 0.0)
            table.add_row(
                f"{res.name}{res.number} {atom.name}",
                self.atom_types.get(bc_idx, "?"),
                f"{charge:+.4f}",
            )
        console.print(table)
        console.print()

        console.print("[bold]Boundary charge treatment options:[/bold]\n")
        console.print(
            "  [cyan]1.[/cyan] [bold]No adjustment[/bold]\n"
            "     Keep all MM charges as-is. Use this for mechanical embedding,\n"
            "     or if you will handle charge scaling outside of proprep.\n"
        )
        console.print(
            "  [cyan]2.[/cyan] [bold]Gaussian ScaleCharge keyword[/bold] (recommended for electronic embedding)\n"
            "     Adds a ScaleCharge keyword to the Gaussian route line. Gaussian will\n"
            "     scale MM charges near the QM region during the QM calculation, without\n"
            "     modifying the charges in the input file.\n"
        )
        console.print(
            "  [cyan]3.[/cyan] [bold]Charge-shift redistribution[/bold]\n"
            "     Zero the link host atom (CA) charges and redistribute them to the\n"
            "     nearest bonded MM atoms. Preserves total system charge.\n"
        )
        console.print(
            "  [cyan]4.[/cyan] [bold]Manual per-atom override[/bold]\n"
            "     Edit each boundary atom charge individually. For advanced users.\n"
        )

        while True:
            choice = prompt_with_context(
                self.processor,
                "Select boundary charge treatment",
                default="1",
                module="ONIOM QM/MM Preparator",
                description="Select boundary charge treatment",
                options_map={"1": "Keep original charges", "2": "ScaleCharge keyword",
                             "3": "Charge-shift redistribution", "4": "Manual per-atom override"},
            )
            if choice in ("1", "2", "3", "4"):
                break
            console.print("[red]Please enter 1, 2, 3, or 4.[/red]")

        if choice == "1":
            console.print("[green]✓ Keeping original boundary charges (no adjustment)[/green]")
        elif choice == "2":
            self._apply_scale_charge_keyword()
        elif choice == "3":
            self._apply_charge_shift()
            self._calculate_layer_charges()
        elif choice == "4":
            self._manual_boundary_charge_override()
            self._calculate_layer_charges()

    def _apply_scale_charge_keyword(self):
        """
        Add Gaussian ScaleCharge keyword to the route line.

        ScaleCharge=ijklmn: each digit × 0.2 gives the scale factor for MM charges
        at increasing bond distances from the QM region.
        """
        console = self.processor.console

        console.print(Panel(
            "[bold]ScaleCharge syntax:[/bold] ScaleCharge=N  where N is one or more digits.\n\n"
            "Each digit is multiplied by 0.2 to get the scale factor for MM charges at\n"
            "that bond distance from the QM region:\n\n"
            "  [cyan]Digit position (right to left):[/cyan]\n"
            "    rightmost digit  → atoms [bold]1 bond[/bold] from QM  (link host atoms)\n"
            "    next digit       → atoms [bold]2 bonds[/bold] from QM\n"
            "    next digit       → atoms [bold]3 bonds[/bold] from QM\n"
            "    ...and so on\n\n"
            "  [cyan]Digit values:[/cyan]\n"
            "    0 → scale factor 0.0 (charge completely zeroed)\n"
            "    1 → scale factor 0.2\n"
            "    3 → scale factor 0.6\n"
            "    5 → scale factor 1.0 (charge unscaled)\n\n"
            "  [cyan]Examples:[/cyan]\n"
            "    [green]500[/green]    → zero charges within 2 bonds, full charges beyond [bold](Gaussian default)[/bold]\n"
            "    [green]5000[/green]   → zero charges within 3 bonds, full charges beyond\n"
            "    [green]53100[/green]  → graduated: 0.0 at 1 bond, 0.0 at 2, 0.2 at 3, 0.6 at 4, full beyond\n\n"
            "[grey50]Note: digits must be monotonically decreasing from left to right.\n"
            "ScaleCharge implies electronic embedding (EmbedCharge).[/grey50]",
            title="Gaussian ScaleCharge Syntax",
            border_style="cyan",
            expand=False,
        ))

        value = prompt_with_context(
            self.processor,
            "Enter ScaleCharge value (digits only, e.g. 500)",
            default="500",
            module="ONIOM QM/MM Preparator",
            description="ScaleCharge value",
        )

        # Validate: must be digits, monotonically decreasing requirement
        # is enforced by Gaussian, but we can do a basic check
        value = value.strip()
        if not value.isdigit():
            console.print("[red]Invalid input — must be digits only. Using default 500.[/red]")
            value = "500"

        keyword = f"ScaleCharge={value}"

        # Add to additional keywords
        if self.additional_keywords:
            self.additional_keywords += f" {keyword}"
        else:
            self.additional_keywords = keyword

        # Describe what this will do
        digits = [int(d) for d in value]
        n_zeroed = sum(1 for d in reversed(digits) if d == 0)

        console.print(f"\n[green]✓ Will add '{keyword}' to Gaussian route line[/green]")
        if n_zeroed > 0:
            console.print(
                f"  Gaussian will zero MM charges within [cyan]{n_zeroed} bond(s)[/cyan] "
                f"of the QM region during the QM calculation."
            )
        console.print(
            "[grey50]  Note: charges in the input file are unchanged; "
            "Gaussian applies scaling at runtime.[/grey50]"
        )

    def _apply_charge_shift(self):
        """
        Zero link host atom charges and redistribute to nearest bonded MM atoms.

        Uses parmed bond topology to find actual bonded MM neighbors
        rather than hardcoded atom names.
        """
        console = self.processor.console
        console.print("\n[bold cyan]Charge-Shift Redistribution[/bold cyan]\n")

        for bc_idx in self.boundary_atom_indices:
            bc_atom = self._parm.atoms[bc_idx]
            res = bc_atom.residue
            ca_charge = self.charges.get(bc_idx, 0.0)
            label = f"{res.name}{res.number} {bc_atom.name}"

            # Find bonded atoms that remain in the LOW layer
            mm_neighbors = []
            for bond in bc_atom.bonds:
                partner = bond.atom1 if bond.atom2 is bc_atom else bond.atom2
                if self.layer_membership.get(partner.idx) == ONIOMLayer.LOW:
                    mm_neighbors.append(partner)

            if not mm_neighbors:
                console.print(
                    f"  [yellow]Warning: No LOW-layer neighbors found for {label}, "
                    f"skipping redistribution[/yellow]"
                )
                continue

            share = ca_charge / len(mm_neighbors)
            for partner in mm_neighbors:
                self.charges[partner.idx] = self.charges.get(partner.idx, 0.0) + share

            self.charges[bc_idx] = 0.0

            neighbor_str = ", ".join(
                f"{p.name} ({self.charges[p.idx]:+.4f})" for p in mm_neighbors
            )
            console.print(
                f"  {label}: [magenta]{ca_charge:+.4f}[/magenta] → [green]0.0000[/green]  "
                f"(redistributed to {neighbor_str})"
            )

        console.print(
            "\n[green]✓ Charge-shift redistribution complete. Total system charge preserved.[/green]"
        )

    def _manual_boundary_charge_override(self):
        """Allow user to manually edit each boundary atom charge."""
        console = self.processor.console
        console.print("\n[bold cyan]Manual Boundary Charge Override[/bold cyan]\n")
        console.print("[grey50]Press Enter to keep the current value for any atom.[/grey50]\n")

        for bc_idx in self.boundary_atom_indices:
            atom = self._parm.atoms[bc_idx]
            res = atom.residue
            current_charge = self.charges.get(bc_idx, 0.0)
            label = f"{res.name}{res.number} {atom.name}"

            new_charge_str = prompt_with_context(
                self.processor,
                f"  {label}  current={current_charge:+.4f}  new charge",
                default="",
                module="ONIOM QM/MM Preparator",
                description=f"Manual charge override for {label}",
            )

            if new_charge_str:
                try:
                    new_charge = float(new_charge_str)
                    self.charges[bc_idx] = new_charge
                    console.print(f"    [green]✓ {current_charge:+.4f} → {new_charge:+.4f}[/green]")
                except ValueError:
                    console.print("    [red]Invalid number, keeping original[/red]")

    def _calculate_suggested_qm_charge(self):
        """
        Suggest the formal charge of the ONIOM model system.

        The model system is what the user selected -- whole HIGH residues plus
        the side chains of CA-CB-trimmed residues -- terminated by capping
        groups (the promoted C=O or N-H of the flanking residue plus a link
        H). The caps are closed-shell neutral groups, so they add no formal
        charge and are not counted, even though their MM partial charges do
        not sum to zero (ff19SB: C+O = +0.03, N+H = -0.14).

        Formal charges are integers per fragment, so each fragment is rounded
        on its own and the integers are summed. Pooling the partial charges
        first and rounding once is wrong for trimmed side chains: an AMBER
        side chain alone carries a fraction (ASP -0.86, GLU -0.88, HIP +0.94)
        and four carboxylates pool to -3.43, which rounds to -3 instead of -4.
        """
        console = self.processor.console

        whole = [[a.idx for a in self._parm.residues[res_idx].atoms]
                 for res_idx in self.high_residue_indices]
        trimmed: Dict[int, List[int]] = {}
        for atom_idx in self.high_atom_indices:
            res_idx = self._parm.atoms[atom_idx].residue.idx
            trimmed.setdefault(res_idx, []).append(atom_idx)

        self.suggested_qm_charge = suggested_model_charge(
            self.charges, whole, list(trimmed.values()))

        console.print(
            f"\n[bold]Suggested QM region charge: {self.suggested_qm_charge}[/bold]"
        )
        console.print(
            "[grey50](Formal charge of the selected residues and trimmed side chains,\n"
            " each rounded to an integer. Capping groups -- the promoted C=O or N-H\n"
            " plus the link H -- are formally neutral and are not counted.)[/grey50]"
        )

        total_system_charge = round(
            self.high_layer_charge + self.medium_layer_charge + self.low_layer_charge
        )
        console.print(
            f"[bold]Total system (real) charge: {total_system_charge}[/bold]"
        )
        console.print(
            f"[grey50]Gaussian charge/multiplicity line: "
            f"{total_system_charge} 1 {self.suggested_qm_charge} 1[/grey50]"
        )

    def _display_mm_summary(self):
        """
        Display final MM configuration summary.
        """
        total_atoms = len(self.atom_types)
        high_atoms = sum(1 for layer in self.layer_membership.values() if layer == ONIOMLayer.HIGH)
        medium_atoms = sum(1 for layer in self.layer_membership.values() if layer == ONIOMLayer.MEDIUM)
        low_atoms = sum(1 for layer in self.layer_membership.values() if layer == ONIOMLayer.LOW)

        # Check charge balance
        total_charge = self.high_layer_charge + self.medium_layer_charge + self.low_layer_charge
        is_balanced = abs(total_charge - round(total_charge)) < 0.01

        self.processor.console.print("\n")
        self.processor.console.print(Panel(
            f"[bold]Force Field:[/bold]         {self.force_field_name}\n"
            f"[bold]Total Atoms Typed:[/bold]   {total_atoms}\n"
            f"\n"
            f"[bold]HIGH Layer:[/bold]          {high_atoms} atoms\n"
            f"  Total Charge:      {self.high_layer_charge:8.4f}\n"
            f"  Boundary Atoms:    {len(self.boundary_atom_indices)} atoms\n" +
            (f"\n[bold]MEDIUM Layer:[/bold]        {medium_atoms} atoms\n"
             f"  Total Charge:      {self.medium_layer_charge:8.4f}\n" if self.n_layers == 3 else "") +
            f"\n"
            f"[bold]LOW Layer:[/bold]           {low_atoms} atoms\n"
            f"  Total Charge:      {self.low_layer_charge:8.4f}\n"
            f"\n"
            f"[bold]Total System Charge:[/bold] {total_charge:8.4f}\n"
            f"{'[green]✓ Charge is balanced[/green]' if is_balanced else '[yellow]⚠ Charge may be unbalanced[/yellow]'}\n"
            f"\n"
            f"[bold]Suggested QM Charge:[/bold] {self.suggested_qm_charge}",
            title="MM Configuration Summary",
            border_style="cyan",
            expand=False,
        ))


    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Job Settings Configuration

    def configure_job_settings(self) -> bool:
        """
        Configure job execution settings interactively.

        Returns:
            True if configuration successful
        """
        self.processor.console.print("\n[bold cyan]Job Execution Settings[/bold cyan]\n")

        # Job type
        job_choices = ["Opt", "SP", "Freq", "Opt Freq"]
        self.job_type = prompt_with_context(
            self.processor,
            "Job type",
            default=self.job_type,
            choices=job_choices,
            module="ONIOM QM/MM Preparator",
            description="Gaussian job type",
        )

        # Processors
        self.n_processors = int_prompt_with_context(
            self.processor,
            "Number of processors",
            default=self.n_processors,
            module="ONIOM QM/MM Preparator",
            description="Number of processors",
        )

        # Memory
        self.memory_gb = int_prompt_with_context(
            self.processor,
            "Memory (GB)",
            default=self.memory_gb,
            module="ONIOM QM/MM Preparator",
            description="Memory in GB",
        )

        # Additional keywords — show route line preview, then ask for extra keywords to append.
        # SoftOnly is embedded in the AMBER spec inside ONIOM(...), matching the writer
        # (oniom_writer._write_route_section); a standalone Amber=SoftOnly would be ignored
        # or treated as a duplicate by Gaussian.
        oniom_methods = f"{self.qm_functional}/{self.qm_basis_set}:AMBER=SoftOnly"
        oniom_options = ["EmbedCharge"]
        preview_remaining = []
        if self.additional_keywords:
            oniom_option_prefixes = ("SCALECHARGE", "ONLYINPUTFILES", "INPUTFILES", "SVALUE")
            for kw in self.additional_keywords.split():
                if kw.upper().startswith(oniom_option_prefixes):
                    oniom_options.append(kw)
                else:
                    preview_remaining.append(kw)

        route_preview = f"#P ONIOM({oniom_methods})=({','.join(oniom_options)}) {self.job_type} Geom=Connect"
        if preview_remaining:
            route_preview += f" {' '.join(preview_remaining)}"

        self.processor.console.print(f"\n[bold]Route line preview:[/bold]")
        self.processor.console.print(f"  [cyan]{route_preview}[/cyan]", highlight=False)
        self.processor.console.print("[grey50]Enter any extra Gaussian keywords to append, or press Enter to skip.[/grey50]")

        kw_input = prompt_with_context(
            self.processor,
            "Extra keywords",
            default="",
            module="ONIOM QM/MM Preparator",
            description="Extra Gaussian keywords to append",
        )
        if kw_input.strip():
            if self.additional_keywords:
                self.additional_keywords += f" {kw_input.strip()}"
            else:
                self.additional_keywords = kw_input.strip()

        self.processor.console.print(
            f"\n[green]✓ Job: {self.job_type}, {self.n_processors} processors, {self.memory_gb}GB[/green]"
        )

        return True

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # ONIOM Preparation

    def prepare_oniom_setup(self) -> bool:
        """
        Prepare complete ONIOM setup using configured parameters.

        Returns:
            True if preparation successful, False otherwise
        """
        # Validate prerequisites
        if self.selected_redox_site is None:
            self.processor.console.print("\n[red]No RedoxSite selected.[/red]")
            self.processor.console.print("Please select a RedoxSite first.")
            return False

        if not self.high_residue_indices:
            self.processor.console.print("\n[red]No HIGH layer residues configured.[/red]")
            self.processor.console.print("Please configure layers first.")
            return False

        # Build layer_assignments from layer_membership using parmed
        layer_assignments_dict = None
        if self.atom_typing_complete and self.layer_membership:
            from proprep.oniom_prep.data_structures import LayerAssignment, FreezeFlag

            layer_assignments_dict = {}
            for atom_idx, layer in self.layer_membership.items():
                atom = self._parm.atoms[atom_idx]
                res = atom.residue
                layer_assignments_dict[atom_idx] = LayerAssignment(
                    atom_idx=atom_idx,
                    coords=(atom.xx, atom.xy, atom.xz),
                    layer=layer,
                    freeze=FreezeFlag.ACTIVE,
                    residue_idx=res.idx,
                    residue_name=res.name,
                    residue_number=res.number,
                    atom_name=atom.name,
                    element=atom.element_name if hasattr(atom, 'element_name') else atom.name[0],
                    assignment_reason="From MM configuration",
                )

        # Create ONIOMPreparator
        self.preparator = ONIOMPreparator(
            redox_site=self.selected_redox_site,
            logger_instance=logger,
            console=self.processor.console
        )

        # Prepare ONIOM setup
        try:
            oniom_setup = self.preparator.prepare_oniom_setup(
                high_residues=self.high_residue_indices,
                medium_residues=self.medium_residue_indices if self.n_layers == 3 else None,
                n_layers=self.n_layers,
                freeze_distance_cutoff=self.freeze_distance_cutoff,
                frozen_residues=self.frozen_residue_indices,
                qm_functional=self.qm_functional,
                qm_basis_set=self.qm_basis_set,
                qm_charge=self.qm_charge,
                qm_multiplicity=self.qm_multiplicity,
                mm_forcefield=self.mm_forcefield,
                force_field_name=self.force_field_name,
                prmtop_path=self._upstream_prmtop if self._using_upstream_frames else None,
                medium_method=self.medium_method if self.n_layers == 3 else None,
                medium_charge=self.medium_charge if self.n_layers == 3 else None,
                medium_multiplicity=self.medium_multiplicity if self.n_layers == 3 else None,
                job_type=self.job_type,
                n_processors=self.n_processors,
                memory_gb=self.memory_gb,
                additional_keywords=self.additional_keywords,
                atom_types=self.atom_types if self.atom_typing_complete else None,
                charges=self.charges if self.atom_typing_complete else None,
                layer_assignments=layer_assignments_dict,
                boundary_atom_indices=self.boundary_atom_indices if self.atom_typing_complete else None,
                parm=self._parm,
            )

            # Store in module
            self.update_workspace("oniom_setup", oniom_setup)

            return True

        except Exception as e:
            self.processor.console.print(f"\n[red]Error during ONIOM preparation: {e}[/red]")
            logger.error(f"ONIOM preparation failed: {e}", exc_info=True)
            return False

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Validation and Reporting

    def view_validation_report(self) -> bool:
        """
        Display detailed validation report.

        Returns:
            True if report displayed successfully
        """
        if self.preparator is None:
            self.processor.console.print("\n[yellow]ONIOM setup not prepared yet.[/yellow]")
            self.processor.console.print("Please run 'Prepare complete ONIOM setup' first.")
            return False

        self.preparator.print_validation_report()
        return True

    def write_oniom_input(self) -> bool:
        """
        Write Gaussian ONIOM input file.

        Returns:
            True if file written successfully
        """
        if self.preparator is None:
            self.processor.console.print("\n[yellow]ONIOM setup not prepared yet.[/yellow]")
            self.processor.console.print("Please run 'Prepare complete ONIOM setup' first.")
            return False

        # Get output directory from workspace
        workspace = self.get_workspace()
        base_dir = workspace.get('working_directory', os.getcwd())

        # Prompt for output filename
        site_label = "_".join(str(i + 1) for i in self.selected_site_indices) if self.selected_site_indices else str(self.selected_site_index + 1)
        default_filename = f"oniom_{site_label}.com"
        output_filename = prompt_with_context(
            self.processor,
            "\nOutput filename",
            default=default_filename,
            module="ONIOM QM/MM Preparator",
            description="ONIOM input output filename",
        )

        output_path = os.path.join(base_dir, "oniom_inputs", output_filename)

        # Write file
        success = self.preparator.write_input_file(
            output_path=output_path,
            include_comments=True,
            include_ff_parameters=True
        )

        if success:
            # Also write summary report
            summary_path = str(Path(output_path).with_suffix('.txt'))
            self.preparator.write_summary_report(summary_path)

            # Update workspace to track completion
            self.update_workspace("oniom_input_file", output_path)

        return success

    def write_diagnostic_input(self) -> bool:
        """
        Write diagnostic Gaussian ONIOM input file.

        Uses ONIOM=OnlyInputFiles so Gaussian writes sub-calculation
        input files and stops without running the calculation.

        Returns:
            True if file written successfully
        """
        if self.preparator is None:
            self.processor.console.print("\n[yellow]ONIOM setup not prepared yet.[/yellow]")
            return False

        workspace = self.get_workspace()
        base_dir = workspace.get('working_directory', os.getcwd())

        site_label = "_".join(str(i + 1) for i in self.selected_site_indices) if self.selected_site_indices else str(self.selected_site_index + 1)
        default_filename = f"oniom_{site_label}_diagnostic.com"
        output_filename = prompt_with_context(
            self.processor,
            "\nDiagnostic output filename",
            default=default_filename,
            module="ONIOM QM/MM Preparator",
            description="Diagnostic input output filename",
        )

        output_path = os.path.join(base_dir, "oniom_inputs", output_filename)

        success = self.preparator.write_diagnostic_file(
            output_path=output_path,
        )

        return success

    def write_summary_report(self) -> bool:
        """
        Write detailed summary report.

        Returns:
            True if report written successfully
        """
        if self.preparator is None:
            self.processor.console.print("\n[yellow]ONIOM setup not prepared yet.[/yellow]")
            return False

        # Get output directory
        workspace = self.get_workspace()
        base_dir = workspace.get('working_directory', os.getcwd())

        site_label = "_".join(str(i + 1) for i in self.selected_site_indices) if self.selected_site_indices else str(self.selected_site_index + 1)
        default_filename = f"oniom_summary_{site_label}.txt"
        output_filename = prompt_with_context(
            self.processor,
            "\nSummary report filename",
            default=default_filename,
            module="ONIOM QM/MM Preparator",
            description="Summary report filename",
        )

        output_path = os.path.join(base_dir, "oniom_inputs", output_filename)

        return self.preparator.write_summary_report(output_path)

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # View Configuration

    def view_current_setup(self) -> bool:
        """
        Display current ONIOM configuration.

        Returns:
            True always
        """
        self.processor.console.print("\n[bold cyan]Current ONIOM Configuration[/bold cyan]\n")

        # Selected site
        if self.selected_redox_site:
            if len(self.selected_site_indices) > 1:
                site_label = ", ".join(str(i + 1) for i in self.selected_site_indices)
                self.processor.console.print(
                    f"[green]Selected RedoxSites:[/green] Sites {site_label} "
                    f"(merged: {len(self.selected_redox_site.atoms)} atoms)"
                )
            else:
                self.processor.console.print(
                    f"[green]Selected RedoxSite:[/green] Site {self.selected_site_index + 1} "
                    f"({len(self.selected_redox_site.atoms)} atoms)"
                )
        else:
            self.processor.console.print("[yellow]Selected RedoxSite: None[/yellow]")

        # Layer configuration
        self.processor.console.print(f"\n[bold]Layer Configuration:[/bold]")
        self.processor.console.print(f"  Number of layers: {self.n_layers}")
        self.processor.console.print(f"  HIGH layer: {len(self.high_residue_indices)} residues")
        if self.n_layers == 3:
            self.processor.console.print(f"  MEDIUM layer: {len(self.medium_residue_indices)} residues")
        if self.freeze_distance_cutoff:
            self.processor.console.print(f"  Freeze cutoff: {self.freeze_distance_cutoff}Å")

        # QM settings
        self.processor.console.print(f"\n[bold]QM Settings:[/bold]")
        self.processor.console.print(f"  Method: {self.qm_functional}/{self.qm_basis_set}")
        self.processor.console.print(f"  Charge: {self.qm_charge}, Multiplicity: {self.qm_multiplicity}")

        # MM settings
        self.processor.console.print(f"\n[bold]MM Settings:[/bold]")
        self.processor.console.print(f"  Force field: {self.force_field_name}")

        # Atom typing status
        if self.atom_typing_complete:
            self.processor.console.print(f"\n[bold]Atom Typing:[/bold]")
            self.processor.console.print(f"  Total atoms typed: {len(self.atom_types)}")
            self.processor.console.print(f"  HIGH layer charge: {self.high_layer_charge:.4f}")
            if self.n_layers == 3:
                self.processor.console.print(f"  MEDIUM layer charge: {self.medium_layer_charge:.4f}")
            self.processor.console.print(f"  LOW layer charge: {self.low_layer_charge:.4f}")
            self.processor.console.print(f"  Boundary atoms: {len(self.boundary_atom_indices)}")
            self.processor.console.print(f"  [green]✓ Atom typing complete[/green]")
        else:
            self.processor.console.print(f"\n[bold]Atom Typing:[/bold] [yellow]Not yet performed[/yellow]")

        # Job settings
        self.processor.console.print(f"\n[bold]Job Settings:[/bold]")
        self.processor.console.print(f"  Type: {self.job_type}")
        self.processor.console.print(f"  Processors: {self.n_processors}, Memory: {self.memory_gb}GB")
        if self.additional_keywords:
            self.processor.console.print(f"  Additional keywords: {self.additional_keywords}")

        # Preparation status
        if self.preparator:
            self.processor.console.print("\n[green]✓ ONIOM setup prepared and ready[/green]")
        else:
            self.processor.console.print("\n[yellow]ONIOM setup not yet prepared[/yellow]")

        return True

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Reset Configuration

    def reset_configuration(self) -> bool:
        """
        Reset all configuration to defaults.

        Returns:
            True if reset successful
        """
        if confirm_with_context(
            self.processor,
            "\n[yellow]Reset all ONIOM configuration to defaults?[/yellow]", default=False,
            module="ONIOM QM/MM Preparator",
            description="Reset ONIOM configuration to defaults",
        ):
            self.__init__()
            self.processor.console.print("[green]✓ Configuration reset to defaults[/green]")
            return True

        return False

    def get_status_info(self, workspace) -> Dict[str, Any]:
        """Get status information for the workspace display."""
        redox_sites = workspace.get("detected_redox_sites", [])
        oniom_input_file = workspace.get("oniom_input_file")

        status = {
            "redox_sites_detected": len(redox_sites),
            "site_selected": self.selected_redox_site is not None,
            "layers_configured": bool(self.high_residue_indices),
            "setup_prepared": self.preparator is not None,
            "input_file_written": bool(oniom_input_file),
        }

        if self.selected_redox_site:
            status["selected_site_indices"] = [i + 1 for i in self.selected_site_indices]
            status["selected_site_index"] = self.selected_site_index + 1
            status["n_layers"] = self.n_layers

        if oniom_input_file:
            status["input_file"] = Path(oniom_input_file).name

        return status
