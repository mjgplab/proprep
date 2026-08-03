"""
ORCA QM/MM Preparator Module

Generates ORCA QM/MM input files from RedoxSite objects and AMBER topology.
Provides a guided, educational workflow for configuring QM/MM calculations
with proper QM region definition, active region selection, and force field conversion.
"""

import logging
import os
import math
import shutil
import subprocess
from typing import Any, Dict, List, Optional, Tuple
from pathlib import Path

import parmed

from rich.panel import Panel
from rich.table import Table
from proprep.utils.prompts import prompt_with_context, confirm_with_context, int_prompt_with_context, float_prompt_with_context

from proprep.utils.module_registry import ProcessingModule
from proprep.utils.workflow_checklist import WorkflowChecklist, WorkflowStep
from proprep.structure_prep.comprehensive_redox_detector import RedoxSite

from proprep.orca_prep.orca_data_structures import (
    ORCAQMMMSetup,
    compress_indices,
    map_redox_atoms_to_prmtop,
)
from proprep.orca_prep.orca_ff_converter import find_orca_mm, convert_prmtop_to_orcaff
from proprep.orca_prep.orca_writer import ORCAWriter

logger = logging.getLogger(__name__)


# =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
# Workflow Steps
# =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#

ORCA_WORKFLOW_STEPS = [
    WorkflowStep(
        id="orca-1", name="Verify AMBER Topology",
        description="Locate and validate prmtop/rst7 files from TLEaP",
        handler="_step_verify_topology",
        section="Prerequisites",
    ),
    WorkflowStep(
        id="orca-2", name="Convert Force Field (orca_mm)",
        description="Run orca_mm -convff -AMBER to generate ORCA force field file",
        handler="_step_convert_forcefield",
        section="Prerequisites",
        dependencies=["orca-1"],
    ),
    WorkflowStep(
        id="orca-3", name="Select QM Region",
        description="Choose RedoxSite atoms for QM treatment; map to prmtop indices",
        handler="_step_select_qm_region",
        section="QM/MM Configuration",
        dependencies=["orca-1"],
    ),
    WorkflowStep(
        id="orca-4", name="Define Active Region",
        description="Set which atoms are free to move during optimization",
        handler="_step_define_active_region",
        section="QM/MM Configuration",
        dependencies=["orca-3"],
    ),
    WorkflowStep(
        id="orca-5", name="Configure QM Settings",
        description="Set QM method, basis set, charge, and multiplicity",
        handler="_step_configure_qm",
        section="QM/MM Configuration",
        dependencies=["orca-3"],
    ),
    WorkflowStep(
        id="orca-6", name="Configure QM/MM Coupling",
        description="Choose QM/MM scheme, embedding, and charge alteration options",
        handler="_step_configure_coupling",
        section="QM/MM Configuration",
        optional=True,
    ),
    WorkflowStep(
        id="orca-7", name="Configure Job Settings",
        description="Set job type, processors, memory, and extra keywords",
        handler="_step_configure_job",
        section="Job Settings",
        optional=True,
    ),
    WorkflowStep(
        id="orca-8", name="Preview & Write ORCA Input",
        description="Preview the ORCA .inp file, then write to disk",
        handler="_step_write_input",
        section="Output",
        dependencies=["orca-2", "orca-3", "orca-4", "orca-5"],
    ),
]


class ORCAQMMMPreparator(ProcessingModule):
    """Module for generating ORCA QM/MM input files from AMBER topology.

    Typically instantiated and driven by the unified QM/MM Preparator,
    which handles topology loading and frame extraction upstream.
    """

    NAME = "ORCA QM/MM Preparator"
    DESCRIPTION = "Generate ORCA QM/MM input files from RedoxSite and AMBER topology"
    CATEGORY = "preparation"
    VERSION = "1.0.0"

    def __init__(self):
        """Initialize the ORCA QM/MM Preparator module"""
        self.processor = None

        # Configuration storage
        self.selected_redox_site: Optional[RedoxSite] = None
        self.selected_site_indices: List[int] = []

        # AMBER topology
        self._parm: Optional[parmed.Structure] = None
        self._prmtop_path: str = ""
        self._rst7_path: str = ""

        # Multiple frames (trajectory mode)
        # First entry is the "primary" frame used for steps 2-7;
        # all frames get input files in step 8.
        self._frame_rst7_paths: List[str] = []
        self._frame_labels: List[str] = []

        # Flag: True when frame data was provided by the unified QM/MM Preparator
        self._using_upstream_frames: bool = False

        # ORCA setup (populated through workflow)
        self.setup: Optional[ORCAQMMMSetup] = None

        # orca_mm path
        self._orca_mm_path: str = ""

    def set_frame_data(self, prmtop_path: str, rst7_paths: List[str],
                       labels: List[str], pdb_paths: Optional[List[str]] = None):
        """Accept frame data from the unified QM/MM Preparator.

        ``pdb_paths`` is optional — passed by the QM/MM dispatcher so
        the viewer hooks have a structure file to show. If omitted,
        viewer hooks silently no-op (the rest of the workflow runs on
        prmtop+rst7 alone).
        """
        self._prmtop_path = prmtop_path
        self._rst7_path = rst7_paths[0]
        self._frame_rst7_paths = list(rst7_paths)
        self._frame_pdb_paths = list(pdb_paths) if pdb_paths else []
        self._frame_labels = list(labels)
        self._using_upstream_frames = True

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Workspace helpers

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

    def get_workspace_requirements(self) -> List[str]:
        return ["parm7_file", "rst7_file", "detected_redox_sites"]

    def get_workspace_outputs(self) -> List[str]:
        return ["orca_qmmm_setup", "orca_input_file", "orca_ff_file", "orca_pdb_file"]

    def get_menu_options(self) -> Dict[str, str]:
        return {
            "run_workflow": "Run ORCA QM/MM workflow (guided checklist)",
            "view_config": "View current ORCA QM/MM configuration",
            "reset": "Reset configuration to defaults",
        }

    def get_menu_suggestion(self, workspace):
        redox_sites = workspace.get("detected_redox_sites", [])
        prmtop = workspace.get("parm7_file")

        if not redox_sites:
            return "No redox sites detected. Run RedoxSite Detector first."
        if not prmtop:
            return "No prmtop file found. Run Topology Generator first."

        orca_input = workspace.get("orca_input_file")
        if orca_input:
            return "ORCA QM/MM input file written. Setup complete."

        return "Run the ORCA QM/MM workflow (option 1) to prepare input files."

    def handle_menu_option(self, option: str) -> bool:
        if option == "run_workflow":
            return self._run_orca_workflow()
        elif option == "view_config":
            return self._view_current_config()
        elif option == "reset":
            return self._reset_configuration()
        return False

    def _run_orca_workflow(self) -> bool:
        """Run the ORCA QM/MM preparation workflow via interactive checklist."""
        workspace = self.get_workspace()
        base_dir = workspace.get('working_directory', os.getcwd())
        state_dir = Path(base_dir) / "orca_inputs"
        state_dir.mkdir(exist_ok=True)

        checklist = WorkflowChecklist(
            steps=ORCA_WORKFLOW_STEPS,
            executor=self,
            processor=self.processor,
            workflow_name="ORCA QM/MM Preparation",
            console=self.processor.console,
            state_dir=state_dir,
        )
        return checklist.run()

    def _view_current_config(self) -> bool:
        """Display current ORCA QM/MM configuration."""
        console = self.processor.console
        if not self.setup:
            console.print("[yellow]No ORCA QM/MM configuration set up yet.[/yellow]")
            return True

        s = self.setup
        lines = [
            f"[bold]QM Method:[/bold] {s.qm_method}/{s.qm_basis_set}",
            f"[bold]Dispersion:[/bold] {s.dispersion or 'None'}",
            f"[bold]Charge/Mult:[/bold] {s.qm_charge} / {s.qm_multiplicity}",
            f"[bold]QM Atoms:[/bold] {len(s.qm_atom_indices)}",
            f"[bold]Active Atoms:[/bold] {len(s.active_atom_indices)}",
            f"[bold]Embedding:[/bold] {s.embedding}",
            f"[bold]Charge Scheme:[/bold] {s.charge_scheme}",
            f"[bold]Job Type:[/bold] {s.job_type}",
            f"[bold]Processors:[/bold] {s.n_processors}",
            f"[bold]Memory/core:[/bold] {s.memory_mb} MB",
        ]
        if s.rigid_mm_water:
            lines.append("[bold]Rigid MM Water:[/bold] Yes")
        if not s.do_nb_fixed_fixed:
            lines.append("[bold]Frozen-Frozen NB:[/bold] Disabled")
        if s.qmmm_setup_only:
            lines.append("[bold]Mode:[/bold] QMMMSetup (verification only)")
        console.print(Panel("\n".join(lines), title="ORCA QM/MM Configuration", border_style="cyan"))
        return True

    def _reset_configuration(self) -> bool:
        """Reset all configuration to defaults."""
        self.setup = None
        self.selected_redox_site = None
        self.selected_site_indices = []
        self._parm = None
        self.processor.console.print("[green]Configuration reset.[/green]")
        return True

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # WorkflowChecklist Handler Methods
    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#

    def _step_verify_topology(self):
        """Step 1: Verify AMBER topology files exist and are valid."""
        console = self.processor.console

        # If upstream frame extraction already provided the data, skip interaction
        if self._using_upstream_frames:
            console.print(f"\n  Prmtop: [cyan]{os.path.basename(self._prmtop_path)}[/cyan]")
            console.print(f"  Frames: {len(self._frame_rst7_paths)}")

            self._parm = parmed.load_file(self._prmtop_path, self._rst7_path)
            n_atoms = len(self._parm.atoms)
            n_residues = len(self._parm.residues)

            self.setup = ORCAQMMMSetup(
                redox_site=None,
                prmtop_path=self._prmtop_path,
                rst7_path=self._rst7_path,
            )
            return {'summary': f"Topology from upstream: {n_atoms} atoms, {n_residues} residues"}

        # Educational panel
        console.print(Panel(
            "[bold]AMBER Topology Files[/bold]\n\n"
            "The prmtop (parameter/topology) file contains all the information needed\n"
            "for molecular mechanics: atom types, partial charges, bond/angle/dihedral\n"
            "parameters, and van der Waals parameters.\n\n"
            "Coordinates can come from either:\n"
            "  [cyan]rst7/inpcrd[/cyan] — A single structure (e.g., from TLEaP or a restart)\n"
            "  [cyan]trajectory[/cyan]  — Extract frame(s) from an MD trajectory (.nc, .mdcrd, etc.)\n\n"
            "ORCA's orca_mm tool reads the prmtop directly to extract all force field\n"
            "parameters, so you do not need to manually format or convert parameters.",
            title="What are prmtop and rst7 files?",
            border_style="grey50",
            expand=False,
        ))

        # Get prmtop path from workspace or auto-detect or browse
        prmtop_path = self.get_from_workspace("parm7_file")
        if not prmtop_path or not os.path.exists(prmtop_path):
            prmtop_path = self._find_or_browse_file(
                "prmtop", [".prmtop", ".parm7"], "AMBER prmtop file"
            )
        if not os.path.exists(prmtop_path):
            raise FileNotFoundError(f"Prmtop file not found: {prmtop_path}")

        console.print(f"\n  Prmtop: [cyan]{os.path.basename(prmtop_path)}[/cyan]")

        # Choose coordinate source
        coord_source = prompt_with_context(
            self.processor,
            "Coordinate source: \\[r]st7/inpcrd or \\[t]rajectory",
            default="r",
            module="ORCA QM/MM Preparator",
            description="Coordinate source type",
            options_map={"r": "rst7/inpcrd file", "t": "Trajectory"},
        ).strip().lower()

        if coord_source.startswith("t"):
            rst7_path = self._load_from_trajectory(prmtop_path)
        else:
            rst7_path = self._load_rst7(prmtop_path)

        # Load with parmed
        console.print(f"\n  Loading topology + coordinates...")
        self._parm = parmed.load_file(prmtop_path, rst7_path)
        self._prmtop_path = prmtop_path
        self._rst7_path = rst7_path

        # Display summary
        n_atoms = len(self._parm.atoms)
        n_residues = len(self._parm.residues)
        n_bonds = len(self._parm.bonds)

        table = Table(title="Topology Summary")
        table.add_column("Property", style="cyan")
        table.add_column("Value", style="green")
        table.add_row("Atoms", str(n_atoms))
        table.add_row("Residues", str(n_residues))
        table.add_row("Bonds", str(n_bonds))
        table.add_row("Prmtop", os.path.basename(prmtop_path))
        table.add_row("Coordinates", os.path.basename(rst7_path))
        console.print(table)

        # Initialize setup object
        self.setup = ORCAQMMMSetup(
            redox_site=None,  # Set later in step 3
            prmtop_path=prmtop_path,
            rst7_path=rst7_path,
        )

        return {'summary': f"Topology validated: {n_atoms} atoms, {n_residues} residues"}

    def _load_rst7(self, prmtop_path: str) -> str:
        """Load coordinates from an rst7/inpcrd file."""
        rst7_path = self.get_from_workspace("rst7_file")
        if not rst7_path or not os.path.exists(rst7_path):
            rst7_path = self._find_or_browse_file(
                "rst7", [".rst7", ".inpcrd", ".crd"], "AMBER rst7/inpcrd file"
            )
        if not os.path.exists(rst7_path):
            raise FileNotFoundError(f"Rst7 file not found: {rst7_path}")
        return rst7_path

    def _load_from_trajectory(self, prmtop_path: str) -> str:
        """
        Load coordinates from an MD trajectory using cpptraj.

        Workflow:
        1. Select trajectory file(s)
        2. Report frame count
        3. Optionally center solute (autoimage)
        4. Select frames (start:end:increment)
        5. Save selected frames as rst7 files
        6. User picks which frame to use

        Returns:
            Path to the selected rst7 file
        """
        console = self.processor.console

        # Find cpptraj
        cpptraj = shutil.which("cpptraj")
        if not cpptraj:
            raise FileNotFoundError(
                "cpptraj not found on PATH. Install AmberTools to use trajectory loading."
            )

        # Select trajectory file(s)
        traj_extensions = [".nc", ".mdcrd", ".crd", ".dcd", ".xtc", ".trr", ".binpos"]
        console.print("\n  [bold]Select trajectory file(s)[/bold]")
        traj_paths = []

        while True:
            traj_path = self._find_or_browse_file(
                "trajectory", traj_extensions, "MD trajectory file"
            )
            traj_paths.append(traj_path)
            console.print(f"  Added: [cyan]{os.path.basename(traj_path)}[/cyan]")

            add_more = confirm_with_context(
                self.processor,
                "Add another trajectory segment?",
                default=False,
                module="ORCA QM/MM Preparator",
                description="Add another trajectory segment",
            )
            if not add_more:
                break

        # Get frame count using cpptraj
        console.print("\n  Reading trajectory info...")
        n_frames = self._cpptraj_frame_count(cpptraj, prmtop_path, traj_paths)
        console.print(f"  Total frames: [cyan]{n_frames}[/cyan]")

        # Centering / autoimage
        do_autoimage = confirm_with_context(
            self.processor,
            "Center solute in box (autoimage)?",
            default=True,
            module="ORCA QM/MM Preparator",
            description="Apply cpptraj autoimage centering",
        )

        # Frame selection
        console.print(Panel(
            "[bold]Frame Selection[/bold]\n\n"
            "Specify frames as [cyan]start:end:increment[/cyan] (1-based).\n\n"
            "Examples:\n"
            f"  [cyan]1:{n_frames}:1[/cyan]       — All frames\n"
            f"  [cyan]{n_frames}:{n_frames}:1[/cyan]   — Last frame only\n"
            f"  [cyan]1:{n_frames}:100[/cyan]     — Every 100th frame\n"
            f"  [cyan]{n_frames - 9 if n_frames > 9 else 1}:{n_frames}:1[/cyan]      — Last 10 frames",
            title="Which frames to extract?",
            border_style="grey50",
            expand=False,
        ))

        frame_spec = prompt_with_context(
            self.processor,
            f"Frame selection (start:end:increment, 1-{n_frames})",
            default=f"{n_frames}:{n_frames}:1",
            module="ORCA QM/MM Preparator",
            description="Trajectory frame selection",
        ).strip()

        # Parse frame spec
        parts = frame_spec.split(":")
        try:
            start = int(parts[0]) if len(parts) > 0 else 1
            end = int(parts[1]) if len(parts) > 1 else start
            increment = int(parts[2]) if len(parts) > 2 else 1
        except ValueError:
            raise ValueError(f"Invalid frame specification: {frame_spec}. Use start:end:increment")

        if start < 1 or end > n_frames or start > end or increment < 1:
            raise ValueError(
                f"Invalid frame range: {start}:{end}:{increment} "
                f"(trajectory has {n_frames} frames)"
            )

        # Calculate how many frames will be extracted
        n_extract = len(range(start, end + 1, increment))
        console.print(f"\n  Extracting {n_extract} frame(s): {start}:{end}:{increment}")

        # Set up output directory
        workspace = self.get_workspace()
        base_dir = workspace.get('working_directory', os.getcwd())
        frames_dir = Path(base_dir) / "orca_inputs" / "frames"
        frames_dir.mkdir(parents=True, exist_ok=True)

        # Run cpptraj to extract frames as rst7
        prmtop_stem = Path(prmtop_path).stem
        rst7_files = self._cpptraj_extract_frames(
            cpptraj, prmtop_path, traj_paths,
            start, end, increment,
            do_autoimage, frames_dir, prmtop_stem
        )

        if not rst7_files:
            raise RuntimeError("cpptraj produced no output frames")

        frame_numbers = list(range(start, end + 1, increment))
        console.print(f"  [green]Extracted {len(rst7_files)} frame(s) to {frames_dir}[/green]")

        if len(rst7_files) == 1:
            # Single frame
            self._frame_rst7_paths = rst7_files
            self._frame_labels = [f"frame{frame_numbers[0]}"]
            console.print(f"  Using: [cyan]{os.path.basename(rst7_files[0])}[/cyan]")
            return rst7_files[0]

        # Multiple frames — let user select which to generate inputs for
        console.print(f"\n  [bold]Select frame(s) for QM/MM input generation:[/bold]\n")
        for i, (rst7, frame_num) in enumerate(zip(rst7_files, frame_numbers), 1):
            console.print(f"    [{i}] Frame {frame_num}: {os.path.basename(rst7)}")

        console.print(Panel(
            "Select frames using the same syntax as RedoxSite selection:\n"
            f"  [cyan]1[/cyan]         — Single frame\n"
            f"  [cyan]1,3,5[/cyan]     — Specific frames\n"
            f"  [cyan]1-{len(rst7_files)}[/cyan]       — Range of frames\n"
            f"  [cyan]all[/cyan]       — All {len(rst7_files)} frames\n\n"
            "The QM/MM settings (steps 2-7) are configured once using the\n"
            "first selected frame. Step 8 writes input files for all selected frames.",
            title="Multi-frame mode",
            border_style="grey50",
            expand=False,
        ))

        choice = prompt_with_context(
            self.processor,
            f"Select frame(s) (1-{len(rst7_files)}, or 'all')",
            default="all",
            module="ORCA QM/MM Preparator",
            description="Select trajectory frames for ORCA input",
        ).strip()

        selected_indices = self._parse_site_selection(choice, len(rst7_files))
        if selected_indices is None:
            selected_indices = [1]

        # Store selected frames (0-based internally)
        self._frame_rst7_paths = [rst7_files[i - 1] for i in selected_indices]
        self._frame_labels = [f"frame{frame_numbers[i - 1]}" for i in selected_indices]

        if len(self._frame_rst7_paths) > 1:
            console.print(
                f"\n  [green]Selected {len(self._frame_rst7_paths)} frames — "
                f"input files will be generated for each in step 8.[/green]"
            )
        primary = self._frame_rst7_paths[0]
        console.print(f"  Primary frame for configuration: [cyan]{os.path.basename(primary)}[/cyan]")
        return primary

    def _cpptraj_frame_count(self, cpptraj: str, prmtop: str, traj_paths: List[str]) -> int:
        """Get total frame count from trajectory file(s) using cpptraj."""
        cpptraj_input = f"parm {prmtop}\n"
        for tp in traj_paths:
            cpptraj_input += f"trajin {tp}\n"
        cpptraj_input += "nframes\ngo\n"

        result = subprocess.run(
            [cpptraj, "-i", "/dev/stdin"],
            input=cpptraj_input, capture_output=True, text=True, timeout=60
        )

        # Parse frame count from output
        for line in result.stdout.splitlines():
            if "frames" in line.lower() and "total" in line.lower():
                # e.g., "Frames: 1000  Total: 1000"
                for part in line.split():
                    if part.isdigit():
                        return int(part)

        # Fallback: look for "X frames" pattern
        import re
        for line in result.stdout.splitlines():
            match = re.search(r'(\d+)\s+frames', line, re.IGNORECASE)
            if match:
                return int(match.group(1))

        # If we still can't find it, try reading trajin summary
        for line in result.stdout.splitlines():
            if "total read" in line.lower() or "Total frames" in line.lower():
                for part in line.split():
                    if part.isdigit():
                        return int(part)

        raise RuntimeError(
            f"Could not determine frame count from cpptraj output:\n{result.stdout}\n{result.stderr}"
        )

    def _cpptraj_extract_frames(
        self, cpptraj: str, prmtop: str, traj_paths: List[str],
        start: int, end: int, increment: int,
        autoimage: bool, output_dir: Path, stem: str
    ) -> List[str]:
        """Extract trajectory frames as rst7 files using cpptraj."""
        console = self.processor.console

        cpptraj_input = f"parm {prmtop}\n"
        for tp in traj_paths:
            cpptraj_input += f"trajin {tp}\n"

        if autoimage:
            cpptraj_input += "autoimage\n"

        # Write each frame as a separate rst7
        frame_numbers = list(range(start, end + 1, increment))

        if len(frame_numbers) == 1:
            # Single frame — write directly
            out_path = str(output_dir / f"{stem}_frame{frame_numbers[0]}.rst7")
            cpptraj_input += f"trajout {out_path} restart onlyframes {frame_numbers[0]}\n"
        else:
            # Multiple frames — use trajout with frame range
            out_pattern = str(output_dir / f"{stem}_frame.rst7")
            cpptraj_input += f"trajout {out_pattern} restart onlyframes "
            cpptraj_input += ",".join(str(f) for f in frame_numbers) + "\n"

        cpptraj_input += "go\n"

        console.print(f"  [grey50]Running cpptraj...[/grey50]")
        result = subprocess.run(
            [cpptraj, "-i", "/dev/stdin"],
            input=cpptraj_input, capture_output=True, text=True, timeout=300
        )

        if result.returncode != 0:
            raise RuntimeError(f"cpptraj failed:\n{result.stderr}")

        # Find the output rst7 files
        rst7_files = sorted(output_dir.glob(f"{stem}_frame*.rst7"))
        return [str(f) for f in rst7_files]

    def _step_convert_forcefield(self):
        """Step 2: Convert AMBER prmtop to ORCA force field format using orca_mm."""
        console = self.processor.console

        # Educational panel
        console.print(Panel(
            "[bold]Force Field Conversion[/bold]\n\n"
            "The orca_mm utility reads all bonded parameters (bonds, angles, dihedrals),\n"
            "non-bonded parameters (van der Waals, partial charges), and topology\n"
            "information from your AMBER prmtop file. It writes them in ORCA's internal\n"
            "ORCAFF format (.prms file) so that the QM/MM engine can evaluate MM energies\n"
            "and compute the QM/MM coupling terms during the calculation.",
            title="What does orca_mm -convff do?",
            border_style="grey50",
            expand=False,
        ))

        # Find orca_mm
        orca_mm = find_orca_mm()
        if not orca_mm:
            console.print("[yellow]orca_mm not found on PATH.[/yellow]")
            orca_dir = prompt_with_context(
                self.processor,
                "Path to ORCA installation directory (containing orca_mm)",
                default="",
                module="ORCA QM/MM Preparator",
                description="ORCA installation directory",
            )
            orca_mm = find_orca_mm(orca_dir)
            if not orca_mm:
                raise FileNotFoundError(
                    f"orca_mm not found in {orca_dir}. "
                    "Please ensure ORCA is installed and orca_mm is accessible."
                )

        self._orca_mm_path = orca_mm
        console.print(f"  Found orca_mm: [cyan]{orca_mm}[/cyan]")

        # Run conversion
        workspace = self.get_workspace()
        output_dir = str(Path(workspace.get('working_directory', os.getcwd())) / "orca_inputs")
        os.makedirs(output_dir, exist_ok=True)

        console.print(f"\n  Running: [grey50]orca_mm -convff -AMBER {os.path.basename(self._prmtop_path)}[/grey50]")
        success, prms_path, output = convert_prmtop_to_orcaff(
            self._prmtop_path, orca_mm, output_dir=output_dir
        )

        if output:
            console.print(Panel(output.strip(), title="orca_mm output", border_style="grey50", expand=False))

        if not success:
            raise RuntimeError(f"Force field conversion failed:\n{output}")

        self.setup.orca_ff_file = prms_path
        self.update_workspace("orca_ff_file", prms_path)

        console.print(f"\n  [green]✓ Force field file: {os.path.basename(prms_path)}[/green]")
        return {'summary': f"Force field converted: {os.path.basename(prms_path)}"}

    def _step_select_qm_region(self):
        """Step 3: Select RedoxSite and map atoms to prmtop indices."""
        console = self.processor.console

        # Educational panel
        console.print(Panel(
            "[bold]QM Region Selection[/bold]\n\n"
            "In QM/MM, a small region of the system (the QM region) is treated with\n"
            "quantum mechanics, while the rest is treated with molecular mechanics.\n"
            "The QM region should include the chemically active site — where bond\n"
            "breaking/forming, electron transfer, or excited states occur.\n\n"
            "Choosing the QM region carefully is important: too small and you miss\n"
            "important electronic effects; too large and the calculation becomes\n"
            "prohibitively expensive. Where the QM/MM boundary crosses covalent bonds,\n"
            "ORCA automatically inserts hydrogen 'link atoms' to cap the QM region.",
            title="What is the QM region?",
            border_style="grey50",
            expand=False,
        ))

        # Select RedoxSite (reuse ONIOM pattern)
        success = self._select_redox_site()
        if not success:
            raise RuntimeError("RedoxSite selection cancelled or failed")

        # Map to prmtop indices by residue ID
        console.print("\n  Mapping RedoxSite atoms to prmtop indices (by residue ID)...")
        qm_indices, unmapped = map_redox_atoms_to_prmtop(self.selected_redox_site, self._parm)

        if unmapped:
            console.print(f"\n  [yellow]Warning: {len(unmapped)} atom(s) could not be mapped:[/yellow]")
            for name in unmapped[:10]:
                console.print(f"    [yellow]{name}[/yellow]")
            if len(unmapped) > 10:
                console.print(f"    [yellow]... and {len(unmapped) - 10} more[/yellow]")

        if not qm_indices:
            raise RuntimeError("No atoms could be mapped to prmtop indices. "
                             "Check that the RedoxSite residue IDs match the prmtop.")

        self.setup.qm_atom_indices = qm_indices
        self.setup.redox_site = self.selected_redox_site

        # Show QM region summary grouped by residue
        from collections import OrderedDict
        residue_summary = OrderedDict()
        for idx in qm_indices:
            atom = self._parm.atoms[idx]
            res_key = (atom.residue.number, atom.residue.name)
            if res_key not in residue_summary:
                residue_summary[res_key] = {'count': 0, 'heavy': 0, 'first_idx': idx}
            residue_summary[res_key]['count'] += 1
            if atom.atomic_number > 1:
                residue_summary[res_key]['heavy'] += 1

        table = Table(title=f"QM Region: {len(qm_indices)} atoms across {len(residue_summary)} residues")
        table.add_column("Residue", style="green")
        table.add_column("ResID", style="yellow", justify="right")
        table.add_column("Atoms", style="cyan", justify="right")
        table.add_column("Heavy", style="magenta", justify="right")
        table.add_column("Index range", style="grey50")

        for (resnum, resname), info in residue_summary.items():
            # Find index range for this residue within QM atoms
            res_indices = [i for i in qm_indices
                          if self._parm.atoms[i].residue.number == resnum
                          and self._parm.atoms[i].residue.name == resname]
            idx_range = f"{min(res_indices)}–{max(res_indices)}"
            table.add_row(resname, str(resnum), str(info['count']),
                         str(info['heavy']), idx_range)

        console.print(table)

        # Analyze and adjust QM/MM boundary
        qm_indices, boundary_bonds = self._analyze_and_adjust_boundary(qm_indices)

        self.setup.qm_atom_indices = qm_indices
        self.setup.boundary_bonds = boundary_bonds

        # Show charge context
        qm_charge_sum = sum(self._parm.atoms[i].charge for i in qm_indices)
        suggested_charge = round(qm_charge_sum)
        console.print(f"\n  Sum of partial charges in QM region: [cyan]{qm_charge_sum:.4f}[/cyan]")
        console.print(f"  Suggested QM charge: [green]{suggested_charge}[/green]")

        if abs(qm_charge_sum - suggested_charge) > 0.1:
            console.print(
                f"  [yellow]Warning: partial charge sum deviates from integer by "
                f"{abs(qm_charge_sum - suggested_charge):.3f}. "
                f"The charge alteration scheme will redistribute charges at the boundary.[/yellow]"
            )

        self.setup.qm_charge = suggested_charge

        # Recount residues after boundary adjustment
        residues_final = set()
        for idx in qm_indices:
            atom = self._parm.atoms[idx]
            residues_final.add((atom.residue.chain, atom.residue.number, atom.residue.name))

        # Hook B + Hook D: show the final QM region in magenta and
        # the boundary atoms in yellow (QM-side and MM-side both
        # rendered, since both are physically present and the user
        # needs to see where the link-atom hydrogens will land).
        # Same colour convention as ONIOM (magenta HIGH, yellow
        # boundary) so users moving between the two interfaces see
        # consistent semantics.
        try:
            from proprep.structure_prep.viewer_coordinator import (
                viewer as _viewer,
            )

            qm_clauses = []
            seen_qm = set()
            for idx in qm_indices:
                atom = self._parm.atoms[idx]
                resnum = atom.residue.idx + 1
                pair = (resnum, atom.name)
                if pair in seen_qm:
                    continue
                seen_qm.add(pair)
                qm_clauses.append(f"({resnum} and .{atom.name})")
            _viewer.unhighlight("orca_qm")
            if qm_clauses:
                _viewer.highlight(
                    " or ".join(qm_clauses),
                    style="ball+stick",
                    color="#ff00ff",
                    label="orca_qm",
                )

            boundary_clauses = []
            seen_b = set()
            for qm_idx, mm_idx in boundary_bonds:
                for idx in (qm_idx, mm_idx):
                    atom = self._parm.atoms[idx]
                    resnum = atom.residue.idx + 1
                    pair = (resnum, atom.name)
                    if pair in seen_b:
                        continue
                    seen_b.add(pair)
                    boundary_clauses.append(f"({resnum} and .{atom.name})")
            _viewer.unhighlight("orca_boundary")
            if boundary_clauses:
                _viewer.highlight(
                    " or ".join(boundary_clauses),
                    style="ball+stick",
                    color="#ffff00",
                    label="orca_boundary",
                )
        except Exception:
            pass

        return {'summary': f"QM region: {len(qm_indices)} atoms across {len(residues_final)} residues"}

    def _step_define_active_region(self):
        """Step 4: Define which atoms are free to move during optimization."""
        console = self.processor.console

        # Educational panel
        console.print(Panel(
            "[bold]Active Region[/bold]\n\n"
            "During geometry optimization, only 'active' atoms are allowed to move.\n"
            "All QM atoms are always active. You can also include nearby MM atoms so\n"
            "the protein environment can relax around the QM region.\n\n"
            "Atoms far from the QM region are frozen in place. This reduces computational\n"
            "cost and prevents artificial drift of distant residues that have no direct\n"
            "interaction with the active site.\n\n"
            "A common choice is to free all atoms within 6–8 Å of any QM atom.",
            title="Why define an active region?",
            border_style="grey50",
            expand=False,
        ))

        # Prompt for cutoff distance
        cutoff = float(prompt_with_context(
            self.processor,
            "Distance cutoff for active MM atoms (Å)",
            default=str(self.setup.active_region_cutoff),
            module="ORCA QM/MM Preparator",
            description="Active region distance cutoff",
        ))
        self.setup.active_region_cutoff = cutoff

        # Find active atoms: QM atoms + MM atoms within cutoff
        qm_set = set(self.setup.qm_atom_indices)
        active_indices = list(qm_set)

        # Get QM atom coordinates for distance calculation
        qm_coords = []
        for idx in self.setup.qm_atom_indices:
            atom = self._parm.atoms[idx]
            qm_coords.append((atom.xx, atom.xy, atom.xz))

        # Find MM atoms within cutoff
        n_mm_active = 0
        for i, atom in enumerate(self._parm.atoms):
            if i in qm_set:
                continue
            ax, ay, az = atom.xx, atom.xy, atom.xz
            for qx, qy, qz in qm_coords:
                dist = math.sqrt((ax - qx)**2 + (ay - qy)**2 + (az - qz)**2)
                if dist <= cutoff:
                    active_indices.append(i)
                    n_mm_active += 1
                    break

        active_indices = sorted(set(active_indices))
        self.setup.active_atom_indices = active_indices

        n_total = len(self._parm.atoms)
        n_frozen = n_total - len(active_indices)
        n_qm = len(self.setup.qm_atom_indices)

        table = Table(title="Active Region Summary", show_header=False, box=None, padding=(0, 2))
        table.add_column("Category", style="cyan")
        table.add_column("Count", style="green", justify="right")
        table.add_row("QM atoms (always active)", str(n_qm))
        table.add_row(f"MM atoms within {cutoff} Å", str(n_mm_active))
        table.add_row("[bold]Total active[/bold]", f"[bold]{len(active_indices)}[/bold]")
        table.add_row("Frozen", str(n_frozen))
        console.print(table)

        # Hook C: highlight active MM residues in element colours so
        # the user sees how far the cutoff reaches around the QM
        # region. The QM region itself stays magenta from Hook B
        # (different label so it survives this overlay). Frozen
        # residues fall back to the default protein cartoon — visible
        # by absence.
        try:
            from proprep.structure_prep.viewer_coordinator import (
                viewer as _viewer,
            )
            qm_set = set(self.setup.qm_atom_indices)
            active_residues = set()
            for idx in active_indices:
                if idx in qm_set:
                    continue
                atom = self._parm.atoms[idx]
                active_residues.add(atom.residue.idx + 1)
            _viewer.unhighlight("orca_active_mm")
            if active_residues:
                sel = "(" + " or ".join(
                    str(r) for r in sorted(active_residues)
                ) + ")"
                _viewer.highlight(
                    sel, style="ball+stick", color="element",
                    label="orca_active_mm",
                )
        except Exception:
            pass

        return {'summary': f"Active: {len(active_indices)} atoms ({n_qm} QM + {n_mm_active} MM within {cutoff}Å)"}

    def _step_configure_qm(self):
        """Step 5: Configure QM method, basis set, charge, and multiplicity."""
        console = self.processor.console

        # Educational panel
        console.print(Panel(
            "[bold]QM Settings[/bold]\n\n"
            "The QM method and basis set determine the level of theory used to describe\n"
            "the electronic structure of the QM region.\n\n"
            "[cyan]Method[/cyan]: The functional (e.g., B3LYP, PBE0, wB97X-D) determines how\n"
            "electron-electron interactions are approximated. Hybrid functionals like B3LYP\n"
            "are a good balance of accuracy and cost for most systems.\n\n"
            "[cyan]Basis set[/cyan]: Describes the mathematical functions used to represent\n"
            "molecular orbitals. Larger basis sets are more accurate but more expensive.\n"
            "def2-SVP is a reasonable starting point; def2-TZVP is more accurate.\n\n"
            "[cyan]Charge & multiplicity[/cyan]: The total charge and spin state of the QM\n"
            "subsystem. These must match the actual electronic state of the atoms in the\n"
            "QM region. Partial charges at the QM/MM boundary are handled by the charge\n"
            "alteration scheme.",
            title="What do QM settings control?",
            border_style="grey50",
            expand=False,
        ))

        # Method
        self.setup.qm_method = prompt_with_context(
            self.processor,
            "QM method (e.g., B3LYP, PBE0, wB97X-D, DLPNO-CCSD(T))",
            default=self.setup.qm_method,
            module="ORCA QM/MM Preparator",
            description="QM method",
        )

        # Basis set
        self.setup.qm_basis_set = prompt_with_context(
            self.processor,
            "Basis set (e.g., def2-SVP, def2-TZVP, cc-pVDZ)",
            default=self.setup.qm_basis_set,
            module="ORCA QM/MM Preparator",
            description="Basis set",
        )

        # Dispersion
        console.print("\n[grey50]Dispersion corrections account for van der Waals interactions\n"
                      "that DFT functionals often underestimate. D3BJ is recommended for proteins.[/grey50]")
        self.setup.dispersion = prompt_with_context(
            self.processor,
            "Dispersion correction (D3BJ, D3, D4, or leave empty for none)",
            default=self.setup.dispersion,
            module="ORCA QM/MM Preparator",
            description="Dispersion correction",
        )

        # Charge
        console.print(f"\n  QM region partial charge sum: [cyan]{self.setup.qm_charge}[/cyan] (from topology)")
        self.setup.qm_charge = int_prompt_with_context(
            self.processor,
            "QM region charge",
            default=self.setup.qm_charge,
            module="ORCA QM/MM Preparator",
            description="QM region charge",
        )

        # Multiplicity
        console.print("\n[grey50]Multiplicity Guide:[/grey50]")
        console.print("[grey50]  1 = singlet (all electrons paired)[/grey50]")
        console.print("[grey50]  2 = doublet (1 unpaired electron, e.g., Cu²⁺)[/grey50]")
        console.print("[grey50]  3 = triplet (2 unpaired electrons)[/grey50]")

        self.setup.qm_multiplicity = int_prompt_with_context(
            self.processor,
            "QM region multiplicity",
            default=self.setup.qm_multiplicity,
            module="ORCA QM/MM Preparator",
            description="QM region multiplicity",
        )

        console.print(
            f"\n  [green]✓ QM: {self.setup.qm_method}/{self.setup.qm_basis_set}"
            f"{' ' + self.setup.dispersion if self.setup.dispersion else ''}, "
            f"charge={self.setup.qm_charge}, mult={self.setup.qm_multiplicity}[/green]"
        )

        return {
            'summary': f"{self.setup.qm_method}/{self.setup.qm_basis_set}, "
                       f"charge={self.setup.qm_charge}, mult={self.setup.qm_multiplicity}"
        }

    def _step_configure_coupling(self):
        """Step 6: Configure embedding, charge alteration, and boundary settings."""
        console = self.processor.console

        # Educational panel
        console.print(Panel(
            "[bold]QM/MM Coupling[/bold]\n\n"
            "The [cyan]!QMMM[/cyan] keyword uses the additive QM/MM scheme: the total energy\n"
            "is the sum of the QM energy for the QM region, the MM energy for the MM region,\n"
            "and coupling terms between them. This step configures how the QM and MM regions\n"
            "interact with each other.",
            title="Additive QM/MM energy partitioning",
            border_style="grey50",
            expand=False,
        ))

        # Embedding
        console.print(Panel(
            "[bold]Electrostatic Embedding[/bold]\n\n"
            "[cyan]Electrostatic embedding[/cyan] (recommended): The MM point charges are\n"
            "included in the QM Hamiltonian, so the QM electron density is polarized by\n"
            "the protein environment. This is more physically accurate because the QM\n"
            "region 'feels' the electrostatic field of the surrounding residues.\n\n"
            "[cyan]Mechanical embedding[/cyan]: Only van der Waals interactions couple the\n"
            "QM and MM regions. The QM calculation does not see MM charges. Simpler and\n"
            "faster, but misses important polarization effects.",
            title="How does the QM region interact with the MM environment?",
            border_style="grey50",
            expand=False,
        ))

        emb = prompt_with_context(
            self.processor,
            "Embedding: \\[e]lectrostatic or \\[m]echanical",
            default="e",
            module="ORCA QM/MM Preparator",
            description="Embedding scheme selection",
        ).strip().lower()
        self.setup.embedding = "Mechanical" if emb.startswith("m") else "Electrostatic"

        # Charge alteration scheme (only for electrostatic embedding)
        if self.setup.embedding == "Electrostatic":
            console.print(Panel(
                "[bold]Charge Alteration at QM/MM Boundary[/bold]\n\n"
                "Where the QM/MM boundary crosses a covalent bond, the partial charges\n"
                "on nearby MM atoms can cause overpolarization of the QM electron density.\n"
                "Charge alteration schemes redistribute or zero out these problematic charges:\n\n"
                "[cyan]CS[/cyan] (Charge Shift, default): Shifts boundary atom charges to\n"
                "    nearby MM atoms. Most robust for general use.\n"
                "[cyan]RCD[/cyan] (Redistributed Charge & Dipole): Preserves both the total\n"
                "    charge and the dipole moment of the redistributed charges.\n"
                "[cyan]Z0-Z3[/cyan]: Zero-charge schemes that set boundary charges to zero\n"
                "    with increasing levels of compensation (Z0=simplest, Z3=most elaborate).",
                title="Why alter charges at the boundary?",
                border_style="grey50",
                expand=False,
            ))

            scheme_choice = prompt_with_context(
                self.processor,
                "Charge alteration scheme (CS, RCD, Z0, Z1, Z2, Z3)",
                default=self.setup.charge_scheme,
                module="ORCA QM/MM Preparator",
                description="Charge alteration scheme",
            ).strip().upper()
            if scheme_choice in ("CS", "RCD", "Z0", "Z1", "Z2", "Z3"):
                self.setup.charge_scheme = scheme_choice

            # Advanced: charge scheme tuning
            if self.setup.charge_scheme == "CS":
                if confirm_with_context(
                    self.processor,
                    "Tune Scale_CS parameter? (default 0.06, controls charge position on bond)",
                    default=False,
                    module="ORCA QM/MM Preparator",
                    description="Tune Scale_CS charge-scheme parameter",
                ):
                    self.setup.scale_cs = float_prompt_with_context(
                        self.processor,
                        "Scale_CS value",
                        default=0.06,
                        module="ORCA QM/MM Preparator",
                        description="Charge shift scaling factor",
                    )
            elif self.setup.charge_scheme == "RCD":
                if confirm_with_context(
                    self.processor,
                    "Tune Scale_RCD parameter? (default 0.5, midpoint of MM1-MM2 bond)",
                    default=False,
                    module="ORCA QM/MM Preparator",
                    description="Tune Scale_RCD charge-scheme parameter",
                ):
                    self.setup.scale_rcd = float_prompt_with_context(
                        self.processor,
                        "Scale_RCD value",
                        default=0.5,
                        module="ORCA QM/MM Preparator",
                        description="RCD scaling factor",
                    )

        # Link atom distance overrides
        if self.setup.boundary_bonds:
            if confirm_with_context(
                self.processor,
                "Customize link atom (H) distances? (default: C=1.09, N=0.99, O=0.98 Angstrom)",
                default=False,
                module="ORCA QM/MM Preparator",
                description="Customize link atom H distances",
            ):
                self.setup.dist_c_hla = float_prompt_with_context(
                    self.processor, "C-H(link) distance (Angstrom)",
                    default=1.09, module="ORCA QM/MM Preparator",
                    description="Link atom C-H distance",
                )
                self.setup.dist_n_hla = float_prompt_with_context(
                    self.processor, "N-H(link) distance (Angstrom)",
                    default=0.99, module="ORCA QM/MM Preparator",
                    description="Link atom N-H distance",
                )
                self.setup.dist_o_hla = float_prompt_with_context(
                    self.processor, "O-H(link) distance (Angstrom)",
                    default=0.98, module="ORCA QM/MM Preparator",
                    description="Link atom O-H distance",
                )

            # Double-counting control
            if confirm_with_context(
                self.processor,
                "Modify link atom double-counting settings? (advanced, defaults usually fine)",
                default=False,
                module="ORCA QM/MM Preparator",
                description="Modify link atom double-counting settings",
            ):
                console.print(Panel(
                    "[bold]Link Atom Double-Counting[/bold]\n\n"
                    "[cyan]DeleteLADoubleCounting[/cyan] (default: true): Removes angle bends\n"
                    "(QM2-QM1-MM1) and torsions (QM3-QM2-QM1-MM1) that would be double-counted\n"
                    "between the QM and MM calculations at the boundary.\n\n"
                    "[cyan]DeleteLABondDoubleCounting[/cyan] (default: true): Removes the bond\n"
                    "stretching term (QM1-MM1) at the boundary. Usually these should both be true.",
                    title="Boundary double-counting removal",
                    border_style="grey50",
                    expand=False,
                ))
                self.setup.delete_la_double_counting = confirm_with_context(
                    self.processor,
                    "DeleteLADoubleCounting (remove angle/torsion double-counting)?",
                    default=True,
                    module="ORCA QM/MM Preparator",
                    description="Remove link-atom angle/torsion double-counting",
                )
                self.setup.delete_la_bond_double_counting = confirm_with_context(
                    self.processor,
                    "DeleteLABondDoubleCounting (remove bond double-counting)?",
                    default=True,
                    module="ORCA QM/MM Preparator",
                    description="Remove link-atom bond double-counting",
                )

        console.print(
            f"\n  [green]✓ Additive QM/MM, {self.setup.embedding} embedding"
            f"{', ' + self.setup.charge_scheme + ' charge scheme' if self.setup.embedding == 'Electrostatic' else ''}[/green]"
        )

        return {
            'summary': f"{self.setup.embedding} embedding"
                       f"{', ' + self.setup.charge_scheme if self.setup.embedding == 'Electrostatic' else ''}"
        }

    def _step_configure_job(self):
        """Step 7: Configure job type, processors, memory, and extra keywords."""
        console = self.processor.console

        # Educational panel
        console.print(Panel(
            "[bold]Job Settings[/bold]\n\n"
            "[cyan]Job types:[/cyan]\n"
            "  Opt      — Geometry optimization (minimize energy w.r.t. nuclear positions)\n"
            "  EnGrad   — Single-point energy + gradient (no optimization)\n"
            "  NumFreq  — Numerical frequency calculation (vibrational analysis)\n"
            "  SP       — Single-point energy only\n\n"
            "[cyan]Memory:[/cyan] ORCA allocates memory per parallel process. 4000 MB/core\n"
            "is a reasonable default. Increase for large basis sets or many atoms.\n\n"
            "[cyan]Extra keywords:[/cyan] Common options include TightSCF (tighter convergence),\n"
            "Grid5 (finer DFT integration grid), RIJCOSX (faster Coulomb/exchange evaluation).",
            title="What do job settings control?",
            border_style="grey50",
            expand=False,
        ))

        # Job type
        job = prompt_with_context(
            self.processor,
            "Job type (Opt, EnGrad, NumFreq, SP)",
            default=self.setup.job_type,
            module="ORCA QM/MM Preparator",
            description="Job type",
        ).strip()
        self.setup.job_type = job

        # Processors
        self.setup.n_processors = int_prompt_with_context(
            self.processor,
            "Number of processors",
            default=self.setup.n_processors,
            module="ORCA QM/MM Preparator",
            description="Number of processors",
        )

        # Memory
        self.setup.memory_mb = int_prompt_with_context(
            self.processor,
            "Memory per core (MB)",
            default=self.setup.memory_mb,
            module="ORCA QM/MM Preparator",
            description="Memory per core in MB",
        )

        # Additional keywords
        self.setup.additional_keywords = prompt_with_context(
            self.processor,
            "Additional ORCA keywords (e.g., TightSCF Grid5 RIJCOSX, or leave empty)",
            default=self.setup.additional_keywords,
            module="ORCA QM/MM Preparator",
            description="Additional ORCA keywords",
        ).strip()

        # Print level
        self.setup.print_level = int_prompt_with_context(
            self.processor,
            "QM/MM print level (0=minimal, 1=default, 2-4=verbose)",
            default=self.setup.print_level,
            module="ORCA QM/MM Preparator",
            description="QM/MM output verbosity",
        )

        # Rigid MM water
        console.print(Panel(
            "[bold]Rigid MM Water[/bold]\n\n"
            "If your system has explicit water molecules in the active MM region,\n"
            "enabling rigid water constraints holds water bond lengths and angles\n"
            "fixed during geometry optimization. This can improve stability and\n"
            "performance when many solvent molecules are in the active region.\n\n"
            "Only relevant if you have explicit solvent (not implicit/GB).",
            title="Should active MM waters be constrained?",
            border_style="grey50",
            expand=False,
        ))
        self.setup.rigid_mm_water = confirm_with_context(
            self.processor,
            "Use rigid MM water constraints?",
            default=False,
            module="ORCA QM/MM Preparator",
            description="Enable rigid MM water constraints",
        )

        # Frozen-frozen non-bonded interactions
        n_total = len(self._parm.atoms) if self._parm else 0
        n_active = len(self.setup.active_atom_indices)
        n_frozen = n_total - n_active
        default_do_nb = n_total < 10000
        if n_frozen > 0:
            console.print(Panel(
                "[bold]Frozen-Frozen Non-Bonded Interactions[/bold]\n\n"
                "By default, ORCA computes non-bonded interactions between all atom pairs,\n"
                "including pairs where both atoms are frozen (outside the active region).\n"
                f"Your system has {n_frozen} frozen atoms out of {n_total} total.\n\n"
                "For large systems, skipping frozen-frozen interactions saves significant\n"
                "time without affecting the optimized geometry or relative energies.",
                title="Performance option for large systems",
                border_style="grey50",
                expand=False,
            ))
            self.setup.do_nb_fixed_fixed = confirm_with_context(
                self.processor,
                f"Compute frozen-frozen non-bonded interactions?",
                default=default_do_nb,
                module="ORCA QM/MM Preparator",
                description="Compute non-bonded interactions between frozen atoms",
            )

        # QMMMSetup verification mode
        console.print(Panel(
            "[bold]QMMMSetup Verification Mode[/bold]\n\n"
            "When enabled, ORCA reads the input, sets up the QM/MM partitioning,\n"
            "and writes a diagnostic file ([cyan]basename_QM.xyz[/cyan]) showing the\n"
            "final QM region with link atoms, then stops without running the actual\n"
            "QM calculation. Useful for verifying your setup before committing\n"
            "to a long calculation.",
            title="Verify before running?",
            border_style="grey50",
            expand=False,
        ))
        self.setup.qmmm_setup_only = confirm_with_context(
            self.processor,
            "Generate a QMMMSetup verification input (no actual calculation)?",
            default=False,
            module="ORCA QM/MM Preparator",
            description="Generate QMMMSetup verification input only",
        )

        # PDB info vs explicit atom lists
        console.print(Panel(
            "[bold]Atom List Format[/bold]\n\n"
            "ORCA can read QM and active atom assignments from:\n"
            "  [cyan]Explicit lists[/cyan] — QMAtoms/ActiveAtoms blocks in the .inp file\n"
            "  [cyan]PDB columns[/cyan]    — Occupancy=1 for QM, B-factor=1 for active\n\n"
            "Using PDB info keeps the input file shorter and lets you edit the\n"
            "PDB directly to change QM/active regions. Both produce identical results.\n"
            "(The PDB already has occupancy/B-factor set correctly either way.)",
            title="How should ORCA find QM/active atoms?",
            border_style="grey50",
            expand=False,
        ))
        use_pdb_info = confirm_with_context(
            self.processor,
            "Use PDB occupancy/B-factor columns instead of explicit atom lists?",
            default=False,
            module="ORCA QM/MM Preparator",
            description="Use PDB columns for atom classification",
        )
        if use_pdb_info:
            self.setup.use_qm_info_from_pdb = True
            self.setup.use_active_info_from_pdb = True

        # MM cutoffs and dielectric (advanced)
        if confirm_with_context(
            self.processor,
            "Configure MM cutoffs and dielectric constant? (advanced, defaults usually fine)",
            default=False,
            module="ORCA QM/MM Preparator",
            description="Configure MM cutoffs and dielectric constant",
        ):
            console.print(Panel(
                "[bold]MM Non-Bonded Settings[/bold]\n\n"
                "ORCA uses cutoff distances for Coulomb and Lennard-Jones interactions:\n"
                "  [cyan]CoulombCutOff[/cyan]  = 12.0 Angstrom (electrostatic cutoff)\n"
                "  [cyan]LJCutOffInner[/cyan]  = 10.0 Angstrom (full LJ interaction)\n"
                "  [cyan]LJCutOffOuter[/cyan]  = 12.0 Angstrom (LJ smoothly switched off)\n"
                "  [cyan]DielecConst[/cyan]    = 1.0 (dielectric constant for MM electrostatics)\n\n"
                "Increase cutoffs for very polar systems. Set DielecConst > 1 to screen\n"
                "long-range electrostatics (e.g., 4.0 for protein interior, 78.0 for water).",
                title="MM non-bonded settings (%mm block)",
                border_style="grey50",
                expand=False,
            ))

            val = float_prompt_with_context(
                self.processor, "Coulomb cutoff (Angstrom)",
                default=12.0, module="ORCA QM/MM Preparator",
                description="Coulomb interaction cutoff",
            )
            if val != 12.0:
                self.setup.coulomb_cutoff = val

            val = float_prompt_with_context(
                self.processor, "LJ inner cutoff (Angstrom)",
                default=10.0, module="ORCA QM/MM Preparator",
                description="LJ inner cutoff",
            )
            if val != 10.0:
                self.setup.lj_cutoff_inner = val

            val = float_prompt_with_context(
                self.processor, "LJ outer cutoff (Angstrom)",
                default=12.0, module="ORCA QM/MM Preparator",
                description="LJ outer cutoff",
            )
            if val != 12.0:
                self.setup.lj_cutoff_outer = val

            val = float_prompt_with_context(
                self.processor, "Dielectric constant",
                default=1.0, module="ORCA QM/MM Preparator",
                description="MM dielectric constant",
            )
            if val != 1.0:
                self.setup.dielec_const = val

        # Live preview
        writer = ORCAWriter(self.setup)
        simple_line = writer.build_simple_input_line()
        console.print(f"\n  [cyan]Simple input line preview:[/cyan]")
        console.print(f"  ! {simple_line}")

        summary = f"{job}, {self.setup.n_processors} cores, {self.setup.memory_mb}MB/core"
        if self.setup.additional_keywords:
            summary += f", {self.setup.additional_keywords}"
        if self.setup.qmmm_setup_only:
            summary += " (QMMMSetup verification)"
        return {'summary': summary}

    def _step_write_input(self):
        """Step 8: Preview and write the ORCA input file(s)."""
        console = self.processor.console

        # Prepare output directory
        workspace = self.get_workspace()
        base_dir = workspace.get('working_directory', os.getcwd())
        output_dir = Path(base_dir) / "orca_inputs"
        output_dir.mkdir(exist_ok=True)

        prmtop_stem = Path(self.setup.prmtop_path).stem
        orca_dir = os.path.dirname(self._orca_mm_path) if self._orca_mm_path else ""

        # Copy .prms file to output dir if not already there
        prms_basename = os.path.basename(self.setup.orca_ff_file)
        prms_dest = str(output_dir / prms_basename)
        if self.setup.orca_ff_file != prms_dest:
            shutil.copy2(self.setup.orca_ff_file, prms_dest)
            console.print(f"  [green]✓ FF file: {prms_dest}[/green]")

        # Build list of (label, rst7_path) for each frame to write
        if len(self._frame_rst7_paths) > 1:
            frames = list(zip(self._frame_labels, self._frame_rst7_paths))
        else:
            # Single frame (rst7 mode or single trajectory frame)
            frames = [("", self._rst7_path)]

        multi_frame = len(frames) > 1

        # Preview using the first frame
        first_label, first_rst7 = frames[0]
        suffix = f"_{first_label}" if first_label else ""
        preview_pdb = f"{prmtop_stem}_qmmm{suffix}.pdb"
        self.setup.pdb_path = str(output_dir / preview_pdb)

        writer = ORCAWriter(self.setup)
        preview_inp = f"{prmtop_stem}_qmmm{suffix}.inp"
        content = writer.generate_input_content(
            title=f"QM/MM calculation for {prmtop_stem}{' ' + first_label if first_label else ''}"
        )

        console.print(Panel(
            content,
            title=f"[bold]{preview_inp}[/bold]"
                  + (f" [grey50](preview — {len(frames)} frames total)[/grey50]" if multi_frame else ""),
            border_style="cyan",
            expand=False,
        ))

        if multi_frame:
            console.print(
                f"\n  [cyan]Multi-frame mode:[/cyan] Will generate {len(frames)} input file set(s) "
                f"with identical QM/MM settings."
            )

        # Confirm write
        if not confirm_with_context(
            self.processor,
            f"Write {'all ' + str(len(frames)) + ' input files' if multi_frame else 'this input file'}?",
            default=True,
        ):
            raise RuntimeError("Write cancelled by user")

        # Write files for each frame
        all_inp_files = []
        all_script_files = []

        for label, rst7_path in frames:
            suffix = f"_{label}" if label else ""
            inp_filename = f"{prmtop_stem}_qmmm{suffix}.inp"
            pdb_filename = f"{prmtop_stem}_qmmm{suffix}.pdb"
            script_filename = f"run_{prmtop_stem}_qmmm{suffix}.sh"

            inp_path = str(output_dir / inp_filename)
            pdb_path = str(output_dir / pdb_filename)
            script_path = str(output_dir / script_filename)

            # Update setup paths for this frame
            self.setup.pdb_path = pdb_path
            self.setup.rst7_path = rst7_path

            frame_writer = ORCAWriter(self.setup)
            frame_title = f"QM/MM calculation for {prmtop_stem}"
            if label:
                frame_title += f" {label}"

            frame_writer.write_input_file(inp_path, title=frame_title)
            frame_writer.prepare_pdb_file(pdb_path, rst7_path=rst7_path)
            frame_writer.write_run_script(script_path, inp_filename, orca_dir=orca_dir)

            all_inp_files.append(inp_filename)
            all_script_files.append(script_filename)

            if multi_frame:
                console.print(f"  [green]✓ {label}:[/green] {inp_filename}")
            else:
                console.print(f"  [green]✓ Input file: {inp_path}[/green]")
                console.print(f"  [green]✓ PDB file:   {pdb_path}[/green]")
                console.print(f"    (occupancy=1.0 for {len(self.setup.qm_atom_indices)} QM atoms, "
                              f"B-factor=1.0 for {len(self.setup.active_atom_indices)} active atoms)")
                console.print(f"  [green]✓ Run script: {script_path}[/green]")

        # Write master run-all script for multi-frame
        if multi_frame:
            master_script = str(output_dir / f"run_all_{prmtop_stem}_qmmm.sh")
            self._write_master_run_script(master_script, all_script_files, output_dir)
            console.print(f"\n  [green]✓ Master script: {master_script}[/green]")

        # Restore primary rst7 path
        self.setup.rst7_path = self._rst7_path

        # Educational summary
        if multi_frame:
            console.print(Panel(
                f"[bold]How to run these calculations[/bold]\n\n"
                f"Run all frames sequentially:\n"
                f"  cd {output_dir}\n"
                f"  ./run_all_{prmtop_stem}_qmmm.sh\n\n"
                f"Or run a single frame:\n"
                f"  ./{all_script_files[0]}\n\n"
                f"[bold]Monitor progress:[/bold]\n"
                f"  tail -f {prmtop_stem}_qmmm_*.out",
                title=f"Next steps ({len(frames)} frames)",
                border_style="green",
                expand=False,
            ))
        else:
            console.print(Panel(
                f"[bold]How to run this calculation[/bold]\n\n"
                f"  cd {output_dir}\n"
                f"  ./{all_script_files[0]}\n\n"
                f"The script runs ORCA in the background (nohup) and writes\n"
                f"the PID to {prmtop_stem}_qmmm.pid for monitoring.\n\n"
                f"[bold]Monitor progress:[/bold]\n"
                f"  tail -f {prmtop_stem}_qmmm.out\n\n"
                f"[bold]Output files to check:[/bold]\n"
                f"  {prmtop_stem}_qmmm.out     — Main output (energies, convergence)\n"
                f"  {prmtop_stem}_qmmm.xyz     — Optimized geometry (if Opt job)\n"
                f"  {prmtop_stem}_qmmm.gbw     — Wavefunction file (for restarts)",
                title="Next steps",
                border_style="green",
                expand=False,
            ))

        # Update workspace
        self.update_workspace("orca_input_file", str(output_dir / all_inp_files[0]))
        self.update_workspace("orca_pdb_file", self.setup.pdb_path)
        self.update_workspace("orca_qmmm_setup", self.setup)

        if multi_frame:
            return {'summary': f"Written {len(frames)} input file set(s)"}
        return {'summary': f"Written: {all_inp_files[0]}"}

    def _write_master_run_script(
        self, output_path: str, script_files: List[str], output_dir: Path
    ) -> None:
        """Write a master script that runs all frame scripts sequentially."""
        lines = [
            "#!/bin/bash",
            "# Master ORCA QM/MM run script — runs all frames sequentially",
            "# Generated by ProPrep",
            "",
            'set -euo pipefail',
            'SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"',
            'cd "$SCRIPT_DIR"',
            "",
            f'echo "Running {len(script_files)} QM/MM calculations sequentially"',
            f'echo "Started: $(date)"',
            'echo ""',
            "",
        ]

        for i, script in enumerate(script_files, 1):
            lines.append(f'echo "=== [{i}/{len(script_files)}] {script} ==="')
            lines.append(f'bash "./{script}"')
            lines.append('echo ""')
            lines.append("")

        lines.extend([
            'echo "======================================"',
            'echo " All calculations complete."',
            'echo " Ended: $(date)"',
            'echo "======================================"',
        ])

        with open(output_path, 'w') as f:
            f.write("\n".join(lines) + "\n")
        os.chmod(output_path, 0o755)

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Helper Methods
    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#

    def _select_redox_site(self) -> bool:
        """Select RedoxSite from workspace. Reuses ONIOM selection pattern."""
        console = self.processor.console
        redox_sites = self.get_from_workspace("detected_redox_sites")

        if not redox_sites:
            console.print("\n[red]No RedoxSite objects found in workspace.[/red]")
            console.print("Please run 'Comprehensive RedoxSite Detector' first.")
            return False

        console.print(f"\n[bold cyan]Found {len(redox_sites)} RedoxSite(s)[/bold cyan]\n")

        # Display available sites
        table = Table(title="Available RedoxSites")
        table.add_column("Index", style="cyan", width=6)
        table.add_column("Redox Centers", style="green")
        table.add_column("Residues", style="yellow")
        table.add_column("Atoms", style="magenta")
        table.add_column("Bonds", style="blue")

        for i, site in enumerate(redox_sites, 1):
            center_description = []
            if hasattr(site, 'centers') and site.centers:
                metal_centers = [c for c in site.centers if c.center_type.value == "metal_ion"]
                organometallic_centers = [c for c in site.centers if c.center_type.value == "organometallic_cofactor"]
                organic_centers = [c for c in site.centers if c.center_type.value == "organic_cofactor"]
                amino_centers = [c for c in site.centers if c.center_type.value == "redox_amino_acid"]

                if metal_centers:
                    elements = ", ".join(sorted(set(c.element for c in metal_centers if c.element)))
                    center_description.append(f"{len(metal_centers)} metal ({elements})")
                if organometallic_centers:
                    names = ", ".join(sorted(set(c.resname for c in organometallic_centers)))
                    center_description.append(f"{len(organometallic_centers)} organometallic ({names})")
                if organic_centers:
                    resnames = ", ".join(sorted(set(c.resname for c in organic_centers)))
                    center_description.append(f"{len(organic_centers)} cofactor ({resnames})")
                if amino_centers:
                    resnames = ", ".join(sorted(set(c.resname for c in amino_centers)))
                    center_description.append(f"{len(amino_centers)} amino ({resnames})")

            center_str = "; ".join(center_description) if center_description else "No centers"
            residues = set((atom.chain, atom.resid, atom.insertion_code) for atom in site.atoms)

            table.add_row(
                str(i),
                center_str,
                str(len(residues)),
                str(len(site.atoms)),
                str(len(site.bonds))
            )

        console.print(table)

        # Hook A: highlight every redox site listed in the table
        # palette-coloured by row number — mirrors the ONIOM picker
        # behaviour. Uses the upstream-extracted PDB so coloured
        # residues line up with what the user types in.
        try:
            from proprep.structure_prep.viewer_coordinator import (
                viewer as _viewer,
            )
            pdb_file = (
                self._frame_pdb_paths[0]
                if getattr(self, '_frame_pdb_paths', None)
                else None
            )
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
                label = f"orca_pick_site_{i}"
                _viewer.highlight(
                    " or ".join(clauses),
                    style="ball+stick",
                    color=f"palette:{i}",
                    label=label,
                )
                applied.append(label)
            self._orca_pick_labels = applied
        except Exception:
            self._orca_pick_labels = []

        # Prompt user to select site(s)
        while True:
            try:
                selection_str = prompt_with_context(
                    self.processor,
                    "\nSelect RedoxSite(s) for QM region (e.g. 1, 1-3, 1,3,5, or all)",
                    default="1",
                    module="ORCA QM/MM Preparator",
                    description="Select RedoxSite(s) for QM region",
                )

                indices = self._parse_site_selection(selection_str, len(redox_sites))
                if indices is None:
                    console.print(f"[red]Invalid selection. Please choose from 1-{len(redox_sites)}[/red]")
                    continue

                self.selected_site_indices = [i - 1 for i in indices]

                if len(indices) == 1:
                    self.selected_redox_site = redox_sites[self.selected_site_indices[0]]
                    console.print(
                        f"\n[green]✓ Selected RedoxSite {indices[0]} with "
                        f"{len(self.selected_redox_site.atoms)} atoms[/green]"
                    )
                else:
                    self.selected_redox_site = self._merge_redox_sites(
                        [redox_sites[i] for i in self.selected_site_indices]
                    )
                    site_labels = ", ".join(str(i) for i in indices)
                    console.print(
                        f"\n[green]✓ Selected RedoxSites {site_labels} "
                        f"(merged: {len(self.selected_redox_site.atoms)} atoms)[/green]"
                    )
                self._clear_orca_pick_labels()
                return True

            except KeyboardInterrupt:
                self._clear_orca_pick_labels()
                console.print("[yellow]Selection cancelled[/yellow]")
                return False

    def _clear_orca_pick_labels(self) -> None:
        """Drop the per-row redox-site pick highlights (Hook A)."""
        labels = getattr(self, "_orca_pick_labels", None) or []
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
        self._orca_pick_labels = []

    @staticmethod
    def _parse_site_selection(selection_str: str, n_sites: int) -> Optional[List[int]]:
        """Parse site selection string into sorted list of unique 1-based indices."""
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
        seen_center_coords = set()
        for site in sites:
            for center in site.centers:
                if center.coords not in seen_center_coords:
                    merged.centers.append(center)
                    seen_center_coords.add(center.coords)
            for atom in site.atoms:
                if atom.coords not in seen_coords:
                    merged.atoms.append(atom)
                    seen_coords.add(atom.coords)
            for bond in site.bonds:
                merged.bonds.append(bond)
            merged.coord_to_pdb.update(site.coord_to_pdb)
            for key, coords_list in site.residue_groups.items():
                if key not in merged.residue_groups:
                    merged.residue_groups[key] = []
                merged.residue_groups[key].extend(coords_list)
        return merged

    def _find_boundary_bonds(self, qm_indices: List[int]) -> List[Tuple[int, int]]:
        """
        Find bonds that cross the QM/MM boundary.

        Returns list of (qm_atom_index, mm_atom_index) tuples.
        """
        qm_set = set(qm_indices)
        boundary = []
        for bond in self._parm.bonds:
            idx1 = bond.atom1.idx
            idx2 = bond.atom2.idx
            if idx1 in qm_set and idx2 not in qm_set:
                boundary.append((idx1, idx2))
            elif idx2 in qm_set and idx1 not in qm_set:
                boundary.append((idx2, idx1))
        return boundary

    def _classify_boundary_bond(self, qm_idx: int, mm_idx: int) -> str:
        """
        Classify a boundary bond for display and adjustment recommendations.

        Returns one of:
            'backbone_amide' - Polar N-C amide bond (should be moved)
            'backbone_ca_c'  - Cα-C bond (acceptable)
            'backbone_n_ca'  - N-Cα bond (acceptable)
            'sidechain'      - Cα-Cβ or similar (acceptable)
            'other'          - Unclassified
        """
        qm_atom = self._parm.atoms[qm_idx]
        mm_atom = self._parm.atoms[mm_idx]

        qm_name = qm_atom.name.strip()
        mm_name = mm_atom.name.strip()

        # Backbone amide: N—C or C—N across residues
        if (qm_name == 'N' and mm_name == 'C') or (qm_name == 'C' and mm_name == 'N'):
            return 'backbone_amide'
        # Cα—C bond
        if (qm_name == 'CA' and mm_name == 'C') or (qm_name == 'C' and mm_name == 'CA'):
            return 'backbone_ca_c'
        # N—Cα bond
        if (qm_name == 'N' and mm_name == 'CA') or (qm_name == 'CA' and mm_name == 'N'):
            return 'backbone_n_ca'
        # Cα—Cβ
        if (qm_name == 'CA' and mm_name == 'CB') or (qm_name == 'CB' and mm_name == 'CA'):
            return 'sidechain'

        return 'other'

    def _get_residue_atom_indices(self, residue) -> List[int]:
        """Get all atom indices for a given residue object."""
        return [atom.idx for atom in residue.atoms]

    def _analyze_and_adjust_boundary(self, qm_indices: List[int]) -> Tuple[List[int], List[Tuple[int, int]]]:
        """
        Analyze QM/MM boundary bonds and let the user adjust boundary placement.

        For each boundary bond, classifies it and offers adjustment options:
        - Backbone amide N-C: flagged as problematic, offer to extend to Cα-C or Cα-N
        - Sidechain Cα-Cβ: acceptable, but user can include the whole residue
        - Cα-C or N-Cα: acceptable for link atom placement

        Returns:
            Tuple of (adjusted qm_indices, final boundary_bonds)
        """
        console = self.processor.console
        qm_set = set(qm_indices)

        boundary_bonds = self._find_boundary_bonds(qm_indices)
        if not boundary_bonds:
            console.print("\n  [green]No QM/MM boundary bonds (entire system or isolated molecules).[/green]")
            return qm_indices, boundary_bonds

        # Classify boundary bonds
        console.print(Panel(
            "[bold]QM/MM Boundary Analysis[/bold]\n\n"
            "Link atoms are placed where covalent bonds cross the QM/MM boundary.\n"
            "The bond that is cut should be non-polar and not involved in conjugation.\n\n"
            "[green]Good cuts:[/green]  Cα—C(=O) or N—Cα (standard backbone cuts)\n"
            "[green]            Cα—Cβ (sidechain truncation)[/green]\n"
            "[red]Bad cuts:[/red]   Amide N—C(=O) bonds (polar, partial double-bond character)",
            title="Where should the QM/MM boundary be?",
            border_style="grey50",
            expand=False,
        ))

        # Display boundary bonds with classification
        table = Table(title=f"QM/MM Boundary: {len(boundary_bonds)} bond(s)")
        table.add_column("#", style="cyan", width=4)
        table.add_column("QM atom", style="green")
        table.add_column("MM atom", style="yellow")
        table.add_column("Bond type", style="magenta")
        table.add_column("Status", style="bold")

        problematic = []
        for i, (qm_idx, mm_idx) in enumerate(boundary_bonds, 1):
            qm_atom = self._parm.atoms[qm_idx]
            mm_atom = self._parm.atoms[mm_idx]
            bond_type = self._classify_boundary_bond(qm_idx, mm_idx)

            qm_label = f"{qm_atom.name}({qm_atom.residue.name}{qm_atom.residue.number})"
            mm_label = f"{mm_atom.name}({mm_atom.residue.name}{mm_atom.residue.number})"

            if bond_type == 'backbone_amide':
                status = "[red]Move boundary[/red]"
                problematic.append((i, qm_idx, mm_idx, bond_type))
            elif bond_type in ('backbone_ca_c', 'backbone_n_ca', 'sidechain'):
                status = "[green]OK[/green]"
            else:
                status = "[yellow]Review[/yellow]"
                problematic.append((i, qm_idx, mm_idx, bond_type))

            type_labels = {
                'backbone_amide': 'Amide N—C',
                'backbone_ca_c': 'Cα—C',
                'backbone_n_ca': 'N—Cα',
                'sidechain': 'Cα—Cβ',
                'other': 'Other',
            }
            table.add_row(str(i), qm_label, mm_label, type_labels[bond_type], status)

        console.print(table)

        if not problematic:
            console.print("\n  [green]All boundary bonds are at acceptable positions.[/green]")
            return qm_indices, boundary_bonds

        # Offer adjustment for problematic bonds
        console.print(
            f"\n  [yellow]{len(problematic)} bond(s) need attention.[/yellow]"
        )

        for num, qm_idx, mm_idx, bond_type in problematic:
            qm_atom = self._parm.atoms[qm_idx]
            mm_atom = self._parm.atoms[mm_idx]
            qm_label = f"{qm_atom.name}({qm_atom.residue.name}{qm_atom.residue.number})"
            mm_label = f"{mm_atom.name}({mm_atom.residue.name}{mm_atom.residue.number})"

            console.print(f"\n  [bold]Bond #{num}:[/bold] {qm_label} — {mm_label}")

            options = self._get_boundary_options(qm_idx, mm_idx, bond_type, qm_set)
            if not options:
                console.print("    [grey50]No automatic adjustment available.[/grey50]")
                continue

            for opt_num, (description, _, _) in enumerate(options, 1):
                console.print(f"    [{opt_num}] {description}")
            console.print(f"    [0] Keep as-is")

            choice = prompt_with_context(
                self.processor,
                f"  Choice for bond #{num}",
                default="1",
                module="ORCA QM/MM Preparator",
                description="Boundary adjustment",
            ).strip()

            try:
                choice_idx = int(choice)
                if 1 <= choice_idx <= len(options):
                    _, add_atoms, remove_atoms = options[choice_idx - 1]
                    qm_set.update(add_atoms)
                    qm_set -= set(remove_atoms)
                    if add_atoms:
                        console.print(f"    [green]Added {len(add_atoms)} atom(s) to QM region[/green]")
                    if remove_atoms:
                        console.print(f"    [yellow]Removed {len(remove_atoms)} atom(s) from QM region[/yellow]")
            except ValueError:
                pass

        # Recompute after adjustments
        qm_indices = sorted(qm_set)
        boundary_bonds = self._find_boundary_bonds(qm_indices)

        # Show final boundary
        console.print(f"\n  [cyan]Final QM/MM boundary:[/cyan] {len(boundary_bonds)} bond(s)")
        for qm_idx, mm_idx in boundary_bonds:
            qm_atom = self._parm.atoms[qm_idx]
            mm_atom = self._parm.atoms[mm_idx]
            bond_type = self._classify_boundary_bond(qm_idx, mm_idx)
            type_labels = {
                'backbone_amide': '[red]Amide N—C[/red]',
                'backbone_ca_c': 'Cα—C',
                'backbone_n_ca': 'N—Cα',
                'sidechain': 'Cα—Cβ',
                'other': 'Other',
            }
            console.print(
                f"    {qm_atom.name}({qm_atom.residue.name}{qm_atom.residue.number}) — "
                f"{mm_atom.name}({mm_atom.residue.name}{mm_atom.residue.number})  "
                f"[grey50]({type_labels.get(bond_type, bond_type)})[/grey50]"
            )

        console.print(f"\n  [green]QM region: {len(qm_indices)} atoms[/green]")
        return qm_indices, boundary_bonds

    def _get_backbone_atoms(self, residue) -> List[int]:
        """Get indices of backbone atoms (N, H, CA, HA, C, O) in a residue."""
        backbone_names = {'N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'HA2', 'HA3', 'C', 'O', 'OXT'}
        return [a.idx for a in residue.atoms if a.name.strip() in backbone_names]

    def _get_sidechain_atoms(self, residue) -> List[int]:
        """Get indices of sidechain atoms (everything except backbone) in a residue."""
        backbone_names = {'N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'HA2', 'HA3', 'C', 'O', 'OXT'}
        return [a.idx for a in residue.atoms if a.name.strip() not in backbone_names]

    def _has_sidechain(self, residue) -> bool:
        """Check if a residue has a sidechain (CB atom)."""
        return any(a.name.strip() == 'CB' for a in residue.atoms)

    def _get_boundary_options(
        self, qm_idx: int, mm_idx: int, bond_type: str, qm_set: set
    ) -> List[Tuple[str, List[int], List[int]]]:
        """
        Generate adjustment options for a boundary bond.

        Each option is a (description, atoms_to_add, atoms_to_remove) tuple.
        """
        qm_atom = self._parm.atoms[qm_idx]
        mm_atom = self._parm.atoms[mm_idx]
        mm_res = mm_atom.residue
        qm_res = qm_atom.residue

        options = []

        if bond_type == 'backbone_amide':
            # Option 1: Extend across the amide to make the cut at Cα-C or N-Cα
            if qm_atom.name.strip() == 'C' and mm_atom.name.strip() == 'N':
                # QM C(resX) — MM N(resX+1): add N,H so cut moves to N—Cα
                atoms_to_add = []
                for a in mm_res.atoms:
                    if a.name.strip() in ('N', 'H', 'H1', 'H2', 'H3'):
                        atoms_to_add.append(a.idx)
                if atoms_to_add:
                    options.append((
                        f"Extend QM: add N,H of {mm_res.name}{mm_res.number} "
                        f"(cut moves to N—Cα of {mm_res.name}{mm_res.number})",
                        atoms_to_add, []
                    ))

            elif qm_atom.name.strip() == 'N' and mm_atom.name.strip() == 'C':
                # QM N(resX) — MM C(resX-1): add C,O so cut moves to Cα—C
                atoms_to_add = []
                for a in mm_res.atoms:
                    if a.name.strip() in ('C', 'O'):
                        atoms_to_add.append(a.idx)
                if atoms_to_add:
                    options.append((
                        f"Extend QM: add C,O of {mm_res.name}{mm_res.number} "
                        f"(cut moves to Cα—C of {mm_res.name}{mm_res.number})",
                        atoms_to_add, []
                    ))

            # Option 2: Cα-Cβ cut — keep only sidechain of the QM residue
            if self._has_sidechain(qm_res):
                backbone_in_qm = [i for i in self._get_backbone_atoms(qm_res) if i in qm_set]
                sidechain = self._get_sidechain_atoms(qm_res)
                if backbone_in_qm and sidechain:
                    options.append((
                        f"Sidechain only: remove backbone of {qm_res.name}{qm_res.number} "
                        f"(cut moves to Cα—Cβ, removes {len(backbone_in_qm)} backbone atoms)",
                        [], backbone_in_qm
                    ))

            # Option 3: Include the entire flanking residue
            flanking_atoms = [a.idx for a in mm_res.atoms if a.idx not in qm_set]
            if flanking_atoms:
                options.append((
                    f"Include entire {mm_res.name}{mm_res.number} in QM "
                    f"({len(flanking_atoms)} atoms, no boundary here)",
                    flanking_atoms, []
                ))

        elif bond_type == 'other':
            if qm_atom.residue == mm_atom.residue:
                # Intra-residue boundary — incomplete residue
                rest = [a.idx for a in mm_atom.residue.atoms if a.idx not in qm_set]
                if rest:
                    options.append((
                        f"Complete residue {mm_res.name}{mm_res.number} "
                        f"({len(rest)} missing atoms)",
                        rest, []
                    ))
            else:
                # Cα-Cβ cut if the QM residue has a sidechain
                if self._has_sidechain(qm_res):
                    backbone_in_qm = [i for i in self._get_backbone_atoms(qm_res) if i in qm_set]
                    sidechain = self._get_sidechain_atoms(qm_res)
                    if backbone_in_qm and sidechain:
                        options.append((
                            f"Sidechain only: remove backbone of {qm_res.name}{qm_res.number} "
                            f"(cut moves to Cα—Cβ, removes {len(backbone_in_qm)} backbone atoms)",
                            [], backbone_in_qm
                        ))

                # Include the entire MM-side residue
                flanking_atoms = [a.idx for a in mm_res.atoms if a.idx not in qm_set]
                if flanking_atoms:
                    options.append((
                        f"Include entire {mm_res.name}{mm_res.number} in QM "
                        f"({len(flanking_atoms)} atoms)",
                        flanking_atoms, []
                    ))

        return options

    def _find_or_browse_file(self, file_type: str, extensions: List[str], description: str) -> str:
        """
        Find a file by auto-detecting in working directory, or offer a file browser.

        1. Scans working directory for files matching the given extensions
        2. If exactly one match, offers it as default
        3. If multiple matches, lets user pick from a numbered list
        4. If no matches, falls back to a directory browser

        Args:
            file_type: Short label (e.g. "prmtop", "rst7")
            extensions: List of extensions to match (e.g. [".prmtop", ".parm7"])
            description: Human-readable description for prompts
        """
        console = self.processor.console
        workspace = self.get_workspace()
        work_dir = Path(workspace.get('working_directory', os.getcwd()))

        # Auto-detect matching files in working directory
        matches = []
        for f in sorted(work_dir.iterdir()):
            if f.is_file() and f.suffix in extensions:
                matches.append(f)

        if len(matches) == 1:
            console.print(f"  Auto-detected {file_type}: [cyan]{matches[0].name}[/cyan]")
            use_it = confirm_with_context(
                self.processor,
                f"Use {matches[0].name}?",
                default=True,
                module="ORCA QM/MM Preparator",
                description=f"Confirm auto-detected {description}",
            )
            if use_it:
                return str(matches[0])

        elif len(matches) > 1:
            console.print(f"\n  Found {len(matches)} {file_type} file(s) in working directory:\n")
            for i, f in enumerate(matches, 1):
                console.print(f"    [{i}] {f.name}")

            choice = prompt_with_context(
                self.processor,
                f"Select {file_type} file (number), or 'b' to browse",
                default="1",
                module="ORCA QM/MM Preparator",
                description=f"Select {description}",
            ).strip()

            if choice != 'b':
                try:
                    idx = int(choice) - 1
                    if 0 <= idx < len(matches):
                        return str(matches[idx])
                except ValueError:
                    pass

        # File browser fallback
        return self._browse_for_file(file_type, extensions, work_dir)

    def _browse_for_file(self, file_type: str, extensions: List[str], start_dir: Path) -> str:
        """Interactive directory browser for selecting a file."""
        from datetime import datetime
        from proprep.utils.file_browser import file_browser

        def _date_detail(p):
            try:
                return datetime.fromtimestamp(os.path.getmtime(p)).strftime("%m/%d/%Y")
            except OSError:
                return ""

        # Shared browser: bare-N, q to cancel, path jump, filename-based replay.
        # Preserve the historical contract: cancelling raises FileNotFoundError.
        result = file_browser(
            directory=str(start_dir),
            extensions=extensions,
            console=self.processor.console,
            processor=self.processor,
            label=f"{file_type} file",
            entry_detail=_date_detail,
            allow_path_jump=True,
            module="ORCA QM/MM Preparator",
        )
        if result is None:
            raise FileNotFoundError(f"No {file_type} file selected")
        return result
