"""
Shared frame extraction and file selection utilities for QM/MM preparation.

Handles topology loading, coordinate source selection (rst7 or trajectory),
and cpptraj-based frame extraction with both rst7 and PDB output.
"""

import logging
import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

from rich.panel import Panel
from rich.table import Table
from proprep.utils.prompts import prompt_with_context, confirm_with_context

logger = logging.getLogger(__name__)

MODULE_NAME = "QM/MM Preparator"


@dataclass
class FrameExtractionResult:
    """Result of the shared frame extraction phase."""
    prmtop_path: str
    frame_rst7_paths: List[str] = field(default_factory=list)
    frame_pdb_paths: List[str] = field(default_factory=list)
    frame_labels: List[str] = field(default_factory=list)


class FrameExtractor:
    """Load AMBER topology, select coordinate source, and extract frames.

    Provides the shared first phase for both Gaussian (ONIOM) and ORCA
    QM/MM preparation workflows.
    """

    def __init__(self, processor, workspace):
        self.processor = processor
        self.console = processor.console
        self.workspace = workspace

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def run(self) -> FrameExtractionResult:
        """Run the full frame extraction workflow.

        1. Load prmtop from workspace or file browser
        2. Choose coordinate source (rst7 or trajectory)
        3. Extract frames, producing both rst7 and PDB per frame
        4. Return FrameExtractionResult
        """
        self.console.print(Panel(
            "[bold]AMBER Topology & Coordinates[/bold]\n\n"
            "QM/MM calculations require an AMBER topology (prmtop) for force\n"
            "field parameters and coordinates from either a single rst7 file\n"
            "or frames extracted from an MD trajectory.\n\n"
            "Both rst7 and PDB files will be generated for each frame,\n"
            "allowing either Gaussian (ONIOM) or ORCA workflows downstream.",
            title="Structure Loading",
            border_style="grey50",
            expand=False,
        ))

        # Get prmtop
        prmtop_path = self._get_prmtop()
        self.console.print(f"\n  Prmtop: [cyan]{os.path.basename(prmtop_path)}[/cyan]")

        # Choose coordinate source
        coord_source = prompt_with_context(
            self.processor,
            "Coordinate source: \\[r]st7/inpcrd or \\[t]rajectory",
            default="r",
            module=MODULE_NAME,
            description="Coordinate source type",
            options_map={"r": "rst7/inpcrd file", "t": "Trajectory"},
        ).strip().lower()

        # Set up output directory
        base_dir = self.workspace.get('working_directory', os.getcwd())
        frames_dir = Path(base_dir) / "qmmm_inputs" / "frames"
        frames_dir.mkdir(parents=True, exist_ok=True)

        if coord_source.startswith("t"):
            result = self._load_from_trajectory(prmtop_path, frames_dir)
        else:
            result = self._load_from_rst7(prmtop_path, frames_dir)

        # Display summary
        table = Table(title="Frame Extraction Summary")
        table.add_column("Property", style="cyan")
        table.add_column("Value", style="green")
        table.add_row("Topology", os.path.basename(result.prmtop_path))
        table.add_row("Frames", str(len(result.frame_rst7_paths)))
        for i, label in enumerate(result.frame_labels):
            table.add_row(
                f"  {label}" if label else "  (single)",
                f"rst7: {os.path.basename(result.frame_rst7_paths[i])}, "
                f"pdb: {os.path.basename(result.frame_pdb_paths[i])}"
            )
        self.console.print(table)

        return result

    # ------------------------------------------------------------------
    # Prmtop loading
    # ------------------------------------------------------------------

    def _get_prmtop(self) -> str:
        """Get prmtop path from workspace or file browser.

        Shows workspace prmtop if available, but always allows browsing.
        """
        prmtop_path = self.workspace.get("parm7_file")
        if prmtop_path and os.path.exists(prmtop_path):
            self.console.print(
                f"  Found topology in workspace: [cyan]{os.path.basename(prmtop_path)}[/cyan]"
            )
            choice = prompt_with_context(
                self.processor,
                "Use this topology? (Enter to accept, or 'b' to browse)",
                default="",
                module=MODULE_NAME,
                description="Select AMBER prmtop file",
            ).strip().lower()
            if choice != 'b':
                return prmtop_path

        prmtop_path = self.find_or_browse_file(
            "prmtop", [".prmtop", ".parm7"], "AMBER prmtop file"
        )
        if not os.path.exists(prmtop_path):
            raise FileNotFoundError(f"Prmtop file not found: {prmtop_path}")
        return prmtop_path

    # ------------------------------------------------------------------
    # Single rst7 loading
    # ------------------------------------------------------------------

    def _load_from_rst7(self, prmtop_path: str, output_dir: Path) -> FrameExtractionResult:
        """Load a single rst7 and generate the corresponding PDB."""
        rst7_path = self.workspace.get("rst7_file")
        if rst7_path and os.path.exists(rst7_path):
            self.console.print(
                f"  Found coordinates in workspace: [cyan]{os.path.basename(rst7_path)}[/cyan]"
            )
            choice = prompt_with_context(
                self.processor,
                "Use these coordinates? (Enter to accept, or 'b' to browse)",
                default="",
                module=MODULE_NAME,
                description="Select AMBER rst7/inpcrd file",
            ).strip().lower()
            if choice == 'b':
                rst7_path = None

        if not rst7_path or not os.path.exists(rst7_path):
            rst7_path = self.find_or_browse_file(
                "rst7", [".rst7", ".inpcrd", ".crd"], "AMBER rst7/inpcrd file"
            )
        if not os.path.exists(rst7_path):
            raise FileNotFoundError(f"Rst7 file not found: {rst7_path}")

        self.console.print(f"  Coordinates: [cyan]{os.path.basename(rst7_path)}[/cyan]")

        # Generate PDB from prmtop + rst7
        pdb_path = self._rst7_to_pdb(prmtop_path, rst7_path, output_dir)

        return FrameExtractionResult(
            prmtop_path=prmtop_path,
            frame_rst7_paths=[rst7_path],
            frame_pdb_paths=[pdb_path],
            frame_labels=[""],
        )

    def _rst7_to_pdb(self, prmtop_path: str, rst7_path: str, output_dir: Path) -> str:
        """Convert rst7 to PDB using parmed."""
        import parmed
        struct = parmed.load_file(prmtop_path, rst7_path)
        stem = Path(rst7_path).stem
        pdb_path = str(output_dir / f"{stem}.pdb")
        struct.save(pdb_path, overwrite=True)
        self.console.print(f"  Generated PDB: [cyan]{os.path.basename(pdb_path)}[/cyan]")
        return pdb_path

    # ------------------------------------------------------------------
    # Trajectory loading
    # ------------------------------------------------------------------

    def _load_from_trajectory(self, prmtop_path: str,
                              output_dir: Path) -> FrameExtractionResult:
        """Load coordinates from MD trajectory via cpptraj.

        1. Select trajectory file(s)
        2. Report frame count
        3. Optionally center solute (autoimage)
        4. Select frames (start:end:increment)
        5. Extract selected frames as rst7 + PDB pairs
        6. Let user pick which frames to use
        """
        cpptraj = shutil.which("cpptraj")
        if not cpptraj:
            raise FileNotFoundError(
                "cpptraj not found on PATH. Install AmberTools to use trajectory loading."
            )

        # Select trajectory file(s)
        traj_extensions = [".nc", ".mdcrd", ".crd", ".dcd", ".xtc", ".trr", ".binpos"]
        traj_paths: List[str] = []

        # Offer pre-loaded trajectory files from the workspace (e.g. loaded
        # via the Structure Loader) before falling back to interactive browse.
        ws_trajs = self.workspace.get("trajectory_files") or []
        ws_trajs = [p for p in ws_trajs if os.path.exists(p)]
        if ws_trajs:
            self.console.print(
                f"\n  Found {len(ws_trajs)} trajectory file(s) in workspace:"
            )
            for p in ws_trajs:
                self.console.print(f"    [cyan]{os.path.basename(p)}[/cyan]")
            if confirm_with_context(
                self.processor,
                "Use these trajectory file(s)?",
                default=True,
                module=MODULE_NAME,
                description="Use workspace-loaded trajectory files",
            ):
                traj_paths = list(ws_trajs)

        if not traj_paths:
            self.console.print("\n  [bold]Select trajectory file(s)[/bold]")

        while not traj_paths or confirm_with_context(
            self.processor,
            "Add another trajectory segment?",
            default=False,
            module=MODULE_NAME,
            description="Add another trajectory segment",
        ):
            traj_path = self.find_or_browse_file(
                "trajectory", traj_extensions, "MD trajectory file"
            )
            traj_paths.append(traj_path)
            self.console.print(f"  Added: [cyan]{os.path.basename(traj_path)}[/cyan]")

        # Frame count
        self.console.print("\n  Reading trajectory info...")
        n_frames = self._cpptraj_frame_count(cpptraj, prmtop_path, traj_paths)
        self.console.print(f"  Total frames: [cyan]{n_frames}[/cyan]")

        # Centering
        do_autoimage = confirm_with_context(
            self.processor,
            "Center solute in box (autoimage)?",
            default=True,
            module=MODULE_NAME,
            description="Apply cpptraj autoimage centering",
        )

        # Frame selection
        ex_all = f"1:{n_frames}:1"
        ex_last = f"{n_frames}:{n_frames}:1"
        ex_every = f"1:{n_frames}:100"
        ex_last10 = f"{max(1, n_frames - 9)}:{n_frames}:1"
        w = max(len(ex_all), len(ex_last), len(ex_every), len(ex_last10))

        self.console.print(Panel(
            "[bold]Frame Selection[/bold]\n\n"
            "Specify frames as [cyan]start:end:increment[/cyan] (1-based).\n\n"
            "Examples:\n"
            f"  [cyan]{ex_all:<{w}}[/cyan]  — All frames\n"
            f"  [cyan]{ex_last:<{w}}[/cyan]  — Last frame only\n"
            f"  [cyan]{ex_every:<{w}}[/cyan]  — Every 100th frame\n"
            f"  [cyan]{ex_last10:<{w}}[/cyan]  — Last 10 frames",
            title="Which frames to extract?",
            border_style="grey50",
            expand=False,
        ))

        frame_spec = prompt_with_context(
            self.processor,
            f"Frame selection (start:end:increment, 1-{n_frames})",
            default=f"{n_frames}:{n_frames}:1",
            module=MODULE_NAME,
            description="Trajectory frame selection",
        ).strip()

        parts = frame_spec.split(":")
        try:
            start = int(parts[0]) if len(parts) > 0 else 1
            end = int(parts[1]) if len(parts) > 1 else start
            increment = int(parts[2]) if len(parts) > 2 else 1
        except ValueError:
            raise ValueError(
                f"Invalid frame specification: {frame_spec}. Use start:end:increment"
            )

        if start < 1 or end > n_frames or start > end or increment < 1:
            raise ValueError(
                f"Invalid frame range: {start}:{end}:{increment} "
                f"(trajectory has {n_frames} frames)"
            )

        n_extract = len(range(start, end + 1, increment))
        self.console.print(f"\n  Extracting {n_extract} frame(s): {start}:{end}:{increment}")

        # Extract frames as both rst7 and PDB
        prmtop_stem = Path(prmtop_path).stem
        rst7_files, pdb_files = self._cpptraj_extract_frames(
            cpptraj, prmtop_path, traj_paths,
            start, end, increment,
            do_autoimage, output_dir, prmtop_stem,
        )

        if not rst7_files:
            raise RuntimeError("cpptraj produced no output frames")

        frame_numbers = list(range(start, end + 1, increment))
        self.console.print(
            f"  [green]Extracted {len(rst7_files)} frame(s) to {output_dir}[/green]"
        )

        if len(rst7_files) == 1:
            return FrameExtractionResult(
                prmtop_path=prmtop_path,
                frame_rst7_paths=rst7_files,
                frame_pdb_paths=pdb_files,
                frame_labels=[f"frame{frame_numbers[0]}"],
            )

        # Multiple frames — let user select which to generate inputs for
        self.console.print(
            f"\n  [bold]Select frame(s) for QM/MM input generation:[/bold]\n"
        )
        for i, (rst7, frame_num) in enumerate(zip(rst7_files, frame_numbers), 1):
            self.console.print(f"    [{i}] Frame {frame_num}: {os.path.basename(rst7)}")

        self.console.print(Panel(
            "Select frames using:\n"
            f"  [cyan]1[/cyan]         — Single frame\n"
            f"  [cyan]1,3,5[/cyan]     — Specific frames\n"
            f"  [cyan]1-{len(rst7_files)}[/cyan]       — Range of frames\n"
            f"  [cyan]all[/cyan]       — All {len(rst7_files)} frames\n\n"
            "The QM/MM settings are configured once using the first selected\n"
            "frame. Input files are written for all selected frames.",
            title="Multi-frame mode",
            border_style="grey50",
            expand=False,
        ))

        choice = prompt_with_context(
            self.processor,
            f"Select frame(s) (1-{len(rst7_files)}, or 'all')",
            default="all",
            module=MODULE_NAME,
            description="Select trajectory frames for QM/MM input",
        ).strip()

        selected = self.parse_selection(choice, len(rst7_files))
        if selected is None:
            selected = [1]

        sel_rst7 = [rst7_files[i - 1] for i in selected]
        sel_pdb = [pdb_files[i - 1] for i in selected]
        sel_labels = [f"frame{frame_numbers[i - 1]}" for i in selected]

        if len(sel_rst7) > 1:
            self.console.print(
                f"\n  [green]Selected {len(sel_rst7)} frames — "
                f"input files will be generated for each.[/green]"
            )
        self.console.print(
            f"  Primary frame: [cyan]{os.path.basename(sel_rst7[0])}[/cyan]"
        )

        return FrameExtractionResult(
            prmtop_path=prmtop_path,
            frame_rst7_paths=sel_rst7,
            frame_pdb_paths=sel_pdb,
            frame_labels=sel_labels,
        )

    # ------------------------------------------------------------------
    # cpptraj helpers
    # ------------------------------------------------------------------

    def _run_cpptraj(self, cpptraj: str, cpptraj_input: str,
                     timeout: int = 300) -> subprocess.CompletedProcess:
        """Run cpptraj with input commands via a temp file.

        Using a temp file instead of stdin piping for cross-platform reliability.
        """
        import tempfile

        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.cpptraj', delete=False
        ) as tmp:
            tmp.write(cpptraj_input)
            tmp_path = tmp.name

        try:
            result = subprocess.run(
                [cpptraj, "-i", tmp_path],
                capture_output=True, text=True, timeout=timeout,
            )
        finally:
            os.unlink(tmp_path)

        return result

    def _cpptraj_frame_count(self, cpptraj: str, prmtop: str,
                             traj_paths: List[str]) -> int:
        """Get total frame count from trajectory file(s) using cpptraj."""
        cpptraj_input = f"parm {prmtop}\n"
        for tp in traj_paths:
            cpptraj_input += f"trajin {tp}\n"
        cpptraj_input += "go\n"

        result = self._run_cpptraj(cpptraj, cpptraj_input, timeout=60)

        # Parse "Coordinate processing will occur on X frames."
        # or "Read X frames and processed X frames."
        for line in result.stdout.splitlines():
            match = re.search(
                r'Coordinate processing will occur on (\d+) frames', line)
            if match:
                return int(match.group(1))

        for line in result.stdout.splitlines():
            match = re.search(r'Read (\d+) frames', line)
            if match:
                return int(match.group(1))

        raise RuntimeError(
            f"Could not determine frame count from cpptraj output:\n"
            f"{result.stdout}\n{result.stderr}"
        )

    def _cpptraj_extract_frames(
        self, cpptraj: str, prmtop: str, traj_paths: List[str],
        start: int, end: int, increment: int,
        autoimage: bool, output_dir: Path, stem: str,
    ) -> tuple:
        """Extract trajectory frames as both rst7 and PDB files using cpptraj.

        Returns:
            (rst7_files, pdb_files) — two parallel lists of file path strings.
        """
        cpptraj_input = f"parm {prmtop}\n"
        for tp in traj_paths:
            cpptraj_input += f"trajin {tp}\n"

        if autoimage:
            cpptraj_input += "autoimage\n"

        frame_numbers = list(range(start, end + 1, increment))
        frame_list = ",".join(str(f) for f in frame_numbers)

        if len(frame_numbers) == 1:
            rst7_out = str(output_dir / f"{stem}_frame{frame_numbers[0]}.rst7")
            pdb_out = str(output_dir / f"{stem}_frame{frame_numbers[0]}.pdb")
            cpptraj_input += f"trajout {rst7_out} restart onlyframes {frame_numbers[0]}\n"
            cpptraj_input += f"trajout {pdb_out} pdb onlyframes {frame_numbers[0]}\n"
        else:
            rst7_pattern = str(output_dir / f"{stem}_frame.rst7")
            pdb_pattern = str(output_dir / f"{stem}_frame.pdb")
            cpptraj_input += f"trajout {rst7_pattern} restart onlyframes {frame_list}\n"
            cpptraj_input += f"trajout {pdb_pattern} pdb onlyframes {frame_list}\n"

        cpptraj_input += "go\n"

        self.console.print(f"  [grey50]Running cpptraj...[/grey50]")
        result = self._run_cpptraj(cpptraj, cpptraj_input, timeout=300)

        if result.returncode != 0:
            raise RuntimeError(f"cpptraj failed:\n{result.stderr}")

        rst7_files = sorted(str(f) for f in output_dir.glob(f"{stem}_frame*.rst7"))
        pdb_files = sorted(str(f) for f in output_dir.glob(f"{stem}_frame*.pdb"))

        return rst7_files, pdb_files

    # ------------------------------------------------------------------
    # File browser utilities
    # ------------------------------------------------------------------

    def find_or_browse_file(self, file_type: str, extensions: List[str],
                            description: str) -> str:
        """Find a file by auto-detecting in working directory, or browse.

        1. Scans working directory for files matching the given extensions
        2. If exactly one match, offers it as default
        3. If multiple matches, lets user pick from a numbered list
        4. If no matches, falls back to a directory browser
        """
        work_dir = Path(self.workspace.get('working_directory', os.getcwd()))

        matches = sorted(
            f for f in work_dir.iterdir()
            if f.is_file() and f.suffix in extensions
        )

        if len(matches) == 1:
            self.console.print(
                f"  Auto-detected {file_type}: [cyan]{matches[0].name}[/cyan]"
            )
            if confirm_with_context(
                self.processor,
                f"Use {matches[0].name}?",
                default=True,
                module=MODULE_NAME,
                description=f"Confirm auto-detected {description}",
            ):
                return str(matches[0])

        elif len(matches) > 1:
            self.console.print(
                f"\n  Found {len(matches)} {file_type} file(s) in working directory:\n"
            )
            for i, f in enumerate(matches, 1):
                self.console.print(f"    [{i}] {f.name}")

            match_options_map = {str(i): f.name for i, f in enumerate(matches, 1)}
            match_options_map["b"] = "Browse for a different file"
            choice = prompt_with_context(
                self.processor,
                f"Select {file_type} file (number), or 'b' to browse",
                default="1",
                module=MODULE_NAME,
                description=f"Select {description}",
                options_map=match_options_map,
            ).strip()

            if choice != 'b':
                try:
                    idx = int(choice) - 1
                    if 0 <= idx < len(matches):
                        return str(matches[idx])
                except ValueError:
                    pass

        return self._browse_for_file(file_type, extensions, work_dir)

    def _browse_for_file(self, file_type: str, extensions: List[str],
                         start_dir: Path) -> str:
        """Interactive directory browser for selecting a file."""
        from datetime import datetime
        from proprep.utils.file_browser import file_browser

        def _date_detail(p):
            try:
                return datetime.fromtimestamp(os.path.getmtime(p)).strftime("%m/%d/%Y")
            except OSError:
                return ""

        # Delegate to the shared browser (bare-N, q to cancel, path jump,
        # filename-based session replay). Preserve this method's historical
        # contract: cancelling raises FileNotFoundError for the caller.
        result = file_browser(
            directory=str(start_dir),
            extensions=extensions,
            console=self.console,
            processor=self.processor,
            label=f"{file_type} file",
            entry_detail=_date_detail,
            allow_path_jump=True,
            module=MODULE_NAME,
        )
        if result is None:
            raise FileNotFoundError(f"No {file_type} file selected")
        return result

    # ------------------------------------------------------------------
    # Selection parsing
    # ------------------------------------------------------------------

    @staticmethod
    def parse_selection(selection_str: str, n_items: int) -> Optional[List[int]]:
        """Parse selection string into sorted list of unique 1-based indices.

        Supports: "1", "1,3,5", "1-5", "all"
        """
        selection_str = selection_str.strip().lower()
        if selection_str == "all":
            return list(range(1, n_items + 1))

        indices: set = set()
        for part in selection_str.split(","):
            part = part.strip()
            if "-" in part:
                bounds = part.split("-", 1)
                try:
                    s, e = int(bounds[0].strip()), int(bounds[1].strip())
                except ValueError:
                    return None
                if s > e:
                    return None
                indices.update(range(s, e + 1))
            else:
                try:
                    indices.add(int(part))
                except ValueError:
                    return None

        if not indices or any(i < 1 or i > n_items for i in indices):
            return None
        return sorted(indices)
