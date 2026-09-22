"""
PACKMOL-Memgen Runner

Wraps the packmol-memgen CLI as a subprocess. Tails the log file in real time
to display progress, warnings, and key data (box dimensions, charges, ion
counts) as they become available.

ProPrep never imports packmol-memgen's Python internals — this keeps the
integration decoupled from its unstable internal API.
"""

import logging
import os
import re
import shutil
import signal
import subprocess
import time
import threading
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .packmol_progress import PackmolProgress

logger = logging.getLogger(__name__)


@dataclass
class PackmolResult:
    """Results from a packmol-memgen run."""
    success: bool
    output_pdb: Optional[str] = None
    box_dimensions: Optional[List[float]] = None  # [X, Y, Z]
    lipid_counts: Optional[Dict[str, int]] = None  # {lipid_name: count}
    ion_counts: Optional[Dict[str, int]] = None  # {ion_name: count}
    total_charge: Optional[int] = None
    water_count: Optional[int] = None
    log_file: Optional[str] = None
    error_message: Optional[str] = None
    cli_command: str = ""
    warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "success": self.success,
            "output_pdb": self.output_pdb,
            "box_dimensions": self.box_dimensions,
            "lipid_counts": self.lipid_counts,
            "ion_counts": self.ion_counts,
            "total_charge": self.total_charge,
            "water_count": self.water_count,
            "log_file": self.log_file,
            "error_message": self.error_message,
            "cli_command": self.cli_command,
            "warnings": self.warnings,
        }


def find_packmol_memgen() -> Optional[str]:
    """Locate the packmol-memgen executable."""
    # Check $AMBERHOME/bin first
    amberhome = os.environ.get("AMBERHOME")
    if amberhome:
        candidate = Path(amberhome) / "bin" / "packmol-memgen"
        if candidate.exists():
            return str(candidate)

    # Check PATH
    found = shutil.which("packmol-memgen")
    if found:
        return found

    return None


def run_packmol_memgen(
    args: List[str],
    working_dir: str,
    console=None,
    time_limit_hours: Optional[float] = None,
) -> PackmolResult:
    """
    Run packmol-memgen as a subprocess with real-time log tailing.

    packmol-memgen writes its output to a log file (packmol-memgen.log)
    via Python's logging module, not to stdout, and goes quiet once it has
    started packmol. packmol writes its own log (packmol.log) for the rest of
    the run, which is nearly all of it. Both files are followed in a
    background thread: the first for the setup, the second for one line per
    packing loop (see packmol_progress).

    Args:
        args: CLI arguments (from MembraneConfig.to_cli_args())
        working_dir: Directory to run in (output files go here)
        console: Optional Rich console for progress display
        time_limit_hours: Stop the run after this long. None (the default)
            sets no limit: how long a packing takes depends on the system
            (about 550,000 atoms needed well over an hour), the progress lines
            say how it is going, and Ctrl-C stops it cleanly.

    Returns:
        PackmolResult with parsed output data.
    """
    executable = find_packmol_memgen()
    if executable is None:
        return PackmolResult(
            success=False,
            error_message=(
                "packmol-memgen not found. Ensure AmberTools is installed and "
                "$AMBERHOME/bin is on your PATH."
            ),
        )

    # What the user is shown, and can run themselves.
    cli_command = " ".join([executable] + args)
    # What is run: the same, started so that its logs are written line by line.
    cmd = _live_log_command(executable) + args
    logger.debug(f"Running: {' '.join(cmd)}")

    if console:
        console.print(f"[grey50]Command: {cli_command}[/grey50]\n")

    log_path = Path(working_dir) / "packmol-memgen.log"
    warnings: List[str] = []

    # Remove stale log files so we only tail fresh output. packmol-memgen
    # overwrites both in any case.
    packmol_log_path = Path(working_dir) / "packmol.log"
    for stale in (log_path, packmol_log_path):
        if stale.exists():
            stale.unlink()

    progress = PackmolProgress()
    leaflets = _LeafletVolumes()

    try:
        process = subprocess.Popen(
            cmd,
            cwd=working_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            # Its own process group: packmol is packmol-memgen's child, and
            # stopping the run has to stop both.
            start_new_session=True,
            # packmol is a gfortran program: without this it block-buffers what it
            # writes to packmol-memgen's pipe, before packmol-memgen's own buffer.
            env=dict(os.environ, GFORTRAN_UNBUFFERED_ALL="1"),
        )

        # Tail the log files in a background thread
        stop_tailing = threading.Event()

        def _tail_log():
            """Follow packmol-memgen.log, then packmol.log beside it, and display both."""
            def on_memgen_line(line):
                _process_log_line(line, console, warnings)
                warning = leaflets.feed(line)
                if warning:
                    warnings.append(warning)
                    console.print(f"\n  [bold red]{warning}[/bold red]\n")

            def on_packmol_lines(lines):
                for shown in progress.feed_many(lines):
                    console.print(f"[grey50]{shown}[/grey50]" if shown.startswith("    ")
                                  else f"[cyan]{shown}[/cyan]")

            # [path, open file, bytes of a line still being written, handler]
            # packmol.log first: packmol has finished before packmol-memgen says so.
            followed = [[packmol_log_path, None, b"", on_packmol_lines],
                        [log_path, None, b"", lambda lines: [on_memgen_line(l) for l in lines]]]
            try:
                while True:
                    exited = process.poll() is not None
                    read_any = False
                    for entry in followed:
                        path, handle, pending, on_line = entry
                        if handle is None:
                            if not path.exists():
                                continue
                            handle = entry[1] = open(path, "rb")
                        chunk = handle.read()
                        if not chunk:
                            continue
                        read_any = True
                        *lines, entry[2] = (pending + chunk).split(b"\n")
                        # Handed over together: of several loops that arrive at once,
                        # the latest is the one worth showing.
                        on_line([line.decode(errors="replace").rstrip() for line in lines])
                    if not read_any:
                        if exited or stop_tailing.is_set():
                            return
                        time.sleep(0.1)
            except Exception as e:
                logger.debug(f"Log tailing error: {e}")
            finally:
                for _, handle, _, _ in followed:
                    if handle is not None:
                        handle.close()

        if console:
            tail_thread = threading.Thread(target=_tail_log, daemon=True)
            tail_thread.start()

        # Wait for process to complete
        limit_s = time_limit_hours * 3600 if time_limit_hours else None
        stdout, stderr = process.communicate(timeout=limit_s)

        # Let tail thread finish reading remaining log lines
        stop_tailing.set()
        if console:
            tail_thread.join(timeout=5)

        if process.returncode != 0:
            # Read full log for error extraction
            log_text = ""
            if log_path.exists():
                log_text = log_path.read_text()

            error_msg = _extract_error_message(
                stdout + "\n" + log_text, stderr,
            )
            if not error_msg:
                error_msg = (
                    f"packmol-memgen exited with code {process.returncode}"
                )
                if stderr.strip():
                    error_msg += f":\n{stderr.strip()}"

            return PackmolResult(
                success=False,
                error_message=error_msg,
                cli_command=cli_command,
                warnings=warnings,
                log_file=str(log_path) if log_path.exists() else None,
            )

        # Parse structured data from log file
        parsed = _parse_output(working_dir, args)
        parsed.cli_command = cli_command
        parsed.warnings = warnings
        return parsed

    except (subprocess.TimeoutExpired, KeyboardInterrupt) as stopped:
        _stop_process_group(process)
        stop_tailing.set()
        if isinstance(stopped, KeyboardInterrupt):
            why = "packmol-memgen was stopped (Ctrl-C)"
        else:
            why = f"packmol-memgen was stopped at the time limit of {time_limit_hours:g} h"
        return PackmolResult(
            success=False,
            error_message=(
                f"{why}. {progress.summary()} "
                f"Its full log is {packmol_log_path}."
            ),
            cli_command=cli_command,
            warnings=warnings,
            log_file=str(log_path) if log_path.exists() else None,
        )
    except Exception as e:
        return PackmolResult(
            success=False,
            error_message=f"Error running packmol-memgen: {e}",
            cli_command=cli_command,
            warnings=warnings,
        )


def _live_log_command(executable: str) -> List[str]:
    """How to start packmol-memgen so that packmol.log is written as packmol runs.

    Through packmol_memgen_live_log.py, under the interpreter named on
    packmol-memgen's own shebang line (it may not be ProPrep's). If
    packmol-memgen is not a Python script with such a line, it is started
    directly, as before: the build is the same, the progress display is only
    less frequent.
    """
    launcher = Path(__file__).with_name("packmol_memgen_live_log.py")
    try:
        with open(executable, "rb") as handle:
            first = handle.readline().decode(errors="replace").strip()
    except OSError:
        return [executable]
    if not first.startswith("#!") or "python" not in first or not launcher.exists():
        return [executable]
    interpreter = first[2:].split()
    if not interpreter or (not os.path.basename(interpreter[0]) == "env" and not os.path.exists(interpreter[0])):
        return [executable]
    return interpreter + [str(launcher), executable]


def _stop_process_group(process) -> None:
    """Stop packmol-memgen and the packmol it started."""
    try:
        os.killpg(os.getpgid(process.pid), signal.SIGTERM)
    except (ProcessLookupError, PermissionError, OSError):
        pass
    try:
        process.wait(timeout=10)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(os.getpgid(process.pid), signal.SIGKILL)
        except (ProcessLookupError, PermissionError, OSError):
            process.kill()


class _LeafletVolumes:
    """Notice, as soon as packmol-memgen reports it, a protein that does not cross the membrane.

    packmol-memgen logs the protein's volume inside each leaflet right after
    orienting it. A protein that crosses the bilayer has volume in both. On
    6R2Q (a transmembrane beta-barrel complex) MEMEMBED left 0 A^3 in the upper
    leaflet and 5% of the protein in the lower one, and the build then packed
    for an hour around a protein lying under the membrane. A peripheral protein
    legitimately touches one leaflet only, so this warns and does not stop.
    """

    _VOLUME = re.compile(r"in (upper|lower) leaflet\s*=\s*([\d.]+)")

    def __init__(self):
        self.volume: Dict[str, float] = {}
        self.warned = False

    def feed(self, line: str) -> Optional[str]:
        m = self._VOLUME.search(line)
        if not m or self.warned:
            return None
        self.volume[m.group(1)] = float(m.group(2))
        if len(self.volume) < 2:
            return None
        empty = [side for side, v in self.volume.items() if v == 0.0]
        if len(empty) != 1:
            return None
        self.warned = True
        other = "lower" if empty[0] == "upper" else "upper"
        return (
            f"The oriented protein has no volume in the {empty[0]} leaflet "
            f"({self.volume[other]:.0f} A^3 in the {other} one): as placed, it does not cross "
            f"the membrane. For a transmembrane protein the orientation has failed; stop "
            f"(Ctrl-C) and choose another method under Protein Orientation (PPM3, or a "
            f"pre-oriented structure). For a peripheral protein this is expected."
        )


# ── Real-time log line processing ────────────────────────────────────────

# Warnings to suppress — these are expected in ProPrep's workflow
_SUPPRESSED_WARNINGS: list = [
    # Add patterns here for warnings that are noise in ProPrep's context
]

# Log line timestamp prefix: "MM/DD/YYYY HH:MM:SS AM/PM:LEVEL:"
_LOG_PREFIX = re.compile(
    r"^\d{2}/\d{2}/\d{4}\s+\d{2}:\d{2}:\d{2}\s+[AP]M:(\w+):\s*"
)


def _process_log_line(line: str, console, warnings: List[str]):
    """Parse a packmol-memgen log line and display it appropriately."""
    if not line.strip() or not console:
        return

    # Strip the timestamp prefix, extract log level
    match = _LOG_PREFIX.match(line)
    if match:
        level = match.group(1).upper()
        content = line[match.end():]
    else:
        level = ""
        content = line

    content = content.strip()
    if not content:
        return

    # Check if this is a suppressed warning
    for pattern in _SUPPRESSED_WARNINGS:
        if pattern.search(content):
            return

    # Warnings — show prominently and collect
    if level == "WARNING" or "WARNING" in content.upper():
        warnings.append(content)
        console.print(f"  [yellow]Warning: {content}[/yellow]")
        return

    # Errors
    if level == "ERROR":
        console.print(f"  [red]{content}[/red]")
        return

    # Key informational lines
    _display_info_line(console, content)


def _display_info_line(console, line: str):
    """Display an informational log line if it's interesting."""

    # Preprocessing / orientation steps
    if re.search(r"Preprocessing\s+\S+", line):
        console.print(f"  [cyan]Preprocessing protein...[/cyan]")
        return
    if re.search(r"Orienting the protein using\s+(\S+)", line):
        m = re.search(r"Orienting the protein using\s+(\S+)", line)
        console.print(f"  [cyan]Orienting protein ({m.group(1)})...[/cyan]")
        return

    # Estimated values block
    if re.search(r"Charge\s+=\s+([-\d]+)", line):
        m = re.search(r"Charge\s+=\s+([-\d]+)", line)
        console.print(f"  Protein charge: {m.group(1)}")
        return
    if re.search(r"Mass\s+=\s+([\d.]+)", line):
        m = re.search(r"Mass\s+=\s+([\d.]+)", line)
        console.print(f"  Protein mass: {float(m.group(1)):.0f} Da")
        return
    if re.search(r"Estimated volume\s+=\s+([\d.]+)", line):
        m = re.search(r"Estimated volume\s+=\s+([\d.]+)", line)
        console.print(f"  Estimated protein volume: {float(m.group(1)):.0f} A^3")
        return

    # Box dimensions
    box_match = re.search(r"(x|y|z)_len\s+=\s+([-\d.]+)", line)
    if box_match:
        axis = box_match.group(1).upper()
        val = float(box_match.group(2))
        console.print(f"  Box {axis}: {val:.1f} A")
        return

    # Lipid charge contribution
    if "lipids contribute a charge" in line.lower():
        m = re.search(r"charge of\s+([-\d]+)", line)
        if m:
            console.print(f"  Lipid charge contribution: {m.group(1)}")
        return

    # Ion concentrations
    if re.search(r"(Positive|Negative) ion concentration:\s+([\d.]+)", line):
        m = re.search(r"(Positive|Negative) ion concentration:\s+([\d.]+)", line)
        console.print(f"  {m.group(1)} ion concentration: {m.group(2)} M")
        return
    if re.search(r"Salt concentration specified:\s+([\d.]+)", line):
        m = re.search(r"Salt concentration specified:\s+([\d.]+)", line)
        console.print(f"  Salt concentration: {m.group(1)} M")
        return

    # Lipid counts in leaflets
    if re.search(r"\d+\s+\w+\s+lipids?\s+in\s+(upper|lower)", line, re.IGNORECASE):
        console.print(f"  [grey50]{line}[/grey50]")
        return

    # Packing progress
    if re.search(r"Information for packing", line, re.IGNORECASE):
        console.print(f"\n  [cyan]Packing...[/cyan]")
        return
    if "running packmol" in line.lower() or "packing" in line.lower():
        console.print(f"  [cyan]{line}[/cyan]")
        return

    # PACKMOL convergence / success
    if "success" in line.lower() and "pack" in line.lower():
        console.print(f"  [green]{line}[/green]")
        return
    if "GENCAN" in line or "function value" in line.lower():
        # PACKMOL optimization lines — show sparingly
        console.print(f"  [grey50]{line}[/grey50]")
        return

    # Output file written
    if re.search(r"\.pdb\s+written", line, re.IGNORECASE):
        console.print(f"  [green]{line}[/green]")
        return

    # Water model selection
    if "water model" in line.lower() and "using" in line.lower():
        console.print(f"  {line}")
        return

    # Input/Output PDB
    if re.search(r"(Input|Output)\s+PDB\s+=", line):
        console.print(f"  [grey50]{line}[/grey50]")
        return

    # Experimental value fallback
    if "experimental value" in line.lower() and "not available" in line.lower():
        console.print(f"  [grey50]{line}[/grey50]")
        return


# ── Error extraction ─────────────────────────────────────────────────────

def _extract_error_message(combined: str, stderr: str) -> Optional[str]:
    """Extract a user-friendly error message from output."""
    text = combined + "\n" + stderr

    # Salt concentration too low
    if "concentration of ions required to neutralize" in text.lower():
        return (
            "The system charge requires more ions than the specified salt "
            "concentration provides. ProPrep should pass --salt_override "
            "automatically -- if you see this, please report it as a bug."
        )

    # PACKMOL convergence failure
    if "could not find" in text.lower() and "packmol" in text.lower():
        return "PACKMOL failed to converge. Try increasing nloop values."

    # MEMEMBED / PPM3 failure
    if "memembed" in text.lower() and "error" in text.lower():
        return "Protein orientation failed (MEMEMBED error). Check the log."

    if "ppm" in text.lower() and ("error" in text.lower() or "fail" in text.lower()):
        return "Protein orientation failed (PPM3 error). Check the log."

    return None


# ── Post-run parsing ─────────────────────────────────────────────────────

def _parse_output(
    working_dir: str,
    args: List[str],
) -> PackmolResult:
    """Parse packmol-memgen output to extract structured data."""

    # Determine output PDB name from args
    output_name = "bilayer.pdb"
    for i, arg in enumerate(args):
        if arg in ("-o", "--output") and i + 1 < len(args):
            output_name = args[i + 1]
            break

    output_pdb = Path(working_dir) / output_name
    if not output_pdb.exists():
        return PackmolResult(
            success=False,
            error_message=f"Expected output file not found: {output_pdb}",
        )

    # Read the log file for structured parsing
    log_path = Path(working_dir) / "packmol-memgen.log"
    log_text = ""
    if log_path.exists():
        try:
            log_text = log_path.read_text()
        except Exception:
            pass

    # Parse box dimensions
    box_dims = _parse_box_dimensions(log_text)

    # Parse protein charge
    total_charge = _parse_charge(log_text)

    # Parse lipid counts from log
    lipid_counts = _parse_lipid_counts(log_text)

    # Count ions and water molecules in the output PDB
    ion_counts, water_count = _count_species_in_pdb(str(output_pdb))

    log_file = str(log_path) if log_path.exists() else None

    return PackmolResult(
        success=True,
        output_pdb=str(output_pdb),
        box_dimensions=box_dims,
        lipid_counts=lipid_counts,
        ion_counts=ion_counts,
        total_charge=total_charge,
        water_count=water_count,
        log_file=log_file,
    )


def _parse_box_dimensions(log_text: str) -> Optional[List[float]]:
    """Extract box dimensions from packmol-memgen log."""
    dims = {}
    for match in re.finditer(r"(x|y|z)_len\s+=\s+([-\d.]+)", log_text):
        dims[match.group(1)] = float(match.group(2))

    if len(dims) == 3:
        return [dims["x"], dims["y"], dims["z"]]

    # Fallback patterns
    patterns = [
        r"Box\s+dimensions?\s*[:\s]+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)",
        r"box\s*=\s*\{?\s*([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)",
        r"setBox.*\{\s*([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s*\}",
    ]
    for pattern in patterns:
        match = re.search(pattern, log_text, re.IGNORECASE)
        if match:
            return [float(match.group(1)), float(match.group(2)), float(match.group(3))]

    return None


def _parse_charge(log_text: str) -> Optional[int]:
    """Extract total protein charge from log."""
    match = re.search(r"Charge\s+=\s+([-\d]+)", log_text)
    if match:
        return int(match.group(1))
    return None


def _parse_lipid_counts(log_text: str) -> Optional[Dict[str, int]]:
    """Extract lipid counts per type from log."""
    counts: Dict[str, int] = {}

    for match in re.finditer(
        r"(\d+)\s+(\w+)\s+lipids?\s+in\s+(upper|lower)",
        log_text,
        re.IGNORECASE,
    ):
        count = int(match.group(1))
        lipid = match.group(2)
        counts[lipid] = counts.get(lipid, 0) + count

    return counts if counts else None


# Common ion residue names in AMBER
ION_RESIDUE_NAMES = {
    "K+", "K", "Na+", "NA", "Li+", "LI", "Rb+", "RB", "Cs+", "CS",
    "Cl-", "CL", "Br-", "BR", "I-", "IOD", "F-",
    "Mg2+", "MG", "Ca2+", "CA", "Zn2+", "ZN", "Fe2+", "FE2",
    "Ba2+", "BA", "Cu2+", "CU",
    "NH4", "NH4+",
}


def _count_species_in_pdb(pdb_path: str) -> Tuple[Dict[str, int], int]:
    """
    Count ions and water molecules in a PDB file.

    Returns:
        Tuple of (ion_counts dict, water_molecule_count).
    """
    ion_counts: Dict[str, int] = {}
    water_count = 0
    seen_water_residues = set()
    seen_ion_atoms = set()

    try:
        with open(pdb_path, "r") as f:
            for line in f:
                if not line.startswith(("ATOM", "HETATM")):
                    continue

                res_name = line[17:20].strip()
                chain_id = line[21]
                res_num = line[22:26].strip()
                residue_key = f"{chain_id}:{res_num}:{res_name}"

                if res_name in ("WAT", "HOH", "TIP3", "TP3", "TP4", "SPC", "OPC"):
                    if residue_key not in seen_water_residues:
                        seen_water_residues.add(residue_key)
                        water_count += 1
                elif res_name in ION_RESIDUE_NAMES:
                    if residue_key not in seen_ion_atoms:
                        seen_ion_atoms.add(residue_key)
                        ion_counts[res_name] = ion_counts.get(res_name, 0) + 1

    except Exception as e:
        logger.warning(f"Error counting species in PDB: {e}")

    return ion_counts, water_count
