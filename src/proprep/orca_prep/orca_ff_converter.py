"""
ORCA Force Field Converter

Wraps the `orca_mm -convff -AMBER` command to convert AMBER prmtop files
into ORCA's native ORCAFF (.prms) format.
"""

import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional, Tuple


def find_orca_mm(orca_dir: Optional[str] = None) -> Optional[str]:
    """
    Locate the orca_mm executable.

    Searches in order:
    1. User-provided ORCA directory
    2. PATH
    3. Common installation directories

    Args:
        orca_dir: Optional path to ORCA installation directory

    Returns:
        Full path to orca_mm executable, or None if not found
    """
    # Check user-provided directory
    if orca_dir:
        candidate = os.path.join(orca_dir, "orca_mm")
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate

    # Check PATH
    found = shutil.which("orca_mm")
    if found:
        return found

    return None


def convert_prmtop_to_orcaff(
    prmtop_path: str,
    orca_mm_path: str,
    output_dir: Optional[str] = None,
    verbose: bool = False,
) -> Tuple[bool, str, str]:
    """
    Run orca_mm -convff -AMBER to convert a prmtop to ORCA force field format.

    The command: orca_mm -convff -AMBER <prmtop_path>
    produces a .prms file in the same directory as the prmtop (or output_dir).

    Args:
        prmtop_path: Path to the AMBER prmtop file
        orca_mm_path: Path to the orca_mm executable
        output_dir: Optional directory for output (defaults to prmtop directory)
        verbose: If True, add -verbose flag

    Returns:
        Tuple of (success, output_prms_path, command_output)
    """
    prmtop = Path(prmtop_path).resolve()
    if not prmtop.exists():
        return False, "", f"Prmtop file not found: {prmtop}"

    # Build command
    cmd = [orca_mm_path, "-convff"]
    if verbose:
        cmd.append("-verbose")
    cmd.extend(["-AMBER", str(prmtop)])

    # Run in the output directory so the .prms file lands there
    work_dir = output_dir if output_dir else str(prmtop.parent)

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=work_dir,
            timeout=120,
        )

        combined_output = ""
        if result.stdout:
            combined_output += result.stdout
        if result.stderr:
            combined_output += result.stderr

        if result.returncode != 0:
            return False, "", f"orca_mm failed (exit code {result.returncode}):\n{combined_output}"

        # orca_mm writes the .prms file next to the input prmtop, not in cwd
        prmtop_dir = prmtop.parent
        search_dirs = [Path(work_dir), prmtop_dir] if str(prmtop_dir) != work_dir else [Path(work_dir)]

        prms_path = None
        candidate_names = [
            prmtop.stem + ".ORCAFF.prms",
            prmtop.stem + ".prms",
            "ORCAFF.prms",
        ]

        for search_dir in search_dirs:
            for candidate_name in candidate_names:
                candidate = search_dir / candidate_name
                if candidate.exists():
                    prms_path = candidate
                    break
            if prms_path:
                break

        if not prms_path:
            # Check for any .prms files in both directories
            for search_dir in search_dirs:
                prms_files = list(search_dir.glob("*.prms"))
                if prms_files:
                    prms_path = max(prms_files, key=lambda p: p.stat().st_mtime)
                    break

        if not prms_path:
            return False, "", f"orca_mm ran successfully but no .prms file was found.\nSearched: {', '.join(str(d) for d in search_dirs)}\nCommand output:\n{combined_output}"

        # Move to output_dir if it was found elsewhere
        if output_dir and prms_path.parent != Path(output_dir):
            dest = Path(output_dir) / prms_path.name
            shutil.move(str(prms_path), str(dest))
            prms_path = dest

        return True, str(prms_path), combined_output

    except subprocess.TimeoutExpired:
        return False, "", "orca_mm timed out after 120 seconds"
    except FileNotFoundError:
        return False, "", f"orca_mm executable not found at: {orca_mm_path}"
    except Exception as e:
        return False, "", f"Error running orca_mm: {e}"
