"""Prepare a trajectory for the structure viewer.

NGL (the viewer's engine) reads Amber NetCDF directly, so no format
conversion is needed. One cpptraj run still earns its keep: ``autoimage``
re-centres the solute so a molecule that crossed the periodic box is not
drawn split, and, optionally, ``strip`` drops water and ions. Because the
structure the viewer loads must have exactly the atoms the frames have,
the same run writes both: a first-frame PDB and the NetCDF, from the same
(possibly stripped) topology.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional, Tuple

SOLVENT_MASK = ":WAT,HOH,Na+,Cl-,K+,Cs+,Rb+,Li+,Mg+,Ca+,Zn+"


def _as_list(trajectory) -> List[str]:
    """One trajectory path, or the list of segments the Structure Loader stores."""
    return [str(trajectory)] if isinstance(trajectory, (str, Path)) else [str(t) for t in trajectory]


def cpptraj_input(prmtop: str, trajectory, out_pdb: str, out_nc: str,
                  strip_solvent: bool) -> str:
    """The cpptraj script: image, optionally strip, write PDB + NetCDF.

    ``trajectory`` is a path or a list of segments, read in order as one
    trajectory. cpptraj reads every format the Structure Loader accepts
    (NetCDF, mdcrd, DCD, XTC, TRR, binpos); the viewer always gets NetCDF.
    """
    lines = [f"parm {prmtop}"]
    lines += [f"trajin {segment}" for segment in _as_list(trajectory)]
    lines.append("autoimage")
    if strip_solvent:
        lines.append(f"strip {SOLVENT_MASK}")
    lines += [
        f"trajout {out_pdb} pdb onlyframes 1",
        f"trajout {out_nc} netcdf",
        "run",
        "quit",
    ]
    return "\n".join(lines) + "\n"


def write_view_files(prmtop: str, trajectory, out_dir: str, *,
                     strip_solvent: bool = True,
                     cpptraj: Optional[str] = None) -> Tuple[Path, Path]:
    """Run cpptraj and return ``(pdb, nc)`` for the viewer.

    Outputs go under ``out_dir`` as ``<trajectory stem>_view[_stripped].pdb``
    and ``.nc`` (the first segment's stem, plus ``_Nseg`` when there are
    several). Raises ``RuntimeError`` with cpptraj's output on failure.
    """
    cpptraj = cpptraj or shutil.which("cpptraj")
    if not cpptraj:
        raise RuntimeError("cpptraj not found on PATH (is AmberTools installed?)")
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    segments = _as_list(trajectory)
    if not segments:
        raise RuntimeError("no trajectory file given")
    stem = (Path(segments[0]).stem + (f"_{len(segments)}seg" if len(segments) > 1 else "")
            + ("_view_stripped" if strip_solvent else "_view"))
    pdb, nc = out / f"{stem}.pdb", out / f"{stem}.nc"
    for f in (pdb, nc):
        if f.exists():
            f.unlink()
    script = out / f"{stem}.cpptraj"
    script.write_text(cpptraj_input(str(prmtop), segments, str(pdb), str(nc), strip_solvent))
    result = subprocess.run([cpptraj, "-i", str(script)], capture_output=True, text=True)
    log = out / f"{stem}.log"
    log.write_text(result.stdout + result.stderr)
    if result.returncode != 0 or not pdb.exists() or not nc.exists():
        tail = "\n".join((result.stdout + result.stderr).splitlines()[-15:])
        raise RuntimeError(f"cpptraj failed (exit {result.returncode}); see {log}\n{tail}")
    return pdb, nc


def frame_count(nc_path: str) -> Optional[int]:
    """Number of frames in an Amber NetCDF trajectory, or None if unreadable."""
    try:
        from scipy.io import netcdf_file  # scipy ships with ProPrep
        with netcdf_file(nc_path, "r", mmap=False) as f:
            return int(f.variables["coordinates"].shape[0])
    except Exception:
        return None


def prepare_and_show(processor, console, prmtop, trajectory, out_dir, *, module: str) -> bool:
    """Ask about solvent, write the viewer files with cpptraj, and open them.

    The one path to a playing trajectory, shared by the Structure Viewer (the
    topology and trajectory the Structure Loader put in the workspace) and the
    MD Manager (picked from a simulation directory). ``trajectory`` is a path
    or a list of segments. Returns False if cpptraj failed.
    """
    from proprep.utils.prompts import confirm_with_context
    from proprep.structure_prep.viewer_coordinator import viewer as _viewer

    strip = confirm_with_context(
        processor,
        "Strip water and ions from the viewed trajectory? (no = keep them visible)",
        default=True, module=module,
        description="Strip solvent for viewing")

    segments = _as_list(trajectory)
    what = Path(segments[0]).name if len(segments) == 1 else f"{len(segments)} trajectory segments"
    console.print(f"[grey50]Running cpptraj (autoimage{', strip solvent' if strip else ''}) on {what}...[/grey50]")
    try:
        pdb, nc = write_view_files(str(prmtop), segments, str(out_dir), strip_solvent=strip)
    except RuntimeError as exc:
        console.print(f"[red]{exc}[/red]")
        return False
    n = frame_count(str(nc))
    console.print(f"[green]✓ {pdb.name} + {nc.name}"
                  f"{f' ({n} frames)' if n else ''} written to {pdb.parent}[/green]")
    _viewer.show_trajectory(str(pdb), str(nc), show_waters=not strip, force=True)
    console.print("[grey50]Use the Trajectory panel in the viewer to play or scrub frames, "
                  "and the Movie panel to save an MP4.[/grey50]")
    return True
