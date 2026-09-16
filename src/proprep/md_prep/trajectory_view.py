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
from typing import Optional, Tuple

SOLVENT_MASK = ":WAT,HOH,Na+,Cl-,K+,Cs+,Rb+,Li+,Mg+,Ca+,Zn+"


def cpptraj_input(prmtop: str, trajectory: str, out_pdb: str, out_nc: str,
                  strip_solvent: bool) -> str:
    """The cpptraj script: image, optionally strip, write PDB + NetCDF."""
    lines = [
        f"parm {prmtop}",
        f"trajin {trajectory}",
        "autoimage",
    ]
    if strip_solvent:
        lines.append(f"strip {SOLVENT_MASK}")
    lines += [
        f"trajout {out_pdb} pdb onlyframes 1",
        f"trajout {out_nc} netcdf",
        "run",
        "quit",
    ]
    return "\n".join(lines) + "\n"


def write_view_files(prmtop: str, trajectory: str, out_dir: str, *,
                     strip_solvent: bool = True,
                     cpptraj: Optional[str] = None) -> Tuple[Path, Path]:
    """Run cpptraj and return ``(pdb, nc)`` for the viewer.

    Outputs go under ``out_dir`` as ``<trajectory stem>_view[_stripped].pdb``
    and ``.nc``. Raises ``RuntimeError`` with cpptraj's output on failure.
    """
    cpptraj = cpptraj or shutil.which("cpptraj")
    if not cpptraj:
        raise RuntimeError("cpptraj not found on PATH (is AmberTools installed?)")
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    stem = Path(trajectory).stem + ("_view_stripped" if strip_solvent else "_view")
    pdb, nc = out / f"{stem}.pdb", out / f"{stem}.nc"
    for f in (pdb, nc):
        if f.exists():
            f.unlink()
    script = out / f"{stem}.cpptraj"
    script.write_text(cpptraj_input(str(prmtop), str(trajectory), str(pdb), str(nc), strip_solvent))
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
