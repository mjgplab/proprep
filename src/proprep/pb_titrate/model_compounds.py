"""On-demand model compound builder via tleap.

The Bashford 4-state cycle compares protein PB energies against a model
compound — a small molecule with the same titrating chemistry as the
protein-embedded site, dissolved in implicit water. This module builds
those compounds on demand and caches the resulting prmtop+rst7 files.

Two flavours:

* Standard amino-acid residues (AS4, GL4, HIP, LYS, CYS, TYR) → an
  ACE-X-NME tripeptide built under leaprc.constph.

* Heme-related titratables (currently just PRN, the heme propionate)
  → a singleton built under leaprc.conste with a methyl cap added on
  CA via parmed, since PRN lacks peptide head/tail connections and
  has dangling valence on its CA when isolated.

Example
-------
    >>> from titrate.model_compounds import build_model
    >>> prm, rst = build_model("AS4")
    >>> prm
    PosixPath('.../titrate/models/ACE_AS4_NME.prmtop')

    >>> prm, rst = build_model("PRN")
    >>> prm
    PosixPath('.../titrate/models/model_PRN.prmtop')

Cache directory defaults to `<package>/models/` but can be overridden.
"""
from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import parmed as pmd

# Prebuilt model compounds shipped with the package. At runtime (installed into
# a system/conda prefix) this directory is typically read-only, so it is only
# ever read from, never written to.
BUNDLED_MODELS_DIR = Path(__file__).resolve().parent / "models"
# Back-compat alias (older name referenced this).
DEFAULT_CACHE_DIR = BUNDLED_MODELS_DIR


def _user_cache_dir() -> Path:
    """Writable fallback cache for on-demand model builds.

    Used when the package ships no prebuilt compound for a residue and the
    caller passed no explicit ``cache_dir``. Honors ``XDG_CACHE_HOME``. Never
    points inside the (possibly read-only) installed package.
    """
    base = os.environ.get("XDG_CACHE_HOME") or str(Path.home() / ".cache")
    return Path(base) / "proprep" / "pb_titrate_models"


# Heme-related titratables that require leaprc.conste and singleton-
# style model compounds rather than ACE-X-NME tripeptides.
HEME_RESNAMES = {"PRN"}


def _model_paths(resname: str, cache_dir: Path) -> Tuple[Path, Path]:
    """Cache filenames for a given residue's model compound."""
    if resname in HEME_RESNAMES:
        return (cache_dir / f"model_{resname}.prmtop",
                cache_dir / f"model_{resname}.rst7")
    return (cache_dir / f"ACE_{resname}_NME.prmtop",
            cache_dir / f"ACE_{resname}_NME.rst7")


def _build_aa(resname: str, prm: Path, rst: Path) -> None:
    """Build an ACE-RES-NME tripeptide via tleap."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        script = td / "build.tleap"
        script.write_text(
            f"source leaprc.constph\n"
            f"m = sequence {{ ACE {resname} NME }}\n"
            f"saveamberparm m {prm} {rst}\n"
            f"quit\n"
        )
        r = subprocess.run(
            ["tleap", "-f", str(script)],
            capture_output=True, text=True, cwd=str(td))
        if r.returncode != 0 or not prm.exists():
            raise RuntimeError(
                f"tleap failed to build ACE-{resname}-NME (exit {r.returncode}):\n"
                f"stdout: {r.stdout[-2000:]}\nstderr: {r.stderr[-1000:]}")


def _build_heme(resname: str, prm: Path, rst: Path) -> None:
    """Build a singleton heme-titratable residue and methyl-cap its dangling CA.

    For PRN, conste.lib defines the propionate as CA(H)(H)–CB(H)(H)–CG(=O1)(=O2)
    with no peptide head/tail connect atoms; CA is bonded only to CB, HA1, HA2
    in isolation, leaving the fourth valence dangling (it bonds to the
    porphyrin's CB carbon in a real heme). To make a chemically sensible
    standalone model — propanoic acid, with PRN's RESP charges intact — we
    add a single HA3 hydrogen on CA in the missing tetrahedral direction.
    The cap is given charge 0 so the residue's per-state net charges remain
    integer, exactly as cpinutil/--describe lists them.
    """
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        script = td / "build.tleap"
        script.write_text(
            f"source leaprc.constph\n"
            f"source leaprc.conste\n"
            f"m = sequence {{ {resname} }}\n"
            f"saveamberparm m {prm} {rst}\n"
            f"quit\n"
        )
        r = subprocess.run(
            ["tleap", "-f", str(script)],
            capture_output=True, text=True, cwd=str(td))
        if r.returncode != 0 or not prm.exists():
            raise RuntimeError(
                f"tleap failed to build {resname} (exit {r.returncode}):\n"
                f"stdout: {r.stdout[-2000:]}\nstderr: {r.stderr[-1000:]}")

    # Methyl-cap CA: load, add HA3, save back over the cached files.
    s = pmd.load_file(str(prm), str(rst))
    res = s.residues[0]
    ca = next((a for a in res.atoms if a.name == "CA"), None)
    if ca is None:
        raise RuntimeError(f"{resname}: no CA atom found in built singleton")

    # Existing bonded neighbours (CB and the two HA hydrogens).
    bonded = [b.atom1 if b.atom2 is ca else b.atom2
              for b in ca.bonds]
    if len(bonded) != 3:
        raise RuntimeError(
            f"{resname}: expected CA to have 3 existing bonds before capping, "
            f"found {len(bonded)} ({[a.name for a in bonded]})")

    coords = s.coordinates  # Nx3 array
    ca_xyz = np.asarray(coords[ca.idx], dtype=float)
    # HA3 placed opposite the centroid of the existing bonded directions, at
    # the typical aliphatic C–H bond length (1.09 Å). This produces a
    # near-tetrahedral geometry without solving any nonlinear minimisation.
    unit_sum = np.zeros(3)
    for nb in bonded:
        v = np.asarray(coords[nb.idx], dtype=float) - ca_xyz
        n = np.linalg.norm(v)
        if n > 1e-6:
            unit_sum += v / n
    direction = -unit_sum
    norm = np.linalg.norm(direction)
    if norm < 1e-6:
        # Degenerate (shouldn't happen for real PRN geometry); pick +x.
        direction = np.array([1.0, 0.0, 0.0])
        norm = 1.0
    ha3_xyz = ca_xyz + (1.09 / norm) * direction

    # Reuse the bond type from an existing CA-HA1 bond so the new bond has
    # k/req parameters and parmed can write the prmtop.
    ha1 = next((a for a in res.atoms if a.name == "HA1"), None)
    if ha1 is None:
        raise RuntimeError(f"{resname}: no HA1 atom to clone bond type from")
    bt_template = next((b.type for b in ca.bonds
                         if b.atom2 is ha1 or b.atom1 is ha1), None)
    if bt_template is None:
        raise RuntimeError(f"{resname}: no CA-HA1 bond type to clone")

    # Build the new atom. Atomic number 1 (H), HC type matches HA1/HA2,
    # mass 1.008, charge 0.0 to preserve integer per-state net charges.
    ha3 = pmd.Atom(name="HA3", type="HC", charge=0.0,
                   mass=1.008, atomic_number=1)
    s.add_atom(ha3, resname=res.name, resnum=res.number,
               chain=res.chain or "")
    s.bonds.append(pmd.Bond(ca, ha3, type=bt_template))
    # Place coordinates: parmed kept the existing coords; append new row.
    new_coords = np.vstack([coords, ha3_xyz[np.newaxis, :]])
    s.coordinates = new_coords

    s.save(str(prm), overwrite=True)
    s.save(str(rst), overwrite=True)


def build_model(resname: str,
                 cache_dir: Optional[Path] = None,
                 rebuild: bool = False) -> Tuple[Path, Path]:
    """Return (prmtop, rst7) paths for the residue's model compound,
    building if not already available.

    Resolution order:
      1. A prebuilt compound shipped in ``BUNDLED_MODELS_DIR`` (read directly;
         works even when the installed package is read-only). Skipped when
         ``rebuild`` is set.
      2. Otherwise build on demand into ``cache_dir`` — an explicit, writable
         directory (e.g. under the run's ``pb_titrate_work/``). When omitted,
         falls back to a per-user cache (``_user_cache_dir()``).

    Never writes inside the package: a read-only install raised
    ``PermissionError`` on the old ``<package>/models`` default.

    Parameters
    ----------
    resname : str
        Residue type. Standard constph names (AS4, GL4, HIP, LYS, CYS, TYR)
        get an ACE-X-NME tripeptide; heme-related names (PRN) get a singleton
        with a methyl-capped CA.
    cache_dir : Path, optional
        Writable directory for on-demand builds. Default: per-user cache.
    rebuild : bool
        Force a rebuild into ``cache_dir`` even if a prebuilt/cached copy
        exists (ignores the bundled compound).
    """
    # 1. Prefer a prebuilt compound shipped with the package (read-only safe).
    if not rebuild:
        b_prm, b_rst = _model_paths(resname, BUNDLED_MODELS_DIR)
        if b_prm.exists() and b_rst.exists():
            return b_prm, b_rst

    # 2. Build on demand into a writable cache (explicit, else per-user).
    cache_dir = Path(cache_dir) if cache_dir is not None else _user_cache_dir()
    cache_dir.mkdir(parents=True, exist_ok=True)
    prm, rst = _model_paths(resname, cache_dir)
    if rebuild or not prm.exists() or not rst.exists():
        if resname in HEME_RESNAMES:
            _build_heme(resname, prm, rst)
        else:
            _build_aa(resname, prm, rst)
    return prm, rst
