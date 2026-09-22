"""
Compare the orientation about to be built with the one in the OPM database.

OPM (Orientations of Proteins in Membranes) holds a membrane placement for
deposited entries, computed with PPM and curated: in particular, which side of
the membrane each face of the protein is on comes from the literature, which
no calculation on the structure alone can know. When the structure has a PDB
ID and OPM has the entry, this reports how far the placement in hand is from
OPM's: the angle between the two membrane normals, where the protein sits
along the normal in each, and whether the protein is the same way up.

It is a second opinion, never a source. Nothing is moved, and no answer
("not in OPM", "no network", "this structure is no longer in the deposited
frame") stops a build.

How the two are compared without naming a single residue. ProPrep-prepared
structures are renumbered and their residues renamed (HIE, HIO, CYO, HCO, ...),
so they cannot be paired with OPM's file by identity. But ProPrep never moves
an atom: on 6R2Q all 12,451 heavy atoms are at exactly their deposited
coordinates after filtering, redox-site transformation, protonation naming and
tLEaP hydrogen addition. The file the structure was loaded from is therefore a
bridge. It pairs with OPM's file by identity (chain, residue number, atom
name), which gives the rigid move from the deposited frame to OPM's; and it
shares its frame with the prepared protein, which is checked, not assumed.
That move applied to the prepared protein is "OPM's placement of this
protein", and it is compared, atom for atom in file order, with the placement
in hand.
"""

import os
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from .orientation import OrientationError, _atoms, rigid_transform, transform_from_oriented_copy

OPM_URL = "https://opm-assets.storage.googleapis.com/pdb/{pdb_id}.pdb"

# Three atoms fix a frame. PDB coordinates have three decimals, so atoms of two
# files that coincide in all of them are the same atoms in the same frame, not a
# coincidence.
MIN_COINCIDENT_ATOMS = 3


class OpmCheckError(RuntimeError):
    """OPM could not be asked, or its answer could not be used; the message says why."""


def fetch_opm_entry(pdb_id: str, out_dir: str, timeout: int = 30) -> Optional[str]:
    """
    Download OPM's oriented file for an entry. A file already there is reused.

    Returns:
        the path, or None if OPM has no such entry (soluble proteins, and
        membrane proteins OPM has not processed)

    Raises:
        OpmCheckError: OPM could not be reached
    """
    import requests

    path = os.path.join(out_dir, f"opm_{pdb_id.lower()}.pdb")
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    try:
        response = requests.get(OPM_URL.format(pdb_id=pdb_id.lower()), timeout=timeout)
    except requests.exceptions.RequestException as e:
        # The exception's own text is a paragraph of connection-pool internals.
        raise OpmCheckError(f"OPM could not be reached ({type(e).__name__}; is there an "
                            f"internet connection?)") from e
    if response.status_code in (403, 404):
        return None
    if response.status_code != 200:
        raise OpmCheckError(f"OPM answered HTTP {response.status_code} for {pdb_id}")
    if b"ATOM" not in response.content:
        return None
    with open(path, "wb") as f:
        f.write(response.content)
    return path


def opm_half_thickness(opm_pdb) -> Optional[float]:
    """OPM's "1/2 of bilayer thickness" remark, in A."""
    with open(opm_pdb, errors="replace") as handle:
        for line in handle:
            if line.startswith("REMARK") and "bilayer thickness" in line:
                try:
                    return float(line.split(":")[-1])
                except ValueError:
                    return None
            if line.startswith(("ATOM", "HETATM")):
                return None
    return None


def _coordinates(pdb_path) -> np.ndarray:
    """Every atom's coordinates, in file order."""
    rows = []
    with open(pdb_path, errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                rows.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return np.array(rows)


def _thousandths(xyz) -> tuple:
    # Integers, so that 70.176 read as a float compares equal to itself.
    return tuple(int(round(v * 1000)) for v in xyz)


def coincident_atoms(pdb_a, pdb_b) -> int:
    """Heavy atoms of ``pdb_a`` found at exactly the coordinates of a heavy atom of ``pdb_b``."""
    there = {_thousandths(xyz) for xyz in _atoms(pdb_b).values()}
    return sum(_thousandths(xyz) in there for xyz in _atoms(pdb_a).values())


def find_loaded_entry(workspace, prepared_pdb) -> Optional[tuple]:
    """
    The PDB entry the prepared protein came from: ``(loaded file, PDB ID)``, or None.

    Candidates are the files the Structure Loader recorded. The one whose atoms
    the prepared protein still sits on is taken (by position, as the viewer's
    electron density does, since ProPrep changes names and numbers along the
    way), and its ID is read from its own HEADER record.
    """
    from proprep.structure_prep.density_map import read_pdb_frame

    if workspace is None:
        return None
    paths = [workspace.get("rcsb_pdb_file"), *(workspace.get("rcsb_pdb_files") or []),
             workspace.get("local_pdb_file")]
    best = None
    for path in dict.fromkeys(p for p in paths if isinstance(p, str) and os.path.exists(p)):
        try:
            pdb_id = read_pdb_frame(path)[3]
            shared = coincident_atoms(prepared_pdb, path) if pdb_id else 0
        except (OSError, ValueError):
            continue
        if shared >= MIN_COINCIDENT_ATOMS and (best is None or shared > best[2]):
            best = (path, pdb_id, shared)
    return best[:2] if best else None


def compare_placements(placed: np.ndarray, reference: np.ndarray) -> Dict[str, float]:
    """
    Two placements of the same atoms (same order), each in a frame whose z axis is
    the membrane normal and whose z = 0 plane is the membrane centre.

    Rotation about the normal and sliding in the membrane plane do not change a
    placement, so only these are reported:

    * ``tilt_deg``: angle between the two membrane normals as seen from the protein
    * ``flipped``: the protein is the other way up (the normals point oppositely)
    * ``centre_z`` / ``reference_centre_z``: where the protein's centre sits along
      the normal in each (with ``placed`` turned over first if it is flipped)
    """
    rotation, _, rmsd = rigid_transform(placed, reference)
    # Our normal (0, 0, 1) as a direction in the reference frame is the third row.
    cosine = float(np.clip(rotation[2, 2], -1.0, 1.0))
    angle = float(np.degrees(np.arccos(cosine)))
    flipped = angle > 90.0
    centre_z = float(placed[:, 2].mean())
    return {
        "tilt_deg": 180.0 - angle if flipped else angle,
        "flipped": flipped,
        "centre_z": -centre_z if flipped else centre_z,
        "reference_centre_z": float(reference[:, 2].mean()),
        "fit_rmsd": rmsd,
    }


def cross_check(placed_pdb, unplaced_pdb, loaded_pdb, opm_pdb) -> Dict[str, object]:
    """
    Compare ``placed_pdb`` with OPM's placement of the same protein.

    Args:
        placed_pdb: the prepared protein in the membrane frame about to be used
        unplaced_pdb: the same file before it was moved (same atoms, same order)
        loaded_pdb: the file the structure was originally loaded from
        opm_pdb: OPM's oriented file for the entry

    Raises:
        OpmCheckError: the comparison cannot be made; the message says why
    """
    placed, unplaced = _coordinates(placed_pdb), _coordinates(unplaced_pdb)
    if placed.shape != unplaced.shape or len(placed) < 3:
        raise OpmCheckError(
            f"{Path(placed_pdb).name} and {Path(unplaced_pdb).name} do not hold the same atoms "
            f"({len(placed)} and {len(unplaced)}).")

    same_frame = coincident_atoms(unplaced_pdb, loaded_pdb)
    if same_frame < MIN_COINCIDENT_ATOMS:
        raise OpmCheckError(
            f"no atom of the prepared protein is still at its coordinates in "
            f"{Path(loaded_pdb).name} ({same_frame} found): the structure has been moved since "
            f"it was loaded (aligned, minimised), so OPM's frame cannot be carried over to it.")

    try:
        rotation, translation, _, paired = transform_from_oriented_copy(loaded_pdb, opm_pdb)
    except OrientationError as exc:
        raise OpmCheckError(
            f"OPM's file is not the same model as {Path(loaded_pdb).name} (another assembly, or "
            f"other chain and residue numbers): {exc}") from exc

    result: Dict[str, object] = compare_placements(placed, unplaced @ rotation + translation)
    result.update({"atoms_in_deposited_frame": same_frame, "atoms_paired_with_opm": paired,
                   "opm_half_thickness": opm_half_thickness(opm_pdb)})
    return result


def report_lines(pdb_id: str, result: Dict[str, object]) -> List[str]:
    """The comparison in words and numbers, for the console."""
    if result["flipped"]:
        side = ("    [bold red]the protein is the OTHER WAY UP: the face OPM puts on the +z side of "
                "the membrane is on the -z side here.[/bold red] OPM's sides come from the "
                "literature; a calculation on the structure alone cannot tell them apart. If "
                "the two sides of your system differ (leaflet composition, ions, a potential), "
                "turn the protein over: Protein Orientation, N-terminus orientation.")
    else:
        side = "    same way up as in OPM"
    lines = [
        f"  OPM has {pdb_id.upper()}. Against its placement "
        f"({result['atoms_paired_with_opm']} atoms paired with OPM's file):",
        side,
        f"    membrane normal: {result['tilt_deg']:.1f}° away from OPM's",
        f"    protein centre along the normal: z = {result['centre_z']:+.1f} Å here, "
        f"{result['reference_centre_z']:+.1f} Å in OPM "
        f"(difference {result['centre_z'] - result['reference_centre_z']:+.1f} Å)",
    ]
    if result["opm_half_thickness"] is not None:
        lines.append(f"    OPM's hydrophobic thickness: {2 * result['opm_half_thickness']:.1f} Å")
    return lines
