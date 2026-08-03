"""
MCPB Force Field Integration Utilities

Pure functions for step 4 (force field integration) of the MCPB workflow:
- Fingerprint parsing for atom types and LINK records
- Unique residue name generation (MCPB convention)
- Mol2 file residue renaming
- PDB residue renaming
- Metadata JSON generation for the FF parameter library
"""

import json
import logging
import shutil
from collections import OrderedDict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

from proprep.utils.tleap_utils import tleap_safe_unit_var

logger = logging.getLogger(__name__)


def parse_fingerprint(fingerprint_path: str,
                      assignments_path: Optional[str] = None) -> Dict[str, Any]:
    """
    Parse a standard.fingerprint file to extract atom types and LINK records.

    The fingerprint has two sections:
    - Atom lines: ``ResID-ResName-AtomName  AtomID  OrigType -> RenamedType``
    - LINK lines: ``LINK AtomID1-AtomName1 AtomID2-AtomName2``

    Args:
        fingerprint_path: Path to standard.fingerprint
        assignments_path: Optional path to atom_type_assignments.json (for element lookup)

    Returns:
        Dictionary with:
        - ``atom_type_entries``: list of addAtomTypes strings like ``'{ "Y1"  "N"  "sp3" }'``
        - ``residues``: OrderedDict of ``{(resid, resname): [{"atom_name", "atom_id", "orig_type", "renamed_type"}, ...]}``
        - ``atom_id_map``: dict ``{atom_id: (resid, resname, atom_name)}``
        - ``link_records``: list of ``{"resid1", "atom1", "resid2", "atom2"}``
    """
    # Load element info from assignments if available
    element_map = {}  # renamed_type -> element
    if assignments_path and Path(assignments_path).exists():
        try:
            with open(assignments_path) as f:
                assignments = json.load(f)
            for _coords, info in assignments.items():
                renamed = info.get("renamed_type", "")
                element = info.get("element", "")
                if renamed and element and renamed not in element_map:
                    element_map[renamed] = element
        except Exception as e:
            logger.warning(f"Could not load atom assignments: {e}")

    residues: Dict[Tuple[int, str], List[Dict]] = OrderedDict()
    atom_id_map: Dict[int, Tuple[int, str, str]] = {}
    seen_renamed_types: Dict[str, str] = {}  # renamed_type -> element

    with open(fingerprint_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("LINK"):
                continue  # Handle LINK lines in second pass

            parts = line.split()
            if len(parts) >= 5 and "->" in line:
                try:
                    # Parse: "ResID-ResName-AtomName  AtomID  OrigType -> RenamedType"
                    label = parts[0]
                    atom_id = int(parts[1])
                    orig_type = parts[2]
                    renamed_type = parts[4]

                    label_parts = label.split("-", 2)
                    if len(label_parts) < 3:
                        continue
                    resid = int(label_parts[0])
                    resname = label_parts[1]
                    atom_name = label_parts[2]

                    key = (resid, resname)
                    if key not in residues:
                        residues[key] = []
                    residues[key].append({
                        "atom_name": atom_name,
                        "atom_id": atom_id,
                        "orig_type": orig_type,
                        "renamed_type": renamed_type,
                    })
                    atom_id_map[atom_id] = (resid, resname, atom_name)

                    # Track renamed types with their elements
                    if renamed_type != orig_type and renamed_type not in seen_renamed_types:
                        element = element_map.get(renamed_type, "")
                        seen_renamed_types[renamed_type] = element

                except (ValueError, IndexError):
                    continue

    # Parse LINK records (second pass)
    link_records = []
    with open(fingerprint_path) as f:
        for line in f:
            line = line.strip()
            if not line or not line.startswith("LINK"):
                continue

            parts = line.split()
            if len(parts) >= 3:
                try:
                    # Format: "LINK AtomID1-AtomName1 AtomID2-AtomName2"
                    id_name1 = parts[1].split("-", 1)
                    id_name2 = parts[2].split("-", 1)
                    atom_id1 = int(id_name1[0])
                    atom_id2 = int(id_name2[0])

                    # Convert atom IDs to residue info
                    if atom_id1 in atom_id_map and atom_id2 in atom_id_map:
                        resid1, resname1, atom_name1 = atom_id_map[atom_id1]
                        resid2, resname2, atom_name2 = atom_id_map[atom_id2]
                        link_records.append({
                            "resid1": resid1,
                            "resname1": resname1,
                            "atom1": atom_name1,
                            "resid2": resid2,
                            "resname2": resname2,
                            "atom2": atom_name2,
                        })
                except (ValueError, IndexError):
                    continue

    # Build addAtomTypes entries for custom (renamed) atom types
    atom_type_entries = []
    for renamed_type, element in sorted(seen_renamed_types.items()):
        if not element:
            # Try to infer element from the type name (M prefix = metal, Y prefix = ligand)
            # This is a fallback; assignments_path should provide accurate elements
            element = _infer_element_from_type(renamed_type)
        # Normalize element symbol to AMBER's case-sensitive convention
        # (first letter upper, remainder lower): a metal stored from a PDB as
        # "MN" must become "Mn", or tleap rejects it ("unknown element MN")
        # and assigns a wrong atomic number / GB radius. Single-letter symbols
        # (O, N, C, H) are unchanged.
        if element:
            element = element[0].upper() + element[1:].lower()
        atom_type_entries.append(f'{{ "{renamed_type}"  "{element}"  "sp3" }}')

    return {
        "atom_type_entries": atom_type_entries,
        "residues": residues,
        "atom_id_map": atom_id_map,
        "link_records": link_records,
    }


def _infer_element_from_type(renamed_type: str) -> str:
    """
    Best-effort element inference from MCPB renamed atom type.

    MCPB uses M-prefix for metals and Y-prefix for ligand atoms,
    but the actual element cannot be reliably determined from the type alone.
    Returns empty string if unknown.
    """
    # This is intentionally conservative — the assignments JSON should be used
    return ""


def generate_unique_residue_names(
    residues: List[Tuple[int, str]],
    existing_names: Optional[Set[str]] = None,
) -> Dict[Tuple[int, str], str]:
    """
    Generate unique 3-character residue names using MCPB convention.

    Convention: first letter + last letter + counter.
    Examples: HID -> HD1, HD2; HIE -> HE1; MNS -> MS1; ZN -> ZN1; CYS -> CS1

    Args:
        residues: List of (resid, resname) tuples
        existing_names: Set of names to avoid (e.g., standard AMBER residues)

    Returns:
        Dict mapping (resid, resname) -> new_name (max 3 chars)
    """
    if existing_names is None:
        existing_names = set()

    # Group residues by base pattern (first+last letter)
    pattern_counters: Dict[str, int] = {}
    result: Dict[Tuple[int, str], str] = {}

    for resid, resname in residues:
        if len(resname) == 0:
            continue

        # Build base pattern: first + last letter
        if len(resname) == 1:
            base = resname[0] + resname[0]
        elif len(resname) == 2:
            base = resname[0] + resname[1]
        else:
            base = resname[0] + resname[-1]

        base = base.upper()

        # Find next available counter
        counter = pattern_counters.get(base, 0) + 1
        while True:
            candidate = f"{base}{counter}"
            if len(candidate) <= 4 and candidate not in existing_names:
                break
            counter += 1
            if counter > 99:
                raise ValueError(
                    f"Cannot generate unique name for {resname} "
                    f"(base '{base}', all counters exhausted)"
                )

        pattern_counters[base] = counter
        result[(resid, resname)] = candidate
        existing_names.add(candidate)

    return result


def rename_mol2_residue(mol2_path: str, new_name: str, output_path: str) -> str:
    """
    Copy a mol2 file with the residue name replaced.

    Updates:
    - @<TRIPOS>MOLECULE section (first data line)
    - @<TRIPOS>ATOM section (residue name column)
    - @<TRIPOS>SUBSTRUCTURE section (substructure name)

    Args:
        mol2_path: Path to source mol2 file
        new_name: New residue name (max 4 chars)
        output_path: Path for output mol2 file

    Returns:
        output_path
    """
    with open(mol2_path) as f:
        lines = f.readlines()

    section = None
    molecule_name_next = False
    output_lines = []

    for line in lines:
        stripped = line.strip()

        # Track section headers
        if stripped.startswith("@<TRIPOS>"):
            section = stripped
            molecule_name_next = (section == "@<TRIPOS>MOLECULE")
            output_lines.append(line)
            continue

        if molecule_name_next:
            # First data line after @<TRIPOS>MOLECULE is the residue name
            output_lines.append(f"{new_name}\n")
            molecule_name_next = False
            continue

        if section == "@<TRIPOS>ATOM" and stripped and not stripped.startswith("@"):
            # Atom line format (from mol2_writer.py):
            # {atom_id:7d} {atom_name:<4s}    {x:10.4f}{y:10.4f}{z:10.4f} {atom_type:<4s} {res_id:6d} {res_name:<4s} {charge:12.6f}
            # Replace the residue name field
            parts = line.split()
            if len(parts) >= 9:
                # Find and replace the residue name (8th field, 0-indexed field 7)
                # The res_name is after res_id (field 7 in 0-indexed, but positions vary)
                # Safest: rebuild the line with the new name
                try:
                    atom_id = int(parts[0])
                    atom_name = parts[1]
                    x = float(parts[2])
                    y = float(parts[3])
                    z = float(parts[4])
                    atom_type = parts[5]
                    res_id = int(parts[6])
                    # parts[7] is old residue name
                    charge = float(parts[8])
                    output_lines.append(
                        f"{atom_id:7d} {atom_name:<4s}    "
                        f"{x:10.4f}{y:10.4f}{z:10.4f} "
                        f"{atom_type:<4s} {res_id:6d} {new_name:<4s} "
                        f"{charge:12.6f}\n"
                    )
                    continue
                except (ValueError, IndexError):
                    pass  # Fall through to write original line

        if section == "@<TRIPOS>SUBSTRUCTURE" and stripped and not stripped.startswith("@"):
            # Format:      1 {residue_name}        1 TEMP              0 ****  ****    0 ROOT
            parts = line.split()
            if len(parts) >= 2:
                parts[1] = new_name
                output_lines.append("     " + "  ".join(parts) + "\n")
                continue

        output_lines.append(line)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        f.writelines(output_lines)

    return output_path


def rename_pdb_residues(
    pdb_path: str,
    residue_map: Dict[Tuple[str, int, str], str],
    output_path: str,
) -> str:
    """
    Rename residues in a PDB file.

    Args:
        pdb_path: Path to input PDB file
        residue_map: Dict of {(chain, resid, old_resname): new_resname}
        output_path: Path for output PDB file

    Returns:
        output_path
    """
    with open(pdb_path) as f:
        lines = f.readlines()

    output_lines = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            # PDB format: cols 18-20 = residue name, col 22 = chain, cols 23-26 = resid
            try:
                old_resname = line[17:20].strip()
                chain = line[21]
                resid = int(line[22:26].strip())

                key = (chain, resid, old_resname)
                if key in residue_map:
                    new_name = residue_map[key]
                    # Pad to 3 chars right-justified in cols 18-20
                    # (standard PDB: right-justified in 3-char field)
                    padded = f"{new_name:>3s}"
                    line = line[:17] + padded + line[20:]
            except (ValueError, IndexError):
                pass

        output_lines.append(line)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        f.writelines(output_lines)

    return output_path


def generate_metadata_json(
    site_type: str,
    description: str,
    atom_type_entries: List[str],
    residue_name_map: Dict[str, str],
    mol2_basenames: List[str],
    frcmod_basename: str,
    redox_state: str = "default",
    spin_state: str = "default",
) -> Dict[str, Any]:
    """
    Generate metadata.json for the FF parameter library entry.

    DEPRECATED / UNUSED: ``create_ff_library`` now writes (and merges) metadata
    via ``proprep.forcefield_params.user_library.promote_state``, which is the
    single source of truth for the on-disk schema. Kept only as a reference for
    the bundled metadata shape; do not wire new callers to it.

    Follows the same schema as existing heme metadata (e.g., bis_his_c_type/metadata.json).

    Args:
        site_type: Site type name (e.g., "zinc_his3_mns")
        description: Human-readable description
        atom_type_entries: List of addAtomTypes strings
        residue_name_map: Dict of {original_label: new_name} for documentation
        mol2_basenames: List of mol2 file basenames (serve as .lib files)
        frcmod_basename: Basename of the frcmod file
        redox_state: Redox state name
        spin_state: Spin state name

    Returns:
        metadata dict ready for JSON serialization
    """
    return {
        "description": description,
        "cofactor_type": site_type,
        "source": "mcpb",
        "references": [],
        "redox_states": {
            redox_state: {
                "description": f"{redox_state} state (MCPB-parameterized)",
                "spin_states": {
                    spin_state: {
                        "description": f"MCPB-generated parameters ({spin_state})",
                        "atom_types": atom_type_entries,
                        "residue_name": residue_name_map,
                        "is_valid": True,
                        "is_default": True,
                        "forcefield_sets": {
                            "mcpb_generated": {
                                "description": f"MCPB-generated parameters for {site_type}",
                                "version": "1.0",
                                "reference": "Generated by ProPrep MCPB workflow",
                                "files": {
                                    "frcmod": frcmod_basename,
                                    "lib": mol2_basenames,
                                },
                                "is_default": True,
                            }
                        },
                    }
                },
            }
        },
    }


def _convert_mol2_to_lib(
    mol2_paths: List[str],
    frcmod_paths: List[str],
    atom_type_entries: List[str],
    output_dir: Path,
) -> List[str]:
    """
    Convert mol2 files to OFF library files via tleap.

    tleap's ``loadmol2`` creates a unit but does NOT register it as a
    residue template for PDB loading.  ``saveoff`` writes an OFF library
    that tleap can later ``loadoff`` to get proper residue matching.

    Returns list of .lib file paths, or empty list on failure.
    """
    import subprocess
    import tempfile

    # Build tleap script
    lines = []

    # Declare custom atom types first
    if atom_type_entries:
        lines.append("addAtomTypes {")
        for entry in atom_type_entries:
            lines.append(f"    {entry}")
        lines.append("}")

    # Load frcmod so tleap knows the parameter types
    for frcmod in frcmod_paths:
        lines.append(f'loadamberparams "{frcmod}"')

    # Load each mol2 and save as OFF library
    lib_paths = []
    for mol2_path in mol2_paths:
        mol2_p = Path(mol2_path)
        unit_name = mol2_p.stem  # e.g., "HD1"
        lib_path = output_dir / f"{unit_name}.lib"
        # tLEaP unit variable must be a valid identifier: a residue code like
        # "9E2" is otherwise lexed as scientific notation. The residue name
        # inside the saved .lib still comes from the mol2, so downstream
        # matching is unaffected. Filenames are quoted for the same reason.
        unit_var = tleap_safe_unit_var(unit_name)
        lines.append(f'{unit_var} = loadmol2 "{mol2_path}"')
        lines.append(f'saveoff {unit_var} "{lib_path}"')
        lib_paths.append(str(lib_path))

    lines.append("quit")

    # Write and run tleap script
    try:
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='_mol2tolib.in', dir=str(output_dir), delete=False
        ) as f:
            f.write('\n'.join(lines))
            script_path = f.name

        result = subprocess.run(
            ["tleap", "-s", "-f", script_path],
            capture_output=True, text=True, cwd=str(output_dir),
            timeout=30
        )

        # Verify all lib files were created
        missing = [p for p in lib_paths if not Path(p).exists()]
        if missing:
            logger.warning(f"mol2→lib conversion: {len(missing)} lib files not created")
            return []

        logger.debug(f"Converted {len(lib_paths)} mol2 files to OFF libraries")
        return lib_paths

    except (subprocess.TimeoutExpired, FileNotFoundError, Exception) as e:
        logger.warning(f"mol2→lib conversion failed: {e}")
        return []


def create_ff_library(
    site_type: str,
    description: str,
    mol2_files: List[str],
    frcmod_files: List[str],
    fingerprint_path: str,
    assignments_path: Optional[str],
    residue_name_map: Dict[Tuple[int, str], str],
    redox_state: str = "default",
    spin_state: str = "default",
    base_dir: Optional[Path] = None,
    atom_type_entries: Optional[List[str]] = None,
    extra_frcmod_files: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Create a complete FF parameter library entry in the user library.

    Deposits (via the merge-safe ``promote_state`` writer) at:
    ``<user library>/metal_sites/{site_type}/{redox_state}/{spin_state}/``

    with renamed mol2/OFF-library files, copied frcmod, and an upserted
    metadata.json leaf. Adding a second redox/spin state for the same site
    merges rather than overwriting; re-running the same state replaces only
    that leaf.

    Args:
        site_type: Site type name
        description: Human-readable description
        mol2_files: List of source mol2 file paths
        frcmod_files: List of source frcmod file paths
        fingerprint_path: Path to standard.fingerprint
        assignments_path: Path to atom_type_assignments.json (optional)
        residue_name_map: Dict of {(resid, resname): new_name}
        redox_state: Redox state name
        spin_state: Spin state name
        base_dir: Base directory for FF params (default: ~/.proprep/forcefield_params/specialized_residues)
        atom_type_entries: Pre-merged addAtomTypes strings. For a multi-site
            protein the caller merges the M*/Y* types across every site's
            fingerprint and passes them here, so the OFF-library conversion and
            metadata cover all sites (not just fingerprint_path's). When None,
            they are parsed from fingerprint_path (single-site behavior).
        extra_frcmod_files: Additional GAFF frcmods that the site's organic
            ligands need (one parmchk2 frcmod per small molecule, e.g.
            ``e4z.frcmod``). Recorded in metadata (files.frcmod becomes a list)
            so the reuse tleap input loads them alongside the bonded frcmod;
            without them the ligand's atoms type but have no vdW/torsion params.

    Returns:
        Dict with:
        - ``library_path``: Path to the library directory
        - ``metadata_path``: Path to metadata.json
        - ``renamed_mol2_files``: list of paths to renamed mol2 files
        - ``frcmod_files``: list of paths to copied frcmod files
    """
    # The MCPB workflow produces parameters for ONE (redox, spin) state at a
    # time, so deposit through the merge-safe user-library writer rather than
    # rewriting metadata.json wholesale (the old path clobbered any other state
    # already recorded for this site). Files are staged, renamed and converted
    # to OFF libraries here; promote_state then copies them into the user
    # library and upserts a single forcefield-set leaf under
    # ``metal_sites/<site_type>/<redox>/<spin>/`` with round-trip validation.
    import tempfile
    from proprep.forcefield_params.user_library import (
        promote_state, PromotionRequest,
    )

    if base_dir is not None:
        logger.debug("create_ff_library: base_dir is ignored; deposits go to "
                     "the loader's user library path.")

    # atom_type_entries is the ONLY thing needed from the fingerprint here
    # (residues come from residue_name_map). For a multi-site protein the caller
    # pre-merges them across all sites; otherwise parse this site's fingerprint.
    if atom_type_entries is None:
        atom_type_entries = parse_fingerprint(fingerprint_path, assignments_path)["atom_type_entries"]

    # Reverse map: old mol2 stem ("HID90") -> new MCPB residue name ("HD1").
    resname_resid_to_new: Dict[str, str] = {
        f"{resname}{resid}": new_name
        for (resid, resname), new_name in residue_name_map.items()
    }

    staging = Path(tempfile.mkdtemp(prefix="proprep_mcpb_promote_"))
    try:
        # Rename mol2 residues into staging.
        renamed_mol2: List[str] = []
        for mol2_path in mol2_files:
            stem = Path(mol2_path).stem  # e.g. "HID90"
            new_name = resname_resid_to_new.get(stem)
            if new_name is None:
                logger.warning(f"No rename mapping for {stem}, copying as-is")
                dest = staging / Path(mol2_path).name
                shutil.copy2(mol2_path, dest)
                renamed_mol2.append(str(dest))
            else:
                out_path = str(staging / f"{new_name}.mol2")
                rename_mol2_residue(mol2_path, new_name, out_path)
                renamed_mol2.append(out_path)
                logger.debug(f"Renamed mol2: {Path(mol2_path).name} -> {new_name}.mol2")

        # Stage frcmod(s): the MCPB bonded frcmod(s) plus each organic ligand's
        # own GAFF frcmod. Dedupe by basename so a ligand frcmod already present
        # in frcmod_files isn't staged twice.
        staged_frcmod = []
        seen_frcmod_names = set()
        for frcmod_path in [*frcmod_files, *(extra_frcmod_files or [])]:
            name = Path(frcmod_path).name
            if name in seen_frcmod_names:
                continue
            seen_frcmod_names.add(name)
            dest = staging / name
            shutil.copy2(frcmod_path, dest)
            staged_frcmod.append(str(dest))
        if not staged_frcmod:
            raise ValueError("create_ff_library: no frcmod files to deposit")

        # Convert mol2 -> OFF library (loadmol2 alone doesn't register a
        # residue template for PDB matching; saveoff does). Fall back to the
        # mol2 files if tleap is unavailable.
        lib_paths = _convert_mol2_to_lib(
            renamed_mol2, staged_frcmod, atom_type_entries, staging
        )
        lib_srcs = lib_paths if lib_paths else renamed_mol2

        # Human-readable residue labels ({"HID90": "HD1", ...}).
        residue_label_map = {
            f"{resname}{resid}": new_name
            for (resid, resname), new_name in residue_name_map.items()
        }

        request = PromotionRequest(
            family="metal_sites",
            type_name=site_type,
            set_name="mcpb_generated",
            frcmod_src=staged_frcmod[0],
            lib_srcs=lib_srcs,
            redox_state=redox_state,
            spin_state=spin_state,
            residue_meta={"description": description, "references": []},
            spin_meta={
                "residue_name": residue_label_map,
                "atom_types": atom_type_entries,
            },
            set_meta={
                "description": description,
                "version": "1.0",
                "reference": "Generated by ProPrep MCPB workflow",
            },
            # A metal site has more than one frcmod: additional per-site bonded
            # frcmods and each organic ligand's own GAFF frcmod. Record them all
            # (files.frcmod becomes a list) so the reuse tleap input emits a
            # loadamberparams for each — otherwise the ligand's atoms type but
            # have no vdW/torsion parameters.
            extra_frcmod_srcs=staged_frcmod[1:],
            # Re-running the same site/redox/spin replaces just that leaf,
            # leaving other states recorded for this site intact.
            on_collision="overwrite",
        )
        result = promote_state(request)
    finally:
        shutil.rmtree(staging, ignore_errors=True)

    copied = result["copied_files"]
    logger.debug(f"Promoted MCPB site '{site_type}' to {result['library_path']}")
    return {
        "library_path": result["library_path"],
        "metadata_path": result["metadata_path"],
        "renamed_mol2_files": [p for p in copied
                               if p.endswith(".lib") or p.endswith(".mol2")],
        "frcmod_files": [p for p in copied if p.endswith(".frcmod")],
        "atom_type_entries": atom_type_entries,
    }
