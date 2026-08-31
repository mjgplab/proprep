"""
User Forcefield Library Writer

The write-side counterpart to :mod:`proprep.forcefield_params.loader`. It
*promotes* freshly-parameterized residues — small molecules, modified amino
acids, and metal sites — into the user library at
``~/.proprep/forcefield_params/specialized_residues/`` so that they extend the
bundled library and become discoverable by the loader (which already walks the
user tree; see ``get_available_cofactor_types``).

Design
------
Each ProPrep parameterizer produces parameters for exactly ONE state at a time
(one redox / protonation / spin state). So the core operation here is an
**upsert of a single leaf**:

    <family>/<type>/<redox_state>/<spin_state>/<set_name>.{frcmod,lib}

into a residue's ``metadata.json`` (creating the file and residue when new,
splicing the leaf in when the residue already exists), guarded by a collision
check on the exact ``(redox_state, spin_state, set_name)`` triple.

Every promotion is **round-tripped**: after writing, the entry is re-read
through the loader API; if it does not come back cleanly the whole operation is
rolled back (copied files removed, ``metadata.json`` restored), so a failed
promotion never leaves the user library in a half-written state.

The metadata schema written here is identical to the bundled schema (see any
``specialized_residues/*/*/metadata.json``), so no loader changes are needed.
Non-redox residues (small molecules, modified AAs) use the degenerate
single-state shape ``single_state/default`` exactly as the bundled
``zinc/cys4`` entry does.
"""

import filecmp
import json
import logging
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from proprep.forcefield_params import loader

logger = logging.getLogger(__name__)

# Default state coordinates for non-redox residues (small molecules, modified
# AAs). Mirrors the bundled zinc/cys4 single-state entry.
DEFAULT_REDOX_STATE = "single_state"
DEFAULT_SPIN_STATE = "default"

# How to handle an existing (redox_state, spin_state, set_name) leaf.
ON_COLLISION_ERROR = "error"          # refuse (default — caller decides)
ON_COLLISION_OVERWRITE = "overwrite"  # replace files + metadata leaf
ON_COLLISION_VERSION = "version_bump"  # write under <set_name>_vN instead


class UserLibraryError(Exception):
    """Base error for user-library promotion failures."""


class LibraryCollisionError(UserLibraryError):
    """Raised when a state/set already exists and on_collision='error'."""


@dataclass
class PromotionRequest:
    """A request to deposit one parameterized state into the user library.

    The four ``*_meta`` dicts map directly onto the metadata schema levels, so
    the writer stays dumb (it only splices and validates) and callers — chiefly
    the interactive wizard — own all human-facing content.

    Levels:
      * ``residue_meta``  → top of metadata.json: ``description``,
        ``cofactor_type`` (defaults to ``type_name``), ``references``,
        ``methodology``, ``set_comparison_guidance``, ``prerequisites``.
      * ``state_meta``    → the redox-state block: ``description``,
        ``formal_charge``.
      * ``spin_meta``     → the spin-state block: ``description``,
        ``residue_name`` (str for a single residue, dict for multi-residue
        metal sites), ``atom_types`` (metal sites only; ``[]`` otherwise),
        ``spin_multiplicity``, ``electronic_config``, plus any ligand-residue
        keys. ``is_default`` defaults True.
      * ``set_meta``      → the forcefield-set block: ``description``,
        ``version``, ``reference``, ``name``, ``methodology``,
        ``ph_treatment``, ``protonation_model``, ``prerequisites``.
        ``is_default`` defaults True.
    """

    family: str                       # e.g. "heme", "small_molecules", "modified_aa"
    type_name: str                    # e.g. "bis_his_c_type", "LIG", "SEP"
    set_name: str                     # e.g. "Guberman_LIG_RESP"
    frcmod_src: str                   # source path to the .frcmod
    lib_srcs: List[str]               # source paths to .lib (and/or .mol2)

    redox_state: str = DEFAULT_REDOX_STATE
    spin_state: str = DEFAULT_SPIN_STATE

    residue_meta: Dict[str, Any] = field(default_factory=dict)
    state_meta: Dict[str, Any] = field(default_factory=dict)
    spin_meta: Dict[str, Any] = field(default_factory=dict)
    set_meta: Dict[str, Any] = field(default_factory=dict)

    extra_file_srcs: List[str] = field(default_factory=list)  # copied verbatim
    # Additional .frcmod sources that must ALSO be recorded in metadata
    # (files.frcmod becomes a list) so the loader/generator emit a
    # `loadamberparams` for each. A metal site has >1 frcmod: the MCPB bonded
    # frcmod plus each organic ligand's own GAFF frcmod. Unlike
    # ``extra_file_srcs`` (copied but not referenced), these are load-bearing.
    extra_frcmod_srcs: List[str] = field(default_factory=list)
    on_collision: str = ON_COLLISION_ERROR


def promote_state(request: PromotionRequest) -> Dict[str, Any]:
    """Deposit one parameterized state into the user library.

    Returns a dict with ``library_path`` (the residue dir), ``metadata_path``,
    ``state_dir``, ``set_name`` (the actually-written name — may differ from
    the requested one under ``version_bump``), and ``copied_files``.

    Raises:
        FileNotFoundError: a source file is missing.
        LibraryCollisionError: the leaf exists and ``on_collision='error'``.
        UserLibraryError: round-trip validation failed (already rolled back).
    """
    _verify_sources_exist(request)

    base = loader.get_user_forcefield_base_path()
    cofactor_path = f"{request.family}/{request.type_name}"
    residue_dir = base / request.family / request.type_name
    metadata_path = residue_dir / "metadata.json"

    # --- read-modify-write: load existing metadata or start a skeleton -------
    original_bytes = metadata_path.read_bytes() if metadata_path.exists() else None
    metadata = json.loads(original_bytes) if original_bytes else _new_metadata_skeleton(request)

    # --- collision resolution on the exact (redox, spin, set) leaf -----------
    set_name = _resolve_set_name(metadata, request)

    state_dir = residue_dir / request.redox_state / request.spin_state

    # Track everything we create so we can undo it on validation failure.
    created_files: List[Path] = []
    created_dirs: List[Path] = []
    replaced: List[tuple] = []
    # Names another set in this same state directory already points at, grown
    # as we go so two sources in THIS deposit cannot collide with each other
    # either.
    reserved = _reserved_basenames(metadata, request, set_name)

    def _copy(src: str) -> str:
        name = _copy_into(src, state_dir, created_files, replaced, reserved)
        reserved.add(name)
        return name

    try:
        _ensure_dirs(state_dir, created_dirs)

        frcmod_basename = _copy(request.frcmod_src)
        extra_frcmod_basenames = [_copy(src) for src in request.extra_frcmod_srcs]
        lib_basenames = [_copy(src) for src in request.lib_srcs]
        for src in request.extra_file_srcs:
            _copy(src)

        # lib metadata field is a bare string for a single file, a list for
        # the MCPB multi-mol2 case (matches the loader's lib_ref handling).
        lib_value: Union[str, List[str]] = (
            lib_basenames[0] if len(lib_basenames) == 1 else lib_basenames
        )
        # frcmod follows the same str-or-list convention: a lone bonded frcmod
        # stays a string; a metal site with ligand GAFF frcmods becomes a list
        # so each gets its own loadamberparams line.
        frcmod_value: Union[str, List[str]] = (
            frcmod_basename if not extra_frcmod_basenames
            else [frcmod_basename, *extra_frcmod_basenames]
        )

        _splice_leaf(metadata, request, set_name, frcmod_value, lib_value)

        metadata_path.write_text(json.dumps(metadata, indent=2))

        _validate_via_loader(cofactor_path, request.redox_state,
                             request.spin_state, set_name)
    except Exception:
        _rollback(metadata_path, original_bytes, created_files, created_dirs,
                  replaced)
        raise

    logger.debug("Promoted %s/%s [%s/%s/%s] into user library at %s",
                request.family, request.type_name, request.redox_state,
                request.spin_state, set_name, residue_dir)

    return {
        "library_path": str(residue_dir),
        "metadata_path": str(metadata_path),
        "state_dir": str(state_dir),
        "set_name": set_name,
        "copied_files": [str(p) for p in created_files],
    }


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #

def _verify_sources_exist(request: PromotionRequest) -> None:
    missing = [p for p in [request.frcmod_src, *request.lib_srcs,
                           *request.extra_frcmod_srcs, *request.extra_file_srcs]
               if not Path(p).is_file()]
    if missing:
        raise FileNotFoundError(
            "Cannot promote — source file(s) not found: " + ", ".join(missing)
        )


def _new_metadata_skeleton(request: PromotionRequest) -> Dict[str, Any]:
    """Top-level metadata for a brand-new residue (no redox states yet)."""
    rm = request.residue_meta
    skeleton: Dict[str, Any] = {
        "description": rm.get("description", f"{request.type_name} (user-parameterized)"),
        "cofactor_type": rm.get("cofactor_type", request.type_name),
        "source": "user",
        "references": list(rm.get("references", [])),
        "redox_states": {},
    }
    # Carry through optional top-level fields only when supplied.
    for key in ("methodology", "set_comparison_guidance", "prerequisites"):
        if key in rm:
            skeleton[key] = rm[key]
    return skeleton


def _resolve_set_name(metadata: Dict[str, Any], request: PromotionRequest) -> str:
    """Apply the collision policy; return the set name to actually write."""
    existing = (metadata.get("redox_states", {})
                .get(request.redox_state, {})
                .get("spin_states", {})
                .get(request.spin_state, {})
                .get("forcefield_sets", {}))

    if request.set_name not in existing:
        return request.set_name

    if request.on_collision == ON_COLLISION_OVERWRITE:
        return request.set_name
    if request.on_collision == ON_COLLISION_VERSION:
        n = 2
        while f"{request.set_name}_v{n}" in existing:
            n += 1
        return f"{request.set_name}_v{n}"
    raise LibraryCollisionError(
        f"Set '{request.set_name}' already exists for "
        f"{request.family}/{request.type_name} "
        f"[{request.redox_state}/{request.spin_state}]. "
        f"Choose overwrite or version_bump."
    )


def _splice_leaf(metadata: Dict[str, Any], request: PromotionRequest,
                 set_name: str, frcmod_value: Union[str, List[str]],
                 lib_value: Union[str, List[str]]) -> None:
    """Insert/replace the redox→spin→set leaf, preserving sibling entries."""
    redox_states = metadata.setdefault("redox_states", {})

    redox_block = redox_states.setdefault(request.redox_state, {})
    redox_block.setdefault(
        "description",
        request.state_meta.get("description", f"{request.redox_state} state"),
    )
    if "formal_charge" in request.state_meta:
        redox_block["formal_charge"] = request.state_meta["formal_charge"]
    spin_states = redox_block.setdefault("spin_states", {})

    spin_block = spin_states.setdefault(request.spin_state, {})
    sm = request.spin_meta
    spin_block.setdefault(
        "description",
        sm.get("description", f"{request.spin_state} ({request.set_name})"),
    )
    # residue_name / atom_types belong to the spin level. Set them on first
    # write; leave existing values alone so a second set added to the same
    # state doesn't silently change them.
    if "residue_name" in sm:
        spin_block.setdefault("residue_name", sm["residue_name"])
    spin_block.setdefault("atom_types", list(sm.get("atom_types", [])))
    for key in ("spin_multiplicity", "electronic_config"):
        if key in sm:
            spin_block.setdefault(key, sm[key])
    # Pass through any ligand-residue keys (e.g. ligand_cys_residues).
    for key, val in sm.items():
        if key.startswith("ligand_"):
            spin_block.setdefault(key, val)
    spin_block.setdefault("is_valid", True)
    spin_block.setdefault("is_default", True)

    ff_sets = spin_block.setdefault("forcefield_sets", {})
    set_block: Dict[str, Any] = {
        "description": request.set_meta.get("description", f"User-parameterized {set_name}"),
        "version": request.set_meta.get("version", "1.0"),
        "reference": request.set_meta.get("reference", ""),
        "files": {"frcmod": frcmod_value, "lib": lib_value},
        "is_default": request.set_meta.get("is_default", True),
    }
    # Optional rich per-set fields, written only when provided.
    for key in ("name", "methodology", "ph_treatment", "protonation_model",
                "prerequisites"):
        if key in request.set_meta:
            set_block[key] = request.set_meta[key]

    # If this set becomes the default, demote the others in this spin state so
    # exactly one default remains.
    if set_block["is_default"]:
        for other in ff_sets.values():
            other["is_default"] = False
    ff_sets[set_name] = set_block


def _validate_via_loader(cofactor_path: str, redox_state: str,
                         spin_state: str, set_name: str) -> None:
    """Re-read the just-written entry through the loader; raise if it's wrong.

    Uses ``load_forcefield_metadata`` + ``discover_forcefield_files`` rather
    than ``validate_forcefield_structure`` because the latter only inspects the
    bundled base path, whereas these two resolve the user base too.
    """
    try:
        loader.load_forcefield_metadata(cofactor_path)
        sets = loader.discover_forcefield_files(cofactor_path, redox_state, spin_state)
    except Exception as exc:  # noqa: BLE001 — surface any loader failure uniformly
        raise UserLibraryError(
            f"Round-trip validation failed for {cofactor_path} "
            f"[{redox_state}/{spin_state}]: {exc}"
        ) from exc

    if not any(s["name"] == set_name for s in sets):
        raise UserLibraryError(
            f"Round-trip validation failed: set '{set_name}' was written but "
            f"the loader did not return it for {cofactor_path} "
            f"[{redox_state}/{spin_state}] (check file paths in metadata)."
        )


def _ensure_dirs(state_dir: Path, created_dirs: List[Path]) -> None:
    """mkdir -p, recording dirs we create (deepest first) for rollback."""
    to_make = []
    p = state_dir
    while not p.exists():
        to_make.append(p)
        p = p.parent
    for d in reversed(to_make):
        d.mkdir()
    created_dirs.extend(to_make)  # already deepest-first


def _copy_into(src: str, dest_dir: Path, created_files: List[Path],
               replaced: List[tuple], reserved: Optional[set] = None) -> str:
    """Copy a source file into dest_dir; return the basename actually written.

    Sets are metadata KEYS, not directories, so every set in one (redox, spin)
    state shares a single folder. ``reserved`` holds the basenames other sets
    there already point at: overwriting one would destroy that set's parameters
    while its metadata still named the file. A colliding source is given a
    disambiguated name instead -- unless it is byte-identical, where sharing
    the file is what the other set already has.

    Files this call overwrites are recorded in ``replaced`` so a rollback can
    put them back; previously they were deleted outright, which turned a failed
    re-deposit into the loss of the parameters that were already there.
    """
    name = Path(src).name
    if reserved and name in reserved:
        existing = dest_dir / name
        identical = (existing.is_file()
                     and filecmp.cmp(src, existing, shallow=False))
        if not identical:
            stem, suffix = Path(name).stem, Path(name).suffix
            n = 2
            while (f"{stem}_{n}{suffix}" in reserved
                   or (dest_dir / f"{stem}_{n}{suffix}").exists()):
                n += 1
            name = f"{stem}_{n}{suffix}"

    dest = dest_dir / name
    if dest.exists():
        try:
            replaced.append((dest, dest.read_bytes()))
        except OSError:
            logger.warning("Could not snapshot %s before overwriting it", dest)
    else:
        created_files.append(dest)
    shutil.copy2(src, dest)
    return dest.name


def _reserved_basenames(metadata: Dict[str, Any], request: "PromotionRequest",
                        set_name: str) -> set:
    """Basenames owned by OTHER sets in the same (redox, spin) state."""
    reserved = set()
    sets = (metadata.get("redox_states", {})
            .get(request.redox_state, {})
            .get("spin_states", {})
            .get(request.spin_state, {})
            .get("forcefield_sets", {}))
    for name, block in sets.items():
        if name == set_name:
            continue
        for value in (block.get("files") or {}).values():
            if isinstance(value, str):
                reserved.add(value)
            elif isinstance(value, list):
                reserved.update(v for v in value if isinstance(v, str))
    return reserved


def _rollback(metadata_path: Path, original_bytes: Optional[bytes],
              created_files: List[Path], created_dirs: List[Path],
              replaced: Optional[List[tuple]] = None) -> None:
    """Undo a failed promotion: restore metadata + overwritten files, drop new
    files and empty dirs."""
    for path, content in (replaced or []):
        try:
            path.write_bytes(content)
        except OSError:
            logger.warning("Rollback: could not restore %s", path)

    for f in created_files:
        try:
            f.unlink(missing_ok=True)
        except OSError:
            logger.warning("Rollback: could not remove %s", f)

    try:
        if original_bytes is None:
            metadata_path.unlink(missing_ok=True)
        else:
            metadata_path.write_bytes(original_bytes)
    except OSError:
        logger.warning("Rollback: could not restore %s", metadata_path)

    # Remove dirs we created, deepest-first, but only if now empty.
    for d in created_dirs:
        try:
            if d.is_dir() and not any(d.iterdir()):
                d.rmdir()
        except OSError:
            pass
