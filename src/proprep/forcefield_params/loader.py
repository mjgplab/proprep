"""
Forcefield Parameter Loader

Utilities for loading and discovering forcefield parameters for redox centers
with support for multiple oxidation states and spin states.
"""

import json
import logging
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any

from proprep.utils.paths import get_package_dir

logger = logging.getLogger(__name__)


class ForcefieldNotFoundError(Exception):
    """Raised when forcefield files or metadata cannot be found."""
    pass


class InvalidForcefieldError(Exception):
    """Raised when forcefield metadata is invalid or incomplete."""
    pass


def get_forcefield_base_path() -> Path:
    """Get the base path for built-in forcefield parameters."""
    return get_package_dir() / "forcefield_params" / "specialized_residues"


def get_user_forcefield_base_path() -> Path:
    """Get the base path for user-supplied forcefield parameters (~/.proprep/)."""
    return Path.home() / ".proprep" / "forcefield_params" / "specialized_residues"


def load_forcefield_metadata(cofactor_path: str) -> Dict[str, Any]:
    """
    Load forcefield metadata for a specific cofactor type.
    
    Args:
        cofactor_path: Relative path like "heme/bis_his_c_type"
        
    Returns:
        Dictionary containing forcefield metadata
        
    Raises:
        ForcefieldNotFoundError: If metadata.json doesn't exist
        InvalidForcefieldError: If metadata is malformed
    """
    base_path = get_forcefield_base_path()
    metadata_file = base_path / cofactor_path / "metadata.json"

    if not metadata_file.exists():
        # Try user-supplied forcefield directory as fallback
        user_base = get_user_forcefield_base_path()
        user_metadata = user_base / cofactor_path / "metadata.json"
        if user_metadata.exists():
            metadata_file = user_metadata
        else:
            raise ForcefieldNotFoundError(
                f"No metadata.json found for cofactor type '{cofactor_path}' "
                f"at {metadata_file}"
            )
    
    try:
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)
        
        # Validate required fields
        required_fields = ['description', 'cofactor_type', 'redox_states']
        for field in required_fields:
            if field not in metadata:
                raise InvalidForcefieldError(
                    f"Missing required field '{field}' in {metadata_file}"
                )
        
        # Validate redox_states structure
        for redox_state, redox_data in metadata['redox_states'].items():
            if 'spin_states' not in redox_data:
                raise InvalidForcefieldError(
                    f"Missing 'spin_states' for redox state '{redox_state}' in {metadata_file}"
                )
            
            for spin_state, spin_data in redox_data['spin_states'].items():
                required_spin_fields = ['description', 'atom_types', 'residue_name']
                for field in required_spin_fields:
                    if field not in spin_data:
                        raise InvalidForcefieldError(
                            f"Missing '{field}' for {redox_state}/{spin_state} in {metadata_file}"
                        )
        
        return metadata
        
    except json.JSONDecodeError as e:
        raise InvalidForcefieldError(f"Invalid JSON in {metadata_file}: {e}")
    except Exception as e:
        raise ForcefieldNotFoundError(f"Error reading {metadata_file}: {e}")


def discover_forcefield_files(cofactor_path: str, redox_state: str, spin_state: str) -> List[Dict[str, str]]:
    """
    Discover available forcefield file sets for a specific cofactor+redox+spin combination.
    
    Args:
        cofactor_path: Relative path like "heme/bis_his_c_type"
        redox_state: Oxidation state like "fe2"
        spin_state: Spin state like "high_spin"
        
    Returns:
        List of dictionaries with forcefield file information:
        [
            {
                "name": "HEH_v1",
                "frcmod": "/path/to/HEH_v1.frcmod", 
                "lib": "/path/to/HEH_v1.lib",
                "description": "Version 1 parameters",
                "version": "1.0",
                "reference": "Smith et al. (2023)",
                "is_default": False
            }
        ]
        
    Raises:
        ForcefieldNotFoundError: If the directory doesn't exist
    """
    base_path = get_forcefield_base_path()
    user_base = get_user_forcefield_base_path()
    search_path = base_path / cofactor_path / redox_state / spin_state

    # Also check user directory
    user_search_path = user_base / cofactor_path / redox_state / spin_state

    # Try to get enhanced metadata first (load_forcefield_metadata checks both dirs)
    try:
        metadata = load_forcefield_metadata(cofactor_path)

        # Check if this redox/spin state exists in metadata
        if (redox_state in metadata['redox_states'] and
            spin_state in metadata['redox_states'][redox_state]['spin_states']):

            spin_data = metadata['redox_states'][redox_state]['spin_states'][spin_state]

            # Check if forcefield_sets is defined
            if 'forcefield_sets' in spin_data:
                forcefield_sets = []

                # Top-level forcefield_set_summary (one block per metadata file,
                # keyed by set name) provides the short, user-facing label +
                # `when to pick which` copy. Falls back to per-set description
                # if absent. See the schema in the cofactor library spec.
                set_summary_block = metadata.get("forcefield_set_summary", {})

                for set_name, set_info in spin_data['forcefield_sets'].items():
                    # Per-set is_valid flag — entries marked False are hidden
                    # from the picker (used for deprecated / pending-refactor
                    # sets whose files are still on disk but whose layout no
                    # longer matches the active transformer).
                    if set_info.get("is_valid", True) is False:
                        continue
                    # Determine which search path has the files
                    active_path = search_path
                    if not search_path.exists() and user_search_path.exists():
                        active_path = user_search_path

                    # frcmod can be a single string or a list (a metal site has
                    # a bonded frcmod plus each organic ligand's GAFF frcmod).
                    frcmod_ref = set_info['files']['frcmod']
                    if isinstance(frcmod_ref, list):
                        frcmod_paths = [active_path / fm for fm in frcmod_ref]
                        all_frcmod_exist = all(fp.exists() for fp in frcmod_paths)
                        frcmod_value = [str(fp.absolute()) for fp in frcmod_paths]
                    else:
                        frcmod_paths = [active_path / frcmod_ref]
                        all_frcmod_exist = frcmod_paths[0].exists()
                        frcmod_value = str(frcmod_paths[0].absolute())

                    # lib can be a single string or a list (MCPB multi-mol2 case)
                    lib_ref = set_info['files']['lib']
                    if isinstance(lib_ref, list):
                        # Multiple lib/mol2 files — verify all exist
                        lib_paths = [active_path / lb for lb in lib_ref]
                        all_lib_exist = all(lp.exists() for lp in lib_paths)
                        lib_value = [str(lp.absolute()) for lp in lib_paths]
                    else:
                        lib_paths = [active_path / lib_ref]
                        all_lib_exist = lib_paths[0].exists()
                        lib_value = str(lib_paths[0].absolute())

                    if all_frcmod_exist and all_lib_exist:
                        # Pull short user-facing copy from the file-level
                        # forcefield_set_summary block when this set has an
                        # entry there. Lets the picker render plain-English
                        # `when to pick which` text without the user wading
                        # through methodology paragraphs first.
                        summary_entry = set_summary_block.get(set_name, {})
                        forcefield_set = {
                            "name": set_name,
                            # display_name falls back to the key if the metadata
                            # entry doesn't carry an explicit user-facing name
                            "display_name": set_info.get("name", set_name),
                            "frcmod": frcmod_value,
                            "lib": lib_value,
                            "description": set_info.get("description", ""),
                            "methodology": set_info.get("methodology", ""),
                            "version": set_info.get("version", ""),
                            "reference": set_info.get("reference", ""),
                            "is_default": set_info.get("is_default", False),
                            "path": str(active_path),
                            # Short user-facing copy (picker prefers these over
                            # `description` when present).
                            "summary_label": summary_entry.get("label", ""),
                            "user_summary": summary_entry.get("user_summary", ""),
                            # Carry cofactor-level fields through to the picker
                            # so the methodology panel can render them without
                            # a second metadata round-trip.
                            "cofactor_path": cofactor_path,
                            "cofactor_methodology": metadata.get("methodology", ""),
                            "cofactor_set_comparison_guidance": metadata.get("set_comparison_guidance", ""),
                            # pH-treatment dimension (constant_pH vs fixed_pH).
                            # Absent ⇒ None: a legacy/treatment-agnostic set that
                            # behaves exactly as before (no protonation fork, no
                            # per-set prereq override). See resolve_residue_names
                            # and get_prerequisite_leaprc_groups(set_name=...).
                            "ph_treatment": set_info.get("ph_treatment"),
                            "protonation_model": set_info.get("protonation_model"),
                            # Per-set prerequisite override (e.g. fixed-pH sets
                            # need ff14SB/ff19SB, not leaprc.constph). None ⇒ fall
                            # back to the cofactor-global prerequisites block.
                            "set_prerequisites": set_info.get("prerequisites"),
                        }
                        forcefield_sets.append(forcefield_set)
                    else:
                        logger.warning(f"Forcefield files missing for {set_name}: frcmod={frcmod_paths}, lib={lib_paths}")

                # Metadata explicitly defines forcefield_sets, so treat it as
                # authoritative — return the surviving sets even when the list
                # is empty (every set is_valid:false, or files are missing). Do
                # NOT fall through to the directory-glob fallback below: that
                # path ignores is_valid and would re-expose sets the metadata
                # deliberately disabled (e.g. a cofactor parked pending
                # re-parameterization, with its lib/frcmod retained in-tree).
                forcefield_sets.sort(key=lambda x: (not x['is_default'], x['name']))
                return forcefield_sets

    except (KeyError, FileNotFoundError, json.JSONDecodeError, ForcefieldNotFoundError):
        # Fallback to file discovery if metadata is missing or invalid
        logger.debug(f"Using file-based discovery for {cofactor_path}/{redox_state}/{spin_state}")

    # Fallback to current file discovery method
    # Check both built-in and user paths
    if not search_path.exists():
        if user_search_path.exists():
            search_path = user_search_path
        else:
            raise ForcefieldNotFoundError(
                f"Forcefield directory not found: {search_path}"
            )
    
    forcefield_sets = []
    
    # Find all .frcmod files and look for matching .lib files
    for frcmod_file in search_path.glob("*.frcmod"):
        base_name = frcmod_file.stem  # e.g., "HEH_v1"
        lib_file = search_path / f"{base_name}.lib"
        
        if lib_file.exists():
            forcefield_set = {
                "name": base_name,
                "frcmod": str(frcmod_file.absolute()),
                "lib": str(lib_file.absolute()),
                "path": str(search_path),
                "description": "",
                "version": "",
                "reference": "",
                "is_default": False
            }
            
            # Look for optional description file
            desc_file = search_path / f"{base_name}.txt"
            if desc_file.exists():
                try:
                    with open(desc_file, 'r') as f:
                        forcefield_set["description"] = f.read().strip()
                except Exception:
                    pass
            
            forcefield_sets.append(forcefield_set)
        else:
            logger.warning(f"Found .frcmod without matching .lib: {frcmod_file}")
    
    return sorted(forcefield_sets, key=lambda x: x["name"])


def validate_forcefield_structure(cofactor_path: str) -> Tuple[bool, List[str]]:
    """
    Validate that a cofactor forcefield directory has proper structure.
    
    Args:
        cofactor_path: Relative path like "heme/bis_his_c_type"
        
    Returns:
        Tuple of (is_valid, list_of_issues)
    """
    issues = []
    base_path = get_forcefield_base_path()
    cofactor_dir = base_path / cofactor_path

    if not cofactor_dir.exists():
        # Try user directory
        user_base = get_user_forcefield_base_path()
        user_dir = user_base / cofactor_path
        if user_dir.exists():
            cofactor_dir = user_dir
        else:
            return False, [f"Cofactor directory does not exist: {cofactor_dir}"]
    
    # Check for metadata.json
    metadata_file = cofactor_dir / "metadata.json"
    if not metadata_file.exists():
        issues.append("Missing metadata.json")
        return False, issues
    
    try:
        metadata = load_forcefield_metadata(cofactor_path)
    except Exception as e:
        issues.append(f"Invalid metadata.json: {e}")
        return False, issues
    
    # Check that all redox/spin combinations have directories and files
    for redox_state, redox_data in metadata['redox_states'].items():
        for spin_state in redox_data['spin_states'].keys():
            state_dir = cofactor_dir / redox_state / spin_state
            
            if not state_dir.exists():
                issues.append(f"Missing directory: {redox_state}/{spin_state}")
                continue
            
            # Check for at least one .frcmod/.lib pair
            try:
                files = discover_forcefield_files(cofactor_path, redox_state, spin_state)
                if not files:
                    issues.append(f"No .frcmod/.lib pairs found in {redox_state}/{spin_state}")
            except Exception as e:
                issues.append(f"Error checking {redox_state}/{spin_state}: {e}")
    
    return len(issues) == 0, issues


def get_available_cofactor_types() -> Dict[str, Dict[str, Any]]:
    """
    Discover all available cofactor types with their metadata.
    
    Returns:
        Dictionary mapping cofactor paths to their metadata:
        {
            "heme/bis_his_c_type": {
                "description": "...",
                "redox_states": {...},
                "valid": True,
                "issues": []
            }
        }
    """
    base_path = get_forcefield_base_path()
    user_base = get_user_forcefield_base_path()
    cofactor_types = {}

    # Collect metadata.json files from both built-in and user directories
    metadata_files = []

    if base_path.exists():
        metadata_files.extend(
            (mf, base_path) for mf in base_path.rglob("metadata.json")
        )
    else:
        logger.warning(f"Forcefield base path does not exist: {base_path}")

    if user_base.exists():
        metadata_files.extend(
            (mf, user_base) for mf in user_base.rglob("metadata.json")
        )

    for metadata_file, source_base in metadata_files:
        # Get relative path from source base to the cofactor directory
        cofactor_dir = metadata_file.parent
        relative_path = cofactor_dir.relative_to(source_base)
        cofactor_path = str(relative_path).replace("\\", "/")  # Ensure forward slashes
        
        try:
            metadata = load_forcefield_metadata(cofactor_path)
            is_valid, issues = validate_forcefield_structure(cofactor_path)
            
            cofactor_types[cofactor_path] = {
                "description": metadata.get("description", ""),
                "cofactor_type": metadata.get("cofactor_type", ""),
                "redox_states": metadata.get("redox_states", {}),
                "valid": is_valid,
                "issues": issues,
                "path": str(cofactor_dir)
            }
        except Exception as e:
            cofactor_types[cofactor_path] = {
                "description": f"Error loading metadata",
                "cofactor_type": "",
                "redox_states": {},
                "valid": False,
                "issues": [str(e)],
                "path": str(cofactor_dir)
            }
    
    return cofactor_types


def get_registered_residue_names() -> Dict[str, str]:
    """Collect every tLEaP unit name already claimed by the forcefield library.

    Walks both the bundled and user libraries and returns the residue/unit names
    they define — the ``residue_name`` at each spin state (a bare string for a
    single residue, or the values of a source->code map for multi-residue metal
    sites and covalent adducts) plus any ``ligand_residue_names``. These are
    exactly the names tLEaP will resolve once the library is loaded, so a newly
    parameterized residue must avoid them: two units saved under the same name
    collide in tLEaP's namespace and the second silently shadows the first.

    Returns:
        Dict mapping each uppercased name -> the cofactor path that owns it
        (first owner wins on the unlikely event of an intra-library duplicate).
        Never raises — a library that cannot be read yields an empty dict, so
        callers can treat "no known names" and "could not look up names" alike.
    """
    names: Dict[str, str] = {}

    def _claim(value: Any, owner: str) -> None:
        # residue_name is a str for a single residue or a {source: code} map for
        # metal sites / covalent adducts; take the code(s) either way.
        candidates = value.values() if isinstance(value, dict) else [value]
        for cand in candidates:
            if isinstance(cand, str) and cand.strip():
                names.setdefault(cand.strip().upper(), owner)

    try:
        cofactor_types = get_available_cofactor_types()
    except Exception as exc:  # noqa: BLE001 — name lookup must never break callers
        logger.warning("Could not enumerate registered residue names: %s", exc)
        return names

    for cofactor_path, meta in cofactor_types.items():
        for redox in (meta.get("redox_states") or {}).values():
            for spin in (redox.get("spin_states") or {}).values():
                if not isinstance(spin, dict):
                    continue
                _claim(spin.get("residue_name"), cofactor_path)
                _claim(spin.get("ligand_residue_names"), cofactor_path)

    return names


def get_atom_types_for_state(cofactor_path: str, redox_state: str, spin_state: str) -> List[str]:
    """
    Get atom types for a specific cofactor+redox+spin combination.
    
    Args:
        cofactor_path: Relative path like "heme/bis_his_c_type"
        redox_state: Oxidation state like "fe2"
        spin_state: Spin state like "high_spin"
        
    Returns:
        List of atom type strings like ["F2H", "N1H", "N2H"]
        
    Raises:
        ForcefieldNotFoundError: If metadata or state not found
        InvalidForcefieldError: If atom_types not properly defined
    """
    metadata = load_forcefield_metadata(cofactor_path)
    
    if redox_state not in metadata['redox_states']:
        raise ForcefieldNotFoundError(
            f"Redox state '{redox_state}' not found for {cofactor_path}"
        )
    
    redox_data = metadata['redox_states'][redox_state]
    if spin_state not in redox_data['spin_states']:
        raise ForcefieldNotFoundError(
            f"Spin state '{spin_state}' not found for {cofactor_path}/{redox_state}"
        )
    
    spin_data = redox_data['spin_states'][spin_state]
    atom_types = spin_data.get('atom_types', [])
    
    if not atom_types:
        raise InvalidForcefieldError(
            f"No atom_types defined for {cofactor_path}/{redox_state}/{spin_state}"
        )
    
    return atom_types


def get_residue_name_for_state(cofactor_path: str, redox_state: str, spin_state: str) -> str:
    """
    Get the transformed residue name for a specific cofactor+redox+spin combination.
    
    Args:
        cofactor_path: Relative path like "heme/bis_his_c_type"
        redox_state: Oxidation state like "fe2"
        spin_state: Spin state like "high_spin"
        
    Returns:
        Residue name string like "HEH"
        
    Raises:
        ForcefieldNotFoundError: If metadata or state not found
    """
    metadata = load_forcefield_metadata(cofactor_path)
    
    if redox_state not in metadata['redox_states']:
        raise ForcefieldNotFoundError(
            f"Redox state '{redox_state}' not found for {cofactor_path}"
        )
    
    redox_data = metadata['redox_states'][redox_state]
    if spin_state not in redox_data['spin_states']:
        raise ForcefieldNotFoundError(
            f"Spin state '{spin_state}' not found for {cofactor_path}/{redox_state}"
        )
    
    spin_data = redox_data['spin_states'][spin_state]
    return spin_data.get('residue_name', '')


def _get_spin_data(metadata: Dict[str, Any], cofactor_path: str,
                   redox_state: str, spin_state: str) -> Dict[str, Any]:
    """Navigate metadata to a spin-state block, raising a clear error if absent."""
    redox_states = metadata.get('redox_states', {})
    if redox_state not in redox_states:
        raise ForcefieldNotFoundError(
            f"Redox state '{redox_state}' not found for {cofactor_path}"
        )
    spin_states = redox_states[redox_state].get('spin_states', {})
    if spin_state not in spin_states:
        raise ForcefieldNotFoundError(
            f"Spin state '{spin_state}' not found for {cofactor_path}/{redox_state}"
        )
    return spin_states[spin_state]


def get_protonation_model(cofactor_path: str, redox_state: str, spin_state: str,
                          set_name: str) -> Optional[Dict[str, Any]]:
    """
    Return the ``protonation_model`` block for one forcefield set, or None.

    The protonation model describes a set's titratable/protonation sites — each
    a ``role`` (matching the transformer's residue-space-plan keys, e.g.
    ``propionate_a``) plus either a ``variants`` map (fixed-pH: static protomer
    choices) or a single ``titratable_residue`` (constant-pH). Sets that don't
    declare one (every current cofactor) return None.

    Args:
        cofactor_path: Relative path like "heme/cys_axial_b_type"
        redox_state, spin_state: The state whose sets to look in
        set_name: The forcefield-set key (e.g. "Guberman_HTO_RESP_FixedpH")

    Returns:
        The ``protonation_model`` dict, or None if the set or block is absent.

    Raises:
        ForcefieldNotFoundError: If metadata / state cannot be found.
    """
    metadata = load_forcefield_metadata(cofactor_path)
    spin_data = _get_spin_data(metadata, cofactor_path, redox_state, spin_state)
    set_info = (spin_data.get('forcefield_sets', {}) or {}).get(set_name)
    if not set_info:
        return None
    return set_info.get('protonation_model')


def resolve_residue_names(cofactor_path: str, redox_state: str, spin_state: str,
                          set_name: str,
                          protonation_choices: Optional[Dict[str, str]] = None
                          ) -> Dict[str, str]:
    """
    Resolve every output residue name a transformer produces, keyed by role.

    This is the single source of truth for the residue codes a transformer
    writes into the structure — so transformers hardcode none of them. Roles
    are drawn from three places in metadata, all merged into one ``{role:
    resname}`` dict:

    - ``center``      ← the spin-state ``residue_name`` (e.g. HEM → HTO)
    - ligand roles    ← the spin-state ``ligand_residue_names`` map
                        (e.g. ``proximal_cys`` → CTO)
    - protonation roles ← the chosen set's ``protonation_model.sites``:
        * fixed-pH site: ``variants[choice]`` where ``choice`` comes from
          ``protonation_choices[role]`` (falling back to the site ``default``)
        * constant-pH site: the single ``titratable_residue``

    Args:
        cofactor_path: Relative path like "heme/cys_axial_b_type"
        redox_state, spin_state: The chosen state
        set_name: The chosen forcefield-set key
        protonation_choices: role → variant-key picks for fixed-pH sites
            (e.g. {"propionate_a": "protonated", "propionate_d": "deprotonated"}).
            Missing entries fall back to each site's ``default``.

    Returns:
        Dict mapping role → residue name, e.g.
        {"center": "HTO", "proximal_cys": "CTO",
         "propionate_a": "PRP", "propionate_d": "PRD"}.

    Raises:
        ForcefieldNotFoundError: If metadata / state cannot be found.
        InvalidForcefieldError: If a requested protonation choice is not a
            valid variant for its site, or a site is malformed.
    """
    protonation_choices = protonation_choices or {}
    metadata = load_forcefield_metadata(cofactor_path)
    spin_data = _get_spin_data(metadata, cofactor_path, redox_state, spin_state)

    resolved: Dict[str, str] = {}

    # center
    center = spin_data.get('residue_name')
    if center:
        resolved['center'] = center

    # ligand roles (proximal_cys, axial_his_stub, ...)
    for role, resname in (spin_data.get('ligand_residue_names', {}) or {}).items():
        resolved[role] = resname

    # protonation-site roles (from the chosen set's protonation_model)
    set_info = (spin_data.get('forcefield_sets', {}) or {}).get(set_name)
    if not set_info:
        raise InvalidForcefieldError(
            f"Forcefield set '{set_name}' not found for "
            f"{cofactor_path}/{redox_state}/{spin_state}"
        )
    model = set_info.get('protonation_model')
    if model:
        for site in model.get('sites', []) or []:
            role = site.get('role')
            if not role:
                raise InvalidForcefieldError(
                    f"protonation_model site missing 'role' in set '{set_name}' "
                    f"({cofactor_path}/{redox_state}/{spin_state})"
                )
            if 'variants' in site:
                variants = site['variants'] or {}
                choice = protonation_choices.get(role, site.get('default'))
                if choice is None:
                    raise InvalidForcefieldError(
                        f"No protonation choice or default for site '{role}' "
                        f"in set '{set_name}'"
                    )
                if choice not in variants:
                    raise InvalidForcefieldError(
                        f"Invalid protonation state '{choice}' for site '{role}' "
                        f"in set '{set_name}'. Valid: {sorted(variants)}"
                    )
                resolved[role] = variants[choice]
            elif 'titratable_residue' in site:
                resolved[role] = site['titratable_residue']
            else:
                raise InvalidForcefieldError(
                    f"protonation_model site '{role}' in set '{set_name}' has "
                    "neither 'variants' nor 'titratable_residue'"
                )

    return resolved


def _parse_leaprc_groups(prereqs: Dict[str, Any]) -> List[List[str]]:
    """Parse a ``prerequisites`` dict into AND-groups of OR'd satisfiers.

    Preferred form ``leaprc_groups`` (each entry a list or a dict with
    ``satisfied_by``); legacy flat ``leaprcs`` becomes one size-1 group per
    entry (strict "every listed leaprc required"). Empty/missing ⇒ ``[]``.
    """
    prereqs = prereqs or {}

    if "leaprc_groups" in prereqs:
        groups: List[List[str]] = []
        for entry in prereqs.get("leaprc_groups", []) or []:
            if isinstance(entry, dict):
                satisfied_by = list(entry.get("satisfied_by", []) or [])
            elif isinstance(entry, list):
                satisfied_by = list(entry)
            else:
                continue
            if satisfied_by:
                groups.append(satisfied_by)
        return groups

    # Legacy flat list: each leaprc is its own required group.
    return [[lr] for lr in (prereqs.get("leaprcs", []) or [])]


def _find_set_prerequisites(metadata: Dict[str, Any], set_name: str) -> Optional[Dict[str, Any]]:
    """Find a forcefield set's own ``prerequisites`` block by set name.

    Walks every redox_state/spin_state, since a set name (e.g.
    ``Guberman_HTO_RESP_FixedpH``) is unique to one state. Returns the first
    matching set's ``prerequisites`` dict, or None if the set isn't found or
    declares none (⇒ caller falls back to the cofactor-global block).
    """
    for redox_data in (metadata.get("redox_states", {}) or {}).values():
        for spin_data in (redox_data.get("spin_states", {}) or {}).values():
            sets = spin_data.get("forcefield_sets", {}) or {}
            if set_name in sets:
                prereqs = sets[set_name].get("prerequisites")
                return prereqs if prereqs else None
    return None


def get_prerequisite_leaprc_groups(cofactor_path: str,
                                   set_name: Optional[str] = None) -> List[List[str]]:
    """
    Return prerequisite leaprcs as AND-groups of OR'd satisfiers.

    Each inner list is one AND-group: at least one of its leaprcs must be
    sourced to satisfy that group. ALL groups must be satisfied (AND across
    groups) for the cofactor's lib/frcmod to load cleanly.

    This is the primary API for the FF picker's "fully-satisfying" rule: a
    picker option fully satisfies a cofactor iff its leaprc(s) hit at least
    one entry in every AND-group.

    Schema: metadata's ``prerequisites.leaprc_groups`` is the preferred form;
    each group is either a list of strings or a dict with a ``satisfied_by``
    list. Legacy ``prerequisites.leaprcs`` (flat list) is converted by
    treating every entry as its own AND-group of size 1 — preserving the
    strict "every listed leaprc is required" semantics that pre-groups code
    had.

    Per-set override: when ``set_name`` is given and that forcefield set
    declares its own ``prerequisites`` block, those groups are returned
    instead of the cofactor-global block. This is what lets a fixed-pH set
    require ff14SB/ff19SB while the constant-pH default set of the same
    cofactor requires leaprc.constph. A set without its own prerequisites
    falls back to the cofactor-global block (current behavior).

    Args:
        cofactor_path: Relative path like "flavin/fmn" or "heme/bis_his_b_type"
        set_name: Optional forcefield-set name to read per-set prerequisites
            from. None ⇒ cofactor-global prerequisites.

    Returns:
        List of AND-groups, each an OR list of satisfying leaprc strings.
        Empty list if no applicable prerequisites are declared.

    Raises:
        ForcefieldNotFoundError: If metadata file cannot be found.
    """
    try:
        metadata = load_forcefield_metadata(cofactor_path)
    except ForcefieldNotFoundError:
        # Caller may pass a path that no longer exists; treat as no prereqs.
        return []

    if set_name is not None:
        set_prereqs = _find_set_prerequisites(metadata, set_name)
        if set_prereqs is not None:
            return _parse_leaprc_groups(set_prereqs)
        # Set has no own prerequisites → fall through to cofactor-global.

    return _parse_leaprc_groups(metadata.get("prerequisites", {}) or {})


_PH_TREATMENT_SUFFIX_RE = re.compile(
    r"_(?:constp[hH]|fixedp[hH])$", re.IGNORECASE)


def _set_base_name(set_name: str) -> str:
    """Strip a trailing pH-treatment suffix from a set name to get its base.

    Companion sets are linked only by naming convention — the constant-pH and
    fixed-pH variants of one parameterization share a base name and differ only
    in the ``_constpH`` / ``_FixedpH`` suffix (e.g. ``Henriques_HCO_RESP_constpH``
    ↔ ``Henriques_HCO_RESP_FixedpH``). There is no explicit companion field in
    the metadata, so the shared base (which also encodes the charge model, e.g.
    ``RESP``) is how we match a set to its opposite-treatment twin.
    """
    return _PH_TREATMENT_SUFFIX_RE.sub("", set_name)


def find_companion_set(cofactor_path: str, redox_state: str, spin_state: str,
                       set_name: str,
                       target_treatment: str = "fixed_pH") -> Optional[str]:
    """Find the ``target_treatment`` companion of a forcefield set.

    Used when a structure's propionates have been switched from the titratable
    PRN to static PRP/PRD (e.g. by the PB-Titrate modern-FF rebuild): the
    recorded constant-pH set must be swapped for its fixed-pH twin so the
    Topology Generator loads the lib that actually defines PRP/PRD and applies
    the ff14SB/ff19SB (not leaprc.constph/conste) prerequisites.

    Resolution, in order of preference:
      1. the set with ``ph_treatment == target_treatment`` whose base name
         (suffix stripped) matches ``set_name``'s base — i.e. the same
         parameterization and charge model;
      2. the ``is_default`` set of the target treatment;
      3. the first set of the target treatment.

    Args:
        cofactor_path: Relative path like "heme/bis_his_c_type".
        redox_state, spin_state: The state whose sets to search.
        set_name: The currently-recorded set (its opposite-treatment twin is
            sought). May be None/empty — then only the default/first fallback
            applies.
        target_treatment: The ph_treatment to switch to (default "fixed_pH").

    Returns:
        The companion set's name, or None if the cofactor defines no set of
        ``target_treatment`` for this state (caller should warn — the structure
        names something the available sets can't build).
    """
    try:
        options = discover_forcefield_files(cofactor_path, redox_state, spin_state)
    except ForcefieldNotFoundError:
        return None

    candidates = [o for o in options
                  if o.get("ph_treatment") == target_treatment]
    if not candidates:
        return None

    if set_name:
        base = _set_base_name(set_name)
        for o in candidates:
            if _set_base_name(o["name"]) == base:
                return o["name"]

    for o in candidates:
        if o.get("is_default"):
            return o["name"]
    return candidates[0]["name"]


def get_prerequisite_leaprcs(cofactor_path: str) -> List[str]:
    """
    Return the flat union of leaprcs that could satisfy any prereq group for
    this cofactor.

    Prefer :py:func:`get_prerequisite_leaprc_groups` for picker logic — this
    flat list collapses AND/OR structure and is only useful for callers that
    just want "which leaprcs are mentioned anywhere in this cofactor's
    prerequisites" (e.g. building a quick "candidate options" set).

    Args:
        cofactor_path: Relative path like "flavin/fmn" or "heme/bis_his_b_type"

    Returns:
        List of leaprc strings, deduplicated while preserving first-seen
        order. Empty list if the cofactor declares no prerequisites.

    Raises:
        ForcefieldNotFoundError: If metadata file cannot be found.
    """
    seen: set = set()
    flat: List[str] = []
    for group in get_prerequisite_leaprc_groups(cofactor_path):
        for lr in group:
            if lr not in seen:
                seen.add(lr)
                flat.append(lr)
    return flat


# Cache: residue name -> bool (does any bundled .lib unit of this name carry
# head/tail connect atoms?). Populated lazily on first call. None until built.
_TEMPLATE_CONNECTIVITY_CACHE: Optional[Dict[str, bool]] = None


def _scan_lib_for_connectivity(lib_path: Path) -> Dict[str, bool]:
    """Parse one AMBER .lib file and return ``{unit_name: has_head_tail}``.

    A unit "has head/tail" iff its ``unit.connect`` block contains two
    non-zero atom seq numbers (i.e. tleap will auto-form peptide bonds to
    neighbours when the residue is loaded as part of a sequence).
    """
    result: Dict[str, bool] = {}
    try:
        text = lib_path.read_text()
    except OSError as e:
        logger.debug("Skipping unreadable lib %s: %s", lib_path, e)
        return result

    lines = text.splitlines()
    n = len(lines)
    i = 0
    while i < n:
        line = lines[i]
        # Match "!entry.<unit>.unit.connect array int"
        if line.startswith("!entry.") and line.rstrip().endswith(".unit.connect array int"):
            # Extract unit name between "!entry." and ".unit.connect"
            try:
                unit = line[len("!entry."):-len(".unit.connect array int")].rstrip()
            except Exception:
                i += 1
                continue
            head = tail = "0"
            if i + 1 < n:
                head = lines[i + 1].strip()
            if i + 2 < n:
                tail = lines[i + 2].strip()
            try:
                has_connectivity = int(head) > 0 and int(tail) > 0
            except ValueError:
                has_connectivity = False
            # First definition wins if a residue appears in multiple libs.
            result.setdefault(unit, has_connectivity)
            i += 3
            continue
        i += 1
    return result


def _build_template_connectivity_cache() -> Dict[str, bool]:
    """Walk all bundled + user .lib files and merge unit -> has-head-tail."""
    merged: Dict[str, bool] = {}
    roots = [get_forcefield_base_path(), get_user_forcefield_base_path()]
    for root in roots:
        if not root.exists():
            continue
        for lib_path in root.rglob("*.lib"):
            for unit, has_ht in _scan_lib_for_connectivity(lib_path).items():
                # Once any lib reports head/tail for this residue, treat it as
                # connected (libs we ship can vary, but the safer bet for
                # avoiding duplicate explicit `bond` statements is to trust the
                # most-connective definition we see).
                if has_ht:
                    merged[unit] = True
                else:
                    merged.setdefault(unit, False)
    return merged


def residue_has_template_connectivity(resname: str) -> bool:
    """Return True if a bundled (or user) .lib defines ``resname`` with both
    head (connect0) and tail (connect1) atoms — meaning tleap will auto-form
    peptide bonds to its sequence neighbours and ProPrep should NOT emit an
    explicit ``bond mol.<i>.C mol.<i+1>.N`` for this residue.

    Result cached for the life of the process.
    """
    global _TEMPLATE_CONNECTIVITY_CACHE
    if _TEMPLATE_CONNECTIVITY_CACHE is None:
        _TEMPLATE_CONNECTIVITY_CACHE = _build_template_connectivity_cache()
    return _TEMPLATE_CONNECTIVITY_CACHE.get(resname, False)


def reset_template_connectivity_cache() -> None:
    """Drop the cached lib scan. Useful in tests or after copying new libs in."""
    global _TEMPLATE_CONNECTIVITY_CACHE
    _TEMPLATE_CONNECTIVITY_CACHE = None