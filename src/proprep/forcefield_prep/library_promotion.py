"""
Interactive promotion of parameterized residues into the user library.

This is the UI layer over :func:`proprep.forcefield_params.promote_state`. It
turns a finished parameterization (or a set of loose, previously-developed
files) into a discoverable entry under
``~/.proprep/forcefield_params/specialized_residues/`` so it extends the
bundled library and can be reused in future projects.

Two entry points:
  * :func:`offer_promotion` — called at the end of a successful
    parameterization run. Opt-in (default no); prefills everything the
    workflow already knows and asks only for what it can't.
  * :func:`run_import_wizard` — a standalone "import existing parameters"
    action for ``.frcmod`` / ``.lib`` files developed at another time or place.

Design choices baked in here (see the user-library plan):
  * Residues are filed by *parameterization route*, not chemistry:
    ``small_molecules/<name>``, ``modified_aa/<name>``, ``metal_sites/<name>``.
  * The parameterizers produce ONE state at a time, so the wizard always
    deposits a single ``single_state/default`` leaf for small molecules and
    modified AAs, and a single (redox, spin) leaf for metal sites — the writer
    merges it into any existing metadata for that residue.
  * Only metal sites can introduce new atom types; small molecules (GAFF) and
    modified AAs (standard / GAFF gap-fill) never do, so ``atom_types`` is left
    empty for them without asking.
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from rich.console import Console
from rich.panel import Panel

from proprep.forcefield_params import (
    promote_state,
    PromotionRequest,
    LibraryCollisionError,
    UserLibraryError,
)
from proprep.forcefield_params import user_library
from proprep.utils.prompts import prompt_with_context, confirm_with_context

logger = logging.getLogger(__name__)

_MODULE = "User Library"

# Parameterization category → user-library family (top-level path segment).
_CATEGORY_FAMILY = {
    "small_molecule": "small_molecules",
    "modified_amino_acid": "modified_aa",
    "modified_protein": "modified_aa",
    "metal_site": "metal_sites",
    "organometallic_cofactor": "metal_sites",
}


def find_lib_artifact(
    search_dir: Union[str, Path],
    residue_name: Optional[str] = None,
    parameter_files: Optional[Dict[str, Any]] = None,
) -> Optional[str]:
    """Locate the OFF library a parameterization left behind.

    The small-molecule / modified-AA workflows ``saveoff`` a ``.lib`` next to
    their other outputs but don't always record its path in
    ``parameter_files``. Prefer an explicit ``lib_file`` entry, then a ``.lib``/
    ``.off`` in ``search_dir`` whose stem matches ``residue_name``, then any.
    """
    if parameter_files and parameter_files.get("lib_file"):
        return parameter_files["lib_file"]
    d = Path(search_dir)
    if not d.is_dir():
        return None
    libs = sorted(d.glob("*.lib")) + sorted(d.glob("*.off"))
    if not libs:
        return None
    if residue_name:
        for p in libs:
            if p.stem.lower() == str(residue_name).lower():
                return str(p.resolve())
    return str(libs[0].resolve())


def offer_promotion(
    console: Console,
    processor,
    *,
    category: str,
    residue_name: Union[str, Dict[str, str]],
    parameter_files: Dict[str, Any],
    atom_types: Optional[List[str]] = None,
) -> Optional[Dict[str, Any]]:
    """Offer to save a just-finished parameterization into the user library.

    ``parameter_files`` is the dict each parameterizer returns, with absolute
    paths under keys ``frcmod_file``, ``lib_file``, ``prep_file`` (any subset).
    Returns the writer result dict on success, or ``None`` if the user declines,
    cancels, or no usable artifacts are present.
    """
    frcmod = parameter_files.get("frcmod_file")
    lib = parameter_files.get("lib_file")
    prep = parameter_files.get("prep_file")

    if not frcmod or not Path(frcmod).is_file():
        # Nothing to promote (e.g. a paused/incomplete run).
        return None

    console.print()
    console.print(Panel(
        "Save these parameters to your personal force-field library?\n"
        "[grey50]They will extend the built-in library and become available "
        "in every future ProPrep project (via the Force Field Explorer and "
        "Topology Generator), not just this one.[/grey50]",
        title="Reuse these parameters?", border_style="green",
    ))
    if not confirm_with_context(
        processor, "Add to your library", default=False,
        module=_MODULE, description="Promote parameters to the user library",
    ):
        return None

    request = _build_request_interactively(
        console, processor,
        category=category, residue_name=residue_name,
        frcmod_file=frcmod, lib_file=lib, prep_file=prep,
        atom_types=atom_types,
    )
    if request is None:
        return None
    return _promote_with_collision_retry(console, processor, request)


def offer_library_promotion(console, processor, *, category, residue_name,
                            frcmod_file, lib_search_dir, prep_file=None,
                            lib_file=None, atom_types=None):
    """Offer to save a finished parameterization into the user library.

    Shared by the final ("Force Field Integration") step of every parameterizer
    checklist, so each one deposits the same way. Resolves the .lib artifact on
    disk first, since parameterizers name it by residue rather than returning it.

    Best-effort and fully guarded: a failure here must never break the
    parameterization, which has already succeeded and stored its results.
    Returns the deposit result dict (with ``library_path``/``metadata_path``)
    when a promotion occurred, else None — callers use it to point a reuse
    transformer at the deposited library.
    """
    try:
        resolved_lib = find_lib_artifact(
            lib_search_dir, residue_name=residue_name,
            parameter_files={"lib_file": lib_file} if lib_file else None,
        )
        return offer_promotion(
            console, processor,
            category=category,
            residue_name=residue_name,
            parameter_files={
                "frcmod_file": frcmod_file,
                "lib_file": resolved_lib,
                "prep_file": prep_file,
            },
            atom_types=atom_types,
        )
    except Exception as e:  # noqa: BLE001 — never let promotion break the workflow
        logger.debug("Library promotion offer skipped: %s", e)
        return None


def run_import_wizard(
    console: Console,
    processor,
) -> Optional[Dict[str, Any]]:
    """Standalone: import previously-developed .frcmod/.lib files.

    Lets the user register parameters they built elsewhere (a prior project, a
    collaborator, a published set) into their library through the same metadata
    wizard the post-parameterization path uses.
    """
    console.print(Panel(
        "Import existing force-field parameters into your personal library.\n"
        "[grey50]Point ProPrep at a finished .frcmod and its .lib/.off (and an "
        "optional .prep) and it will be registered for reuse in any "
        "project.[/grey50]",
        title="Import parameters", border_style="cyan",
    ))

    frcmod = _prompt_existing_file(
        console, processor, "Path to the .frcmod file",
        suffixes=(".frcmod",), required=True,
    )
    if frcmod is None:
        return None
    lib = _prompt_existing_file(
        console, processor, "Path to the .lib / .off library file",
        suffixes=(".lib", ".off"), required=True,
    )
    if lib is None:
        return None
    prep = _prompt_existing_file(
        console, processor, "Path to a .prep file (optional, press Enter to skip)",
        suffixes=(".prep", ".prepi", ".prepin"), required=False,
    )

    # Category for an imported set is whatever the user says it is.
    category = _prompt_category(console, processor)
    default_name = Path(lib).stem

    request = _build_request_interactively(
        console, processor,
        category=category, residue_name=default_name,
        frcmod_file=frcmod, lib_file=lib, prep_file=prep,
        atom_types=None,
    )
    if request is None:
        return None
    return _promote_with_collision_retry(console, processor, request)


# --------------------------------------------------------------------------- #
# request building
# --------------------------------------------------------------------------- #

def _build_request_interactively(
    console: Console,
    processor,
    *,
    category: str,
    residue_name: Union[str, Dict[str, str]],
    frcmod_file: str,
    lib_file: Optional[str],
    prep_file: Optional[str],
    atom_types: Optional[List[str]],
) -> Optional[PromotionRequest]:
    family = _CATEGORY_FAMILY.get(category, "small_molecules")
    is_metal = family == "metal_sites"

    # A single residue name is a string; metal sites may hand us a rename map.
    primary_name = (next(iter(residue_name.values())) if isinstance(residue_name, dict)
                    else residue_name)

    # type_name: the path segment + identity of the entry.
    type_label = "site identifier" if is_metal else "residue name"
    type_name = prompt_with_context(
        processor, f"Library entry name ({type_label})",
        default=str(primary_name),
        module=_MODULE, description=f"User-library {type_label}",
    ).strip()
    if not type_name:
        type_name = str(primary_name)

    description = prompt_with_context(
        processor, "One-line description",
        default=f"{primary_name} ({_default_method_blurb(family)})",
        module=_MODULE, description="Residue description",
    ).strip()

    set_name = prompt_with_context(
        processor, "Parameter set name",
        default=f"user_{type_name}_RESP" if is_metal else f"user_{type_name}",
        module=_MODULE, description="Force-field set name",
    ).strip()

    references = _prompt_references(console, processor)

    # State coordinates: degenerate for non-metals, real dimensions for metals.
    if is_metal:
        redox_state = prompt_with_context(
            processor, "Redox state label", default="as_parameterized",
            module=_MODULE, description="Redox state",
        ).strip() or "as_parameterized"
        spin_state = prompt_with_context(
            processor, "Spin state label", default="default",
            module=_MODULE, description="Spin state",
        ).strip() or "default"
        state_meta = {}
        charge = _prompt_optional_int(console, processor, "Formal charge of the site")
        if charge is not None:
            state_meta["formal_charge"] = charge
        resolved_atom_types = list(atom_types or [])
    else:
        redox_state = user_library.DEFAULT_REDOX_STATE
        spin_state = user_library.DEFAULT_SPIN_STATE
        state_meta = {}
        resolved_atom_types = []  # GAFF / standard types — nothing new to declare

    lib_srcs = [lib_file] if lib_file else []
    if not lib_srcs:
        console.print("[yellow]No .lib/.off library file is available — a residue "
                      "library is required for the loader to use these parameters. "
                      "Cannot promote.[/yellow]")
        return None

    extra_files = [prep_file] if prep_file and prep_file not in lib_srcs else []

    spin_meta: Dict[str, Any] = {
        "residue_name": residue_name,   # str or dict, preserved as-is
        "atom_types": resolved_atom_types,
    }

    return PromotionRequest(
        family=family,
        type_name=type_name,
        set_name=set_name,
        frcmod_src=frcmod_file,
        lib_srcs=lib_srcs,
        redox_state=redox_state,
        spin_state=spin_state,
        residue_meta={"description": description, "references": references},
        state_meta=state_meta,
        spin_meta=spin_meta,
        set_meta={"description": description, "version": "1.0"},
        extra_file_srcs=extra_files,
    )


def _promote_with_collision_retry(
    console: Console,
    processor,
    request: PromotionRequest,
) -> Optional[Dict[str, Any]]:
    """Call the writer; on a collision, let the user pick how to resolve it."""
    while True:
        try:
            result = promote_state(request)
        except LibraryCollisionError as exc:
            console.print(f"[yellow]{exc}[/yellow]")
            choice = prompt_with_context(
                processor, "How would you like to proceed?",
                choices=["1", "2", "3"], default="3",
                module=_MODULE, description="Resolve library collision",
                options_map={
                    "1": "Save as a new version (keep the existing one)",
                    "2": "Overwrite the existing set",
                    "3": "Cancel",
                },
            )
            if choice == "1":
                request.on_collision = user_library.ON_COLLISION_VERSION
                continue
            if choice == "2":
                request.on_collision = user_library.ON_COLLISION_OVERWRITE
                continue
            console.print("[grey50]Promotion cancelled.[/grey50]")
            return None
        except (UserLibraryError, FileNotFoundError) as exc:
            console.print(f"[red]Could not save to library: {exc}[/red]")
            return None

        console.print(Panel(
            f"[green]Saved to your library.[/green]\n"
            f"[grey50]{result['library_path']}\n"
            f"set: {result['set_name']}[/grey50]\n"
            "It will now appear in the Force Field Explorer and Topology "
            "Generator in any project.",
            title="Parameters saved", border_style="green",
        ))
        return result


# --------------------------------------------------------------------------- #
# small prompt helpers
# --------------------------------------------------------------------------- #

def _default_method_blurb(family: str) -> str:
    return {
        "small_molecules": "GAFF, AM1-BCC/RESP charges, user-parameterized",
        "modified_aa": "modified amino acid, user-parameterized",
        "metal_sites": "MCPB metal site, user-parameterized",
    }.get(family, "user-parameterized")


def _prompt_references(console: Console, processor) -> List[str]:
    raw = prompt_with_context(
        processor,
        "Reference(s) / citation (optional; separate multiple with ' | ')",
        default="",
        module=_MODULE, description="Literature references",
    ).strip()
    if not raw:
        return []
    return [r.strip() for r in raw.split("|") if r.strip()]


def _prompt_optional_int(console: Console, processor, label: str) -> Optional[int]:
    raw = prompt_with_context(
        processor, f"{label} (optional, press Enter to skip)", default="",
        module=_MODULE, description=label,
    ).strip()
    if not raw:
        return None
    try:
        return int(raw)
    except ValueError:
        console.print(f"[yellow]'{raw}' is not an integer — skipping.[/yellow]")
        return None


def _prompt_category(console: Console, processor) -> str:
    choice = prompt_with_context(
        processor, "What kind of parameters are these?",
        choices=["1", "2", "3"], default="1",
        module=_MODULE, description="Parameter category",
        options_map={
            "1": "Small molecule / ligand",
            "2": "Modified amino acid",
            "3": "Metal site",
        },
    )
    return {"1": "small_molecule", "2": "modified_amino_acid",
            "3": "metal_site"}[choice]


def _prompt_existing_file(
    console: Console,
    processor,
    label: str,
    *,
    suffixes: tuple,
    required: bool,
) -> Optional[str]:
    """Prompt for a file path, validating existence and (loosely) suffix."""
    while True:
        raw = prompt_with_context(
            processor, label, default="",
            module=_MODULE, description=label,
        ).strip()
        if not raw:
            if required:
                console.print("[yellow]A path is required.[/yellow]")
                if not confirm_with_context(processor, "Try again?", default=True,
                                            module=_MODULE,
                                            description="Retry file entry"):
                    return None
                continue
            return None
        path = Path(raw).expanduser()
        if not path.is_file():
            console.print(f"[yellow]No such file: {path}[/yellow]")
            if not confirm_with_context(processor, "Try again?", default=True,
                                        module=_MODULE,
                                        description="Retry file entry"):
                return None
            continue
        if suffixes and path.suffix.lower() not in suffixes:
            if not confirm_with_context(
                processor,
                f"{path.name} doesn't end in {'/'.join(suffixes)} — use it anyway?",
                default=False, module=_MODULE, description="Confirm unusual suffix",
            ):
                continue
        return str(path.resolve())
