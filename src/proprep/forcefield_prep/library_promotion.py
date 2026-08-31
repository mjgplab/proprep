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
import os
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from rich.console import Console
from rich.markup import escape
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
        console, processor, ".frcmod file",
        suffixes=(".frcmod",), required=True,
    )
    if frcmod is None:
        return None
    # The companion files are almost always beside the frcmod, so start there
    # rather than sending the user back to the working directory.
    beside = str(Path(frcmod).parent)
    lib = _prompt_existing_file(
        console, processor, ".lib / .off library file",
        suffixes=(".lib", ".off"), required=True, start_dir=beside,
    )
    if lib is None:
        return None
    prep = _prompt_existing_file(
        console, processor, ".prep file (optional)",
        suffixes=(".prep", ".prepi", ".prepin"), required=False, start_dir=beside,
    )

    # Category for an imported set is whatever the user says it is.
    category = _prompt_category(console, processor)

    # The residue name recorded for this entry has to be the name tLEaP will
    # match, which is the unit inside the library -- not the filename stem.
    # They diverge often enough to matter (a GDP.lib holding unit "gdp"), and
    # the stem-derived name is then one no structure will ever carry.
    default_name = Path(lib).stem
    units = library_unit_names(lib)
    if units:
        if units[0] != default_name:
            console.print(
                f"\n[grey50]The library defines unit [cyan]{units[0]}[/cyan]; the "
                f"filename says [cyan]{default_name}[/cyan]. tLEaP matches the unit "
                f"name and the match is case-sensitive, so [cyan]{units[0]}[/cyan] is "
                f"recorded as the residue name.[/grey50]", highlight=False)
        if len(units) > 1:
            console.print(
                f"[grey50]It declares {len(units)} units "
                f"({', '.join(units)}); the first is used as the residue "
                f"name.[/grey50]", highlight=False)
        default_name = units[0]

    # Imported files can declare atom types; generated ones cannot. See
    # _prompt_imported_atom_types.
    imported_types = _prompt_imported_atom_types(console, processor, frcmod)

    request = _build_request_interactively(
        console, processor,
        category=category, residue_name=default_name,
        frcmod_file=frcmod, lib_file=lib, prep_file=prep,
        atom_types=imported_types or None,
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

    # The next two prompts name two DIFFERENT levels of the library, and the
    # difference is invisible from the prompt text alone: the entry is a
    # directory, the set is a key inside that directory's metadata. Asked back
    # to back with similar wording, they read as the same question twice.
    _print_naming_help(console, family, is_metal)

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
        # GENERATED small molecules / modified AAs declare nothing new -- they
        # are typed with existing GAFF or standard types -- so callers pass
        # None and this stays empty. An IMPORT can supply types read from the
        # frcmod's MASS section, and dropping them here would leave tleap with
        # no addAtomTypes entry for them.
        resolved_atom_types = list(atom_types or [])

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
        residue_meta={
            "description": description,
            "references": references,
            # What must be sourced for these files to load. Recorded here
            # because nothing else will: the Topology Generator enforces
            # declared prerequisites but almost no entry declared any, so the
            # requirement lived only in the depositor's head.
            **_inferred_prerequisites([frcmod_file, *lib_srcs]),
        },
        state_meta=state_meta,
        spin_meta=spin_meta,
        set_meta={"description": description, "version": "1.0"},
        extra_file_srcs=extra_files,
    )


# key, label, what it does. Printed before the prompt because options_map
# does not render -- see the comment at the call site. The meaning is a
# template so the version-bump line can name the set it would actually write
# rather than a placeholder.
_COLLISION_CHOICES = (
    ("1", "Save as a new version",
     "keep the existing set; write this one as {next_name}"),
    ("2", "Overwrite the existing set",
     "replace its files and metadata with these"),
    ("3", "Cancel", "leave the library untouched -- the default"),
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
            # escape(): the message carries the state coordinates in square
            # brackets ("[single_state/default]"), which Rich would read as
            # markup and swallow -- taking with it the one detail that says
            # WHICH leaf collided.
            console.print(f"\n[yellow]{escape(str(exc))}[/yellow]")
            # options_map only feeds the session recorder, and passing it
            # SUPPRESSES the inline choice list, so the options have to be
            # printed here. See _prompt_category, which hit the same trap.
            next_name = f"{request.set_name}_v2"
            for key, label, meaning in _COLLISION_CHOICES:
                # highlight=False: Rich's repr highlighter otherwise recolours
                # bare words inside brackets and parentheses in prose.
                console.print(f"  {key}. [green]{label}[/green] — "
                              f"{meaning.format(next_name=next_name)}",
                              highlight=False)
            choice = prompt_with_context(
                processor, "How would you like to proceed?",
                choices=[k for k, _, _ in _COLLISION_CHOICES], default="3",
                module=_MODULE, description="Resolve library collision",
                options_map={k: label for k, label, _ in _COLLISION_CHOICES},
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
            console.print(f"[red]Could not save to library: {escape(str(exc))}[/red]")
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

def _print_naming_help(console: Console, family: str, is_metal: bool) -> None:
    """Explain entry-vs-set before asking for both."""
    subject = "site" if is_metal else "molecule"
    console.print(Panel(
        f"[bold]Entry name[/bold]   which {subject} — one directory per {subject}\n"
        f"[grey50]specialized_residues/{family}/[/grey50][cyan]<entry>[/cyan][grey50]/[/grey50]\n"
        f"\n"
        f"[bold]Set name[/bold]     whose numbers — one entry can hold several\n"
        f"[grey50]  ...   /metadata.json  ->  forcefield_sets  ->  [/grey50]"
        f"[cyan]<set>[/cyan]\n"
        f"[grey50]A later refit, or a published set for the same {subject}, sits "
        f"beside this one as a second set under the same entry. Exactly one is "
        f"the default.[/grey50]",
        title="Two names, two levels", border_style="cyan",
    ))


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


# Offered category -> (label, what it means, library family it lands in).
_CATEGORY_CHOICES = {
    "1": ("Small molecule / ligand",
          "a cofactor or ligand parameterized as its own residue (FAD, NAD, a drug)",
          "small_molecules"),
    "2": ("Modified amino acid",
          "a residue that stays part of the chain (phosphoserine, a crosslink)",
          "modified_aa"),
    "3": ("Metal site",
          "a metal centre and its coordinating residues (MCPB output)",
          "metal_sites"),
}


# --------------------------------------------------------------------------- #
# atom types declared by an imported frcmod
# --------------------------------------------------------------------------- #

# Masses that identify an element unambiguously. Used to propose the element
# for an addAtomTypes entry; the user confirms before it is stored.
_MASS_TO_ELEMENT = {
    1.008: "H", 4.003: "He", 6.941: "Li", 9.012: "Be", 10.81: "B", 12.01: "C",
    14.01: "N", 16.00: "O", 19.00: "F", 22.99: "Na", 24.305: "Mg",
    26.98: "Al", 28.09: "Si", 30.97: "P", 32.06: "S", 35.45: "Cl",
    39.10: "K", 40.08: "Ca", 47.87: "Ti", 50.94: "V", 52.00: "Cr",
    54.94: "Mn", 55.85: "Fe", 58.93: "Co", 58.69: "Ni", 63.55: "Cu",
    65.38: "Zn", 65.4: "Zn", 79.90: "Br", 95.96: "Mo", 95.94: "Mo",
    126.9: "I", 183.84: "W",
}


def _inferred_prerequisites(paths) -> dict:
    """``{"prerequisites": {...}}`` for a deposit, or ``{}`` when undeterminable.

    Derived from the atom types the files use; see
    ``forcefield_params.prerequisites``. Returns an empty dict rather than an
    empty block so metadata gains nothing when nothing is known.
    """
    try:
        from proprep.forcefield_params.prerequisites import infer_leaprc_groups
        groups = infer_leaprc_groups(list(paths))
    except Exception as exc:  # noqa: BLE001 - a deposit must not fail on this
        logger.debug("Could not infer prerequisites: %s", exc)
        return {}
    if not groups:
        return {}
    return {"prerequisites": {"leaprc_groups": groups}}


def library_unit_names(lib_path) -> List[str]:
    """Unit names an OFF/lib file declares, in file order.

    The unit name is what tLEaP matches a PDB residue against, and the match is
    case-sensitive: a library whose unit is ``gdp`` will not bind a residue
    named ``GDP``. The FILENAME is not authoritative -- ``saveoff`` names the
    entry after the tLEaP variable it was given, which need not resemble the
    file it was written to -- so deriving the residue name from the stem
    records a name no structure may contain.
    """
    names: List[str] = []
    try:
        text = Path(lib_path).read_text(errors="ignore")
    except OSError:
        return names

    in_index = False
    for line in text.splitlines():
        if line.startswith("!!index array str"):
            in_index = True
            continue
        if in_index:
            if line.startswith("!"):
                break
            match = re.match(r'\s*"([^"]+)"', line)
            if match:
                names.append(match.group(1).strip())
    return names


def library_atom_names(lib_path) -> set:
    """Non-hydrogen atom names declared by an OFF/mol2 library unit.

    Hydrogens are excluded because a crystal structure usually lacks them, so
    including them would depress every comparison against a PDB residue.
    """
    names = set()
    try:
        text = Path(lib_path).read_text(errors="ignore")
    except OSError:
        return names

    in_atoms = False
    for line in text.splitlines():
        if line.startswith("!entry.") and ".unit.atoms table" in line:
            in_atoms = True
            continue
        if in_atoms:
            if line.startswith("!"):
                break
            match = re.match(r'\s*"([^"]+)"\s+"([^"]+)"', line)
            if match:
                name = match.group(1).strip()
                if name and not name.upper().startswith("H"):
                    names.add(name.upper())
    return names


def parse_frcmod_mass_types(frcmod_path) -> List[tuple]:
    """``(atom_type, mass)`` for each entry in a frcmod's MASS section.

    The section runs from the ``MASS`` header to the first blank line. An
    frcmod that introduces no new types has an empty one, which is the common
    case for GAFF output.
    """
    types = []
    try:
        in_mass = False
        with open(frcmod_path, errors="ignore") as fh:
            for line in fh:
                stripped = line.strip()
                if not in_mass:
                    if stripped.upper() == "MASS":
                        in_mass = True
                    continue
                if not stripped:
                    break
                parts = stripped.split()
                if len(parts) < 2:
                    continue
                try:
                    types.append((parts[0], float(parts[1])))
                except ValueError:
                    continue
    except OSError:
        return []
    return types


def _known_atom_types() -> set:
    """Atom types the shipped Amber parameter files already define.

    A type already in parm19/parm10/gaff2 needs no addAtomTypes entry, and
    redeclaring one risks changing an existing definition.
    """
    known = set()
    amberhome = os.environ.get("AMBERHOME")
    candidates = []
    if amberhome:
        candidates.append(Path(amberhome) / "dat" / "leap" / "parm")
    candidates.append(Path(sys.prefix) / "dat" / "leap" / "parm")

    for parm_dir in candidates:
        if not parm_dir.is_dir():
            continue
        # The .dat files open with their MASS block; a frcmod labels it. Both
        # are needed: ff19SB's 2C/3C are declared in frcmod.ff19SB, not in
        # parm19.dat, so scanning only the .dat files reports them as new.
        sources = [(parm_dir / n, False)
                   for n in ("parm19.dat", "parm10.dat", "gaff2.dat", "gaff.dat")]
        # Only the frcmods a leaprc actually sources. Scanning all of them is
        # over-broad: frcmod.tumuc, a niche nucleic-acid set, declares YA as a
        # hydrogen, which would vouch for a type ProPrep uses for a metal-bound
        # oxygen and suppress a prompt the user needs.
        for pattern in ("frcmod.ff*", "frcmod.tip*", "frcmod.opc*",
                        "frcmod.spce*", "frcmod.ions*"):
            sources += [(f, True) for f in sorted(parm_dir.glob(pattern))]

        for path, labelled in sources:
            if not path.is_file():
                continue
            try:
                with open(path, errors="ignore") as fh:
                    in_mass = not labelled
                    if not labelled:
                        next(fh, None)      # title line
                    for line in fh:
                        stripped = line.strip()
                        if labelled and not in_mass:
                            if stripped.upper() == "MASS":
                                in_mass = True
                            continue
                        if not stripped:
                            break           # MASS block ends at the first blank
                        parts = stripped.split()
                        if len(parts) >= 2:
                            try:
                                float(parts[1])
                            except ValueError:
                                continue
                            known.add(parts[0])
            except OSError:
                continue
        if known:
            break
    return known


def _prompt_imported_atom_types(console: Console, processor, frcmod_path) -> List[str]:
    """addAtomTypes entries for types an imported frcmod introduces.

    Generated parameters cannot introduce types -- antechamber assigns existing
    GAFF ones -- which is why the post-parameterization path never asks. That
    reasoning does not carry to an IMPORT: those files came from a collaborator
    or a paper and may declare anything, and a type with no addAtomTypes entry
    fails at tleap with nothing in this wizard having mentioned it.

    Nothing is asked when the frcmod declares no new types, which is the common
    case, so this is silent for ordinary GAFF sets.
    """
    declared = parse_frcmod_mass_types(frcmod_path)
    if not declared:
        return []

    known = _known_atom_types()
    if not known:
        # Falling back to "everything is new" is the safe direction -- assuming
        # they are all known would drop the addAtomTypes entry a genuinely new
        # type needs and fail later inside tLEaP. But widening SILENTLY reads as
        # "this frcmod declares unusual types" when it actually means "I could
        # not find Amber's parameter files", which is nearly always ProPrep
        # running outside the environment AmberTools is installed in.
        console.print(Panel(
            "[bold]Could not check these types against Amber's.[/bold]\n"
            "[grey50]parm10/parm19/gaff were not found, so ProPrep cannot tell "
            "which of this frcmod's types are already defined — every one is "
            "listed below, standard types like CT and OS included.\n\n"
            f"Looked in: $AMBERHOME/dat/leap/parm "
            f"(AMBERHOME={os.environ.get('AMBERHOME') or 'unset'})\n"
            f"           {Path(sys.prefix) / 'dat' / 'leap' / 'parm'}\n\n"
            "Declaring a type Amber already defines writes a redundant "
            "addAtomTypes entry. If these should not all be new, cancel and "
            "re-run from the environment AmberTools is installed in.[/grey50]",
            title="Atom types unverified", border_style="yellow"))
    novel = [(t, m) for t, m in declared if t not in known] if known else declared
    if not novel:
        console.print(
            f"[grey50]The frcmod declares {len(declared)} atom type(s), all "
            f"already defined by the standard parameter files.[/grey50]")
        return []

    console.print(
        f"\n[bold]This frcmod introduces {len(novel)} atom type(s)[/bold]")
    console.print(
        "[grey50]tLEaP needs an addAtomTypes entry for each, giving its element "
        "and hybridization. The element is proposed from the mass.[/grey50]")

    entries = []
    for atom_type, mass in novel:
        element = _MASS_TO_ELEMENT.get(round(mass, 3)) or _MASS_TO_ELEMENT.get(round(mass, 2))
        if element is None:
            element = prompt_with_context(
                processor, f"Element for atom type {atom_type} (mass {mass})",
                default="", module=_MODULE,
                description="Element for imported atom type").strip()
        else:
            element = prompt_with_context(
                processor, f"Element for atom type {atom_type} (mass {mass})",
                default=element, module=_MODULE,
                description="Element for imported atom type").strip()
        if not element:
            console.print(f"[yellow]  Skipping {atom_type} — no element given.[/yellow]")
            continue
        hybridization = prompt_with_context(
            processor, f"Hybridization for {atom_type}",
            default="sp3", module=_MODULE,
            description="Hybridization for imported atom type").strip() or "sp3"
        entries.append(f'{{ "{atom_type}" "{element}" "{hybridization}" }}')
        console.print(f"[grey50]  {entries[-1]}[/grey50]")

    return entries


def _prompt_category(console: Console, processor) -> str:
    # options_map only feeds the session recorder, and passing it SUPPRESSES
    # the inline choice list (show_choices defaults to False on the assumption
    # they are printed separately). They were not, so this read as a bare
    # "What kind of parameters are these? (1):" with nothing to go on.
    console.print("\n[bold]What kind of parameters are these?[/bold]")
    for key, (label, meaning, family) in _CATEGORY_CHOICES.items():
        console.print(f"  {key}. [green]{label}[/green] — {meaning}")
        console.print(f"     [grey50]stored under {family}/[/grey50]")

    choice = prompt_with_context(
        processor, "Category",
        choices=list(_CATEGORY_CHOICES), default="1",
        module=_MODULE, description="Parameter category",
        options_map={k: v[0] for k, v in _CATEGORY_CHOICES.items()},
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
    start_dir: Optional[str] = None,
) -> Optional[str]:
    """Browse for an existing file, filtered to ``suffixes``.

    Uses the shared file browser rather than asking for a typed path. Importing
    means pointing at files the user did not produce in this session -- a prior
    project, a collaborator, a published set -- so these are the paths they are
    least likely to have memorised. Typing one still works: the browser's path
    jump accepts a full path, and a directory navigates there.
    """
    from proprep.utils.file_browser import SKIP, default_size_detail, file_browser

    while True:
        selection = file_browser(
            directory=start_dir or os.getcwd(),
            extensions=list(suffixes) if suffixes else None,
            console=console,
            processor=processor,
            label=label,
            entry_detail=default_size_detail,
            allow_path_jump=True,
            optional=not required,
            module=_MODULE,
        )

        if selection is SKIP:
            return None
        if selection is None:
            if not required:
                return None
            console.print("[yellow]A file is required.[/yellow]")
            if not confirm_with_context(processor, "Try again?", default=True,
                                        module=_MODULE,
                                        description="Retry file selection"):
                return None
            continue

        path = Path(str(selection)).expanduser()
        if not path.is_file():
            console.print(f"[yellow]No such file: {path}[/yellow]")
            continue

        # A typed path bypasses the browser's extension filter, so keep the
        # loose suffix check for that route.
        if suffixes and path.suffix.lower() not in suffixes:
            if not confirm_with_context(
                processor,
                f"{path.name} doesn't end in {'/'.join(suffixes)} — use it anyway?",
                default=False, module=_MODULE, description="Confirm unusual suffix",
            ):
                continue

        return str(path.resolve())
