#!/usr/bin/env python3
"""
Comprehensive Redox Site Detector - ProPrep Utility

Complete implementation of the corrected redox site detection workflow:
1. Configuration & Setup
2. Inventory with relational information  
3. Center selection & grouping
4. Iterative site refinement with all 9 search method combinations

ProPrep-style utility class for redox site detection integration.
Author: Claude Code Implementation
"""

import logging
import math
import numpy as np
from typing import Dict, List, Set, Tuple, Optional, Any, Union
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
import copy

# Rich imports for ProPrep-style UI
from rich.console import Console
from proprep.utils.prompts import prompt_with_context, confirm_with_context
from rich.panel import Panel
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn

# BioPython imports
try:
    from Bio.PDB import Structure, Model, Chain, Residue, Atom
    from Bio.PDB.NeighborSearch import NeighborSearch
    from Bio.PDB.PDBParser import PDBParser
    from Bio.PDB.vectors import Vector
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False
    # Fallback type hints
    Structure = Any
    Model = Any
    Chain = Any
    Residue = Any
    Atom = Any
    NeighborSearch = Any
    Vector = Any

# Module-level helpers for routing residue/center highlights through the
# ViewerCoordinator. The 7 sites below all used to spin up their own
# InteractiveStructureViewer (~30 LOC each) and call _launch_viewer, which
# stopped+restarted the underlying server on every "view" press. The
# coordinator's update_annotations path is a live /version-poll bump
# instead — no flicker, no log noise.
def _find_workspace_structure(processor) -> Optional[str]:
    """Return the current PDB path from workspace, or None.

    NOTE: legacy 5-key priority loop. Flagged for migration to
    ``get_priority_pdb_file`` in the legacy-priority-loops backlog; kept
    here verbatim so the viewer-coordinator work doesn't conflate two
    cleanups.
    """
    if processor is None:
        return None
    workspace = processor._get_workspace()
    for key in ['processed_pdb_file', 'local_pdb_file', 'rcsb_pdb_file',
                'alphafold_pdb_file', 'alphafill_pdb_file']:
        structure_file = workspace.get(key)
        if structure_file:
            return structure_file
    return None


def _show_residues_in_viewer(
    processor,
    console: Console,
    ngl_selection: str,
    *,
    label: str,
    focused: bool = True,
) -> bool:
    """Show ``ngl_selection`` highlighted on the current workspace structure.

    Used by user-initiated ``view`` commands — ``force=True`` ensures the
    viewer launches in CLI mode too (where the auto-launch policy is off
    by default).
    """
    if not ngl_selection:
        console.print("[yellow]No residues to view[/yellow]")
        return False

    structure_file = _find_workspace_structure(processor)
    if not structure_file:
        console.print("[yellow]No structure file found in workspace[/yellow]")
        console.print("[grey50]Viewer requires a PDB file path from workspace[/grey50]")
        return False

    from proprep.structure_prep.viewer_coordinator import viewer as _viewer
    _viewer.show_structure(structure_file, force=True)
    _viewer.highlight(
        ngl_selection,
        style="ball+stick",
        color="element",
        label=label,
        focused=focused,
        force=True,
    )
    return True


def _auto_show_palette_highlights(
    processor,
    items: list,
    *,
    replace: bool = True,
    style: str = "ball+stick",
    focused: bool = False,
) -> None:
    """Auto-render palette-coloured highlights at a workflow waypoint.

    ``items`` is a list of ``{"selection": str, "label": str, ...}`` dicts.
    Each item gets a palette colour from its 1-based list position by
    default. Optional per-item overrides:

    - ``"color_index"`` — use a specific palette index instead of the
      list position (lets two reps share one colour, e.g. an inventory
      center and its nearby residues).
    - ``"style"`` — override the default representation style for this
      item only (e.g. ``"line"`` for context rep, ``"ball+stick"`` for
      the center).

    ``replace=True`` clears prior coordinator-managed annotations before
    drawing — used at stage transitions where the set of highlighted
    things changes.

    No ``force=True``: this is auto-update, so it's silent in CLI mode
    when no viewer is running. The coordinator queues the state and it
    will appear when the user explicitly launches the viewer from the menu.
    """
    if not items:
        return
    structure_file = _find_workspace_structure(processor)
    if not structure_file:
        return
    from proprep.structure_prep.viewer_coordinator import viewer as _viewer
    _viewer.show_structure(structure_file)
    if replace:
        _viewer.clear_annotations()
    for i, item in enumerate(items, start=1):
        sel = item.get("selection")
        if not sel:
            continue
        idx = item.get("color_index", i)
        _viewer.highlight(
            sel,
            color=f"palette:{idx}",
            label=item.get("label") or f"item_{i}",
            style=item.get("style", style),
            focused=focused,
        )


def _center_to_ngl(center) -> str:
    """RedoxCenter -> NGL selection string, with altloc when present."""
    sel = f":{center.chain} and {center.resid}"
    if getattr(center, 'altloc', None) and center.altloc.strip():
        sel += f" and %{center.altloc.strip()}"
    return sel


_INVENTORY_PREVIEW_LABEL = "inventory_preview"


def _element_to_ngl(elem: str) -> str:
    """Pick an NGL selection clause that actually selects atoms of ``elem``.

    NGL's ``_X`` element selector parses correctly for single-letter
    element symbols (``_C``, ``_N``, ``_O``, ``_S``) but silently
    breaks for two-letter symbols (``_FE``, ``_SE``, ``_CU``) — the
    selection falls through to "match everything", which blew out the
    inventory coordination halo when the user filtered by FE.

    The clause therefore matches on residue name OR atom name: ``[FE]``
    catches free metal ions modelled as their own residue, ``.FE``
    catches metals embedded in cofactors (HEM iron, FES sulfurs, MSE
    selenium). But the atom-name half collides with protein atom names
    for several metals whose symbol is also a standard PDB atom name —
    ``.CA`` (calcium) matches every alpha-carbon, ``.CD`` (cadmium)
    every Cdelta, ``.CE`` (cerium) the Lys/Met Cepsilon, etc. — so
    filtering "all metals" used to light up most of the protein.

    Guarding the whole clause with the ``hetero`` keyword fixes that:
    free ions and metals embedded in cofactors are all HETATM records,
    while the colliding carbons live in polymer residues, so they fall
    outside ``hetero`` and drop out. Applied uniformly (single- and
    two-letter) since every metal of interest is a hetero atom anyway.
    """
    e = elem.strip().upper()
    if not e:
        return "none"
    return f"(hetero and ([{e}] or .{e}))"


def _resolve_inventory_filter(value: str, *, all_default) -> list[str]:
    """Translate the metals input string to a concrete element list.

    ``ALL``/empty → ``all_default`` (the full METALS set). ``NONE`` →
    empty list. Anything else → the comma-separated tokens.
    """
    v = (value or '').strip().upper()
    if v == 'NONE':
        return []
    if v == 'ALL' or v == '':
        return list(all_default)
    return [t.strip() for t in v.split(',') if t.strip()]


def _show_all_sites_overview(processor, sites) -> None:
    """At redox-detector exit, replace lingering halos with one rep per site.

    The Site Refinement workflow ends with the viewer narrowed to the
    last-processed site (yellow halos on bonded atoms, bond-line rep
    drawn). Carrying that single-site state out of the module is
    misleading — the user just finished detecting *all* sites and is
    moving back to the main menu, where the viewer should reflect the
    full set.

    This helper clears all coordinator-managed annotations and lays
    down one ball+stick rep per detected site, each in a distinct
    palette colour, covering every residue (centers + atoms) the site
    contains. The result is a clean overview: "here are the N redox
    sites you just detected, each in its own colour."

    Empty ``sites`` clears the annotation slate without re-laying anything.
    """
    if not sites:
        from proprep.structure_prep.viewer_coordinator import viewer as _viewer
        _viewer.clear_annotations()
        return

    items = []
    for i, site in enumerate(sites, start=1):
        residues = sorted({
            (c.chain, c.resid) for c in getattr(site, 'centers', [])
        } | {
            (a.chain, a.resid) for a in getattr(site, 'atoms', [])
        })
        if not residues:
            continue
        sel = " or ".join(f"(:{ch} and {rid})" for ch, rid in residues)
        items.append({
            "selection": sel,
            "label": getattr(site, 'site_id', f"site_{i}"),
        })

    if not items:
        from proprep.structure_prep.viewer_coordinator import viewer as _viewer
        _viewer.clear_annotations()
        return

    _auto_show_palette_highlights(processor, items)


def _show_inventory_metals_preview(processor, answers) -> None:
    """At the inventory metals prompt, halo every metal matching the filter.

    Confirms "your metal filter caught these atoms" right when the user
    typed it. We deliberately don't try to preview coordinating non-
    metals: actually identifying coordination partners (a metal-by-metal
    distance search within the bond radius) is precisely what the
    inventory scan does next, and Chunk A already halos those real
    results when the inventory display lands. A preview that re-runs
    the same distance walk just to colour the same atoms a moment
    earlier is duplicate work for no extra signal — and the equivalent
    NGL ``within X of Y`` selection-string syntax doesn't exist anyway.

    Persists through the non-metals (Q2) and search-radius (Q3) prompts
    since those don't trigger an update — the metal halo stays visible
    for the whole atom-filtering section. Replaced when section 2's
    disulfide halo lands. Single label ``_INVENTORY_PREVIEW_LABEL``.
    """
    metal_elems = _resolve_inventory_filter(
        answers.get('metals'), all_default=METALS,
    )

    structure_file = _find_workspace_structure(processor)
    if not structure_file:
        return
    from proprep.structure_prep.viewer_coordinator import viewer as _viewer
    _viewer.show_structure(structure_file)

    if not metal_elems:
        _viewer.unhighlight(_INVENTORY_PREVIEW_LABEL)
        return

    metal_clause = " or ".join(_element_to_ngl(e) for e in metal_elems)
    _viewer.highlight(
        f"({metal_clause})",
        style="halo",
        color="#ffff00",
        label=_INVENTORY_PREVIEW_LABEL,
    )


def _show_inventory_disulfide_preview(processor) -> None:
    """At the disulfide-threshold prompt, halo every CYS SG in the structure.

    The user is calibrating a distance threshold — the natural context
    is "here are all the SG atoms whose pairwise distances will be
    compared against your threshold." Replaces the previous section's
    halo via the shared ``_INVENTORY_PREVIEW_LABEL``.
    """
    structure_file = _find_workspace_structure(processor)
    if not structure_file:
        return
    from proprep.structure_prep.viewer_coordinator import viewer as _viewer
    _viewer.show_structure(structure_file)
    _viewer.highlight(
        "([CYS] and .SG)",
        style="halo",
        color="#ffff00",
        label=_INVENTORY_PREVIEW_LABEL,
    )


def _show_inventory_residue_preview(processor, answers) -> None:
    """At the residue-inclusion prompts, halo the residue classes the user opted into.

    Builds a halo combining whichever of the two residue-inclusion
    answers are currently ``True``:

    - ``include_nonstandard`` → all non-standard residues
      (``not protein and not nucleic and not water and not ion``)
    - ``include_redox_aa`` → ``TYR``, ``TRP``, ``PHE``, ``MET``

    Re-fires after each yes/no prompt so the halo accumulates or
    shrinks as the user toggles. If neither is selected, clears via
    ``unhighlight``.
    """
    structure_file = _find_workspace_structure(processor)
    if not structure_file:
        return
    from proprep.structure_prep.viewer_coordinator import viewer as _viewer
    _viewer.show_structure(structure_file)

    clauses = []
    if answers.get('include_nonstandard'):
        clauses.append("(not protein and not nucleic and not water and not ion)")
    if answers.get('include_redox_aa'):
        redox_clause = " or ".join(
            f"[{r}]" for r in REDOX_ACTIVE_AMINO_ACIDS_VIZ
        )
        clauses.append(f"({redox_clause})")

    if not clauses:
        _viewer.unhighlight(_INVENTORY_PREVIEW_LABEL)
        return

    _viewer.highlight(
        " or ".join(clauses),
        style="halo",
        color="#ffff00",
        label=_INVENTORY_PREVIEW_LABEL,
    )


# Mirror of SimpleCCDQuerier.REDOX_ACTIVE_AMINO_ACIDS but as an ordered
# list for stable NGL selection rendering. CYS is excluded here because
# the inventory's redox-AA toggle covers TYR/TRP/PHE/MET — see Q6 prompt
# text "Show redox-active amino acids (TYR,TRP,PHE,MET)?".
REDOX_ACTIVE_AMINO_ACIDS_VIZ = ['TYR', 'TRP', 'PHE', 'MET']


def _show_boundary_indicator(
    processor,
    site,
    boundary_choice: str,
    *,
    custom_indices=None,
) -> None:
    """Visualize the search-origin geometry implied by the boundary mode.

    Each mode of the boundary prompt (Question 1 in the distance-search
    config; Question 3 in the count-based config) gets a distinct
    overlay so the user can see where the search will radiate from
    *before* they commit to a radius:

    - ``"1"`` (geometric center): a translucent sphere at the centroid
      of all site centers + atoms — the literal point distances will
      be measured from.
    - ``"2"`` (min distance from any site atoms): a yellow halo on
      every center + atom in the site — distances will be measured to
      whichever atom is closest, so the user sees the full envelope.
    - ``"3"`` (min distance from custom atoms): a yellow halo on just
      the user-picked atoms (passed via ``custom_indices``) — same
      semantics as mode 2 but restricted to the chosen subset.

    Re-firing replaces any prior boundary indicator (single stable
    label ``boundary_indicator``).
    """
    structure_file = _find_workspace_structure(processor)
    if not structure_file:
        return

    from proprep.structure_prep.viewer_coordinator import viewer as _viewer
    _viewer.show_structure(structure_file)

    # Always clear first — we may be swapping from a sphere (mode 1) to
    # a halo (mode 2/3) or vice versa, and the two storage backends
    # (annotation_config vs shape_config) don't auto-replace each other.
    _viewer.unhighlight("boundary_indicator")

    if boundary_choice == "1":
        all_coords = (
            [c.coords for c in getattr(site, 'centers', [])]
            + [a.coords for a in getattr(site, 'atoms', [])]
        )
        if not all_coords:
            return
        cx = sum(c[0] for c in all_coords) / len(all_coords)
        cy = sum(c[1] for c in all_coords) / len(all_coords)
        cz = sum(c[2] for c in all_coords) / len(all_coords)
        _viewer.show_sphere(
            (cx, cy, cz),
            radius=1.5,
            label="boundary_indicator",
            color="#ffaa00",
            opacity=0.5,
        )
        return

    if boundary_choice == "2":
        residues = sorted({
            (c.chain, c.resid) for c in getattr(site, 'centers', [])
        } | {
            (a.chain, a.resid) for a in getattr(site, 'atoms', [])
        })
        if not residues:
            return
        sel = " or ".join(f"(:{ch} and {rid})" for ch, rid in residues)
        _viewer.highlight(
            sel,
            style="halo",
            color="#ffff00",
            label="boundary_indicator",
        )
        return

    if boundary_choice == "3" and custom_indices:
        atoms = getattr(site, 'atoms', [])
        by_residue: dict[tuple[str, int], set[str]] = {}
        for idx in custom_indices:
            if not (0 <= idx < len(atoms)):
                continue
            atom = atoms[idx]
            key = (atom.chain, atom.resid)
            by_residue.setdefault(key, set()).add(atom.atom_name)
        clauses = []
        for (ch, rid), names in by_residue.items():
            sorted_names = sorted(n for n in names if n)
            if not sorted_names:
                continue
            names_clause = " or ".join(f".{n}" for n in sorted_names)
            clauses.append(f"(:{ch} and {rid} and ({names_clause}))")
        if not clauses:
            return
        _viewer.highlight(
            " or ".join(clauses),
            style="halo",
            color="#ffff00",
            label="boundary_indicator",
        )


def _atom_to_ngl_selection(info) -> str:
    """coord_to_pdb entry -> NGL selection that picks exactly one atom."""
    sel = f":{info.get('chain', '')} and {info.get('resid', 0)} and .{info.get('atom_name', '')}"
    altloc = (info.get('altloc') or '').strip()
    if altloc:
        sel += f" and %{altloc}"
    return sel


def _highlight_site_bonds(processor, site) -> None:
    """Draw a yellow halo + bond-line overlay for every bond in ``site``.

    Two overlays at one moment:

    - ``bond_def_bonded`` (halo): translucent spacefill around every atom
      that participates in any defined bond, across the whole site.
    - ``bond_def_lines`` (bond): an NGL distance representation drawing
      a line per bond between its two endpoint atoms. Surfaces the
      coordinate / covalent bonds the detector defined for the site —
      these aren't in the structure's CONECT records so without this
      rep they're invisible in the viewer.

    Iterates ``site.bonds`` and resolves each endpoint via ``coord_to_pdb``
    to find the (chain, resid, atom_name) of every bonded atom. So
    creating an FE↔NE2 bond lights up the heme iron and the HIS
    coordinating nitrogen at once with a halo, plus draws a yellow line
    between them.

    Called from interactive bond definition (after each new bond) and
    from ``apply_template_to_site`` (once, after the template's bonds
    are stamped). If no bonds are defined yet, both reps are removed
    via ``unhighlight`` so stale overlays don't linger.
    """
    structure_file = _find_workspace_structure(processor)
    if not structure_file:
        return

    from proprep.structure_prep.viewer_coordinator import viewer as _viewer
    _viewer.show_structure(structure_file)

    # Group bonded atom names by (chain, resid) so the NGL selection has
    # one residue clause per residue with all its bonded atom names OR-ed
    # together. Skip endpoints whose coords aren't in coord_to_pdb (should
    # never happen, but guard so a corrupt site doesn't wedge the viewer).
    #
    # Every bond gets a yellow distance line, even ones that the existing
    # ball+stick rep already draws (intra-residue CA-CB etc.). The line
    # may be hidden behind the underlying cylinder in those cases — the
    # halo on both endpoints still indicates the atoms are part of a
    # defined bond, and not filtering keeps "what the user defined"
    # exactly equal to "what the bond rep represents".
    by_residue: dict[tuple[str, int], set[str]] = {}
    bond_pairs: list[tuple[str, str]] = []
    for bond in site.bonds:
        info_a = site.coord_to_pdb.get(bond.atom1_coords)
        info_b = site.coord_to_pdb.get(bond.atom2_coords)
        if not info_a or not info_b:
            continue
        for info in (info_a, info_b):
            key = (info.get('chain', ''), info.get('resid', 0))
            by_residue.setdefault(key, set()).add(info.get('atom_name', ''))
        bond_pairs.append(
            (_atom_to_ngl_selection(info_a), _atom_to_ngl_selection(info_b))
        )

    if not by_residue:
        _viewer.unhighlight("bond_def_bonded")
        _viewer.unhighlight("bond_def_lines")
        return

    residue_clauses = []
    for (ch, rid), atom_names in by_residue.items():
        names = sorted(n for n in atom_names if n)
        if not names:
            continue
        names_clause = " or ".join(f".{n}" for n in names)
        residue_clauses.append(f"(:{ch} and {rid} and ({names_clause}))")

    if not residue_clauses:
        _viewer.unhighlight("bond_def_bonded")
        _viewer.unhighlight("bond_def_lines")
        return

    selection = " or ".join(residue_clauses)
    _viewer.highlight(
        selection,
        style="halo",
        color="#ffff00",
        label="bond_def_bonded",
    )
    _viewer.show_bonds(
        bond_pairs,
        label="bond_def_lines",
        color="#ffff00",
        show_labels=False,
    )


def _narrow_viewer_to_site(processor, site) -> None:
    """Narrow the docked viewer to a single ``RedoxSite``'s contents.

    Called from both refinement entry points — ``refine_site_interactively``
    (the first site of each template type, where the user configures the
    template interactively) and ``apply_template_to_site`` (the auto-
    applied subsequent sites). Both should give the same visual cue
    of "this is the site we're working on now". Palette index 1 is used
    as a fixed "current site" marker (no group color, since only one
    site is visible at a time).

    The protein cartoon stays visible for context; the camera auto-
    centers on the site selection so the visual change is obvious even
    when the un-narrowed centers were close to the kept one.
    """
    site_center_residues = sorted({
        (c.chain, c.resid) for c in getattr(site, 'centers', [])
    })
    site_other_residues = sorted({
        (a.chain, a.resid) for a in getattr(site, 'atoms', [])
        if (a.chain, a.resid) not in site_center_residues
    })
    items = []
    site_selection_parts = []
    if site_center_residues:
        sel = " or ".join(f":{ch} and {rid}" for ch, rid in site_center_residues)
        items.append({
            "selection": sel, "label": "current_site_center", "color_index": 1,
        })
        site_selection_parts.append(sel)
    if site_other_residues:
        sel = " or ".join(f":{ch} and {rid}" for ch, rid in site_other_residues)
        items.append({
            "selection": sel, "label": "current_site_members", "color_index": 1,
        })
        site_selection_parts.append(sel)
    _auto_show_palette_highlights(processor, items)
    # Auto-zoom on the site so the camera move makes the narrowing
    # visually obvious. Without this, the camera stays where it was and
    # narrowing is barely perceptible if both sites were close.
    if site_selection_parts:
        from proprep.structure_prep.viewer_coordinator import viewer as _viewer
        _viewer.focus_on(" or ".join(site_selection_parts))


from proprep.utils.prompts import (
    prompt_with_context,
    confirm_with_context,
    prompt_float_with_retry,
    prompt_int_with_retry,
    NavigationException,
)

logger = logging.getLogger(__name__)


def with_navigation_loop(question_sequence):
    """
    Decorator/wrapper to add back navigation to a sequence of questions.

    Args:
        question_sequence: List of tuples (prompt_func, args_dict) where:
            - prompt_func is one of the validated prompt functions
            - args_dict contains the arguments to pass to that function

    Returns:
        List of answers in order, or None if user exits

    Example:
        answers = with_navigation_loop([
            (prompt_float_with_retry, {'prompt': 'Enter radius', 'default': 5.0}),
            (prompt_int_with_retry, {'prompt': 'Enter count', 'default': 10}),
            (confirm_with_context, {'prompt': 'Confirm?', 'default': True}),
        ])
    """
    answers = []
    current_step = 0

    while current_step < len(question_sequence):
        prompt_func, args = question_sequence[current_step]

        # Ensure allow_back is set for all but the first question
        if current_step > 0:
            args = dict(args)  # Copy to avoid modifying original
            args['allow_back'] = True
        else:
            # First question - no back navigation
            args = dict(args)
            args['allow_back'] = False

        try:
            # Ask the question
            answer = prompt_func(**args)

            # Store answer and move forward
            if current_step < len(answers):
                answers[current_step] = answer  # Updating existing answer
            else:
                answers.append(answer)

            current_step += 1

        except NavigationException:
            # User wants to go back
            if current_step > 0:
                current_step -= 1
                # Remove the answer we're going back from
                if current_step < len(answers):
                    answers = answers[:current_step]
            # If current_step is 0, we can't go back further (shouldn't happen due to allow_back=False)

    return answers


def prompt_with_navigation(processor, prompt, choices=None, default=None,
                           module=None, description=None, options_map=None,
                           allow_back=True):
    """
    Enhanced prompt that allows 'back' navigation.

    Args:
        processor: The processor object (has session manager)
        prompt: The question to ask
        choices: Valid answers
        default: Default answer if user just presses Enter
        module: Name of this part of ProPrep
        description: What this question is asking
        options_map: Dictionary mapping choices to descriptions
        allow_back: If True, user can type 'back' or 'b' to go back

    Returns:
        The user's answer

    Raises:
        NavigationException: If user types 'back' and allow_back=True
    """
    # Add back to choices if allowed
    extended_choices = choices
    if allow_back and choices:
        extended_choices = list(choices) + ['back', 'b']

    # Add hint about back navigation
    enhanced_prompt = prompt
    if allow_back:
        enhanced_prompt = f"{prompt}\n[grey50]  (Type 'back' or 'b' to return to previous question)[/grey50]"

    response = prompt_with_context(
        processor=processor,
        prompt=enhanced_prompt,
        choices=extended_choices if extended_choices != choices else choices,
        default=default,
        module=module,
        description=description,
        options_map=options_map
    )

    # Check for back navigation
    if allow_back and response and response.lower() in ['back', 'b']:
        raise NavigationException("User requested to go back")

    return response


# Lazy import of transformer registry to avoid circular import
# (transformers import CenterType from this module)
_redox_transformer_registry = None
_transformer_registry_import_attempted = False

def _get_transformer_registry():
    """Lazy import of transformer registry to avoid circular import"""
    global _redox_transformer_registry, _transformer_registry_import_attempted

    if not _transformer_registry_import_attempted:
        _transformer_registry_import_attempted = True
        try:
            from proprep.redoxsite_prep.transformation import redox_transformer_registry
            _redox_transformer_registry = redox_transformer_registry
        except ImportError as e:
            logger.debug(f"Transformer registry not available: {e}")
            _redox_transformer_registry = None

    # Register user-developed transformers (auto-emitted .py + JSON specs under
    # ~/.proprep/transformers) so the detector lists them alongside the built-ins
    # when the user assigns site-type/template names. Called on every access:
    # load_user_transformers is idempotent (skips already-loaded files) but still
    # picks up transformers created earlier this session (e.g. via the Redox Site
    # Preparer's transformer creator).
    if _redox_transformer_registry is not None:
        try:
            from proprep.redoxsite_prep.transformation.auto_rename import (
                load_user_transformers,
            )
            load_user_transformers()
        except Exception as e:  # pragma: no cover - defensive
            logger.debug(f"Could not load user transformers into detector: {e}")

    return _redox_transformer_registry

# ===== COMPREHENSIVE METAL CLASSIFICATION =====

# Complete periodic table metals classification
METALS = {
    # Group 1: Alkali metals
    'LI', 'NA', 'K', 'RB', 'CS', 'FR',
    
    # Group 2: Alkaline earth metals  
    'BE', 'MG', 'CA', 'SR', 'BA', 'RA',
    
    # Group 3-12: Transition metals
    'SC', 'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN',
    'Y', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD',
    'LA', 'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG',
    'AC', 'RF', 'DB', 'SG', 'BH', 'HS', 'MT', 'DS', 'RG', 'CN',
    
    # Lanthanides
    'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU',
    
    # Actinides
    'TH', 'PA', 'U', 'NP', 'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM', 'MD', 'NO', 'LR',
    
    # Post-transition metals
    'AL', 'GA', 'IN', 'SN', 'TL', 'PB', 'BI', 'PO',
    
    # Metalloids that often behave as metals in coordination
    'SB', 'TE'
}

# ===== CENTER TYPE CLASSIFICATION =====

class CenterType(Enum):
    """Classification of redox center types.

    - METAL_ION: Isolated metal ions (single-atom, e.g., Zn coordinated to protein)
    - ORGANOMETALLIC_COFACTOR: Multi-atom cofactors with embedded metal (heme, Fe-S clusters)
    - ORGANIC_COFACTOR: Purely organic cofactors without metal (flavin, quinones)
    - REDOX_AMINO_ACID: Standard redox-active amino acids (TYR, TRP, PHE, MET, CYS)
    """
    METAL_ION = "metal_ion"
    ORGANOMETALLIC_COFACTOR = "organometallic_cofactor"
    ORGANIC_COFACTOR = "organic_cofactor"
    REDOX_AMINO_ACID = "redox_amino_acid"

# ===== SEARCH METHOD ENUMS =====

class SearchConstraint(Enum):
    """Search constraint types"""
    DISTANCE_CUTOFF = "distance_cutoff"
    COUNT_CUTOFF = "count_cutoff"

class DistanceMethod(Enum):
    """Distance-based search methods"""
    FIXED = "fixed"
    ADAPTIVE = "adaptive"

class BoundaryDefinition(Enum):
    """Site boundary definition methods"""
    GEOMETRIC_CENTER = "geometric_center"
    MIN_DISTANCE_ALL = "min_distance_all"
    MIN_DISTANCE_CUSTOM = "min_distance_custom"

# ===== CORE DATA STRUCTURES =====

@dataclass
class RedoxCenter:
    """A redox-active center identified from structure analysis"""
    # PDB identification
    chain: str
    resname: str
    resid: int
    coords: Tuple[float, float, float]  # Immutable coordinates (permanent identifier)
    center_type: CenterType
    atom_name: Optional[str] = None  # None if whole residue is the center
    insertion_code: str = ""
    altloc: str = ""  # Alternate location indicator (e.g., 'A', 'B', or '' for none)
    element: Optional[str] = None
    properties: Dict[str, Any] = field(default_factory=dict)

def sync_redox_sites_from_pdb(pdb_file, sites) -> bool:
    """Update RedoxSite objects to match a PDB via coordinate-based matching.

    After a residue/atom RENAME that preserves atomic coordinates (e.g. the
    modified-AA Route B integration renaming CYF->CF1 / RBF->RF1 and their
    atoms to library names), the workspace RedoxSite objects still carry the
    OLD residue names, atom names, and ``coord_to_pdb`` mapping. Downstream
    consumers (PDB Filter, Topology Generator) then emit stale tLEaP `bond`
    commands (`bond mol.44.SG ...` when the renamed lib atom is now `S1`),
    which tLEaP rejects.

    Because coordinates are invariant under a pure rename, we re-key every
    RedoxSite field off the atom coordinates: parse the PDB into a
    coord -> (chain, resname, resid, atom_name, insertion_code) lookup and
    stamp the fresh identifiers onto each site's atoms, centers,
    ``coord_to_pdb``, bonds, and ``residue_groups``.

    This is the module-level, reusable form of
    ``RedoxTransformationManager._update_redox_sites_from_existing_pdb`` so
    the modified-AA parameterizer can sync without a full manager instance.

    Args:
        pdb_file: Path to the renamed/prepared PDB.
        sites: List of RedoxSite objects to update in place.

    Returns:
        True if any atom was remapped, False otherwise.
    """
    pdb_path = Path(pdb_file)
    if not pdb_path.exists():
        logger.warning(f"PDB file not found: {pdb_file}")
        return False

    try:
        with open(pdb_path, 'r') as f:
            pdb_lines = f.readlines()
    except Exception as e:
        logger.error(f"Error reading PDB file {pdb_file}: {e}")
        return False

    # coord (3-decimal-rounded) -> atom identity, matching RedoxSiteAtom coords
    coord_lookup = {}
    for line in pdb_lines:
        if line.startswith(("ATOM", "HETATM")):
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                coord_lookup[(round(x, 3), round(y, 3), round(z, 3))] = {
                    'chain': line[21],
                    'resname': line[17:20].strip(),
                    'resid': int(line[22:26]),
                    'atom_name': line[12:16].strip(),
                    'insertion_code': line[26].strip() if len(line) > 26 else '',
                }
            except (ValueError, IndexError):
                continue

    if not coord_lookup:
        logger.warning(f"No atoms found in PDB file: {pdb_file}")
        return False

    if not sites:
        return False

    def find_matching_atom(coords):
        key = (round(float(coords[0]), 3), round(float(coords[1]), 3), round(float(coords[2]), 3))
        return coord_lookup.get(key)

    updates_made = 0
    for site in sites:
        for atom in site.atoms:
            new_info = find_matching_atom(atom.coords)
            if new_info:
                atom.chain = new_info['chain']
                atom.resname = new_info['resname']
                atom.resid = new_info['resid']
                atom.atom_name = new_info['atom_name']
                atom.insertion_code = new_info.get('insertion_code', '')
                updates_made += 1

        for center in site.centers:
            new_info = find_matching_atom(center.coords)
            if new_info:
                center.chain = new_info['chain']
                center.resname = new_info['resname']
                center.resid = new_info['resid']
                center.atom_name = new_info['atom_name']
                center.insertion_code = new_info.get('insertion_code', '')

        new_coord_to_pdb = {}
        for coords, pdb_info in site.coord_to_pdb.items():
            new_info = find_matching_atom(coords)
            if new_info:
                updated_info = copy.copy(pdb_info)
                updated_info.update(new_info)
                new_coord_to_pdb[coords] = updated_info
            else:
                new_coord_to_pdb[coords] = pdb_info
        site.coord_to_pdb = new_coord_to_pdb

        for bond in site.bonds:
            if hasattr(bond, 'atom1_coords'):
                new_info = find_matching_atom(bond.atom1_coords)
                if new_info and hasattr(bond, 'atom1_residue_info'):
                    bond.atom1_residue_info.update(new_info)
            if hasattr(bond, 'atom2_coords'):
                new_info = find_matching_atom(bond.atom2_coords)
                if new_info and hasattr(bond, 'atom2_residue_info'):
                    bond.atom2_residue_info.update(new_info)

        new_residue_groups = {}
        for (chain, resid, icode), coords_list in site.residue_groups.items():
            if coords_list:
                new_info = find_matching_atom(coords_list[0])
                if new_info:
                    new_key = (new_info['chain'], new_info['resid'], new_info.get('insertion_code', ''))
                    new_residue_groups.setdefault(new_key, []).extend(coords_list)
                else:
                    new_residue_groups[(chain, resid, icode)] = coords_list
            else:
                new_residue_groups[(chain, resid, icode)] = coords_list
        site.residue_groups = new_residue_groups

    logger.debug(f"Updated {updates_made} atoms in {len(sites)} RedoxSites from {pdb_file}")
    return updates_made > 0


@dataclass
class RedoxSiteAtom:
    """Individual atom in a redox site with PDB metadata + immutable coordinates"""
    # Current PDB identification (changes during transformations)
    chain: str
    resname: str  
    resid: int
    atom_name: str
    coords: Tuple[float, float, float]  # NEVER changes (permanent identifier)
    element: str
    altloc: str = ""
    insertion_code: str = ""
    occupancy: Optional[float] = None
    bfactor: Optional[float] = None
    properties: Dict[str, Any] = field(default_factory=dict)  # For serial_number, etc.

@dataclass
class RedoxSiteBond:
    """Bond between two atoms identified by coordinates"""
    atom1_coords: Tuple[float, float, float]
    atom2_coords: Tuple[float, float, float]
    bond_type: str  # "intraresidue" or "interresidue"
    chemical_type: str  # "coordinate", "covalent", "disulfide", "unknown"
    distance: float
    atom1_element: str
    atom2_element: str
    atom1_residue_info: Dict[str, Any] = field(default_factory=dict)
    atom2_residue_info: Dict[str, Any] = field(default_factory=dict)
    # How the contact is REALIZED in the force field, independent of what it IS
    # chemically (chemical_type). "bonded" = MCPB bonded model (metal-ligand
    # bond term, Y-type renaming, RESP repartition, tLEaP bond). "restrained" =
    # kept as a nonbonded residue (e.g. TIP3P water), held near the metal by an
    # MD distance restraint; the contact still keeps the coordinating residue in
    # the QM models for correct electronics but emits no bonded term. Default
    # "bonded" preserves all existing behavior.
    treatment: str = "bonded"

# ===== BOND CLASSIFICATION FUNCTION =====

def classify_bond_types(atom1_element: str, atom2_element: str, 
                       distance: float, atom1_residue: str, atom2_residue: str,
                       atom1_resid: int, atom2_resid: int, 
                       atom1_chain: str, atom2_chain: str,
                       atom1_name: str = None, atom2_name: str = None) -> Tuple[str, str]:
    """
    Classify bond type and residue relationship.
    
    Returns:
        Tuple of (bond_type, chemical_type)
        bond_type: "intraresidue" or "interresidue"  
        chemical_type: "covalent", "coordinate", "metal-metal", "disulfide"
    """
    elem1, elem2 = atom1_element.upper(), atom2_element.upper()
    is_metal1 = elem1 in METALS
    is_metal2 = elem2 in METALS
    
    # Determine bond_type (intra vs inter-residue)
    same_residue = (atom1_chain == atom2_chain and 
                   atom1_resid == atom2_resid and 
                   atom1_residue == atom2_residue)
    bond_type = "intraresidue" if same_residue else "interresidue"
    
    # Check for disulfide bonds first (SG-SG between CYS residues)
    if (atom1_name == 'SG' and atom1_residue == 'CYS' and
        atom2_name == 'SG' and atom2_residue == 'CYS'):
        chemical_type = "disulfide"
    # Determine chemical_type based on metal content
    elif not is_metal1 and not is_metal2:
        # 0 metal atoms = covalent bond
        chemical_type = "covalent"
    elif (is_metal1 and not is_metal2) or (is_metal2 and not is_metal1):
        # 1 metal atom and 1 non-metal atom = coordinate bond
        chemical_type = "coordinate"
    elif is_metal1 and is_metal2:
        # 2 metal atoms = metal-metal bond
        chemical_type = "metal-metal"
    
    return bond_type, chemical_type


def resolve_bond_residue_token(token, residue_list, on_error=None):
    """Map a bond-definition token to a 0-based index into ``residue_list``.

    ``residue_list`` is the ordered list of ``(chain, resname, resid,
    insertion_code)`` keys shown in the bond-definition table. A token may be a
    residue ID — bare resid (``"202"``), resname+resid (``"MN202"``, ``"E4Z201"``)
    or chain:resid (``"A:202"``) — so the user can read IDs straight off the
    viewer. If it matches no residue ID it falls back to the 1-based table row
    number, so older session logs that recorded row numbers still replay
    (real resids rarely collide with the small row numbers). Residue-ID matches
    take precedence over row numbers.

    Returns the int index, or ``None`` (calling ``on_error(msg)`` if provided).
    """
    import re as _re

    def _fail(msg):
        if on_error:
            on_error(msg)
        return None

    t = (token or "").strip()
    if not t:
        return None
    chain_filter = None
    body = t
    if ':' in body:
        cpart, body = body.split(':', 1)
        chain_filter = cpart.strip().upper() or None
    m = _re.fullmatch(r'([A-Za-z0-9]*?)(\d+)', body.strip())
    if not m:
        return _fail(f"Could not parse residue '{token}'. Try 202, MN202, or A:202.")
    name_part = m.group(1).upper()
    num = int(m.group(2))
    matches = [
        idx for idx, (chain, resname, resid, icode) in enumerate(residue_list)
        if resid == num
        and (not name_part or (resname or "").upper() == name_part)
        and (not chain_filter or (chain or "").upper() == chain_filter)
    ]
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        opts = ", ".join(
            f"{residue_list[i][1]} {residue_list[i][0]}:{residue_list[i][2]}"
            for i in matches
        )
        return _fail(f"Residue '{token}' is ambiguous ({opts}); qualify with chain, e.g. A:{num}.")
    # No residue-ID match: fall back to the 1-based table row index.
    if not name_part and chain_filter is None and 1 <= num <= len(residue_list):
        return num - 1
    example = f"{residue_list[0][1]}{residue_list[0][2]}" if residue_list else "MN202"
    return _fail(f"No residue matches '{token}'. Use a residue ID (e.g. {example}) or row 1-{len(residue_list)}.")


# ===== REDOX SITE CLASS =====

class RedoxSite:
    """Complete redox site definition with coordinate tracking for transformations"""
    
    def __init__(self, site_id: str, structure_id: str):
        self.site_id = site_id
        self.structure_id = structure_id
        
        # Site composition
        self.centers: List[RedoxCenter] = []  # The redox-active centers
        self.atoms: List[RedoxSiteAtom] = []  # All atoms in site
        self.bonds: List[RedoxSiteBond] = []  # All bonds detected
        
        # Coordinate mapping (key for transformations)
        self.coord_to_pdb: Dict[Tuple[float, float, float], Dict[str, Any]] = {}
        
        # Site organization
        self.residue_groups: Dict[Tuple[str, int, str], List[Tuple[float, float, float]]] = {}
        
        # Detection metadata
        self.detection_method: str = ""
        self.detection_parameters: Dict[str, Any] = {}
        self.boundary_definition: str = ""
        self.search_history: List[Dict[str, Any]] = []  # Track iterative searches
        self.site_type: str = ""  # User-specified site type from template categorization
    
    def add_center(self, center: RedoxCenter):
        """Add a redox center to the site"""
        self.centers.append(center)
        
        # Add to coordinate mapping
        self.coord_to_pdb[center.coords] = {
            'chain': center.chain,
            'resname': center.resname,
            'resid': center.resid,
            'atom_name': center.atom_name,
            'element': center.element,
            'insertion_code': center.insertion_code,
            'center_type': center.center_type.value
        }
    
    def add_atom(self, atom: RedoxSiteAtom):
        """Add an atom to the site"""
        self.atoms.append(atom)
        
        # Add to coordinate mapping
        self.coord_to_pdb[atom.coords] = {
            'chain': atom.chain,
            'resname': atom.resname,
            'resid': atom.resid,
            'atom_name': atom.atom_name,
            'element': atom.element,
            'altloc': atom.altloc,
            'insertion_code': atom.insertion_code,
            'occupancy': atom.occupancy,
            'bfactor': atom.bfactor
        }
        
        # Update residue groupings
        res_key = (atom.chain, atom.resid, atom.insertion_code)
        if res_key not in self.residue_groups:
            self.residue_groups[res_key] = []
        self.residue_groups[res_key].append(atom.coords)
    
    def add_bond_with_classification(self, atom1_coords: Tuple[float, float, float],
                                   atom2_coords: Tuple[float, float, float],
                                   distance: float,
                                   treatment: str = "bonded"):
        """Add bond with automatic type classification.

        treatment: "bonded" (default, MCPB bonded model) or "restrained"
        (keep as nonbonded residue held by an MD distance restraint).
        """
        
        # Get atom information from site
        atom1_info = self.coord_to_pdb.get(atom1_coords)
        atom2_info = self.coord_to_pdb.get(atom2_coords)
        
        if not atom1_info or not atom2_info:
            logger.warning(f"Cannot classify bond - missing atom info")
            return
        
        # Classify bond
        bond_type, chemical_type = classify_bond_types(
            atom1_info['element'], atom2_info['element'],
            distance,
            atom1_info['resname'], atom2_info['resname'],
            atom1_info['resid'], atom2_info['resid'],
            atom1_info['chain'], atom2_info['chain'],
            atom1_info.get('atom_name'), atom2_info.get('atom_name')
        )
        
        bond = RedoxSiteBond(
            atom1_coords=atom1_coords,
            atom2_coords=atom2_coords,
            bond_type=bond_type,
            chemical_type=chemical_type,
            distance=distance,
            atom1_element=atom1_info['element'],
            atom2_element=atom2_info['element'],
            atom1_residue_info=atom1_info,
            atom2_residue_info=atom2_info,
            treatment=treatment
        )

        self.bonds.append(bond)
    
    def update_atom_metadata(self, coord_to_new_metadata: Dict[Tuple[float, float, float], Dict[str, Any]]):
        """Update PDB metadata after transformations - coordinates stay same"""
        for coords, new_meta in coord_to_new_metadata.items():
            # DEBUG: Check metadata
            if 'chain' not in new_meta:
                print(f"ERROR: Missing 'chain' in metadata for coords {coords}")
                print(f"Metadata keys: {new_meta.keys()}")
                print(f"Metadata: {new_meta}")
                raise ValueError(f"Missing 'chain' in metadata for coords {coords}")
            # Update coord_to_pdb mapping
            if coords in self.coord_to_pdb:
                self.coord_to_pdb[coords].update(new_meta)
            
            # Update atom objects
            for atom in self.atoms:
                if atom.coords == coords:
                    atom.chain = new_meta.get('chain', atom.chain)
                    atom.resname = new_meta.get('resname', atom.resname)
                    atom.resid = new_meta.get('resid', atom.resid)
                    atom.atom_name = new_meta.get('atom_name', atom.atom_name)
                    break
            
            # Update center objects
            for center in self.centers:
                if center.coords == coords:
                    center.chain = new_meta.get('chain', center.chain)
                    center.resname = new_meta.get('resname', center.resname)
                    center.resid = new_meta.get('resid', center.resid)
                    center.atom_name = new_meta.get('atom_name', center.atom_name)
                    break
            
            # Update residue groupings
            old_groups = list(self.residue_groups.keys())
            for old_group in old_groups:
                if coords in self.residue_groups[old_group]:
                    # Remove from old group
                    self.residue_groups[old_group].remove(coords)
                    if not self.residue_groups[old_group]:
                        del self.residue_groups[old_group]
                    
                    # Add to new group - ensure resid is an integer
                    new_resid = new_meta.get('resid')
                    if isinstance(new_resid, str):
                        new_resid = int(new_resid)
                    new_group = (new_meta.get('chain'), new_resid, new_meta.get('insertion_code', ''))
                    if new_group not in self.residue_groups:
                        self.residue_groups[new_group] = []
                    self.residue_groups[new_group].append(coords)
                    break
    
    def get_atoms_by_residue(self, chain: str, resid: int, insertion_code: str = "") -> List[RedoxSiteAtom]:
        """Get atoms by current residue identification"""
        # For "no insertion" cases, check both empty string and space variants
        if insertion_code in ('', ' '):
            # Try both possible "no insertion" representations
            primary_key = (chain, resid, insertion_code)
            alternate_key = (chain, resid, ' ' if insertion_code == '' else '')
            
            # Collect coordinates from both possible keys
            coord_list = []
            if primary_key in self.residue_groups:
                coord_list.extend(self.residue_groups[primary_key])
            if alternate_key in self.residue_groups:
                coord_list.extend(self.residue_groups[alternate_key])
            
            if coord_list:
                return [atom for atom in self.atoms if atom.coords in coord_list]
            else:
                return []
        else:
            # For actual insertions, use exact match
            key = (chain, resid, insertion_code)
            
            if key in self.residue_groups:
                coord_list = self.residue_groups[key]
                return [atom for atom in self.atoms if atom.coords in coord_list]
            else:
                return []
    
    def get_current_pdb_info(self, coords: Tuple[float, float, float]) -> Optional[Dict[str, Any]]:
        """Get current PDB info for an atom by coordinates"""
        return self.coord_to_pdb.get(coords)
    
    def get_site_boundary_coords(self, boundary_def: BoundaryDefinition, 
                               custom_coords: List[Tuple[float, float, float]] = None) -> List[Tuple[float, float, float]]:
        """Get boundary coordinates based on boundary definition"""
        if boundary_def == BoundaryDefinition.GEOMETRIC_CENTER:
            # Calculate geometric center
            all_coords = [center.coords for center in self.centers] + [atom.coords for atom in self.atoms]
            if not all_coords:
                return []
            center = tuple(np.mean(all_coords, axis=0))
            return [center]
        elif boundary_def == BoundaryDefinition.MIN_DISTANCE_ALL:
            # Use all atoms in site
            return [center.coords for center in self.centers] + [atom.coords for atom in self.atoms]
        elif boundary_def == BoundaryDefinition.MIN_DISTANCE_CUSTOM:
            # Use custom selected coordinates
            return custom_coords if custom_coords else []
        else:
            return []

# ===== CCD QUERIER =====

class SimpleCCDQuerier:
    """Simplified CCD querier for identifying non-standard residues"""
    
    # Standard amino acids (non-redox-active except where noted)
    STANDARD_AMINO_ACIDS = {
        'ALA', 'ARG', 'ASN', 'ASP', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU',
        'LYS', 'PRO', 'SER', 'THR', 'VAL'
    }
    
    # Redox-active standard amino acids
    REDOX_ACTIVE_AMINO_ACIDS = {
        'PHE', 'TYR', 'TRP',  # Aromatic - can participate in electron transfer
        'MET',                # Sulfur - can be oxidized
        'CYS'                 # Cysteine - can form multi-center redox sites, metal coordination
    }
    
    # Common water/ion residues to ignore
    COMMON_SOLVENTS = {
        'HOH', 'WAT', 'H2O', 'SO4', 'PO4', 'CL', 'NA', 'MG', 'CA', 'K'
    }
    
    def __init__(self, use_cache: bool = True):
        self.use_cache = use_cache
        self._cache = {}
    
    def is_non_standard_residue(self, residue_name: str) -> bool:
        """Check if residue is non-standard (not a standard amino acid)"""
        residue_name = residue_name.strip().upper()
        
        return (residue_name not in self.STANDARD_AMINO_ACIDS and 
                residue_name not in self.REDOX_ACTIVE_AMINO_ACIDS and
                residue_name not in self.COMMON_SOLVENTS)
    
    def is_redox_active_amino_acid(self, residue_name: str) -> bool:
        """Check if residue is a redox-active amino acid"""
        return residue_name.strip().upper() in self.REDOX_ACTIVE_AMINO_ACIDS
    
    def get_basic_residue_info(self, residue_name: str, residue_obj=None) -> Dict[str, Any]:
        """Get basic information about a residue from structure analysis"""
        residue_name = residue_name.strip().upper()
        
        if residue_name in self._cache:
            return self._cache[residue_name]
        
        info = {
            'residue_name': residue_name,
            'is_standard_amino_acid': residue_name in self.STANDARD_AMINO_ACIDS,
            'is_redox_active_amino_acid': residue_name in self.REDOX_ACTIVE_AMINO_ACIDS,
            'is_solvent': residue_name in self.COMMON_SOLVENTS,
            'is_non_standard': self.is_non_standard_residue(residue_name)
        }
        
        # Add basic composition analysis if residue object provided
        if residue_obj:
            atom_count = len(list(residue_obj.get_atoms()))
            elements = set(atom.element for atom in residue_obj.get_atoms())
            
            info.update({
                'atom_count': atom_count,
                'elements': sorted(list(elements)),
                'contains_metals': bool(elements & METALS),
                'contains_sulfur': 'S' in elements,
                'contains_nitrogen': 'N' in elements,
                'heavy_atom_count': len([a for a in residue_obj.get_atoms() if a.element != 'H'])
            })
            
            # Simple heuristic for potential redox activity
            info['potentially_redox_active'] = (
                atom_count > 10 and  # Substantial molecule
                len(elements) > 2 and  # Multiple element types
                ('N' in elements or 'S' in elements or bool(elements & METALS))  # Likely redox atoms
            )
        
        if self.use_cache:
            self._cache[residue_name] = info
            
        return info

# ===== CONFIGURATION CLASS =====

@dataclass
class DetectionConfig:
    """Configuration for redox site detection"""
    # Search radius for inventory phase (finding bonds around redox centers)
    bond_search_distance: float = 4.0

    # Disulfide bond detection threshold (SG-SG distance in Å)
    disulfide_distance_threshold: float = 2.5

    # What to include in inventory
    include_non_standard_residues: bool = True
    include_redox_amino_acids: bool = False

    # Atom selection criteria for inventory display
    inventory_include_metals: List[str] = field(default_factory=list)  # Empty list means all metals
    inventory_include_nonmetals: List[str] = field(default_factory=list)  # Empty list means all non-metals

    # Interactive settings
    interactive_mode: bool = True
    verbose: bool = True

class DetectionConfigInterface:
    """Interface for configuring detection parameters"""
    
    def __init__(self, console: Console = None, processor=None):
        self.config = DetectionConfig()
        self.console = console if console else Console()
        self.processor = processor
    
    def configure_interactively(self) -> DetectionConfig:
        """Interactive configuration interface with genuine back navigation"""
        self.console.print("\n[bold underline]Redox Site Detection Configuration[/bold underline]")

        # State tracking for navigation
        answers = {}
        step = 0

        while True:
            try:
                if step == 0:
                    # Question 1: Metals
                    self.console.print("\n[bold]Metal coordination search:[/bold]")
                    self.console.print(
                        "Detect each metal and its coordinating non-metal atoms "
                        "within a cutoff distance."
                    )
                    self._display_recognized_metals()

                    metals_input = prompt_with_navigation(
                        processor=self.processor,
                        prompt="Metal elements to scan (comma-separated list, 'none', or 'all')",
                        default=answers.get('metals', "all"),
                        module="Redox Detector Config",
                        description="Select metals to include in inventory",
                        allow_back=False  # First question
                    ).strip().upper()
                    answers['metals'] = metals_input
                    # Halo metals immediately so the user sees what
                    # their filter caught. Persists through Q2/Q3 and
                    # is replaced by section 2's halo after Q4.
                    _show_inventory_metals_preview(self.processor, answers)
                    step = 1

                elif step == 1:
                    # Question 2: Non-metals
                    nonmetals_input = prompt_with_navigation(
                        processor=self.processor,
                        prompt="Coordinating non-metal elements (comma-separated list, 'none', or 'all')",
                        default=answers.get('nonmetals', "all"),
                        module="Redox Detector Config",
                        description="Select non-metals to include in inventory",
                        allow_back=True
                    ).strip().upper()
                    answers['nonmetals'] = nonmetals_input
                    step = 2

                elif step == 2:
                    # Question 3: Bond search radius
                    bond_distance = prompt_float_with_retry(
                        processor=self.processor,
                        prompt="Search radius around each metal (Å)",
                        default=answers.get('bond_distance', self.config.bond_search_distance),
                        min_value=0.1,
                        max_value=20.0,
                        module="Redox Detector Config",
                        description="Set inventory search radius"
                    )
                    answers['bond_distance'] = bond_distance
                    step = 3

                elif step == 3:
                    # Question 4: Disulfide threshold
                    self.console.print("\n[bold]Disulfide bond cutoff:[/bold]")
                    disulfide_threshold = prompt_float_with_retry(
                        processor=self.processor,
                        prompt="Detect CYS SG–SG pairs within this distance as candidate disulfides (Å)",
                        default=answers.get('disulfide_threshold', self.config.disulfide_distance_threshold),
                        min_value=1.0,
                        max_value=5.0,
                        module="Redox Detector Config",
                        description="Set disulfide bond detection threshold"
                    )
                    answers['disulfide_threshold'] = disulfide_threshold
                    # Section transition AFTER the answer: replace
                    # section-1's coordination halo with section-2's
                    # CYS SG halo. Setting it before the prompt would
                    # blow away section-1's halo before the user gets
                    # a chance to see it.
                    _show_inventory_disulfide_preview(self.processor)
                    step = 4

                elif step == 4:
                    # Question 5: Non-standard residues
                    self.console.print("\n[bold]Additional centers to detect:[/bold]")
                    include_nonstandard = confirm_with_context(
                        processor=self.processor,
                        prompt="Non-standard residues (HEM, FAD, etc.)?",
                        default=answers.get('include_nonstandard', self.config.include_non_standard_residues),
                        module="Redox Detector Config",
                        description="Include non-standard residues in inventory",
                        allow_back=True
                    )
                    answers['include_nonstandard'] = include_nonstandard
                    # Section transition: drop the disulfide halo and
                    # show whichever residue classes are currently
                    # toggled on. Q6 will refresh again to reflect
                    # the redox-AA answer alongside this one.
                    _show_inventory_residue_preview(self.processor, answers)
                    step = 5

                elif step == 5:
                    # Question 6: Redox amino acids
                    include_redox_aa = confirm_with_context(
                        processor=self.processor,
                        prompt="Redox-active amino acids (TYR, TRP, PHE, MET)?",
                        default=answers.get('include_redox_aa', self.config.include_redox_amino_acids),
                        module="Redox Detector Config",
                        description="Include redox-active amino acids in inventory",
                        allow_back=True
                    )
                    answers['include_redox_aa'] = include_redox_aa
                    _show_inventory_residue_preview(self.processor, answers)
                    # All questions answered - break out of loop
                    break

            except NavigationException:
                # User pressed back - go to previous question
                if step > 0:
                    step -= 1
                    self.console.print(f"[grey50]← Going back to previous question...[/grey50]\n")
                # If step is 0, we can't go further back

        # Apply all answers to config
        if answers['metals'] == 'NONE':
            self.config.inventory_include_metals = ['NONE']
        elif answers['metals'] == 'ALL' or answers['metals'] == '':
            self.config.inventory_include_metals = []
        else:
            self.config.inventory_include_metals = [m.strip() for m in answers['metals'].split(',') if m.strip()]

        if answers['nonmetals'] == 'NONE':
            self.config.inventory_include_nonmetals = ['NONE']
        elif answers['nonmetals'] == 'ALL' or answers['nonmetals'] == '':
            self.config.inventory_include_nonmetals = []
        else:
            self.config.inventory_include_nonmetals = [nm.strip() for nm in answers['nonmetals'].split(',') if nm.strip()]

        self.config.bond_search_distance = answers['bond_distance']
        self.config.disulfide_distance_threshold = answers['disulfide_threshold']
        self.config.include_non_standard_residues = answers['include_nonstandard']
        self.config.include_redox_amino_acids = answers['include_redox_aa']
        
        # Display configuration summary
        self.console.print("\n[bold green]Configuration complete![/bold green]")
        
        config_summary_table = Table(title="Final Configuration")
        config_summary_table.add_column("Setting", style="cyan")
        config_summary_table.add_column("Value", style="green")
        
        config_summary_table.add_row("Inventory search radius", f"{self.config.bond_search_distance}Å")
        config_summary_table.add_row("Disulfide detection threshold", f"{self.config.disulfide_distance_threshold}Å")

        # Display inventory atom filtering
        if not self.config.inventory_include_metals:
            metals_value = "all"
        elif self.config.inventory_include_metals == ['NONE']:
            metals_value = "none"
        else:
            metals_value = ','.join(self.config.inventory_include_metals)
        config_summary_table.add_row("Inventory metals", metals_value)
        
        if not self.config.inventory_include_nonmetals:
            nonmetals_value = "all"
        elif self.config.inventory_include_nonmetals == ['NONE']:
            nonmetals_value = "none"
        else:
            nonmetals_value = ','.join(self.config.inventory_include_nonmetals)
        config_summary_table.add_row("Inventory non-metals", nonmetals_value)
        
        config_summary_table.add_row("Include non-standard residues", str(self.config.include_non_standard_residues))
        config_summary_table.add_row("Include redox amino acids", str(self.config.include_redox_amino_acids))
        
        self.console.print(config_summary_table)
        
        return self.config
    
    def _display_recognized_metals(self):
        """Display all recognized metals in an organized table"""
        metals_table = Table(title="Recognized metal elements", show_header=True, header_style="bold magenta")
        metals_table.add_column("Group", style="cyan", width=20)
        metals_table.add_column("Elements", style="green")
        
        # Group 1: Alkali metals
        alkali = ['LI', 'NA', 'K', 'RB', 'CS', 'FR']
        metals_table.add_row("Alkali metals", ", ".join(alkali))
        
        # Group 2: Alkaline earth metals  
        alkaline_earth = ['BE', 'MG', 'CA', 'SR', 'BA', 'RA']
        metals_table.add_row("Alkaline earth", ", ".join(alkaline_earth))
        
        # Transition metals (split into rows for readability)
        transition_1 = ['SC', 'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN']
        transition_2 = ['Y', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD']
        transition_3 = ['LA', 'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG']
        transition_4 = ['AC', 'RF', 'DB', 'SG', 'BH', 'HS', 'MT', 'DS', 'RG', 'CN']
        
        metals_table.add_row("Transition (3d)", ", ".join(transition_1))
        metals_table.add_row("Transition (4d)", ", ".join(transition_2))
        metals_table.add_row("Transition (5d)", ", ".join(transition_3))
        metals_table.add_row("Transition (6d)", ", ".join(transition_4))
        
        # Lanthanides
        lanthanides = ['CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU']
        metals_table.add_row("Lanthanides", ", ".join(lanthanides))
        
        # Actinides
        actinides = ['TH', 'PA', 'U', 'NP', 'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM', 'MD', 'NO', 'LR']
        metals_table.add_row("Actinides", ", ".join(actinides))
        
        # Post-transition metals
        post_transition = ['AL', 'GA', 'IN', 'SN', 'TL', 'PB', 'BI', 'PO']
        metals_table.add_row("Post-transition", ", ".join(post_transition))
        
        # Metalloids
        metalloids = ['SB', 'TE']
        metals_table.add_row("Metalloids", ", ".join(metalloids))
        
        self.console.print(metals_table)
        self.console.print("[grey50]Tip: Common examples include CA, MG, FE, ZN, CU, MN[/grey50]")

# ===== INVENTORY SCANNERS =====

class MetalIonScanner:
    """Scanner for all metal ions in structure"""
    
    def __init__(self, config: DetectionConfig):
        self.config = config
        self.ccd_querier = SimpleCCDQuerier()
    
    def scan_structure(self, structure: Structure, existing_centers: List[RedoxCenter] = None,
                      selected_chains: List[str] = None) -> List[RedoxCenter]:
        """Find ALL metal ions in structure, skipping residues already detected as organic cofactors"""
        metal_centers = []

        # Build set of residues already detected (e.g., as organic cofactors)
        existing_residues = set()
        if existing_centers:
            existing_residues = {(c.chain, c.resid, c.insertion_code) for c in existing_centers}

        logger.debug(f"Scanning for metal ions...")
        if existing_residues:
            logger.debug(f"Skipping {len(existing_residues)} residues already detected as organic cofactors")

        for model in structure:
            for chain in model:
                if selected_chains and chain.id not in selected_chains:
                    continue

                for residue in chain:
                    res_key = (chain.id, residue.id[1], residue.id[2])

                    # Skip residues already detected as organic cofactors
                    if res_key in existing_residues:
                        logger.debug(f"Skipping residue {residue.resname} {chain.id}:{residue.id[1]} - already detected as organic cofactor")
                        continue

                    # Check each atom for metal character (process all metal atoms individually)
                    # Use get_unpacked_list() to get ALL atoms including all alternate locations
                    for atom in residue.get_unpacked_list():
                        element = atom.element.upper()
                        if element in METALS:
                            # Apply inventory metal filtering
                            if not self.config.inventory_include_metals:
                                # Empty list means include all metals
                                pass
                            elif self.config.inventory_include_metals == ['NONE']:
                                # Special marker for no metals
                                continue
                            elif element not in [m.upper() for m in self.config.inventory_include_metals]:
                                # This metal is not in the include list, skip it
                                continue

                            center = RedoxCenter(
                                chain=chain.id,
                                resname=residue.resname,
                                resid=residue.id[1],
                                atom_name=atom.name,
                                insertion_code=residue.id[2],
                                altloc=atom.altloc,
                                coords=tuple(round(x, 3) for x in atom.coord),
                                element=atom.element,
                                center_type=CenterType.METAL_ION,
                                properties={
                                    'is_metal': True,
                                    'metal_type': atom.element.upper(),
                                    'serial_number': atom.serial_number,
                                    'occupancy': atom.occupancy,
                                    'bfactor': atom.bfactor
                                }
                            )
                            metal_centers.append(center)

                            logger.debug(f"Found metal: {atom.element} {atom.name} in {residue.resname} {chain.id}:{residue.id[1]}")

        logger.debug(f"Found {len(metal_centers)} metal centers")
        return metal_centers


# Map: input PDB residue name -> cofactor-family hint. Consumed by the
# transformer selector to short-circuit per-transformer evaluation. Keys
# include both the PDB-CCD canonical names users will load AND the
# ProPrep-emitted state codes (so re-loaded post-ProPrep PDBs still hint
# correctly).
_COFACTOR_FAMILY_MAP = {
    # Flavin: FMN side
    "FMN": "flavin_fmn", "FMS": "flavin_fmn", "FMR": "flavin_fmn",
    "FMH": "flavin_fmn", "FMQ": "flavin_fmn",
    # Flavin: FAD side
    "FAD": "flavin_fad", "FAS": "flavin_fad", "FAR": "flavin_fad",
    "FAH": "flavin_fad", "FAQ": "flavin_fad",
    # Nicotinamide: NADP side
    "NAP": "nicotinamide_nadp", "NDP": "nicotinamide_nadp", "NP2": "nicotinamide_nadp",
    # Nicotinamide: NAD side ("NAI" is PDB-CCD NADH; "NDH" is our library code)
    "NAD": "nicotinamide_nad", "NAI": "nicotinamide_nad", "NDH": "nicotinamide_nad",
    # Pterin: biopterin family
    "H4B": "pterin_biopterin", "H4C": "pterin_biopterin", "H3R": "pterin_biopterin",
    "H2Q": "pterin_biopterin", "H2B": "pterin_biopterin", "BIO": "pterin_biopterin",
}


def _organic_cofactor_family(resname: str) -> Optional[str]:
    """Return the cofactor-family hint (matches transformer TRANSFORMER_NAME)
    for a known organic-cofactor residue code, or None if the residue is not
    in the cofactor library."""
    return _COFACTOR_FAMILY_MAP.get(resname)


class NonStandardResidueScanner:
    """Scanner for non-standard residues (organic cofactors)"""
    
    def __init__(self, config: DetectionConfig):
        self.config = config
        self.ccd_querier = SimpleCCDQuerier()
    
    def scan_structure(self, structure: Structure, existing_centers: List[RedoxCenter],
                      selected_chains: List[str] = None) -> List[RedoxCenter]:
        """Find non-standard residues as potential organic cofactors"""
        if not self.config.include_non_standard_residues:
            return []
        
        organic_centers = []
        existing_residues = {(c.chain, c.resid, c.insertion_code) for c in existing_centers}
        
        logger.debug("Scanning for non-standard residues...")
        
        for model in structure:
            for chain in model:
                if selected_chains and chain.id not in selected_chains:
                    continue
                    
                for residue in chain:
                    res_key = (chain.id, residue.id[1], residue.id[2])

                    # Skip if already included as metal center
                    if res_key in existing_residues:
                        continue

                    res_info = self.ccd_querier.get_basic_residue_info(residue.resname, residue)

                    # Only include non-standard residues that are NOT single-atom metal ions
                    # Single-atom metals (CU, ZN, FE, etc.) should be handled by MetalIonScanner
                    # Multi-atom metal-containing residues (hemes, Fe-S clusters) ARE organic cofactors
                    if res_info['is_non_standard']:
                        # Skip single-atom metal ions
                        if res_info.get('atom_count', 0) == 1 and res_info.get('contains_metals', False):
                            logger.debug(f"Skipping single-atom metal ion {residue.resname} - will be handled by MetalIonScanner")
                            continue

                        # Calculate centroid coordinates for whole residue
                        coords = np.array([atom.coord for atom in residue.get_atoms()])
                        centroid = tuple(round(x, 3) for x in np.mean(coords, axis=0))

                        # Distinguish organometallic (with metal) from organic (without)
                        contains_metal = res_info.get('contains_metals', False)
                        if contains_metal:
                            center_type = CenterType.ORGANOMETALLIC_COFACTOR
                        else:
                            center_type = CenterType.ORGANIC_COFACTOR

                        # Hint which cofactor-family transformer covers this
                        # residue, so the transformer selector can short-circuit
                        # the per-transformer evaluate_redox_site loop. Pure
                        # property — does not affect site grouping.
                        cofactor_family = _organic_cofactor_family(residue.resname)

                        center = RedoxCenter(
                            chain=chain.id,
                            resname=residue.resname,
                            resid=residue.id[1],
                            atom_name=None,  # Whole residue
                            insertion_code=residue.id[2],
                            coords=centroid,
                            element=None,  # Not a single element
                            center_type=center_type,
                            properties={
                                'is_organic_cofactor': not contains_metal,
                                'is_organometallic': contains_metal,
                                'atom_count': res_info.get('atom_count', 0),
                                'elements': res_info.get('elements', []),
                                'potentially_redox_active': res_info.get('potentially_redox_active', False),
                                'cofactor_family': cofactor_family,
                            }
                        )
                        organic_centers.append(center)
                        
                        logger.debug(f"Found non-standard residue: {residue.resname} {chain.id}:{residue.id[1]}")
        
        logger.debug(f"Found {len(organic_centers)} non-standard residues")
        return organic_centers

class RedoxAminoAcidScanner:
    """Scanner for redox-active amino acids"""
    
    def __init__(self, config: DetectionConfig):
        self.config = config
        self.ccd_querier = SimpleCCDQuerier()
        self.redox_amino_acids = {'TYR', 'TRP', 'PHE', 'MET', 'CYS'}
    
    def scan_structure(self, structure: Structure, existing_centers: List[RedoxCenter],
                      selected_chains: List[str] = None) -> List[RedoxCenter]:
        """Find redox-active amino acids"""
        if not self.config.include_redox_amino_acids:
            return []
        
        amino_centers = []
        existing_residues = {(c.chain, c.resid, c.insertion_code) for c in existing_centers}
        
        logger.debug("Scanning for redox-active amino acids...")
        
        for model in structure:
            for chain in model:
                if selected_chains and chain.id not in selected_chains:
                    continue
                    
                for residue in chain:
                    res_key = (chain.id, residue.id[1], residue.id[2])
                    
                    # Skip if already included
                    if res_key in existing_residues:
                        continue
                    
                    if residue.resname in self.redox_amino_acids:
                        # For amino acids, use specific atoms as centers
                        # Determine target atom name based on residue type
                        target_atom_name = None

                        if residue.resname == 'CYS':
                            target_atom_name = 'SG'
                        elif residue.resname == 'TYR':
                            target_atom_name = 'CZ'
                        elif residue.resname in ['TRP', 'PHE']:
                            target_atom_name = 'CD1'
                        elif residue.resname == 'MET':
                            target_atom_name = 'SD'

                        # Find all altlocs of the target atom (handles alternate locations)
                        if target_atom_name:
                            target_atoms = [atom for atom in residue.get_unpacked_list()
                                          if atom.name == target_atom_name]

                            for atom in target_atoms:
                                properties = {
                                    'is_redox_amino_acid': True,
                                    'amino_acid_type': residue.resname
                                }

                                # Special properties for cysteine
                                if residue.resname == 'CYS':
                                    properties['can_form_disulfide'] = True
                                    properties['can_coordinate_metal'] = True

                                center = RedoxCenter(
                                    chain=chain.id,
                                    resname=residue.resname,
                                    resid=residue.id[1],
                                    atom_name=atom.name,
                                    insertion_code=residue.id[2],
                                    altloc=atom.altloc,
                                    coords=tuple(round(x, 3) for x in atom.coord),
                                    element=atom.element,
                                    center_type=CenterType.REDOX_AMINO_ACID,
                                    properties=properties
                                )
                                amino_centers.append(center)

                                logger.debug(f"Found redox amino acid: {residue.resname} {chain.id}:{residue.id[1]} {atom.name} altloc={atom.altloc}")
        
        logger.debug(f"Found {len(amino_centers)} redox-active amino acids")
        return amino_centers


class DisulfideBondScanner:
    """Scanner for disulfide bonds between CYS residues"""

    def __init__(self, config: DetectionConfig):
        self.config = config
        self.distance_threshold = config.disulfide_distance_threshold  # Å for SG-SG bond
        self.detected_bonds = []  # Store (chain1, res1, chain2, res2, distance)

    def scan_structure(self, structure: Structure, existing_centers: List[RedoxCenter],
                      selected_chains: List[str] = None,
                      source_pdb_file: str = None) -> List[RedoxCenter]:
        """
        Detect disulfide bonds and return as paired RedoxCenters

        Strategy:
        1. Try SSBOND records from PDB file (if available)
        2. Verify SSBOND records against structure coordinates
        3. If no SSBOND or verification fails, do distance-based detection
        4. Return centers as pairs (each bond = 2 centers flagged as disulfide)
        """
        # Build set of CYS residues already detected as individual redox amino acids
        existing_cys_residues = set()
        for center in existing_centers:
            if center.resname == 'CYS' and center.center_type == CenterType.REDOX_AMINO_ACID:
                existing_cys_residues.add((center.chain, center.resid, center.insertion_code))

        logger.debug("Scanning for disulfide bonds...")

        # Step 1: Try SSBOND records first
        ssbond_records = []
        if source_pdb_file:
            ssbond_records = self._parse_ssbond_records(source_pdb_file, selected_chains)
            if ssbond_records:
                logger.debug(f"Found {len(ssbond_records)} SSBOND record(s) in PDB file")

        # Step 2: Verify or detect bonds
        if ssbond_records:
            # Verify SSBOND records against structure
            verified_bonds = self._verify_ssbond_records(structure, ssbond_records)
            self.detected_bonds = verified_bonds
            logger.debug(f"Verified {len(verified_bonds)} of {len(ssbond_records)} SSBOND bonds")
        else:
            # No SSBOND records, use distance-based detection
            detected_bonds = self._detect_by_distance(structure, selected_chains)
            self.detected_bonds = detected_bonds
            logger.debug(f"Distance-based detection found {len(detected_bonds)} disulfide bond(s)")

        # Step 3: Create RedoxCenter objects for disulfide-bonded CYS
        disulfide_centers = []
        for bond_data in self.detected_bonds:
            # Handle both old format (5 elements) and new format (6 elements with altloc)
            if len(bond_data) == 6:
                chain1, res1, chain2, res2, distance, altloc = bond_data
            else:
                chain1, res1, chain2, res2, distance = bond_data
                altloc = ''

            # Create center for first CYS (if not already in existing_centers)
            if (chain1, res1, '') not in existing_cys_residues:
                center1 = self._create_cys_center(structure, chain1, res1, distance,
                                                  partner_chain=chain2, partner_res=res2, altloc=altloc)
                if center1:
                    disulfide_centers.append(center1)

            # Create center for second CYS (if not already in existing_centers)
            if (chain2, res2, '') not in existing_cys_residues:
                center2 = self._create_cys_center(structure, chain2, res2, distance,
                                                  partner_chain=chain1, partner_res=res1, altloc=altloc)
                if center2:
                    disulfide_centers.append(center2)

        logger.debug(f"Created {len(disulfide_centers)} disulfide bond center(s)")
        return disulfide_centers

    def _parse_ssbond_records(self, pdb_file: str, selected_chains: List[str] = None) -> List[Tuple]:
        """Parse SSBOND records from PDB file header"""
        ssbond_records = []

        try:
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith('SSBOND'):
                        # SSBOND format (PDB specification):
                        # Columns 16    : Chain 1
                        # Columns 18-21 : Residue 1 sequence number
                        # Columns 30    : Chain 2
                        # Columns 32-35 : Residue 2 sequence number
                        chain1 = line[15:16].strip()
                        res1 = int(line[17:21].strip())
                        chain2 = line[29:30].strip()
                        res2 = int(line[31:35].strip())

                        # Apply chain filtering
                        if selected_chains:
                            if chain1 not in selected_chains or chain2 not in selected_chains:
                                continue

                        ssbond_records.append((chain1, res1, chain2, res2))
                        logger.debug(f"SSBOND: {chain1}:{res1} ←→ {chain2}:{res2}")

                    # Stop at first ATOM/MODEL record
                    if line.startswith(('ATOM', 'MODEL')):
                        break
        except Exception as e:
            logger.warning(f"Could not parse SSBOND records: {e}")

        return ssbond_records

    def _verify_ssbond_records(self, structure: Structure,
                               ssbond_records: List[Tuple]) -> List[Tuple]:
        """Verify SSBOND records by checking CYS exists and measuring SG-SG distance"""
        verified_bonds = []

        for chain1, res1, chain2, res2 in ssbond_records:
            try:
                # Find both CYS residues
                cys1 = None
                cys2 = None

                for model in structure:
                    if chain1 in model and chain2 in model:
                        for residue in model[chain1]:
                            if residue.id[1] == res1 and residue.resname in ['CYS', 'CYX']:
                                cys1 = residue
                                break
                        for residue in model[chain2]:
                            if residue.id[1] == res2 and residue.resname in ['CYS', 'CYX']:
                                cys2 = residue
                                break
                        break

                if cys1 and cys2 and 'SG' in cys1 and 'SG' in cys2:
                    # Measure actual SG-SG distance
                    distance = cys1['SG'] - cys2['SG']

                    # Verify distance is reasonable for disulfide bond
                    if distance <= 3.0:  # Allow slightly more than threshold for SSBOND verification
                        verified_bonds.append((chain1, res1, chain2, res2, distance))
                        logger.debug(f"✓ Verified SSBOND {chain1}:{res1} ←→ {chain2}:{res2} ({distance:.2f}Å)")
                    else:
                        logger.warning(f"✗ SSBOND {chain1}:{res1} ←→ {chain2}:{res2} distance too large ({distance:.2f}Å)")
                else:
                    logger.warning(f"✗ Could not find CYS residues for SSBOND {chain1}:{res1} ←→ {chain2}:{res2}")

            except Exception as e:
                logger.warning(f"Error verifying SSBOND {chain1}:{res1} ←→ {chain2}:{res2}: {e}")

        return verified_bonds

    def _detect_by_distance(self, structure: Structure,
                           selected_chains: List[str] = None) -> List[Tuple]:
        """Detect disulfide bonds by measuring SG-SG distances

        Handles alternate locations by only pairing atoms with matching altlocs.
        """
        bonds = []
        cys_list = []  # (chain_id, res_id, SG_atom, residue)

        # Collect all CYS/CYX residues with SG atoms (including all altlocs)
        for model in structure:
            for chain in model:
                if selected_chains and chain.id not in selected_chains:
                    continue

                for residue in chain:
                    if residue.resname in ['CYS', 'CYX']:
                        # Get all SG atoms including alternate locations
                        sg_atoms = [atom for atom in residue.get_unpacked_list() if atom.name == 'SG']
                        for sg_atom in sg_atoms:
                            cys_list.append((chain.id, residue.id[1], sg_atom, residue))

        # Pairwise distance check - ONLY pair atoms with matching altlocs
        for i, (chain1, res1, sg1, residue1) in enumerate(cys_list):
            for chain2, res2, sg2, residue2 in cys_list[i+1:]:
                # Check if altlocs are compatible (both same, or one/both blank)
                altloc1 = sg1.altloc.strip()
                altloc2 = sg2.altloc.strip()

                # Allow pairing if:
                # 1. Both have same altloc (A-A, B-B)
                # 2. One or both are blank (default conformation)
                if altloc1 == altloc2 or altloc1 == '' or altloc2 == '':
                    distance = sg1 - sg2  # BioPython calculates distance

                    if distance <= self.distance_threshold:
                        # Use the non-blank altloc if one is blank
                        bond_altloc = altloc1 if altloc1 else altloc2
                        bonds.append((chain1, res1, chain2, res2, distance, bond_altloc))
                        logger.debug(f"Disulfide bond detected: {chain1}:{res1}{bond_altloc} ←→ {chain2}:{res2}{bond_altloc} ({distance:.2f}Å)")

        return bonds

    def _create_cys_center(self, structure: Structure, chain_id: str, res_id: int,
                          distance: float, partner_chain: str, partner_res: int, altloc: str = '') -> RedoxCenter:
        """Create a RedoxCenter for a CYS involved in disulfide bond

        Args:
            altloc: Specific alternate location to use (empty string for default)
        """
        try:
            for model in structure:
                if chain_id in model:
                    for residue in model[chain_id]:
                        if residue.id[1] == res_id and residue.resname in ['CYS', 'CYX']:
                            # Get all SG atoms and find the one matching the altloc
                            sg_atoms = [atom for atom in residue.get_unpacked_list() if atom.name == 'SG']

                            # Find the SG atom with matching altloc
                            sg_atom = None
                            for atom in sg_atoms:
                                if atom.altloc.strip() == altloc.strip():
                                    sg_atom = atom
                                    break

                            # If no exact match, try blank altloc (default conformation)
                            if not sg_atom:
                                for atom in sg_atoms:
                                    if atom.altloc.strip() == '':
                                        sg_atom = atom
                                        break

                            # If still no match, use first one
                            if not sg_atom and sg_atoms:
                                sg_atom = sg_atoms[0]

                            if sg_atom:
                                center = RedoxCenter(
                                    chain=chain_id,
                                    resname=residue.resname,
                                    resid=res_id,
                                    atom_name='SG',
                                    insertion_code=residue.id[2],
                                    altloc=sg_atom.altloc,
                                    coords=tuple(round(x, 3) for x in sg_atom.coord),
                                    element='S',
                                    center_type=CenterType.REDOX_AMINO_ACID,
                                    properties={
                                        'is_disulfide_bonded': True,
                                        'disulfide_partner_chain': partner_chain,
                                        'disulfide_partner_res': partner_res,
                                        'disulfide_bond_distance': round(distance, 2),
                                        'can_form_disulfide': True,
                                        'amino_acid_type': 'CYS'
                                    }
                                )
                                return center
        except Exception as e:
            logger.warning(f"Could not create center for {chain_id}:{res_id}: {e}")

        return None

    def get_all_cys_residues(self, structure: Structure,
                            selected_chains: List[str] = None) -> List[Tuple[str, int, str]]:
        """Get all CYS/CYX residues for manual specification interface"""
        cys_list = []

        for model in structure:
            for chain in model:
                if selected_chains and chain.id not in selected_chains:
                    continue

                for residue in chain:
                    if residue.resname in ['CYS', 'CYX'] and 'SG' in residue:
                        cys_list.append((chain.id, residue.id[1], residue.resname))

        return sorted(cys_list)  # Sort by chain then residue number


# ===== RELATIONAL ANALYZER =====

class RelationalAnalyzer:
    """Analyze relationships between centers (bonds, clusters, etc.)"""
    
    def __init__(self, config: DetectionConfig):
        self.config = config
    
    def find_nearby_bonds(self, structure: Structure, centers: List[RedoxCenter],
                         selected_chains: List[str] = None) -> Dict[Tuple[str, int, str], List[Dict[str, Any]]]:
        """Find nearby atoms that could form bonds with each center"""
        nearby_bonds = {}
        
        # Build NeighborSearch for efficient distance queries
        atoms = list(structure.get_atoms())
        ns = NeighborSearch(atoms)
        
        for center in centers:
            center_key = (center.chain, center.resid, center.insertion_code)
            bonds = []
            
            # Find atoms within bond search distance
            nearby_atoms = ns.search(center.coords, self.config.bond_search_distance)
            
            for atom in nearby_atoms:
                # Skip the center atom itself
                if (atom.parent.parent.id == center.chain and
                    atom.parent.id[1] == center.resid and
                    atom.parent.id[2] == center.insertion_code and
                    atom.name == center.atom_name):
                    continue

                # Skip other atoms from the same residue as the center
                if (atom.parent.parent.id == center.chain and
                    atom.parent.id[1] == center.resid and
                    atom.parent.id[2] == center.insertion_code):
                    continue
                
                # Apply inventory atom filtering
                if not self._atom_meets_inventory_criteria(atom):
                    continue
                
                distance = np.linalg.norm(np.array(center.coords) - np.array(atom.coord))
                
                # Classify the potential bond
                bond_type, chemical_type = classify_bond_types(
                    center.element or 'X', atom.element,
                    distance,
                    center.resname, atom.parent.resname,
                    center.resid, atom.parent.id[1],
                    center.chain, atom.parent.parent.id,
                    center.atom_name, atom.name
                )
                
                bond_info = {
                    'chain': atom.parent.parent.id,
                    'resname': atom.parent.resname,
                    'resid': atom.parent.id[1],
                    'atom_name': atom.name,
                    'element': atom.element,
                    'distance': distance,
                    'bond_type': bond_type,
                    'chemical_type': chemical_type,
                    'coords': tuple(round(x, 3) for x in atom.coord)
                }
                bonds.append(bond_info)
            
            # Sort by distance
            bonds.sort(key=lambda x: x['distance'])
            nearby_bonds[center_key] = bonds
        
        return nearby_bonds
    
    def find_metal_clusters(self, centers: List[RedoxCenter]) -> Dict[int, List[RedoxCenter]]:
        """No automatic metal clustering - user decides grouping"""
        # Return empty dict - user will group centers manually in Phase 3
        return {}
    
    def _atom_meets_inventory_criteria(self, atom) -> bool:
        """Check if atom meets the inventory filtering criteria"""
        element = atom.element.upper()
        atom_name = atom.name.upper()
        
        # Check if atom is metal
        if element in METALS:
            # Metal atom - check against inventory_include_metals criteria
            if not self.config.inventory_include_metals:
                # Empty list means include all metals
                return True
            elif self.config.inventory_include_metals == ['NONE']:
                # Special marker for no metals
                return False
            else:
                # Check if this metal is in the include list
                return element in [m.upper() for m in self.config.inventory_include_metals]
        else:
            # Non-metal atom - check against inventory_include_nonmetals criteria
            if not self.config.inventory_include_nonmetals:
                # Empty list means include all non-metals
                return True
            elif self.config.inventory_include_nonmetals == ['NONE']:
                # Special marker for no non-metals
                return False
            else:
                # Check if element name OR atom name is in the include list
                nonmetal_criteria = [nm.upper() for nm in self.config.inventory_include_nonmetals]
                return element in nonmetal_criteria or atom_name in nonmetal_criteria


# ===== INVENTORY DISPLAY =====

class InventoryDisplay:
    """Display comprehensive inventory with relational information"""
    
    def __init__(self, config: DetectionConfig, console: Console = None):
        self.config = config
        self.console = console if console else Console()
    
    def display_inventory(self, centers: List[RedoxCenter], 
                         nearby_bonds: Dict[Tuple[str, int, str], List[Dict[str, Any]]],
                         metal_clusters: Dict[int, List[RedoxCenter]]) -> None:
        """Display comprehensive inventory with all relational information"""
        
        if not centers:
            self.console.print("[yellow]No redox centers found.[/yellow]")
            return
        
        # Create main inventory table
        inventory_table = Table(title=f"Redox Center Inventory ({len(centers)} centers found)")
        inventory_table.add_column("Center ID", style="cyan", width=10)
        inventory_table.add_column("Center Desc.", style="yellow", width=20)
        inventory_table.add_column("Nearby Atom", style="green", width=20)
        inventory_table.add_column("Distance", style="blue", width=10)
        
        # Group centers by type for organized display
        metal_centers = [c for c in centers if c.center_type == CenterType.METAL_ION]
        organometallic_centers = [c for c in centers if c.center_type == CenterType.ORGANOMETALLIC_COFACTOR]
        organic_centers = [c for c in centers if c.center_type == CenterType.ORGANIC_COFACTOR]
        disulfide_centers = [c for c in centers if c.properties.get('is_disulfide_bonded', False)]
        other_amino_centers = [c for c in centers if c.center_type == CenterType.REDOX_AMINO_ACID and not c.properties.get('is_disulfide_bonded', False)]
        # Combine organic and organometallic for display purposes (both are "cofactor" category)
        cofactor_centers = organometallic_centers + organic_centers

        center_id = 1

        # Display disulfide bonds first (highlighted)
        if disulfide_centers:
            # Group into pairs for better display
            displayed_pairs = set()
            for center in disulfide_centers:
                center_key = (center.chain, center.resid)
                partner_key = (center.properties['disulfide_partner_chain'],
                              center.properties['disulfide_partner_res'])
                pair_key = tuple(sorted([center_key, partner_key]))

                if pair_key not in displayed_pairs:
                    displayed_pairs.add(pair_key)
                    # Find partner center
                    partner = None
                    for c in disulfide_centers:
                        if (c.chain == partner_key[0] and c.resid == partner_key[1]) or \
                           (c.chain == partner_key[1] and c.resid == partner_key[1]):
                            partner = c
                            break

                    # Display both CYS in the pair with special formatting
                    distance = center.properties['disulfide_bond_distance']
                    self._add_disulfide_center_to_table(inventory_table, center, center_id, nearby_bonds, distance)
                    center_id += 1

                    if partner:
                        self._add_disulfide_center_to_table(inventory_table, partner, center_id, nearby_bonds, distance)
                        center_id += 1

        # Display metal centers
        if metal_centers:
            for center in metal_centers:
                self._add_center_to_table(inventory_table, center, center_id, nearby_bonds)
                center_id += 1

        # Display cofactor centers (organometallic + organic)
        if cofactor_centers:
            for center in cofactor_centers:
                self._add_center_to_table(inventory_table, center, center_id, nearby_bonds)
                center_id += 1

        # Display other redox amino acid centers (non-disulfide)
        if other_amino_centers:
            for center in other_amino_centers:
                self._add_center_to_table(inventory_table, center, center_id, nearby_bonds)
                center_id += 1

        # Display summary by type first
        self._display_inventory_summary(metal_centers, organometallic_centers, organic_centers, disulfide_centers, other_amino_centers)

        # Then display detailed inventory table
        self.console.print(inventory_table)
    
    def _add_disulfide_center_to_table(self, table: Table, center: 'RedoxCenter', center_id: int, nearby_bonds: Dict, distance: float):
        """Add a disulfide-bonded CYS center with special highlighting"""
        # Create center description with disulfide indicator
        center_desc = f"[bold yellow]{center.element} {center.resname} {center.chain}:{center.resid}"
        if center.altloc:
            center_desc += center.altloc
        center_desc += "[/bold yellow]"

        # Show disulfide bond partner info
        partner_chain = center.properties['disulfide_partner_chain']
        partner_res = center.properties['disulfide_partner_res']
        bond_info = f"[bold yellow]⟷[/bold yellow] {partner_chain}:{partner_res}"
        distance_str = f"[bold yellow]{distance:.2f}Å[/bold yellow]"

        table.add_row(f"[{center_id}]", center_desc, bond_info, distance_str)
        table.add_section()

    def _add_center_to_table(self, table: Table, center: 'RedoxCenter', center_id: int, nearby_bonds: Dict):
        """Add a redox center and its nearby atoms to the inventory table"""
        # Create center description
        # For residue-level centers (element is None), show just the residue
        if center.element is None or center.atom_name is None:
            center_desc = f"{center.resname} {center.chain}:{center.resid}"
        else:
            center_desc = f"{center.element} {center.resname} {center.chain}:{center.resid}"

        if center.altloc:
            center_desc += center.altloc  # e.g., "CU A:154A"

        # Get center key for nearby bonds lookup - need to match the key format used in nearby_bonds
        center_key = (center.chain, center.resid, center.insertion_code)
        nearby_atoms = nearby_bonds.get(center_key, [])

        if not nearby_atoms:
            # Center with no nearby atoms
            table.add_row(f"[{center_id}]", center_desc, "No nearby atoms found", "—")
        else:
            # First row shows center ID and description
            first_atom = nearby_atoms[0]
            atom_desc = f"{first_atom.get('resname', '?')} {first_atom.get('chain', '?')}:{first_atom.get('resid', '?')} {first_atom.get('atom_name', '?')}"
            distance = f"{first_atom.get('distance', 0):.1f}Å"
            table.add_row(f"[{center_id}]", center_desc, atom_desc, distance)

            # Subsequent rows show only nearby atoms (empty center columns)
            for atom in nearby_atoms[1:]:
                atom_desc = f"{atom.get('resname', '?')} {atom.get('chain', '?')}:{atom.get('resid', '?')} {atom.get('atom_name', '?')}"
                distance = f"{atom.get('distance', 0):.1f}Å"
                table.add_row("", "", atom_desc, distance)

        # Add separator line after each center (except the last one will be handled by Rich)
        table.add_section()

    def _display_inventory_summary(self, metal_centers: List, organometallic_centers: List,
                                   organic_centers: List, disulfide_centers: List,
                                   other_amino_centers: List):
        """Display a summary of the inventory by center type"""
        summary_table = Table(title="Inventory Summary", show_header=True)
        summary_table.add_column("Type", style="bold cyan")
        summary_table.add_column("Count", style="bold green")
        summary_table.add_column("Examples", style="grey50")

        # Display disulfide bonds first (with warning if present)
        if disulfide_centers:
            disulfide_bond_count = len(disulfide_centers) // 2  # Pairs
            example_str = f"{disulfide_bond_count} disulfide bond(s)"
            summary_table.add_row("[bold yellow]Disulfide Bonds[/bold yellow]",
                                 f"[bold yellow]{disulfide_bond_count}[/bold yellow]",
                                 f"[bold yellow]{example_str}[/bold yellow]")

        if metal_centers:
            examples = [f"{c.element}" for c in metal_centers[:3]]
            example_str = ", ".join(examples)
            if len(metal_centers) > 3:
                example_str += f" (+ {len(metal_centers) - 3} more)"
            summary_table.add_row("Metal Ions", str(len(metal_centers)), example_str)

        if organometallic_centers:
            examples = [f"{c.resname}" for c in organometallic_centers[:3]]
            example_str = ", ".join(examples)
            if len(organometallic_centers) > 3:
                example_str += f" (+ {len(organometallic_centers) - 3} more)"
            summary_table.add_row("Organometallic Cofactors", str(len(organometallic_centers)), example_str)

        if organic_centers:
            examples = [f"{c.resname}" for c in organic_centers[:3]]
            example_str = ", ".join(examples)
            if len(organic_centers) > 3:
                example_str += f" (+ {len(organic_centers) - 3} more)"
            summary_table.add_row("Organic Cofactors", str(len(organic_centers)), example_str)

        if other_amino_centers:
            examples = [f"{c.resname}" for c in other_amino_centers[:3]]
            example_str = ", ".join(examples)
            if len(other_amino_centers) > 3:
                example_str += f" (+ {len(other_amino_centers) - 3} more)"
            summary_table.add_row("Other Redox Amino Acids", str(len(other_amino_centers)), example_str)

        self.console.print("\n")
        self.console.print(summary_table)

        # Add warning panel if disulfides found
        if disulfide_centers:
            disulfide_count = len(disulfide_centers) // 2
            warning_text = f"[bold yellow]⚠ Found {disulfide_count} disulfide bond(s)[/bold yellow]\n"
            warning_text += "[grey50]Disulfide bonds are often functionally important for protein stability and redox activity.[/grey50]"
            from rich.panel import Panel
            self.console.print(Panel(warning_text, border_style="yellow", padding=(1, 2), expand=False))

# ===== SEARCH PARAMETERS FOR PHASE 4 =====

@dataclass
class SearchParameters:
    """Parameters for all 9 search method combinations in Phase 4"""
    # Search constraint type
    constraint: SearchConstraint
    
    # Distance-based parameters (if constraint == DISTANCE_CUTOFF)
    distance_method: Optional[DistanceMethod] = None
    boundary_definition: Optional[BoundaryDefinition] = None
    
    # Fixed distance parameters
    radius: Optional[float] = None
    
    # Adaptive distance parameters
    min_radius: Optional[float] = None
    max_radius: Optional[float] = None
    radius_step: Optional[float] = None
    target_coordination: Optional[int] = None
    
    # Count-based parameters (if constraint == COUNT_CUTOFF)
    target_residue_types: Optional[List[str]] = None
    target_residue_count: Optional[int] = None  # Total count (for backward compatibility)
    target_residue_count_map: Optional[Dict[str, int]] = None  # Per-type counts
    
    # Custom boundary coordinates (if boundary_definition == MIN_DISTANCE_CUSTOM)
    custom_boundary_coords: Optional[List[Tuple[float, float, float]]] = None
    
    # Atom filtering for this specific search
    search_include_metals: List[str] = field(default_factory=list)  # Empty list means all metals
    search_include_nonmetals: List[str] = field(default_factory=list)  # Empty list means all non-metals
    
    def get_method_description(self) -> str:
        """Get human-readable description of search method"""
        if self.constraint == SearchConstraint.DISTANCE_CUTOFF:
            distance_desc = "Fixed" if self.distance_method == DistanceMethod.FIXED else "Adaptive"
            boundary_desc = {
                BoundaryDefinition.GEOMETRIC_CENTER: "geometric center",
                BoundaryDefinition.MIN_DISTANCE_ALL: "all atoms",
                BoundaryDefinition.MIN_DISTANCE_CUSTOM: "custom atoms"
            }.get(self.boundary_definition, "unknown")
            return f"{distance_desc} distance from {boundary_desc}"
        else:
            if self.target_residue_count_map:
                type_counts = [f"{count} {res_type}" for res_type, count in self.target_residue_count_map.items()]
                return f"Find {', '.join(type_counts)} residues"
            else:
                return f"Find {self.target_residue_count} closest {self.target_residue_types} residues"

@dataclass
class SiteTemplate:
    """Template for automated site refinement of similar redox sites"""
    site_type: str

    # Search configuration
    residue_types: List[str]
    residue_counts: List[int]
    atom_filtering_metals: List[str]
    atom_filtering_nonmetals: List[str]
    boundary_definition: BoundaryDefinition

    # User selections
    residue_selection_choice: str  # e.g., "1-4"
    bond_pairs: List[str]  # e.g., ["1-2", "1-3", "1-4", "1-5"]
    bond_selections: Dict[str, List[Dict[str, str]]]  # "1-2" -> [{"source_atom": "NE2", "target_atom": "SG"}]
    continue_searching: bool

    # Optional fields (must come after required fields)
    custom_boundary_atom_indices: Optional[List[int]] = None  # Atom indices for custom boundary (relative to site)

    # Distance-based search parameters (if search used distance cutoff)
    search_constraint: Optional[str] = None  # "distance" or "count"
    distance_method: Optional[str] = None  # "fixed" or "adaptive"
    distance_radius: Optional[float] = None  # Fixed radius
    distance_min_radius: Optional[float] = None  # Adaptive min
    distance_max_radius: Optional[float] = None  # Adaptive max
    distance_radius_step: Optional[float] = None  # Adaptive step
    distance_target_coordination: Optional[int] = None  # Adaptive target

@dataclass
class SearchResult:
    """Result from a site search operation"""
    detected_atoms: List[Dict[str, Any]]
    detected_residues: List[Dict[str, Any]]
    search_parameters: SearchParameters
    boundary_coords: List[Tuple[float, float, float]]
    search_radius_used: Optional[float] = None
    total_atoms_found: int = 0
    total_residues_found: int = 0
    search_timestamp: str = ""

    def _add_center_to_table(self, table: Table, center: 'RedoxCenter', center_id: int, nearby_bonds: Dict):
        """Add a redox center and its nearby atoms to the inventory table"""
        # Create center description
        # For residue-level centers (element is None), show just the residue
        if center.element is None or center.atom_name is None:
            center_desc = f"{center.resname} {center.chain}:{center.resid}"
        else:
            center_desc = f"{center.element} {center.resname} {center.chain}:{center.resid}"

        if center.altloc:
            center_desc += center.altloc  # e.g., "CU A:154A"

        # Get center key for nearby bonds lookup
        center_key = (center.chain, center.resid, center.atom_name)
        nearby_atoms = nearby_bonds.get(center_key, [])
        
        if not nearby_atoms:
            # Center with no nearby atoms
            table.add_row(f"[{center_id}]", center_desc, "No nearby atoms found", "—")
        else:
            # First row shows center ID and description
            first_atom = nearby_atoms[0]
            atom_desc = f"{first_atom.get('resname', '?')} {first_atom.get('chain', '?')}:{first_atom.get('resid', '?')} {first_atom.get('atom_name', '?')}"
            distance = f"{first_atom.get('distance', 0):.1f}Å"
            table.add_row(f"[{center_id}]", center_desc, atom_desc, distance)
            
            # Subsequent rows show only nearby atoms (empty center columns)
            for atom in nearby_atoms[1:]:
                atom_desc = f"{atom.get('resname', '?')} {atom.get('chain', '?')}:{atom.get('resid', '?')} {atom.get('atom_name', '?')}"
                distance = f"{atom.get('distance', 0):.1f}Å"
                table.add_row("", "", atom_desc, distance)
            
        # Add separator line after each center (except the last one will be handled by Rich)
        table.add_section()
    
    def _display_inventory_summary(self, metal_centers: List, cofactor_centers: List, amino_centers: List):
        """Display a summary of the inventory by center type.

        Note: cofactor_centers includes both organometallic and organic cofactors combined.
        """
        summary_table = Table(title="Inventory Summary", show_header=True)
        summary_table.add_column("Type", style="bold cyan")
        summary_table.add_column("Count", style="bold green")
        summary_table.add_column("Examples", style="grey50")

        if metal_centers:
            examples = [f"{c.element}" for c in metal_centers[:3]]
            example_str = ", ".join(examples)
            if len(metal_centers) > 3:
                example_str += f" (+ {len(metal_centers) - 3} more)"
            summary_table.add_row("Metal Ions", str(len(metal_centers)), example_str)

        if cofactor_centers:
            # Separate organometallic from organic for display
            organometallic = [c for c in cofactor_centers if c.center_type == CenterType.ORGANOMETALLIC_COFACTOR]
            organic = [c for c in cofactor_centers if c.center_type == CenterType.ORGANIC_COFACTOR]
            if organometallic:
                examples = [f"{c.resname}" for c in organometallic[:3]]
                example_str = ", ".join(examples)
                if len(organometallic) > 3:
                    example_str += f" (+ {len(organometallic) - 3} more)"
                summary_table.add_row("Organometallic Cofactors", str(len(organometallic)), example_str)
            if organic:
                examples = [f"{c.resname}" for c in organic[:3]]
                example_str = ", ".join(examples)
                if len(organic) > 3:
                    example_str += f" (+ {len(organic) - 3} more)"
                summary_table.add_row("Organic Cofactors", str(len(organic)), example_str)

        if amino_centers:
            examples = [f"{c.resname}" for c in amino_centers[:3]]
            example_str = ", ".join(examples)
            if len(amino_centers) > 3:
                example_str += f" (+ {len(amino_centers) - 3} more)"
            summary_table.add_row("Redox Amino Acids", str(len(amino_centers)), example_str)

        self.console.print("\n")
        self.console.print(summary_table)

# ===== CENTER SELECTION AND GROUPING =====

class CenterSelectionInterface:
    """Interface for selecting which centers to include in redox sites"""

    def __init__(self, config: DetectionConfig, console: Console = None, processor=None):
        self.config = config
        self.console = console if console else Console()
        self.processor = processor
    
    def select_centers(self, centers: List[RedoxCenter]) -> tuple:
        """
        Two-stage center selection with disulfide bond review

        Returns:
            Tuple of (selected_centers, confirmed_disulfide_pairs)
        """
        if not centers:
            self.console.print("[yellow]No centers found to select from.[/yellow]")
            return [], []

        # ============================================================
        # STAGE 1: DISULFIDE BOND REVIEW (if any detected)
        # ============================================================
        potential_disulfide_pairs = self._identify_disulfide_pairs(centers)
        confirmed_disulfide_pairs = []

        if potential_disulfide_pairs:
            self.console.print("\n[bold cyan]═══ STAGE 1: Disulfide Bond Review ═══[/bold cyan]\n")
            confirmed_disulfide_pairs = self._review_disulfide_candidates(
                potential_disulfide_pairs, centers
            )
            # Update center properties based on confirmation
            self._update_disulfide_confirmations(centers, confirmed_disulfide_pairs)

        # ============================================================
        # STAGE 2: CENTER SELECTION (all centers with disulfide status)
        # ============================================================
        self.console.print("\n[bold cyan]═══ STAGE 2: Center Selection ═══[/bold cyan]\n")
        selected_centers = self._select_centers_with_disulfide_status(
            centers, confirmed_disulfide_pairs
        )

        return selected_centers, confirmed_disulfide_pairs

    def _identify_disulfide_pairs(self, centers: List[RedoxCenter]) -> List[tuple]:
        """Identify pairs of disulfide-bonded centers"""
        pairs = []
        processed = set()

        for center in centers:
            if center.properties.get('is_disulfide_bonded', False):
                center_key = (center.chain, center.resid)

                if center_key in processed:
                    continue

                # Find partner
                partner_chain = center.properties['disulfide_partner_chain']
                partner_res = center.properties['disulfide_partner_res']

                for other in centers:
                    if (other.chain == partner_chain and other.resid == partner_res and
                        other.properties.get('is_disulfide_bonded', False)):

                        pairs.append((center, other))
                        processed.add(center_key)
                        processed.add((other.chain, other.resid))
                        break

        return pairs

    def _review_disulfide_candidates(self, potential_pairs: List[tuple],
                                     all_centers: List[RedoxCenter]) -> List[tuple]:
        """Let user review and confirm which disulfide bonds to keep"""

        self.console.print(f"[yellow]Detected {len(potential_pairs)} potential disulfide bond(s)[/yellow]\n")
        self.console.print("[grey50]Review each bond and decide whether to treat it as a disulfide bond[/grey50]")
        self.console.print("[grey50]or as two independent CYS residues.[/grey50]\n")

        # Display all candidates in a table
        candidate_table = Table(title="Potential Disulfide Bonds")
        candidate_table.add_column("Bond ID", style="cyan", width=8)
        candidate_table.add_column("CYS 1", style="green", width=15)
        candidate_table.add_column("CYS 2", style="yellow", width=15)
        candidate_table.add_column("SG-SG Dist", style="magenta", width=12)
        candidate_table.add_column("Source", style="grey50", width=15)

        for i, (c1, c2) in enumerate(potential_pairs, 1):
            dist = c1.properties.get('disulfide_bond_distance', 0)
            source = "SSBOND" if c1.properties.get('from_ssbond_record', False) else "Distance"

            candidate_table.add_row(
                f"SS{i}",
                f"{c1.chain}:{c1.resid}",
                f"{c2.chain}:{c2.resid}",
                f"{dist:.2f}Å",
                source
            )

        self.console.print(candidate_table)
        self.console.print()

        # Selection options
        self.console.print("[bold]Disulfide Bond Selection Options:[/bold]")
        self.console.print("  [green]all[/green]       - Confirm all detected disulfide bonds")
        self.console.print("  [green]none[/green]      - Treat all CYS as independent (no disulfides)")
        self.console.print("  [green]review[/green]    - Review each bond individually")
        self.console.print("  [green]view[/green]      - View disulfide bonds in 3D structure viewer")
        self.console.print("  [green]SS1,SS3[/green]   - Select specific bonds (e.g., 'SS1,SS3,SS5')")
        self.console.print("  [green]SS1-SS3[/green]   - Select bond ranges (e.g., 'SS1-SS3' or 'SS1-3')")
        self.console.print("  [green]SS1,SS3-5[/green] - Mix individual and ranges")
        self.console.print()

        while True:  # Loop to allow viewing and returning to prompt
            choice = prompt_with_context(
                self.processor,
                "How to handle disulfide bonds?",
                default="all",
                module="Redox Site Detector",
                description="Disulfide bond selection strategy"
            )

            # Process choice
            if choice.lower() == 'all':
                self.console.print(f"[green]✓ Confirmed all {len(potential_pairs)} disulfide bond(s)[/green]")
                return potential_pairs
            elif choice.lower() == 'none':
                self.console.print(f"[yellow]All CYS will be treated as independent residues[/yellow]")
                return []
            elif choice.lower() == 'review':
                return self._review_each_disulfide(potential_pairs)
            elif choice.lower() == 'view':
                # Launch viewer with all CYS centers involved in disulfides
                centers_to_view = []
                for c1, c2 in potential_pairs:
                    if c1 not in centers_to_view:
                        centers_to_view.append(c1)
                    if c2 not in centers_to_view:
                        centers_to_view.append(c2)
                self._launch_quick_structure_viewer(centers_to_view, "Disulfide Bond Candidates")
                continue  # Return to prompt after viewing
            else:
                # Try to parse as specific bond selection
                result = self._parse_disulfide_selection(choice, potential_pairs)
                if result is not None:
                    return result
                # If parsing failed, loop will repeat prompt

    def _review_each_disulfide(self, potential_pairs: List[tuple]) -> List[tuple]:
        """Review each disulfide bond individually"""

        self.console.print("\n[bold]Reviewing each disulfide bond individually...[/bold]\n")

        confirmed_pairs = []

        for i, (c1, c2) in enumerate(potential_pairs, 1):
            dist = c1.properties.get('disulfide_bond_distance', 0)
            source = "SSBOND record" if c1.properties.get('from_ssbond_record', False) else "distance detection"

            self.console.print(f"[bold]Bond SS{i}:[/bold] {c1.chain}:{c1.resid} ←→ {c2.chain}:{c2.resid}")
            self.console.print(f"  Distance: {dist:.2f}Å")
            self.console.print(f"  Source: {source}")

            if dist > 2.5:
                self.console.print(f"  [yellow]⚠ Note: Distance exceeds typical threshold (2.5Å)[/yellow]")

            keep = confirm_with_context(
                self.processor,
                f"Treat as disulfide bond?",
                default=True,
                module="Redox Site Detector",
                description=f"Confirm disulfide bond SS{i}"
            )

            if keep:
                confirmed_pairs.append((c1, c2))
                self.console.print(f"[green]✓ Confirmed as disulfide bond[/green]\n")
            else:
                self.console.print(f"[grey50]○ Will be treated as independent CYS residues[/grey50]\n")

        self.console.print(f"[bold]Summary:[/bold] Confirmed {len(confirmed_pairs)} of {len(potential_pairs)} disulfide bond(s)\n")

        return confirmed_pairs

    def _parse_disulfide_selection(self, selection: str,
                                   potential_pairs: List[tuple]) -> List[tuple]:
        """Parse disulfide bond selection string like 'SS1,SS3,SS5' or 'SS1-SS3'"""
        confirmed_pairs = []

        try:
            # Parse comma-separated bond IDs (can include ranges)
            bond_tokens = [s.strip().upper() for s in selection.split(',')]

            for token in bond_tokens:
                # Check for range syntax (SS1-SS3)
                if '-' in token and token.startswith('SS'):
                    # Split range
                    range_parts = token.split('-')
                    if len(range_parts) == 2:
                        start_part = range_parts[0].strip()
                        end_part = range_parts[1].strip()

                        # Parse start
                        if start_part.startswith('SS'):
                            try:
                                start_idx = int(start_part[2:]) - 1
                            except ValueError:
                                self.console.print(f"[yellow]Warning: Invalid range start '{start_part}'[/yellow]")
                                continue
                        else:
                            self.console.print(f"[yellow]Warning: Range start should be in format 'SS1'[/yellow]")
                            continue

                        # Parse end (allow just number like SS1-3 or full SS1-SS3)
                        if end_part.startswith('SS'):
                            try:
                                end_idx = int(end_part[2:]) - 1
                            except ValueError:
                                self.console.print(f"[yellow]Warning: Invalid range end '{end_part}'[/yellow]")
                                continue
                        else:
                            # Just a number
                            try:
                                end_idx = int(end_part) - 1
                            except ValueError:
                                self.console.print(f"[yellow]Warning: Invalid range end '{end_part}'[/yellow]")
                                continue

                        # Add all bonds in range
                        if start_idx > end_idx:
                            self.console.print(f"[yellow]Warning: Invalid range {token} (start > end)[/yellow]")
                            continue

                        for idx in range(start_idx, end_idx + 1):
                            if 0 <= idx < len(potential_pairs):
                                if potential_pairs[idx] not in confirmed_pairs:
                                    confirmed_pairs.append(potential_pairs[idx])
                            else:
                                self.console.print(f"[yellow]Warning: SS{idx+1} is out of range (1-{len(potential_pairs)})[/yellow]")
                    else:
                        self.console.print(f"[yellow]Warning: Invalid range format '{token}'[/yellow]")

                # Single bond ID
                elif token.startswith('SS'):
                    try:
                        idx = int(token[2:]) - 1  # Convert SS1 → index 0

                        if 0 <= idx < len(potential_pairs):
                            if potential_pairs[idx] not in confirmed_pairs:
                                confirmed_pairs.append(potential_pairs[idx])
                        else:
                            self.console.print(f"[yellow]Warning: {token} is out of range (1-{len(potential_pairs)})[/yellow]")
                    except ValueError:
                        self.console.print(f"[yellow]Warning: Invalid bond ID format '{token}'[/yellow]")
                else:
                    self.console.print(f"[yellow]Warning: Bond IDs should be in format 'SS1', 'SS2', etc.[/yellow]")

            if confirmed_pairs:
                bond_list = ', '.join([f'SS{potential_pairs.index(p)+1}' for p in confirmed_pairs])
                self.console.print(f"[green]✓ Selected {len(confirmed_pairs)} disulfide bond(s): {bond_list}[/green]")
            else:
                self.console.print("[yellow]No valid disulfide bonds selected[/yellow]")

        except Exception as e:
            self.console.print(f"[red]Error parsing selection: {e}[/red]")
            self.console.print("[yellow]Treating as 'none' - no disulfides selected[/yellow]")
            return []

        return confirmed_pairs

    def _update_disulfide_confirmations(self, all_centers: List[RedoxCenter],
                                       confirmed_pairs: List[tuple]):
        """Update center properties based on which disulfides were confirmed"""
        confirmed_center_keys = set()
        for c1, c2 in confirmed_pairs:
            confirmed_center_keys.add((c1.chain, c1.resid))
            confirmed_center_keys.add((c2.chain, c2.resid))

        # Update all centers
        for center in all_centers:
            if center.properties.get('is_disulfide_bonded', False):
                center_key = (center.chain, center.resid)

                if center_key not in confirmed_center_keys:
                    # This was a candidate but user rejected it
                    center.properties['is_disulfide_bonded'] = False
                    center.properties['was_candidate_disulfide'] = True
                    center.properties['disulfide_rejected'] = True

                    # Keep the partner info for reference
                    center.properties['candidate_disulfide_partner_chain'] = center.properties.pop('disulfide_partner_chain', None)
                    center.properties['candidate_disulfide_partner_res'] = center.properties.pop('disulfide_partner_res', None)

    def _select_centers_with_disulfide_status(self, all_centers: List[RedoxCenter],
                                              confirmed_disulfide_pairs: List[tuple]) -> List[RedoxCenter]:
        """Select centers with awareness of confirmed disulfide bonds"""

        # AUTOMATICALLY select all confirmed disulfide centers
        auto_selected_centers = []
        for c1, c2 in confirmed_disulfide_pairs:
            if c1 not in auto_selected_centers:
                auto_selected_centers.append(c1)
            if c2 not in auto_selected_centers:
                auto_selected_centers.append(c2)

        # Get remaining centers (not part of confirmed disulfides)
        remaining_centers = [c for c in all_centers if c not in auto_selected_centers]

        # Count different types
        disulfide_count = len(confirmed_disulfide_pairs)
        disulfide_center_count = len(auto_selected_centers)
        rejected_disulfide_count = len([c for c in remaining_centers if c.properties.get('disulfide_rejected', False)])
        other_center_count = len(remaining_centers) - rejected_disulfide_count

        # Display what's automatically selected
        if auto_selected_centers:
            auto_panel = Panel(
                f"[bold green]✓ Automatically Selected:[/bold green]\n\n"
                f"  {disulfide_count} confirmed disulfide bond(s) = {disulfide_center_count} CYS centers\n"
                f"  [grey50]These centers are included automatically[/grey50]",
                border_style="green",
                padding=(1, 2),
                expand=False
            )
            self.console.print(auto_panel)
            self.console.print()

        # If there are remaining centers, let user select from them
        if remaining_centers:
            # Display remaining centers summary
            summary_panel = Panel(
                f"[bold]Remaining Centers to Select:[/bold]\n\n"
                f"  [yellow]○ {rejected_disulfide_count} independent CYS[/yellow] (rejected disulfides)\n"
                f"  [cyan]• {other_center_count} other redox center(s)[/cyan]\n"
                f"\n[grey50]Total remaining: {len(remaining_centers)} centers[/grey50]",
                border_style="cyan",
                padding=(1, 2),
                expand=False
            )
            self.console.print(summary_panel)
            self.console.print()

            # Display remaining centers table
            display_table = Table(title="Remaining Centers to Select")
            display_table.add_column("ID", style="cyan", width=6)
            display_table.add_column("Location", style="yellow", width=20)
            display_table.add_column("Element", style="green", width=8)
            display_table.add_column("Type", style="magenta", width=25)

            # Add rows for remaining centers only
            for i, center in enumerate(remaining_centers, 1):
                location = f"{center.resname} {center.chain}:{center.resid}"
                if center.altloc:
                    location += center.altloc
                if center.atom_name:
                    location += f" {center.atom_name}"

                element = center.element or "—"
                center_type = center.center_type.value.replace('_', ' ') if center.center_type else "Unknown"

                display_table.add_row(str(i), location, element, center_type)

            self.console.print(display_table)
            self.console.print()

            # Selection options for remaining centers
            self.console.print("[bold]Selection Options:[/bold]")
            self.console.print("  [grey50]•[/grey50] [green]all[/green]      - Select all remaining centers")
            self.console.print("  [grey50]•[/grey50] [green]none[/green]     - Select only the disulfide bonds (no others)")
            self.console.print("  [grey50]•[/grey50] [green]view[/green]     - View remaining centers in 3D structure viewer")
            self.console.print("  [grey50]•[/grey50] [green]1-10[/green]     - Select by ID range")
            self.console.print("  [grey50]•[/grey50] [green]1,3,5-8[/green]  - Select specific IDs")
            self.console.print()

            while True:  # Loop to allow viewing and returning to prompt
                selection = prompt_with_context(
                    self.processor,
                    "Select additional centers to include",
                    default="all",
                    module="Redox Site Detector",
                    description="Select additional centers beyond disulfides"
                )

                # Parse selection for remaining centers
                if selection.lower() == 'all':
                    selected_from_remaining = remaining_centers.copy()
                    break
                elif selection.lower() == 'none':
                    selected_from_remaining = []
                    break
                elif selection.lower() == 'view':
                    # Launch viewer with remaining centers
                    self._launch_quick_structure_viewer(remaining_centers, "Remaining Centers")
                    continue  # Return to prompt after viewing
                else:
                    # Parse ID selection
                    indices = self._parse_selection_input(selection, len(remaining_centers))
                    if indices is not None:
                        selected_from_remaining = [remaining_centers[i] for i in indices]
                        break
                    else:
                        self.console.print("[yellow]Invalid selection. Please try again.[/yellow]")
                        continue  # Repeat prompt for valid input

        else:
            # No remaining centers
            self.console.print("[grey50]No additional centers to select.[/grey50]\n")
            selected_from_remaining = []

        # Combine auto-selected disulfides with user-selected remaining
        final_selected = auto_selected_centers + selected_from_remaining

        # Display final selection summary
        selected_disulfide_count = len(confirmed_disulfide_pairs)
        selected_other_count = len(selected_from_remaining)

        self.console.print(f"\n[bold green]✓ Final Selection: {len(final_selected)} centers total[/bold green]")
        if selected_disulfide_count > 0:
            self.console.print(f"  • {selected_disulfide_count} disulfide bond(s) = {len(auto_selected_centers)} centers [green](auto-selected)[/green]")
        if selected_other_count > 0:
            self.console.print(f"  • {selected_other_count} additional center(s)")

        # Re-emit highlights so the viewer drops the un-selected centers
        # and shows only what the user kept. Each surviving center keeps
        # its original inventory color (looked up by position in
        # ``all_centers``) so the table-to-viewer mapping the user just
        # used to make the choice stays consistent. Nearby residues are
        # preserved for the kept centers, dropped for the rest.
        nearby_bonds = getattr(self, 'inventory_nearby_bonds', {}) or {}
        selected_items = []
        for c in final_selected:
            try:
                inv_idx = all_centers.index(c) + 1
            except ValueError:
                continue
            selected_items.append({
                "selection": _center_to_ngl(c),
                "label": f"selected_{inv_idx}_center",
                "color_index": inv_idx,
            })
            nearby_residues = sorted({
                (atom['chain'], atom['resid'])
                for atom in nearby_bonds.get((c.chain, c.resid, c.insertion_code), [])
            })
            if nearby_residues:
                selected_items.append({
                    "selection": " or ".join(f":{ch} and {rid}" for ch, rid in nearby_residues),
                    "label": f"selected_{inv_idx}_nearby",
                    "color_index": inv_idx,
                })
        _auto_show_palette_highlights(getattr(self, 'processor', None), selected_items)

        return final_selected

    def _parse_selection_input(self, input_str: str, max_num: int) -> Optional[List[int]]:
        """Parse selection input supporting ranges and comma-separated values
        
        Examples:
        - "1-10" -> [0,1,2,3,4,5,6,7,8,9]
        - "1,3,5-8" -> [0,2,4,5,6,7]
        - "1-5 10-15" -> [0,1,2,3,4,9,10,11,12,13,14]
        """
        try:
            indices = set()
            
            # Split by both spaces and commas
            parts = input_str.replace(',', ' ').split()
            
            for part in parts:
                if '-' in part:
                    # Handle range (e.g., "1-10")
                    range_parts = part.split('-')
                    if len(range_parts) == 2:
                        start = int(range_parts[0])
                        end = int(range_parts[1])
                        if start <= end:
                            for i in range(start, end + 1):
                                if 1 <= i <= max_num:
                                    indices.add(i - 1)  # Convert to 0-based
                                else:
                                    return None  # Invalid range
                        else:
                            return None  # Invalid range (start > end)
                    else:
                        return None  # Invalid range format
                else:
                    # Handle single number
                    num = int(part)
                    if 1 <= num <= max_num:
                        indices.add(num - 1)  # Convert to 0-based
                    else:
                        return None  # Invalid number
            
            return sorted(list(indices))

        except ValueError:
            return None

    # ========================================================================
    # Structure Viewer Integration
    # ========================================================================

    def _launch_quick_structure_viewer(self, centers_to_highlight: List[RedoxCenter],
                                      title: str = "Selected Centers") -> bool:
        """
        Launch structure viewer with specific centers highlighted.

        Args:
            centers_to_highlight: Centers to highlight in the viewer
            title: Title for the annotation

        Returns:
            True if viewer launched successfully
        """
        ngl_selection = self._centers_to_ngl_selection(centers_to_highlight)
        return _show_residues_in_viewer(
            getattr(self, 'processor', None),
            self.console,
            ngl_selection,
            label='selected_centers',
        )

    def _centers_to_ngl_selection(self, centers: List[RedoxCenter]) -> str:
        """
        Convert RedoxCenter list to NGL selection string.

        Args:
            centers: List of RedoxCenter objects

        Returns:
            NGL selection string like ":A and 154 or :A and 96"
        """
        selections = []
        for center in centers:
            selection = f":{center.chain} and {center.resid}"

            # Add altloc if present (and not blank/empty)
            if center.altloc and center.altloc.strip():
                selection += f" and %{center.altloc.strip()}"

            selections.append(selection)

        return " or ".join(selections) if selections else ""

class CenterGroupingInterface:
    """Interface for grouping centers into redox sites"""

    def __init__(self, config: DetectionConfig, console: Console = None, processor=None):
        self.config = config
        self.console = console if console else Console()
        self.processor = processor
    
    def group_centers_into_sites(self, selected_centers: List[RedoxCenter],
                                confirmed_disulfide_pairs: List[tuple] = None) -> List[List[RedoxCenter]]:
        """
        Group selected centers into redox sites with automatic disulfide pairing

        Args:
            selected_centers: Centers selected by user
            confirmed_disulfide_pairs: Disulfide pairs confirmed in selection stage

        Returns:
            List of site groups (each group is a list of RedoxCenters)
        """
        if not selected_centers:
            return []

        sites = []
        remaining_centers = selected_centers.copy()
        site_number = 1

        self.console.print("\n[bold underline]Grouping Into Redox Sites[/bold underline]")

        # ===== AUTOMATIC DISULFIDE GROUPING (no prompt) =====
        if confirmed_disulfide_pairs:
            # Filter confirmed pairs to only those in selected_centers
            selected_disulfide_pairs = []
            for c1, c2 in confirmed_disulfide_pairs:
                if c1 in selected_centers and c2 in selected_centers:
                    selected_disulfide_pairs.append((c1, c2))

            if selected_disulfide_pairs:
                self.console.print(f"\n[bold cyan]═══ Disulfide Bond Grouping ═══[/bold cyan]\n")

                # Automatically create one site per disulfide pair
                for i, (c1, c2) in enumerate(selected_disulfide_pairs, 1):
                    sites.append([c1, c2])
                    remaining_centers.remove(c1)
                    remaining_centers.remove(c2)

                self.console.print(f"[green]✓ Automatically created {len(selected_disulfide_pairs)} disulfide site(s)[/green]")
                self.console.print("[grey50]  Each disulfide bond grouped as a separate site[/grey50]\n")
                site_number += len(selected_disulfide_pairs)

        # ===== CONTINUE WITH STANDARD AUTO-GROUPING =====
        # Check for auto-grouping opportunities
        auto_groups = self._find_auto_grouping_candidates(remaining_centers)
        if auto_groups:
            # Display auto-grouping candidates in table
            auto_table = Table(title=f"Found {len(auto_groups)} Potential Auto-Groupings")
            auto_table.add_column("ID", style="cyan", width=8)
            auto_table.add_column("Location", style="green", width=20)
            auto_table.add_column("Centers", style="yellow", width=15)
            auto_table.add_column("Elements", style="magenta")
            
            for i, (key, centers) in enumerate(auto_groups.items(), 1):
                resname, chain, resid = key
                elements = []
                for center in centers:
                    if center.atom_name:
                        elements.append(f"{center.element} {center.atom_name}")
                    else:
                        elements.append(f"{center.element}")
                
                auto_table.add_row(
                    f"[{i}]",
                    f"{resname} {chain}:{resid}",
                    f"{len(centers)} centers",
                    ", ".join(elements)
                )
            
            self.console.print(auto_table)
            
            auto_choice = confirm_with_context(
                processor=self.processor,
                prompt="\n[bold cyan]Auto-group centers with same chain:resid:resname?[/bold cyan]",
                default=True,
                module="Redox Detector",
                description="Auto-group centers by location"
            )
            if auto_choice:
                # Create sites from auto-groups with Rich display
                created_table = Table(title="Auto-Created Sites")
                created_table.add_column("Site", style="cyan")
                created_table.add_column("Location", style="green")
                created_table.add_column("Centers", style="yellow")
                
                for key, centers in auto_groups.items():
                    sites.append(centers)
                    for center in centers:
                        remaining_centers.remove(center)
                    resname, chain, resid = key
                    
                    created_table.add_row(
                        f"Site {site_number}",
                        f"{resname} {chain}:{resid}",
                        f"{len(centers)} centers"
                    )
                    site_number += 1
                
                self.console.print(created_table)
        
        # New comprehensive grouping interface
        if remaining_centers:
            self.console.print(f"\n[bold cyan]═══ Manual Site Grouping ═══[/bold cyan]\n")

            # Display all remaining centers for grouping
            centers_table = Table(title=f"Centers to Group ({len(remaining_centers)} total)")
            centers_table.add_column("ID", style="cyan", width=8)
            centers_table.add_column("Location", style="green", width=25)
            centers_table.add_column("Element", style="yellow", width=10)
            centers_table.add_column("Type", style="magenta", no_wrap=True)
            
            for i, center in enumerate(remaining_centers, 1):
                type_str = center.center_type.value.replace('_', ' ')
                location = f"{center.resname} {center.chain}:{center.resid}"
                if center.altloc:
                    location += center.altloc
                if center.atom_name:
                    location += f" {center.atom_name}"
                
                centers_table.add_row(
                    f"[{i}]",
                    location,
                    center.element or "?",
                    type_str
                )
            
            self.console.print(f"\n")
            self.console.print(centers_table)
            
            # Comprehensive site assignment input
            while True:
                from rich.text import Text

                # Create format text without markup interpretation
                format_text = Text()
                format_text.append("\n")
                format_text.append("Enter site assignments", style="bold cyan")
                format_text.append(" (format: site:centers site:centers)")
                self.console.print(format_text)

                self.console.print("[bold]Examples:[/bold]")
                self.console.print("  [grey50]•[/grey50] [green]1:1 2:2 3:3 4:4[/green]     [grey50](individual sites)[/grey50]")
                self.console.print("  [grey50]•[/grey50] [green]1:1,2 2:3,4[/green]        [grey50](paired sites)[/grey50]")
                self.console.print("  [grey50]•[/grey50] [green]1:1-4[/green]              [grey50](single large site)[/grey50]")
                self.console.print("  [grey50]•[/grey50] [green]1:1,3 2:2,4-6[/green]      [grey50](mixed grouping)[/grey50]")
                self.console.print("  [grey50]•[/grey50] [green]\\[Enter][/green]            [grey50](auto-assign individual sites)[/grey50]")
                
                assignment_input = prompt_with_context(
                    processor=self.processor,
                    prompt="[green]Site assignments[/green]",
                    default="",
                    module="Redox Detector",
                    description="Enter site assignments (e.g., 1:1,2 2:3,4)"
                ).strip()
                
                if assignment_input == '':
                    # Auto-assign each center to its own site (1:1 2:2 3:3 etc.)
                    num_centers = len(remaining_centers)
                    auto_assignments = ' '.join([f"{i+1}:{i+1}" for i in range(num_centers)])
                    self.console.print(f"[cyan]Auto-assigning individual sites: {auto_assignments}[/cyan]")
                    assignment_input = auto_assignments
                
                # Parse and validate assignments
                parsed_sites = self._parse_site_assignments(assignment_input, len(remaining_centers))
                
                if parsed_sites is None:
                    self.console.print("[red]Invalid assignment format or values.[/red]")
                    continue
                
                # Validate complete 1-to-1 mapping
                validation_result = self._validate_site_assignments(parsed_sites, len(remaining_centers))
                
                if validation_result is not True:
                    self.console.print(f"[red]Validation error: {validation_result}[/red]")
                    continue
                
                # Create preview sites from assignments
                preview_sites = self._create_sites_from_assignments(parsed_sites, remaining_centers)
                
                # Display preview and get confirmation
                preview_table = Table(title=f"Proposed Site Grouping ({len(preview_sites)} sites)")
                preview_table.add_column("Site", style="cyan", width=10)
                preview_table.add_column("Centers", style="green")
                preview_table.add_column("Type", style="yellow")
                
                for i, site_centers in enumerate(preview_sites, 1):
                    site_type = self._get_site_type_description(site_centers)
                    
                    if len(site_centers) == 1:
                        center = site_centers[0]
                        location = f"{center.resname} {center.chain}:{center.resid}"
                        if center.atom_name:
                            location += f" {center.atom_name}"
                        centers_str = location
                    else:
                        center_strs = []
                        for center in site_centers:
                            location = f"{center.resname} {center.chain}:{center.resid}"
                            if center.atom_name:
                                location += f" {center.atom_name}"
                            center_strs.append(location)
                        centers_str = ' + '.join(center_strs)
                    
                    preview_table.add_row(
                        f"Site {i}",
                        centers_str,
                        site_type
                    )
                
                self.console.print("\n")
                self.console.print(preview_table)

                sites.extend(preview_sites)
                remaining_centers.clear()  # All centers are now assigned
                break
        
        # Display final confirmed site grouping
        if sites:
            final_table = Table(title=f"Final Site Grouping ({len(sites)} sites)")
            final_table.add_column("Site", style="cyan", width=10)
            final_table.add_column("Centers", style="green")
            final_table.add_column("Type", style="yellow")
            
            for i, site_centers in enumerate(sites, 1):
                site_type = self._get_site_type_description(site_centers)
                
                if len(site_centers) == 1:
                    center = site_centers[0]
                    location = f"{center.resname} {center.chain}:{center.resid}"
                    if center.atom_name:
                        location += f" {center.atom_name}"
                    centers_str = location
                else:
                    center_strs = []
                    for center in site_centers:
                        location = f"{center.resname} {center.chain}:{center.resid}"
                        if center.atom_name:
                            location += f" {center.atom_name}"
                        center_strs.append(location)
                    centers_str = ' + '.join(center_strs)
                
                final_table.add_row(
                    f"Site {i}",
                    centers_str,
                    site_type
                )
            
            self.console.print("\n")
            self.console.print(final_table)
        
        return sites
    
    def _find_auto_grouping_candidates(self, centers: List[RedoxCenter]) -> Dict[Tuple[str, str, int], List[RedoxCenter]]:
        """Find centers that can be auto-grouped by same chain:resid:resname"""
        groups = {}
        
        for center in centers:
            key = (center.resname, center.chain, center.resid)
            if key not in groups:
                groups[key] = []
            groups[key].append(center)
        
        # Only return groups with multiple centers
        auto_groups = {k: v for k, v in groups.items() if len(v) > 1}
        return auto_groups
    
    def _get_site_type_description(self, site_centers: List[RedoxCenter]) -> str:
        """Get a descriptive name for the site type"""
        if len(site_centers) == 1:
            center = site_centers[0]
            if center.center_type == CenterType.METAL_ION:
                return f"mono-nuclear {center.element}"
            elif center.center_type == CenterType.ORGANOMETALLIC_COFACTOR:
                return f"organometallic cofactor ({center.resname})"
            elif center.center_type == CenterType.ORGANIC_COFACTOR:
                return "organic cofactor"
            elif center.center_type == CenterType.REDOX_AMINO_ACID:
                if center.resname == 'CYS':
                    return "free cysteine"
                else:
                    return f"redox-active {center.resname}"
        else:
            # Multi-center site
            metals = [c for c in site_centers if c.center_type == CenterType.METAL_ION]
            organometallics = [c for c in site_centers if c.center_type == CenterType.ORGANOMETALLIC_COFACTOR]
            organics = [c for c in site_centers if c.center_type == CenterType.ORGANIC_COFACTOR]
            amino_acids = [c for c in site_centers if c.center_type == CenterType.REDOX_AMINO_ACID]

            if len(metals) > 1:
                if all(m.element == metals[0].element for m in metals):
                    return f"multi-nuclear {metals[0].element} cluster"
                else:
                    return "multi-metal cluster"
            elif organometallics:
                return f"organometallic site ({organometallics[0].resname})"
            elif len(metals) == 1 and organics:
                return f"{metals[0].element} in {organics[0].resname}"
            elif len(amino_acids) == 2 and all(c.resname == 'CYS' for c in amino_acids):
                return "disulfide pair"
            else:
                return "multi-center site"
    
    def _parse_selection_input(self, input_str: str, max_num: int) -> Optional[List[int]]:
        """Parse selection input supporting ranges and comma-separated values
        
        Examples:
        - "1-10" -> [0,1,2,3,4,5,6,7,8,9]
        - "1,3,5-8" -> [0,2,4,5,6,7]
        - "1-5 10-15" -> [0,1,2,3,4,9,10,11,12,13,14]
        """
        try:
            indices = set()
            
            # Split by both spaces and commas
            parts = input_str.replace(',', ' ').split()
            
            for part in parts:
                if '-' in part:
                    # Handle range (e.g., "1-10")
                    range_parts = part.split('-')
                    if len(range_parts) == 2:
                        start = int(range_parts[0])
                        end = int(range_parts[1])
                        if start <= end:
                            for i in range(start, end + 1):
                                if 1 <= i <= max_num:
                                    indices.add(i - 1)  # Convert to 0-based
                                else:
                                    return None  # Invalid range
                        else:
                            return None  # Invalid range (start > end)
                    else:
                        return None  # Invalid range format
                else:
                    # Handle single number
                    num = int(part)
                    if 1 <= num <= max_num:
                        indices.add(num - 1)  # Convert to 0-based
                    else:
                        return None  # Invalid number
            
            return sorted(list(indices))
            
        except ValueError:
            return None
    
    def _parse_site_assignments(self, assignment_str: str, max_centers: int) -> Optional[Dict[int, List[int]]]:
        """
        Parse site assignment string into site_id -> center_indices mapping.
        
        Format: "site:centers site:centers"
        Examples: "1:1,2 2:3,4" or "1:1-4" or "1:1 2:2 3:3"
        
        Returns:
            Dict mapping site_id to list of 0-based center indices, or None if invalid
        """
        try:
            sites = {}
            
            # Split by spaces to get individual site assignments
            assignments = assignment_str.split()
            
            for assignment in assignments:
                if ':' not in assignment:
                    return None
                
                # Split site_id:center_spec
                parts = assignment.split(':', 1)
                if len(parts) != 2:
                    return None
                
                site_id_str, center_spec = parts
                
                # Parse site ID
                try:
                    site_id = int(site_id_str)
                    if site_id < 1:
                        return None
                except ValueError:
                    return None
                
                # Parse center specification (same logic as _parse_selection_input)
                center_indices = self._parse_selection_input(center_spec, max_centers)
                if center_indices is None:
                    return None
                
                # Check for duplicate site ID
                if site_id in sites:
                    return None  # Duplicate site ID
                
                sites[site_id] = center_indices
            
            return sites
            
        except Exception:
            return None
    
    def _validate_site_assignments(self, sites: Dict[int, List[int]], total_centers: int) -> Union[bool, str]:
        """
        Validate that site assignments form a complete 1-to-1 mapping.
        
        Returns:
            True if valid, or error message string if invalid
        """
        if not sites:
            return "No site assignments provided"
        
        # Check site IDs are sequential starting from 1
        site_ids = sorted(sites.keys())
        expected_ids = list(range(1, len(site_ids) + 1))
        if site_ids != expected_ids:
            return f"Site IDs must be sequential starting from 1. Got: {site_ids}, expected: {expected_ids}"
        
        # Collect all assigned center indices
        all_assigned = []
        for site_id, center_indices in sites.items():
            if not center_indices:
                return f"Site {site_id} has no centers assigned"
            all_assigned.extend(center_indices)
        
        # Check for duplicates
        if len(all_assigned) != len(set(all_assigned)):
            duplicates = [x for x in set(all_assigned) if all_assigned.count(x) > 1]
            duplicate_centers = [str(x + 1) for x in duplicates]  # Convert back to 1-based
            return f"Centers assigned to multiple sites: {', '.join(duplicate_centers)}"
        
        # Check all centers are assigned
        expected_centers = set(range(total_centers))
        assigned_centers = set(all_assigned)
        
        if assigned_centers != expected_centers:
            missing = expected_centers - assigned_centers
            extra = assigned_centers - expected_centers
            
            if missing:
                missing_1based = [str(x + 1) for x in sorted(missing)]
                return f"Centers not assigned to any site: {', '.join(missing_1based)}"
            
            if extra:
                extra_1based = [str(x + 1) for x in sorted(extra)]
                return f"Invalid center IDs: {', '.join(extra_1based)}"
        
        return True
    
    def _create_sites_from_assignments(self, sites: Dict[int, List[int]], centers: List[RedoxCenter]) -> List[List[RedoxCenter]]:
        """
        Create site groups from validated assignments.
        
        Returns:
            List of site groups (each group is a list of RedoxCenter objects)
        """
        # Sort by site ID to ensure consistent ordering
        sorted_sites = sorted(sites.items())
        
        site_groups = []
        for site_id, center_indices in sorted_sites:
            site_centers = [centers[i] for i in center_indices]
            site_groups.append(site_centers)
        
        return site_groups
    
    def _get_site_type_description(self, site_centers: List[RedoxCenter]) -> str:
        """Get a descriptive name for the site type"""
        if len(site_centers) == 1:
            center = site_centers[0]
            if center.center_type == CenterType.METAL_ION:
                return f"mono-nuclear {center.element}"
            elif center.center_type == CenterType.ORGANOMETALLIC_COFACTOR:
                return f"organometallic cofactor ({center.resname})"
            elif center.center_type == CenterType.ORGANIC_COFACTOR:
                return "organic cofactor"
            elif center.center_type == CenterType.REDOX_AMINO_ACID:
                if center.resname == 'CYS':
                    return "free cysteine"
                else:
                    return f"redox-active {center.resname}"
        else:
            # Multi-center site
            metals = [c for c in site_centers if c.center_type == CenterType.METAL_ION]
            organometallics = [c for c in site_centers if c.center_type == CenterType.ORGANOMETALLIC_COFACTOR]
            organics = [c for c in site_centers if c.center_type == CenterType.ORGANIC_COFACTOR]
            amino_acids = [c for c in site_centers if c.center_type == CenterType.REDOX_AMINO_ACID]

            if len(metals) > 1:
                if all(m.element == metals[0].element for m in metals):
                    return f"multi-nuclear {metals[0].element} cluster"
                else:
                    return "multi-metal cluster"
            elif organometallics:
                return f"organometallic site ({organometallics[0].resname})"
            elif len(metals) == 1 and organics:
                return f"{metals[0].element} in {organics[0].resname}"
            elif len(amino_acids) == 2 and all(c.resname == 'CYS' for c in amino_acids):
                return "disulfide pair"
            else:
                return "multi-center site"

# ===== SITE REFINEMENT INTERFACE =====

class SiteRefinementInterface:
    """Complete interface for iterative site refinement with all 9 search method combinations"""

    def __init__(self, config: DetectionConfig, console: Console = None, template_mode: bool = False, transformer_mappings: dict = None, processor=None):
        self.config = config
        self.console = console if console else Console()
        self.template_mode = template_mode
        self.current_template = None
        self.processor = processor
        self.transformer_mappings = transformer_mappings or {}
    
    def refine_site_interactively(self, site: RedoxSite, structure: Structure) -> RedoxSite:
        """
        Complete iterative refinement interface for a redox site.
        Implements ALL 9 search method combinations:
        - Distance cutoff: Fixed radius + Adaptive radius  
        - Boundary definitions: Geometric center + Min distance (all atoms) + Min distance (custom atoms)
        - Count cutoff: Residue-based (find N closest residues of specified types)
        """
        
        self.console.print(f"\n[bold underline]REFINING SITE: {site.site_id}[/bold underline]")
        self.console.print(f"Current site contains [cyan]{len(site.centers)}[/cyan] center(s), [cyan]{len(site.atoms)}[/cyan] atoms")

        # Narrow the docked viewer to just this site. The same hook fires
        # for templated subsequent sites in ``apply_template_to_site``.
        _narrow_viewer_to_site(getattr(self, 'processor', None), site)
        
        # Show current centers in table
        centers_table = Table(title="Current Centers")
        centers_table.add_column("#", style="cyan", width=6)
        centers_table.add_column("Location", style="green")
        centers_table.add_column("Type", style="yellow")
        
        for i, center in enumerate(site.centers, 1):
            location = f"{center.resname} {center.chain}:{center.resid}"
            if center.atom_name:
                location += f" {center.atom_name}"
            type_str = center.center_type.value.replace('_', ' ')
            centers_table.add_row(str(i), location, type_str)
        
        self.console.print(centers_table)
        
        refined_site = copy.deepcopy(site)
        search_iteration = 1
        
        # Main refinement loop with parameter adjustment capability
        current_search_method = None
        current_search_params = None
        
        while True:
            # Check if we need to configure new search method or adjust existing parameters
            if current_search_method is None:
                self.console.print(f"\n[bold cyan]Search Iteration {search_iteration}[/bold cyan]")
                
                self.console.print("[grey50]1.[/grey50] Distance cutoff (search within radius)")
                self.console.print("[grey50]2.[/grey50] Residue count cutoff (find N closest residues)")
                self.console.print("[grey50]3.[/grey50] Finish refinement")
                if self.template_mode:
                    self.console.print(
                        "[yellow]Caution: Distance cutoffs may not be satisfied at every site "
                        "due to geometric variations between sites. Residue count is generally "
                        "more reliable for templates.[/yellow]"
                    )

                constraint_choice = prompt_with_context(
                    processor=self.processor,
                    prompt="[green]Choose search constraint[/green]",
                    choices=["1", "2", "3"],
                    default="2" if self.template_mode else "1",
                    module="Redox Detector - Refinement",
                    description="Select search constraint method",
                    options_map={
                        "1": "Distance cutoff (search within radius)",
                        "2": "Residue count cutoff (find N closest residues)",
                        "3": "Finish refinement"
                    }
                ).strip()
                
                if constraint_choice == '3':
                    self._finalize_refinement(refined_site)
                    break
                elif constraint_choice == '1':
                    current_search_method = "distance"
                elif constraint_choice == '2':
                    current_search_method = "count"
                # Note: Invalid choices handled by Prompt.ask with choices parameter
            
            # Configure or reconfigure parameters
            if current_search_method == "distance":
                search_params = self._configure_distance_search(refined_site)
                if not search_params:
                    current_search_method = None  # Reset to method selection
                    continue

                current_search_params = search_params
                
                # Perform search
                search_result = self._perform_distance_search(refined_site, structure, search_params)
                
            elif current_search_method == "count":
                search_params = self._configure_count_search(refined_site)
                if not search_params:
                    current_search_method = None  # Reset to method selection
                    continue
                current_search_params = search_params

                # Perform search
                search_result = self._perform_count_search(refined_site, structure, search_params)
            
            # Process results with action choice
            refined_site, action = self._process_search_results(refined_site, search_result, structure)
            
            if action == "finish":
                break
            elif action == "finish_refine":
                # User chose "finish" at the residue-selection prompt without
                # adding any found residues. Finalize like the "3. Finish
                # refinement" menu choice (bond definition) and exit the loop.
                self._finalize_refinement(refined_site)
                break
            elif action == "adjust_parameters":
                # Stay with same search method, reconfigure parameters in next iteration
                self.console.print(f"\n[yellow]Adjusting {current_search_method} search parameters...[/yellow]")
                # current_search_method stays the same, will reconfigure parameters
            elif action == "continue":
                # Continue with new search method selection
                current_search_method = None
                search_iteration += 1
            else:
                # Fallback
                current_search_method = None
                search_iteration += 1
        
        self.console.print(f"\n[bold green]Refinement complete![/bold green]")
        self.console.print(f"Final site contains [cyan]{len(refined_site.centers)}[/cyan] center(s), [cyan]{len(refined_site.atoms)}[/cyan] atoms, [cyan]{len(refined_site.bonds)}[/cyan] bonds")
        
        return refined_site

    def _finalize_refinement(self, refined_site: RedoxSite) -> None:
        """Finalize a refinement session and prepare to exit the loop.

        Shared by the "3. Finish refinement" menu choice and the "finish"
        command in the residue-selection prompt, so a site finishes identically
        no matter which way the user exits. In template mode, ensure a (minimal)
        template exists. Then route through the same define/skip-bond prompt a
        completed distance/count search reaches via _process_search_results: a
        site that needs no further searching -- e.g. a disulfide pair, where
        both SG atoms are already present -- otherwise exits without ever being
        offered bond definition, so the SG-SG `bond` line never makes it into
        the tleap script. Guard on having no bonds yet so finishing after an
        earlier iteration that already defined bonds doesn't clear/clobber them
        (the prompt clears site.bonds on entry).
        """
        if self.template_mode and self.current_template is None:
            self.current_template = SiteTemplate(
                site_type="",  # Will be set later
                residue_types=[],  # No additional residues
                residue_counts=[],  # No additional residues
                atom_filtering_metals=[],
                atom_filtering_nonmetals=[],
                boundary_definition=BoundaryDefinition.GEOMETRIC_CENTER,
                custom_boundary_atom_indices=None,
                residue_selection_choice="none",  # No residues added
                bond_pairs=[],  # No bonds
                bond_selections={},  # No bonds
                continue_searching=False  # Refinement finished
            )
        if not refined_site.bonds:
            self._define_bonds_interactively(refined_site)

    def _configure_distance_search(self, site: RedoxSite) -> Optional[SearchParameters]:
        """Configure distance-based search parameters with back navigation"""
        self.console.print(f"\n[bold underline]Distance-Based Search Configuration[/bold underline]")

        # Get previous search parameters if available
        previous_search_metals = None
        previous_search_nonmetals = None
        if site.search_history:
            last_search = site.search_history[-1]['search_parameters']
            previous_search_metals = last_search.search_include_metals
            previous_search_nonmetals = last_search.search_include_nonmetals

        # State tracking for navigation with conditional branching
        answers = {}
        step = 0

        while True:
            try:
                if step == 0:
                    # Question 1: Boundary definition
                    self.console.print("\n[bold]Boundary definition:[/bold]")
                    self.console.print("[grey50]1.[/grey50] Geometric center (calculate centroid and search from there)")
                    self.console.print("[grey50]2.[/grey50] Minimum distance from any site atoms")
                    self.console.print("[grey50]3.[/grey50] Minimum distance from custom selected atoms")

                    boundary_choice = prompt_with_navigation(
                        processor=self.processor,
                        prompt="[green]Select boundary definition[/green]",
                        choices=["1", "2", "3"],
                        default=answers.get('boundary_choice', "1"),
                        module="Redox Detector - Refinement",
                        description="Select boundary definition for distance search",
                        options_map={
                            "1": "Geometric center",
                            "2": "Minimum distance from any site atoms",
                            "3": "Minimum distance from custom selected atoms"
                        },
                        allow_back=False  # First question
                    ).strip()
                    answers['boundary_choice'] = boundary_choice

                    # Modes 1/2 get their indicator immediately; mode 3
                    # waits until the user picks atoms (handled in the
                    # 'custom' step) so the halo lands on the actual
                    # picked subset rather than nothing.
                    if boundary_choice in ('1', '2'):
                        _show_boundary_indicator(
                            getattr(self, 'processor', None), site, boundary_choice,
                        )

                    # If custom boundary, immediately ask for atom selection
                    if boundary_choice == '3':
                        step = 'custom'
                    else:
                        step = 1

                elif step == 'custom':
                    # Custom atom selection - immediately after choosing custom boundary
                    result = self._select_custom_boundary_atoms(site)
                    if not result:
                        self.console.print("[yellow]No custom atoms selected. Returning to boundary selection.[/yellow]")
                        step = 0
                        continue
                    answers['custom_coords'], answers['custom_indices'] = result
                    _show_boundary_indicator(
                        getattr(self, 'processor', None), site, '3',
                        custom_indices=answers['custom_indices'],
                    )
                    step = 1  # Proceed to distance method

                elif step == 1:
                    # Question 2: Distance method
                    self.console.print("\n[bold]Distance method:[/bold]")
                    self.console.print("[grey50]1.[/grey50] Fixed radius (search within specified distance)")
                    self.console.print("[grey50]2.[/grey50] Adaptive radius (adjust radius to achieve target coordination)")

                    method_choice = prompt_with_navigation(
                        processor=self.processor,
                        prompt="[green]Select distance method[/green]",
                        choices=["1", "2"],
                        default=answers.get('method_choice', "1"),
                        module="Redox Detector - Refinement",
                        description="Select distance search method",
                        options_map={
                            "1": "Fixed radius",
                            "2": "Adaptive radius"
                        },
                        allow_back=True
                    ).strip()
                    answers['method_choice'] = method_choice

                    # Determine next step based on method choice
                    if method_choice == '1':  # Fixed
                        step = 2  # Go to fixed radius prompt
                    else:  # Adaptive
                        step = 5  # Go to adaptive parameters

                elif step == 2:
                    # Question 3: Fixed radius
                    radius = prompt_float_with_retry(
                        processor=self.processor,
                        prompt="[green]Enter search radius[/green] (Å)",
                        default=answers.get('radius', 5.0),
                        min_value=0.1,
                        max_value=50.0,
                        module="Redox Detector - Refinement",
                        description="Enter fixed search radius"
                    )
                    answers['radius'] = radius
                    step = 10  # Go to atom filtering

                elif step == 5:
                    # Question 4: Adaptive min radius
                    min_radius = prompt_float_with_retry(
                        processor=self.processor,
                        prompt="[green]Enter minimum radius[/green] (Å)",
                        default=answers.get('min_radius', 3.0),
                        min_value=0.1,
                        max_value=20.0,
                        module="Redox Detector - Refinement",
                        description="Enter minimum radius for adaptive search"
                    )
                    answers['min_radius'] = min_radius
                    step = 6

                elif step == 6:
                    # Question 5: Adaptive max radius
                    max_radius = prompt_float_with_retry(
                        processor=self.processor,
                        prompt="[green]Enter maximum radius[/green] (Å)",
                        default=answers.get('max_radius', 8.0),
                        min_value=answers['min_radius'],
                        max_value=50.0,
                        module="Redox Detector - Refinement",
                        description="Enter maximum radius for adaptive search"
                    )
                    answers['max_radius'] = max_radius
                    step = 7

                elif step == 7:
                    # Question 6: Radius step
                    radius_step = prompt_float_with_retry(
                        processor=self.processor,
                        prompt="[green]Enter radius step size[/green] (Å)",
                        default=answers.get('radius_step', 0.5),
                        min_value=0.01,
                        max_value=5.0,
                        module="Redox Detector - Refinement",
                        description="Enter radius step size for adaptive search"
                    )
                    answers['radius_step'] = radius_step
                    step = 8

                elif step == 8:
                    # Question 7: Target coordination
                    target_coordination = prompt_int_with_retry(
                        processor=self.processor,
                        prompt="[green]Enter target coordination number[/green]",
                        default=answers.get('target_coordination', 6),
                        min_value=1,
                        max_value=20,
                        module="Redox Detector - Refinement",
                        description="Enter target coordination number"
                    )
                    answers['target_coordination'] = target_coordination
                    step = 10  # Go to atom filtering

                elif step == 10:
                    # Final step: Atom filtering (handled by separate method, not navigable here)
                    break  # Exit navigation loop

            except NavigationException:
                # User pressed back - determine where to go
                if step != 0:  # Can go back from any step except first
                    if step == 10:
                        # Coming back from atom filtering
                        # Return to appropriate step based on method choice
                        if answers['method_choice'] == '1':  # Fixed
                            step = 2  # Back to radius
                        else:  # Adaptive
                            step = 8  # Back to target coordination
                    elif step in [6, 7, 8]:
                        # In adaptive flow
                        step -= 1
                    elif step == 5:
                        # From first adaptive question, back to method choice
                        step = 1
                    elif step == 2:
                        # From fixed radius, back to method choice
                        step = 1
                    elif step == 1:
                        # From method choice, back to boundary or custom atom selection
                        if answers.get('boundary_choice') == '3':
                            step = 'custom'  # Back to custom atom selection
                        else:
                            step = 0  # Back to boundary selection
                    elif step == 'custom':
                        # From custom atom selection, back to boundary
                        step = 0

                    self.console.print(f"[grey50]← Going back...[/grey50]\n")

        # Build boundary definition from answers
        boundary_choice = answers['boundary_choice']
        if boundary_choice == '1':
            boundary_def = BoundaryDefinition.GEOMETRIC_CENTER
        elif boundary_choice == '2':
            boundary_def = BoundaryDefinition.MIN_DISTANCE_ALL
        else:  # '3'
            boundary_def = BoundaryDefinition.MIN_DISTANCE_CUSTOM

        # Build distance method from answers
        method_choice = answers['method_choice']
        if method_choice == '1':
            distance_method = DistanceMethod.FIXED
        else:
            distance_method = DistanceMethod.ADAPTIVE

        # Get custom boundary data from answers (collected during navigation loop)
        custom_coords = answers.get('custom_coords')
        custom_indices = answers.get('custom_indices')

        # Create template if in template mode (distance-based path)
        if self.template_mode and self.current_template is None:
            self.current_template = SiteTemplate(
                site_type="",  # Will be set later
                residue_types=[],  # Not used for distance search
                residue_counts=[],  # Not used for distance search
                atom_filtering_metals=[],  # Will be captured below
                atom_filtering_nonmetals=[],  # Will be captured below
                boundary_definition=boundary_def,
                custom_boundary_atom_indices=None,
                residue_selection_choice="",  # Will be captured in refinement
                bond_pairs=[],  # Will be captured in bond definition
                bond_selections={},  # Will be captured in bond definition
                continue_searching=False,  # Will be captured at end
                search_constraint="distance",
                distance_method="fixed" if distance_method == DistanceMethod.FIXED else "adaptive",
                distance_radius=answers.get('radius'),
                distance_min_radius=answers.get('min_radius'),
                distance_max_radius=answers.get('max_radius'),
                distance_radius_step=answers.get('radius_step'),
                distance_target_coordination=answers.get('target_coordination'),
            )

        # Store custom boundary indices in template if in template mode
        if boundary_choice == '3' and self.template_mode and self.current_template:
            self.current_template.custom_boundary_atom_indices = custom_indices

        # Get atom filtering (not part of navigation loop since it's a complex sub-workflow)
        search_metals, search_nonmetals = self._configure_search_atom_filtering(previous_search_metals, previous_search_nonmetals)

        # Build and return SearchParameters based on collected answers
        if distance_method == DistanceMethod.FIXED:
            return SearchParameters(
                constraint=SearchConstraint.DISTANCE_CUTOFF,
                distance_method=distance_method,
                boundary_definition=boundary_def,
                radius=answers['radius'],
                custom_boundary_coords=custom_coords,
                search_include_metals=search_metals,
                search_include_nonmetals=search_nonmetals
            )
        else:  # Adaptive
            return SearchParameters(
                constraint=SearchConstraint.DISTANCE_CUTOFF,
                distance_method=distance_method,
                boundary_definition=boundary_def,
                min_radius=answers['min_radius'],
                max_radius=answers['max_radius'],
                radius_step=answers['radius_step'],
                target_coordination=answers['target_coordination'],
                custom_boundary_coords=custom_coords,
                search_include_metals=search_metals,
                search_include_nonmetals=search_nonmetals
            )
    
    def _configure_count_search(self, site: RedoxSite) -> Optional[SearchParameters]:
        """Configure count-based search parameters with back navigation"""
        self.console.print(f"\n[bold underline]Count-Based Search Configuration[/bold underline]")
        self.console.print("Find N closest residues of specified types (regardless of distance)")

        # State tracking for navigation
        answers = {}
        step = 0

        while True:
            try:
                if step == 0:
                    # Question 1: Residue types
                    residue_types_input = prompt_with_navigation(
                        processor=self.processor,
                        prompt="[green]Enter residue types to find[/green] (comma-separated, e.g., HIS,CYS,ASP)",
                        default=answers.get('residue_types_input', "HIS,CYS"),
                        module="Redox Detector - Refinement",
                        description="Enter residue types for count search",
                        allow_back=False  # First question
                    ).strip()

                    if not residue_types_input:
                        self.console.print("[red]No residue types specified. Please try again.[/red]")
                        continue  # Ask again

                    answers['residue_types_input'] = residue_types_input
                    answers['residue_types'] = [t.strip().upper() for t in residue_types_input.split(',')]
                    step = 1

                elif step == 1:
                    # Question 2: Counts per residue type
                    residue_types = answers['residue_types']
                    self.console.print(f"\n[bold]Found residue types:[/bold] {', '.join(residue_types)}")

                    counts_input = prompt_with_navigation(
                        processor=self.processor,
                        prompt=f"[green]Enter number of each type to find[/green] (comma-separated, e.g., 2,3 for {len(residue_types)} types)",
                        default=answers.get('counts_input', "2"),
                        module="Redox Detector - Refinement",
                        description="Enter count of each residue type to find",
                        allow_back=True
                    ).strip()

                    try:
                        target_counts = [int(c.strip()) for c in counts_input.split(',')]
                        if len(target_counts) != len(residue_types):
                            self.console.print(f"[red]Number of counts ({len(target_counts)}) must match number of residue types ({len(residue_types)}). Please try again.[/red]")
                            continue  # Ask again
                    except ValueError:
                        self.console.print("[red]Invalid count values. Please enter comma-separated integers.[/red]")
                        continue  # Ask again

                    answers['counts_input'] = counts_input
                    answers['target_counts'] = target_counts
                    step = 2

                elif step == 2:
                    # Question 3: Boundary definition
                    self.console.print("\n[bold]Boundary definition (distance measurement origin):[/bold]")
                    self.console.print("[grey50]1.[/grey50] Geometric center (calculate centroid and measure distances from there)")
                    self.console.print("[grey50]2.[/grey50] Minimum distance from any site atoms")
                    self.console.print("[grey50]3.[/grey50] Minimum distance from custom selected atoms")

                    boundary_choice = prompt_with_navigation(
                        processor=self.processor,
                        prompt="[green]Select boundary definition[/green]",
                        choices=["1", "2", "3"],
                        default=answers.get('boundary_choice', "1"),
                        module="Redox Detector - Refinement",
                        description="Select boundary definition for count search",
                        options_map={
                            "1": "Geometric center",
                            "2": "Minimum distance from any site atoms",
                            "3": "Minimum distance from custom selected atoms"
                        },
                        allow_back=True
                    ).strip()

                    answers['boundary_choice'] = boundary_choice
                    # Modes 1/2 get their indicator immediately; mode 3
                    # waits until ``_select_custom_boundary_atoms`` runs
                    # below, so the halo lands on the actual picked
                    # subset rather than nothing.
                    if boundary_choice in ('1', '2'):
                        _show_boundary_indicator(
                            getattr(self, 'processor', None), site, boundary_choice,
                        )
                    break  # Done with navigation

            except NavigationException:
                # User pressed back
                if step > 0:
                    step -= 1
                    self.console.print(f"[grey50]← Going back...[/grey50]\n")

        # Extract answers
        residue_types = answers['residue_types']
        target_counts = answers['target_counts']
        boundary_choice = answers['boundary_choice']

        # Create mapping of residue type to count
        residue_count_map = dict(zip(residue_types, target_counts))
        total_count = sum(target_counts)

        # Capture template data if in template mode
        if self.template_mode and self.current_template is None:
            self.current_template = SiteTemplate(
                site_type="",  # Will be set later
                residue_types=residue_types,
                residue_counts=target_counts,
                atom_filtering_metals=[],  # Will be captured below
                atom_filtering_nonmetals=[],  # Will be captured below
                boundary_definition=BoundaryDefinition.GEOMETRIC_CENTER,  # Will be captured below
                custom_boundary_atom_indices=None,  # Will be captured if custom boundary selected
                residue_selection_choice="",  # Will be captured in refinement
                bond_pairs=[],  # Will be captured in bond definition
                bond_selections={},  # Will be captured in bond definition
                continue_searching=False  # Will be captured at end
            )

        # Handle custom boundary atom selection BEFORE search atom filtering
        # This ensures better UX - users select boundary atoms right after choosing custom boundary
        custom_coords = None
        custom_indices = None
        if boundary_choice == '3':
            result = self._select_custom_boundary_atoms(site)
            if not result:
                print("No custom atoms selected. Please reconfigure.")
                return None
            custom_coords, custom_indices = result

            _show_boundary_indicator(
                getattr(self, 'processor', None), site, '3',
                custom_indices=custom_indices,
            )

            # Store custom boundary indices in template if in template mode
            if self.template_mode and self.current_template:
                self.current_template.custom_boundary_atom_indices = custom_indices

        # Get previous search parameters if available
        previous_search_metals = None
        previous_search_nonmetals = None
        if site.search_history:
            last_search = site.search_history[-1]['search_parameters']
            previous_search_metals = last_search.search_include_metals
            previous_search_nonmetals = last_search.search_include_nonmetals

        # Ask about atom filtering (not part of navigation since it's a complex sub-workflow)
        search_metals, search_nonmetals = self._configure_search_atom_filtering(previous_search_metals, previous_search_nonmetals)

        # Capture atom filtering in template
        if self.template_mode and self.current_template:
            self.current_template.atom_filtering_metals = search_metals.copy()
            self.current_template.atom_filtering_nonmetals = search_nonmetals.copy()

        # Build boundary definition from answers
        if boundary_choice == '1':
            boundary_def = BoundaryDefinition.GEOMETRIC_CENTER
        elif boundary_choice == '2':
            boundary_def = BoundaryDefinition.MIN_DISTANCE_ALL
        else:  # '3'
            boundary_def = BoundaryDefinition.MIN_DISTANCE_CUSTOM

        # Capture boundary definition in template
        if self.template_mode and self.current_template:
            self.current_template.boundary_definition = boundary_def
        
        return SearchParameters(
            constraint=SearchConstraint.COUNT_CUTOFF,
            target_residue_types=residue_types,
            target_residue_count=total_count,
            target_residue_count_map=residue_count_map,
            boundary_definition=boundary_def,
            custom_boundary_coords=custom_coords,
            search_include_metals=search_metals,
            search_include_nonmetals=search_nonmetals
        )
    
    def _configure_search_atom_filtering(self, previous_search_metals=None, previous_search_nonmetals=None) -> Tuple[List[str], List[str]]:
        """Configure atom filtering for search phase with back navigation"""
        print("\n=== Search Atom Filtering Configuration ===")

        # State tracking for navigation
        answers = {}
        step = 0

        while True:
            try:
                if step == 0:
                    # Question 1: Use previous/current filtering?
                    if previous_search_metals is not None and previous_search_nonmetals is not None:
                        self.console.print("\n[bold]Previous search filtering:[/bold]")
                        self.console.print(f"  Include metals: [cyan]{previous_search_metals}[/cyan]")
                        self.console.print(f"  Include non-metals: [cyan]{previous_search_nonmetals}[/cyan]")
                        choice = confirm_with_context(
                            processor=self.processor,
                            prompt="[green]Use same atom filtering for search?[/green]",
                            default=answers.get('use_previous', True),
                            module="Redox Detector - Refinement",
                            description="Use previous search atom filtering",
                            allow_back=False  # First question
                        )
                        answers['use_previous'] = choice
                        if choice:
                            # User wants to use previous - return immediately
                            return previous_search_metals, previous_search_nonmetals
                        else:
                            step = 1  # Continue to custom configuration
                    else:
                        self.console.print("\n[bold]Current filtering:[/bold]")
                        self.console.print(f"  Include metals: [cyan]{self.config.inventory_include_metals}[/cyan]")
                        self.console.print(f"  Include non-metals: [cyan]{self.config.inventory_include_nonmetals}[/cyan]")
                        choice = confirm_with_context(
                            processor=self.processor,
                            prompt="[green]Use same atom filtering for search?[/green]",
                            default=answers.get('use_current', True),
                            module="Redox Detector - Refinement",
                            description="Use current atom filtering for search",
                            allow_back=False  # First question
                        )
                        answers['use_current'] = choice
                        if choice:
                            # User wants to use current - return immediately
                            return self.config.inventory_include_metals, self.config.inventory_include_nonmetals
                        else:
                            step = 1  # Continue to custom configuration

                elif step == 1:
                    # Question 2: Configure metals
                    self.console.print("\n[bold underline]Configure search-specific atom filtering:[/bold underline]")

                    metals_input = prompt_with_navigation(
                        processor=self.processor,
                        prompt="[green]Include metals[/green] (comma-separated list, 'none', or 'all')",
                        default=answers.get('metals_input', "all"),
                        module="Redox Detector - Site Refinement",
                        description="Configure search metals filter",
                        options_map={"all": "All metals", "none": "No metals", "custom": "Custom list (comma-separated)"},
                        allow_back=True
                    ).strip().upper()
                    answers['metals_input'] = metals_input
                    step = 2

                elif step == 2:
                    # Question 3: Configure non-metals
                    nonmetals_input = prompt_with_navigation(
                        processor=self.processor,
                        prompt="[green]Include non-metals[/green] (comma-separated list, 'none', or 'all')",
                        default=answers.get('nonmetals_input', "all"),
                        module="Redox Detector - Site Refinement",
                        description="Configure search nonmetals filter",
                        options_map={"all": "All non-metals", "none": "No non-metals", "custom": "Custom list (comma-separated)"},
                        allow_back=True
                    ).strip().upper()
                    answers['nonmetals_input'] = nonmetals_input
                    break  # Done

            except NavigationException:
                # User pressed back
                if step > 0:
                    if step == 2:
                        step = 1  # From nonmetals back to metals
                    elif step == 1:
                        step = 0  # From metals back to use previous/current question
                    self.console.print(f"[grey50]← Going back...[/grey50]\n")

        # Build results from answers
        metals_input = answers['metals_input']
        if metals_input == 'NONE':
            search_include_metals = ['NONE']
        elif metals_input == 'ALL' or metals_input == '':
            search_include_metals = []
        else:
            search_include_metals = [m.strip() for m in metals_input.split(',') if m.strip()]

        nonmetals_input = answers['nonmetals_input']
        if nonmetals_input == 'NONE':
            search_include_nonmetals = ['NONE']
        elif nonmetals_input == 'ALL' or nonmetals_input == '':
            search_include_nonmetals = []
        else:
            search_include_nonmetals = [nm.strip() for nm in nonmetals_input.split(',') if nm.strip()]

        return search_include_metals, search_include_nonmetals
    
    def _select_custom_boundary_atoms(self, site: RedoxSite) -> Optional[Tuple[List[Tuple[float, float, float]], List[int]]]:
        """
        Let user select custom atoms for boundary definition.

        Returns:
            Tuple of (coordinates, indices) where indices are 0-based positions in site.atoms
            or None if cancelled
        """
        print("\n=== Custom Boundary Atom Selection ===")

        if not site.atoms:
            self.console.print("[red]No atoms in current site to select from.[/red]")
            return None
        
        # Display all atoms in the site using Rich table
        atom_table = Table(title="Available atoms in current site")
        atom_table.add_column("#", style="cyan", width=6)
        atom_table.add_column("Residue", style="green", width=12)
        atom_table.add_column("Atom", style="yellow", width=8)
        atom_table.add_column("Element", style="magenta", width=8)

        atom_list = []
        for i, atom in enumerate(site.atoms, 1):
            atom_list.append(atom)
            residue_info = f"{atom.resname} {atom.chain}:{atom.resid}"
            atom_table.add_row(f"[{i}]", residue_info, atom.atom_name, atom.element)

        self.console.print(atom_table)

        while True:
            selection_input = prompt_with_context(
                processor=self.processor,
                prompt=f"[green]Select atoms for boundary[/green] (1-{len(atom_list)}, atom names, comma-separated for multiple, ranges, or 'all'/'cancel')",
                module="Redox Detector - Site Refinement",
                description="Select custom boundary atoms",
                options_map={"all": "All atoms", "cancel": "Cancel selection", "custom": "Numbers, names, or ranges (e.g., 1-4, FE,NA,NC, 1,2,5-8)"}
            ).strip()

            if selection_input.lower() == 'cancel':
                return None
            elif selection_input.lower() == 'all':
                all_indices = list(range(len(atom_list)))
                return ([atom.coords for atom in atom_list], all_indices)
            else:
                # Try to parse as numbers/ranges first
                indices = self._parse_selection_input(selection_input, len(atom_list))

                # If that fails, try to parse as atom names
                if indices is None:
                    indices = self._parse_atom_name_selection(selection_input, atom_list)

                if indices is not None:
                    selected_coords = [atom_list[i].coords for i in indices]

                    # Show selection confirmation with Rich table
                    selection_table = Table(title=f"Selected {len(selected_coords)} atoms")
                    selection_table.add_column("Residue", style="green")
                    selection_table.add_column("Atom", style="yellow")
                    selection_table.add_column("Element", style="magenta")

                    for i in indices:
                        atom = atom_list[i]
                        residue_info = f"{atom.resname} {atom.chain}:{atom.resid}"
                        selection_table.add_row(residue_info, atom.atom_name, atom.element)

                    self.console.print(selection_table)

                    return (selected_coords, indices)
                else:
                    self.console.print(f"[red]Invalid selection. Please use numbers (1-{len(atom_list)}), ranges (1-4), atom names (FE,NA,NC), or comma-separated combinations.[/red]")
                    continue

    def _parse_atom_name_selection(self, input_str: str, atom_list: List) -> Optional[List[int]]:
        """Parse selection input as atom names (comma-separated)

        Examples:
        - "FE" -> indices of all atoms named FE
        - "FE,NA,NC" -> indices of atoms named FE, NA, or NC
        - If multiple atoms share a name, prompts for clarification

        Args:
            input_str: User input string with atom names
            atom_list: List of RedoxSiteAtom objects

        Returns:
            List of 0-based indices, or None if invalid/cancelled
        """
        try:
            atom_names = [name.strip().upper() for name in input_str.split(',') if name.strip()]
            if not atom_names:
                return None

            selected_indices = []

            for atom_name in atom_names:
                # Find all atoms with this name
                matching_atoms = [(i, atom) for i, atom in enumerate(atom_list)
                                if atom.atom_name.upper() == atom_name]

                if len(matching_atoms) == 0:
                    self.console.print(f"[red]No atoms named '{atom_name}' found.[/red]")
                    return None
                elif len(matching_atoms) == 1:
                    selected_indices.append(matching_atoms[0][0])
                else:
                    # Multiple atoms with same name - ask for clarification
                    self.console.print(f"[yellow]Multiple atoms named '{atom_name}' found:[/yellow]")
                    for idx, (atom_idx, atom) in enumerate(matching_atoms, 1):
                        original_number = atom_idx + 1
                        residue_info = f"{atom.resname} {atom.chain}:{atom.resid}"
                        self.console.print(f"  {idx}. [{original_number}] {atom_name} in {residue_info}")

                    try:
                        sub_choice = prompt_with_context(
                            processor=self.processor,
                            prompt=f"[green]Select which '{atom_name}' atom(s)[/green] (1-{len(matching_atoms)}, comma-separated for multiple)",
                            module="Redox Detector - Site Refinement",
                            description=f"Select which '{atom_name}' atom",
                            options_map={str(i+1): f"Atom {i+1}" for i in range(len(matching_atoms))}
                        ).strip()

                        # Parse sub-selection (can be comma-separated)
                        sub_choices = [c.strip() for c in sub_choice.split(',') if c.strip()]
                        for choice in sub_choices:
                            sub_idx = int(choice) - 1
                            if 0 <= sub_idx < len(matching_atoms):
                                selected_indices.append(matching_atoms[sub_idx][0])
                            else:
                                self.console.print(f"[red]Invalid choice: {choice}[/red]")
                                return None
                    except (ValueError, IndexError):
                        self.console.print(f"[red]Invalid selection for '{atom_name}'[/red]")
                        return None

            return sorted(list(set(selected_indices)))  # Remove duplicates and sort

        except Exception as e:
            logger.debug(f"Error parsing atom name selection: {e}")
            return None

    def _perform_distance_search(self, site: RedoxSite, structure: Structure,
                                params: SearchParameters) -> SearchResult:
        """Perform distance-based search with specified parameters"""
        
        # Get boundary coordinates
        boundary_coords = site.get_site_boundary_coords(
            params.boundary_definition, 
            params.custom_boundary_coords
        )
        
        if not boundary_coords:
            self.console.print("[red]No boundary coordinates available for search.[/red]")
            return SearchResult([], [], params, [])
        
        if params.distance_method == DistanceMethod.FIXED:
            return self._fixed_distance_search(site, structure, boundary_coords, params)
        else:
            return self._adaptive_distance_search(site, structure, boundary_coords, params)
    
    def _perform_count_search(self, site: RedoxSite, structure: Structure,
                             params: SearchParameters) -> SearchResult:
        """Perform count-based search for closest residues"""

        # Calculate boundary coordinates based on boundary definition
        boundary_coords = site.get_site_boundary_coords(
            params.boundary_definition,
            params.custom_boundary_coords
        )

        # Use boundary coords if available, otherwise fall back to all site atoms
        if boundary_coords:
            search_coords = boundary_coords
        else:
            search_coords = [center.coords for center in site.centers] + [atom.coords for atom in site.atoms]

        if not search_coords:
            return SearchResult([], [], params, [])

        # Find all residues of target types
        target_residues = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.resname in params.target_residue_types:
                        # Skip if residue is already in site
                        res_key = (chain.id, residue.id[1], residue.id[2])
                        if res_key in site.residue_groups:
                            continue

                        # Calculate minimum distance from boundary coords to any residue atom
                        min_distance = float('inf')
                        for residue_atom in residue.get_atoms():
                            # Apply search atom filtering
                            if not self._atom_meets_search_criteria(residue_atom, params):
                                continue

                            for boundary_coord in search_coords:
                                distance = np.linalg.norm(np.array(boundary_coord) - np.array(residue_atom.coord))
                                min_distance = min(min_distance, distance)

                        # Skip residues with no atoms meeting criteria
                        if min_distance == float('inf'):
                            continue

                        target_residues.append({
                            'chain': chain.id,
                            'resname': residue.resname,
                            'resid': residue.id[1],
                            'insertion_code': residue.id[2],
                            'distance': min_distance,
                            'residue_obj': residue
                        })
        
        # Group residues by type and select specified number of each
        if params.target_residue_count_map:
            # Use per-type counts
            closest_residues = []
            residues_by_type = {}
            
            # Group residues by type
            for res in target_residues:
                res_type = res['resname']
                if res_type not in residues_by_type:
                    residues_by_type[res_type] = []
                residues_by_type[res_type].append(res)
            
            # Sort each type by distance and take the specified count
            for res_type, count in params.target_residue_count_map.items():
                if res_type in residues_by_type:
                    sorted_type_residues = sorted(residues_by_type[res_type], key=lambda x: x['distance'])
                    closest_residues.extend(sorted_type_residues[:count])
        else:
            # Fall back to original behavior (total count)
            target_residues.sort(key=lambda x: x['distance'])
            closest_residues = target_residues[:params.target_residue_count]
        
        # Convert to detected atoms format
        detected_atoms = []
        detected_residues = []
        
        for res_info in closest_residues:
            residue = res_info['residue_obj']
            
            # Add residue info
            detected_residues.append({
                'chain': res_info['chain'],
                'resname': res_info['resname'],
                'resid': res_info['resid'],
                'insertion_code': res_info['insertion_code'],
                'distance': res_info['distance']
            })
            
            # Add all atoms in residue that meet criteria
            for atom in residue.get_atoms():
                # Apply search-specific atom selection criteria
                if not self._atom_meets_search_criteria(atom, params):
                    continue

                # Calculate minimum distance from this atom to boundary coords
                atom_distance = min(
                    np.linalg.norm(np.array(boundary_coord) - np.array(atom.coord))
                    for boundary_coord in search_coords
                )
                
                # Classify potential bonds
                bond_type, chemical_type = self._classify_atom_bond_potential(
                    atom, site, atom_distance
                )
                
                detected_atoms.append({
                    'chain': res_info['chain'],
                    'resname': res_info['resname'],
                    'resid': res_info['resid'],
                    'atom_name': atom.name,
                    'element': atom.element,
                    'coords': tuple(round(x, 3) for x in atom.coord),
                    'distance': atom_distance,
                    'bond_type': bond_type,
                    'chemical_type': chemical_type,
                    'insertion_code': res_info['insertion_code']
                })

        return SearchResult(
            detected_atoms=detected_atoms,
            detected_residues=detected_residues,
            search_parameters=params,
            boundary_coords=boundary_coords,
            total_atoms_found=len(detected_atoms),
            total_residues_found=len(detected_residues)
        )
    
    def _fixed_distance_search(self, site: RedoxSite, structure: Structure, boundary_coords: List[Tuple[float, float, float]],
                              params: SearchParameters) -> SearchResult:
        """Perform fixed radius distance search"""
        
        # Build atom list for efficient searching
        atoms = list(structure.get_atoms())
        ns = NeighborSearch(atoms)
        
        detected_atoms = []
        detected_residues = set()
        
        # Search from each boundary point
        for boundary_coord in boundary_coords:
            nearby_atoms = ns.search(boundary_coord, params.radius)
            
            for atom in nearby_atoms:
                # Apply search-specific atom selection criteria
                if not self._atom_meets_search_criteria(atom, params):
                    continue
                
                distance = np.linalg.norm(np.array(boundary_coord) - np.array(atom.coord))
                
                # Find metal center in site for proper bond classification
                metal_element = None
                for redox_center in site.centers:
                    if redox_center.element and redox_center.element.upper() in METALS:
                        metal_element = redox_center.element
                        break
                
                # Classify bond properly using metal content
                if metal_element:
                    bond_type, chemical_type = classify_bond_types(
                        atom1_element=metal_element,
                        atom2_element=atom.element,
                        distance=distance,
                        atom1_residue="",
                        atom2_residue="",
                        atom1_resid=0,
                        atom2_resid=0,
                        atom1_chain="",
                        atom2_chain=""
                    )
                else:
                    # No metal center found, treat as covalent
                    bond_type = "interresidue"
                    chemical_type = "covalent"
                
                atom_info = {
                    'chain': atom.parent.parent.id,
                    'resname': atom.parent.resname,
                    'resid': atom.parent.id[1],
                    'atom_name': atom.name,
                    'element': atom.element,
                    'coords': tuple(round(x, 3) for x in atom.coord),
                    'distance': distance,
                    'bond_type': bond_type,
                    'chemical_type': chemical_type,
                    'insertion_code': atom.parent.id[2]
                }
                detected_atoms.append(atom_info)
                
                # Track unique residues
                res_key = (atom.parent.parent.id, atom.parent.resname, atom.parent.id[1], atom.parent.id[2])
                detected_residues.add(res_key)
        
        # Convert residue set to list format
        residue_list = []
        for chain, resname, resid, icode in detected_residues:
            residue_list.append({
                'chain': chain,
                'resname': resname,
                'resid': resid,
                'insertion_code': icode
            })
        
        return SearchResult(
            detected_atoms=detected_atoms,
            detected_residues=residue_list,
            search_parameters=params,
            boundary_coords=boundary_coords,
            search_radius_used=params.radius,
            total_atoms_found=len(detected_atoms),
            total_residues_found=len(residue_list)
        )
    
    def _adaptive_distance_search(self, site: RedoxSite, structure: Structure, boundary_coords: List[Tuple[float, float, float]],
                                 params: SearchParameters) -> SearchResult:
        """Perform adaptive radius search to achieve target coordination"""
        
        atoms = list(structure.get_atoms())
        ns = NeighborSearch(atoms)
        
        # Try different radii to achieve target coordination
        current_radius = params.min_radius
        best_result = None
        
        while current_radius <= params.max_radius:
            detected_atoms = []
            detected_residues = set()
            
            # Search with current radius
            for boundary_coord in boundary_coords:
                nearby_atoms = ns.search(boundary_coord, current_radius)
                
                for atom in nearby_atoms:
                    # Apply search-specific atom selection criteria
                    if not self._atom_meets_search_criteria(atom, params):
                        continue
                    
                    distance = np.linalg.norm(np.array(boundary_coord) - np.array(atom.coord))
                    
                    # Only count meaningful coordination atoms (non-hydrogen, reasonable distance)
                    if atom.element != 'H' and distance <= current_radius:
                        # Find metal center in site for proper bond classification
                        metal_element = None
                        for redox_center in site.centers:
                            if redox_center.element and redox_center.element.upper() in METALS:
                                metal_element = redox_center.element
                                break
                        
                        # Classify bond properly using metal content
                        if metal_element:
                            bond_type, chemical_type = classify_bond_types(
                                atom1_element=metal_element,
                                atom2_element=atom.element,
                                distance=distance,
                                atom1_residue="",
                                atom2_residue="",
                                atom1_resid=0,
                                atom2_resid=0,
                                atom1_chain="",
                                atom2_chain=""
                            )
                        else:
                            # No metal center found, treat as covalent
                            bond_type = "interresidue"
                            chemical_type = "covalent"
                        
                        atom_info = {
                            'chain': atom.parent.parent.id,
                            'resname': atom.parent.resname,
                            'resid': atom.parent.id[1],
                            'atom_name': atom.name,
                            'element': atom.element,
                            'coords': tuple(round(x, 3) for x in atom.coord),
                            'distance': distance,
                            'bond_type': bond_type,
                            'chemical_type': chemical_type,
                            'insertion_code': atom.parent.id[2]
                        }
                        detected_atoms.append(atom_info)
                        
                        res_key = (atom.parent.parent.id, atom.parent.resname, atom.parent.id[1], atom.parent.id[2])
                        detected_residues.add(res_key)
            
            # Check if we've achieved target coordination
            coordination_count = len([a for a in detected_atoms if a['chemical_type'] in ['coordinate', 'covalent']])
            
            if coordination_count >= params.target_coordination:
                # Found target coordination
                residue_list = []
                for chain, resname, resid, icode in detected_residues:
                    residue_list.append({
                        'chain': chain,
                        'resname': resname,
                        'resid': resid,
                        'insertion_code': icode
                    })
                
                return SearchResult(
                    detected_atoms=detected_atoms,
                    detected_residues=residue_list,
                    search_parameters=params,
                    boundary_coords=boundary_coords,
                    search_radius_used=current_radius,
                    total_atoms_found=len(detected_atoms),
                    total_residues_found=len(residue_list)
                )
            
            current_radius += params.radius_step
        
        # If we didn't achieve target coordination, return the best result from max radius
        print(f"Warning: Could not achieve target coordination of {params.target_coordination}")
        print(f"Returning results from maximum radius {params.max_radius}Å")
        
        return self._fixed_distance_search(site, structure, boundary_coords, 
                                         SearchParameters(
                                             constraint=SearchConstraint.DISTANCE_CUTOFF,
                                             distance_method=DistanceMethod.FIXED,
                                             boundary_definition=params.boundary_definition,
                                             radius=params.max_radius,
                                             search_include_metals=params.search_include_metals,
                                             search_include_nonmetals=params.search_include_nonmetals
                                         ))
    
    def _classify_atom_bond_potential(self, atom, site: RedoxSite, distance: float) -> Tuple[str, str]:
        """Classify the potential bond type for a detected atom"""
        # Find a metal center in the site to determine bond classification
        metal_element = None
        for redox_center in site.centers:
            if redox_center.element and redox_center.element.upper() in METALS:
                metal_element = redox_center.element
                break
        
        if metal_element:
            # Use the proper classification function with metal center
            bond_type, chemical_type = classify_bond_types(
                atom1_element=metal_element,
                atom2_element=atom.element,
                distance=distance,
                atom1_residue="",  # Not needed for chemical_type classification
                atom2_residue="",
                atom1_resid=0,
                atom2_resid=0,
                atom1_chain="",
                atom2_chain=""
            )
        else:
            # Fallback for sites without metal centers
            bond_type = "interresidue"
            chemical_type = "covalent"  # No metals means covalent
        
        return bond_type, chemical_type

    def _get_site_residue_keys(self, site: RedoxSite) -> Set[Tuple[str, str, int, str]]:
        """Get set of (chain, resname, resid, icode) tuples for all residues in the site.

        This includes residues from both atoms and centers to ensure complete coverage.
        """
        existing_residues = set()

        # Add residues from all atoms in site
        for atom in site.atoms:
            res_key = (atom.chain, atom.resname, atom.resid, atom.insertion_code)
            existing_residues.add(res_key)

        # Add residues from centers (they may not have corresponding atoms yet)
        for center in site.centers:
            res_key = (center.chain, center.resname, center.resid, center.insertion_code)
            existing_residues.add(res_key)

        return existing_residues

    def _filter_search_result_for_existing_residues(self, site: RedoxSite, search_result: SearchResult) -> SearchResult:
        """Filter search results to exclude residues already in the site.

        Returns a new SearchResult with atoms/residues from existing site residues removed.
        """
        existing_residues = self._get_site_residue_keys(site)

        # Filter detected atoms
        filtered_atoms = []
        for atom_info in search_result.detected_atoms:
            res_key = (atom_info['chain'], atom_info['resname'], atom_info['resid'], atom_info['insertion_code'])
            if res_key not in existing_residues:
                filtered_atoms.append(atom_info)

        # Filter detected residues
        filtered_residues = []
        for res_info in search_result.detected_residues:
            res_key = (res_info['chain'], res_info['resname'], res_info['resid'], res_info['insertion_code'])
            if res_key not in existing_residues:
                filtered_residues.append(res_info)

        # Return new SearchResult with filtered data
        return SearchResult(
            detected_atoms=filtered_atoms,
            detected_residues=filtered_residues,
            search_parameters=search_result.search_parameters,
            boundary_coords=search_result.boundary_coords,
            search_radius_used=search_result.search_radius_used,
            total_atoms_found=len(filtered_atoms),
            total_residues_found=len(filtered_residues)
        )

    def _process_search_results(self, site: RedoxSite, search_result: SearchResult,
                               structure: Structure) -> Tuple[RedoxSite, str]:
        """Process search results with user action choice"""

        # Filter out residues already in the site to prevent duplication
        filtered_result = self._filter_search_result_for_existing_residues(site, search_result)

        if not filtered_result.detected_atoms:
            if search_result.detected_atoms:
                # Atoms were found but all belong to residues already in the site
                self.console.print("[yellow]All detected residues are already in the site. "
                                 "Try a larger search radius to find additional residues.[/yellow]")
            else:
                self.console.print("[yellow]No atoms detected in search.[/yellow]")
            return site, "continue"
        
        # Display search results (filtered to exclude existing residues)
        self._display_search_results(filtered_result)

        # Direct residue selection with option to adjust search parameters
        selected_residues = self._select_residues_to_add(site, filtered_result, structure)

        if selected_residues == "finish":
            # User chose to finish refinement without adding any found
            # residues. Signal the loop to finalize (bonds) and exit.
            return site, "finish_refine"

        if selected_residues is None:
            # User chose to adjust search parameters (entered "0")
            site.search_history.append({
                'search_parameters': filtered_result.search_parameters,
                'atoms_found': filtered_result.total_atoms_found,
                'residues_found': filtered_result.total_residues_found,
                'residues_added': 0,  # No residues added since we're adjusting
                'atoms_added': 0,     # No atoms added since we're adjusting
                'search_radius': filtered_result.search_radius_used,
                'adjustment_only': True  # Flag to indicate this was just for parameter adjustment
            })
            return site, "adjust_parameters"
        elif not selected_residues:
            self.console.print("[yellow]No residues selected for addition.[/yellow]")
            return site, "continue"

        # Add complete residues to site
        updated_site = copy.deepcopy(site)
        atoms_added = self._add_complete_residues_to_site(updated_site, selected_residues, structure)

        if atoms_added == 0:
            self.console.print("[yellow]No atoms were added.[/yellow]")
            return site, "continue"

        self.console.print(f"[green]Added {atoms_added} atoms from {len(selected_residues)} complete residues to site.[/green]")

        # User-defined bond specification
        self._define_bonds_interactively(updated_site)

        # Update search history
        updated_site.search_history.append({
            'search_parameters': filtered_result.search_parameters,
            'atoms_found': filtered_result.total_atoms_found,
            'residues_found': filtered_result.total_residues_found,
            'residues_added': len(selected_residues),
            'atoms_added': atoms_added,
            'search_radius': filtered_result.search_radius_used
        })
        
        self.console.print(f"Site now contains [cyan]{len(updated_site.atoms)}[/cyan] atoms, [cyan]{len(updated_site.bonds)}[/cyan] bonds")

        # Ask if user wants to continue searching with expanded site
        continue_choice = confirm_with_context(
            processor=self.processor,
            prompt="[green]Search for additional residues to add to this site?[/green]",
            default=True,
            module="Redox Detector - Site Refinement",
            description="Search for additional residues to add to this site"
        )
        if continue_choice:
            return updated_site, "continue"
        else:
            return updated_site, "finish"
    
    def _display_search_results(self, search_result: SearchResult):
        """Display search results showing complete residues that will be added"""
        self.console.print(f"\n[bold underline]SEARCH RESULTS[/bold underline]")
        self.console.print(f"Search method: [cyan]{search_result.search_parameters.get_method_description()}[/cyan]")
        if search_result.search_radius_used:
            self.console.print(f"Search radius used: [cyan]{search_result.search_radius_used:.2f}Å[/cyan]")
        self.console.print(f"Total atoms found: [cyan]{search_result.total_atoms_found}[/cyan]")
        self.console.print(f"Total residues found: [cyan]{search_result.total_residues_found}[/cyan]")

        if not search_result.detected_atoms:
            self.console.print("[yellow]No atoms detected.[/yellow]")
            return

        # Group atoms by residue for display
        residue_atoms = {}
        for atom_info in search_result.detected_atoms:
            res_key = (atom_info['chain'], atom_info['resname'], atom_info['resid'], atom_info['insertion_code'])
            if res_key not in residue_atoms:
                residue_atoms[res_key] = []
            residue_atoms[res_key].append(atom_info)

        # Auto-highlight each candidate residue in the viewer with the
        # same 1-based palette ID the table uses below. Drawn before the
        # table prints so the user sees both views land together.
        _auto_show_palette_highlights(
            getattr(self, 'processor', None),
            [
                {
                    "selection": f":{chain} and {resid}",
                    "label": f"candidate_{i}",
                }
                for i, (chain, resname, resid, icode) in enumerate(residue_atoms.keys(), start=1)
            ],
        )
        
        # Create Rich table for results
        results_table = Table(title="Detected residues that have atoms meeting search criteria")
        results_table.add_column("#", style="cyan", width=6)
        results_table.add_column("Residue", style="green")
        results_table.add_column("Atoms", style="yellow")
        results_table.add_column("Details", style="grey50")
        
        for i, (res_key, atoms) in enumerate(residue_atoms.items(), 1):
            chain, resname, resid, icode = res_key
            residue_str = f"{resname} {chain}:{resid}"
            
            # Use appropriate text based on search type
            if search_result.search_parameters.constraint == SearchConstraint.COUNT_CUTOFF:
                atoms_str = f"{len(atoms)} atoms meeting criteria"
            else:
                atoms_str = f"{len(atoms)} atoms within range"
            
            # Create details string with atom info
            details_parts = []
            for atom in atoms:
                details_parts.append(f"{atom['atom_name']} {atom['element']} ({atom['distance']:.2f}Å)")
            details_str = ", ".join(details_parts[:3])  # Show first 3 atoms
            if len(atoms) > 3:
                details_str += f" + {len(atoms) - 3} more"
            
            results_table.add_row(f"[{i}]", residue_str, atoms_str, details_str)
        
        self.console.print(results_table)
    
    def _select_residues_to_add(self, site: RedoxSite, search_result: SearchResult, structure: Structure) -> Union[List[Tuple[str, int, str]], str, None]:
        """Interactive selection of which detected residues to add completely to the site

        Returns a list of selected residues, or one of two sentinels:
        ``None`` to adjust search parameters and retry (0 / none), or the string
        ``"finish"`` to finish refinement without adding anything (finish/done/exit).

        Args:
            site: The current RedoxSite (for viewing original centers as context)
            search_result: Search results containing detected residues
            structure: BioPython structure object

        Returns:
            List of (chain, resid, icode) tuples for selected residues,
            empty list if 'none' selected, or None if user wants to adjust search parameters
        """
        
        if not search_result.detected_atoms:
            return []
        
        # Group by residue for selection
        residue_atoms = {}
        for atom_info in search_result.detected_atoms:
            res_key = (atom_info['chain'], atom_info['resname'], atom_info['resid'], atom_info['insertion_code'])
            if res_key not in residue_atoms:
                residue_atoms[res_key] = []
            residue_atoms[res_key].append(atom_info)
        
        self.console.print(f"\n[bold underline]RESIDUE SELECTION[/bold underline]")
        self.console.print("Select residues to add (ALL atoms from selected residues will be added)")
        
        residue_list = list(residue_atoms.keys())
        
        # Get total atom counts for each residue from structure
        for i, res_key in enumerate(residue_list, 1):
            chain, resname, resid, icode = res_key
            
            # Find the actual residue in structure to get total atom count
            total_atoms = 0
            for model in structure:
                for struct_chain in model:
                    if struct_chain.id == chain:
                        for residue in struct_chain:
                            if (residue.id[1] == resid and 
                                residue.resname == resname and 
                                residue.id[2] == icode):
                                total_atoms = len(list(residue.get_atoms()))
                                break
            
            criteria_count = len(residue_atoms[res_key])
            # Building table row data for later display
        
        # Create table for residue selection
        selection_table = Table(title="Residues Available for Addition")
        selection_table.add_column("#", style="cyan", width=6)
        selection_table.add_column("Residue", style="green")
        selection_table.add_column("Will Add", style="yellow")
        selection_table.add_column("Met Criteria", style="magenta")
        
        for i, res_key in enumerate(residue_list, 1):
            chain, resname, resid, icode = res_key
            
            # Find the actual residue in structure to get total atom count
            total_atoms = 0
            for model in structure:
                for struct_chain in model:
                    if struct_chain.id == chain:
                        for residue in struct_chain:
                            if (residue.id[1] == resid and 
                                residue.resname == resname and 
                                residue.id[2] == icode):
                                total_atoms = len(list(residue.get_atoms()))
                                break
            
            criteria_count = len(residue_atoms[res_key])
            selection_table.add_row(
                f"[{i}]",
                f"{resname} {chain}:{resid}",
                f"all {total_atoms} atoms",
                f"{criteria_count} atoms"
            )
        
        self.console.print(selection_table)

        self.console.print("\n[grey50]Numbers, ranges, or 'all'/'none'/'view'/'finish'[/grey50]")
        self.console.print("[grey50]Examples: 1-3, 1,3, 1-2 4, view, 0, or finish[/grey50]")
        self.console.print("[grey50]Note: 0 or none = Adjust search parameters and repeat analysis[/grey50]")
        self.console.print("[grey50]Note: finish = Finish refinement now (add nothing more, keep the site as-is)[/grey50]")
        self.console.print("[grey50]Note: view = View available residues in 3D structure viewer[/grey50]")

        while True:
            selection_input = prompt_with_context(
                processor=self.processor,
                prompt="[green]Select residues to add[/green]",
                module="Redox Detector - Site Refinement",
                description="Select residues to add to site",
                options_map={"all": "All residues", "none": "None (adjust search)", "0": "Adjust search parameters", "finish": "Finish refinement (add nothing more)", "view": "View in 3D viewer", "custom": "Custom selection (e.g., 1-3, 1,3)"}
            ).strip()

            if selection_input == '0' or selection_input.lower() == 'none':
                return None  # Signal to adjust search parameters
            elif selection_input.lower() in ('finish', 'done', 'exit'):
                return "finish"  # Signal to finish refinement and exit the loop
            elif selection_input.lower() == 'all':
                return [(chain, resid, icode) for chain, resname, resid, icode in residue_list]
            elif selection_input.lower() == 'view':
                # Launch viewer with available residues AND original centers (for context)
                # Create a temporary site-like object with centers + found residues
                self._launch_viewer_for_residues_with_centers(site, residue_list, "Available Residues + Original Centers")
                continue  # Return to prompt after viewing
            else:
                indices = self._parse_selection_input(selection_input, len(residue_list))
                if indices is not None:
                    selected_residues = []
                    for i in indices:
                        chain, resname, resid, icode = residue_list[i]
                        selected_residues.append((chain, resid, icode))
                    return selected_residues
                else:
                    self.console.print(f"[red]Invalid selection. Please use numbers 1-{len(residue_list)}, ranges (1-3), comma-separated (1,2,4), 'view', 0 to adjust, or 'finish' to exit.[/red]")
                    continue
    
    def _atom_meets_inventory_criteria(self, atom) -> bool:
        """Check if atom meets the inventory filtering criteria"""
        element = atom.element.upper()
        atom_name = atom.name.upper()
        
        # Check if atom is metal
        if element in METALS:
            # Metal atom - check against inventory_include_metals criteria
            if not self.config.inventory_include_metals:
                # Empty list means include all metals
                return True
            elif self.config.inventory_include_metals == ['NONE']:
                # Special marker for no metals
                return False
            else:
                # Check if this metal is in the include list
                return element in [m.upper() for m in self.config.inventory_include_metals]
        else:
            # Non-metal atom - check against inventory_include_nonmetals criteria
            if not self.config.inventory_include_nonmetals:
                # Empty list means include all non-metals
                return True
            elif self.config.inventory_include_nonmetals == ['NONE']:
                # Special marker for no non-metals
                return False
            else:
                # Check if element name OR atom name is in the include list
                nonmetal_criteria = [nm.upper() for nm in self.config.inventory_include_nonmetals]
                return element in nonmetal_criteria or atom_name in nonmetal_criteria
    
    def _parse_selection_input(self, input_str: str, max_num: int) -> Optional[List[int]]:
        """Parse selection input supporting ranges and comma-separated values
        
        Examples:
        - "1-10" -> [0,1,2,3,4,5,6,7,8,9]
        - "1,3,5-8" -> [0,2,4,5,6,7]
        - "1-5 10-15" -> [0,1,2,3,4,9,10,11,12,13,14]
        """
        try:
            indices = set()
            
            # Split by both spaces and commas
            parts = input_str.replace(',', ' ').split()
            
            for part in parts:
                if '-' in part:
                    # Handle range (e.g., "1-10")
                    range_parts = part.split('-')
                    if len(range_parts) == 2:
                        start = int(range_parts[0])
                        end = int(range_parts[1])
                        if start <= end:
                            for i in range(start, end + 1):
                                if 1 <= i <= max_num:
                                    indices.add(i - 1)  # Convert to 0-based
                                else:
                                    return None  # Invalid range
                        else:
                            return None  # Invalid range (start > end)
                    else:
                        return None  # Invalid range format
                else:
                    # Handle single number
                    num = int(part)
                    if 1 <= num <= max_num:
                        indices.add(num - 1)  # Convert to 0-based
                    else:
                        return None  # Invalid number
            
            return sorted(list(indices))
            
        except ValueError:
            return None
    
    def _atom_meets_search_criteria(self, atom, search_params: SearchParameters) -> bool:
        """Check if atom meets the search-specific filtering criteria"""
        element = atom.element.upper()
        atom_name = atom.name.upper()
        
        # Check if atom is metal
        if element in METALS:
            # Metal atom - check against search_include_metals criteria
            if not search_params.search_include_metals:
                # Empty list means include all metals
                return True
            elif search_params.search_include_metals == ['NONE']:
                # Special marker for no metals
                return False
            else:
                # Check if this metal is in the include list
                return element in [m.upper() for m in search_params.search_include_metals]
        else:
            # Non-metal atom - check against search_include_nonmetals criteria
            if not search_params.search_include_nonmetals:
                # Empty list means include all non-metals
                return True
            elif search_params.search_include_nonmetals == ['NONE']:
                # Special marker for no non-metals
                return False
            else:
                # Check if element name OR atom name is in the include list
                nonmetal_criteria = [nm.upper() for nm in search_params.search_include_nonmetals]
                return element in nonmetal_criteria or atom_name in nonmetal_criteria
    
    def _add_complete_residues_to_site(self, site: RedoxSite, selected_residues: List[Tuple[str, int, str]], structure: Structure) -> int:
        """Add all atoms from complete residues to the site"""
        atoms_added = 0
        
        for chain_id, resid, icode in selected_residues:
            # Find the residue in structure
            for model in structure:
                for chain in model:
                    if chain.id == chain_id:
                        for residue in chain:
                            if residue.id[1] == resid and residue.id[2] == icode:
                                # Add all atoms from this residue
                                for atom in residue.get_atoms():
                                    site_atom = RedoxSiteAtom(
                                        chain=chain_id,
                                        resname=residue.resname,
                                        resid=resid,
                                        atom_name=atom.name,
                                        coords=tuple(round(x, 3) for x in atom.coord),
                                        element=atom.element,
                                        altloc=atom.altloc,
                                        insertion_code=icode,
                                        occupancy=atom.occupancy,
                                        bfactor=atom.bfactor,
                                        properties={'serial_number': atom.serial_number}
                                    )
                                    site.add_atom(site_atom)
                                    atoms_added += 1
                                break
                        break
                    
        return atoms_added
    
    def _define_bonds_interactively(self, site: RedoxSite):
        """Interactive bond definition between residue pairs"""
        print(f"\n=== BOND DEFINITION ===")
        
        # Clear any existing bonds - user defines all bonds manually
        site.bonds.clear()
        
        # Group atoms by residue
        residue_groups = {}
        for atom in site.atoms:
            res_key = (atom.chain, atom.resname, atom.resid, atom.insertion_code)
            if res_key not in residue_groups:
                residue_groups[res_key] = []
            residue_groups[res_key].append(atom)
        
        residue_list = list(residue_groups.keys())
        
        if len(residue_list) < 2:
            print("Need at least 2 residues to define inter-residue bonds.")
            return
        
        print("Current site residues:")
        # Build table data for display later
        
        # Create table for residues first
        residues_table = Table(title="Current site residues")
        residues_table.add_column("#", style="cyan", width=6)
        residues_table.add_column("Residue", style="green")
        residues_table.add_column("Atoms", style="yellow")
        
        for i, res_key in enumerate(residue_list, 1):
            chain, resname, resid, icode = res_key
            atom_count = len(residue_groups[res_key])
            atom_types = sorted(set(atom.element for atom in residue_groups[res_key]))
            residues_table.add_row(
                f"[{i}]",
                f"{resname} {chain}:{resid}",
                f"{atom_count} atoms: {','.join(atom_types)}"
            )
        
        self.console.print(residues_table)

        self.console.print("\n[bold]Define bonds between residues:[/bold]")
        self.console.print("[grey50]1.[/grey50] Define bonds between residue pairs")
        self.console.print("[grey50]2.[/grey50] Skip bond definition")
        self.console.print("[grey50]v.[/grey50] View site in 3D structure viewer")

        while True:  # Loop to allow viewing and returning
            choice = prompt_with_context(
                processor=self.processor,
                prompt="[green]Choose option[/green]",
                default="1",
                module="Redox Detector - Bond Management",
                description="Define bonds between residues",
                options_map={"1": "Define bonds between residue pairs", "2": "Skip bond definition", "v": "View in 3D viewer", "view": "View in 3D viewer"}
            ).strip()

            if choice == "2":
                return
            elif choice.lower() == "v" or choice.lower() == "view":
                # Launch viewer showing both centers (with altloc) and added residues
                # We need to show centers + all other residues
                self._launch_viewer_for_site(site, "Current Site")
                continue  # Return to prompt after viewing
            elif choice == "1":
                break  # Continue to bond definition
        
        # Show transformer bond requirements if site has a mapped transformer
        self._display_transformer_bond_requirements(site)

        # Get residue pairs from user
        while True:
            self.console.print(f"\n[bold underline]RESIDUE PAIR BOND DEFINITION[/bold underline]")
            
            # Show available residues again
            self.console.print(residues_table)
            
            self.console.print(
                f"\n[bold]Enter residue pairs to define bonds[/bold] by residue ID "
                f"(e.g., [cyan]202-80[/cyan], [cyan]MN202-GLU80[/cyan], or [cyan]A:202-A:203[/cyan])")
            self.console.print(
                f"[grey50]Table row numbers (1-{len(residue_list)}) also work.[/grey50]")

            pairs_input = prompt_with_context(
                processor=self.processor,
                prompt="[green]Enter residue pairs[/green] (comma or space-separated, or 'done')",
                module="Redox Detector - Bond Management",
                description="Enter residue pairs for bonds",
                options_map={"done": "Finish defining bonds", "custom": "Residue pairs (e.g., 1-2, 1-3, 2-3)"}
            ).strip()
            
            if pairs_input.upper() == 'DONE':
                return
            
            # Parse residue pairs (handle both commas and spaces)
            residue_pairs = []
            try:
                # Split by both commas and spaces, remove empty strings
                pair_strings = [p.strip() for p in pairs_input.replace(',', ' ').split() if p.strip()]
                
                for pair_str in pair_strings:
                    if '-' not in pair_str:
                        print(f"Invalid pair format: {pair_str}. Use format 1-2.")
                        continue
                    
                    parts = pair_str.split('-')
                    if len(parts) != 2:
                        print(f"Invalid pair format: {pair_str}. Use format 1-2.")
                        continue
                    
                    first_token, second_token = parts[0].strip(), parts[1].strip()

                    # Resolve each side to a residue-list index. Accepts residue
                    # IDs (202, MN202, A:202) with a fallback to table row number.
                    first_idx = resolve_bond_residue_token(first_token, residue_list, on_error=print)
                    second_idx = resolve_bond_residue_token(second_token, residue_list, on_error=print)
                    if first_idx is None or second_idx is None:
                        continue

                    # Avoid duplicate pairs (1-2 same as 2-1)
                    pair_tuple = tuple(sorted([first_idx, second_idx]))
                    if pair_tuple not in [tuple(sorted([p[0], p[1]])) for p in residue_pairs]:
                        residue_pairs.append((first_idx, second_idx))
                
                if residue_pairs:
                    # Capture bond pairs in template mode
                    if self.template_mode and self.current_template is not None:
                        bond_pair_strings = []
                        for first_idx, second_idx in residue_pairs:
                            bond_pair_strings.append(f"{first_idx + 1}-{second_idx + 1}")  # Convert to 1-based
                        self.current_template.bond_pairs = bond_pair_strings
                    break
                else:
                    print("No valid pairs entered. Please try again.")
                    
            except Exception as e:
                print(f"Error parsing input: {e}. Please use format: A-B C-D")
        
        # Group pairs by source residue for efficient bond definition
        source_to_targets = {}
        for first_idx, second_idx in residue_pairs:
            # Treat the first residue in each pair as source, second as target (unidirectional)
            if first_idx not in source_to_targets:
                source_to_targets[first_idx] = []
            source_to_targets[first_idx].append(second_idx)
        
        # Process each source residue
        for source_idx, target_indices in source_to_targets.items():
            self._define_bonds_from_source_residue(site, residue_groups, residue_list, source_idx, target_indices)
        
        total_bonds = len(site.bonds)
        self.console.print(f"[green]Bond definition complete. Defined {total_bonds} bonds.[/green]")
        
        # Check if user wants to add more bonds
        while True:
#           add_more = prompt_with_context(None,
#               "\nWould you like to add more bonds or edit existing ones?",
#               choices=["yes", "no", "y", "n"],
#               default="no"
#           )
            
            add_more = "n"
            if add_more.lower() in ["yes", "y"]:
                # Show current bonds
                if site.bonds:
                    self.console.print("\n[bold]Current bonds:[/bold]")
                    for i, bond in enumerate(site.bonds, 1):
                        atom1_info = site.coord_to_pdb.get(bond.atom1_coords, {})
                        atom2_info = site.coord_to_pdb.get(bond.atom2_coords, {})
                        self.console.print(f"  {i}. {atom1_info.get('resname', '?')} {atom1_info.get('chain', '?')}:{atom1_info.get('resid', '?')} {atom1_info.get('atom_name', '?')} ↔ {atom2_info.get('resname', '?')} {atom2_info.get('chain', '?')}:{atom2_info.get('resid', '?')} {atom2_info.get('atom_name', '?')} ({bond.chemical_type})")
                
                # Allow defining more bonds using the same process
                self.console.print(f"\n[bold cyan]=== ADDITIONAL BOND DEFINITION ===[/bold cyan]")
                self._display_current_site_residues(residue_groups, residue_list)
                pairs_input = prompt_with_context(
                    processor=self.processor,
                    prompt="\nEnter additional residue pairs (e.g., 1-2, 1-3) or 'cancel'",
                    module="Redox Detector - Bond Management",
                    description="Enter additional residue pairs",
                    options_map={"cancel": "Cancel additional bonds", "custom": "Residue pairs (e.g., 1-2, 1-3)"}
                ).strip()
                
                if pairs_input.lower() == 'cancel':
                    break
                    
                # Parse and process additional pairs
                source_to_targets = self._parse_residue_pairs(pairs_input, residue_list)
                if source_to_targets:
                    for source_idx, target_indices in source_to_targets.items():
                        self._define_bonds_from_source_residue(site, residue_groups, residue_list, source_idx, target_indices)
                    
                    total_bonds = len(site.bonds)
                    self.console.print(f"[green]Now have {total_bonds} total bonds.[/green]")
                else:
                    break
            else:
                break
        
        # Display comprehensive site summary
        self._display_site_summary(site)
    
    def _define_bonds_from_source_residue(self, site: RedoxSite, residue_groups, residue_list, source_idx, target_indices):
        """Define bonds from a source residue to all its target residues"""
        source_key = residue_list[source_idx]
        source_atoms = residue_groups[source_key]
        chain_src, resname_src, resid_src, icode_src = source_key
        
        # Collect all target atoms
        target_atoms = []
        for target_idx in target_indices:
            target_key = residue_list[target_idx]
            target_residue_atoms = residue_groups[target_key]
            chain_tgt, resname_tgt, resid_tgt, icode_tgt = target_key
            
            for atom in target_residue_atoms:
                target_atoms.append({
                    'atom': atom,
                    'residue_info': f"{resname_tgt} {chain_tgt}:{resid_tgt}",
                    'target_idx': target_idx
                })
        
        self.console.print(f"\n[bold underline]DEFINING BONDS FROM {resname_src} {chain_src}:{resid_src}[/bold underline]")
        target_residue_names = [residue_list[idx][1] + ' ' + residue_list[idx][0] + ':' + str(residue_list[idx][2]) for idx in target_indices]
        self.console.print(f"Target residues: [cyan]{', '.join(target_residue_names)}[/cyan]")
        
        # Bond definition loop for this source residue
        while True:
            # Show atoms with bond status in Rich table
            source_table = Table(title=f"Atoms in {resname_src} {chain_src}:{resid_src}")
            source_table.add_column("#", style="cyan", width=6)
            source_table.add_column("Atom", style="green", width=12)
            source_table.add_column("Status", style="yellow")
            
            for i, atom in enumerate(source_atoms, 1):
                # Check if this atom already has bonds
                existing_bonds = [b for b in site.bonds if 
                                (b.atom1_coords == atom.coords or b.atom2_coords == atom.coords)]
                if existing_bonds:
                    bond_partners = []
                    for bond in existing_bonds:
                        # Determine which atom is the partner
                        if bond.atom1_coords == atom.coords:
                            partner_info = bond.atom2_residue_info
                        else:
                            partner_info = bond.atom1_residue_info
                        bond_partners.append(f"{partner_info['chain']}:{partner_info['resname']}{partner_info['resid']}:{partner_info['atom_name']}")
                    status = f"bonded to {', '.join(bond_partners)}"
                else:
                    status = "available"
                
                source_table.add_row(f"[{i}]", atom.atom_name, status)

            self.console.print(source_table)

            # Refresh the yellow halo on every bonded atom across the
            # whole site — both ends of each defined bond. Re-fires every
            # iteration of this inner loop, so the halo grows in sync
            # with the table's Status column as new bonds are created.
            _highlight_site_bonds(getattr(self, 'processor', None), site)
            
            # Select source atom
            try:
                choice = prompt_with_context(
                    processor=self.processor,
                    prompt=f"[green]Select source atom from {resname_src} {chain_src}:{resid_src}[/green] (1-{len(source_atoms)}, atom name, or 'done')",
                    module="Redox Detector - Bond Management",
                    description=f"Select source atom from {resname_src}",
                    options_map={"done": "Finish selecting atoms", "custom": "Atom number or name"}
                ).strip()
                if choice.lower() == 'done':
                    return
                
                selected_source_atom = None
                
                # Try to parse as number first
                try:
                    source_atom_idx = int(choice) - 1
                    if 0 <= source_atom_idx < len(source_atoms):
                        selected_source_atom = source_atoms[source_atom_idx]
                except ValueError:
                    pass
                
                # If not a valid number, try to find by atom name
                if selected_source_atom is None:
                    matching_atoms = [(i, atom) for i, atom in enumerate(source_atoms) if atom.atom_name.upper() == choice.upper()]
                    if len(matching_atoms) == 1:
                        selected_source_atom = matching_atoms[0][1]
                    elif len(matching_atoms) > 1:
                        self.console.print(f"[yellow]Multiple atoms named '{choice}' found:[/yellow]")
                        for idx, (atom_idx, atom) in enumerate(matching_atoms, 1):
                            original_number = atom_idx + 1
                            self.console.print(f"  {idx}. [{original_number}] {atom.atom_name}")
                        
                        try:
                            sub_choice = prompt_with_context(
                                processor=self.processor,
                                prompt=f"[green]Select which '{choice}' atom[/green] (1-{len(matching_atoms)})",
                                module="Redox Detector - Bond Management",
                                description=f"Select which '{choice}' atom",
                                options_map={str(i+1): f"Atom {i+1}" for i in range(len(matching_atoms))}
                            ).strip()
                            sub_idx = int(sub_choice) - 1
                            if 0 <= sub_idx < len(matching_atoms):
                                selected_source_atom = matching_atoms[sub_idx][1]
                            else:
                                self.console.print(f"[red]Invalid choice. Please use 1-{len(matching_atoms)}.[/red]")
                                continue
                        except ValueError:
                            self.console.print(f"[red]Invalid input. Please enter a number.[/red]")
                            continue
                    else:
                        self.console.print(f"[red]Invalid choice '{choice}'. Please use 1-{len(source_atoms)} or atom name.[/red]")
                        continue
                
                # Check if this source atom already has bonds defined
                existing_bonds = []
                for bond in site.bonds:
                    # Check if source atom matches either end of existing bond
                    if (bond.atom1_residue_info['chain'] == selected_source_atom.chain and
                        bond.atom1_residue_info['resname'] == selected_source_atom.resname and
                        bond.atom1_residue_info['resid'] == selected_source_atom.resid and
                        bond.atom1_residue_info['atom_name'] == selected_source_atom.atom_name):
                        existing_bonds.append(bond)
                    elif (bond.atom2_residue_info['chain'] == selected_source_atom.chain and
                          bond.atom2_residue_info['resname'] == selected_source_atom.resname and
                          bond.atom2_residue_info['resid'] == selected_source_atom.resid and
                          bond.atom2_residue_info['atom_name'] == selected_source_atom.atom_name):
                        existing_bonds.append(bond)
                
                # If bonds exist, ask user what to do
                if existing_bonds:
                    self.console.print(f"\n[yellow]Source atom {selected_source_atom.atom_name} already has {len(existing_bonds)} bond(s):[/yellow]")
                    for i, bond in enumerate(existing_bonds, 1):
                        # Determine which atom is the partner
                        if (bond.atom1_residue_info['chain'] == selected_source_atom.chain and
                            bond.atom1_residue_info['atom_name'] == selected_source_atom.atom_name):
                            partner = bond.atom2_residue_info
                        else:
                            partner = bond.atom1_residue_info
                        self.console.print(f"  {i}. {partner['resname']} {partner['chain']}:{partner['resid']} {partner['atom_name']}")
                    
                    self.console.print("\n[bold]Options:[/bold]")
                    self.console.print("[grey50]1.[/grey50] Delete existing bonds and create new one")
                    self.console.print("[grey50]2.[/grey50] Replace one existing bond")
                    self.console.print("[grey50]3.[/grey50] Add another bond (keep existing)")
                    self.console.print("[grey50]4.[/grey50] Cancel (keep existing bonds unchanged)")

                    action_choice = prompt_with_context(
                        processor=self.processor,
                        prompt="[green]Choose action[/green]",
                        choices=["1", "2", "3", "4"],
                        default="4",
                        module="Redox Detector - Bond Management",
                        description="Handle existing bonds",
                        options_map={
                            "1": "Delete existing bonds and create new one",
                            "2": "Replace one existing bond",
                            "3": "Add another bond (keep existing)",
                            "4": "Cancel (keep existing bonds unchanged)"
                        }
                    ).strip()
                    
                    if action_choice == "1":
                        # Delete all existing bonds for this atom and return to source selection
                        for bond in existing_bonds:
                            site.bonds.remove(bond)
                        self.console.print(f"[yellow]Deleted {len(existing_bonds)} existing bond(s). Select a new source atom.[/yellow]")
                        continue
                    elif action_choice == "2":
                        # Replace specific bond
                        if len(existing_bonds) == 1:
                            site.bonds.remove(existing_bonds[0])
                            self.console.print("[yellow]Replaced existing bond[/yellow]")
                        else:
                            bond_choice = prompt_with_context(
                                processor=self.processor,
                                prompt=f"[green]Which bond to replace?[/green] (1-{len(existing_bonds)})",
                                choices=[str(i) for i in range(1, len(existing_bonds)+1)],
                                module="Redox Detector - Bond Management",
                                description="Select bond to replace",
                                options_map={str(i): f"Bond {i}" for i in range(1, len(existing_bonds)+1)}
                            )
                            bond_to_replace = existing_bonds[int(bond_choice) - 1]
                            site.bonds.remove(bond_to_replace)
                            self.console.print(f"[yellow]Replaced bond #{bond_choice}[/yellow]")
                    elif action_choice == "3":
                        # Add another bond - no action needed, will proceed to target selection
                        self.console.print("[cyan]Adding additional bond...[/cyan]")
                    elif action_choice == "4":
                        # Cancel - go back to source atom selection
                        self.console.print("[yellow]Cancelled - keeping existing bonds[/yellow]")
                        continue
                
                # Calculate distances from selected source atom to all target atoms
                target_atoms_with_distances = []
                for target_info in target_atoms:
                    target_atom = target_info['atom']
                    distance = np.linalg.norm(
                        np.array(selected_source_atom.coords) - np.array(target_atom.coords)
                    )
                    target_atoms_with_distances.append({
                        **target_info,
                        'distance': distance
                    })
                
                # Sort by distance
                target_atoms_with_distances.sort(key=lambda x: x['distance'])
                
                # Group atoms by residue while preserving distance order
                residue_groups = {}
                residue_order = []
                for target_info in target_atoms_with_distances:
                    residue_info = target_info['residue_info']
                    if residue_info not in residue_groups:
                        residue_groups[residue_info] = []
                        residue_order.append(residue_info)
                    residue_groups[residue_info].append(target_info)
                
                # Create indexed list matching the display order
                display_ordered_atoms = []
                
                # Display target atoms in Rich table grouped by residue but ordered by closest atom distance
                target_table = Table(title=f"Target atoms (ordered by distance to {selected_source_atom.atom_name} {selected_source_atom.element})")
                target_table.add_column("#", style="cyan", width=6)
                target_table.add_column("Residue", style="green", width=12)
                target_table.add_column("Atom", style="yellow", width=8)
                target_table.add_column("Element", style="magenta", width=8)
                target_table.add_column("Distance", style="blue", width=10)
                
                i = 1
                for residue_idx, residue_info in enumerate(residue_order):
                    atoms_in_residue = residue_groups[residue_info]
                    
                    # Add separator line before each residue (except the first one)
                    if residue_idx > 0:
                        target_table.add_section()
                    
                    for target_info in atoms_in_residue:
                        atom = target_info['atom']
                        distance = target_info['distance']
                        target_table.add_row(
                            f"[{i}]",
                            residue_info,
                            atom.atom_name,
                            atom.element,
                            f"{distance:.2f}Å"
                        )
                        display_ordered_atoms.append(target_info)
                        i += 1
                
                self.console.print(target_table)
                
                # Select target atom
                try:
                    target_choice = prompt_with_context(
                        processor=self.processor,
                        prompt=f"[green]Select target atom(s)[/green] (1-{len(display_ordered_atoms)}, atom names, comma-separated for multiple, or 'cancel')",
                        module="Redox Detector - Bond Management",
                        description="Select target atom(s) for bonding",
                        options_map={"cancel": "Cancel bond creation", "custom": "Atom numbers or names (comma-separated)"}
                    ).strip()
                    if target_choice.lower() == 'cancel':
                        continue
                        
                    # Parse multiple selections (numbers or atom names, but not mixed)
                    target_indices = []
                    choices = [choice.strip() for choice in target_choice.split(',') if choice.strip()]
                    
                    if not choices:
                        continue
                    
                    # Determine if input is all numbers or all atom names
                    is_numeric = []
                    for choice in choices:
                        try:
                            int(choice)
                            is_numeric.append(True)
                        except ValueError:
                            is_numeric.append(False)
                    
                    # Check for mixed format
                    if any(is_numeric) and not all(is_numeric):
                        self.console.print(f"[red]Please use either all numbers (1,2,3) or all atom names (FE,ND1,SG), not mixed.[/red]")
                        continue
                    
                    # Process based on format
                    if all(is_numeric):
                        # Handle numeric input
                        for choice in choices:
                            target_atom_idx = int(choice) - 1
                            if 0 <= target_atom_idx < len(display_ordered_atoms):
                                target_indices.append(target_atom_idx)
                            else:
                                self.console.print(f"[red]Invalid number: {choice}. Please use 1-{len(display_ordered_atoms)}.[/red]")
                                target_indices = []
                                break
                    else:
                        # Handle atom name input
                        for choice in choices:
                            matching_atoms = [(i, info) for i, info in enumerate(display_ordered_atoms) 
                                            if info['atom'].atom_name.upper() == choice.upper()]
                            
                            if len(matching_atoms) == 1:
                                target_indices.append(matching_atoms[0][0])
                            elif len(matching_atoms) > 1:
                                self.console.print(f"[yellow]Multiple atoms named '{choice}' found:[/yellow]")
                                for idx, (atom_idx, info) in enumerate(matching_atoms, 1):
                                    original_number = atom_idx + 1
                                    residue_info = info['residue_info']
                                    distance = info['distance']
                                    self.console.print(f"  {idx}. [{original_number}] {info['atom'].atom_name} in {residue_info} ({distance:.2f}Å)")
                                
                                n_match = len(matching_atoms)
                                sub_choice = prompt_with_context(
                                    processor=self.processor,
                                    prompt=f"[green]Select which '{choice}' atom(s)[/green] (1-{n_match}; comma/space-separated and/or ranges like 1-{n_match}, or 'all')",
                                    module="Redox Detector - Bond Management",
                                    description=f"Select which '{choice}' atom(s)",
                                    options_map={str(i+1): f"Atom {i+1}" for i in range(n_match)}
                                ).strip()
                                # Accept numbers, comma/space-separated lists, ranges
                                # (e.g. 1-4), and 'all'. A range like '1-4' previously hit
                                # int() and was rejected as invalid even though the prompt
                                # showed '1-4' as the valid bound.
                                sub_nums, ok = [], True
                                if sub_choice.lower() == 'all':
                                    sub_nums = list(range(1, n_match + 1))
                                else:
                                    for tok in sub_choice.replace(',', ' ').split():
                                        if '-' in tok[1:]:  # range like 1-4 (not a leading minus)
                                            a, _, b = tok.partition('-')
                                            if a.isdigit() and b.isdigit():
                                                lo, hi = int(a), int(b)
                                                sub_nums.extend(range(lo, hi + 1) if lo <= hi else range(hi, lo + 1))
                                            else:
                                                ok = False
                                                break
                                        elif tok.isdigit():
                                            sub_nums.append(int(tok))
                                        else:
                                            ok = False
                                            break
                                if not ok or not sub_nums:
                                    self.console.print(f"[red]Invalid input. Use numbers, ranges (e.g. 1-{n_match}), or 'all'.[/red]")
                                    target_indices = []
                                    break
                                for s in sub_nums:
                                    if 1 <= s <= n_match:
                                        idx0 = matching_atoms[s - 1][0]
                                        if idx0 not in target_indices:
                                            target_indices.append(idx0)
                                    else:
                                        self.console.print(f"[red]Invalid choice: {s}. Please use 1-{n_match}.[/red]")
                                        target_indices = []
                                        break
                                if not target_indices:
                                    break
                            else:
                                self.console.print(f"[red]No atom named '{choice}' found.[/red]")
                                target_indices = []
                                break
                    
                    if not target_indices:
                        continue  # Skip if there were invalid selections
                    
                    # Create bonds for each selected target atom
                    bonds_created = 0
                    for target_atom_idx in target_indices:
                        selected_target_info = display_ordered_atoms[target_atom_idx]
                        selected_target_atom = selected_target_info['atom']
                        
                        # Capture template data if in template mode
                        if self.template_mode and self.current_template is not None:
                            # Find which residue pair this bond belongs to
                            source_res_idx = next(i for i, key in enumerate(residue_list) if key == source_key)
                            target_res_idx = next(i for i, key in enumerate(residue_list) if key == (selected_target_atom.chain, selected_target_atom.resname, selected_target_atom.resid, selected_target_atom.insertion_code))
                            
                            bond_pair_key = f"{source_res_idx + 1}-{target_res_idx + 1}"  # Convert to 1-based
                            
                            if bond_pair_key not in self.current_template.bond_selections:
                                self.current_template.bond_selections[bond_pair_key] = []
                            self.current_template.bond_selections[bond_pair_key].append({
                                'source_atom': selected_source_atom.atom_name,
                                'target_atom': selected_target_atom.atom_name
                            })
                        
                        # Classify bond properly. Pass atom names so SG-SG
                        # between CYS is recognised as a disulfide (the
                        # disulfide test in classify_bond_types keys off
                        # atom1_name/atom2_name == 'SG'); omitting them made
                        # interactively-defined disulfides fall through to
                        # "covalent", inconsistent with the template-apply path.
                        bond_type, chemical_type = classify_bond_types(
                            atom1_element=selected_source_atom.element,
                            atom2_element=selected_target_atom.element,
                            distance=selected_target_info['distance'],
                            atom1_residue=selected_source_atom.resname,
                            atom2_residue=selected_target_atom.resname,
                            atom1_resid=selected_source_atom.resid,
                            atom2_resid=selected_target_atom.resid,
                            atom1_chain=selected_source_atom.chain,
                            atom2_chain=selected_target_atom.chain,
                            atom1_name=selected_source_atom.atom_name,
                            atom2_name=selected_target_atom.atom_name
                        )
                        
                        # Create the bond
                        bond = RedoxSiteBond(
                            atom1_coords=selected_source_atom.coords,
                            atom2_coords=selected_target_atom.coords,
                            bond_type=bond_type,
                            chemical_type=chemical_type,
                            distance=selected_target_info['distance'],
                            atom1_element=selected_source_atom.element,
                            atom2_element=selected_target_atom.element,
                            atom1_residue_info={
                                'chain': selected_source_atom.chain,
                                'resname': selected_source_atom.resname,
                                'resid': selected_source_atom.resid,
                                'atom_name': selected_source_atom.atom_name
                            },
                            atom2_residue_info={
                                'chain': selected_target_atom.chain,
                                'resname': selected_target_atom.resname,
                                'resid': selected_target_atom.resid,
                                'atom_name': selected_target_atom.atom_name
                            }
                        )
                        
                        site.bonds.append(bond)
                        self.console.print(f"[green]Created bond: {selected_source_atom.atom_name} {selected_source_atom.element} ↔ {selected_target_atom.atom_name} {selected_target_atom.element} ({selected_target_info['distance']:.2f}Å)[/green]")
                        bonds_created += 1
                    
                    if bonds_created > 1:
                        self.console.print(f"[green]Created {bonds_created} bonds total.[/green]")
                    
                except ValueError:
                    self.console.print("[red]Invalid input. Please enter a number.[/red]")
                    continue
                    
            except ValueError:
                self.console.print("[red]Invalid input. Please enter a number.[/red]")
                continue
    
    def _define_bonds_between_residues(self, site: RedoxSite, residue_groups, residue_list, first_idx, second_idx):
        """Define bonds between atoms in two specific residues"""
        res1_key = residue_list[first_idx]
        res2_key = residue_list[second_idx]
        res1_atoms = residue_groups[res1_key]
        res2_atoms = residue_groups[res2_key]
        
        chain1, resname1, resid1, icode1 = res1_key
        chain2, resname2, resid2, icode2 = res2_key
        
        print(f"\nResidue pair: {resname1} {chain1}:{resid1} ↔ {resname2} {chain2}:{resid2}")
        
        # Display atoms in each residue
        print(f"\nAtoms in {resname1} {chain1}:{resid1}:")
        for i, atom in enumerate(res1_atoms, 1):
            print(f"  [{i}] {atom.atom_name} {atom.element}")
        
        print(f"\nAtoms in {resname2} {chain2}:{resid2}:")
        for i, atom in enumerate(res2_atoms, 1):
            print(f"  [{i}] {atom.atom_name} {atom.element}")
        
        # Bond definition loop for this pair
        while True:
            # Select atom from first residue
            while True:
                try:
                    _a1_map = {str(i): f"{a.atom_name} {a.element}" for i, a in enumerate(res1_atoms, 1)}
                    _a1_map["done"] = "Done (finish this pair)"
                    atom1_choice = prompt_with_context(
                        self.processor,
                        f"Select atom from {resname1} {chain1}:{resid1} (1-{len(res1_atoms)} or 'done')",
                        module="Redox Bond Editor",
                        description="Select first atom for cross-residue bond",
                        options_map=_a1_map,
                    ).strip()
                    if atom1_choice.lower() == 'done':
                        return
                    atom1_idx = int(atom1_choice) - 1
                    if 0 <= atom1_idx < len(res1_atoms):
                        break
                    print(f"Invalid choice. Please use 1-{len(res1_atoms)}.")
                except ValueError:
                    print("Invalid input. Please enter a number.")
            
            # Select atom from second residue
            while True:
                try:
                    _a2_map = {str(i): f"{a.atom_name} {a.element}" for i, a in enumerate(res2_atoms, 1)}
                    _a2_map["cancel"] = "Cancel this bond"
                    atom2_choice = prompt_with_context(
                        self.processor,
                        f"Select atom from {resname2} {chain2}:{resid2} (1-{len(res2_atoms)} or 'cancel')",
                        module="Redox Bond Editor",
                        description="Select second atom for cross-residue bond",
                        options_map=_a2_map,
                    ).strip()
                    if atom2_choice.lower() == 'cancel':
                        break
                    atom2_idx = int(atom2_choice) - 1
                    if 0 <= atom2_idx < len(res2_atoms):
                        break
                    print(f"Invalid choice. Please use 1-{len(res2_atoms)}.")
                except ValueError:
                    print("Invalid input. Please enter a number.")
            
            if atom2_choice.lower() == 'cancel':
                continue
            
            # Create the bond
            atom1 = res1_atoms[atom1_idx]
            atom2 = res2_atoms[atom2_idx]
            distance = np.linalg.norm(np.array(atom1.coords) - np.array(atom2.coords))
            
            # Auto-determine bond types
            bond_type = "intraresidue" if (chain1 == chain2 and resid1 == resid2 and icode1 == icode2) else "interresidue"
            
            # Auto-determine chemical type using METALS
            elem1, elem2 = atom1.element.upper(), atom2.element.upper()
            is_metal1 = elem1 in METALS
            is_metal2 = elem2 in METALS
            
            # Use the proper classification function
            bond_type, chemical_type = classify_bond_types(
                atom1_element=atom1.element,
                atom2_element=atom2.element,
                distance=distance,
                atom1_residue=resname1,
                atom2_residue=resname2,
                atom1_resid=resid1,
                atom2_resid=resid2,
                atom1_chain=chain1,
                atom2_chain=chain2,
                atom1_name=atom1.atom_name,
                atom2_name=atom2.atom_name
            )
            
            print(f"\nBond: {resname1} {chain1}:{resid1} {atom1.atom_name} ↔ {resname2} {chain2}:{resid2} {atom2.atom_name} ({distance:.2f}Å)")
            print(f"Chemical type: {chemical_type} (auto-detected: {'metal-nonmetal' if chemical_type=='coordinate' else 'nonmetal-nonmetal' if chemical_type=='covalent' else 'disulfide bridge' if chemical_type=='disulfide' else 'metal-metal' if chemical_type=='metal-metal' else 'other'})")
            print(f"Bond type: {bond_type} (auto-detected: {'same residue' if bond_type=='intraresidue' else 'different residues'})")
            
            # Create and add bond
            bond = RedoxSiteBond(
                atom1_coords=atom1.coords,
                atom2_coords=atom2.coords,
                bond_type=bond_type,
                chemical_type=chemical_type,
                distance=distance,
                atom1_element=atom1.element,
                atom2_element=atom2.element,
                atom1_residue_info=site.coord_to_pdb[atom1.coords],
                atom2_residue_info=site.coord_to_pdb[atom2.coords]
            )
            site.bonds.append(bond)
            
            # Ask if user wants to add another bond for this pair
            another_bond = prompt_with_context(
                self.processor,
                "Add another bond for this residue pair? [n]",
                default="n",
                module="Redox Bond Editor",
                description="Add another cross-residue bond?",
            ).strip().lower()
            if not another_bond.startswith('y'):
                break
    
    def _recompute_site_bonds(self, site: RedoxSite, structure: Structure):
        """Clear bonds - user will define bonds manually"""
        # Clear all bonds - no automatic bond detection
        site.bonds.clear()
        # Note: Bonds are now defined interactively by user
    
    def _display_site_summary(self, site: RedoxSite):
        """Display comprehensive site summary table"""
        
        # Create the main summary table
        summary_table = Table(title=f"Site Summary: {site.site_id}")
        summary_table.add_column("Content", style="white", no_wrap=False)
        
        # Section 1: Redox Centers
        summary_table.add_row(f"[bold cyan]Redox Centers ({len(site.centers)})[/bold cyan]")
        if site.centers:
            for center in site.centers:
                location = f"{center.resname} {center.chain}:{center.resid}"
                if center.altloc:
                    location += center.altloc
                if center.atom_name:
                    location += f" {center.atom_name}"
                center_type = center.center_type.value.replace('_', ' ')
                summary_table.add_row(f"  {location} ({center_type})")
        else:
            summary_table.add_row("  No redox centers defined")
        
        # Add separator
        summary_table.add_section()
        
        # Section 2: Site Residues
        residue_groups = {}
        for atom in site.atoms:
            res_key = (atom.chain, atom.resname, atom.resid)
            if res_key not in residue_groups:
                residue_groups[res_key] = []
            residue_groups[res_key].append(atom)
        
        summary_table.add_row(f"[bold green]Site Residues ({len(residue_groups)})[/bold green]")
        if residue_groups:
            for (chain, resname, resid), atoms in residue_groups.items():
                atom_count = len(atoms)
                elements = sorted(set(atom.element for atom in atoms))
                elements_str = ", ".join(elements)
                summary_table.add_row(f"  {resname} {chain}:{resid} ({atom_count} atoms: {elements_str})")
        else:
            summary_table.add_row("  No residues in site")
        
        # Add separator
        summary_table.add_section()
        
        # Section 3: Site Bonds
        summary_table.add_row(f"[bold yellow]Site Bonds ({len(site.bonds)})[/bold yellow]")
        if site.bonds:
            for bond in site.bonds:
                atom1_info = bond.atom1_residue_info
                atom2_info = bond.atom2_residue_info
                bond_desc = (f"  {atom1_info['resname']} {atom1_info['chain']}:{atom1_info['resid']} {atom1_info['atom_name']} "
                           f"↔ {atom2_info['resname']} {atom2_info['chain']}:{atom2_info['resid']} {atom2_info['atom_name']} "
                           f"({bond.chemical_type}, {bond.distance:.2f}Å)")
                summary_table.add_row(bond_desc)
        else:
            summary_table.add_row("  No bonds defined")
        
        self.console.print("\n")
        self.console.print(summary_table)
    
    def refine_site_for_template(self, site: RedoxSite, structure: Structure, site_type: str) -> Tuple[RedoxSite, Optional[SiteTemplate]]:
        """Refine site and capture template for batch processing"""
        # Perform normal interactive refinement but capture template
        if not self.template_mode:
            self.console.print("[yellow]Warning: refine_site_for_template called without template_mode=True[/yellow]")
            
        # Set template capture mode
        self.current_template = None
        template_capture_mode = True
        
        # Refine site normally but capture choices
        refined_site = self.refine_site_interactively(site, structure)
        
        # Template should have been captured during refinement
        template = self.current_template
        if template:
            template.site_type = site_type
            
        return refined_site, template
        
    def _apply_template_bonds(self, site: RedoxSite, template: SiteTemplate):
        """Apply bond patterns from template to site"""
        
        # Group atoms by residue for bond matching
        residue_groups = {}
        for atom in site.atoms:
            res_key = (atom.chain, atom.resname, atom.resid, atom.insertion_code)
            if res_key not in residue_groups:
                residue_groups[res_key] = []
            residue_groups[res_key].append(atom)
        
        residue_list = list(residue_groups.keys())
        bonds_created = 0
        
        # Apply each captured bond selection
        for bond_pair, bond_list in template.bond_selections.items():
            try:
                # Parse bond pair (e.g., "1-2")
                res1_idx, res2_idx = map(int, bond_pair.split('-'))
                res1_idx -= 1  # Convert to 0-based
                res2_idx -= 1
                
                if res1_idx >= len(residue_list) or res2_idx >= len(residue_list):
                    self.console.print(f"[yellow]Skipping bond pair {bond_pair}: residue index out of range[/yellow]")
                    continue
                    
                source_key = residue_list[res1_idx]
                target_key = residue_list[res2_idx]
                
                # Apply each bond in this bond pair
                for bond_info in bond_list:
                    # Find source atom by name
                    source_atoms = residue_groups[source_key]
                    source_atom = next((a for a in source_atoms if a.atom_name == bond_info['source_atom']), None)
                    
                    # Find target atom by name
                    target_atoms = residue_groups[target_key]
                    target_atom = next((a for a in target_atoms if a.atom_name == bond_info['target_atom']), None)
                    
                    if source_atom and target_atom:
                        # Calculate distance
                        distance = np.linalg.norm(
                            np.array(source_atom.coords) - np.array(target_atom.coords)
                        )
                        
                        # Classify bond properly
                        bond_type, chemical_type = classify_bond_types(
                            atom1_element=source_atom.element,
                            atom2_element=target_atom.element,
                            distance=distance,
                            atom1_residue=source_atom.resname,
                            atom2_residue=target_atom.resname,
                            atom1_resid=source_atom.resid,
                            atom2_resid=target_atom.resid,
                            atom1_chain=source_atom.chain,
                            atom2_chain=target_atom.chain,
                            atom1_name=source_atom.atom_name,
                            atom2_name=target_atom.atom_name
                        )
                        
                        # Create bond using same format as interactive mode
                        bond = RedoxSiteBond(
                            atom1_coords=source_atom.coords,
                            atom2_coords=target_atom.coords,
                            bond_type=bond_type,
                            chemical_type=chemical_type,
                            distance=distance,
                            atom1_element=source_atom.element,
                            atom2_element=target_atom.element,
                            atom1_residue_info={
                                'chain': source_atom.chain,
                                'resname': source_atom.resname,
                                'resid': source_atom.resid,
                                'atom_name': source_atom.atom_name
                            },
                            atom2_residue_info={
                                'chain': target_atom.chain,
                                'resname': target_atom.resname,
                                'resid': target_atom.resid,
                                'atom_name': target_atom.atom_name
                            }
                        )
                        
                        site.bonds.append(bond)
                        bonds_created += 1
                        self.console.print(f"[green]Applied template bond: {source_atom.atom_name} ↔ {target_atom.atom_name} ({distance:.2f}Å)[/green]")
                    else:
                        missing_atoms = []
                        if not source_atom:
                            missing_atoms.append(f"source {bond_info['source_atom']}")
                        if not target_atom:
                            missing_atoms.append(f"target {bond_info['target_atom']}")
                        self.console.print(f"[yellow]Skipping bond {bond_pair}: {', '.join(missing_atoms)} not found[/yellow]")
                    
            except (ValueError, IndexError) as e:
                self.console.print(f"[yellow]Skipping invalid bond pattern {bond_pair}: {e}[/yellow]")
                continue
        
        self.console.print(f"[green]Template bonds applied: {bonds_created} bonds created[/green]")

    def apply_template_to_site(self, site: RedoxSite, structure: Structure, template: SiteTemplate) -> RedoxSite:
        """Apply template to site automatically"""
        if not template:
            raise ValueError("No template provided")

        self.console.print(f"[cyan]Applying template for {template.site_type} to {site.site_id}[/cyan]")

        # Narrow the viewer to this site too. ``refine_site_interactively``
        # has the same hook for the first site of each template type;
        # this matches it for the auto-applied subsequent sites.
        _narrow_viewer_to_site(getattr(self, 'processor', None), site)

        # Check if template is minimal (no residues to add and no distance search)
        if (not template.residue_types or sum(template.residue_counts) == 0) and template.search_constraint != "distance":
            self.console.print(f"[green]Template applied: minimal site (center only, no additional residues)[/green]")
            # Apply bond patterns from template (if any)
            if template.bond_pairs:
                self._apply_template_bonds(site, template)
            _highlight_site_bonds(getattr(self, 'processor', None), site)
            return site

        # Reconstruct custom boundary coords if needed
        custom_coords = None
        if template.boundary_definition == BoundaryDefinition.MIN_DISTANCE_CUSTOM and template.custom_boundary_atom_indices:
            # Use the stored indices to get coordinates from this site's atoms
            custom_coords = []
            for idx in template.custom_boundary_atom_indices:
                if idx < len(site.atoms):
                    custom_coords.append(site.atoms[idx].coords)
                else:
                    self.console.print(f"[yellow]Warning: Template boundary atom index {idx} out of range for site {site.site_id}[/yellow]")

        # Create SearchParameters from template
        if template.search_constraint == "distance":
            # Distance-based search
            if template.distance_method == "fixed":
                search_params = SearchParameters(
                    constraint=SearchConstraint.DISTANCE_CUTOFF,
                    distance_method=DistanceMethod.FIXED,
                    boundary_definition=template.boundary_definition,
                    radius=template.distance_radius,
                    custom_boundary_coords=custom_coords,
                    search_include_metals=template.atom_filtering_metals,
                    search_include_nonmetals=template.atom_filtering_nonmetals
                )
            else:
                search_params = SearchParameters(
                    constraint=SearchConstraint.DISTANCE_CUTOFF,
                    distance_method=DistanceMethod.ADAPTIVE,
                    boundary_definition=template.boundary_definition,
                    min_radius=template.distance_min_radius,
                    max_radius=template.distance_max_radius,
                    radius_step=template.distance_radius_step,
                    target_coordination=template.distance_target_coordination,
                    custom_boundary_coords=custom_coords,
                    search_include_metals=template.atom_filtering_metals,
                    search_include_nonmetals=template.atom_filtering_nonmetals
                )
            search_result = self._perform_distance_search(site, structure, search_params)
        else:
            # Count-based search
            search_params = SearchParameters(
                constraint=SearchConstraint.COUNT_CUTOFF,
                target_residue_types=template.residue_types,
                target_residue_count=sum(template.residue_counts),
                target_residue_count_map=dict(zip(template.residue_types, template.residue_counts)),
                boundary_definition=template.boundary_definition,
                custom_boundary_coords=custom_coords,
                search_include_metals=template.atom_filtering_metals,
                search_include_nonmetals=template.atom_filtering_nonmetals
            )
            search_result = self._perform_count_search(site, structure, search_params)

        if not search_result or search_result.total_residues_found == 0:
            raise ValueError("Template search found no residues")

        # Auto-select residues (use template's selection pattern)
        # Convert detected_residues to format expected by _add_complete_residues_to_site
        selected_residues = [(res['chain'], res['resid'], res['insertion_code'])
                           for res in search_result.detected_residues]

        # Add complete residues to site using existing method
        atoms_added = self._add_complete_residues_to_site(site, selected_residues, structure)

        # Apply bond patterns from template
        if template.bond_pairs:
            self._apply_template_bonds(site, template)

        self.console.print(f"[green]Template applied: added {atoms_added} atoms from {len(selected_residues)} residues[/green]")

        # Re-narrow now that the template has populated ``site.atoms`` —
        # the early narrow at the top of this method only saw the centers,
        # so the added residues weren't yet visible. Re-narrow then add
        # the bonded-atom halo on top.
        _narrow_viewer_to_site(getattr(self, 'processor', None), site)
        _highlight_site_bonds(getattr(self, 'processor', None), site)
        return site
    
    def _display_transformer_bond_requirements(self, site):
        """Display bond requirements for the transformer mapped to this site's type"""
        if not self.transformer_mappings:
            return
            
        site_type = getattr(site, 'site_type', None)
        if not site_type or site_type not in self.transformer_mappings:
            return
            
        transformer_name = self.transformer_mappings[site_type]

        registry = _get_transformer_registry()
        if not registry:
            return

        transformers = registry.get_all_transformers()
        transformer_class = transformers.get(transformer_name)
        
        if not transformer_class:
            return
            
        try:
            requirements = transformer_class.get_site_requirements()
            if 'bonds' not in requirements:
                return
                
            bond_info = requirements['bonds']
            self.console.print(f"\n[bold cyan]Bond Requirements for {transformer_name}:[/bold cyan]")
            
            # Handle new required_bond_groups structure
            if 'required_bond_groups' in bond_info:
                bond_groups = bond_info['required_bond_groups']
                require_one_group = bond_info.get('require_one_group', False)
                
                if require_one_group:
                    self.console.print(f"[grey50](At least one variant must be satisfied)[/grey50]")
                
                for group in bond_groups:
                    desc = group.get('description', 'Bond group')
                    min_count = group.get('min_count', 0)
                    bond_types = group.get('bond_types', {})
                    
                    self.console.print(f"• {desc} (min: {min_count} bonds):", markup=False)
                    
                    for bond_type, atom_pairs in bond_types.items():
                        self.console.print(f"  {bond_type.title()}:", markup=False)
                        for pair in atom_pairs:
                            (res1, atom1), (res2, atom2) = pair
                            bond_text = f"    - {res1}:{atom1} ↔ {res2}:{atom2}"
                            self.console.print(bond_text, markup=False, highlight=False)
            
            # Handle legacy required_bond_types structure for compatibility
            elif 'required_bond_types' in bond_info:
                for bond_type, details in bond_info.get('required_bond_types', {}).items():
                    desc = details.get('description', f'{bond_type} bonds')
                    min_count = details.get('min_count', 0)
                    
                    atom_pairs = details.get('atom_pairs', [])
                    if atom_pairs:
                        self.console.print(f"• {desc} (min: {min_count}):", markup=False)
                        for pair in atom_pairs:
                            (res1, atom1), (res2, atom2) = pair
                            bond_text = f"  - {res1}:{atom1} ↔ {res2}:{atom2}"
                            self.console.print(bond_text, markup=False, highlight=False)
                    else:
                        self.console.print(f"• {desc} (min: {min_count})", markup=False)
                        
        except (NotImplementedError, AttributeError):
            # Transformer doesn't implement get_site_requirements
            pass

    # ========================================================================
    # Structure Viewer Integration
    # ========================================================================

    def _launch_viewer_for_site(self, site: 'RedoxSite', title: str = "Site") -> bool:
        """
        Launch structure viewer showing entire site (centers + all residues).

        Args:
            site: RedoxSite object to visualize
            title: Title for the annotation

        Returns:
            True if viewer launched successfully
        """
        # Build NGL selection: centers (with altloc) + non-center residues.
        selections = []
        center_residues = set()
        if hasattr(site, 'centers') and site.centers:
            for center in site.centers:
                selection = f":{center.chain} and {center.resid}"
                if center.altloc and center.altloc.strip():
                    selection += f" and %{center.altloc.strip()}"
                selections.append(selection)
                center_residues.add((center.chain, center.resid))
        if hasattr(site, 'atoms') and site.atoms:
            seen_residues = set()
            for atom in site.atoms:
                res_key = (atom.chain, atom.resid)
                if res_key not in center_residues and res_key not in seen_residues:
                    selections.append(f":{atom.chain} and {atom.resid}")
                    seen_residues.add(res_key)

        ngl_selection = " or ".join(selections) if selections else ""
        return _show_residues_in_viewer(
            getattr(self, 'processor', None),
            self.console,
            ngl_selection,
            label='site_view',
        )

    def _launch_viewer_with_centers(self, centers: List[RedoxCenter], title: str = "Selected Centers") -> bool:
        """
        Launch structure viewer with specific centers highlighted.

        Args:
            centers: Centers to highlight in the viewer
            title: Title for the annotation

        Returns:
            True if viewer launched successfully
        """
        selections = []
        for center in centers:
            selection = f":{center.chain} and {center.resid}"
            if center.altloc and center.altloc.strip():
                selection += f" and %{center.altloc.strip()}"
            selections.append(selection)
        ngl_selection = " or ".join(selections) if selections else ""
        return _show_residues_in_viewer(
            getattr(self, 'processor', None),
            self.console,
            ngl_selection,
            label='selected_centers',
        )

    def _launch_viewer_for_residues(self, residue_list: List[Tuple], title: str = "Residues") -> bool:
        """
        Launch structure viewer with specific residues highlighted.

        Args:
            residue_list: List of (chain, resname, resid, icode) tuples
            title: Title for the annotation

        Returns:
            True if viewer launched successfully
        """
        selections = [f":{chain} and {resid}" for chain, resname, resid, icode in residue_list]
        ngl_selection = " or ".join(selections) if selections else ""
        return _show_residues_in_viewer(
            getattr(self, 'processor', None),
            self.console,
            ngl_selection,
            label='selected_residues',
        )

    def _launch_viewer_for_residues_with_centers(self, site: RedoxSite, residue_list: List[Tuple],
                                                 title: str = "Residues + Centers") -> bool:
        """
        Launch structure viewer showing both the site centers (for context) and found residues.
        This is used during residue selection to help users see spatial relationships.

        Args:
            site: RedoxSite containing the original centers
            residue_list: List of (chain, resname, resid, icode) tuples for found residues
            title: Title for the annotation

        Returns:
            True if viewer launched successfully
        """
        # Build NGL selection: centers (with altloc) + found residues not
        # already covered by centers.
        selections = []
        center_residues = set()
        if hasattr(site, 'centers') and site.centers:
            for center in site.centers:
                selection = f":{center.chain} and {center.resid}"
                if center.altloc and center.altloc.strip():
                    selection += f" and %{center.altloc.strip()}"
                selections.append(selection)
                center_residues.add((center.chain, center.resid))
        for chain, resname, resid, icode in residue_list:
            if (chain, resid) not in center_residues:
                selections.append(f":{chain} and {resid}")

        ngl_selection = " or ".join(selections) if selections else ""
        return _show_residues_in_viewer(
            getattr(self, 'processor', None),
            self.console,
            ngl_selection,
            label='found_residues_and_centers',
        )

# ===== MAIN WORKFLOW ORCHESTRATOR =====

class ComprehensiveRedoxDetector:
    """Main orchestrator for the complete redox site detection workflow"""
    
    def __init__(self, console: Console = None, processor=None):
        """
        Initialize the detector.

        Args:
            console: Rich console for output. If None, creates a new one.
            processor: Optional processor instance for session recording
        """
        self.console = console if console else Console()
        self.processor = processor
        self.config = DetectionConfig()
        self.config_interface = DetectionConfigInterface(console=self.console, processor=processor)
        self.all_centers = []
        self.selected_centers = []
        self.site_groups = []
        self.final_sites = []
        
        # Template-based automation
        self.use_templates = False
        self.site_types = {}  # site_idx -> site_type_name
        self.templates = {}  # site_type_name -> SiteTemplate
        self.template_results = {}  # site_type_name -> List[site_results]
        self.transformer_mappings = {}  # site_type_name -> transformer_name
        
        # Set up logging
        self.logger = logging.getLogger(__name__)
    
    def detect_redox_sites(self, structure_file: str = None, structure: Structure = None, selected_chains: List[str] = None, interactive: bool = True) -> List[RedoxSite]:
        """
        Complete redox site detection workflow
        
        Phase 1: Configuration & Setup
        Phase 2: Inventory with relational information  
        Phase 3: Center selection & grouping
        Phase 4: Site refinement (simplified for this implementation)
        """
        
        if not HAS_BIOPYTHON:
            self.console.print("[red]Error: BioPython not available. Please install BioPython to run structure analysis.[/red]")
            return []
        
        # Load or use provided structure
        if structure is not None:
            # Use provided structure object
            self.console.print(f"\n[bold cyan]Analyzing provided structure[/bold cyan]")
            structure_name = getattr(structure, 'id', 'provided_structure')
        elif structure_file:
            # Load structure from file
            try:
                with Progress(
                    SpinnerColumn(),
                    TextColumn("[progress.description]{task.description}"),
                    console=self.console
                ) as progress:
                    task = progress.add_task("Loading structure...", total=None)
                    parser = PDBParser(QUIET=True)
                    structure = parser.get_structure("structure", structure_file)
                    progress.update(task, description="Structure loaded successfully")
                    
                # Display structure info
                self.console.print(f"\n[bold cyan]Analyzing structure:[/bold cyan] {structure_file}")
                structure_name = structure_file
                    
            except Exception as e:
                self.console.print(f"[red]Error loading structure: {e}[/red]")
                self.logger.error(f"Structure loading failed: {e}", exc_info=True)
                return []
        else:
            self.console.print("[red]Error: No structure provided. Must supply either structure_file or structure object.[/red]")
            return []
            
        # Display chain selection info
        if selected_chains:
            self.console.print(f"[bold]Selected chains:[/bold] {', '.join(selected_chains)}")
        else:
            self.console.print("[bold]Selected chains:[/bold] all")
        
        # Phase 1: Configuration
        if interactive:
            self.console.print("\n[bold underline]Configuration Setup[/bold underline]")
            self.config = self.config_interface.configure_interactively()
        else:
            config_table = Table(title="Default Configuration")
            config_table.add_column("Parameter", style="cyan")
            config_table.add_column("Value", style="green")
            config_table.add_row("Inventory search radius", f"{self.config.bond_search_distance}Å")
            self.console.print(config_table)
        
        # Phase 2: Inventory Building
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=self.console
        ) as progress:
            task = progress.add_task("Scanning structure for redox centers...", total=None)
            self.all_centers = self._build_comprehensive_inventory(structure, selected_chains)
            progress.update(task, description=f"Found {len(self.all_centers)} potential redox centers")
        
        if not self.all_centers:
            self.console.print("[yellow]No redox centers found in structure.[/yellow]")
            return []
        
        # Build relational information
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=self.console
        ) as progress:
            task = progress.add_task("Analyzing structural relationships...", total=None)
            
            relational_analyzer = RelationalAnalyzer(self.config)
            inventory_display = InventoryDisplay(self.config, console=self.console)
            
            nearby_bonds = relational_analyzer.find_nearby_bonds(structure, self.all_centers, selected_chains)
            metal_clusters = relational_analyzer.find_metal_clusters(self.all_centers)
            
            progress.update(task, description="Structural analysis complete")
        
        # Display comprehensive inventory
        self.console.print("\n[bold underline]Redox Center Inventory[/bold underline]")
        inventory_display.display_inventory(self.all_centers, nearby_bonds, metal_clusters)

        # Auto-highlight every detected center with its 1-based ID color so
        # the user can read the inventory table and the 3D viewer side by
        # side. For each center we add TWO labeled reps sharing the
        # palette color: the center atoms (ball+stick) and the nearby
        # residues from the table's "Nearby Atom" column (ball+stick).
        # Same color = same row in the table; separate labels mean each
        # can be toggled in NGL's Representation Manager.
        inventory_items = []
        for i, c in enumerate(self.all_centers, start=1):
            inventory_items.append({
                "selection": _center_to_ngl(c),
                "label": f"inventory_{i}_center",
                "color_index": i,
            })
            # Collect unique residues from this center's nearby_bonds entry.
            center_key = (c.chain, c.resid, c.insertion_code)
            nearby_residues = sorted({
                (atom['chain'], atom['resid'])
                for atom in nearby_bonds.get(center_key, [])
            })
            if nearby_residues:
                nearby_sel = " or ".join(
                    f":{ch} and {rid}" for ch, rid in nearby_residues
                )
                inventory_items.append({
                    "selection": nearby_sel,
                    "label": f"inventory_{i}_nearby",
                    "color_index": i,
                })
        _auto_show_palette_highlights(self.processor, inventory_items)

        # Manual disulfide specification prompt (interactive mode only)
        if interactive:
            disulfide_count = len([c for c in self.all_centers
                                  if c.properties.get('is_disulfide_bonded', False)]) // 2

            if disulfide_count == 0:
                # No auto-detected disulfides — go straight into the manual
                # helper. It shows the full SG–SG distance ranking and accepts
                # 'none' if the user just wants to skip after seeing it, so a
                # separate y/n gate would only add a step.
                self.console.print(
                    f"\n[yellow]⚠ No SG–SG pairs were within the "
                    f"{self.config.disulfide_distance_threshold:.2f} Å auto-detection cutoff.[/yellow]"
                )
                manually_added = self._prompt_manual_disulfide_specification(structure, selected_chains)
                if manually_added:
                    self.all_centers.extend(manually_added)
                    self.console.print(f"[green]✓ Added {len(manually_added)} CYS center(s) from {len(manually_added)//2} manual disulfide bond(s)[/green]")
            else:
                # Some disulfides found - offer to add more
                if confirm_with_context(
                    self.processor,
                    f"\nFound {disulfide_count} disulfide bond(s). Specify additional disulfides?",
                    default=False,  # Default No when already found some
                    module="Redox Site Detector",
                    description="Additional disulfide specification"
                ):
                    manually_added = self._prompt_manual_disulfide_specification(structure, selected_chains)
                    if manually_added:
                        self.all_centers.extend(manually_added)
                        self.console.print(f"[green]✓ Added {len(manually_added)} CYS center(s) from {len(manually_added)//2} additional disulfide bond(s)[/green]")

        if not interactive:
            # In non-interactive mode, each center becomes its own site
            self.console.print("\n[bold cyan]Non-interactive mode:[/bold cyan] Creating individual sites for each center...")
            self.selected_centers = self.all_centers.copy()
            
            # Each center becomes its own site - no automatic grouping
            self.site_groups = [[center] for center in self.selected_centers]
            
        else:
            # Phase 3: Interactive Center Selection & Grouping
            self.console.print("\n[bold underline]Center Selection & Grouping[/bold underline]")
            selection_interface = CenterSelectionInterface(self.config, console=self.console, processor=self.processor)
            # Hand nearby_bonds to the selection interface so its post-
            # selection viewer hook can keep nearby-residue context for
            # the centers the user kept (matches the inventory hook's
            # center+nearby color grouping).
            selection_interface.inventory_nearby_bonds = nearby_bonds
            self.selected_centers, confirmed_disulfide_pairs = selection_interface.select_centers(self.all_centers)

            if not self.selected_centers:
                self.console.print("[yellow]No centers selected.[/yellow]")
                return []

            grouping_interface = CenterGroupingInterface(self.config, console=self.console, processor=self.processor)
            self.site_groups = grouping_interface.group_centers_into_sites(
                self.selected_centers,
                confirmed_disulfide_pairs=confirmed_disulfide_pairs
            )

            # Re-colour the viewer by site assignment: every center within
            # a site shares one palette colour, distinct from the next
            # site's. Lets the user visually verify the grouping is what
            # they intended.
            _auto_show_palette_highlights(
                self.processor,
                [
                    {
                        "selection": " or ".join(_center_to_ngl(c) for c in group),
                        "label": f"site_{i}",
                    }
                    for i, group in enumerate(self.site_groups, start=1)
                ],
            )

        # Auto-assign disulfide site types and transformer mappings
        self._auto_assign_disulfide_site_types()
        self._auto_assign_zinc_cys4_site_types()

        # Site Type Categorization for Template-Based Automation
        if interactive and len(self.site_groups) > 1:
            self._configure_site_types()
        
        # Phase 4: Complete Site Refinement with ALL 9 search method combinations
        self.console.print("\n[bold underline]Site Refinement[/bold underline]")
        
        if self.use_templates:
            self._process_sites_with_templates(structure, interactive)
        else:
            self._process_sites_standard(structure, interactive)
        
        # Show final results
        if self.final_sites:
            self.console.print(f"\n[bold green]✓ Detection complete![/bold green] Created {len(self.final_sites)} fully refined redox sites.")

            # Offer site review and manual refinement options
            self._offer_site_review_options(structure)
        else:
            self.console.print("\n[yellow]Detection complete - no redox sites identified.[/yellow]")

        return self.final_sites
    
    def _map_site_types_to_transformers(self, assignments: List[Tuple[List[int], str]]) -> Dict[str, str]:
        """Map site types to transformer names if they match"""
        mappings = {}

        registry = _get_transformer_registry()
        if not registry:
            return mappings

        available_transformers = registry.get_all_transformers()
        
        for site_indices, site_type in assignments:
            # Check if site_type matches any transformer name
            for transformer_name, transformer_class in available_transformers.items():
                if (site_type.lower() == transformer_class.TRANSFORMER_NAME.lower() or
                    site_type.lower() in [st.lower() for st in transformer_class.SUPPORTED_SITE_TYPES]):
                    mappings[site_type] = transformer_class.TRANSFORMER_NAME
                    break
        
        return mappings
    
    def _display_available_transformers(self):
        """Display available transformers and their bond requirements"""
        registry = _get_transformer_registry()
        if not registry:
            self.console.print("\n[yellow]Note: Transformer registry not available. Using custom site types.[/yellow]")
            return

        transformers = registry.get_all_transformers()
        if not transformers:
            self.console.print("\n[yellow]No transformers available.[/yellow]")
            return
            
        # Group transformers by chemistry category for easier visual scanning.
        # Categories render in this order; any transformer not in the map goes
        # under "Other".
        CATEGORY_ORDER = [
            "Hemes",
            "Iron-sulfur clusters",
            "Organic cofactors",
            "Metal sites",
            "Utility",
        ]
        CATEGORY_FOR_TRANSFORMER = {
            "heme_bis_his_c_type":   "Hemes",
            "heme_bis_his_b_type":   "Hemes",
            "heme_cys_axial_b_type": "Hemes",
            "heme_his_met_axial_c_type": "Hemes",
            "iron_sulfur_4fe4s":     "Iron-sulfur clusters",
            "flavin_fmn":            "Organic cofactors",
            "flavin_fad":            "Organic cofactors",
            "nicotinamide_nadp":     "Organic cofactors",
            "nicotinamide_nad":      "Organic cofactors",
            "pterin_biopterin":      "Organic cofactors",
            "zinc_cys4":             "Metal sites",
            "disulfide":             "Utility",
            "no_transformation":     "Utility",
        }

        # Category tags that show up in SUPPORTED_SITE_TYPES as a category
        # marker but aren't meant to be typed by the user as a template name.
        # Filtered out of the alias display.
        CATEGORY_TAG_ALIASES = {"organic_cofactor"}

        from collections import OrderedDict
        buckets = OrderedDict((cat, []) for cat in CATEGORY_ORDER)
        buckets["Other"] = []
        for transformer_class in transformers.values():
            name = transformer_class.TRANSFORMER_NAME
            cat = CATEGORY_FOR_TRANSFORMER.get(name, "Other")
            buckets[cat].append(transformer_class)

        # Use one borderless Rich table per category so descriptions wrap with
        # proper column alignment. (Plain f-string formatting breaks the
        # alignment when descriptions wrap past the terminal width.) A unified
        # name-column width keeps the description-column left edge aligned
        # across categories.
        all_classes = [tc for cls in buckets.values() for tc in cls]
        if not all_classes:
            return
        name_col_w = max(len(tc.TRANSFORMER_NAME) for tc in all_classes)

        self.console.print(
            "\n[bold cyan]Available transformers — use any name shown in bold:[/bold cyan]"
        )

        for cat in buckets:
            classes = buckets[cat]
            if not classes:
                continue
            self.console.print(f"\n[bold cyan]{cat}[/bold cyan]")
            table = Table(
                show_header=False, box=None, pad_edge=False,
                show_edge=False, padding=(0, 2, 0, 2),
            )
            table.add_column("Name", style="bold", no_wrap=True,
                             width=name_col_w, min_width=name_col_w)
            table.add_column("Description", overflow="fold")
            for tc in classes:
                # Build a description that includes any meaningful aliases
                # (skipping the category-tag noise) on a dim continuation line.
                desc = tc.DESCRIPTION
                if hasattr(tc, "SUPPORTED_SITE_TYPES") and tc.SUPPORTED_SITE_TYPES:
                    if tc.SUPPORTED_SITE_TYPES == ["*"]:
                        desc += "\n[grey50]· compatible with any site type[/grey50]"
                    else:
                        aliases = [
                            a for a in tc.SUPPORTED_SITE_TYPES
                            if a != tc.TRANSFORMER_NAME
                            and a not in CATEGORY_TAG_ALIASES
                        ]
                        if aliases:
                            desc += f"\n[grey50]· also: {', '.join(aliases)}[/grey50]"
                table.add_row(tc.TRANSFORMER_NAME, desc)
            self.console.print(table)
    
    
    def _auto_assign_disulfide_site_types(self):
        """
        Automatically assign site_type = "disulfide" to disulfide bond sites.

        This method identifies sites that are disulfide bonds (two CYS centers with
        is_disulfide_bonded property) and automatically assigns them the "disulfide"
        site type and transformer mapping.

        This happens before user-interactive site type configuration, so:
        - Disulfide sites are always auto-assigned
        - User can still override if using template mode
        - The disulfide transformer will be automatically selected during transformation
        """
        if not self.site_groups:
            return

        disulfide_count = 0

        for site_idx, site_centers in enumerate(self.site_groups, 1):
            # Check if this is a disulfide site:
            # 1. Exactly 2 centers
            # 2. Both are CYS residues
            # 3. Both have is_disulfide_bonded property set to True

            if len(site_centers) != 2:
                continue

            # Check if both centers are CYS with disulfide bond property
            is_disulfide_site = True
            for center in site_centers:
                # Check if center is CYS
                if not (hasattr(center, 'resname') and center.resname in ['CYS', 'CYX']):
                    is_disulfide_site = False
                    break

                # Check if center has is_disulfide_bonded property
                if not (hasattr(center, 'properties') and
                        center.properties.get('is_disulfide_bonded', False)):
                    is_disulfide_site = False
                    break

            if is_disulfide_site:
                # Auto-assign site_type = "disulfide"
                self.site_types[site_idx] = "disulfide"
                disulfide_count += 1

        # Add transformer mapping for disulfide sites (if any were found)
        if disulfide_count > 0:
            # Only add if not already present (don't overwrite user settings)
            if "disulfide" not in self.transformer_mappings:
                self.transformer_mappings["disulfide"] = "disulfide"

            # Log the auto-assignment
            self.logger.debug(f"Auto-assigned 'disulfide' site type to {disulfide_count} site(s)")
            self.console.print(f"[grey50]Auto-assigned 'disulfide' transformer to {disulfide_count} disulfide bond site(s)[/grey50]")

    def _auto_assign_zinc_cys4_site_types(self):
        """
        Automatically assign site_type = "zinc_cys4" to tetrahedral Zn(Cys)4 sites.

        A Zn(Cys)4 site is recognized as: 1 Zn METAL_ION center + 4 CYS
        REDOX_AMINO_ACID centers grouped into the same site by the spatial
        grouping pass. Auto-assignment routes these sites through the
        bundled zinc_cys4 transformer + Guberman_ZnCys4 forcefield set.
        """
        if not self.site_groups:
            return

        zinc_cys4_count = 0

        for site_idx, site_centers in enumerate(self.site_groups, 1):
            # Skip sites already auto-assigned (e.g., disulfide).
            if site_idx in self.site_types:
                continue

            if len(site_centers) != 5:
                continue

            metals = [c for c in site_centers
                      if c.center_type == CenterType.METAL_ION
                      and getattr(c, "element", "").upper() == "ZN"]
            cys_ligands = [c for c in site_centers
                           if c.center_type == CenterType.REDOX_AMINO_ACID
                           and c.resname == "CYS"
                           and c.atom_name == "SG"]

            if len(metals) == 1 and len(cys_ligands) == 4:
                self.site_types[site_idx] = "zinc_cys4"
                zinc_cys4_count += 1

        if zinc_cys4_count > 0:
            if "zinc_cys4" not in self.transformer_mappings:
                self.transformer_mappings["zinc_cys4"] = "zinc_cys4"

            self.logger.debug(f"Auto-assigned 'zinc_cys4' site type to {zinc_cys4_count} site(s)")
            self.console.print(f"[grey50]Auto-assigned 'zinc_cys4' transformer to {zinc_cys4_count} tetrahedral Zn(Cys)4 site(s)[/grey50]")

    def _configure_site_types(self):
        """Configure site types for template-based automation"""
        from proprep.utils.prompts import prompt_with_context, confirm_with_context
        
        self.console.print("\n[bold cyan]Template-Based Automation[/bold cyan]")
        self.console.print()
        self.console.print("[bold]Concepts:[/bold]")
        self.console.print("  • [cyan]Templates[/cyan]: User-specified tags for detection recipes")
        self.console.print("    Define what to search for (coordination, ligands, geometry)")
        self.console.print()
        self.console.print("  • [cyan]Transformers[/cyan]: Pre-coded instructions for forcefield compatibility")
        self.console.print("    Define what to change (residue names, atom names, etc.)")
        self.console.print()
        self.console.print("[bold]How they work together:[/bold]")
        self.console.print("  1. Assign distinct template names for each unique grouping")
        self.console.print("     (e.g., 'ca_site_type1', 'ca_site_type2' for Ca sites with different ligands)")
        self.console.print("  2. Templates automate refinement after configuring the first site of each type")
        self.console.print("  3. During Redox Site Preparation: If template name matches a transformer, it's automatically applied")
        self.console.print("  4. During Redox Site Preparation: Otherwise, you'll assign transformers manually")
        self.console.print("  5. During Redox Site Preparation: Use 'no_transformation' if sites are already forcefield-compliant")

        # Ask if user wants to use templates
        self.use_templates = confirm_with_context(
            processor=self.processor,
            prompt="\n[green]Use template-based automation for systematic redox sites?[/green]",
            default=False,
            module="Redox Detector - Template Configuration",
            description="Use template-based automation"
        )
        
        if not self.use_templates:
            return
            
        # Show available sites
        self.console.print(f"\n[bold]Available sites:[/bold] 1-{len(self.site_groups)}")
        
        # Display available transformers and their information
        self._display_available_transformers()
        
        self.console.print("\n[bold]Site type assignment examples:[/bold]")
        registry = _get_transformer_registry()
        if registry:
            # Get first available transformer for example
            transformers = registry.get_all_transformers()
            if transformers:
                first_transformer = list(transformers.values())[0]
                transformer_name = first_transformer.TRANSFORMER_NAME
                self.console.print(f'1-18:"{transformer_name}"')
                self.console.print(f'19,21:"{transformer_name}", 20:"custom_site_type"')
            else:
                self.console.print('1-18:"custom_site_type"')
                self.console.print('19,21:"custom_site_type", 20:"another_type"')
        else:
            self.console.print('1-18:"bis-his c-type heme"')
            self.console.print('19,21:"Fe4S4 cluster", 20:"his-met b-type heme"')
        
        # Get site type assignments
        while True:
            assignments_input = prompt_with_context(
                processor=self.processor,
                prompt="\n[green]Enter site type assignments[/green]",
                default="",
                module="Redox Detector - Template Configuration",
                description="Enter site type assignments",
                options_map={"custom": "Site assignments (e.g., 1-18:\"Fe4S4 cluster\", 19:\"heme\")"}
            ).strip()
            
            if not assignments_input:
                self.console.print("[yellow]No assignments entered. Template automation disabled.[/yellow]")
                self.use_templates = False
                return
                
            try:
                # Parse assignments like: 1-18:"bis-his c-type heme", 19,21:"Fe4S4 cluster"
                assignments = self._parse_site_type_assignments(assignments_input)
                
                # Validate all sites are assigned
                assigned_sites = set()
                for site_indices, site_type in assignments:
                    assigned_sites.update(site_indices)
                
                all_sites = set(range(1, len(self.site_groups) + 1))
                if assigned_sites != all_sites:
                    missing = all_sites - assigned_sites
                    extra = assigned_sites - all_sites
                    if missing:
                        self.console.print(f"[red]Missing site assignments: {sorted(missing)}[/red]")
                    if extra:
                        self.console.print(f"[red]Invalid site numbers: {sorted(extra)}[/red]")
                    continue
                    
                # Store assignments
                for site_indices, site_type in assignments:
                    for site_idx in site_indices:
                        self.site_types[site_idx] = site_type
                
                # Map site types to transformers if available
                self.transformer_mappings = self._map_site_types_to_transformers(assignments)
                        
                # Show confirmation
                type_counts = {}
                for site_type in self.site_types.values():
                    type_counts[site_type] = type_counts.get(site_type, 0) + 1
                    
                self.console.print("\n[bold green]Site type assignments confirmed:[/bold green]")
                for site_type, count in type_counts.items():
                    self.console.print(f"  {site_type}: {count} sites")
                
                # Display transformer mappings if any were found
                if self.transformer_mappings:
                    self.console.print("\n[bold cyan]Transformer Mappings:[/bold cyan]")
                    for site_type, transformer_name in self.transformer_mappings.items():
                        self.console.print(f"  {site_type} → {transformer_name}")
                    self.console.print("[grey50]These sites will be pre-configured for the appropriate transformer.[/grey50]")

                break
                
            except ValueError as e:
                self.console.print(f"[red]Error parsing assignments: {e}[/red]")
                self.console.print("Please use format: site_range:\"site_type\", ...")
                continue
    
    def _parse_site_type_assignments(self, assignments_input: str) -> List[Tuple[List[int], str]]:
        """Parse site type assignments like: 1-18:\"bis-his c-type heme\", 19,21:\"Fe4S4 cluster\""""
        assignments = []

        # Split by commas that separate assignments (after closing quote)
        # This is more complex because commas can appear in:
        # 1. Site ranges (e.g., "19,21:...")
        # 2. Between assignments (e.g., ...site1", site2:...)
        parts = []
        current = ""
        in_quotes = False
        last_quote_pos = -1

        for i, char in enumerate(assignments_input):
            if char == '"' and (not current or current[-1] != '\\'):
                in_quotes = not in_quotes
                current += char
                if not in_quotes:  # Closing quote
                    last_quote_pos = len(current) - 1
            elif char == ',' and not in_quotes:
                # This comma is outside quotes
                # Check if we've seen a closing quote and a colon before this comma
                # If yes, this is an assignment separator
                # If no, this is part of a site range
                if ':' in current and last_quote_pos > current.rfind(':'):
                    # We have a colon and the last quote came after it
                    # This means we completed an assignment, so split here
                    if current.strip():
                        parts.append(current.strip())
                    current = ""
                    last_quote_pos = -1
                else:
                    # This comma is part of the site range (before the colon)
                    current += char
            else:
                current += char

        if current.strip():
            parts.append(current.strip())

        for part in parts:
            if ':' not in part:
                raise ValueError(f"Missing ':' in assignment: {part}")

            site_range, site_type = part.split(':', 1)
            site_range = site_range.strip()
            site_type = site_type.strip().strip('"')

            if not site_type:
                raise ValueError(f"Empty site type in: {part}")

            # Parse site range (e.g., "1-18" or "19,21" or "1,3,5")
            site_indices = []
            for range_part in site_range.split(','):
                range_part = range_part.strip()
                if '-' in range_part:
                    start, end = range_part.split('-', 1)
                    site_indices.extend(range(int(start), int(end) + 1))
                else:
                    site_indices.append(int(range_part))

            assignments.append((site_indices, site_type))

        return assignments
    
    def _process_sites_standard(self, structure, interactive):
        """Standard site processing without templates"""
        for i, site_centers in enumerate(self.site_groups, 1):
            site_id = f"site_{i}"
            site = RedoxSite(site_id, structure.id or "structure")
            
            # Add centers to site
            for center in site_centers:
                site.add_center(center)
            
            # Add atoms from all residues in the site (initial)
            self._add_atoms_to_site(site, site_centers, structure)
            
            # No automatic bond detection - bonds are user-defined
            site.bonds.clear()
            
            site.detection_method = "comprehensive_workflow_with_refinement"
            site.detection_parameters = {
                'bond_search_distance': self.config.bond_search_distance
            }
            
            # Log site creation at debug level
            self.logger.debug(f"Created initial {site_id} with {len(site.centers)} center(s), {len(site.atoms)} atoms, {len(site.bonds)} bonds")
            
            # FULL INTERACTIVE REFINEMENT - implement ALL 9 search method combinations
            if interactive:
                refinement_interface = SiteRefinementInterface(self.config, console=self.console, transformer_mappings=getattr(self, 'transformer_mappings', {}), processor=self.processor)
                refined_site = refinement_interface.refine_site_interactively(site, structure)
                self.final_sites.append(refined_site)
            else:
                # Non-interactive mode - perform basic expansion with default parameters
                refined_site = self._auto_refine_site(site, structure)
                self.final_sites.append(refined_site)
    
    def _process_sites_with_templates(self, structure, interactive):
        """Template-based site processing"""
        # Group sites by type and process each type
        sites_by_type = {}
        for i, site_centers in enumerate(self.site_groups, 1):
            site_type = self.site_types.get(i, "unknown")
            if site_type not in sites_by_type:
                sites_by_type[site_type] = []
            sites_by_type[site_type].append((i, site_centers))
        
        # Process each site type
        for site_type, site_list in sites_by_type.items():
            self.console.print(f"\n[bold cyan]Processing {len(site_list)} sites of type: {site_type}[/bold cyan]")
            
            # Process first site of this type to create template
            first_site_idx, first_site_centers = site_list[0]
            first_site = self._create_initial_site(first_site_idx, first_site_centers, structure)
            
            self.console.print(f"\n[bold]Configuring template using site_{first_site_idx}[/bold]")
            refinement_interface = SiteRefinementInterface(self.config, console=self.console, template_mode=True, transformer_mappings=getattr(self, 'transformer_mappings', {}), processor=self.processor)
            refined_first_site, template = refinement_interface.refine_site_for_template(first_site, structure, site_type)
            
            if template:
                self.templates[site_type] = template
                self.final_sites.append(refined_first_site)
                
                # Ask user if they want to apply template to remaining sites
                if len(site_list) > 1:
                    apply_template = confirm_with_context(
                        processor=self.processor,
                        prompt=f"\n[green]Apply these settings to all other {site_type} sites ({len(site_list)-1} remaining)?[/green]",
                        default=True,
                        module="Redox Detector - Template Application",
                        description=f"Apply template to {len(site_list)-1} {site_type} sites"
                    )
                    
                    if apply_template:
                        # Apply template to remaining sites
                        remaining_sites = site_list[1:]
                        template_results = self._apply_template_to_sites(template, remaining_sites, structure)
                        
                        # Show results and handle failures
                        self._review_template_results(site_type, template_results, structure)
                    else:
                        # Process remaining sites manually
                        for site_idx, site_centers in site_list[1:]:
                            site = self._create_initial_site(site_idx, site_centers, structure)
                            refined_site = refinement_interface.refine_site_interactively(site, structure)
                            self.final_sites.append(refined_site)
                else:
                    self.console.print(f"[green]✓ Template created for {site_type} (only 1 site of this type)[/green]")
            else:
                # Template creation failed, process manually
                self.console.print(f"[yellow]Template creation failed for {site_type}. Processing remaining sites manually.[/yellow]")
                for site_idx, site_centers in site_list[1:]:
                    site = self._create_initial_site(site_idx, site_centers, structure)
                    refined_site = refinement_interface.refine_site_interactively(site, structure)
                    self.final_sites.append(refined_site)
    
    def _create_initial_site(self, site_idx, site_centers, structure):
        """Create initial site with centers and atoms"""
        site_id = f"site_{site_idx}"
        site = RedoxSite(site_id, structure.id or "structure")
        
        # Set site type from user categorization
        site.site_type = self.site_types.get(site_idx, "unknown")
        
        # Add centers to site
        for center in site_centers:
            site.add_center(center)
        
        # Add atoms from all residues in the site (initial)
        self._add_atoms_to_site(site, site_centers, structure)
        
        # No automatic bond detection - bonds are user-defined
        site.bonds.clear()
        
        site.detection_method = "comprehensive_workflow_with_refinement"
        site.detection_parameters = {
            'bond_search_distance': self.config.bond_search_distance
        }
        
        # Log site creation (debug level)
        self.logger.debug(f"Created initial {site_id} with {len(site.centers)} center(s), {len(site.atoms)} atoms, {len(site.bonds)} bonds")
        
        return site
    
    def _apply_template_to_sites(self, template: SiteTemplate, sites_list, structure):
        """Apply template to multiple sites and return results"""
        results = []
        
        for site_idx, site_centers in sites_list:
            site = self._create_initial_site(site_idx, site_centers, structure)
            
            try:
                # Apply template using automated refinement
                refinement_interface = SiteRefinementInterface(self.config, console=self.console, template_mode=True, transformer_mappings=getattr(self, 'transformer_mappings', {}), processor=self.processor)
                refined_site = refinement_interface.apply_template_to_site(site, structure, template)
                
                results.append({
                    'site_idx': site_idx,
                    'site': refined_site,
                    'status': 'success',
                    'error': None
                })
                
            except Exception as e:
                results.append({
                    'site_idx': site_idx,
                    'site': site,
                    'status': 'failed',
                    'error': str(e)
                })
        
        return results
    
    def _review_template_results(self, site_type, results, structure):
        """Review template application results and handle failures"""
        from rich.table import Table
        from proprep.utils.prompts import prompt_with_context
        
        # Show results summary including all sites of this type
        results_table = Table(title=f"Site Results for \"{site_type}\"")
        results_table.add_column("Site", style="cyan", width=8)
        results_table.add_column("Status", style="green", width=8)
        results_table.add_column("Details", style="white", width=20)
        results_table.add_column("Warnings", style="yellow", width=15)

        successful_sites = []
        failed_sites = []

        # First, add any existing sites of this type that are already in final_sites
        existing_sites_of_type = [site for site in self.final_sites if hasattr(site, 'site_type') and site.site_type == site_type]
        for site in existing_sites_of_type:
            # Extract site number from site_id (e.g., "site_1" -> 1)
            site_num = site.site_id.split('_')[-1] if hasattr(site, 'site_id') else "?"

            # Check for long bonds (>4Å)
            long_bonds = [bond for bond in site.bonds if bond.distance > 4.0]
            warnings = f"{len(long_bonds)} bond(s) >4Å" if long_bonds else "—"

            results_table.add_row(
                f"Site {site_num}",
                "✓",
                f"{len(site.atoms)} atoms, {len(site.bonds)} bonds",
                warnings
            )

        # Then add the new template results
        for result in results:
            site_idx = result['site_idx']
            status = result['status']

            if status == 'success':
                site = result['site']
                successful_sites.append(site)

                # Check for long bonds (>4Å)
                long_bonds = [bond for bond in site.bonds if bond.distance > 4.0]
                warnings = f"{len(long_bonds)} bond(s) >4Å" if long_bonds else "—"

                results_table.add_row(
                    f"Site {site_idx}",
                    "✓",
                    f"{len(site.atoms)} atoms, {len(site.bonds)} bonds",
                    warnings
                )
            else:
                failed_sites.append(result)
                results_table.add_row(
                    f"Site {site_idx}",
                    "✗ FAILED",
                    result['error'][:50] + "..." if len(result['error']) > 50 else result['error'],
                    "—"
                )
        
        self.console.print(results_table)
        
        # Add successful sites to final results
        self.final_sites.extend(successful_sites)
        
        # Handle failed sites
        if failed_sites:
            failed_indices = [str(result['site_idx']) for result in failed_sites]
            rerun_input = prompt_with_context(
                processor=self.processor,
                prompt=f"\n[yellow]Re-run which sites manually (comma-separated list)[/yellow]",
                default="",
                module="Redox Detector - Template Application",
                description="Select failed sites to re-run manually",
                options_map={"custom": f"Site numbers (e.g., {', '.join(failed_indices[:3])})"}
            ).strip()
            
            if rerun_input:
                try:
                    rerun_indices = [int(idx.strip()) for idx in rerun_input.split(',') if idx.strip()]
                    
                    for result in failed_sites:
                        if result['site_idx'] in rerun_indices:
                            self.console.print(f"\n[bold]Manually refining site_{result['site_idx']}[/bold]")
                            refinement_interface = SiteRefinementInterface(self.config, console=self.console, transformer_mappings=getattr(self, 'transformer_mappings', {}), processor=self.processor)
                            refined_site = refinement_interface.refine_site_interactively(result['site'], structure)
                            self.final_sites.append(refined_site)
                            
                except ValueError:
                    self.console.print("[red]Invalid site numbers entered. Skipping manual refinement.[/red]")
        
        self.console.print(f"[green]✓ Completed processing {site_type}: {len(successful_sites)} sites automatically processed using template, {len(failed_sites)} failed[/green]")
    
    def _offer_site_review_options(self, structure):
        """Offer options to review sites and manually refine selected ones"""
        from proprep.utils.prompts import prompt_with_context, confirm_with_context
        from rich.table import Table
        
        while True:
            self.console.print("\n[bold]Site Review Options:[/bold]")
            self.console.print("[grey50]1.[/grey50] View site summary for specific sites")
            self.console.print("[grey50]2.[/grey50] Edit detected sites")
            self.console.print("[grey50]3.[/grey50] Manually re-run detection on specific sites")
            self.console.print("[grey50]4.[/grey50] Export sites")
            self.console.print("[grey50]5.[/grey50] Continue")

            choice = prompt_with_context(
                processor=self.processor,
                prompt="[green]Choose option[/green]",
                choices=["1", "2", "3", "4", "5"],
                default="5",
                module="Redox Detector - Site Review",
                description="Select site review action",
                options_map={
                    "1": "View site summary for specific sites",
                    "2": "Edit detected sites",
                    "3": "Manually re-run detection on specific sites",
                    "4": "Export sites",
                    "5": "Continue"
                }
            ).strip()

            if choice == "5":
                # Module exit: drop the last-processed-site halos and
                # lay down a clean overview of every detected site so
                # the viewer reflects the full result set rather than
                # the in-flight refinement state.
                _show_all_sites_overview(
                    getattr(self, 'processor', None), self.final_sites,
                )
                return
            elif choice == "1":
                self._show_site_summaries()
            elif choice == "2":
                self._edit_detected_sites()
            elif choice == "3":
                self._manually_refine_selected_sites(structure)
            elif choice == "4":
                # Export and then return to menu
                self._offer_export_options()

    def _show_site_summaries(self):
        """Show detailed summaries for selected sites"""
        from proprep.utils.prompts import prompt_with_context
        
        # Show available sites table
        sites_table = Table(title="Available Sites")
        sites_table.add_column("#", style="cyan", width=6)
        sites_table.add_column("Site ID", style="green")
        sites_table.add_column("Atoms", style="yellow")
        sites_table.add_column("Bonds", style="blue")
        sites_table.add_column("Centers", style="magenta")
        
        for i, site in enumerate(self.final_sites, 1):
            sites_table.add_row(
                f"[{i}]",
                site.site_id,
                str(len(site.atoms)),
                str(len(site.bonds)),
                str(len(site.centers))
            )
        
        self.console.print(sites_table)

        selection = prompt_with_context(
            processor=self.processor,
            prompt=f"[green]Select sites to view[/green] (1-{len(self.final_sites)}, comma-separated, 'all', or 'back')",
            module="Redox Detector - Site Review",
            description="Select sites to view details",
            options_map={"back": "Return to menu", "all": "View all sites", "custom": "Site numbers (e.g., 1-3, 5, 7)"}
        ).strip()

        if selection.lower() == 'back':
            return

        if selection.lower() == 'all':
            indices = list(range(len(self.final_sites)))
        else:
            try:
                # Use existing selection parser that handles ranges and comma-separated values
                indices = self._parse_selection_input(selection, len(self.final_sites))
                if indices is None:
                    self.console.print(f"[red]Invalid selection. Valid range is 1-{len(self.final_sites)}. Use numbers, ranges (1-3), or comma-separated (1,2,4).[/red]")
                    return
            except ValueError:
                self.console.print("[red]Invalid selection format[/red]")
                return

        # Show detailed summaries for selected sites
        for site_idx in indices:
            site = self.final_sites[site_idx]
            self.console.print(f"\n[bold underline]Site Summary: {site.site_id}[/bold underline]")

            # Import the method from SiteRefinementInterface
            refinement_interface = SiteRefinementInterface(self.config, console=self.console, transformer_mappings=getattr(self, 'transformer_mappings', {}), processor=self.processor)
            refinement_interface._display_site_summary(site)

    def _manually_refine_selected_sites(self, structure):
        """Allow manual refinement of selected sites"""
        from proprep.utils.prompts import prompt_with_context, confirm_with_context
        from rich.table import Table
        
        # Show available sites table
        sites_table = Table(title="Available Sites for Manual Refinement")
        sites_table.add_column("#", style="cyan", width=6)
        sites_table.add_column("Site ID", style="green")
        sites_table.add_column("Atoms", style="yellow")
        sites_table.add_column("Bonds", style="blue")
        sites_table.add_column("Centers", style="magenta")
        
        for i, site in enumerate(self.final_sites, 1):
            sites_table.add_row(
                f"[{i}]",
                site.site_id,
                str(len(site.atoms)),
                str(len(site.bonds)),
                str(len(site.centers))
            )
        
        self.console.print(sites_table)

        selection = prompt_with_context(
            processor=self.processor,
            prompt=f"[green]Select sites to manually refine[/green] (1-{len(self.final_sites)}, comma-separated, 'all', or 'back')",
            module="Redox Detector - Site Review",
            description="Select sites to manually refine",
            options_map={"back": "Return to menu", "all": "Refine all sites", "custom": "Site numbers (e.g., 1-3, 5, 7)"}
        ).strip()

        if selection.lower() == 'back':
            return

        if selection.lower() == 'all':
            indices = list(range(len(self.final_sites)))
        else:
            try:
                # Use existing selection parser that handles ranges and comma-separated values
                indices = self._parse_selection_input(selection, len(self.final_sites))
                if indices is None:
                    self.console.print(f"[red]Invalid selection. Valid range is 1-{len(self.final_sites)}. Use numbers, ranges (1-3), or comma-separated (1,2,4).[/red]")
                    return
            except ValueError:
                self.console.print("[red]Invalid selection format[/red]")
                return

        if confirm_with_context(
            processor=self.processor,
            prompt=f"[yellow]This will replace {len(indices)} site(s) with manually refined versions. Continue?[/yellow]",
            default=True,
            module="Redox Detector - Site Review",
            description=f"Confirm manual refinement of {len(indices)} sites"
        ):
            # Refine selected sites
            for site_idx in indices:
                site = self.final_sites[site_idx]
                self.console.print(f"\n[bold]Manually refining {site.site_id}[/bold]")

                # Create a fresh site with centers and their residue atoms (no bonds)
                fresh_site = RedoxSite(site.site_id, site.structure_id)
                fresh_site.site_type = site.site_type  # Preserve site type
                for center in site.centers:
                    fresh_site.add_center(center)

                # Add atoms from center residues to the fresh site
                self._add_atoms_to_site(fresh_site, list(site.centers), structure)

                refinement_interface = SiteRefinementInterface(self.config, console=self.console, transformer_mappings=getattr(self, 'transformer_mappings', {}), processor=self.processor)
                refined_site = refinement_interface.refine_site_interactively(fresh_site, structure)
                self.final_sites[site_idx] = refined_site

            self.console.print(f"[green]✓ Manual refinement complete for {len(indices)} site(s)[/green]")

    def _edit_detected_sites(self):
        """Launch the comprehensive RedoxSite editor"""
        from .redox_site_editor import RedoxSiteEditor

        if not self.final_sites:
            self.console.print("[yellow]No redox sites have been detected yet[/yellow]")
            return

        workspace = self.processor.workspace if self.processor else {}

        # Launch editor with workspace for structure loading
        editor = RedoxSiteEditor(
            sites=self.final_sites,
            console=self.console,
            processor=self.processor,
            workspace=workspace
        )

        modified_sites = editor.run()

        if modified_sites is not None:
            # User saved changes
            self.final_sites = modified_sites
            # Also update workspace if it exists
            if workspace:
                workspace.set("detected_redox_sites", modified_sites)
            self.console.print("\n[green]✓ Changes saved[/green]")
        else:
            # User cancelled
            self.console.print("\n[grey50]No changes made[/grey50]")

    def _offer_export_options(self):
        """Offer export options for detected redox sites"""
        from proprep.utils.prompts import prompt_with_context, confirm_with_context

        if not self.final_sites:
            return

        # Show export format options
        self.console.print("\n[bold]Export format options:[/bold]")
        self.console.print("[grey50]1.[/grey50] JSON (complete site data)")
        self.console.print("[grey50]2.[/grey50] PDB - all sites in one file")
        self.console.print("[grey50]3.[/grey50] PDB - each site in separate files")

        format_choice = prompt_with_context(
            processor=self.processor,
            prompt="[green]Choose export format[/green] (comma-separated list)",
            default="1",
            module="Redox Detector - Export",
            description="Choose export format",
            options_map={
                "1": "JSON (complete site data)",
                "2": "PDB - all sites in one file",
                "3": "PDB - each site in separate files"
            }
        ).strip()
        
        # Parse comma-separated choices
        choices = [choice.strip() for choice in format_choice.split(',') if choice.strip()]
        
        # Validate choices
        valid_choices = {'1', '2', '3'}
        invalid_choices = [c for c in choices if c not in valid_choices]
        if invalid_choices:
            self.console.print(f"[red]Invalid choices: {', '.join(invalid_choices)}. Valid options are 1, 2, 3[/red]")
            return
            
        if not choices:
            self.console.print("[yellow]No export formats selected[/yellow]")
            return
        
        # Get the source PDB filename for export functions
        source_pdb = getattr(self, 'source_pdb_file', 'structure.pdb')
        
        try:
            if "1" in choices:
                _export_to_json(
                    self.final_sites,
                    source_pdb,
                    self.console,
                    transformer_mappings=getattr(self, 'transformer_mappings', {}),
                )
                self.console.print("[green]✓ JSON export complete[/green]")
                
            if "2" in choices:
                _export_to_pdb_single(self.final_sites, source_pdb, self.console)
                self.console.print("[green]✓ Single PDB export complete[/green]")
                
            if "3" in choices:
                _export_to_pdb_separate(self.final_sites, source_pdb, self.console)
                self.console.print("[green]✓ Separate PDB export complete[/green]")
                
        except Exception as e:
            self.console.print(f"[red]Export failed: {str(e)}[/red]")
    
    def _auto_refine_site(self, site: RedoxSite, structure: Structure) -> RedoxSite:
        """Automatic site refinement for non-interactive mode"""
        self.logger.info(f"Auto-refining {site.site_id}...")
        
        # Create simple refinement interface for automatic expansion
        refinement_interface = SiteRefinementInterface(self.config, console=self.console, transformer_mappings=getattr(self, 'transformer_mappings', {}), processor=self.processor)
        
        # Perform one round of fixed distance search from geometric center
        search_params = SearchParameters(
            constraint=SearchConstraint.DISTANCE_CUTOFF,
            distance_method=DistanceMethod.FIXED,
            boundary_definition=BoundaryDefinition.GEOMETRIC_CENTER,
            radius=self.config.bond_search_distance
        )
        
        search_result = refinement_interface._perform_distance_search(site, structure, search_params)
        
        # Auto-select all coordination bonds for addition
        if search_result.detected_atoms:
            coordination_atoms = [
                atom for atom in search_result.detected_atoms 
                if atom['chemical_type'] in ['coordinate', 'covalent']
            ]
            
            if coordination_atoms:
                print(f"Auto-adding {len(coordination_atoms)} coordination atoms to {site.site_id}")
                
                updated_site = copy.deepcopy(site)
                for atom_info in coordination_atoms:
                    site_atom = RedoxSiteAtom(
                        chain=atom_info['chain'],
                        resname=atom_info['resname'],
                        resid=atom_info['resid'],
                        atom_name=atom_info['atom_name'],
                        coords=atom_info['coords'],
                        element=atom_info['element'],
                        insertion_code=atom_info['insertion_code']
                    )
                    updated_site.add_atom(site_atom)
                
                # Clear bonds - no automatic bond detection in auto-mode
                updated_site.bonds.clear()
                
                # Update search history
                updated_site.search_history.append({
                    'search_parameters': search_params,
                    'atoms_found': search_result.total_atoms_found,
                    'atoms_added': len(coordination_atoms),
                    'search_radius': search_params.radius,
                    'auto_refinement': True
                })
                
                return updated_site
        
        return site
    
    def _build_comprehensive_inventory(self, structure: Structure, selected_chains: List[str] = None) -> List[RedoxCenter]:
        """Build comprehensive inventory of all redox centers"""
        all_centers = []

        # Scan for non-standard residues FIRST (highest priority)
        # This ensures metal-containing cofactors (like hemes) are classified as
        # ORGANOMETALLIC_COFACTOR rather than having their metal atoms extracted as METAL_ION
        organic_scanner = NonStandardResidueScanner(self.config)
        organic_centers = organic_scanner.scan_structure(structure, all_centers, selected_chains)
        all_centers.extend(organic_centers)

        # Scan for metals SECOND, skipping residues already detected as cofactors
        metal_scanner = MetalIonScanner(self.config)
        metal_centers = metal_scanner.scan_structure(structure, existing_centers=all_centers,
                                                     selected_chains=selected_chains)
        all_centers.extend(metal_centers)

        # Scan for redox-active amino acids THIRD (individual CYS, TYR, etc.)
        amino_scanner = RedoxAminoAcidScanner(self.config)
        amino_centers = amino_scanner.scan_structure(structure, all_centers, selected_chains)
        all_centers.extend(amino_centers)

        # Scan for disulfide bonds FOURTH
        # This runs after individual CYS detection so we can:
        # 1. Skip CYS already detected as individual centers
        # 2. Flag CYS centers that are disulfide-bonded
        disulfide_scanner = DisulfideBondScanner(self.config)
        disulfide_centers = disulfide_scanner.scan_structure(
            structure,
            all_centers,
            selected_chains,
            source_pdb_file=self.source_pdb_file
        )
        all_centers.extend(disulfide_centers)

        # Store disulfide scanner for later manual specification
        self.disulfide_scanner = disulfide_scanner

        return all_centers

    def _prompt_manual_disulfide_specification(self, structure: Structure,
                                               selected_chains: List[str]) -> List[RedoxCenter]:
        """Interactive interface for manually specifying disulfide bonds."""

        all_cys = self.disulfide_scanner.get_all_cys_residues(structure, selected_chains)

        if not all_cys:
            self.console.print("[yellow]No CYS/CYX residues found in structure.[/yellow]")
            return []

        # CYS index table — indices here are what the user types in the pair prompt.
        self.console.print("\n[bold cyan]CYS/CYX residues:[/bold cyan]\n")
        cys_table = Table()
        cys_table.add_column("Index", style="cyan", width=6)
        cys_table.add_column("Chain", style="green", width=6)
        cys_table.add_column("Residue", style="yellow", width=8)
        cys_table.add_column("Type", style="magenta", width=8)
        idx_to_cys = {}
        for idx, (chain, resid, resname) in enumerate(all_cys, 1):
            cys_table.add_row(str(idx), chain, str(resid), resname)
            idx_to_cys[idx] = (chain, resid, resname)
        self.console.print(cys_table)

        # A disulfide requires two cysteines. With only one CYS/CYX there are no
        # pairs to rank or specify, so skip the (necessarily empty) distance
        # table and the pair prompt -- just note why and return.
        if len(all_cys) < 2:
            self.console.print(
                "\n[grey50]Only one CYS/CYX residue is present -- a disulfide bond "
                "needs two, so there are no pairs to specify.[/grey50]"
            )
            return []

        # Pairwise SG–SG distances sorted ascending. Lets the user see at a glance
        # which pairs are plausibly disulfide-bonded vs. pure background, including
        # any near-misses that fell outside the auto-detection cutoff. Capped at
        # 10 Å because past that, no χ-angle rotation could bring the SGs into
        # bonding range without rebuilding the backbone — so listing those pairs
        # is just noise.
        DISPLAY_LIMIT_A = 10.0
        cutoff = self.config.disulfide_distance_threshold
        pair_rows = []
        for i in range(1, len(all_cys) + 1):
            for j in range(i + 1, len(all_cys) + 1):
                ch1, r1, _ = idx_to_cys[i]
                ch2, r2, _ = idx_to_cys[j]
                try:
                    d = self._measure_sg_distance(structure, ch1, r1, ch2, r2)
                except ValueError:
                    d = None
                pair_rows.append((i, j, ch1, r1, ch2, r2, d))
        # Measured pairs first (ascending), unmeasurable pairs at the end.
        pair_rows.sort(key=lambda r: (r[6] is None, r[6] if r[6] is not None else 0))

        # Always show pairs within the display limit and any unmeasurable pairs
        # (they signal a structural problem worth flagging). If every measured
        # pair exceeds the limit, fall back to showing the single closest one
        # so the user isn't left staring at an empty table.
        shown_rows = [r for r in pair_rows if r[6] is None or r[6] <= DISPLAY_LIMIT_A]
        hidden_rows = [r for r in pair_rows if r[6] is not None and r[6] > DISPLAY_LIMIT_A]
        if not any(r[6] is not None for r in shown_rows) and hidden_rows:
            shown_rows.insert(0, hidden_rows.pop(0))

        self.console.print(
            f"\n[bold cyan]Pairwise SG–SG distances[/bold cyan] "
            f"(auto-detect cutoff: {cutoff:.2f} Å, shown: ≤ {DISPLAY_LIMIT_A:.0f} Å):\n"
        )
        dist_table = Table()
        dist_table.add_column("Pair", style="cyan", width=8)
        dist_table.add_column("Residues", style="yellow")
        dist_table.add_column("SG–SG", justify="right", style="green")
        dist_table.add_column("Note", style="magenta")
        for i, j, ch1, r1, ch2, r2, d in shown_rows:
            pair_str = f"{i}-{j}"
            res_str = f"{ch1}:{r1} – {ch2}:{r2}"
            if d is None:
                dist_str = "n/a"
                note = "missing SG atom"
            else:
                dist_str = f"{d:.2f} Å"
                if d <= cutoff:
                    note = "≤ cutoff"
                elif d <= cutoff + 1.0:
                    note = "near cutoff"
                elif d > DISPLAY_LIMIT_A:
                    note = "closest pair (all > 10 Å)"
                else:
                    note = ""
            dist_table.add_row(pair_str, res_str, dist_str, note)
        self.console.print(dist_table)
        if hidden_rows:
            self.console.print(
                f"[grey50]({len(hidden_rows)} more pair(s) with SG–SG > {DISPLAY_LIMIT_A:.0f} Å hidden — "
                f"backbone would need rebuilding for a disulfide to form there.)[/grey50]"
            )

        # Single combined prompt — accepts 'none' (or empty/default) to skip,
        # or a comma-separated list of pair indices like '1-2, 3-4'.
        response = prompt_with_context(
            self.processor,
            "\nSpecify disulfide bonds as index pairs (e.g. '1-2, 3-4') or 'none'",
            default="none",
            module="Redox Site Detector",
            description="Manual disulfide bond specification"
        ).strip()

        if not response or response.lower() == "none":
            return []

        added_centers = []
        added_bonds = []  # dedupe within this entry

        for token in response.split(","):
            token = token.strip()
            if not token:
                continue
            parts = token.split("-")
            if len(parts) != 2:
                self.console.print(f"[yellow]Skipping '{token}': expected format 'i-j'[/yellow]")
                continue
            try:
                idx1, idx2 = int(parts[0]), int(parts[1])
            except ValueError:
                self.console.print(f"[yellow]Skipping '{token}': non-integer index[/yellow]")
                continue
            if idx1 not in idx_to_cys or idx2 not in idx_to_cys:
                self.console.print(f"[yellow]Skipping '{token}': index out of range 1-{len(all_cys)}[/yellow]")
                continue
            if idx1 == idx2:
                self.console.print(f"[yellow]Skipping '{token}': cannot bond a CYS to itself[/yellow]")
                continue

            chain1, res1, _ = idx_to_cys[idx1]
            chain2, res2, _ = idx_to_cys[idx2]
            bond_key = tuple(sorted([(chain1, res1), (chain2, res2)]))
            if bond_key in added_bonds:
                self.console.print(f"[yellow]Skipping '{token}': duplicate of an already-specified pair[/yellow]")
                continue

            try:
                distance = self._measure_sg_distance(structure, chain1, res1, chain2, res2)
            except Exception as e:
                self.console.print(f"[red]Skipping '{token}': could not measure SG–SG distance ({e})[/red]")
                continue

            center1 = self.disulfide_scanner._create_cys_center(
                structure, chain1, res1, distance, chain2, res2
            )
            center2 = self.disulfide_scanner._create_cys_center(
                structure, chain2, res2, distance, chain1, res1
            )
            if not (center1 and center2):
                self.console.print(f"[red]Skipping '{token}': could not create centers[/red]")
                continue

            added_centers.extend([center1, center2])
            added_bonds.append(bond_key)
            warn = "  [yellow](>3.0 Å)[/yellow]" if distance > 3.0 else ""
            self.console.print(
                f"[green]✓ Added disulfide:[/green] "
                f"{chain1}:{res1} ↔ {chain2}:{res2} ({distance:.2f} Å){warn}"
            )

        return added_centers

    def _measure_sg_distance(self, structure: Structure,
                            chain1: str, res1: int,
                            chain2: str, res2: int) -> float:
        """Measure SG-SG distance between two CYS residues"""
        for model in structure:
            if chain1 in model and chain2 in model:
                cys1 = None
                cys2 = None

                for residue in model[chain1]:
                    if residue.id[1] == res1 and residue.resname in ['CYS', 'CYX']:
                        cys1 = residue
                        break

                for residue in model[chain2]:
                    if residue.id[1] == res2 and residue.resname in ['CYS', 'CYX']:
                        cys2 = residue
                        break

                if cys1 and cys2 and 'SG' in cys1 and 'SG' in cys2:
                    return cys1['SG'] - cys2['SG']

        raise ValueError(f"Could not measure distance for {chain1}:{res1} - {chain2}:{res2}")


    def _add_atoms_to_site(self, site: RedoxSite, site_centers: List[RedoxCenter], structure: Structure):
        """Add atoms from site residues to the RedoxSite object"""
        residue_keys = set()
        for center in site_centers:
            residue_keys.add((center.chain, center.resid, center.insertion_code))
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    res_key = (chain.id, residue.id[1], residue.id[2])
                    if res_key in residue_keys:
                        for atom in residue:
                            site_atom = RedoxSiteAtom(
                                chain=chain.id,
                                resname=residue.resname,
                                resid=residue.id[1],
                                atom_name=atom.name,
                                coords=tuple(round(x, 3) for x in atom.coord),
                                element=atom.element,
                                insertion_code=residue.id[2],
                                occupancy=atom.occupancy,
                                bfactor=atom.bfactor,
                                properties={'serial_number': atom.serial_number}
                            )
                            site.add_atom(site_atom)
    
    def _detect_site_bonds(self, site: RedoxSite, structure: Structure):
        """Bond detection for all atoms using configurable distance cutoff"""
        atoms = site.atoms
        
        logger.info(f"Detecting bonds for {len(atoms)} atoms in site {site.site_id}")
        
        bonds_detected = 0
        for i, atom1 in enumerate(atoms):
            for j, atom2 in enumerate(atoms):
                if i >= j:
                    continue
                
                distance = np.linalg.norm(np.array(atom1.coords) - np.array(atom2.coords))
                
                # Use configurable bond search distance
                if distance <= self.config.bond_search_distance:
                    # Classify the bond based on atom types
                    bond_type, chemical_type = classify_bond_types(
                        atom1.element, atom2.element, distance,
                        atom1.resname, atom2.resname,
                        atom1.resid, atom2.resid,
                        atom1.chain, atom2.chain,
                        atom1.atom_name, atom2.atom_name
                    )
                    site.add_bond_with_classification(atom1.coords, atom2.coords, distance, chemical_type)
                    bonds_detected += 1
        
        logger.info(f"Detected {bonds_detected} bonds in site {site.site_id}")
    
    def _parse_selection_input(self, input_str: str, max_num: int) -> Optional[List[int]]:
        """Parse selection input supporting ranges and comma-separated values
        
        Examples:
        - "1-10" -> [0,1,2,3,4,5,6,7,8,9]
        - "1,3,5-8" -> [0,2,4,5,6,7]
        - "1-5 10-15" -> [0,1,2,3,4,9,10,11,12,13,14]
        """
        try:
            indices = set()
            
            # Split by both spaces and commas
            parts = input_str.replace(',', ' ').split()
            
            for part in parts:
                if '-' in part:
                    # Handle range (e.g., "1-10")
                    range_parts = part.split('-')
                    if len(range_parts) == 2:
                        start = int(range_parts[0])
                        end = int(range_parts[1])
                        if start <= end:
                            for i in range(start, end + 1):
                                if 1 <= i <= max_num:
                                    indices.add(i - 1)  # Convert to 0-based
                                else:
                                    return None  # Invalid range
                        else:
                            return None  # Invalid range (start > end)
                    else:
                        return None  # Invalid range format
                else:
                    # Handle single number
                    num = int(part)
                    if 1 <= num <= max_num:
                        indices.add(num - 1)  # Convert to 0-based
                    else:
                        return None  # Invalid number
            
            return sorted(list(indices))
            
        except ValueError:
            return None

    # ========================================================================
    # Structure Viewer Integration
    # ========================================================================

    def _launch_quick_structure_viewer(self, centers_to_highlight: List[RedoxCenter],
                                      title: str = "Selected Centers") -> bool:
        """
        Launch structure viewer with specific centers highlighted.

        Args:
            centers_to_highlight: Centers to highlight in the viewer
            title: Title for the annotation

        Returns:
            True if viewer launched successfully
        """
        ngl_selection = self._centers_to_ngl_selection(centers_to_highlight)
        return _show_residues_in_viewer(
            getattr(self, 'processor', None),
            self.console,
            ngl_selection,
            label='selected_centers',
        )

    def _centers_to_ngl_selection(self, centers: List[RedoxCenter]) -> str:
        """
        Convert RedoxCenter list to NGL selection string.

        Args:
            centers: List of RedoxCenter objects

        Returns:
            NGL selection string like ":A and 154 or :A and 96"
        """
        selections = []
        for center in centers:
            selection = f":{center.chain} and {center.resid}"

            # Add altloc if present (and not blank/empty)
            if center.altloc and center.altloc.strip():
                selection += f" and %{center.altloc.strip()}"

            selections.append(selection)

        return " or ".join(selections) if selections else ""

    def _launch_viewer_for_residues(self, residue_list: List[Tuple], title: str = "Residues") -> bool:
        """
        Launch structure viewer with specific residues highlighted.

        Args:
            residue_list: List of (chain, resname, resid, icode) tuples
            title: Title for the annotation

        Returns:
            True if viewer launched successfully
        """
        # Note: this variant deliberately keeps protein/ligand/ion/water
        # defaults visible (focused=False) so the highlighted residues are
        # shown in their structural context — matches the previous behavior
        # of this method (no viewer_config flags set before migration).
        selections = [f":{chain} and {resid}" for chain, resname, resid, icode in residue_list]
        ngl_selection = " or ".join(selections) if selections else ""
        return _show_residues_in_viewer(
            getattr(self, 'processor', None),
            self.console,
            ngl_selection,
            label='selected_residues',
            focused=False,
        )


# ===== COMMAND LINE INTERFACE =====

def list_json_files(directory="."):
    """
    List all JSON files in the given directory.

    Args:
        directory: Directory path to search for JSON files

    Returns:
        List of absolute paths to JSON files
    """
    import os
    import glob
    from pathlib import Path

    json_files = []
    search_dir = Path(directory).absolute()

    # Search for JSON files with various patterns
    for pattern in ["*redox*.json", "*_redox_sites.json", "*.json"]:
        json_files.extend(glob.glob(str(search_dir / pattern)))

    # Remove duplicates, sort, and return absolute paths
    return sorted(list(set([os.path.abspath(f) for f in json_files])))


def display_json_file_menu(directory=".", console=None, processor=None):
    """
    Display a menu of available JSON files with directory navigation.
    Uses consistent interface matching display_pdb_file_menu.

    Args:
        directory: Directory path to search for JSON files (default: current directory)
        console: Rich console object for display (optional)
        processor: Processor object for session recording (optional)

    Returns:
        str or None: Selected file path or None if canceled
    """
    # Thin wrapper over the shared interactive browser: unified bare-N / q UX,
    # comma+range multi disabled here (single JSON pick), filename-based replay.
    from proprep.utils.file_browser import file_browser, default_size_detail

    return file_browser(
        directory=directory,
        extensions=[".json"],
        console=console,
        processor=processor,
        label="JSON file",
        entry_detail=default_size_detail,
        allow_path_jump=True,
        module="Redox Site Detector - JSON Import",
    )


def _prompt_for_json_import(console=None, processor=None) -> Optional[str]:
    """
    Prompt user to import JSON file using interactive file browser.

    Args:
        console: Rich console object for display (optional)
        processor: Processor object for session recording (optional)

    Returns:
        str or None: Selected JSON file path or None if canceled
    """
    if console is None:
        console = Console()

    # Ask if user wants to import
    if processor:
        choice = prompt_with_context(
            processor=processor,
            prompt="Do you want to import redox sites from a JSON file?",
            choices=["y", "n"],
            default="n",
            module="Redox Site Detector",
            description="Import redox sites from JSON file",
            options_map={"y": "Yes, import from JSON", "n": "No, run detection"}
        )
    else:
        choice = confirm_with_context(None, "Do you want to import redox sites from a JSON file?", default=False)
        choice = "y" if choice else "n"

    if choice.lower() in ['n', 'no']:
        return None

    # Use the file browser to select JSON file
    return display_json_file_menu(directory=".", console=console, processor=processor)


def _import_from_json(json_file: str) -> Tuple[List[RedoxSite], Dict[str, str]]:
    """Import redox sites from a JSON file.

    Returns:
        (sites, transformer_mappings) — the second is the site-type → transformer-name
        map persisted by _export_to_json. Missing in JSONs written by older ProPrep
        versions, in which case an empty dict is returned.
    """
    import json
    import os

    if not os.path.exists(json_file):
        raise FileNotFoundError(f"JSON file not found: {json_file}")

    with open(json_file, 'r') as f:
        data = json.load(f)

    sites = [dict_to_redox_site(site_data) for site_data in data.get("sites", [])]
    transformer_mappings = data.get("transformer_mappings", {}) or {}

    print(f"Successfully imported {len(sites)} redox site(s) from {json_file}")
    return sites, transformer_mappings


def dict_to_redox_site(site_data):
    """Reconstruct a RedoxSite from its JSON/dict serialization.

    Inverse of the per-site dict written by _export_to_json /
    redox_transformation_manager.export_redox_sites_to_json. Idempotent: a
    RedoxSite (or any non-dict) is returned unchanged, so it doubles as a
    normalizer for consumers that may receive either form — detected_redox_sites
    round-trips through JSON workspace state, so a resumed session can hand
    object-expecting code dict-form sites.
    """
    if not isinstance(site_data, dict):
        return site_data

    site = RedoxSite(site_data["site_id"], site_data.get("structure_id", "imported"))
    site.site_type = site_data.get("site_type", "unknown")

    for center_data in site_data.get("centers", []):
        center = RedoxCenter(
            chain=center_data["chain"],
            resname=center_data["resname"],
            resid=center_data["resid"],
            atom_name=center_data.get("atom_name"),
            insertion_code=center_data.get("insertion_code", ""),
            altloc=center_data.get("altloc", ""),
            coords=tuple(round(x, 3) for x in center_data["coordinates"]),
            center_type=CenterType(center_data["center_type"]),
            element=center_data.get("element")
        )
        site.add_center(center)

    for atom_data in site_data.get("atoms", []):
        atom = RedoxSiteAtom(
            chain=atom_data["chain"],
            resname=atom_data["resname"],
            resid=atom_data["resid"],
            atom_name=atom_data["atom_name"],
            coords=tuple(round(x, 3) for x in atom_data["coordinates"]),
            element=atom_data["element"]
        )
        site.add_atom(atom)

    for bond_data in site_data.get("bonds", []):
        bond = RedoxSiteBond(
            atom1_coords=tuple(round(x, 3) for x in bond_data["atom1_coordinates"]),
            atom2_coords=tuple(round(x, 3) for x in bond_data["atom2_coordinates"]),
            bond_type=bond_data["bond_type"],
            chemical_type=bond_data.get("chemical_type", "unknown"),
            distance=bond_data["distance"],
            atom1_element=bond_data.get("atom1_element", ""),
            atom2_element=bond_data.get("atom2_element", ""),
            atom1_residue_info=bond_data["atom1"],
            atom2_residue_info=bond_data["atom2"],
            treatment=bond_data.get("treatment", "bonded")
        )
        site.bonds.append(bond)

    return site


def _export_to_json(sites: List[RedoxSite], pdb_file: str, console=None,
                    transformer_mappings: Optional[Dict[str, str]] = None) -> None:
    """Export redox sites to JSON format.

    transformer_mappings (site-type → transformer-name) is persisted at the
    top level so that a fresh ProPrep session importing this JSON can
    re-populate workspace['redox_transformer_mappings'] and let the
    transformation manager auto-assign without re-running detection.
    """
    import json
    import os

    # Create export data
    export_data = {
        "source_pdb": os.path.basename(pdb_file),
        "transformer_mappings": transformer_mappings or {},
        "sites": []
    }
    
    for site in sites:
        site_data = {
            "site_id": site.site_id,
            "structure_id": site.structure_id,
            "site_type": site.site_type,
            "centers": [],
            "atoms": [],
            "bonds": []
        }
        
        # Export centers
        for center in site.centers:
            center_data = {
                "chain": center.chain,
                "resname": center.resname,
                "resid": center.resid,
                "atom_name": center.atom_name,
                "element": center.element,
                "center_type": center.center_type.value,
                "coordinates": [round(float(x), 3) for x in center.coords]
            }
            site_data["centers"].append(center_data)
        
        # Export atoms
        for atom in site.atoms:
            atom_data = {
                "chain": atom.chain,
                "resname": atom.resname,
                "resid": atom.resid,
                "atom_name": atom.atom_name,
                "element": atom.element,
                "coordinates": [round(float(x), 3) for x in atom.coords]
            }
            site_data["atoms"].append(atom_data)
        
        # Export bonds
        for bond in site.bonds:
            bond_data = {
                "atom1": bond.atom1_residue_info,
                "atom2": bond.atom2_residue_info,
                "bond_type": bond.bond_type,
                "chemical_type": bond.chemical_type,
                "distance": float(bond.distance),
                "atom1_element": bond.atom1_element,
                "atom2_element": bond.atom2_element,
                "atom1_coordinates": [round(float(x), 3) for x in bond.atom1_coords],
                "atom2_coordinates": [round(float(x), 3) for x in bond.atom2_coords],
                "treatment": getattr(bond, "treatment", "bonded")
            }
            site_data["bonds"].append(bond_data)
        
        export_data["sites"].append(site_data)
    
    # Write to file
    base_name = os.path.splitext(os.path.basename(pdb_file))[0]
    output_file = f"{base_name}_redox_sites.json"
    
    with open(output_file, 'w') as f:
        json.dump(export_data, f, indent=2)
    
    if console:
        console.print(f"[green]Exported redox sites to {output_file}[/green]")
    else:
        print(f"Exported redox sites to {output_file}")


def _export_to_pdb_single(sites: List[RedoxSite], pdb_file: str, console=None) -> None:
    """Export all redox sites to a single PDB file"""
    import os
    
    base_name = os.path.splitext(os.path.basename(pdb_file))[0]
    output_file = f"{base_name}_redox_sites.pdb"
    
    with open(output_file, 'w') as f:
        f.write("HEADER    REDOX SITES EXTRACTED\n")
        f.write(f"REMARK   1 Source: {os.path.basename(pdb_file)}\n")
        f.write(f"REMARK   2 Total sites: {len(sites)}\n")
        
        atom_num = 1
        for site_idx, site in enumerate(sites, 1):
            f.write(f"REMARK {site_idx + 2:3d} Site {site_idx}: {site.site_id}\n")
            
            # Write atoms for this site
            for atom in site.atoms:
                f.write(f"ATOM  {atom_num:5d}  {atom.atom_name:<4s}{atom.resname:>3s} {atom.chain}{atom.resid:>4d}    {atom.coords[0]:8.3f}{atom.coords[1]:8.3f}{atom.coords[2]:8.3f}  1.00 99.00           {atom.element:>2s}\n")
                atom_num += 1
            
            f.write(f"TER\n")
        
        f.write("END\n")
    
    if console:
        console.print(f"[green]Exported all redox sites to {output_file}[/green]")
    else:
        print(f"Exported all redox sites to {output_file}")


def _export_to_pdb_separate(sites: List[RedoxSite], pdb_file: str, console=None) -> None:
    """Export each redox site to a separate PDB file"""
    import os
    
    base_name = os.path.splitext(os.path.basename(pdb_file))[0]
    
    for site_idx, site in enumerate(sites, 1):
        output_file = f"{base_name}_site_{site_idx}.pdb"
        
        with open(output_file, 'w') as f:
            f.write("HEADER    REDOX SITE\n")
            f.write(f"REMARK   1 Source: {os.path.basename(pdb_file)}\n")
            f.write(f"REMARK   2 Site ID: {site.site_id}\n")
            f.write(f"REMARK   3 Centers: {len(site.centers)}\n")
            f.write(f"REMARK   4 Atoms: {len(site.atoms)}\n")
            f.write(f"REMARK   5 Bonds: {len(site.bonds)}\n")
            
            # Write atoms
            atom_num = 1
            for atom in site.atoms:
                f.write(f"ATOM  {atom_num:5d}  {atom.atom_name:<4s}{atom.resname:>3s} {atom.chain}{atom.resid:>4d}    {atom.coords[0]:8.3f}{atom.coords[1]:8.3f}{atom.coords[2]:8.3f}  1.00 99.00           {atom.element:>2s}\n")
                atom_num += 1
            
            f.write("END\n")
        
        if console:
            console.print(f"[green]Exported site {site_idx} to {output_file}[/green]")
        else:
            print(f"Exported site {site_idx} to {output_file}")


def main():
    """Main command line interface"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Comprehensive Redox Site Detector")
    parser.add_argument("pdb_file", help="PDB file to analyze")
    parser.add_argument("--chains", help="Specific chains to analyze (comma-separated)")
    parser.add_argument("--no-interactive", action="store_true", help="Run in non-interactive mode")
    parser.add_argument("--verbose", action="store_true", help="Verbose output")
    parser.add_argument("--import-json", help="Import redox sites from JSON file instead of detecting")
    
    args = parser.parse_args()
    
    # Set up logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    
    # Check if importing from JSON
    json_file_to_import = args.import_json
    
    # If no JSON file specified but in interactive mode, prompt for import
    if not json_file_to_import and not args.no_interactive:
        json_file_to_import = _prompt_for_json_import()
    
    if json_file_to_import:
        try:
            sites = _import_from_json(json_file_to_import)
        except Exception as e:
            print(f"Error importing JSON: {e}")
            return 0
    else:
        # Parse chain selection
        selected_chains = None
        if args.chains:
            selected_chains = [c.strip() for c in args.chains.split(',')]
        
        # Run detection
        detector = ComprehensiveRedoxDetector()
        sites = detector.detect_redox_sites(
            args.pdb_file, 
            selected_chains=selected_chains,
            interactive=not args.no_interactive
        )
    
    if sites:
        print(f"\n=== FINAL RESULTS ===")
        if json_file_to_import:
            print(f"Imported {len(sites)} redox site(s)")
        else:
            print(f"Detected {len(sites)} redox site(s)")
        
        for i, site in enumerate(sites, 1):
            print(f"\n--- Site {i} ({site.site_id}) ---")
            
            # Show redox centers
            print(f"Redox Centers ({len(site.centers)}):")
            for center in site.centers:
                if center.atom_name:
                    print(f"  • {center.resname} {center.chain}:{center.resid} {center.atom_name} ({center.center_type.value})")
                else:
                    print(f"  • {center.resname} {center.chain}:{center.resid} ({center.center_type.value})")
            
            # Show all atoms in site
            print(f"Site Atoms ({len(site.atoms)}):")
            by_residue = {}
            for atom in site.atoms:
                res_key = f"{atom.chain}:{atom.resname}{atom.resid}"
                if res_key not in by_residue:
                    by_residue[res_key] = []
                by_residue[res_key].append(atom.atom_name)
            
            for res_key in sorted(by_residue.keys()):
                atoms_str = ", ".join(sorted(by_residue[res_key]))
                print(f"  • {res_key}: {atoms_str}")
            
            # Show bonds
            print(f"Bonds ({len(site.bonds)}):")
            if site.bonds:
                for bond in site.bonds:
                    atom1_info = bond.atom1_residue_info
                    atom2_info = bond.atom2_residue_info
                    print(f"  • {atom1_info['chain']}:{atom1_info['resname']}{atom1_info['resid']} {atom1_info['atom_name']} — {atom2_info['chain']}:{atom2_info['resname']}{atom2_info['resid']} {atom2_info['atom_name']} ({bond.bond_type}, {bond.distance:.2f}Å)")
            else:
                print("  • No bonds defined")
        
        # Export options
        print(f"\n=== EXPORT OPTIONS ===")
        while True:
            print("1. Export to JSON")
            print("2. Export to PDB (single file)")
            print("3. Export to separate PDB files")
            print("4. Continue without export")
            
            choice = prompt_with_context(
                self.processor,
                "Select export option (1-4)",
                module="Redox Site Detector",
                description="Select export option",
                options_map={
                    "1": "Export to JSON",
                    "2": "Export to PDB (single file)",
                    "3": "Export to separate PDB files",
                    "4": "Continue without export",
                },
            ).strip()
            
            if choice == "1":
                _export_to_json(sites, args.pdb_file, None)
            elif choice == "2":
                _export_to_pdb_single(sites, args.pdb_file, None)
            elif choice == "3":
                _export_to_pdb_separate(sites, args.pdb_file, None)
            elif choice == "4":
                break
            else:
                print("Invalid choice. Please select 1-4.")
    
    return len(sites)

if __name__ == "__main__":
    main()