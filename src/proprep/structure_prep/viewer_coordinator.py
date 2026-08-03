"""ViewerCoordinator — terminal-as-driver for the Interactive Structure Viewer.

Workflow modules call this coordinator's intent-level methods (``show_structure``,
``highlight``, ``focus_on``) at natural moments — load complete, redox site
detected, residue picked, etc. The coordinator decides whether to launch a
viewer, push a live update to an existing one, or silently no-op based on
the current process state and execution mode.

Auto-launch policy (mode-dependent default):

- **Web shell mode** (``PROPREP_WEB_SHELL=1`` set by ``proprep-web``): the
  iframe is right there, so the first relevant event auto-launches the
  viewer and the iframe sources itself via the existing announce channel.
- **CLI mode**: the first relevant event is a silent no-op — popping a
  browser tab unbidden during a CLI workflow is intrusive. The user is
  expected to launch the viewer explicitly via the Structure Viewer menu.
- **Both modes**: events that fire while a viewer is *already running*
  push a live update to it via the existing /version-poll channel.

This module never raises from public methods. A failed update is logged at
DEBUG and dropped — viewer problems must not stop a module from completing.

Usage::

    from proprep.structure_prep.viewer_coordinator import viewer

    # In a structure loader, after a successful download:
    viewer.show_structure(file_path)

    # In a redox detector, after detecting a site:
    viewer.highlight(":A and 56", style="ball+stick", color="element",
                     label="redox_site_1")

    # When a module wants to drop all annotations (e.g., new workflow stage):
    viewer.clear_annotations()
"""

from __future__ import annotations

import logging
import os
import threading
from typing import Optional

logger = logging.getLogger(__name__)


def _is_web_shell_mode() -> bool:
    """True if running under ``proprep-web`` (PROPREP_WEB_SHELL is set)."""
    return bool(os.environ.get("PROPREP_WEB_SHELL"))


# 12-colour qualitative palette for indexed highlights (Stage-2 inventory,
# Stage-7 search results, etc.). ColorBrewer "Paired" — distinguishable
# colour pairs that preserve identity across hue families. Indexed
# 1-based to match the ID columns in the detector's tables.
_PALETTE_PAIRED12 = [
    "#a6cee3", "#1f78b4",  # 1, 2  — light/dark blue
    "#b2df8a", "#33a02c",  # 3, 4  — light/dark green
    "#fb9a99", "#e31a1c",  # 5, 6  — light/dark red
    "#fdbf6f", "#ff7f00",  # 7, 8  — light/dark orange
    "#cab2d6", "#6a3d9a",  # 9, 10 — light/dark purple
    "#ffff99", "#b15928",  # 11,12 — yellow/brown
]


def palette_color(index: int) -> str:
    """Return a hex color for ``index`` (1-based) from the qualitative palette.

    Wraps modulo the palette length, so callers can pass an arbitrary
    number of indexed items without index-out-of-range concerns.
    """
    if index < 1:
        index = 1
    return _PALETTE_PAIRED12[(index - 1) % len(_PALETTE_PAIRED12)]


# Common CSS colour words → hex. NGL reads a representation ``color`` string
# first as a colour-*scheme* id (``element``, ``resname``, ...); a bare colour
# word matches no scheme and is silently dropped (the rep never renders), so
# named colours must be converted to hex to show up. Scheme ids are NOT in this
# map, so they pass through untouched. Mirrors the hex used by the viewer's own
# colour picker (interactive_structure_viewer color_map).
_NAMED_COLOR_HEX = {
    "red": "#ff0000", "green": "#00ff00", "blue": "#0000ff",
    "yellow": "#ffff00", "orange": "#ff8c00", "cyan": "#00ffff",
    "magenta": "#ff00ff", "white": "#ffffff", "black": "#000000",
    "purple": "#800080", "pink": "#ff69b4", "gray": "#808080",
    "grey": "#808080",
}


def _resolve_color(color: str) -> str:
    """Normalise a colour for NGL.

    ``palette:N`` → hex from the qualitative palette; a bare CSS colour word →
    its hex (NGL would otherwise treat it as an unknown scheme and drop the
    rep); scheme ids (``element``, ``resname``, ...) and hex pass through.
    """
    if not isinstance(color, str):
        return color
    if color.startswith("palette:"):
        try:
            return palette_color(int(color.split(":", 1)[1]))
        except (ValueError, IndexError):
            return color  # malformed; let NGL try
    return _NAMED_COLOR_HEX.get(color.lower(), color)


class ViewerCoordinator:
    """Process-wide singleton facade over ``InteractiveStructureViewer``.

    Holds one persistent ``InteractiveStructureViewer`` instance for the
    duration of the process so module calls accumulate state on the same
    viewer (current structure, current annotations) rather than each
    creating a fresh viewer that stomps the previous one.
    """

    def __init__(self) -> None:
        self._lock = threading.RLock()
        self._viewer = None  # InteractiveStructureViewer, created lazily

    # ----- public intent methods -------------------------------------------

    def show_structure(
        self, file_path: str, *, name: Optional[str] = None, force: bool = False,
    ) -> None:
        """Display ``file_path`` as the current structure.

        If a viewer is already running, swaps to this structure (relaunching
        the underlying server is necessary because the server is bound to
        a fixed structure-file list at construction time). If no viewer is
        running, launches one only in web shell mode — *unless* ``force``
        is True, which makes the call always launch (used by explicit
        user-initiated commands like the redox detector's ``view`` prompt
        that must work in CLI mode too).

        Repeating the same path is a no-op so module-level idempotent
        callers don't churn the iframe.
        """
        if not file_path:
            return
        with self._lock:
            self._safely(lambda: self._show_structure_impl(file_path, name=name, force=force))

    def show_structures(
        self, file_paths, *, force: bool = False,
    ) -> None:
        """Display multiple structures overlaid in the viewer.

        Renders each path as its own NGL StructureComponent on the same
        stage — useful for alignment workflows where the user wants to
        see reference and aligned target(s) together. Replaces any
        prior single- or multi-structure view (and any annotations,
        which are stale once the structure list changes).

        ``force=True`` launches even in CLI mode if no viewer is running
        (mirrors ``show_structure``'s force flag for explicit user-
        initiated views). Empty list is a no-op rather than a clear —
        use ``clear_annotations`` for that.
        """
        if not file_paths:
            return
        paths = [str(p) for p in file_paths if p]
        if not paths:
            return
        with self._lock:
            self._safely(lambda: self._show_structures_impl(paths, force=force))

    def highlight(
        self,
        selection: str,
        *,
        style: str = "ball+stick",
        color: str = "element",
        label: Optional[str] = None,
        focused: bool = False,
        force: bool = False,
        opacity: Optional[float] = None,
        scale: Optional[float] = None,
    ) -> None:
        """Add an annotation overlay to the currently shown structure.

        ``selection`` is an NGL selection string (e.g., ``":A and 56"``).
        ``label`` lets later calls replace the same overlay; default labels
        are auto-generated from the selection.

        ``focused=True`` hides the default protein/ligand/ion/water reps so
        the overlay is the only thing visible — used when the caller is
        zooming in on a few specific residues. ``force=True`` launches the
        viewer even in CLI mode if none is running (mirrors
        ``show_structure``'s force flag for explicit user-initiated views).

        ``style="halo"`` renders a translucent yellow-ish spacefill sphere
        around each atom, similar to the measurement tool's pick marker. The
        halo is scaled off each atom's van der Waals radius, so it tracks
        atom size; the pick marker uses a fixed radius instead. Combine with a
        hex ``color`` for a coloured halo. ``opacity`` and ``scale`` override
        the style's defaults if supplied.
        """
        if not selection:
            return
        with self._lock:
            self._safely(lambda: self._highlight_impl(
                selection, style=style, color=color, label=label,
                focused=focused, force=force, opacity=opacity, scale=scale,
            ))

    def set_opacity(self, target: str, opacity: float, *, force: bool = False) -> None:
        """Make one of the default representations translucent.

        ``target`` is a default rep name — ``"protein"``, ``"ligands"``,
        ``"ions"``, ``"waters"``, ``"nonstandard"`` — or the full rep id
        (``"default_protein"``). ``opacity`` is 0.0 (invisible) to 1.0
        (solid).

        Use this when the thing worth looking at is *inside* the protein: a
        buried metal site or a channel is hidden by a solid cartoon, but
        reads clearly through one at ~0.3. Unlike ``highlight``, this
        changes an existing default rep rather than adding an overlay, and
        it leaves everything visible — so the site keeps its structural
        context instead of floating on its own the way ``focused=True``
        leaves it.

        The browser exposes the same control as a per-representation
        slider, so this only sets the starting state; the user can drag it
        from there.
        """
        with self._lock:
            self._safely(lambda: self._set_opacity_impl(
                target, opacity, force=force,
            ))

    def show_sphere(
        self,
        coords,
        *,
        radius: float,
        label: str,
        color: str = "#ffaa00",
        opacity: float = 0.5,
        force: bool = False,
    ) -> None:
        """Draw a translucent sphere at an arbitrary 3D coordinate.

        ``coords`` is an ``(x, y, z)`` triple in the structure's frame.
        The sphere is rendered via NGL's Shape API as its own component
        — independent of any structure atom — so it can mark a centroid
        or other geometric feature that isn't an atom in the loaded PDB.

        ``label`` is the unique identifier for this sphere; subsequent
        calls with the same label replace the previous sphere. Passing
        the label to ``unhighlight`` removes it.
        """
        if coords is None:
            self.unhighlight(label)
            return
        with self._lock:
            self._safely(lambda: self._show_sphere_impl(
                coords, radius=radius, label=label,
                color=color, opacity=opacity, force=force,
            ))

    def show_bonds(
        self,
        pairs,
        *,
        label: str,
        color: str = "#ffff00",
        show_labels: bool = False,
        force: bool = False,
    ) -> None:
        """Draw a bond-network overlay between pairs of NGL selections.

        ``pairs`` is a list of ``(selection_a, selection_b)`` tuples — each
        tuple resolves to a single bond drawn as a line between two
        atoms. Backed by NGL's ``distance`` representation.

        Useful for surfacing bonds that aren't in the structure file's
        CONECT records (e.g. coordinate bonds in a redox site that the
        detector defined explicitly). ``show_labels=True`` adds a bond-
        length label per line; off by default since 6+ labels per site
        are visually busy.
        """
        if not pairs:
            self.unhighlight(label)
            return
        with self._lock:
            self._safely(lambda: self._show_bonds_impl(
                pairs, label=label, color=color,
                show_labels=show_labels, force=force,
            ))

    def focus_on(self, selection: str) -> None:
        """Re-centre the viewer on the given NGL selection (no overlay)."""
        # Focusing is implemented as an annotation with focused_mode=True
        # so the existing per-residue refocus path (alt-loc picker) handles it.
        if not selection:
            return
        with self._lock:
            self._safely(lambda: self._focus_on_impl(selection))

    def clear_annotations(self) -> None:
        """Drop all coordinator-managed overlays and push a refresh."""
        with self._lock:
            self._safely(self._clear_annotations_impl)

    def unhighlight(self, label: str) -> None:
        """Remove a single labelled annotation, leaving others in place.

        Symmetric to a labelled ``highlight(...)``: pass the same label and
        the rep disappears. Used for "stop showing X" without disturbing
        any other coordinator-managed reps.
        """
        if not label:
            return
        with self._lock:
            self._safely(lambda: self._unhighlight_impl(label))

    def is_running(self) -> bool:
        """True if a viewer server is currently up (regardless of who started it)."""
        from proprep.structure_prep import interactive_structure_viewer as ivm
        srv = ivm._active_viewer_server
        return bool(srv is not None and srv.is_running())

    def set_processor(self, processor) -> None:
        """Plumb a processor into the underlying viewer for session recording.

        The viewer's interactive commands (LaunchViewerCommand,
        ShowAnnotationInfoCommand, etc.) record themselves through this
        processor reference. Callers that previously instantiated
        ``InteractiveStructureViewer`` directly used to call
        ``viewer.set_processor(...)`` to wire this up; this method exposes
        the same hook on the coordinator so the conversion is one-line.
        Silent no-op if anything goes wrong — viewer-state plumbing must
        never bubble into the caller's flow.
        """
        with self._lock:
            self._safely(lambda: self._ensure_viewer().set_processor(processor))

    # ----- internals --------------------------------------------------------

    def _ensure_viewer(self):
        """Lazily build the persistent ``InteractiveStructureViewer``."""
        if self._viewer is None:
            from proprep.structure_prep.interactive_structure_viewer import (
                InteractiveStructureViewer,
            )
            self._viewer = InteractiveStructureViewer()
            # Initialise the empty annotation/viewer-config/shape dicts so
            # later mutate-then-update_annotations calls don't AttributeError.
            self._viewer.annotation_config = {}
            self._viewer.viewer_config = {}
            self._viewer.shape_config = {}
        elif not hasattr(self._viewer, 'shape_config'):
            # Backstop for viewers constructed before shape_config existed.
            self._viewer.shape_config = {}
        return self._viewer

    def _show_structures_impl(
        self, file_paths, *, force: bool,
    ) -> None:
        v = self._ensure_viewer()
        same = (list(v.selected_structures) == list(file_paths))
        running = self.is_running()

        if same and running:
            if force:
                # Explicit user-initiated view: relaunch to guarantee a live
                # tab even when the structure list is unchanged (see the
                # matching note in _show_structure_impl).
                v._launch_viewer(open_browser=True)
            return

        v.selected_structures = list(file_paths)
        # New structure list — annotations / shapes are stale.
        v.annotation_config = {}
        v.viewer_config = {}
        v.shape_config = {}

        if running:
            # Relaunch only to swap the structure list. Open a fresh browser
            # tab solely when this swap was an explicit user view request
            # (force=True); a passive workflow waypoint reuses the existing
            # tab, which re-sources itself from the reused port.
            v._launch_viewer(open_browser=force)
            logger.debug(
                "ViewerCoordinator: swapped to overlay of %d structures", len(file_paths)
            )
            return
        if _is_web_shell_mode() or force:
            v._launch_viewer()
            logger.debug(
                "ViewerCoordinator: launched overlay of %d structures (web=%s, force=%s)",
                len(file_paths), _is_web_shell_mode(), force,
            )

    def _show_structure_impl(
        self, file_path: str, *, name: Optional[str], force: bool,
    ) -> None:
        v = self._ensure_viewer()
        same_structure = (v.selected_structures == [file_path])
        running = self.is_running()

        if same_structure and running:
            if force:
                # Explicit user-initiated view (e.g. the checklist "view
                # structure" command). The structure is unchanged, but the
                # user may have closed the browser tab -- or is_running() may
                # be reporting a stale/dead server -- so relaunch to guarantee
                # a live tab instead of silently no-op'ing.
                v._launch_viewer(open_browser=True)
            return  # already showing it; nothing to swap

        v.selected_structures = [file_path]
        # New structure means coordinator-managed annotations / shapes /
        # focused-mode flags are stale — caller must re-issue any
        # highlights / spheres.
        v.annotation_config = {}
        v.viewer_config = {}
        v.shape_config = {}

        if running:
            # Different structure on a live viewer — relaunch is the only
            # way to swap structure files in the underlying http.server
            # (its ``ViewerHTTPRequestHandler.structure_files`` is set at
            # construction). The iframe re-sources via the cache-buster in
            # app.js when the new server announces.
            # Relaunch only to swap the bound structure file. Pop a fresh
            # browser tab only for an explicit user view request (force=True);
            # a passive waypoint reuses the existing tab, which re-sources
            # itself from the reused port via app.js's cache-buster.
            v._launch_viewer(open_browser=force)
            logger.debug("ViewerCoordinator: swapped structure to %s", file_path)
            return

        if _is_web_shell_mode() or force:
            v._launch_viewer()
            logger.debug(
                "ViewerCoordinator: launched for %s (web=%s, force=%s)",
                file_path, _is_web_shell_mode(), force,
            )
            return

        # CLI mode + no running viewer + no force: silent.
        logger.debug(
            "ViewerCoordinator: deferring show_structure(%s) — CLI mode and no viewer running",
            file_path,
        )

    def _highlight_impl(
        self,
        selection: str,
        *,
        style: str,
        color: str,
        label: Optional[str],
        focused: bool,
        force: bool,
        opacity: Optional[float] = None,
        scale: Optional[float] = None,
    ) -> None:
        v = self._ensure_viewer()
        if not v.selected_structures:
            logger.debug("ViewerCoordinator.highlight: no structure shown yet — ignoring")
            return

        key = label or f"hl_{abs(hash(selection)) & 0xFFFF:04x}"
        resolved_color = _resolve_color(color)
        v.annotation_config = dict(v.annotation_config or {})
        cfg = {
            "count": 1,
            "ngl_selection": selection,
            "colors_per_structure": {0: resolved_color},
            "style": style,
            "priority": "HIGH",
            "structure_indices": [0],
        }
        if opacity is not None:
            cfg["opacity"] = opacity
        if scale is not None:
            cfg["scale"] = scale
        v.annotation_config[key] = cfg

        if focused:
            # Hide the default protein/ligand/ion/water reps so only the
            # overlay is visible. Reset on the next show_structure(...)
            # swap (handled in _show_structure_impl).
            v.viewer_config = dict(v.viewer_config or {})
            v.viewer_config.update({
                "show_protein_ribbon": False,
                "show_ligands": False,
                "show_ions": False,
                "show_waters": False,
            })

        if self.is_running():
            v.update_annotations()
        elif _is_web_shell_mode() or force:
            v._launch_viewer()
        # else: CLI silent — annotation persists in coordinator state and
        # will appear when the user explicitly launches the viewer.

    def _set_opacity_impl(
        self,
        target: str,
        opacity: float,
        *,
        force: bool,
    ) -> None:
        v = self._ensure_viewer()
        if not v.selected_structures:
            logger.debug(
                "ViewerCoordinator.set_opacity: no structure shown yet — ignoring"
            )
            return

        rep_id = str(target)
        if not rep_id.startswith("default_"):
            rep_id = f"default_{rep_id}"

        v.viewer_config = dict(v.viewer_config or {})
        rep_opacity = dict(v.viewer_config.get("rep_opacity") or {})
        rep_opacity[rep_id] = max(0.0, min(1.0, float(opacity)))
        v.viewer_config["rep_opacity"] = rep_opacity

        if self.is_running():
            v.update_annotations()
        elif _is_web_shell_mode() or force:
            v._launch_viewer()

    def _show_sphere_impl(
        self,
        coords,
        *,
        radius: float,
        label: str,
        color: str,
        opacity: float,
        force: bool,
    ) -> None:
        v = self._ensure_viewer()
        if not v.selected_structures:
            logger.debug("ViewerCoordinator.show_sphere: no structure shown yet — ignoring")
            return

        v.shape_config = dict(v.shape_config or {})
        v.shape_config[label] = {
            "type": "sphere",
            "coords": [float(coords[0]), float(coords[1]), float(coords[2])],
            "radius": float(radius),
            "color": color,
            "opacity": float(opacity),
        }

        if self.is_running():
            v.update_annotations()
        elif _is_web_shell_mode() or force:
            v._launch_viewer()

    def _show_bonds_impl(
        self,
        pairs,
        *,
        label: str,
        color: str,
        show_labels: bool,
        force: bool,
    ) -> None:
        v = self._ensure_viewer()
        if not v.selected_structures:
            logger.debug("ViewerCoordinator.show_bonds: no structure shown yet — ignoring")
            return

        atom_pairs = [[str(a), str(b)] for (a, b) in pairs]
        v.annotation_config = dict(v.annotation_config or {})
        v.annotation_config[label] = {
            "count": len(atom_pairs),
            "ngl_selection": "",  # unused for bond style
            "atom_pairs": atom_pairs,
            "show_labels": show_labels,
            "colors_per_structure": {0: color},
            "style": "bond",
            "priority": "HIGH",
            "structure_indices": [0],
        }

        if self.is_running():
            v.update_annotations()
        elif _is_web_shell_mode() or force:
            v._launch_viewer()

    def _focus_on_impl(self, selection: str) -> None:
        v = self._ensure_viewer()
        v.viewer_config = dict(v.viewer_config or {})
        # ``focused_mode=True`` makes the iframe call NGL's
        # ``component.autoView(<first visible ann_ rep selection>)``
        # whenever a config update lands. We also keep the original
        # selection in viewer_config for any future caller that wants
        # an explicit per-selection autoView (currently the iframe just
        # uses the first visible annotation, so this is informational).
        v.viewer_config["focused_mode"] = True
        v.viewer_config["focus_selection"] = selection
        if self.is_running():
            v.update_annotations()
        elif _is_web_shell_mode():
            v._launch_viewer()

    def _clear_annotations_impl(self) -> None:
        v = self._ensure_viewer()
        v.annotation_config = {}
        v.shape_config = {}
        if self.is_running():
            v.update_annotations()

    def _unhighlight_impl(self, label: str) -> None:
        v = self._ensure_viewer()
        changed = False
        if v.annotation_config and label in v.annotation_config:
            v.annotation_config = {
                k: cfg for k, cfg in v.annotation_config.items() if k != label
            }
            changed = True
        # The same label may identify a shape (e.g. a centroid sphere)
        # rather than an annotation; nuke from both stores so callers
        # don't have to know which kind of overlay carries the label.
        if getattr(v, 'shape_config', None) and label in v.shape_config:
            v.shape_config = {
                k: cfg for k, cfg in v.shape_config.items() if k != label
            }
            changed = True
        if changed and self.is_running():
            v.update_annotations()

    def _safely(self, fn) -> None:
        """Run ``fn`` and swallow exceptions at DEBUG level.

        Public coordinator calls must never raise into module workflows —
        a wedged viewer can't be allowed to abort an unrelated step.
        """
        try:
            fn()
        except Exception as exc:
            logger.debug("ViewerCoordinator silenced exception: %s", exc, exc_info=True)


# Module-level singleton. Callers do ``from ... import viewer`` and use
# ``viewer.show_structure(...)``. Tests can construct fresh instances.
viewer = ViewerCoordinator()


__all__ = ["ViewerCoordinator", "viewer", "palette_color"]
