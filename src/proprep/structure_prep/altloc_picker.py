"""
Choosing one alternate location per residue, residue by residue.

Shared by the Structure Fixer (Structure Completeness) and Molecular Docking.
Each residue with alternates is shown with every alternate's occupancy, how
many of the residue's alternate-carrying atoms it models, and, with the live
viewer open, the colour it is drawn in. An alternate that models only part of
the residue is flagged: choosing it keeps the atoms it lacks from the
alternate with the highest occupancy that has them (the fill plan), so the
residue stays complete. That is the case of PDB 8ZST's C-terminal LEU 129,
whose alternate C models only the side chain.

The callers own their structure formats: the Structure Fixer reads Biopython
structures, Docking reads PDB records. Both describe a residue's alternates as
a ``ResidueAlternates`` and the viewer's environment shell as plain points and
their residues, so the display, prompt, colours and fill rules are one code.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from proprep.utils import prompts

logger = logging.getLogger(__name__)

# Per-altloc colours used in the 3D viewer, with the name the prompt prints
# next to each alternate so the user can tell which occupancy is which.
PALETTE = {
    "A": ("#e74c3c", "red"),
    "B": ("#3498db", "blue"),
    "C": ("#2ecc71", "green"),
    "D": ("#f39c12", "orange"),
    "E": ("#9b59b6", "purple"),
    "F": ("#1abc9c", "teal"),
}
FALLBACK_PALETTE = [
    ("#e67e22", "dark orange"), ("#34495e", "slate"),
    ("#c0392b", "dark red"), ("#16a085", "sea green"),
]
# Radius (Angstrom) of the environment shell drawn around each residue: the
# tight first-contact shell (clashes, H-bond partners, packing).
ENV_DISTANCE = 5.0

VIEWER_PROMPT = "Launch the 3D viewer to help pick alternate locations?"
PICK_PROMPT = "Select alternate to keep"

ResidueKey = Tuple[str, int, str]          # (chain, residue number, insertion code)


def altloc_color(alt: str, index: int) -> Tuple[str, str]:
    """(hex, name) for an altloc letter; letters beyond F cycle the fallback palette."""
    return PALETTE.get(alt.upper(), FALLBACK_PALETTE[index % len(FALLBACK_PALETTE)])


def fill_source(atom_letters: Dict[str, Dict[str, Any]]) -> Dict[str, str]:
    """For each atom carrying alternates, the letter to fall back on.

    Highest occupancy wins; ties go to the first letter alphabetically so the
    choice is reproducible across runs and replays.
    """
    source = {}
    for name, letters in atom_letters.items():
        if not letters:
            continue
        source[name] = min(letters, key=lambda l: (-(letters[l] if letters[l] is not None else -1.0), l))
    return source


def fill_plan(atom_letters: Dict[str, Dict[str, Any]], letter: str) -> Dict[str, str]:
    """Atoms ``letter`` does not model, and the alternate each comes from.

    Empty when the chosen alternate models every atom that carries alternates,
    which is the ordinary case.
    """
    source = fill_source(atom_letters)
    return {name: source[name] for name, letters in atom_letters.items()
            if letter not in letters and name in source}


@dataclass
class ResidueAlternates:
    """One residue's alternate locations: which atoms each letter models, at what occupancy."""

    chain: str
    resseq: int
    icode: str
    resname: str
    atom_letters: Dict[str, Dict[str, Optional[float]]] = field(default_factory=dict)   # atom -> {letter: occ}
    occupancies: Dict[str, List[float]] = field(default_factory=dict)                   # letter -> occupancies

    def add(self, atom_name: str, letter: str, occupancy: Optional[float]) -> None:
        self.occupancies.setdefault(letter, [])
        if occupancy is not None:
            self.occupancies[letter].append(occupancy)
        self.atom_letters.setdefault(atom_name, {})[letter] = occupancy

    @property
    def key(self) -> ResidueKey:
        return (self.chain, self.resseq, self.icode)

    @property
    def number(self) -> str:
        """Residue number with any insertion code, as the PDB writes it (40, 40A)."""
        return f"{self.resseq}{self.icode}"

    @property
    def letters(self) -> List[str]:
        return sorted(self.occupancies)

    def average_occupancy(self, letter: str) -> str:
        values = self.occupancies.get(letter) or []
        return f"{sum(values) / len(values):.2f}" if values else "?"

    def missing(self, letter: str) -> List[str]:
        """Alternate-carrying atoms that ``letter`` does not model."""
        return sorted(n for n, letters in self.atom_letters.items() if letter not in letters)

    def fills(self, letter: str) -> Dict[str, str]:
        return fill_plan(self.atom_letters, letter) if self.missing(letter) else {}


# ------------------------------------------------------------------ viewer ---

def ask_to_open_viewer(processor, module: str) -> bool:
    """Ask before launching: in CLI mode the picker must not pop a browser tab unbidden."""
    return prompts.confirm_with_context(
        processor=processor,
        prompt=VIEWER_PROMPT,
        default=False,
        module=module,
        description="Optionally launch the structure viewer with per-residue refocus to aid alt-loc selection",
    )


def open_viewer(console, pdb_path: str, points: Sequence[Sequence[float]],
                owners: Sequence[ResidueKey]) -> Optional[Dict[str, Any]]:
    """Show ``pdb_path`` (which must still carry its alternates) in the live viewer.

    ``points`` and ``owners`` are every atom's coordinates and residue, used
    to draw the environment shell around each residue. Returns the viewer
    state for ``focus_viewer``, or None when the viewer cannot be launched;
    the picker works the same without it.
    """
    try:
        import numpy as np
        from proprep.structure_prep.viewer_coordinator import viewer as _viewer

        # force=True: the user has just opted in, so this is a user-initiated view.
        _viewer.show_structure(pdb_path, force=True)
        console.print("[grey50]Live 3D viewer is open in your browser; it will refocus "
                      "on each residue as you choose.[/grey50]")
        return {
            "pdb_path": pdb_path,
            "prev_labels": [],
            "env_distance": ENV_DISTANCE,
            "points": np.asarray(points, dtype=float).reshape(-1, 3),
            "owners": list(owners),
        }
    except Exception as e:
        logger.debug(f"Live alt-loc viewer unavailable: {e}")
        return None


def _ngl_residue(chain: str, resseq: int, icode: str = "") -> str:
    return f":{chain} and {resseq}" + (f"^{icode}" if icode else "")


def environment_selection(viewer_state: Dict[str, Any], key: ResidueKey) -> Optional[str]:
    """NGL selection of the residues within ``env_distance`` of residue ``key``, or None.

    Measured from every atom of the residue, all alternates included, and
    grouped by chain, e.g. ``(:A and (54 or 55 or 90)) or (:B and (12))``.
    """
    env_distance = viewer_state.get("env_distance") or 0.0
    points = viewer_state.get("points")
    owners = viewer_state.get("owners") or []
    if not env_distance or points is None or not len(owners):
        return None
    try:
        import numpy as np
        mine = np.array([owner == key for owner in owners])
        if not mine.any():
            return None
        target = points[mine]
        near = np.zeros(len(owners), dtype=bool)
        for xyz in target:
            near |= np.einsum("ij,ij->i", points - xyz, points - xyz) <= env_distance ** 2
        by_chain: Dict[str, set] = {}
        for owner, close in zip(owners, near):
            if close and owner != key:
                by_chain.setdefault(owner[0], set()).add(owner[1:])
        if not by_chain:
            return None
        groups = []
        for chain in sorted(by_chain):
            numbers = " or ".join(f"{n}^{i}" if i else str(n) for n, i in sorted(by_chain[chain]))
            groups.append(f"(:{chain} and ({numbers}))")
        return " or ".join(groups)
    except Exception as e:
        logger.debug(f"Could not compute alt-loc environment: {e}")
        return None


def focus_viewer(viewer_state: Dict[str, Any], key: ResidueKey, letters: Sequence[str],
                 occupancies: Optional[Dict[str, str]] = None) -> None:
    """Refocus the viewer on one residue's alternates.

    Clears the previous residue's reps, then draws a grey licorice scaffold of
    the whole residue (focused, so the camera centres on it), one ball+stick
    rep per alternate in its colour (NGL ``%A`` selects altloc A), and the
    environment shell as a faint non-focused line overlay.
    """
    try:
        from proprep.structure_prep.viewer_coordinator import viewer as _viewer

        for stale_label in viewer_state.get("prev_labels", []):
            _viewer.unhighlight(stale_label)

        base = _ngl_residue(*key)
        new_labels = ["altloc_scaffold"]
        _viewer.highlight(base, style="licorice", color="#bdc3c7", label="altloc_scaffold", focused=True)
        for i, alt in enumerate(letters):
            color, color_name = altloc_color(alt, i)
            label = f"altloc_{alt}"
            new_labels.append(label)
            occ = (occupancies or {}).get(alt)
            display = f"Alt {alt} ({color_name})" + (f", occ {occ}" if occ else "")
            _viewer.highlight(f"{base} and %{alt}", style="ball+stick", color=color, label=label,
                              display_label=display)

        env_selection = environment_selection(viewer_state, key)
        if env_selection:
            _viewer.highlight(env_selection, style="line", color="#7f8c8d", label="altloc_environment",
                              opacity=0.6)
            new_labels.append("altloc_environment")
        viewer_state["prev_labels"] = new_labels
    except Exception as e:
        logger.debug(f"Could not refocus alt-loc viewer: {e}")


def close_viewer(viewer_state: Optional[Dict[str, Any]]) -> None:
    """Clear the per-altloc reps so they do not leak into later views; the structure stays loaded."""
    if not viewer_state:
        return
    try:
        from proprep.structure_prep.viewer_coordinator import viewer as _viewer
        for stale_label in viewer_state.get("prev_labels", []):
            _viewer.unhighlight(stale_label)
        viewer_state["prev_labels"] = []
    except Exception as e:
        logger.debug(f"Could not tear down alt-loc viewer: {e}")


# ------------------------------------------------------------------ picker ---

def pick(processor, console, residue: ResidueAlternates, viewer_state: Optional[Dict[str, Any]], *,
         module: str, default_letter: Optional[str] = None) -> Tuple[str, Dict[str, str]]:
    """Show one residue's alternates, ask which to keep; return (letter, fill plan).

    The default is the first alternate, or ``default_letter`` when the caller
    has an earlier choice to offer again.
    """
    letters = residue.letters
    n_lettered = len(residue.atom_letters)
    source = fill_source(residue.atom_letters)
    occupancies: Dict[str, str] = {}

    console.print(f"[bold]Chain {residue.chain}, {residue.resname} {residue.number}:[/bold]")
    for i, letter in enumerate(letters, 1):
        occupancies[letter] = residue.average_occupancy(letter)
        missing = residue.missing(letter)
        detail = f"occupancy: {occupancies[letter]}, {n_lettered - len(missing)} of {n_lettered} atoms"
        if viewer_state is not None:
            hex_color, color_name = altloc_color(letter, i - 1)
            console.print(f"  {i}. Alternate {letter} ({detail})  [{hex_color}]■ {color_name} in viewer[/{hex_color}]")
        else:
            console.print(f"  {i}. Alternate {letter} ({detail})")
        if missing:
            filled = ", ".join(f"{n} from {source[n]}" for n in missing if n in source)
            console.print(f"     [yellow]partial: does not model {', '.join(missing)}[/yellow]")
            if filled:
                console.print(f"     [grey50]choosing it takes {filled}[/grey50]")

    if viewer_state is not None:
        focus_viewer(viewer_state, residue.key, letters, occupancies=occupancies)

    default = str(letters.index(default_letter) + 1) if default_letter in letters else "1"
    choice = prompts.prompt_with_context(
        processor=processor,
        prompt=PICK_PROMPT,
        choices=[str(i) for i in range(1, len(letters) + 1)],
        default=default,
        module=module,
        description=f"Select alternate for {residue.resname} {residue.chain}:{residue.number}",
        options_map={str(i + 1): f"Alternate {letter}" for i, letter in enumerate(letters)},
    )
    chosen = letters[int(choice) - 1]
    gaps = residue.missing(chosen)
    plan = residue.fills(chosen)
    if plan:
        detail = ", ".join(f"{n} from {l}" for n, l in sorted(plan.items()))
        console.print(f"[yellow]  Alternate {chosen} does not model {', '.join(gaps)}; keeping {detail} "
                      f"so the residue stays complete.[/yellow]")
    console.print(f"[green]✓ Will keep alternate {chosen}[/green]\n")
    return chosen, plan


def keeps_atom(atom_name: str, altloc: str, chosen: Optional[str], fills: Optional[Dict[str, str]]) -> bool:
    """Whether an atom survives a residue's choice: unlabelled, the chosen alternate, or a fill from its source."""
    if not altloc:
        return True
    if altloc == chosen:
        return True
    return (fills or {}).get(atom_name) == altloc


def residues_from_records(records: Iterable[Tuple[str, int, str, str, str, str, Optional[float]]]
                          ) -> Dict[ResidueKey, ResidueAlternates]:
    """ResidueAlternates for every residue with labelled atoms.

    ``records`` are (chain, resseq, icode, resname, atom name, altloc, occupancy)
    for each atom; unlabelled atoms are skipped. Residues with a single letter
    are included: callers decide what a lone label means for them.
    """
    out: Dict[ResidueKey, ResidueAlternates] = {}
    for chain, resseq, icode, resname, name, altloc, occupancy in records:
        if not altloc:
            continue
        key = (chain, resseq, icode)
        if key not in out:
            out[key] = ResidueAlternates(chain, resseq, icode, resname)
        out[key].add(name, altloc, occupancy)
    return out
