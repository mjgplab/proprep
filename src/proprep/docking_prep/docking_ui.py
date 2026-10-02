"""
Prompts, tables and typed chemistry edits for the Molecular Docking menus.

Every question the docking menus ask goes through the helpers here, so the
wording lives in one place and stays fixed. Prompt text, table cells and notes
are escaped for Rich, which otherwise reads "[s]" in "[s] skip" as its
strikethrough tag and prints nothing of it (found in the running program): session replay matches recorded
answers to the exact prompt text. Numeric prompts show their default in the
text through ProPrep's retry helpers; defaults are constants, never computed
from the structure, for the same reason.

Chemistry edits are typed as short commands (also the keyboard alternative to
picking in the viewer):

    remove H2A          remove atom H2A
    charge O2A 0        set the formal charge of O2A
    bond C1 C2 2        set the C1-C2 bond order (1, 2 or 3)
    add-h O2A           add a hydrogen to O2A
"""

from __future__ import annotations

import re
from typing import Dict, Iterable, List, Optional, Sequence

from rich.markup import escape
from rich.table import Table

from proprep.utils import prompts

MODULE_NAME = "Molecular Docking"
BLUE = "#1f6feb"          # option keys, as in the other modules' dashboards


def ask(processor, text: str, *, description: str, choices: Optional[Sequence[str]] = None,
        default: Optional[str] = None, options_map: Optional[Dict[str, str]] = None) -> str:
    kwargs = dict(module=MODULE_NAME, description=description)
    if choices is not None:
        kwargs["choices"] = list(choices)
    if default is not None:
        kwargs["default"] = default
    if options_map is not None:
        kwargs["options_map"] = options_map
    return prompts.prompt_with_context(processor, escape(text), **kwargs).strip()


def ask_float(processor, text: str, *, default: float, description: str,
              min_value: Optional[float] = None, max_value: Optional[float] = None) -> float:
    return prompts.prompt_float_with_retry(processor, escape(text), default=default, module=MODULE_NAME,
                                           description=description, min_value=min_value, max_value=max_value)


def ask_int(processor, text: str, *, default: int, description: str,
            min_value: Optional[int] = None, max_value: Optional[int] = None) -> int:
    return prompts.prompt_int_with_retry(processor, escape(text), default=default, module=MODULE_NAME,
                                         description=description, min_value=min_value, max_value=max_value)


def ask_required_int(processor, console, text: str, *, description: str) -> int:
    """An integer with NO default: for choices the program must not make (a metal's charge)."""
    while True:
        answer = ask(processor, text, description=description)
        try:
            return int(answer)
        except ValueError:
            console.print(f"[red]'{answer}' is not a whole number.[/red]")


def confirm(processor, text: str, *, default: bool, description: str) -> bool:
    return prompts.confirm_with_context(processor, escape(text), default=default, module=MODULE_NAME,
                                        description=description)


# ---------------------------------------------------------------- parsing ---

def parse_list(text: str) -> List[str]:
    """'A:45, A:88' or 'A:45 A:88' -> ['A:45', 'A:88']."""
    return [item for item in re.split(r"[,\s]+", text.strip()) if item]


def parse_numbers(text: str, count: int) -> List[int]:
    """'1,3-5' -> [1, 3, 4, 5], each between 1 and ``count``; raises ValueError otherwise."""
    numbers: List[int] = []
    for part in parse_list(text):
        if "-" in part:
            low, high = part.split("-", 1)
            numbers.extend(range(int(low), int(high) + 1))
        else:
            numbers.append(int(part))
    bad = [n for n in numbers if not 1 <= n <= count]
    if bad:
        raise ValueError(f"{', '.join(map(str, bad))} not between 1 and {count}")
    return sorted(set(numbers))


EDIT_HELP = ("remove NAME | charge NAME VALUE | bond NAME1 NAME2 ORDER | add-h NAME")
# Takes back the last chemistry edit. Not "undo": typing undo at any prompt
# rewinds the whole session (utils.session_rewind), so it never reaches here.
DROP = "drop"
DROP_HELP = "'drop' to take back the last edit"


def parse_edit(text: str):
    """One typed edit command (see the module docstring) -> ChemistryEdit; ValueError if malformed."""
    from .chemistry_edits import ChemistryEdit
    words = text.split()
    if not words:
        raise ValueError("empty edit")
    verb = words[0].lower()
    if verb == "remove" and len(words) == 2:
        return ChemistryEdit("remove_atom", (words[1],))
    if verb == "charge" and len(words) == 3:
        return ChemistryEdit("set_formal_charge", (words[1],), int(words[2]))
    if verb == "bond" and len(words) == 4:
        return ChemistryEdit("set_bond_order", (words[1], words[2]), int(words[3]))
    if verb in ("add-h", "addh") and len(words) == 2:
        return ChemistryEdit("add_hydrogen", (words[1],))
    raise ValueError(f"'{text}' is not an edit; use: {EDIT_HELP}")


# ------------------------------------------------------------------- picks ---

def pick(m, kind: str, prompt: str, structure_index: int):
    """One click in the open viewer (see ViewerCoordinator.pick); None, with the reason said,
    when no viewer is open or the pick was cancelled. Typed input always remains."""
    from proprep.structure_prep.viewer_coordinator import viewer
    m.console.print(f"  Waiting for a click in the viewer: {prompt}. Cancel there, or press Ctrl-C here.",
                    highlight=False)
    result = viewer.pick(kind, prompt, structure_index=structure_index)
    if result is None:
        m.console.print("  [dark_orange3]No pick: the viewer is not open (open it from the Structure Viewer) "
                        "or the pick was cancelled. Type it instead.[/dark_orange3]", highlight=False)
    return result


# ----------------------------------------------------------------- output ---

def table(title: str, columns: Sequence[str], rows: Iterable[Sequence[object]]) -> Table:
    out = Table(title=title, header_style="bold", show_lines=False, title_justify="left")
    for column in columns:
        out.add_column(column)
    for row in rows:
        out.add_row(*[escape(str(value)) for value in row])
    return out


def print_notes(console, notes: Iterable[str], style: str = "") -> None:
    for note in (escape(n) for n in notes):
        console.print(f"  {'[' + style + ']' if style else ''}- {note}{'[/' + style + ']' if style else ''}",
                      highlight=False)
