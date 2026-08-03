"""Shared force-field selection menu rendering.

ProPrep offers force-field *choices* (which protein FF / water model to build
with) in two places: the Topology Generator and the component parameterizers.
Both draw from the same curated catalog (``forcefield_catalog.FORCEFIELD_OPTIONS``)
but historically each hand-wrote its own Rich table, so they diverged visually
and the parameterizer even kept a stale hardcoded copy of the option list.

This module is the single renderer for that curated catalog. It is deliberately
NOT used by the Force Field Explorer, which browses the full leaprc universe
(DNA/RNA/lipid/GAFF/water/...) for parameter inspection rather than offering a
buildable-protein choice.

Style (matches the Topology Generator): a bold-blue "━━━ title ━━━" divider, a
grey50 description line, and a headerless table (bold choice # / green name /
description) with an optional "None" row and inline recommendation markers.
"""

from typing import Any, Dict, List, Optional

from rich.console import Console
from rich.table import Table


def render_forcefield_category(
    console: Console,
    title: str,
    description: Optional[str],
    options: List[Dict[str, Any]],
    *,
    allow_none: bool = False,
    none_text: str = "Skip this category",
    start_index: int = 1,
) -> int:
    """Render one force-field category in the shared ProPrep style.

    Each option dict may carry ``name``, ``description``, ``recommended``,
    ``recommendation_reason`` and ``cofactor_satisfier`` (the last two drive the
    inline recommendation markers). "None" is always the first row when allowed
    so the user never has to count to a category-dependent final index.

    Returns the ``none_offset`` (1 if a "None" row was rendered, else 0) so the
    caller can map a typed 1-based choice back into ``options``.
    """
    console.print(f"\n[bold blue]━━━ {title} ━━━[/bold blue]")
    if description:
        console.print(f"[grey50]{description}[/grey50]")

    none_offset = 1 if allow_none else 0

    table = Table(show_header=False, box=None, padding=(0, 2))
    table.add_column("Choice", style="bold")
    table.add_column("Name", style="green")
    table.add_column("Description")

    if allow_none:
        table.add_row(str(start_index), "None", none_text)

    for i, opt in enumerate(options, start=start_index + none_offset):
        name = opt["name"]
        if opt.get("recommended"):
            reason = opt.get("recommendation_reason")
            # Cofactor-driven satisfiers render prominently (green ←) to show the
            # option covers a loaded cofactor's prerequisites; softer picker
            # recommendations use the yellow marker.
            if opt.get("cofactor_satisfier") and reason:
                name += f"  [bold green]← {reason}[/bold green]"
            elif reason:
                name += f" [yellow](Recommended — {reason})[/yellow]"
            else:
                name += " [yellow](Recommended)[/yellow]"
        table.add_row(str(i), name, opt.get("description", ""))

    console.print(table)
    return none_offset


def _catalog_options(category: str) -> List[Dict[str, Any]]:
    """Fetch a category's options from the curated catalog."""
    from proprep.forcefield_params.forcefield_catalog import FORCEFIELD_OPTIONS
    return FORCEFIELD_OPTIONS.get(category, {}).get("options", [])


def _default_choice(options: List[Dict[str, Any]], none_offset: int) -> str:
    """First recommended option's 1-based index, else the first real option."""
    for i, opt in enumerate(options, start=1 + none_offset):
        if opt.get("recommended"):
            return str(i)
    return str(1 + none_offset)


def select_single_forcefield(
    console: Console,
    processor: Any,
    category: str,
    *,
    module: str,
    prompt_label: str,
    title: Optional[str] = None,
    description: Optional[str] = None,
    options: Optional[List[Dict[str, Any]]] = None,
    allow_none: bool = False,
    none_text: str = "None / skip",
) -> Optional[Dict[str, Any]]:
    """Render a curated-catalog category and prompt for ONE selection.

    Pulls ``options`` from ``FORCEFIELD_OPTIONS[category]`` unless supplied.
    Defaults to the first recommended option. Returns the chosen option dict,
    or ``None`` when the user selects the "None" row (only offered when
    ``allow_none`` is True).
    """
    from proprep.utils.prompts import prompt_with_context

    if options is None:
        options = _catalog_options(category)
    if title is None:
        from proprep.forcefield_params.forcefield_catalog import FORCEFIELD_OPTIONS
        title = FORCEFIELD_OPTIONS.get(category, {}).get("title", category.upper())

    none_offset = render_forcefield_category(
        console, title, description, options,
        allow_none=allow_none, none_text=none_text,
    )

    default = _default_choice(options, none_offset)
    choice_str = prompt_with_context(
        processor, prompt_label, default=default,
        module=module, description=prompt_label,
    )

    try:
        choice = int(choice_str)
    except (ValueError, TypeError):
        return options[0] if options else None

    if allow_none and choice == 1:
        return None
    idx = choice - 1 - none_offset
    if 0 <= idx < len(options):
        return options[idx]
    return options[0] if options else None
