"""
AMBER Parameter Display Utilities
===================================

Formatting and display layer for the educational parameter database.
Renders parameter hints at configurable verbosity levels and provides
a reference guide printer.

Depends on Rich for display and the parameter database for content.
"""

from enum import Enum
from typing import Any, Dict, List, Optional

from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich import box

from .amber_parameter_database import Parameter, get_parameters_by_category


class Verbosity(Enum):
    """Educational content verbosity levels."""
    TERSE = "terse"
    COMPACT = "compact"
    VERBOSE = "verbose"


def format_commonality_tag(param: Parameter) -> str:
    """Return Rich markup string for commonly/rarely changed tag."""
    if param.common:
        return "[bold green]\\[commonly changed][/bold green]"
    return "[grey50]\\[rarely changed][/grey50]"


def _format_options_line(param: Parameter, recommended: Any = None) -> str:
    """Format options or range into a compact display string.

    Marks the recommended value with '-->' prefix.
    """
    if recommended is None:
        recommended = param.default

    if param.options:
        parts = []
        for val, desc in param.options.items():
            if val == recommended:
                parts.append(f"[bold cyan]-->{val}={desc}[/bold cyan]")
            else:
                parts.append(f"{val}={desc}")
        return "  ".join(parts)

    # No options dict -- show range or default
    range_parts = []
    if param.min_val is not None:
        range_parts.append(f"min={param.min_val}")
    if param.max_val is not None:
        range_parts.append(f"max={param.max_val}")
    default_str = f"[bold cyan]default={recommended}[/bold cyan]"
    if range_parts:
        return f"{default_str}  ({', '.join(range_parts)})"
    return default_str


def format_parameter_hint(
    param: Parameter,
    verbosity: Verbosity,
    console: Console,
    recommended: Any = None,
) -> None:
    """Print educational content for a parameter at the given verbosity level.

    Terse:   Commonality tag, param name, brief, options with recommendation.
    Compact: Terse + what_it_does explanation.
    Verbose: Compact + manual_notes (when the manual says anything) + related.

    Verbose used to add why_default and when_to_change. Both fields are gone;
    they could only be filled with advice the manual does not give. Most
    parameters now show nothing extra at Verbose beyond `related`, which is the
    honest outcome -- see the module docstring of amber_parameter_database.
    """
    tag = format_commonality_tag(param)
    options_line = _format_options_line(param, recommended)

    # --- Terse: single compact block ---
    console.print(f"\n  {tag} [bold]{param.name}[/bold]: {param.brief}")
    console.print(f"  {options_line}")

    if verbosity == Verbosity.TERSE:
        return

    # --- Compact: add what_it_does ---
    what = param.what_it_does.strip()
    console.print(f"  [grey50]{what}[/grey50]")

    if verbosity == Verbosity.COMPACT:
        return

    # --- Verbose: add whatever the manual itself says ---
    if param.manual_notes:
        console.print(f"  [italic]Amber manual:[/italic] {param.manual_notes}")
    if param.related:
        console.print(f"  [grey50]Related: {', '.join(param.related)}[/grey50]")


def print_reference_guide(
    db: Dict[str, Parameter],
    console: Console,
    category_filter: Optional[str] = None,
) -> None:
    """Print the full parameter reference guide organized by category.

    Always shows full verbose content regardless of the wizard's verbosity setting.

    Args:
        db: The parameter database.
        console: Rich console for output.
        category_filter: Optional category name to show only one category.
    """
    grouped = get_parameters_by_category(db)

    with console.pager(styles=True):
        console.print(Panel(
            "[bold]AMBER MDIN Parameter Reference Guide[/bold]\n"
            "Complete educational reference for all AMBER input parameters.\n"
            "Content sourced from the AMBER Manual (Chapters 22-23).",
            border_style="cyan",
            expand=False,
        ))

        for category, params in grouped.items():
            if category_filter and category != category_filter:
                continue

            table = Table(
                title=f"[bold yellow]{category}[/bold yellow]",
                box=box.ROUNDED,
                expand=False,
                show_lines=True,
            )
            table.add_column("Parameter", style="cyan", width=16)
            table.add_column("Default", style="green", width=10)
            table.add_column("Description", style="white", width=60)
            table.add_column("Commonality", width=18)

            for param in sorted(params, key=lambda p: p.name):
                tag = "commonly changed" if param.common else "rarely changed"
                tag_str = f"[green]{tag}[/green]" if param.common else tag

                # Build description cell
                desc_parts = [f"[bold]{param.brief}[/bold]"]
                desc_parts.append(param.what_it_does.strip())
                if param.manual_notes:
                    desc_parts.append(f"[italic]Amber manual:[/italic] {param.manual_notes}")
                if param.options:
                    opts = "; ".join(f"{k}={v}" for k, v in param.options.items())
                    desc_parts.append(f"Options: {opts}")
                if param.related:
                    desc_parts.append(f"Related: {', '.join(param.related)}")
                description = "\n".join(desc_parts)

                default_str = f"{param.default:g}" if isinstance(param.default, float) else str(param.default)

                table.add_row(param.name, default_str, description, tag_str)

            console.print(table)
            console.print()
