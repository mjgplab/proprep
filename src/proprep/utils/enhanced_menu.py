"""
Enhanced Menu Display Framework

Provides color-coded status indicators and dependency tracking for module menus.
"""

from typing import Dict, List, Tuple, Optional, Any
from enum import Enum


class OptionStatus(Enum):
    """Status of a menu option."""
    AVAILABLE = "available"      # Can be executed now (green)
    READY = "ready"              # Prerequisites met, not yet done (green)
    COMPLETED = "completed"      # Already executed (dim green)
    BLOCKED = "blocked"          # Prerequisites not met (dim)
    RECOMMENDED = "recommended"   # Suggested next step (yellow/bold)


# Leading status glyph + colour per state. Mirrors the top-level menu:
# a colour-coded glyph at the START of the row (not the end), in
# white-robust colours (no ``dim``), each state visually distinct.
_STATE_GLYPH = {
    OptionStatus.AVAILABLE:   ("✓", "green"),
    OptionStatus.READY:       ("✓", "green"),
    OptionStatus.COMPLETED:   ("●", "green"),
    OptionStatus.BLOCKED:     ("○", "grey50"),
    OptionStatus.RECOMMENDED: ("→", "dark_orange"),
}

# Name + note styling, matched to the top-level menu (menu_commands.py:
# _MENU_NAME_STYLE / _MENU_NOTE_STYLE). dark_orange3 (not dark_orange) for
# the note keeps it legible on a white background; see the rationale there.
_OPTION_NAME_STYLE = "bold blue"
_OPTION_NOTE_STYLE = "bold dark_orange3"


def _clean_dependency_text(text: str) -> str:
    """Normalise a legacy ``dependency_text`` such as
    ``"[Need to setup simulations first] ○"`` into a plain reason
    (``"Need to setup simulations first"``) for the ⚠ note line —
    stripping the bracket wrapper and any trailing status glyph the old
    end-of-line format baked in."""
    if not text:
        return ""
    import re
    t = text.strip()
    m = re.match(r"^\[(.*?)\]\s*", t)
    if m:
        t = (m.group(1) + " " + t[m.end():]).strip()
    return t.rstrip(" ○●✓→✗•").strip()


class MenuOption:
    """Represents a single menu option with status and dependencies."""

    def __init__(
        self,
        key: str,
        description: str,
        status: OptionStatus = OptionStatus.AVAILABLE,
        dependency_text: str = "",
        is_separator: bool = False
    ):
        """
        Initialize a menu option.

        Args:
            key: Option key (e.g., "1", "2", "load")
            description: Human-readable description
            status: Current status of the option
            dependency_text: Text explaining dependencies (e.g., "Need to load first")
            is_separator: Whether this is a section separator
        """
        self.key = key
        self.description = description
        self.status = status
        self.dependency_text = dependency_text  # Fixed: was missing 'self.'
        self.is_separator = is_separator

    def status_glyph(self) -> Tuple[str, str]:
        """Return ``(glyph, colour)`` for this option's status — a
        leading, colour-coded indicator like the top-level menu."""
        return _STATE_GLYPH.get(self.status, ("•", "white"))

    def clean_note(self) -> str:
        """Plain-text blocked reason for the ⚠ note line, or ''."""
        return _clean_dependency_text(self.dependency_text)


class EnhancedMenuDisplay:
    """Helper class for displaying menus with color-coded status indicators."""

    def __init__(self, console):
        """
        Initialize the enhanced menu display.

        Args:
            console: Rich Console instance
        """
        self.console = console

    def print_option(self, option: MenuOption):
        """
        Print a single menu option, matching the top-level menu layout:
        a leading colour-coded status glyph, the option name in bold
        blue, and — for blocked options — an indented ⚠ reason line.

        Args:
            option: MenuOption to display
        """
        if option.is_separator:
            self.console.print(f"\n{option.description}")
            return

        from rich.padding import Padding

        glyph, gcolor = option.status_glyph()
        prefix = f"  {option.key}. "
        self.console.print(
            f"{prefix}[{gcolor}]{glyph}[/{gcolor}] "
            f"[{_OPTION_NAME_STYLE}]{option.description}[/{_OPTION_NAME_STYLE}]",
            highlight=False,
        )

        note = option.clean_note()
        if note:
            # Hanging-indented to align under the name (prefix width +2
            # for the status glyph and its trailing space).
            self.console.print(
                Padding(
                    f"[{_OPTION_NOTE_STYLE}]⚠ {note}[/{_OPTION_NOTE_STYLE}]",
                    (0, 0, 0, len(prefix) + 2),
                ),
                highlight=False,
            )

    def print_menu(
        self,
        title: str,
        options: List[MenuOption],
        suggestion: Optional[str] = None
    ):
        """
        Print a complete menu with title, options, and optional suggestion.

        Args:
            title: Menu title
            options: List of MenuOption objects
            suggestion: Optional suggestion text to display at bottom
        """
        # Print title
        self.console.print(f"\n[bold]===== {title} =====[/bold]")

        # Print each option
        for option in options:
            self.print_option(option)

        # Print suggestion if provided
        if suggestion:
            self.console.print(f"\n[bright_blue]→ Suggestion: {suggestion}[/bright_blue]")


def create_dependency_checker(workspace):
    """
    Factory function to create dependency checking functions.

    Args:
        workspace: Current workspace state

    Returns:
        Function that checks if a workspace key exists and has a value
    """
    def has_data(key: str) -> bool:
        """Check if workspace has data for the given key."""
        value = workspace.get(key)
        return value is not None

    return has_data


def determine_option_status(
    workspace,
    completion_key: Optional[str] = None,
    required_keys: Optional[List[str]] = None,
    optional_keys: Optional[List[str]] = None
) -> Tuple[OptionStatus, str]:
    """
    Determine the status of a menu option based on workspace state.

    Args:
        workspace: Current workspace
        completion_key: Key that tracks if this option was completed
        required_keys: List of workspace keys that must exist (AND logic)
        optional_keys: List of workspace keys where at least one must exist (OR logic)

    Returns:
        Tuple of (OptionStatus, dependency_text)
    """
    # Check if already completed
    if completion_key and workspace.get(completion_key) is not None:
        return (OptionStatus.COMPLETED, "")

    # Friendly phrasing for raw workspace keys (e.g. "parm7_file" ->
    # "a topology", structure keys -> "a loaded structure"). Reuses the
    # same map the top-level menu uses so wording stays consistent.
    from proprep.utils.module_registry import _friendly_requirement

    # Check required keys (AND logic)
    if required_keys:
        missing = [key for key in required_keys if workspace.get(key) is None]
        if missing:
            return (OptionStatus.BLOCKED, f"Need {_friendly_requirement(missing)}")

    # Check optional keys (OR logic) — the whole group is one pipe-OR
    # requirement, satisfied if any one key is present.
    if optional_keys:
        has_any = any(workspace.get(key) is not None for key in optional_keys)
        if not has_any:
            reason = _friendly_requirement([" | ".join(optional_keys)])
            return (OptionStatus.BLOCKED, f"Need {reason}")

    # If we got here, prerequisites are met
    if completion_key:
        # Prerequisites met but not done yet
        return (OptionStatus.READY, "")
    else:
        # No completion tracking, always available
        return (OptionStatus.AVAILABLE, "")


def suggest_next_action(options: List[MenuOption]) -> Optional[str]:
    """
    Suggest the next recommended action based on option statuses.

    Args:
        options: List of menu options

    Returns:
        Suggestion text or None
    """
    # Find first available or ready option
    for option in options:
        if option.status in [OptionStatus.AVAILABLE, OptionStatus.READY]:
            if option.key.isdigit():
                return f"Try option {option.key}: {option.description}"

    # Check if everything is done
    all_completed = all(
        opt.status == OptionStatus.COMPLETED or opt.is_separator
        for opt in options
    )
    if all_completed:
        return "All options completed! Continue to next module"

    return None
