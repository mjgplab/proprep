"""
Command History Module

Command history management for MPSA processor commands.
"""

from typing import List, Optional
from .command import Command


class CommandHistory:
    """Manages command execution history and breadcrumbs."""

    def __init__(self):
        self.history: List["Command"] = []
        self.breadcrumbs: List[str] = []
        self.max_history = 100  # Limit history size

    def add_command(self, command: "Command"):
        """Add a command to history."""
        self.history.append(command)
        if len(self.history) > self.max_history:
            self.history.pop(0)

    def add_breadcrumb(self, breadcrumb: str):
        """Add a breadcrumb to the navigation trail."""
        # Deduplicate - don't add the same breadcrumb twice in succession
        if not self.breadcrumbs or self.breadcrumbs[-1] != breadcrumb:
            self.breadcrumbs.append(breadcrumb)

    def pop_breadcrumb(self) -> Optional[str]:
        """Remove and return the last breadcrumb."""
        return self.breadcrumbs.pop() if self.breadcrumbs else None

    def get_breadcrumb_trail(self) -> str:
        """Get formatted breadcrumb trail."""
        if not self.breadcrumbs:
            return "Home"
        return " > ".join(self.breadcrumbs)

    def clear_breadcrumbs(self):
        """Clear all breadcrumbs."""
        self.breadcrumbs.clear()

    def get_recent_commands(self, count: int = 10) -> List[Command]:
        """Get recent commands from history."""
        return self.history[-count:]
