"""
Workspace Management for MPSA

This module provides a centralized workspace for sharing data between modules.
"""

import logging
from typing import Any, Dict, List, Optional, Set, Tuple

from rich.console import Console
from rich.table import Table

from .workspace import Workspace

logger = logging.getLogger(__name__)


class MPSAWorkspace(Workspace):
    """
    Centralized workspace for sharing data between modules

    The Workspace class provides a dictionary-like interface for storing
    and retrieving data shared between modules. It includes features for
    validating data types, logging changes, and displaying workspace contents.
    """

    def __init__(self, debug=False):
        """
        Initialize the workspace

        Args:
            debug: Enable debug logging
        """
        self._data = {"debug": debug}

        # Configure logger
        self.logger = logging.getLogger(__name__)
        if debug:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.INFO)

    def get(self, key: str, default: Any = None) -> Any:
        """
        Get a workspace value with a default

        Args:
            key: Key to retrieve
            default: Default value if key doesn't exist

        Returns:
            The value for the given key, or the default if not found
        """
        return self._data.get(key, default)

    def set(self, key: str, value: Any) -> None:
        """
        Set a workspace value

        Args:
            key: Key to set
            value: Value to store
        """
        old_value = self._data.get(key)
        self._data[key] = value

        if self._data.get("debug", False):
            self._log_update(key, value, old_value)

    def has(self, key: str) -> bool:
        """
        Check if a key exists in the workspace

        Args:
            key: Key to check

        Returns:
            True if key exists, False otherwise
        """
        return key in self._data and key != "debug"

    def update(self, updates: Dict[str, Any]) -> None:
        """
        Update multiple workspace values

        Args:
            updates: Dictionary of key-value pairs to update
        """
        for key, value in updates.items():
            self.set(key, value)

    def dump(self) -> None:
        """
        Dump workspace contents to stdout (simple format)
        """
        print("Workspace Contents:")
        for key, value in sorted(self._data.items()):
            if key == "debug":
                continue

            value_type = type(value).__name__

            # Size info
            size_info = ""
            if hasattr(value, "__len__") and not isinstance(value, (str, bytes)):
                size_info = f"({len(value)} items)"
            elif isinstance(value, str):
                size_info = f"({len(value)} chars)"

            print(f"  {key}: {value_type} {size_info}")

    def dump_rich(self, console: Console, detailed: bool = False) -> None:
        """
        Dump workspace contents to Rich console

        Args:
            console: Rich console for output
            detailed: Whether to show detailed information
        """
        table = Table(title="Workspace Contents")
        table.add_column("Key", style="cyan")
        table.add_column("Type", style="green")
        table.add_column("Size", style="blue")

        for key, value in sorted(self._data.items()):
            if key == "debug":
                continue

            value_type = type(value).__name__

            # Size info
            size_info = ""
            if hasattr(value, "__len__") and not isinstance(value, (str, bytes)):
                size_info = str(len(value))
            elif isinstance(value, str):
                size_info = f"{len(value)} chars"

            table.add_row(key, value_type, size_info)

        console.print(table)

    def clear(self) -> None:
        """Clear all workspace data except debug flag"""
        debug = self._data.get("debug", False)
        self._data = {"debug": debug}

    def get_keys(self) -> List[str]:
        """Get all keys in the workspace"""
        return list(self._data.keys())

    def items(self) -> List[Tuple[str, Any]]:
        """
        Get all key-value pairs in the workspace

        Returns:
            List of tuples (key, value)
        """
        return list([item for item in self._data.items() if item[0] != "debug"])

    def __iter__(self):
        """
        Allow iteration over key-value pairs

        Returns:
            Iterator over (key, value) pairs
        """
        return iter(self.items())

    def __setitem__(self, key: str, value: Any) -> None:
        """
        Dictionary-style item assignment

        Args:
            key: Key to set
            value: Value to store
        """
        logger.warning("This method is deprecated. Use set() instead.")
        self.set(key, value)

    def __getitem__(self, key: str) -> Any:
        """
        Dictionary-style item access

        Args:
            key: Key to retrieve

        Returns:
            The value for the given key

        Raises:
            KeyError: If key doesn't exist
        """
        logger.warning("This method is deprecated. Use get() instead.")
        if key not in self._data:
            raise KeyError(f"Key '{key}' not found in workspace")
        return self._data[key]

    def __contains__(self, key: str) -> bool:
        """
        Check if key exists in workspace

        Args:
            key: Key to check

        Returns:
            True if key exists, False otherwise
        """
        logger.warning("This method is deprecated. Use has() instead.")
        return key in self._data

    def _log_update(self, key, value, old_value):
        """
        Log a workspace update

        Args:
            key: Key being updated
            value: New value
            old_value: Previous value
        """
        value_type = type(value).__name__
        value_info = f"{value_type}"
        if hasattr(value, "__len__") and not isinstance(value, (str, bytes)):
            value_info = f"{value_type}({len(value)})"

        if old_value is None:
            self.logger.debug(f"Added new value for '{key}': {value_info}")
        else:
            old_type = type(old_value).__name__
            old_info = f"{old_type}"
            if hasattr(old_value, "__len__") and not isinstance(
                old_value, (str, bytes)
            ):
                old_info = f"{old_type}({len(old_value)})"

            self.logger.debug(f"Updated '{key}' from {old_info} to {value_info}")
