"""
Settings Manager for ProPrep

Manages user preferences and settings stored in ~/.proprep/settings.json
"""

import json
import os
from pathlib import Path
from typing import Literal, Optional
from datetime import datetime


MenuMode = Literal["workflow", "full-menu", "ask"]
MenuLayout = Literal["list", "grid"]


class SettingsManager:
    """Manages user settings and preferences."""

    def __init__(self, settings_dir: Optional[Path] = None):
        """
        Initialize settings manager.

        Args:
            settings_dir: Optional custom directory for settings.
                         Defaults to ~/.proprep
        """
        if settings_dir is None:
            self.settings_dir = Path.home() / ".proprep"
        else:
            self.settings_dir = Path(settings_dir)

        self.settings_file = self.settings_dir / "settings.json"
        self._ensure_settings_dir()
        self._settings = self._load_settings()

    def _ensure_settings_dir(self):
        """Ensure settings directory exists."""
        self.settings_dir.mkdir(parents=True, exist_ok=True)

    def _load_settings(self) -> dict:
        """Load settings from file or create default settings."""
        if self.settings_file.exists():
            try:
                with open(self.settings_file, 'r') as f:
                    return json.load(f)
            except (json.JSONDecodeError, IOError):
                # If file is corrupted, use defaults
                return self._get_default_settings()
        else:
            return self._get_default_settings()

    def _get_default_settings(self) -> dict:
        """Get default settings."""
        return {
            "default_menu_mode": "ask",
            "menu_layout": "list",
            "github_token": None,
            "github_repo": None,
            "last_updated": datetime.now().isoformat()
        }

    def _save_settings(self):
        """Save settings to file."""
        self._settings["last_updated"] = datetime.now().isoformat()
        try:
            with open(self.settings_file, 'w') as f:
                json.dump(self._settings, f, indent=2)
        except IOError:
            # Fail silently - settings won't persist but app continues
            return

        # This file can hold a GitHub personal access token, so restrict it to
        # owner read/write only (0o600) to match the "stored securely" promise
        # shown in the feedback setup. Best-effort: the save already succeeded,
        # so a chmod failure on an exotic filesystem must not be treated as
        # fatal.
        try:
            os.chmod(self.settings_file, 0o600)
        except OSError:
            pass

    def get_menu_mode(self) -> MenuMode:
        """
        Get the default menu mode preference.

        Returns:
            Menu mode: "workflow", "full-menu", or "ask"
        """
        mode = self._settings.get("default_menu_mode", "ask")
        if mode not in ["workflow", "full-menu", "ask"]:
            return "ask"
        return mode

    def set_menu_mode(self, mode: MenuMode):
        """
        Set the default menu mode preference.

        Args:
            mode: Menu mode to set ("workflow", "full-menu", or "ask")
        """
        if mode not in ["workflow", "full-menu", "ask"]:
            raise ValueError(f"Invalid menu mode: {mode}")

        self._settings["default_menu_mode"] = mode
        self._save_settings()

    def get_menu_layout(self) -> MenuLayout:
        """
        Get the full-menu layout preference.

        Returns:
            Layout: "list" (single column, default) or "grid"
            (2-column panel grid, good for screenshots/figures)
        """
        layout = self._settings.get("menu_layout", "list")
        if layout not in ["list", "grid"]:
            return "list"
        return layout

    def set_menu_layout(self, layout: MenuLayout):
        """
        Set the full-menu layout preference.

        Args:
            layout: Layout to set ("list" or "grid")
        """
        if layout not in ["list", "grid"]:
            raise ValueError(f"Invalid menu layout: {layout}")

        self._settings["menu_layout"] = layout
        self._save_settings()

    def reset_to_defaults(self):
        """Reset all settings to defaults."""
        self._settings = self._get_default_settings()
        self._save_settings()

    def get_github_token(self) -> str:
        """
        Get the GitHub Personal Access Token for feedback submission.

        Returns:
            GitHub PAT or None if not configured
        """
        return self._settings.get("github_token")

    def set_github_token(self, token: str):
        """
        Set the GitHub Personal Access Token.

        Args:
            token: GitHub PAT with repo:issues scope
        """
        self._settings["github_token"] = token
        self._save_settings()

    def get_github_repo(self) -> str:
        """
        Get the GitHub repository for feedback submission.

        Returns:
            Repository in format "owner/repo" or None if not configured
        """
        return self._settings.get("github_repo")

    def set_github_repo(self, repo: str):
        """
        Set the GitHub repository for feedback submission.

        Args:
            repo: Repository in format "owner/repo" (e.g., "mjgp/proprep")
        """
        self._settings["github_repo"] = repo
        self._save_settings()
