"""
Command pattern implementation for MPSA processor.

This package provides command classes for all processor operations,
enabling better organization, undo functionality, and breadcrumb tracking.
"""

from .command import Command
from .processor_command import MenuCommand, ProcessorCommand, ModuleCommand
from .module_commands import RunModuleCommand, RunModuleMenuCommand
from .processor import Processor
from .pdbprocessor import PDBProcessor
from .workspace_commands import (
    ShowWorkspaceStatusCommand,
    ShowDetailedWorkspaceCommand,
    SaveWorkspaceCommand,
    LoadWorkspaceCommand,
    ResetWorkspaceCommand,
    ToggleDebugCommand,
    ShowWorkspaceHistoryCommand,
    SaveCommandHistoryCommand,
    LoadCommandHistoryCommand,
    ReplayCommandHistoryCommand,
)
from .menu_commands import MainMenuCommand, WorkspaceOptionsMenuCommand

__all__ = [
    "Command",
    "MenuCommand",
    "ProcessorCommand",
    "PDBProcessor",
    "RunModuleCommand",
    "RunModuleMenuCommand",
    "ShowWorkspaceStatusCommand",
    "ShowDetailedWorkspaceCommand",
    "SaveWorkspaceCommand",
    "LoadWorkspaceCommand",
    "ResetWorkspaceCommand",
    "ToggleDebugCommand",
    "ShowWorkspaceHistoryCommand",
    "SaveCommandHistoryCommand",
    "LoadCommandHistoryCommand",
    "ReplayCommandHistoryCommand",
    "ResetCommandHistoryCommand",
    "MainMenuCommand",
    "WorkspaceOptionsMenuCommand",
]
