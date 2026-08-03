"""
Workspace-related command implementations.

This module contains commands for workspace operations like
saving, loading, resetting, and displaying workspace status.
"""

import json
import os
from typing import Any, Dict, List, Optional

import yaml
from proprep.utils.prompts import prompt_with_context, confirm_with_context

from .processor_command import WorkspaceCommand, ModuleActionCommand
from .processor import Processor
from proprep.utils.prompts import prompt_with_context, confirm_with_context


class ShowWorkspaceStatusCommand(WorkspaceCommand):
    """Command to display workspace status."""

    def __init__(self, processor: Processor, detailed: bool = False):
        super().__init__(processor, "Show workspace status")
        self.detailed = detailed

    def execute(self) -> None:
        """Display workspace status."""
        self.console.print("\n[bold]===== Workspace Status =====[/bold]")
        self.workspace.dump_rich(self.console, detailed=self.detailed)


class ShowDetailedWorkspaceCommand(ShowWorkspaceStatusCommand):
    """Command to display detailed workspace information."""

    def __init__(self, processor: Processor):
        super().__init__(processor, detailed=True)
        self.description = "Show detailed workspace information"


class SaveWorkspaceCommand(WorkspaceCommand):
    """Command to save workspace to a file."""

    def __init__(self, processor: Processor, filename: str = None):
        super().__init__(processor, "Save workspace")
        self.filename = filename

    def execute(self) -> bool:
        """Save workspace to file."""
        if not self.filename:
            self.filename = prompt_with_context(
                self.processor,
                "Enter filename to save workspace",
                default="workspace.json",
                module="Workspace",
                description="Save workspace filename"
            )
        return self.processor.save_workspace(self.filename)


class LoadWorkspaceCommand(WorkspaceCommand):
    """Command to load workspace from a file."""

    def __init__(self, processor: Processor, filename: str = None):
        super().__init__(processor, "Load workspace")
        self.filename = filename

    def execute(self) -> bool:
        """Load workspace from file."""
        if not self.filename:
            self.filename = prompt_with_context(
                self.processor,
                "Enter filename to load workspace",
                default="workspace.json",
                module="Workspace",
                description="Load workspace filename"
            )

        if not os.path.exists(self.filename):
            self.console.print(f"[red]File {self.filename} not found[/red]")
            return False

        return self.processor.load_workspace(self.filename)


class ResetWorkspaceCommand(WorkspaceCommand):
    """Command to reset workspace to initial state."""

    def __init__(self, processor: Processor):
        super().__init__(processor, "Reset workspace")

    def execute(self) -> bool:
        """Reset workspace."""
        return self.processor.reset_workspace()


class ToggleDebugCommand(WorkspaceCommand):
    """Command to toggle debug mode."""

    def __init__(self, processor: Processor):
        super().__init__(processor, "Toggle debug mode")

    def execute(self) -> None:
        """Toggle debug mode."""
        workspace = self.workspace
        debug_enabled = workspace.get("debug", False)
        new_debug_state = not debug_enabled
        workspace.set("debug", new_debug_state)

        status = "enabled" if new_debug_state else "disabled"
        self.console.print(f"[green]Debug mode {status}[/green]")


class ClearWorkspaceCommand(WorkspaceCommand):
    """Command to clear all workspace data."""

    def __init__(self, processor: Processor):
        super().__init__(processor, "Clear workspace")

    def execute(self) -> bool:
        """Clear workspace data."""
        if confirm_with_context(
            self.processor,
            "[yellow]Are you sure you want to clear all workspace data?[/yellow]",
            default=False,
            module="Workspace",
            description="Confirm clear workspace"
        ):
            workspace = self.workspace
            debug_setting = workspace.get("debug", False)

            workspace.clear()
            workspace.set("debug", debug_setting)

            self.console.print("[green]Workspace data cleared[/green]")
            return True
        return False


class ShowBreadcrumbsCommand(WorkspaceCommand):
    """Command to show current navigation breadcrumbs."""

    def __init__(self, processor: Processor):
        super().__init__(processor, "Show navigation breadcrumbs")

    def execute(self) -> None:
        """Show current breadcrumb trail."""
        trail = self.processor.command_history.get_breadcrumb_trail()
        self.console.print(f"\n[bold]Current Location:[/bold] {trail}")
        self.console.print(f"\n[bold]Current Location:[/bold] {trail}")


class ExportWorkspaceCommand(WorkspaceCommand):
    """Command to export workspace in different formats."""

    def __init__(self, processor: Processor, format_type: str = "json"):
        super().__init__(processor, f"Export workspace as {format_type}")
        self.format_type = format_type

    def execute(self) -> bool:
        """Export workspace in specified format."""
        filename = prompt_with_context(
            self.processor,
            f"Enter filename for {self.format_type} export",
            default=f"workspace_export.{self.format_type}",
            module="Workspace",
            description="Export workspace filename"
        )

        try:
            workspace_data = dict(self.workspace.items())

            if self.format_type == "json":
                return self._export_json(filename, workspace_data)
            elif self.format_type == "yaml":
                return self._export_yaml(filename, workspace_data)
            else:
                self.console.print(f"[red]Unsupported format: {self.format_type}[/red]")
                return False

        except Exception as e:
            self.console.print(f"[red]Export failed: {str(e)}[/red]")
            return False

    def _export_json(self, filename: str, data: dict) -> bool:
        """Export as JSON."""
        with open(filename, "w") as f:
            json.dump(data, f, indent=2, default=str)
        self.console.print(f"[green]Workspace exported to {filename}[/green]")
        return True

    def _export_yaml(self, filename: str, data: dict) -> bool:
        """Export as YAML."""
        try:
            with open(filename, "w") as f:
                yaml.dump(data, f, default_flow_style=False)
            self.console.print(f"[green]Workspace exported to {filename}[/green]")
            return True
        except ImportError:
            self.console.print("[red]YAML export requires PyYAML package[/red]")
            return False


class ShowWorkspaceHistoryCommand(WorkspaceCommand):
    """Command to display command execution history."""

    def __init__(self, processor: Processor, count: int = 10):
        super().__init__(processor, "Show command history")
        self.count = count

    def execute(self) -> None:
        """Display recent command history."""
        self.processor.show_command_history(self.count)


class SaveCommandHistoryCommand(WorkspaceCommand):
    """Command to save command history to a JSON file."""

    def __init__(self, processor: Processor, filename: str = None):
        super().__init__(processor, "Save command history")
        self.filename = filename

    def execute(self) -> bool:
        """Save command history to JSON file."""
        if not self.filename:
            self.filename = prompt_with_context(
                self.processor,
                "Enter filename to save command history",
                default="history.json",
                module="Workspace",
                description="Save command history filename"
            )

        try:
            history_data = []
            for cmd in self.processor.command_history.history:
                cmd_data = self._serialize_command(cmd)
                history_data.append(cmd_data)

            export_data = {
                "command_history": history_data,
                "breadcrumb_trail": self.processor.command_history.breadcrumbs.copy(),
                "export_timestamp": (
                    self.processor.command_history.history[-1].timestamp.isoformat()
                    if self.processor.command_history.history
                    else None
                ),
            }

            with open(self.filename, "w") as f:
                json.dump(export_data, f, indent=2)

            self.console.print(
                f"[green]Command history saved to {self.filename}[/green]"
            )
            return True

        except Exception as e:
            self.console.print(f"[red]Failed to save command history: {str(e)}[/red]")
            return False

    def _serialize_command(self, cmd) -> Dict[str, Any]:
        """Create a serializable representation of a command."""
        cmd_data = {
            "description": cmd.description,
            "timestamp": cmd.timestamp.isoformat(),
            "executed": cmd.executed,
            "error": str(cmd.error) if cmd.error else None,
            "class_name": cmd.__class__.__name__,
            "breadcrumb": (
                cmd.get_breadcrumb()
                if callable(getattr(cmd, "get_breadcrumb", None))
                else cmd.description
            ),
        }

        # Add ModuleActionCommand specific attributes
        if isinstance(cmd, ModuleActionCommand):
            cmd_data.update(
                {
                    "module_name": cmd.module_name,
                    "action_name": cmd._action_to_run,
                    "args": cmd.args,
                    "captured_state": self._filter_serializable_state(
                        cmd.captured_state
                    ),
                }
            )

        # Add other command-specific attributes
        for attr in ["filename", "format_type", "count", "option_key"]:
            value = getattr(cmd, attr, None)
            if value is not None:
                cmd_data[attr] = value

        return cmd_data

    def _filter_serializable_state(self, captured_state):
        """Filter out non-serializable objects from captured state."""
        if not captured_state:
            return captured_state

        filtered_state = {}

        for key, value in captured_state.items():
            try:
                # Test if the value is JSON serializable
                json.dumps(value)
                filtered_state[key] = value
            except (TypeError, ValueError):
                # If not serializable, store a description instead
                if hasattr(value, "__class__"):
                    filtered_state[key] = (
                        f"<Non-serializable {value.__class__.__module__}.{value.__class__.__name__} object>"
                    )
                else:
                    filtered_state[key] = (
                        f"<Non-serializable {type(value).__name__} object>"
                    )

        return filtered_state


class LoadCommandHistoryCommand(WorkspaceCommand):
    """Command to load command history from a JSON file."""

    def __init__(self, processor: Processor, filename: str = None):
        super().__init__(processor, "Load command history")
        self.filename = filename

    def execute(self) -> bool:
        """Load command history from JSON file."""
        if not self.filename:
            self.filename = prompt_with_context(
                self.processor,
                "Enter filename to load command history",
                default="history.json",
                module="Workspace",
                description="Load command history filename"
            )

        if not os.path.exists(self.filename):
            self.console.print(f"[red]File {self.filename} not found[/red]")
            return False

        try:
            with open(self.filename, "r") as f:
                data = json.load(f)

            if "command_history" not in data:
                self.console.print(f"[red]Invalid history file format[/red]")
                return False

            if self.processor.command_history.history:
                if not confirm_with_context(
                    self.processor,
                    "[yellow]This will replace current command history. Continue?[/yellow]",
                    default=False,
                    module="Workspace",
                    description="Confirm replace command history"
                ):
                    return False

            self.processor.command_history.history.clear()
            self.processor.command_history.breadcrumbs.clear()

            if "breadcrumb_trail" in data:
                self.processor.command_history.breadcrumbs.extend(
                    data["breadcrumb_trail"]
                )

            self.console.print(
                f"[green]Command history loaded from {self.filename}[/green]"
            )
            self.console.print(
                f"[bright_blue]Loaded {len(data['command_history'])} command records[/bright_blue]"
            )

            if data.get("export_timestamp"):
                self.console.print(
                    f"[bright_blue]History exported at: {data['export_timestamp']}[/bright_blue]"
                )

            return True

        except Exception as e:
            self.console.print(f"[red]Failed to load command history: {str(e)}[/red]")
            return False


class ReplayCommandHistoryCommand(WorkspaceCommand):
    """Command to replay command history with full ModuleActionCommand support."""

    def __init__(self, processor: Processor, filename: str = None):
        super().__init__(processor, "Replay command history")
        self.filename = filename
        self.history_data = None

    def execute(self) -> None:
        """Replay command history with full execution support."""
        if not self._load_history_data():
            return

        self._show_replay_menu()

    def _load_history_data(self) -> bool:
        """Load history data from file."""
        if not self.filename:
            self.filename = prompt_with_context(
                self.processor,
                "Enter filename to replay command history",
                default="history.json",
                module="Workspace",
                description="Replay command history filename"
            )

        if not os.path.exists(self.filename):
            self.console.print(f"[red]File {self.filename} not found[/red]")
            return False

        try:
            with open(self.filename, "r") as f:
                self.history_data = json.load(f)

            if "command_history" not in self.history_data:
                self.console.print(f"[red]Invalid history file format[/red]")
                return False

            return True

        except Exception as e:
            self.console.print(f"[red]Failed to load history file: {str(e)}[/red]")
            return False

    def _show_replay_menu(self) -> None:
        """Show the command history replay menu."""
        commands = self.history_data["command_history"]

        if not commands:
            self.console.print("[yellow]No commands found in history[/yellow]")
            return

        while True:
            self.console.print(f"\n[bold]===== Command History Replay =====[/bold]")
            self.console.print(f"[bright_blue]File: {self.filename}[/bright_blue]")
            self.console.print(f"[bright_blue]Commands: {len(commands)}[/bright_blue]")

            if self.history_data.get("export_timestamp"):
                self.console.print(
                    f"[bright_blue]Exported: {self.history_data['export_timestamp']}[/bright_blue]"
                )

            self.console.print("\n[bold]Options:[/bold]")
            self.console.print("1. View all commands", highlight=False)
            self.console.print("2. Replay all functional commands", highlight=False)
            self.console.print("3. Restore workspace state", highlight=False)
            self.console.print("4. Restore breadcrumb trail", highlight=False)
            self.console.print("5. Back to workspace menu", highlight=False)

            choice = prompt_with_context(
                self.processor,
                "Enter your choice",
                choices=["1", "2", "3", "4", "5"],
                default="1",
                module="Workspace",
                description="Replay menu choice",
                options_map={"1": "View all", "2": "Replay functional", "3": "Restore state", "4": "Restore breadcrumbs", "5": "Back"}
            )

            if choice == "1":
                self._view_all_commands(commands)
            elif choice == "2":
                self._replay_functional_commands(commands)
            elif choice == "3":
                self._restore_workspace_state(commands)
            elif choice == "4":
                self._restore_breadcrumbs()
            elif choice == "5":
                break

    def _view_all_commands(self, commands: List[Dict]) -> None:
        """View all commands in chronological order."""
        self.console.print(f"\n[bold]===== All Commands ({len(commands)}) =====[/bold]")

        for i, cmd in enumerate(commands, 1):
            status = "✓" if cmd["executed"] else "✗"
            error_info = f" [red](Error: {cmd['error']})[/red]" if cmd["error"] else ""
            module_info = (
                f" [{cmd.get('module_name', 'General')}]"
                if cmd.get("module_name")
                else ""
            )
            self.console.print(
                f"{i:3d}. {status} [{cmd['timestamp'][:19]}]{module_info} {cmd['description']}{error_info}"
            )

        prompt_with_context(
            self.processor,
            "\nPress Enter to continue",
            default="",
            module="Workspace",
            description="Pause after command history view",
        )

    def _replay_functional_commands(self, commands: List[Dict]) -> None:
        """Replay all functional commands (module actions + workspace operations)."""
        functional_commands = [
            cmd for cmd in commands if self._is_functional_command(cmd)
        ]

        if not functional_commands:
            self.console.print(
                "[yellow]No functional commands found to replay[/yellow]"
            )
            return

        self.console.print(
            f"\n[bold]===== Functional Commands ({len(functional_commands)}) =====[/bold]"
        )
        for i, cmd in enumerate(functional_commands, 1):
            status = "✓" if cmd["executed"] else "✗"
            module_info = (
                f" [{cmd.get('module_name', 'General')}]"
                if cmd.get("module_name")
                else ""
            )
            self.console.print(f"{i:3d}. {status}{module_info} {cmd['description']}")

        if confirm_with_context(
            self.processor,
            f"[yellow]Execute {len(functional_commands)} functional command(s)?[/yellow]",
            default=False,
            module="Workspace",
            description="Confirm execute functional commands"
        ):
            self._execute_commands_batch(functional_commands)

    def _replay_all_non_navigation_commands(self, commands: List[Dict]) -> None:
        """Replay all non-navigation commands (excludes pure navigation but includes functional commands)."""
        non_navigation_commands = [
            cmd for cmd in commands if not self._is_pure_navigation_command(cmd)
        ]

        if not non_navigation_commands:
            self.console.print(
                "[yellow]No non-navigation commands found to replay[/yellow]"
            )
            return

        self.console.print(
            f"\n[bold]===== Non-Navigation Commands ({len(non_navigation_commands)}) =====[/bold]"
        )
        for i, cmd in enumerate(non_navigation_commands, 1):
            status = "✓" if cmd["executed"] else "✗"
            module_info = (
                f" [{cmd.get('module_name', 'General')}]"
                if cmd.get("module_name")
                else ""
            )
            self.console.print(f"{i:3d}. {status}{module_info} {cmd['description']}")

        if confirm_with_context(
            self.processor,
            f"[yellow]Execute {len(non_navigation_commands)} non-navigation command(s)?[/yellow]",
            default=False,
            module="Workspace",
            description="Confirm execute non-navigation commands"
        ):
            self._execute_commands_batch(non_navigation_commands)

    def _restore_workspace_state(self, commands: List[Dict]) -> None:
        """Restore workspace state by replaying ModuleActionCommands with captured state."""
        stateful_commands = [
            cmd
            for cmd in commands
            if cmd.get("class_name") == "ModuleActionCommand"
            and cmd.get("captured_state")
        ]

        if not stateful_commands:
            self.console.print(
                "[yellow]No commands with captured workspace state found[/yellow]"
            )
            return

        self.console.print(
            f"\n[bold]===== Stateful Commands ({len(stateful_commands)}) =====[/bold]"
        )
        for i, cmd in enumerate(stateful_commands, 1):
            module_name = cmd.get("module_name", "Unknown")
            action_name = cmd.get("action_name", "Unknown")
            state_keys = list(cmd.get("captured_state", {}).keys())
            self.console.print(
                f"{i:3d}. {module_name}.{action_name} (State: {', '.join(state_keys)})"
            )

        if confirm_with_context(
            self.processor,
            "[yellow]Restore workspace state from these commands?[/yellow]",
            default=False,
            module="Workspace",
            description="Confirm restore workspace state"
        ):
            self._restore_workspace_from_commands(stateful_commands)

    def _restore_breadcrumbs(self) -> None:
        """Restore the breadcrumb trail from the history."""
        if "breadcrumb_trail" not in self.history_data:
            self.console.print("[yellow]No breadcrumb trail found in history[/yellow]")
            return

        breadcrumbs = self.history_data["breadcrumb_trail"]

        if not breadcrumbs:
            self.console.print("[yellow]No breadcrumbs to restore[/yellow]")
            return

        if confirm_with_context(
            self.processor,
            f"[yellow]Restore breadcrumb trail: {' > '.join(breadcrumbs)}?[/yellow]",
            default=True,
            module="Workspace",
            description="Confirm restore breadcrumb trail"
        ):
            self.processor.command_history.breadcrumbs.clear()
            self.processor.command_history.breadcrumbs.extend(breadcrumbs)
            self.console.print("[green]Breadcrumb trail restored[/green]")

    def _restore_workspace_from_commands(self, commands_list: List[Dict]) -> None:
        """Restore workspace state from captured states in commands."""
        restored_keys = set()

        for cmd_data in commands_list:
            captured_state = cmd_data.get("captured_state", {})
            if captured_state:
                workspace = self.processor._get_workspace()
                for key, value in captured_state.items():
                    if value is not None:
                        workspace.set(key, value)
                        restored_keys.add(key)

        if restored_keys:
            self.console.print(
                f"[green]Restored workspace keys: {', '.join(sorted(restored_keys))}[/green]"
            )
        else:
            self.console.print("[yellow]No workspace state to restore[/yellow]")

    def _is_functional_command(self, cmd_data: Dict) -> bool:
        """Check if a command performs functional work."""
        class_name = cmd_data.get("class_name", "")

        # Known non-functional/navigation commands
        navigation_classes = [
            "MainMenuCommand",
            "WorkspaceOptionsMenuCommand",
            "SaveCommandHistoryCommand",
            "LoadCommandHistoryCommand",
            "ReplayCommandHistoryCommand",
            "ResetCommandHistoryCommand",
        ]

        # Known functional workspace commands
        functional_workspace_classes = [
            "ModuleActionCommand",
            "SaveWorkspaceCommand",
            "LoadWorkspaceCommand",
            "ExportWorkspaceCommand",
            "ShowWorkspaceStatusCommand",
            "ShowDetailedWorkspaceCommand",
            "ToggleDebugCommand",
            "ClearWorkspaceCommand",
        ]

        # OptionCommand handling
        if class_name == "OptionCommand":
            # OptionCommand with module_name is functional
            return bool(cmd_data.get("module_name"))

        # Navigation commands are not functional
        if class_name in navigation_classes:
            return False

        # Known functional workspace commands
        if class_name in functional_workspace_classes:
            return True

        # Generic approach: Any command with module_name and action_name is treated as functional
        # This covers all custom module commands (DetectMetalSitesCommand, ImportMetalSitesCommand, etc.)
        if cmd_data.get("module_name") and cmd_data.get("action_name"):
            return True

        # Default to non-functional for unknown commands
        return False

    def _is_pure_navigation_command(self, cmd_data: Dict) -> bool:
        """Check if a command is pure navigation (should be skipped during replay)."""
        class_name = cmd_data.get("class_name", "")

        navigation_classes = [
            "MainMenuCommand",
            "WorkspaceOptionsMenuCommand",
            "SaveCommandHistoryCommand",
            "LoadCommandHistoryCommand",
            "ReplayCommandHistoryCommand",
            "ResetCommandHistoryCommand",
        ]

        # OptionCommand without module_name is navigation
        if class_name == "OptionCommand" and not cmd_data.get("module_name"):
            return True

        return class_name in navigation_classes

    def _execute_commands_batch(self, commands_list: List[Dict]) -> None:
        """Execute a batch of commands with error handling."""
        if not commands_list:
            self.console.print("[yellow]No commands to execute[/yellow]")
            return

        executed_count = 0
        failed_count = 0

        self.console.print(
            f"\n[bold]Starting batch execution of {len(commands_list)} commands...[/bold]"
        )

        for i, cmd_data in enumerate(commands_list, 1):
            self.console.print(
                f"\n[bright_blue]({i}/{len(commands_list)}) {cmd_data['description']}[/bright_blue]"
            )

            try:
                command_obj = self._recreate_command_from_data(cmd_data)
                if command_obj:
                    result = command_obj.execute()
                    self.console.print(f"[green]✓ Success[/green]")
                    executed_count += 1
                else:
                    self.console.print(
                        f"[yellow]⚠ Skipped (navigation command)[/yellow]"
                    )

            except Exception as e:
                self.console.print(f"[red]✗ Failed: {str(e)}[/red]")
                failed_count += 1

                if i < len(commands_list):
                    if not confirm_with_context(
                        self.processor,
                        "Continue with remaining commands?",
                        default=True,
                        module="Workspace",
                        description="Confirm continue batch execution"
                    ):
                        break

        self.console.print(f"\n[bold]Batch execution complete:[/bold]")
        self.console.print(f"[green]Executed: {executed_count}[/green]")
        self.console.print(f"[red]Failed: {failed_count}[/red]")

    def _recreate_command_from_data(self, cmd_data: Dict) -> Optional[WorkspaceCommand]:
        """Recreate command objects from history data."""
        class_name = cmd_data.get("class_name", "")

        # Skip navigation commands
        if self._is_pure_navigation_command(cmd_data):
            return None

        # Handle workspace commands
        if class_name == "ShowWorkspaceStatusCommand":
            return ShowWorkspaceStatusCommand(self.processor)
        elif class_name == "ShowDetailedWorkspaceCommand":
            return ShowDetailedWorkspaceCommand(self.processor)
        elif class_name == "SaveWorkspaceCommand":
            return SaveWorkspaceCommand(self.processor, cmd_data.get("filename"))
        elif class_name == "LoadWorkspaceCommand":
            return LoadWorkspaceCommand(self.processor, cmd_data.get("filename"))
        elif class_name == "ExportWorkspaceCommand":
            return ExportWorkspaceCommand(
                self.processor, cmd_data.get("format_type", "json")
            )
        elif class_name == "ShowWorkspaceHistoryCommand":
            return ShowWorkspaceHistoryCommand(
                self.processor, cmd_data.get("count", 10)
            )
        elif class_name == "OptionCommand" and cmd_data.get("module_name"):
            # Convert OptionCommand with module to RunModuleCommand
            from .module_commands import RunModuleCommand

            return RunModuleCommand(self.processor, cmd_data.get("module_name"))
        elif class_name == "ModuleActionCommand":
            return ModuleActionCommand(
                self.processor,
                cmd_data.get("module_name"),
                cmd_data.get("action_name"),
                cmd_data.get("args", {}),
            )
        elif cmd_data.get("module_name") and cmd_data.get("action_name"):
            # Generic approach: Any command with module_name and action_name
            # is treated as a ModuleActionCommand equivalent
            # This handles all custom module commands (DetectMetalSitesCommand, ImportMetalSitesCommand, etc.)
            return ModuleActionCommand(
                self.processor,
                cmd_data.get("module_name"),
                cmd_data.get("action_name"),
                cmd_data.get("args", {}),
            )

        return None


class ResetCommandHistoryCommand(WorkspaceCommand):
    """Command to reset/clear command history."""

    def __init__(self, processor: Processor):
        super().__init__(processor, "Reset command history")

    def execute(self) -> bool:
        """Reset command history to empty state."""
        if not self.processor.command_history.history:
            self.console.print("[yellow]Command history is already empty[/yellow]")
            return True

        current_count = len(self.processor.command_history.history)
        breadcrumb_count = len(self.processor.command_history.breadcrumbs)

        self.console.print(
            f"[yellow]Current command history: {current_count} commands[/yellow]"
        )
        if breadcrumb_count > 0:
            self.console.print(
                f"[yellow]Current breadcrumb trail: {breadcrumb_count} items[/yellow]"
            )

        if confirm_with_context(
            self.processor,
            "[yellow]Are you sure you want to clear all command history and breadcrumbs?[/yellow]",
            default=False,
            module="Workspace",
            description="Confirm clear command history"
        ):
            self.processor.command_history.history.clear()
            self.processor.command_history.breadcrumbs.clear()

            self.console.print("[green]Command history and breadcrumbs cleared[/green]")
            return True
        else:
            self.console.print("[bright_blue]Command history reset cancelled[/bright_blue]")
            return False
