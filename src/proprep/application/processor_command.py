"""
Base command classes for the MPSA processor command pattern implementation.

This module defines the core command interfaces and base implementations
that provide breadcrumb functionality, and common command utilities.
"""

from abc import abstractmethod
from typing import Any, Dict, Optional, TYPE_CHECKING

from .command import Command

if TYPE_CHECKING:
    from .processor import Processor
    from proprep.utils.workspace import Workspace


class ProcessorCommand(Command):
    """
    Base class for commands that operate on the processor.

    Provides common functionality for commands that need access to
    processor state, workspace, and console output.
    """

    def __init__(self, processor: "Processor", description: str = ""):
        super().__init__(processor, description)

    @property
    def console(self):
        """Get console from processor."""
        return self.processor.console

    @property
    def workspace(self):
        """Get workspace from processor."""
        return self.processor._get_workspace()

    def execute_with_error_handling(self) -> Any:
        """Execute command with error handling and history tracking."""
        try:
            if not self.can_execute():
                raise RuntimeError(f"Cannot execute command: {self.description}")

            if self.processor.command_history:
                self.processor.command_history.add_command(self)

            # Set flag to prevent duplicate history recording in execute()
            self._in_error_handling = True
            try:
                self.result = self.execute()
            finally:
                self._in_error_handling = False

            self.executed = True
            return self.result

        except Exception as e:
            self.error = e
            self.console.print(
                f"[red]Error executing {self.description}: {str(e)}[/red]"
            )
            raise


class MenuCommand(ProcessorCommand):
    """
    Base class for menu-related commands.

    Handles menu navigation, breadcrumb management, and common menu operations.
    Enhanced with navigation escape hatch system.
    """

    def __init__(self, processor: "Processor", menu_name: str, description: str = "", 
                 menu_category: str = None):
        super().__init__(processor, description or f"Show {menu_name} menu")
        self.menu_name = menu_name
        self.menu_category = menu_category  # e.g., "Structure Preparation", "Energetic Evaluation"
        self.options: Dict[str, Command] = {}
        self.navigation_options: Dict[str, str] = {}

    def add_option(self, key: str, command: Command):
        """Add a menu option."""
        self.options[key] = command

    def get_breadcrumb(self) -> str:
        """Get breadcrumb for this menu."""
        return self.menu_name

    def show_breadcrumbs(self):
        """Display current breadcrumb trail."""
        trail = self.processor.command_history.get_breadcrumb_trail()
        self.console.print(f"[grey50]Navigation: {trail}[/grey50]\n")

    def enter_menu(self):
        """Enter this menu (add breadcrumb)."""
        self.processor.command_history.add_breadcrumb(self.get_breadcrumb())

    def exit_menu(self):
        """Exit this menu (remove breadcrumb)."""
        self.processor.command_history.pop_breadcrumb()

    def add_navigation_options(self, options_map: Dict[str, str], option_counter: int = None) -> int:
        """
        Add standard navigation options to the menu.

        Args:
            options_map: Dictionary to add navigation options to
            option_counter: Current option number (for numbered menus)

        Returns:
            Updated option counter (if applicable)
        """
        # Build compact navigation footer
        options_map["b"] = "back"
        options_map["m"] = "main"
        options_map["x"] = "exit"

        menu_mode = self.processor.workspace.get("menu_mode", "full-menu")
        menu_label = "menu" if menu_mode == "workflow" else "main"

        self.console.print(f"\n[grey50]\\[b]ack · \\[m]{menu_label} · \\[x]exit[/grey50]")

        return option_counter

    def handle_navigation_choice(self, choice: str, options_map: Dict[str, str]) -> str:
        """
        Handle navigation choices. Returns the action to take:
        'back', 'section', 'main', 'exit', or None if not a navigation choice.
        """
        if choice in options_map:
            nav_action = options_map[choice]
            if nav_action in ["back", "section", "main", "exit"]:
                return nav_action
        return None

    def execute_navigation(self, nav_action: str) -> bool:
        """
        Execute navigation action. Returns True if navigation was handled.
        """
        if nav_action == "back":
            return False  # Exit current menu to go back
        elif nav_action == "section":
            # Navigate to category menu
            self._navigate_to_section()
            return True
        elif nav_action == "main":
            # Navigate to main menu
            self._navigate_to_main()
            return True
        elif nav_action == "exit":
            # Exit the application
            self._navigate_to_exit()
            return True
        return False

    def _navigate_to_section(self):
        """Navigate to the category/section menu."""
        if not self.menu_category:
            self.console.print("[yellow]No section menu available[/yellow]")
            return
            
        # Clear current breadcrumbs and execute appropriate category menu
        self.processor.command_history.clear_breadcrumbs()
        
        if self.menu_category == "Structure Preparation":
            from .menu_commands import StructurePreparationMenuCommand
            menu_cmd = StructurePreparationMenuCommand(self.processor)
            menu_cmd.execute()
        elif self.menu_category == "Energetic Evaluation":
            from .menu_commands import EnergeticEvaluationMenuCommand
            menu_cmd = EnergeticEvaluationMenuCommand(self.processor)
            menu_cmd.execute()
        elif self.menu_category == "Kinetic Analysis":
            from .menu_commands import KineticAnalysisMenuCommand
            menu_cmd = KineticAnalysisMenuCommand(self.processor)
            menu_cmd.execute()

    def _navigate_to_main(self):
        """Navigate to the main/workflow menu based on current mode."""
        # Clear breadcrumbs and execute appropriate menu based on mode
        self.processor.command_history.clear_breadcrumbs()

        # Check which menu mode we're in
        menu_mode = self.processor.workspace.get("menu_mode", "full-menu")

        if menu_mode == "workflow":
            from .menu_commands import WorkflowMenuCommand
            menu = WorkflowMenuCommand(self.processor)
        else:
            from .menu_commands import MainMenuCommand
            menu = MainMenuCommand(self.processor)

        menu.execute()
        
    def _navigate_to_exit(self):
        """Exit the application with confirmation."""
        from proprep.utils.prompts import confirm_with_context
        if confirm_with_context(
            self.processor,
            "Are you sure you want to exit?",
            default=False,
            module="Navigation",
            description="Confirm exit application"
        ):
            self.console.print("[green]Goodbye![/green]")
            # Signal to exit the entire application
            import sys
            sys.exit(0)


class UndoableCommand(ProcessorCommand):
    """
    Base class for commands that can be undone.

    Stores state before execution to enable undo functionality.
    """

    def __init__(self, processor: "Processor", description: str = ""):
        super().__init__(processor, description)
        self.previous_state: Optional[Dict[str, Any]] = None
        self.undone = False

    @abstractmethod
    def undo(self):
        """Undo the command."""
        pass

    def can_undo(self) -> bool:
        """Check if command can be undone."""
        return self.executed and not self.undone and self.previous_state is not None

    def save_state(self):
        """Save current state for undo."""
        # This would save relevant processor/workspace state
        # Implementation depends on what needs to be undoable
        pass


class WorkspaceCommand(ProcessorCommand):
    """
    Base class for workspace-related commands.

    Provides common functionality for commands that modify workspace state.
    """

    def __init__(self, processor: "Processor", description: str = ""):
        super().__init__(processor, description)

    def can_execute(self) -> bool:
        """Check if workspace operations are available."""
        return self.workspace is not None

    def sync_workspace(self):
        """Sync workspace state with processor."""
        self.processor._sync_workspace_to_attributes()
        self.processor._sync_attributes_to_workspace()


class ModuleCommand(ProcessorCommand):
    """
    Base class for module-related commands.

    Provides common functionality for commands that interact with processing modules.
    """

    def __init__(self, processor: "Processor", module_name: str, description: str = ""):
        super().__init__(processor, description or f"Run {module_name}")
        self.module_name = module_name

    def get_module(self):
        """Get the module instance."""
        return self.processor.get_module_instance(self.module_name)

    def can_execute(self) -> bool:
        """Check if module can be executed."""
        module = self.get_module()
        if not module:
            return False
        return module.can_process(self.workspace)

    def get_breadcrumb(self) -> str:
        """Get breadcrumb for this module."""
        return self.module_name


class ModuleActionCommand(ModuleCommand):
    """
    Command that wraps a call to a module's action method, capturing and restoring workspace state for replay.
    """

    def __init__(
        self,
        processor: "Processor",
        module_name: str,
        action_name: str,
        args: Dict[str, Any],
    ):
        # Use description to reflect the action for breadcrumbs/history
        description = f"{module_name}.{action_name}"
        super().__init__(
            processor=processor, module_name=module_name, description=description
        )
        self._action_to_run = action_name
        self.args = args or {}
        self.captured_state: Dict[str, Any] = {}  # Workspace state before execution

    def execute_with_error_handling(self) -> Any:
        """Execute command with better error handling for missing requirements."""
        try:
            module = self.get_module()
            if not module:
                raise RuntimeError(
                    f"Module '{self.module_name}' not found or not available."
                )

            if not self.can_execute():
                # Provide specific error message about missing requirements
                requirements = module.get_workspace_requirements()
                workspace = self.workspace
                missing_reqs = [
                    req for req in requirements if workspace.get(req) is None
                ]

                if missing_reqs:
                    missing_str = ", ".join(missing_reqs)
                    # Provide user-friendly error messages for common requirements
                    if "structure" in missing_reqs:
                        error_msg = "No PDB structure loaded. Please load a PDB file first using the PDB Loader module."
                    else:
                        error_msg = f"Missing required workspace data: {missing_str}"
                    raise RuntimeError(error_msg)

                raise RuntimeError(f"Cannot execute command: {self.description}")

            if self.processor.command_history:
                self.processor.command_history.add_command(self)

            # Set flag to prevent duplicate history recording in execute()
            self._in_error_handling = True
            try:
                self.result = self.execute()
            finally:
                self._in_error_handling = False

            self.executed = True
            return self.result

        except Exception as e:
            self.error = e
            self.console.print(
                f"[red]Error executing {self.description}: {str(e)}[/red]"
            )
            raise

    def execute(self) -> Any:
        """
        1. Look up the module instance
        2. Capture required workspace keys before running
        3. Call the module's method
        4. Record this command in module and global history
        5. Return the result
        """
        module_inst = self.get_module()
        if module_inst is None:
            raise RuntimeError(
                f"Cannot find module '{self.module_name}' for execution."
            )

        # 1. Capture required workspace state before running
        workspace = self.workspace
        requirements = module_inst.get_workspace_requirements()
        for key in requirements:
            self.captured_state[key] = (
                workspace.get(key) if workspace.has(key) else None
            )

        # 2. Call the module's method
        method = getattr(module_inst, self._action_to_run, None)
        if method is None:
            raise AttributeError(
                f"Module '{self.module_name}' does not have an action called '{self._action_to_run}'."
            )
        result = method(**self.args)

        # 3. Record in module's local history if available
        if hasattr(module_inst, "_record_local_command"):
            module_inst._record_local_command(self)

        # 4. Add to global history (only if not called via execute_with_error_handling)
        if self.processor.command_history and not getattr(
            self, "_in_error_handling", False
        ):
            self.processor.command_history.add_command(self)

        self.executed = True
        self.result = result
        return result

    def restore_state(self, workspace):
        """
        Restore the captured workspace state before replaying this command.
        """
        for key, val in self.captured_state.items():
            workspace.set(key, val)

    def get_breadcrumb(self) -> str:
        """
        Return a breadcrumb string for this command (module.action).
        """
        return f"{self.module_name}.{self._action_to_run}"
