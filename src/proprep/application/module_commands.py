"""
Module-related command implementations.

This module contains commands for running processing modules,
displaying module menus, and handling module workflows.
"""

from rich.prompt import Prompt, Confirm
from rich.status import Status

from Bio.PDB import PDBParser


from .processor_command import ModuleCommand, MenuCommand, ProcessorCommand
from .processor import Processor
from proprep.utils.prompts import confirm_with_context


class RunModuleCommand(ModuleCommand):
    """Command to run a processing module workflow step."""

    def __init__(self, processor: Processor, module_name: str):
        super().__init__(processor, module_name, f"Run {module_name} workflow")

    def execute(self) -> bool:
        """Execute the module workflow step."""
        module = self.get_module()
        if not module:
            self.console.print(f"[red]Module {self.module_name} not found[/red]")
            return False

        workspace = self.workspace

        # Special case for PDB Loader with pdb_id
        if self.module_name == "PDB Loader" and workspace.get("pdb_id"):
            return self._handle_pdb_loader_with_id()

        # Check if module can process current workspace
        if not module.can_process(workspace):
            self._show_requirements_error(module)
            return False

        # Process workspace using the module
        with Status(f"Running {self.module_name}..."):
            updated_workspace = module.process(workspace)

        self.console.print(f"[green]Successfully ran {self.module_name}[/green]")

        # Sync workspace to processor attributes
        self.processor._sync_workspace_to_attributes()

        return True

    def _handle_pdb_loader_with_id(self) -> bool:
        """Handle PDB Loader with pdb_id case."""
        try:
            from ..structure_prep.pdb_loader import download_pdb
            from ..structure_prep.pdb_loader import PDBMetadataExtractor
        except ImportError:
            # Fallback imports
            try:
                from proprep.structure_prep.pdb_loader import (
                    download_pdb,
                    PDBMetadataExtractor,
                )
            except ImportError:
                self.console.print(
                    "[red]Error: Could not import PDB loader modules[/red]"
                )
                return False

        pdb_id = self.workspace.get("pdb_id")

        # Download the PDB file
        pdb_file = download_pdb(pdb_id)
        if pdb_file:
            self.console.print(f"[green]Downloaded PDB file: {pdb_file}[/green]")
            self.processor.pdb_file = pdb_file
            self.workspace.set("pdb_file", pdb_file)

            # Load metadata
            self.processor.metadata = PDBMetadataExtractor(pdb_file)
            self.workspace.set("metadata", self.processor.metadata)

            # Try to load structure
            try:
                parser = PDBParser(QUIET=True)
                self.processor.structure = parser.get_structure("structure", pdb_file)
                self.workspace.set("structure", self.processor.structure)

            except Exception as e:
                self.console.print(
                    f"[yellow]Warning: Could not load structure: {str(e)}[/yellow]"
                )

            self.console.print(f"[green]Successfully loaded PDB {pdb_id}[/green]")
            return True

        return False

    def _show_requirements_error(self, module):
        """Show requirements error message."""
        requirements = module.get_workspace_requirements()
        self.console.print(
            f"[yellow]Cannot run {self.module_name} - missing requirements:[/yellow]"
        )
        for req in requirements:
            # Check if requirement is an OR condition (contains '|')
            if '|' in req:
                # Split on '|' and check if ANY of them are present
                or_parts = [part.strip() for part in req.split('|')]
                any_present = any(self.workspace.get(part) is not None for part in or_parts)
                status = "✓" if any_present else "✗"
                # Show as "key1 OR key2 OR key3"
                display_req = " OR ".join(or_parts)
                self.console.print(f"  {status} {display_req}")
            else:
                # Simple requirement
                status = "✓" if self.workspace.get(req) is not None else "✗"
                self.console.print(f"  {status} {req}")


class RunModuleMenuCommand(MenuCommand):
    """Command to display and handle a module's menu using the command pattern and breadcrumbs."""

    def __init__(self, processor: Processor, module_name: str, menu_category: str = None):
        super().__init__(processor, module_name, f"Show {module_name} menu", menu_category)
        self.module_name = module_name
        self._option_commands = None

    def _build_option_commands(self, module, menu_options):
        """Builds a mapping of option keys to Command objects for the menu."""
        option_commands = {}
        for key, (option, description) in enumerate(menu_options.items()):
            # Each option becomes a Command with a breadcrumb
            class OptionCommand(MenuCommand):
                def __init__(self, processor, module, option_key, description):
                    super().__init__(processor, description)
                    self.module = module
                    self.option_key = option_key

                def execute(self):
                    self.enter_menu()
                    try:
                        self.module.handle_menu_option(self.option_key)
                        if hasattr(self.module, "processor") and self.module.processor:
                            self.module.processor._sync_attributes_to_workspace()
                    finally:
                        self.exit_menu()

                def get_breadcrumb(self):
                    return description

            option_commands[str(key + 1)] = OptionCommand(
                self.processor, module, option, description
            )
        return option_commands

    def execute(self) -> None:
        """Display and handle the module menu using the command pattern and breadcrumbs."""
        module = self.processor.get_module_instance(self.module_name)
        if not module:
            self.console.print(f"[red]Module {self.module_name} not found[/red]")
            return

        self.enter_menu()

        if getattr(module, 'EXPERIMENTAL', False):
            self.console.print(
                "\n[yellow]Note: This module is experimental and may change in future releases.[/yellow]"
            )

        menu_options = module.get_menu_options()
        workspace = self.workspace
        debug_enabled = workspace.get("debug", False)
        if debug_enabled:
            self.console.print(
                f"[yellow]DEBUG: Module menu options for {self.module_name}: {list(menu_options.keys())}[/yellow]"
            )

        if not menu_options:
            self.console.print(
                f"[yellow]Module {self.module_name} has no menu options[/yellow]"
            )
            self.exit_menu()
            return

        # Build command objects for each menu option
        self._option_commands = self._build_option_commands(module, menu_options)

        while True:
            self.show_breadcrumbs()
            self.console.print(f"\n[bold blue]══ {self.module_name} Menu ══[/bold blue]")
            can_process = module.can_process(workspace)
            if not can_process:
                # Human-readable reason (not raw workspace keys), matching
                # the top-level menu's unmet-prerequisite note.
                try:
                    note = module.availability_note(workspace)
                except Exception:
                    note = None
                self.console.print(
                    f"[bold dark_orange3]⚠ {note or 'Requirements not met'}"
                    "[/bold dark_orange3]\n"
                )

            # Check if module supports enhanced menu display
            if hasattr(module, 'get_enhanced_menu_options'):
                # Use enhanced display with color coding and status
                from proprep.utils.enhanced_menu import EnhancedMenuDisplay
                menu_display = EnhancedMenuDisplay(self.console)

                enhanced_options = module.get_enhanced_menu_options(workspace)
                for option in enhanced_options:
                    menu_display.print_option(option)

                # Get suggestion if module provides one
                if hasattr(module, 'get_menu_suggestion'):
                    suggestion = module.get_menu_suggestion(workspace)
                    if suggestion:
                        self.console.print(f"\n[bright_blue]→ Suggestion: {suggestion}[/bright_blue]")
            else:
                # Fall back to original display
                for i, (option, description) in enumerate(menu_options.items()):
                    self.console.print(f"{i+1}. {description}")
            
            # Add navigation options using the base class method
            options_map = {}
            for i in range(len(menu_options)):
                options_map[str(i + 1)] = list(menu_options.keys())[i]
            
            current_option_num = len(menu_options) + 1
            final_option_num = self.add_navigation_options(options_map, current_option_num)

            # Create choices list for all valid options
            choices = [str(i + 1) for i in range(len(menu_options))]
            # Add letter-based navigation options
            choices.extend(["b", "s", "m", "x"])

            # Build options_map for rich context (includes menu options + navigation)
            full_options_map = {}
            for i, (option, description) in enumerate(menu_options.items()):
                full_options_map[str(i + 1)] = description
            full_options_map["b"] = "Back"
            full_options_map["s"] = "Section Menu"
            full_options_map["m"] = "Main Menu"
            full_options_map["x"] = "Exit"

            # Import prompt_with_context from menu_commands
            from .menu_commands import prompt_with_context

            choice = prompt_with_context(
                processor=self.processor,
                prompt="Enter your choice",
                choices=choices,
                default="1",
                module=f"{self.module_name} Menu",
                description=f"Select {self.module_name} option",
                options_map=full_options_map
            )

            # Handle navigation
            nav_action = self.handle_navigation_choice(choice, options_map)
            if nav_action:
                if nav_action == "back":
                    self.exit_menu()
                    break
                else:
                    # Execute navigation (section/main)
                    self.execute_navigation(nav_action)
                    self.exit_menu()  # Exit after navigation
                    break

            # Execute the selected command (option)
            option_cmd = self._option_commands.get(choice)
            if option_cmd:
                # Check this specific option's enhanced status (enhanced menu).
                if hasattr(module, 'get_enhanced_menu_options'):
                    from proprep.utils.enhanced_menu import OptionStatus
                    enhanced_options = module.get_enhanced_menu_options(workspace)
                    selected_opt = next(
                        (o for o in enhanced_options if o.key == choice), None
                    )
                    if selected_opt and selected_opt.status == OptionStatus.BLOCKED:
                        # Prerequisites genuinely unmet — explain and don't run,
                        # using the same ⚠ note style/cleaning as the menu's own
                        # per-option prereq notes (strips the legacy "[...] ○"
                        # wrapper) so the message reads consistently.
                        from proprep.utils.enhanced_menu import _OPTION_NOTE_STYLE
                        note = selected_opt.clean_note() or "Prerequisites not met"
                        self.console.print(
                            f"[{_OPTION_NOTE_STYLE}]⚠ {note}[/{_OPTION_NOTE_STYLE}]"
                        )
                        continue
                    if selected_opt and selected_opt.status == OptionStatus.AVAILABLE:
                        # An option explicitly marked AVAILABLE is structure-
                        # independent (e.g. Help, or "Import existing
                        # parameters") — run it directly, without the
                        # module-level can_process "try anyway?" gate that would
                        # otherwise block it inside a module whose *main*
                        # workflow needs a structure.
                        option_cmd.execute_with_error_handling()
                        continue
                    # Any other non-blocked option keeps the existing behaviour.
                    if not can_process:
                        if not confirm_with_context(
                            self.processor,
                            "[yellow]Module requirements not met. Try anyway?[/yellow]",
                            default=False,
                            module=self.module_name,
                            description="Confirm run despite requirements not met"
                        ):
                            continue
                    option_cmd.execute_with_error_handling()
                    continue
                if not can_process:
                    if not confirm_with_context(
                        self.processor,
                        "[yellow]Module requirements not met. Try anyway?[/yellow]",
                        default=False,
                        module=self.module_name,
                        description="Confirm run despite requirements not met"
                    ):
                        continue
                option_cmd.execute_with_error_handling()

        self.processor._sync_workspace_to_attributes()


class ListModulesCommand(ProcessorCommand):
    """Command to list all available modules."""

    def __init__(self, processor: Processor):
        super().__init__(processor, "List available modules")

    def execute(self) -> None:
        """List all available modules with their status."""
        self.console.print("\n[bold blue]══ Available Modules ══[/bold blue]")

        workspace = self.workspace

        for name, module_class in self.processor.registry.modules.items():
            module_instance = self.processor.get_module_instance(name)
            can_process = (
                module_instance.can_process(workspace) if module_instance else False
            )
            status = "[green]✓[/green]" if can_process else "[yellow]○[/yellow]"

            self.console.print(f"{status} {name} - {module_class.DESCRIPTION}")


class ShowModuleInfoCommand(ModuleCommand):
    """Command to show detailed information about a module."""

    def __init__(self, processor: Processor, module_name: str):
        super().__init__(processor, module_name, f"Show {module_name} information")

    def can_execute(self) -> bool:
        """This command can always execute as it just shows information."""
        return self.get_module() is not None

    def execute(self) -> None:
        """Show detailed module information."""
        module = self.get_module()
        if not module:
            self.console.print(f"[red]Module {self.module_name} not found[/red]")
            return

        self.console.print(f"\n[bold blue]══ {self.module_name} Information ══[/bold blue]")
        self.console.print(f"[bright_blue]Description:[/bright_blue] {module.DESCRIPTION}")

        # Show requirements
        requirements = module.get_workspace_requirements()
        if requirements:
            self.console.print(f"[bright_blue]Requirements:[/bright_blue]")
            workspace = self.workspace
            for req in requirements:
                status = "✓" if workspace.get(req) is not None else "✗"
                self.console.print(f"  {status} {req}")

        # Show menu options
        menu_options = module.get_menu_options()
        if menu_options:
            self.console.print(f"[bright_blue]Menu Options:[/bright_blue]")
            for option, description in menu_options.items():
                self.console.print(f"  • {description}")

        # Show processing status
        can_process = module.can_process(self.workspace)
        status_text = (
            "[green]Ready[/green]" if can_process else "[yellow]Not Ready[/yellow]"
        )
        self.console.print(f"[bright_blue]Status:[/bright_blue] {status_text}")
