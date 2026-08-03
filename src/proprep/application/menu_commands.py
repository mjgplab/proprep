"""
Menu command implementations.

This module contains commands for displaying and handling
the main menu and various submenus.
"""

from typing import Dict, TYPE_CHECKING

from proprep.utils.prompts import prompt_with_context, confirm_with_context

from .processor_command import MenuCommand
from proprep.utils.prompts import prompt_with_context, confirm_with_context

from .module_commands import RunModuleMenuCommand
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
    ResetCommandHistoryCommand,
)

from .processor import Processor


# Menu item styling. Kept full-brightness (no ``dim``) so the menu
# reproduces cleanly in screenshots/manuscript figures, and chosen to
# stay legible on BOTH light and dark backgrounds:
#   - ``bold blue`` is a mid-tone accent that reads on white (dark
#     blue) and on dark terminals (bright blue) alike.
#   - ``default`` follows the terminal's foreground, so the
#     description is black on a white background and white on a dark
#     one — never the invisible white-on-white that a fixed ``white``
#     would produce.
# Swap these constants to retheme the menu. ``_MENU_NOTE_STYLE`` is the
# unmet-prerequisite note shown under unavailable (○) items —
# ``dark_orange3`` (#d75f00) not ``yellow``/``dark_orange``: the lighter
# oranges only reach ~2.4-2.9:1 on white, while dark_orange3 clears the
# bold-text contrast floor (3.8:1 white, 5.5:1 black, so it survives the
# light AND dark web-shell themes) while still reading as a "warning" hue
# rather than an error (red stays reserved for errors/invalid input).
_MENU_NAME_STYLE = "bold blue"
_MENU_DESC_STYLE = "default"
_MENU_NOTE_STYLE = "bold dark_orange3"


def _print_menu_entry(console, prefix, status, name, description, tag="", note=None):
    """Render one menu item as a two-line entry (three when ``note``).

    Line 1 is the interactive row: the option key/number, an
    availability status glyph, and the module NAME in
    ``_MENU_NAME_STYLE`` — the bright scan column the eye follows down
    the menu. An optional ``tag`` (e.g. ``[exp]`` or a plugin version
    badge) trails the name.

    Line 2 is the explainer: the description in ``_MENU_DESC_STYLE``,
    hanging-indented to align under the name, so long descriptions wrap
    into a tidy block instead of stretching the row to the terminal
    edge. ``prefix`` must be plain text (no Rich markup) so its visible
    width equals its ``len``; the ``+2`` accounts for the status glyph
    and the space after it.

    Line 3 (optional) is ``note`` — a "why this is unavailable" reason
    shown only for ○ items, in ``_MENU_NOTE_STYLE`` with a ⚠ glyph,
    hanging-indented to match the description.
    """
    from rich.padding import Padding

    console.print(
        f"{prefix}{status} [{_MENU_NAME_STYLE}]{name}[/{_MENU_NAME_STYLE}]{tag}",
        highlight=False,
    )
    if description:
        console.print(
            Padding(
                f"[{_MENU_DESC_STYLE}]{description}[/{_MENU_DESC_STYLE}]",
                (0, 0, 0, len(prefix) + 2),
            ),
            highlight=False,
        )
    if note:
        console.print(
            Padding(
                f"[{_MENU_NOTE_STYLE}]⚠ {note}[/{_MENU_NOTE_STYLE}]",
                (0, 0, 0, len(prefix) + 2),
            ),
            highlight=False,
        )


def _plugin_display_label(tool_key: str) -> str:
    """Display label for a plugin sentinel (``__plugin::<name>``): the
    plugin's declared display_name, falling back to the bare name."""
    pname = tool_key.removeprefix("__plugin::")
    pmeta = WorkflowMenuCommand._plugin_metadata.get(pname)
    if pmeta and pmeta.display_name:
        return pmeta.display_name
    return pname


def _menu_entry_renderable(entry):
    """One menu item as a stacked grid cell: status glyph + key + name
    (+ trailing tag), the description, and the ⚠ note when blocked — the
    same content the list renderer shows via ``_print_menu_entry``. The
    ``tag`` is reused verbatim as Rich markup ([exp] badge / plugin chip)."""
    from rich.console import Group
    from rich.padding import Padding
    from rich.text import Text

    glyph = (
        Text("✓", style="bold green") if entry["available"]
        else Text("○", style="bold yellow")
    )
    head = Text()
    head.append(f"{entry['key']:>3}. ")
    head.append_text(glyph)
    head.append(" ")
    head.append(entry["label"], style=_MENU_NAME_STYLE)
    if entry["tag"]:
        head.append_text(Text.from_markup(entry["tag"]))
    parts = [head]
    if entry["desc"]:
        parts.append(Padding(Text(entry["desc"], style=_MENU_DESC_STYLE), (0, 0, 0, 6)))
    if entry["note"]:
        parts.append(Padding(Text(f"⚠ {entry['note']}", style=_MENU_NOTE_STYLE), (0, 0, 0, 6)))
    return Group(*parts)


def _menu_section_panel(section):
    """Bordered panel for one menu section, titled with the section
    heading and containing its entries (empty when none qualify)."""
    from rich.console import Group
    from rich.panel import Panel
    from rich.text import Text

    if section["entries"]:
        body = Group(*[_menu_entry_renderable(e) for e in section["entries"]])
    else:
        body = Text("")
    return Panel(
        body, title=Text(section["title"], style="bold blue"),
        title_align="left", border_style="blue", padding=(0, 1),
    )


def _menu_footer_panel(footer):
    """Bordered WORKSPACE & SETTINGS panel: footer actions laid out two
    per row so labels never wrap mid-word in a half-width column."""
    from rich.panel import Panel
    from rich.table import Table
    from rich.text import Text

    cells = []
    for key, label in footer:
        t = Text()
        t.append(f"{key}. ", style="bold blue")
        t.append(label, style="default")
        cells.append(t)
    grid = Table.grid(padding=(0, 3))
    grid.add_column()
    grid.add_column()
    for i in range(0, len(cells), 2):
        pair = cells[i:i + 2]
        if len(pair) == 1:
            pair.append(Text(""))
        grid.add_row(*pair)
    return Panel(
        grid, title=Text("WORKSPACE & SETTINGS", style="bold blue"),
        title_align="left", border_style="blue", padding=(0, 1),
    )


def _dispatch_plugin_in_menu(processor, console, plugin_name: str) -> None:
    """Build a HostContext, call the plugin's ``launch``, isolate
    exceptions. Shared by both ``MainMenuCommand`` (full-menu mode)
    and ``WorkflowMenuCommand`` (guided mode) so plugin handling
    behaves identically across menu modes.

    The two menu commands maintain plugin lookup tables on
    ``WorkflowMenuCommand`` (``_plugin_instances`` / ``_plugin_metadata``)
    populated by ``register_plugin_tool`` at processor startup;
    those are read here regardless of which menu is dispatching.
    """
    import logging
    from pathlib import Path
    from rich.panel import Panel

    plugin = WorkflowMenuCommand._plugin_instances.get(plugin_name)
    meta = WorkflowMenuCommand._plugin_metadata.get(plugin_name)
    if plugin is None or meta is None:
        console.print(
            f"[red]Plugin {plugin_name!r} is registered in the menu "
            "but missing from the dispatch table.[/red]"
        )
        return

    # Lazy import — keeps proprep.plugins out of menu_commands'
    # import-time graph for users who never instantiate the menu
    # commands (e.g. CLI-only batch flows).
    from proprep.plugins import HostContext, PLUGIN_API_VERSION

    seed = processor.build_plugin_seed(meta)

    host_recorder = None
    sm = getattr(processor, "session_manager", None)
    if sm is not None and getattr(sm, "is_recording", lambda: False)():
        host_recorder = sm.recorder

    host = HostContext(
        working_dir=Path.cwd(),
        seed_state=seed,
        session_recorder=host_recorder,
        console=console,
        logger=logging.getLogger(f"proprep.plugin.{plugin_name}"),
        api_version=PLUGIN_API_VERSION,
    )

    try:
        plugin.launch(host)
    except Exception as e:  # noqa: BLE001 — plugin isolation
        console.print(Panel(
            f"[red]{type(e).__name__}: {e}[/red]\n"
            f"[grey50]Plugin: {plugin_name} (v{meta.version})[/grey50]",
            title="Plugin error",
            border_style="red", expand=False,
        ))


def _format_stage_shortcut_label(stage_name: str, shortcut: str) -> str:
    """Bracket the shortcut letter inside the stage name for the
    workflow menu's nav footer.

    Mirrors the native shortcut style ``\\[l]oad`` — the letter is
    bracketed at its first case-insensitive occurrence in the stage
    name. When the shortcut letter doesn't appear in the name we
    fall back to ``<name> \\[<letter>]`` so the footer still tells
    the user which key to press.

    Examples:
      ``("Analyze", "z")``  → ``"analy\\[z]e"``  (first 'z' in
                                   "Analyze" lowercase is index 5)
      ``("Analyze", "a")``  → ``"\\[a]nalyze"``
      ``("Cooper",  "x")``  → ``"cooper \\[x]"``  (no 'x' in name)
    """
    lower_name = stage_name.lower()
    letter = shortcut.lower()
    idx = lower_name.find(letter)
    if idx < 0:
        return f"{stage_name.lower()} \\[{letter}]"
    head = stage_name[:idx].lower()
    tail = stage_name[idx + 1:].lower()
    return f"{head}\\[{letter}]{tail}"


class StructurePreparationMenuCommand(MenuCommand):
    """Command to display and handle the structure preparation menu."""

    def __init__(self, processor: Processor):
        super().__init__(processor, "Structure Preparation", "Protein structure preparation modules")

        # Fixed order of modules for display
        self.module_order = [
            "Structure Loader",  # New unified structure loader
            # "PDB Loader",  # Kept for now, will remove after testing
            # "AlphaFold Predictor",  # Unregistered - replaced by Structure Loader
            "Homology Searcher",
            "PDB Filter",
            "Structure Aligner",
            "Amino Acid Mutator",
            "Structure Fixer",
            "Redox Site Preparer",
            "MD Restraint Manager",
            "Protonation State Analyzer",
            "Membrane Builder",
            "Topology Generator",
            "Molecular Dynamics Manager",
            "QM/MM Preparator",
            # Utility modules moved to end
            "Force Field Parameterizer",
            "Structure Viewer",
        ]

    def execute(self) -> None:
        """Display and handle the structure preparation menu."""
        self.enter_menu()

        try:
            while True:
                self.show_breadcrumbs()
                self._show_menu()

                choice = self._get_user_choice()
                if not self._handle_choice(choice):
                    break

        finally:
            self.exit_menu()
            
    def _show_menu(self):
        """Display the structure preparation menu options."""
        self.console.print(
            "\n[bold blue]══ Structure Preparation Modules ══[/bold blue]"
        )
        self.console.print(
            "[italic]Prepare protein structures for MD simulations[/italic]\n"
        )

        # Find the longest module name for proper alignment
        longest_name = max(
            len(name)
            for name in self.module_order
            if name in self.processor.registry.modules
        )

        # Calculate the max number width for right alignment
        max_option_num = len(self.module_order) + 1  # +1 for back option
        num_width = len(str(max_option_num))

        # Display modules in fixed order with aligned descriptions
        option_num = 1
        self.options_map = {}

        workspace = self.workspace

        for name in self.module_order:
            if name in self.processor.registry.modules:
                module_class = self.processor.registry.modules[name]
                module_desc = module_class.DESCRIPTION

                # Check if requirements are met
                module_instance = self.processor.get_module_instance(name)
                can_process = module_instance.can_process(workspace)
                status = "[green]✓[/green]" if can_process else "[yellow]○[/yellow]"

                # Use a fixed width for the module name to align all descriptions
                padded_name = name.ljust(longest_name)

                # Right-align the option number with consistent width
                formatted_num = f"{option_num}.".rjust(num_width + 1)

                self.console.print(
                    f"{formatted_num} {status} {padded_name} - {module_desc}",
                    highlight=False
                )
                self.options_map[str(option_num)] = name
                option_num += 1

        # Add navigation options using the base class method
        option_num = self.add_navigation_options(self.options_map, option_num)

    def _get_user_choice(self) -> str:
        """Get user menu choice."""
        return prompt_with_context(
            processor=self.processor,
            prompt="\nEnter your choice",
            choices=list(self.options_map.keys()),
            default="1",
            module="Structure Preparation Menu",
            description="Select structure preparation module",
            options_map=self.options_map
        )

    def _handle_choice(self, choice: str) -> bool:
        """Handle user menu choice. Returns True to continue, False to exit."""
        # Check if it's a navigation choice first
        nav_action = self.handle_navigation_choice(choice, self.options_map)
        if nav_action:
            if nav_action == "back":
                return False  # Exit to main menu
            else:
                # Execute navigation (section/main) - this will handle the navigation
                self.execute_navigation(nav_action)
                return False  # Exit this menu since navigation was executed

        # Handle regular module selection
        selected = self.options_map[choice]
        module_menu = RunModuleMenuCommand(self.processor, selected, "Structure Preparation")
        module_menu.execute()

        return True


class EnergeticEvaluationMenuCommand(MenuCommand):
    """Command to display and handle the energetic evaluation menu."""

    def __init__(self, processor: Processor):
        super().__init__(processor, "Energetic Evaluation", "Energy analysis and evaluation modules")

        # Module order for energetic evaluation
        self.module_order = [
            "Pathway Finder",
            "Interaction Energy Calculator"
        ]

    def execute(self) -> None:
        """Display and handle the energetic evaluation menu."""
        self.enter_menu()

        try:
            while True:
                self.show_breadcrumbs()
                self._show_menu()

                choice = self._get_user_choice()
                if not self._handle_choice(choice):
                    break

        finally:
            self.exit_menu()
            
    def _show_menu(self):
        """Display the energetic evaluation menu options."""
        self.console.print(
            "\n[bold blue]══ Energetic Evaluation Modules ══[/bold blue]"
        )
        self.console.print(
            "[italic]Analyze and evaluate molecular energetics[/italic]\n"
        )

        if not self.module_order:
            self.console.print("[yellow]No energetic evaluation modules available yet.[/yellow]")
            self.console.print("[italic]Modules will be added here as they are integrated.[/italic]\n")
            self.options_map = {"1": "back"}
            self.console.print("1. Back to Main Menu", highlight=False)
        else:
            # Find the longest module name for proper alignment
            longest_name = max(
                len(name)
                for name in self.module_order
                if name in self.processor.registry.modules
            )

            # Calculate the max number width for right alignment
            max_option_num = len(self.module_order) + 1  # +1 for back option
            num_width = len(str(max_option_num))

            # Display modules in fixed order with aligned descriptions
            option_num = 1
            self.options_map = {}

            workspace = self.workspace

            for name in self.module_order:
                if name in self.processor.registry.modules:
                    module_class = self.processor.registry.modules[name]
                    module_desc = module_class.DESCRIPTION

                    # Check if requirements are met
                    module_instance = self.processor.get_module_instance(name)
                    can_process = module_instance.can_process(workspace)
                    status = "[green]✓[/green]" if can_process else "[yellow]○[/yellow]"

                    # Use a fixed width for the module name to align all descriptions
                    padded_name = name.ljust(longest_name)

                    # Right-align the option number with consistent width
                    formatted_num = f"{option_num}.".rjust(num_width + 1)

                    self.console.print(
                        f"{formatted_num} {status} {padded_name} - {module_desc}",
                        highlight=False
                    )
                    self.options_map[str(option_num)] = name
                    option_num += 1

            # Add navigation options using the base class method
            option_num = self.add_navigation_options(self.options_map, option_num)

    def _get_user_choice(self) -> str:
        """Get user menu choice."""
        return prompt_with_context(
            processor=self.processor,
            prompt="\nEnter your choice",
            choices=list(self.options_map.keys()),
            default="1",
            module="Energetic Evaluation Menu",
            description="Select energetic evaluation module",
            options_map=self.options_map
        )

    def _handle_choice(self, choice: str) -> bool:
        """Handle user menu choice. Returns True to continue, False to exit."""
        # Check if it's a navigation choice first
        nav_action = self.handle_navigation_choice(choice, self.options_map)
        if nav_action:
            if nav_action == "back":
                return False  # Exit to main menu
            else:
                # Execute navigation (section/main)
                self.execute_navigation(nav_action)
                return False  # Exit this menu since navigation was executed

        # Handle regular module selection
        selected = self.options_map[choice]
        module_menu = RunModuleMenuCommand(self.processor, selected, "Energetic Evaluation")
        module_menu.execute()

        return True


class KineticAnalysisMenuCommand(MenuCommand):
    """Command to display and handle the kinetic analysis menu."""

    def __init__(self, processor: Processor):
        super().__init__(processor, "Kinetic Analysis", "Kinetic analysis and simulation modules")

        # Module order for kinetic analysis
        self.module_order = [
            "Diffusion Analysis",
            "Flux Analysis", 
            "Parameter Optimization"
        ]

    def execute(self) -> None:
        """Display and handle the kinetic analysis menu."""
        self.enter_menu()

        try:
            while True:
                self.show_breadcrumbs()
                self._show_menu()

                choice = self._get_user_choice()
                if not self._handle_choice(choice):
                    break

        finally:
            self.exit_menu()
            
    def _show_menu(self):
        """Display the kinetic analysis menu options."""
        self.console.print(
            "\n[bold blue]══ Kinetic Analysis Modules ══[/bold blue]"
        )
        self.console.print(
            "[italic]Analyze molecular kinetics and dynamics[/italic]\n"
        )

        if not self.module_order:
            self.console.print("[yellow]No kinetic analysis modules available yet.[/yellow]")
            self.console.print("[italic]Modules will be added here as they are integrated.[/italic]\n")
            self.options_map = {"1": "back"}
            self.console.print("1. Back to Main Menu", highlight=False)
        else:
            # Find the longest module name for proper alignment
            longest_name = max(
                len(name)
                for name in self.module_order
                if name in self.processor.registry.modules
            )

            # Calculate the max number width for right alignment
            max_option_num = len(self.module_order) + 1  # +1 for back option
            num_width = len(str(max_option_num))

            # Display modules in fixed order with aligned descriptions
            option_num = 1
            self.options_map = {}

            workspace = self.workspace

            for name in self.module_order:
                if name in self.processor.registry.modules:
                    module_class = self.processor.registry.modules[name]
                    module_desc = module_class.DESCRIPTION

                    # Check if requirements are met
                    module_instance = self.processor.get_module_instance(name)
                    can_process = module_instance.can_process(workspace)
                    status = "[green]✓[/green]" if can_process else "[yellow]○[/yellow]"

                    # Use a fixed width for the module name to align all descriptions
                    padded_name = name.ljust(longest_name)

                    # Right-align the option number with consistent width
                    formatted_num = f"{option_num}.".rjust(num_width + 1)

                    self.console.print(
                        f"{formatted_num} {status} {padded_name} - {module_desc}",
                        highlight=False
                    )
                    self.options_map[str(option_num)] = name
                    option_num += 1

            # Add navigation options using the base class method
            option_num = self.add_navigation_options(self.options_map, option_num)

    def _get_user_choice(self) -> str:
        """Get user menu choice."""
        return prompt_with_context(
            processor=self.processor,
            prompt="\nEnter your choice",
            choices=list(self.options_map.keys()),
            default="1",
            module="Kinetic Analysis Menu",
            description="Select kinetic analysis module",
            options_map=self.options_map
        )

    def _handle_choice(self, choice: str) -> bool:
        """Handle user menu choice. Returns True to continue, False to exit."""
        # Check if it's a navigation choice first
        nav_action = self.handle_navigation_choice(choice, self.options_map)
        if nav_action:
            if nav_action == "back":
                return False  # Exit to main menu
            else:
                # Execute navigation (section/main)
                self.execute_navigation(nav_action)
                return False  # Exit this menu since navigation was executed

        # Handle regular module selection
        selected = self.options_map[choice]
        module_menu = RunModuleMenuCommand(self.processor, selected, "Kinetic Analysis")
        module_menu.execute()

        return True


class MainMenuCommand(MenuCommand):
    """Command to display and handle the main menu."""

    def __init__(self, processor: Processor):
        super().__init__(processor, "Main Menu", "Show main menu")

    def execute(self) -> None:
        """Display and handle the main menu."""
        self.enter_menu()

        try:
            while True:
                # If menu mode was changed elsewhere (e.g. by a nested
                # menu instance created via _navigate_to_main), exit so
                # run_main_menu can reload with the correct menu class.
                if self.workspace.get("menu_mode", "full-menu") != "full-menu":
                    break

                self.show_breadcrumbs()
                self._show_menu()

                choice = self._get_user_choice()
                if not self._handle_choice(choice):
                    break

        finally:
            self.exit_menu()

    def _show_menu(self):
        """Display the main menu with structure preparation modules."""
        workspace = self.workspace

        # Define module groups with section-based numbering. Plugin
        # sections are appended after the native ones (just before
        # UTILITIES) so the user sees the same numbering for native
        # tools regardless of which plugins happen to be installed.
        # Within the plugin block, sections are ordered by
        # stage_order, ties broken alphabetically by plugin name.
        module_groups = [
            ("1. STRUCTURE INPUT", [
                "Structure Loader",
                "Biological Assembly Generator",
            ]),
            ("2. REDOX SITE DETECTION", [
                "Redox Site Detector",
                "Force Field Explorer",
                "Force Field Parameterizer",
            ]),
            ("3. STRUCTURE COMPARISON", [
                "Homology Searcher",
                "Structure Aligner",
            ]),
            ("4. STRUCTURE PREPARATION", [
                "PDB Filter",
                "Amino Acid Mutator",
                "Structure Fixer",
                "Redox Site Preparer",
                "Protonation State Analyzer",
                "Structure Orientator",
                "Membrane Builder",
            ]),
            ("5. SIMULATION SETUP & EXECUTION", [
                "Topology Generator",
                "Molecular Dynamics Manager",
                "QM/MM Preparator",
            ]),
        ]
        # Plugin-contributed sections — one per plugin stage, with
        # all plugins targeting the same stage grouped together.
        # Section number begins after the last native section so the
        # native numbering stays stable.
        next_section_num = len(module_groups) + 1
        plugin_sections_by_stage: Dict[str, List[str]] = {}
        plugin_stage_order: Dict[str, int] = {}
        for pname, meta in WorkflowMenuCommand._plugin_metadata.items():
            sentinel = f"__plugin::{pname}"
            plugin_sections_by_stage.setdefault(meta.stage, []).append(sentinel)
            # First-registered plugin's stage_order wins for the
            # section's display position; matches the precedent in
            # WorkflowStateManager.register_plugin_stage.
            plugin_stage_order.setdefault(meta.stage, meta.stage_order)
        # Sort plugin stages by stage_order then alphabetically.
        # Lazy import — keeps startup paths that never instantiate
        # MainMenuCommand from paying for the workflow_state_manager
        # module (matches the pattern WorkflowMenuCommand already uses
        # for the same import).
        from proprep.utils.workflow_state_manager import WorkflowStateManager

        for stage in sorted(
            plugin_sections_by_stage,
            key=lambda s: (plugin_stage_order[s], s),
        ):
            stage_label = (
                WorkflowStateManager.STAGE_NAMES.get(stage, stage).upper()
            )
            module_groups.append(
                (f"{next_section_num}. {stage_label}",
                 plugin_sections_by_stage[stage]),
            )
            next_section_num += 1
        module_groups.append(
            ("UTILITIES", ["Structure Viewer"]),
        )

        self.options_map = {}

        # Build a structured model from ``module_groups`` and populate
        # ``self.options_map``. Rendering is a pure function of this
        # model, so the same data drives either the single-column list
        # (default) or the 2-column panel grid (``menu_layout == "grid"``).
        sections = []
        for section_title, module_names in module_groups:
            section_num = (
                section_title.split(".")[0]
                if section_title.split(".")[0].isdigit() else None
            )
            entries = []
            subsection_letter = 'a'
            for name in module_names:
                option_key = (
                    f"{section_num}{subsection_letter}" if section_num
                    else f"u{subsection_letter}"
                )
                entry = self._build_menu_entry(workspace, option_key, name)
                if entry is not None:
                    entries.append(entry)
                    subsection_letter = chr(ord(subsection_letter) + 1)
            sections.append({"title": section_title, "entries": entries})

        # Workspace & settings footer — same keys/actions in both layouts.
        footer = [
            ("w", "Workspace Options"),
            ("p", "Preferences"),
            ("h", "Help & Workflows"),
            ("f", "Feedback & Support"),
            ("g", "Switch to Workflow Mode"),
            ("x", "Exit"),
        ]
        for key, action in (
            ("w", "workspace_options"), ("p", "preferences"), ("h", "help"),
            ("f", "feedback"), ("g", "workflow_mode"), ("x", "exit"),
        ):
            self.options_map[key] = action

        if self._menu_layout() == "grid":
            self._render_menu_grid(sections, footer)
        else:
            self._render_menu_list(sections, footer)

    def _menu_layout(self) -> str:
        """Resolve the full-menu layout: ``"grid"`` or ``"list"``.

        A per-run workspace override (set from the ``--menu-grid`` /
        ``--menu-list`` CLI flags) wins; otherwise the persisted
        preference from ``SettingsManager`` applies (default ``list``).
        """
        override = self.workspace.get("menu_layout")
        if override in ("list", "grid"):
            return override
        from proprep.utils.settings_manager import SettingsManager
        return SettingsManager().get_menu_layout()

    def _build_menu_entry(self, workspace, option_key, name):
        """Resolve one menu item to a display record, or ``None`` when it
        is not shown in this build (unregistered module name).

        Centralizes the three item kinds — the Force Field Explorer
        pseudo-tool, plugin sentinels (``__plugin::<name>``), and
        registered modules — so the list and grid renderers stay in
        lock-step and ``options_map`` is populated identically. ``tag`` is
        trailing Rich markup (the ``[exp]`` badge or a plugin version
        chip) consumed verbatim by both renderers.
        """
        if name == "Force Field Explorer":
            self.options_map[option_key] = "forcefield_explorer"
            return {
                "key": option_key, "label": name,
                "desc": "Browse force field parameters",
                "available": True, "tag": "", "note": None,
            }

        if name.startswith("__plugin::"):
            plugin_name = name.removeprefix("__plugin::")
            meta = WorkflowMenuCommand._plugin_metadata.get(plugin_name)
            self.options_map[option_key] = name
            return {
                "key": option_key, "label": _plugin_display_label(name),
                "desc": meta.description if meta else "(plugin)",
                "available": True,
                "tag": f"  [grey50]\\[plugin v{meta.version if meta else '?'}][/grey50]",
                "note": None,
            }

        if name in self.processor.registry.modules:
            module_class = self.processor.registry.modules[name]
            module_instance = self.processor.get_module_instance(name)
            can_process = module_instance.can_process(workspace)
            exp_tag = (
                " [bold dark_orange3]\\[exp][/bold dark_orange3]"
                if getattr(module_class, 'EXPERIMENTAL', False) else ""
            )
            # Explain WHY an unavailable (○) module is blocked.
            note = None
            if not can_process:
                try:
                    note = module_instance.availability_note(workspace)
                except Exception:
                    note = None
            self.options_map[option_key] = name
            return {
                "key": option_key, "label": name,
                "desc": module_class.DESCRIPTION,
                "available": can_process, "tag": exp_tag, "note": note,
            }

        return None

    def _render_menu_list(self, sections, footer):
        """Single-column list renderer (default): one heading per
        section, two/three lines per entry via ``_print_menu_entry``."""
        for section in sections:
            self.console.print(
                f"\n[bold blue]══ {section['title']} ══[/bold blue]",
                highlight=False,
            )
            for e in section["entries"]:
                prefix = f" {e['key']:>3}. "
                status = (
                    "[bold green]✓[/bold green]" if e["available"]
                    else "[bold yellow]○[/bold yellow]"
                )
                _print_menu_entry(
                    self.console, prefix, status,
                    e["label"], e["desc"], e["tag"], e["note"],
                )

        self.console.print("\n[bold blue]══ WORKSPACE & SETTINGS ══[/bold blue]")
        for key, label in footer:
            self.console.print(f"  {key}. {label}")

    def _render_menu_grid(self, sections, footer):
        """2-column panel-grid renderer (``menu_layout == "grid"``).

        Each section becomes a bordered panel; sections 1,3,5… stack in
        the left column and 2,4,6… in the right. Paired rows get their
        top edges aligned by measuring the rendered grid and padding the
        higher panel down, while a trailing unpaired panel flows up tight.
        UTILITIES and the workspace footer sit across the bottom. Status
        glyphs, ``[exp]`` tags, and ⚠ notes are all preserved.
        """
        from rich.console import Console, Group
        from rich.table import Table
        from rich.text import Text

        width = self.console.width or 118

        # Render UTILITIES with the workspace footer along the bottom;
        # the numbered/plugin sections form the 2-column body above it.
        util_section = None
        body_sections = []
        for s in sections:
            if s["title"] == "UTILITIES":
                util_section = s
            else:
                body_sections.append(s)

        panels = [_menu_section_panel(s) for s in body_sections]

        def _two_col(left, right):
            g = Table.grid(expand=True, padding=(0, 1))
            g.add_column(ratio=1)
            g.add_column(ratio=1)
            g.add_row(left, right)
            return g

        left_panels = [p for i, p in enumerate(panels) if i % 2 == 0]
        right_panels = [p for i, p in enumerate(panels) if i % 2 == 1]
        left_titles = [s["title"] for i, s in enumerate(body_sections) if i % 2 == 0]
        right_titles = [s["title"] for i, s in enumerate(body_sections) if i % 2 == 1]
        left_pad = [0] * len(left_panels)
        right_pad = [0] * len(right_panels)

        def _column(col_panels, pads):
            parts = []
            for i, p in enumerate(col_panels):
                if i:
                    parts.append(Text(""))
                parts.extend(Text("") for _ in range(pads[i]))
                parts.append(p)
            return Group(*parts) if parts else Text("")

        def _title_row(grid, title):
            con = Console(width=width)
            lines = con.render_lines(grid, con.options.update_width(width))
            for idx, segs in enumerate(lines):
                if title in "".join(s.text for s in segs):
                    return idx
            return -1

        # Top-down: align each paired row's top edges by pushing the
        # higher panel down. Re-measure each step so lower rows account
        # for padding added above them. Width-independent (reads the
        # real rendered grid rather than assuming a column width).
        for r in range(min(len(left_panels), len(right_panels))):
            grid = _two_col(_column(left_panels, left_pad), _column(right_panels, right_pad))
            lt = _title_row(grid, left_titles[r])
            rt = _title_row(grid, right_titles[r])
            if lt >= 0 and rt >= 0:
                if lt < rt:
                    left_pad[r] += rt - lt
                elif rt < lt:
                    right_pad[r] += lt - rt

        top = _two_col(_column(left_panels, left_pad), _column(right_panels, right_pad))

        util_panel = _menu_section_panel(util_section) if util_section else Text("")
        bottom = _two_col(util_panel, _menu_footer_panel(footer))

        self.console.print(top)
        self.console.print()
        self.console.print(bottom)

    def _get_user_choice(self) -> str:
        """Get user menu choice."""
        # Don't use Prompt.ask with choices parameter as it shows the entire list
        # Instead, validate manually
        while True:
            response = prompt_with_context(
                self.processor,
                "\nEnter your choice",
                default="1a",
                module="Main Menu",
                description="Select menu option",
                options_map=self.options_map
            )

            # Validate the choice
            if response in self.options_map:
                return response
            else:
                self.console.print(f"[red]Invalid choice: {response}. Please try again.[/red]")
                self.console.print("[grey50]Hint: Use format like '1a', '3c', 'w', 'h', or 'x'[/grey50]")

    def _handle_choice(self, choice: str) -> bool:
        """Handle user menu choice. Returns True to continue, False to exit."""
        selected = self.options_map[choice]

        # Plugin tools are encoded as ``__plugin::<name>`` sentinels.
        # Dispatch them before the generic module-selection path so a
        # plugin name can never collide with a registered module.
        # Same dispatcher as WorkflowMenuCommand uses, so plugin
        # behaviour is mode-invariant.
        if selected.startswith("__plugin::"):
            _dispatch_plugin_in_menu(
                self.processor, self.console,
                selected.removeprefix("__plugin::"),
            )
            return True

        if selected == "exit":
            if confirm_with_context(
                self.processor,
                "Are you sure you want to exit?",
                default=False,
                module="Main Menu",
                description="Confirm exit"
            ):
                self.console.print("[green]Goodbye![/green]")
                import sys
                sys.exit(0)  # Exit the application
            return True  # Continue the menu loop
        elif selected == "workspace_options":
            workspace_menu = WorkspaceOptionsMenuCommand(self.processor)
            workspace_menu.execute()
        elif selected == "preferences":
            preferences_menu = PreferencesMenuCommand(self.processor)
            preferences_menu.execute()
        elif selected == "help":
            self._show_help_and_workflows()
        elif selected == "feedback":
            from .feedback_command import FeedbackCommand
            feedback_cmd = FeedbackCommand(self.processor)
            feedback_cmd.execute()
        elif selected == "forcefield_explorer":
            from .forcefield_explorer import ForceFieldExplorerCommand
            explorer = ForceFieldExplorerCommand(self.processor)
            explorer.execute()
        elif selected == "workflow_mode":
            # Switch to workflow mode
            self.processor.workspace.set("menu_mode", "workflow")
            self.console.print("[bright_blue]Switching to Workflow Mode...[/bright_blue]")
            return False  # Exit to trigger menu reload
        else:
            # Handle module selection directly
            module_menu = RunModuleMenuCommand(self.processor, selected, "ProPrep Main")
            module_menu.execute()

        return True

    def _show_help_and_workflows(self):
        """Display help information and common workflows."""
        self.console.print("\n[bold bright_blue]═══════════════════════════════════════════════════════════════[/bold bright_blue]")
        self.console.print("[bold bright_blue]                    Common Workflows[/bold bright_blue]")
        self.console.print("[bold bright_blue]═══════════════════════════════════════════════════════════════[/bold bright_blue]\n")

        workflows = [
            ("1️⃣  BASIC PROTEIN MD SIMULATION",
             "1a.Load → 3a.Filter → 3c.Fix → 3d.Protonate → 5b.tLEaP → 5c.MD"),

            ("2️⃣  METALLOPROTEIN MD SIMULATION",
             "1a.Load → 3a.Filter → 3c.Fix → 4a.RedoxSite → 5a.FF-Params → 5b.tLEaP → 5c.MD"),

            ("3️⃣  QM/MM CALCULATION",
             "1a.Load → 3a.Filter → 3c.Fix → 4a.RedoxSite → 5b.tLEaP → 5d.QM/MM"),

            ("4️⃣  STRUCTURE COMPARISON",
             "1a.Load (multiple) → 2b.Aligner → ua.Viewer"),

            ("5️⃣  CONTINUE EXISTING MD",
             "5c.MD (load existing topology/coordinates)"),
        ]

        for title, workflow in workflows:
            self.console.print(f"[bold]{title}[/bold]")
            self.console.print(f"   {workflow}\n")

        self.console.print("\n[bold]Tool Categories:[/bold]\n")

        categories = [
            ("STRUCTURE INPUT", "Load PDB files from database or local files"),
            ("REDOX SITE DETECTION", "Detect and parameterize non-standard residues (optional)"),
            ("STRUCTURE COMPARISON", "Sequence homology and structural alignment (optional)"),
            ("STRUCTURE PREPARATION", "Clean, filter, fix, and prepare structures"),
            ("SIMULATION SETUP & EXECUTION", "Create topologies and run MD/QM simulations"),
            ("UTILITIES", "Visualization and workspace management"),
        ]

        for category, description in categories:
            self.console.print(f"[bright_blue]{category}[/bright_blue]: {description}")

        self.console.print("\n[grey50]Press Enter to return to main menu...[/grey50]")
        prompt_with_context(
            self.processor,
            "",
            default="",
            module="Main Menu",
            description="Pause after category view",
        )


class WorkflowMenuCommand(MenuCommand):
    """Command to display workflow mode menu with progressive disclosure."""

    # Map workflow stages to tool categories
    STAGE_TOOLS = {
        "load": {
            "section": "1. STRUCTURE INPUT",
            "tools": ["Structure Loader", "Biological Assembly Generator"]
        },
        "detect": {
            "section": "2. REDOX SITE DETECTION",
            "tools": ["Redox Site Detector", "Force Field Explorer", "Force Field Parameterizer"]
        },
        "compare": {
            "section": "3. STRUCTURE COMPARISON",
            "tools": ["Homology Searcher", "Structure Aligner"]
        },
        "fix": {
            "section": "4. STRUCTURE PREPARATION",
            "tools": ["PDB Filter", "Amino Acid Mutator", "Structure Fixer", "Redox Site Preparer", "Protonation State Analyzer", "Structure Orientator", "Membrane Builder"]
        },
        "simulate": {
            "section": "5. SIMULATION SETUP & EXECUTION",
            "tools": ["Topology Generator", "Molecular Dynamics Manager", "QM/MM Preparator"]
        },
    }

    # Plugin lookup tables, populated by ``register_plugin_tool`` at
    # PDBProcessor init time. Keyed by the plugin's metadata.name so
    # the menu dispatcher can route a sentinel ``__plugin::<name>``
    # selection back to the plugin instance + its metadata.
    _plugin_instances = {}        # Dict[str, ProPrepPlugin]
    _plugin_metadata = {}         # Dict[str, PluginMetadata]

    @classmethod
    def register_plugin_tool(cls, meta, plugin) -> None:
        """Splice a plugin tool into the STAGE_TOOLS entry for its stage.

        Creates the stage entry if it doesn't yet exist (the plugin's
        stage was just added to ``WorkflowStateManager.STAGES`` by
        ``register_plugin_stage``, but no native tools live there).
        Tool name encodes the plugin's identity as
        ``__plugin::<name>`` — the dispatcher checks for that prefix
        and routes accordingly. Idempotent: re-registering the same
        plugin replaces its prior entry rather than duplicating it.
        """
        sentinel = f"__plugin::{meta.name}"
        stage_entry = cls.STAGE_TOOLS.setdefault(
            meta.stage,
            {
                "section": (
                    f"{meta.stage_order // 10 + 1}. {meta.stage_name.upper()}"
                ),
                "tools": [],
            },
        )
        # Drop any prior entry for this plugin first (idempotency).
        stage_entry["tools"] = [
            t for t in stage_entry["tools"] if t != sentinel
        ]
        # Insert respecting tool_order: collect existing plugin tools'
        # order values, sort with this one inserted.
        ordered_pairs = []
        for t in stage_entry["tools"]:
            if t.startswith("__plugin::"):
                pname = t.removeprefix("__plugin::")
                pmeta = cls._plugin_metadata.get(pname)
                ordered_pairs.append((pmeta.tool_order if pmeta else 0, t))
            else:
                # Native modules sort first — they're not plugins.
                ordered_pairs.append((-1, t))
        ordered_pairs.append((meta.tool_order, sentinel))
        ordered_pairs.sort(key=lambda kv: (kv[0], kv[1]))
        stage_entry["tools"] = [t for _, t in ordered_pairs]

        cls._plugin_instances[meta.name] = plugin
        cls._plugin_metadata[meta.name] = meta

    def __init__(self, processor: Processor):
        super().__init__(processor, "Main Menu (Workflow Mode)", "Workflow mode menu")

        # Create workflow state manager
        from proprep.utils.workflow_state_manager import WorkflowStateManager
        self.workflow_state = WorkflowStateManager(processor.workspace)

    def execute(self) -> None:
        """Display and handle workflow mode menu."""
        self.enter_menu()

        try:
            while True:
                # If menu mode was changed elsewhere (e.g. by a nested
                # menu instance created via _navigate_to_main), exit so
                # run_main_menu can reload with the correct menu class.
                if self.workspace.get("menu_mode", "full-menu") != "workflow":
                    break

                self.show_breadcrumbs()
                self._show_menu()

                choice = self._get_user_choice()
                if not self._handle_choice(choice):
                    break

        finally:
            self.exit_menu()

    def _show_menu(self):
        """Display workflow menu showing only current stage tools."""
        from rich.text import Text
        from rich.console import Console

        workspace = self.workspace

        # Show updated workflow guide panel
        self._show_workflow_guide_v2(workspace)

        # Get current stage and its tools
        current_stage = self.workflow_state.get_current_stage()
        stage_config = self.STAGE_TOOLS[current_stage]

        self.options_map = {}

        # Display current stage section
        section_title = stage_config["section"]
        module_names = stage_config["tools"]

        self.console.print(f"\n[bold blue]══ {section_title} ══[/bold blue]", highlight=False)

        # Find longest name for alignment. Plugin tools are encoded as
        # ``__plugin::<name>`` sentinels — for layout we use the
        # plugin's display label (display_name when the plugin set
        # one, else falling back to name) instead of the sentinel.
        def _display_label(tool_key: str) -> str:
            if tool_key.startswith("__plugin::"):
                pname = tool_key.removeprefix("__plugin::")
                meta = WorkflowMenuCommand._plugin_metadata.get(pname)
                if meta and meta.display_name:
                    return meta.display_name
                return pname
            return tool_key

        displayable = [
            name for name in module_names
            if name in self.processor.registry.modules
            or name == "Force Field Explorer"
            or name.startswith("__plugin::")
        ]
        if displayable:
            # Display modules in this section with simple numbering
            option_num = 1
            for name in module_names:
                prefix = f"  {option_num}. "

                if name == "Force Field Explorer":
                    ff_desc = "Browse force field parameters"
                    _print_menu_entry(
                        self.console, prefix, "[bold green]✓[/bold green]",
                        name, ff_desc,
                    )
                    self.options_map[str(option_num)] = "forcefield_explorer"
                    option_num += 1
                elif name.startswith("__plugin::"):
                    # Plugin tool — read description + version from
                    # registered metadata; status is always available
                    # because is_available was True at discovery (a
                    # later runtime regression would surface as an
                    # exception in launch, caught by the dispatcher).
                    plugin_name = name.removeprefix("__plugin::")
                    meta = self._plugin_metadata.get(plugin_name)
                    desc = meta.description if meta else "(plugin)"
                    # Prefer display_name (plugin's branded label,
                    # e.g. "ETAnalyze") over the lowercase identifier.
                    label = _display_label(name)
                    _print_menu_entry(
                        self.console, prefix, "[bold green]✓[/bold green]",
                        label, desc,
                        f"  [grey50]\\[plugin v{meta.version if meta else '?'}][/grey50]",
                    )
                    self.options_map[str(option_num)] = name
                    option_num += 1
                elif name in self.processor.registry.modules:
                    module_class = self.processor.registry.modules[name]
                    module_desc = module_class.DESCRIPTION

                    # Check if requirements are met
                    module_instance = self.processor.get_module_instance(name)
                    can_process = module_instance.can_process(workspace)
                    status = "[bold green]✓[/bold green]" if can_process else "[bold yellow]○[/bold yellow]"

                    exp_tag = " [bold dark_orange3]\\[exp][/bold dark_orange3]" if getattr(module_class, 'EXPERIMENTAL', False) else ""

                    # Explain WHY an unavailable (○) module is blocked.
                    note = None
                    if not can_process:
                        try:
                            note = module_instance.availability_note(workspace)
                        except Exception:
                            note = None

                    _print_menu_entry(
                        self.console, prefix, status,
                        name, module_desc, exp_tag, note,
                    )
                    self.options_map[str(option_num)] = name
                    option_num += 1

        # Quick access and utilities
        self.console.print(f"\n  v. Structure Viewer  |  u. Utilities...")
        self.options_map["v"] = "Structure Viewer"
        self.options_map["u"] = "utilities"

        # Compact navigation footer
        nav_parts = []
        if self.workflow_state.can_go_back():
            nav_parts.append("\\[b]ack")
            self.options_map["b"] = "back"
        if self.workflow_state.can_advance():
            nav_parts.append("\\[n]ext")
            self.options_map["n"] = "next"

        # Stage jump shortcuts. Native ones are hardcoded; plugin
        # shortcuts (declared via PluginMetadata.stage_shortcut) get
        # appended below if they don't collide.
        stage_shortcuts = {
            "l": "load", "d": "detect", "c": "compare",
            "f": "fix", "s": "simulate",
        }
        stage_label_parts = [
            "\\[l]oad", "\\[d]etect", "\\[c]ompare", "\\[f]ix", "\\[s]imulate",
        ]
        # Reserved keys that the menu uses for things other than
        # stage jumps — plugin shortcuts may not steal any of these.
        # Single source of truth so the next person adding a top-level
        # shortcut doesn't have to remember to update collision logic.
        _reserved_nonstage_keys = {"a", "b", "n", "u", "v", "x"}

        for plugin_name, meta in WorkflowMenuCommand._plugin_metadata.items():
            if not meta.stage_shortcut:
                continue
            key = meta.stage_shortcut.lower()
            if (
                key in stage_shortcuts
                or key in _reserved_nonstage_keys
                or len(key) != 1
            ):
                # Silent skip — the stage is still reachable via
                # [n]ext / [b]ack. A logger.warning at registration
                # time would be the place to surface the collision
                # to plugin authors; for now keep the menu clean.
                continue
            stage_shortcuts[key] = meta.stage
            stage_label_parts.append(
                _format_stage_shortcut_label(meta.stage_name, key)
            )

        stage_labels = " ".join(stage_label_parts)
        nav_parts.append(stage_labels)
        for key, stage in stage_shortcuts.items():
            self.options_map[key] = f"jump:{stage}"

        nav_parts.append("\\[a]ll tools")
        self.options_map["a"] = "all"
        nav_parts.append("e\\[x]it")
        self.options_map["x"] = "exit"

        # highlight=False so Rich's ReprHighlighter doesn't style the
        # "[b]ack" brackets as braces — that styling renders them
        # invisibly on a light/white terminal background.
        self.console.print(f"\n{' | '.join(nav_parts)}", highlight=False)

    def _show_workflow_guide_v2(self, workspace):
        """Display workflow progress indicator with proper symbols using Rich Panel."""
        from rich.panel import Panel

        # Build progress text with colors. Read names from
        # WorkflowStateManager so plugin-contributed stages (added
        # via register_plugin_stage) render with their declared
        # display label instead of IndexError-ing on a hardcoded list.
        stages = self.workflow_state.get_all_stages()

        progress_parts = []
        for stage in stages:
            stage_name = self.workflow_state.STAGE_NAMES.get(stage, stage)
            symbol = self.workflow_state.get_stage_symbol(stage)
            color = self.workflow_state.get_stage_color(stage)
            # Bold so the lighter status colors (green/yellow) keep
            # enough contrast to read on a light/white terminal too.
            progress_parts.append(f"[bold {color}]{stage_name}:{symbol}[/bold {color}]")

        progress_line = " → ".join(progress_parts)

        # Generate smart suggestion
        suggestion = self._get_smart_suggestion_v2(workspace)

        # Build the guide content
        guide_content = f"{progress_line}\n\n{suggestion}"

        # Fixed-width panel so the suggestion text wraps inside a stable
        # box instead of stretching the panel to the full terminal width
        # (expand=False would size to the longest line, i.e. the
        # single-line suggestion). 76 cols comfortably holds the
        # progress line on one row while wrapping the longer suggestion.
        self.console.print(Panel(
            guide_content,
            title="WORKFLOW GUIDE",
            border_style="bright_blue",
            width=76,
            expand=False
        ))

    def _get_smart_suggestion_v2(self, workspace):
        """Generate context-aware suggestions based on workflow state."""
        from proprep.utils.structure_selector import StructureSelector

        current_stage = self.workflow_state.get_current_stage()
        selector = StructureSelector(workspace, self.processor.console)
        status = selector.get_structure_status()
        has_structure = status.get("has_any", False)

        # Return plain text without markup - styling will be applied at display time
        suggestions = {
            "load": "Start by loading a structure from databases or local files and then optionally applying symmetry transformations to generate a biological assembly.",
            "detect": "Detect redox-active sites, explore the parameter space of AMBER FFs, and if needed, parameterize small organic molecules, modified amino acids, or metal sites.",
            "compare": "Compare sequences or align structures",
            "fix": "Filter, mutate, repair, prepare redox sites, assign protonation states, orient your structure, or embed it in a membrane",
            "simulate": "Build topologies, configure, execute, or analyze MD simulations, or setup QM/MM calculations on MD snapshots.",
        }

        return suggestions.get(current_stage, "")

    def _get_user_choice(self) -> str:
        """Get user menu choice (case-insensitive)."""
        while True:
            response = prompt_with_context(
                processor=self.processor,
                prompt="\nEnter your choice",
                default="1",
                module="Main Menu (Workflow Mode)",
                description="Select option",
                options_map=self.options_map
            )
            normalized = response.lower()
            if normalized in self.options_map:
                return normalized
            else:
                self.console.print(f"[red]Invalid choice: {response}. Please try again.[/red]")

    def _handle_choice(self, choice: str) -> bool:
        """Handle user menu choice. Returns True to continue, False to exit."""
        selected = self.options_map[choice]

        # Plugin tools are encoded as ``__plugin::<name>`` sentinels.
        # Dispatch them before the generic module-selection path so a
        # plugin name can never collide with a registered module.
        if selected.startswith("__plugin::"):
            self._dispatch_plugin(selected.removeprefix("__plugin::"))
            return True

        if selected == "exit":
            if confirm_with_context(
                self.processor,
                "Are you sure you want to exit?",
                default=False,
                module="Workflow Menu",
                description="Confirm exit"
            ):
                self.console.print("[green]Goodbye![/green]")
                import sys
                sys.exit(0)
            return True
        elif selected == "utilities":
            self._show_utilities_submenu()
        elif selected == "back":
            self.workflow_state.go_back()
            self.console.print(f"[yellow]← Moved back to {self.workflow_state.STAGE_NAMES[self.workflow_state.get_current_stage()]} stage[/yellow]")
        elif selected == "next":
            old_stage = self.workflow_state.get_current_stage()
            self.workflow_state.mark_completed(old_stage)
            self.workflow_state.advance_stage()
            new_stage = self.workflow_state.get_current_stage()
            self.console.print(f"[green]→ Moved to {self.workflow_state.STAGE_NAMES[new_stage]} stage[/green]")
        elif selected.startswith("jump:"):
            target_stage = selected.split(":")[1]
            target_name = self.workflow_state.STAGE_NAMES[target_stage]
            if target_stage == self.workflow_state.get_current_stage():
                self.console.print(f"[grey50]Already on {target_name} stage[/grey50]")
            else:
                self.workflow_state.jump_to_stage(target_stage)
                self.console.print(f"[bright_blue]Jumped to {target_name} stage[/bright_blue]")
        elif selected == "forcefield_explorer":
            from .forcefield_explorer import ForceFieldExplorerCommand
            explorer = ForceFieldExplorerCommand(self.processor)
            explorer.execute()
        elif selected == "all":
            # Switch to full menu mode
            self.processor.workspace.set("menu_mode", "full-menu")
            self.console.print("[bright_blue]Switching to Full Menu mode...[/bright_blue]")
            return False  # Exit to trigger menu reload
        else:
            # Handle module selection
            # Auto-advance only after completing a tool in Load stage
            current_stage = self.workflow_state.get_current_stage()

            module_menu = RunModuleMenuCommand(self.processor, selected, "Workflow Menu")
            module_menu.execute()

            # Auto-advance after Load stage if structure was loaded
            from proprep.utils.structure_selector import StructureSelector
            selector = StructureSelector(self.workspace, self.console)
            status = selector.get_structure_status()
            if current_stage == "load" and status.get("has_any", False):
                self.workflow_state.mark_completed("load")
                self.workflow_state.advance_stage()
                self.console.print("[green]✓ Structure loaded! Advanced to Detect stage[/green]")

        return True

    def _dispatch_plugin(self, plugin_name: str) -> None:
        """Thin wrapper around the shared module-level dispatcher
        so workflow-mode and full-menu mode handle plugins identically."""
        _dispatch_plugin_in_menu(self.processor, self.console, plugin_name)

    def _show_utilities_submenu(self):
        """Display utilities submenu with workspace, preferences, help, and feedback options."""
        self.console.print("\n[bold bright_blue]═══ Utilities ═══[/bold bright_blue]")
        self.console.print("  w. Workspace Options")
        self.console.print("  p. Preferences")
        self.console.print("  h. Help & Workflows")
        self.console.print("  f. Feedback & Support")
        self.console.print("  c. Cancel")

        choice = prompt_with_context(
            self.processor,
            "\nSelect option",
            choices=["w", "p", "h", "f", "c"],
            default="c",
            module="Utilities Menu",
            description="Select utility option",
            options_map={"w": "Workspace Options", "p": "Preferences", "h": "Help", "f": "Feedback", "c": "Cancel"}
        )

        if choice == "w":
            workspace_menu = WorkspaceOptionsMenuCommand(self.processor)
            workspace_menu.execute()
        elif choice == "p":
            preferences_menu = PreferencesMenuCommand(self.processor)
            preferences_menu.execute()
        elif choice == "h":
            self._show_help_and_workflows()
        elif choice == "f":
            from .feedback_command import FeedbackCommand
            feedback_cmd = FeedbackCommand(self.processor)
            feedback_cmd.execute()
        # 'c' or any other choice just returns to main menu

    def _show_help_and_workflows(self):
        """Display help information and common workflows."""
        self.console.print("\n[bold bright_blue]═══════════════════════════════════════════════════════════════[/bold bright_blue]")
        self.console.print("[bold bright_blue]                    Workflow Mode Help[/bold bright_blue]")
        self.console.print("[bold bright_blue]═══════════════════════════════════════════════════════════════[/bold bright_blue]\n")

        self.console.print("[bold]How Workflow Mode Works:[/bold]")
        self.console.print("  • Shows only tools relevant to your current stage")
        self.console.print("  • Guides you step-by-step through the preparation process")
        self.console.print("  • Use 'n' to advance, 'b' to go back, 's' to skip optional stages")
        self.console.print("  • Use 'a' to switch to Full Menu mode if you need all tools\n")

        self.console.print("[bold]Workflow Stages:[/bold]")
        stages_info = [
            ("Load", "Load PDB structure from database or file"),
            ("Detect", "Optional: Detect redox-active sites for alignment or analysis"),
            ("Compare", "Optional: Sequence/structure comparison and alignment"),
            ("Fix", "Filter, repair, and prepare structure"),
            ("Setup", "Specialized setup: redox sites, restraints, force fields, topology"),
            ("Simulate", "Run MD simulations or QM/MM calculations"),
        ]
        for stage, desc in stages_info:
            self.console.print(f"  [bright_blue]{stage}:[/bright_blue] {desc}")

        self.console.print("\n[bold]Common Workflows:[/bold]")
        workflows = [
            ("Basic Protein MD", "Load → Fix → Setup (tLEaP) → Simulate (MD)"),
            ("Metalloprotein MD", "Load → Detect → Fix → Setup (RedoxSite, FF-Params, tLEaP) → Simulate (MD)"),
            ("QM/MM Calculation", "Load → Detect → Fix → Setup (RedoxSite) → tLEaP → QM/MM"),
            ("With Alignment", "Load → Detect → Compare (Align on redox sites) → Fix → ..."),
        ]
        for title, workflow in workflows:
            self.console.print(f"  [bold]{title}:[/bold] {workflow}")

        self.console.print("\n[grey50]Press Enter to return to menu...[/grey50]")
        prompt_with_context(
            self.processor,
            "",
            default="",
            module="Workflow Menu",
            description="Pause after workflow overview",
        )


class WorkspaceOptionsMenuCommand(MenuCommand):
    """Command to display and handle workspace options menu."""

    def __init__(self, processor: Processor):
        super().__init__(processor, "Workspace Options", "Workspace options menu")

    def execute(self) -> None:
        """Display and handle workspace options menu."""
        self.enter_menu()

        try:
            while True:
                self.show_breadcrumbs()
                self._show_menu()

                choice = self._get_user_choice()
                if not self._handle_choice(choice):
                    break

        finally:
            self.exit_menu()

    def _show_menu(self):
        """Display workspace options menu."""
        self.console.print("\n[bold blue]══ Workspace Options ══[/bold blue]")

        workspace = self.workspace
        debug_enabled = workspace.get("debug", False)

        self.options = {
            "1": ("Show workspace status", ShowWorkspaceStatusCommand(self.processor)),
            "2": (
                "Show detailed workspace",
                ShowDetailedWorkspaceCommand(self.processor),
            ),
            "3": ("Save workspace", SaveWorkspaceCommand(self.processor)),
            "4": ("Load workspace", LoadWorkspaceCommand(self.processor)),
            "5": ("Reset workspace", ResetWorkspaceCommand(self.processor)),
            "6": (
                f"Toggle debug mode (currently {'ON' if debug_enabled else 'OFF'})",
                ToggleDebugCommand(self.processor),
            ),
            "7": ("Show command history", ShowWorkspaceHistoryCommand(self.processor)),
            "8": ("Save command history", SaveCommandHistoryCommand(self.processor)),
            "9": ("Load command history", LoadCommandHistoryCommand(self.processor)),
            "10": (
                "Replay command history",
                ReplayCommandHistoryCommand(self.processor),
            ),
            "11": (
                "Reset command history",
                ResetCommandHistoryCommand(self.processor),
            ),
            "12": ("Back to main menu", None),
        }

        # Display options
        for key, (description, _) in self.options.items():
            self.console.print(f"{key}. {description}")

    def _get_user_choice(self) -> str:
        """Get user menu choice."""
        # Create options_map for rich context
        options_map = {key: description for key, (description, _) in self.options.items()}

        return prompt_with_context(
            processor=self.processor,
            prompt="\nEnter your choice",
            choices=list(self.options.keys()),
            default="1",
            module="Workspace Options Menu",
            description="Select workspace option",
            options_map=options_map
        )

    def _handle_choice(self, choice: str) -> bool:
        """Handle user menu choice. Returns True to continue, False to exit."""
        description, command = self.options[choice]

        if command is None:  # Back to main menu
            return False

        try:
            command.execute_with_error_handling()
        except Exception:
            # Error already handled in execute_with_error_handling
            pass

        return True


class PreferencesMenuCommand(MenuCommand):
    """Command to display and handle user preferences menu."""

    def __init__(self, processor: Processor):
        super().__init__(processor, "Preferences", "User preferences menu")

    def execute(self) -> None:
        """Display and handle preferences menu."""
        self.enter_menu()

        try:
            while True:
                self.show_breadcrumbs()
                self._show_menu()

                choice = self._get_user_choice()
                if not self._handle_choice(choice):
                    break

        finally:
            self.exit_menu()

    def _show_menu(self):
        """Display preferences menu."""
        from proprep.utils.settings_manager import SettingsManager

        self.console.print("\n[bold blue]══ Preferences ══[/bold blue]")

        settings_mgr = SettingsManager()
        current_mode = settings_mgr.get_menu_mode()
        current_layout = settings_mgr.get_menu_layout()

        # Map internal values to display names
        mode_display = {
            "workflow": "Workflow Mode (guided)",
            "full-menu": "Full Menu Mode (all tools)",
            "ask": "Ask each time"
        }
        layout_display = {
            "list": "List (single column)",
            "grid": "Grid (2-column panels)",
        }

        # Check GitHub configuration status
        github_token = settings_mgr.get_github_token()
        github_repo = settings_mgr.get_github_repo()
        github_status = "configured" if (github_token and github_repo) else "not configured"

        self.options = {
            "1": (
                f"Change default menu mode (currently: {mode_display[current_mode]})",
                "change_menu_mode"
            ),
            "2": (
                f"Change full-menu layout (currently: {layout_display[current_layout]})",
                "change_menu_layout"
            ),
            "3": (f"Configure GitHub feedback (currently: {github_status})", "configure_github"),
            "4": ("Reset all preferences to defaults", "reset_preferences"),
            "5": ("Back to main menu", None),
        }

        # Display options
        for key, (description, _) in self.options.items():
            self.console.print(f"{key}. {description}")

    def _get_user_choice(self) -> str:
        """Get user menu choice."""
        # Create options_map for rich context
        options_map = {key: description for key, (description, _) in self.options.items()}

        return prompt_with_context(
            processor=self.processor,
            prompt="\nEnter your choice",
            choices=list(self.options.keys()),
            default="1",
            module="Preferences Menu",
            description="Select preference option",
            options_map=options_map
        )

    def _handle_choice(self, choice: str) -> bool:
        """Handle user menu choice. Returns True to continue, False to exit."""
        from proprep.utils.settings_manager import SettingsManager
        from proprep.utils.prompts import prompt_with_context

        description, action = self.options[choice]

        if action is None:  # Back to main menu
            return False

        try:
            if action == "change_menu_mode":
                self._change_menu_mode()
            elif action == "change_menu_layout":
                self._change_menu_layout()
            elif action == "configure_github":
                self._configure_github()
            elif action == "reset_preferences":
                self._reset_preferences()
        except Exception as e:
            self.console.print(f"[red]Error: {e}[/red]")

        return True

    def _change_menu_mode(self):
        """Change the default menu mode preference."""
        from proprep.utils.settings_manager import SettingsManager

        settings_mgr = SettingsManager()
        current_mode = settings_mgr.get_menu_mode()

        self.console.print("\n[bold bright_blue]Select Default Menu Mode[/bold bright_blue]")
        self.console.print("  1. Workflow Mode (guided, step-by-step)")
        self.console.print("  2. Full Menu Mode (all tools visible)")
        self.console.print("  3. Ask each time")

        mode_map = {
            "1": "workflow",
            "2": "full-menu",
            "3": "ask"
        }

        choice = prompt_with_context(
            self.processor,
            "\nSelect mode",
            choices=["1", "2", "3"],
            default="1",
            module="Preferences",
            description="Select default menu mode",
            options_map={"1": "Workflow Mode", "2": "Full Menu Mode", "3": "Ask each time"}
        )

        new_mode = mode_map[choice]

        if new_mode != current_mode:
            settings_mgr.set_menu_mode(new_mode)
            self.console.print(f"[green]✓ Default menu mode set to: {new_mode}[/green]")
            self.console.print("[yellow]Note: Restart ProPrep for the change to take effect[/yellow]")
        else:
            self.console.print("[grey50]Mode unchanged[/grey50]")

    def _change_menu_layout(self):
        """Change the full-menu layout preference (list vs panel grid)."""
        from proprep.utils.settings_manager import SettingsManager

        settings_mgr = SettingsManager()
        current_layout = settings_mgr.get_menu_layout()

        self.console.print("\n[bold bright_blue]Select Full-Menu Layout[/bold bright_blue]")
        self.console.print("  1. List (single column, one tool per row)")
        self.console.print("  2. Grid (2-column panels, good for screenshots)")

        layout_map = {"1": "list", "2": "grid"}

        choice = prompt_with_context(
            self.processor,
            "\nSelect layout",
            choices=["1", "2"],
            default="1",
            module="Preferences",
            description="Select full-menu layout",
            options_map={"1": "List", "2": "Grid"}
        )

        new_layout = layout_map[choice]

        if new_layout != current_layout:
            settings_mgr.set_menu_layout(new_layout)
            self.console.print(f"[green]✓ Full-menu layout set to: {new_layout}[/green]")
        else:
            self.console.print("[grey50]Layout unchanged[/grey50]")

    def _configure_github(self):
        """Configure GitHub integration for feedback submission."""
        from proprep.utils.settings_manager import SettingsManager
        from rich.panel import Panel

        settings_mgr = SettingsManager()

        # Show introduction
        intro_text = """[bold bright_blue]GitHub Feedback Configuration[/bold bright_blue]

To submit feedback directly to GitHub Issues, you need:

1. [bold]GitHub Personal Access Token (PAT)[/bold]
   - Classic token with 'repo' or 'public_repo' scope
   - Create at: https://github.com/settings/tokens

2. [bold]Repository[/bold]
   - Format: owner/repo (e.g., "Mag14011/ProPrep")
   - Must have Issues enabled

[grey50]Your token is stored locally in ~/.proprep/settings.json
Not sharing with anyone. For feedback only.[/grey50]"""

        self.console.print(Panel(intro_text, border_style="bright_blue"))
        self.console.print()

        # Get current settings
        current_repo = settings_mgr.get_github_repo()
        current_token = settings_mgr.get_github_token()

        if current_repo:
            self.console.print(f"[grey50]Current repository: {current_repo}[/grey50]")
        if current_token:
            self.console.print(f"[grey50]Current token: {'*' * 20} (hidden)[/grey50]")
        self.console.print()

        # Configure repository
        repo = prompt_with_context(
            self.processor,
            "GitHub repository (owner/repo)",
            default=current_repo if current_repo else "Mag14011/ProPrep",
            module="Preferences",
            description="Configure GitHub repository"
        )

        # Configure token
        self.console.print("\n[yellow]Enter your GitHub Personal Access Token[/yellow]")
        self.console.print("[grey50](Input will be hidden)[/grey50]")

        # Use getpass for secure input
        import getpass
        token = getpass.getpass("Token: ")

        if not token and current_token:
            self.console.print("[grey50]Keeping existing token[/grey50]")
            token = current_token

        # Validate format
        if not repo or "/" not in repo:
            self.console.print("[red]Invalid repository format. Expected: owner/repo[/red]")
            return

        if not token:
            self.console.print("[red]Token is required[/red]")
            return

        # Save settings
        settings_mgr.set_github_repo(repo)
        settings_mgr.set_github_token(token)

        self.console.print("\n[green]✓ GitHub feedback configuration saved[/green]")
        self.console.print(f"[green]  Repository: {repo}[/green]")
        self.console.print("[grey50]  Token: (hidden)[/grey50]")
        self.console.print("\n[bright_blue]You can now use the Feedback menu to submit issues![/bright_blue]")

    def _reset_preferences(self):
        """Reset all preferences to defaults."""
        from proprep.utils.settings_manager import SettingsManager

        if confirm_with_context(
            self.processor,
            "\n[yellow]This will reset all preferences to defaults. Continue?[/yellow]",
            default=False,
            module="Preferences",
            description="Confirm reset preferences"
        ):
            settings_mgr = SettingsManager()
            settings_mgr.reset_to_defaults()
            self.console.print("[green]✓ All preferences reset to defaults[/green]")
        else:
            self.console.print("[grey50]Reset cancelled[/grey50]")


class ExitCommand(MenuCommand):
    """Command to exit the application."""

    def __init__(self, processor: Processor):
        super().__init__(processor, "Exit", "Exit application")

    def execute(self) -> bool:
        """Exit the application with confirmation."""
        if confirm_with_context(
            self.processor,
            "Are you sure you want to exit?",
            default=False,
            module="Exit",
            description="Confirm exit"
        ):
            self.console.print("[green]Goodbye![/green]")
            return True
        return False
