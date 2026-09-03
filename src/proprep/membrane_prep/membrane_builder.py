"""
Membrane Builder Module

Builds membrane-protein systems using packmol-memgen. Provides a comprehensive
UI for configuring lipid composition, geometry, ions, and packing parameters.

Parametrization (tLEaP) and minimization are delegated to ProPrep's existing
Topology Generator and MD Manager modules, respectively. This ensures
ProPrep retains full control over force field sourcing, bond directives, and
custom parameters for redox-active proteins.
"""

import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from proprep.utils.module_registry import ProcessingModule, register_module
from proprep.utils.prompts import (
    confirm_with_context,
    float_prompt_with_context,
    int_prompt_with_context,
    prompt_float_with_retry,
    prompt_int_with_retry,
    prompt_with_context,
)

from .lipid_library import LipidLibrary
from .membrane_config import MembraneConfig, SoluteConfig

logger = logging.getLogger(__name__)

MODULE_NAME = "Membrane Builder"


@register_module
class MembraneBuilderModule(ProcessingModule):
    """Module for building membrane-protein systems via packmol-memgen."""

    NAME = MODULE_NAME
    DESCRIPTION = "Build membrane-protein systems using packmol-memgen"
    VERSION = "1.0.0"
    CATEGORY = "preparation"
    PRIORITY = 55  # After Protonation State Analyzer (53) and Structure Orientator (50)

    def __init__(self):
        super().__init__()
        self.config = MembraneConfig()
        self.lipid_library = LipidLibrary()
        self._last_result = None

    @property
    def console(self) -> Console:
        if self.processor and hasattr(self.processor, "console"):
            return self.processor.console
        return Console()

    # ── Workspace interface ──────────────────────────────────────────────

    def get_workspace_requirements(self) -> list:
        # No hard requirements — can build protein-free bilayers
        return []

    def get_workspace_outputs(self) -> list:
        return [
            "membrane_packed_pdb",
            "membrane_config",
            "membrane_leaprc_requirements",
            "membrane_box_dimensions",
            "membrane_ion_summary",
            "membrane_solutes",
            "is_membrane_system",
        ]

    def can_process(self, workspace) -> bool:
        # Always available — protein is optional
        return True

    def get_menu_options(self) -> Dict[str, str]:
        return {
            "build_membrane": "Build membrane-protein system",
        }

    def get_enhanced_menu_options(self, workspace):
        from proprep.utils.enhanced_menu import MenuOption, OptionStatus

        options = []

        has_membrane = workspace.get("membrane_packed_pdb") is not None

        if has_membrane:
            status = OptionStatus.COMPLETED
            desc = "Build membrane-protein system"
        else:
            status = OptionStatus.AVAILABLE
            desc = "Build membrane-protein system"

        options.append(
            MenuOption(
                key="1",
                description=desc,
                status=status,
            )
        )

        return options

    def handle_menu_option(self, option: str) -> bool:
        if option == "build_membrane":
            return self.process(self.processor.workspace)
        return False

    # ── Main entry point ─────────────────────────────────────────────────

    def process(self, workspace) -> bool:
        """Run the membrane builder interactive workflow."""
        self.config = MembraneConfig()
        self._auto_configure_from_workspace(workspace)

        while True:
            self._show_main_menu(workspace)

            choice = prompt_with_context(
                self.processor,
                "\nEnter choice",
                choices=["1", "2", "3", "4", "5", "6", "7", "8", "9", "r", "g", "b"],
                default="2",
                module=MODULE_NAME,
                description="Main menu selection",
                options_map={
                    "1": "Protein Selection",
                    "2": "Lipid Composition",
                    "3": "Solvent & Ions",
                    "4": "Box & Membrane Dimensions",
                    "5": "Protein Orientation",
                    "6": "Specialized Geometry",
                    "7": "PACKMOL Settings",
                    "8": "Force Fields",
                    "9": "Custom Parameters & Solutes",
                    "r": "Review Full Configuration",
                    "g": "Generate & Run",
                    "b": "Back",
                },
            )

            if choice == "b":
                return False
            elif choice == "1":
                self._menu_protein_selection(workspace)
            elif choice == "2":
                self._menu_lipid_composition()
            elif choice == "3":
                self._menu_solvent_ions()
            elif choice == "4":
                self._menu_box_dimensions()
            elif choice == "5":
                self._menu_protein_orientation()
            elif choice == "6":
                self._menu_specialized_geometry()
            elif choice == "7":
                self._menu_packmol_settings()
            elif choice == "8":
                self._menu_force_fields(workspace)
            elif choice == "9":
                self._menu_custom_parameters()
            elif choice == "r":
                self._show_review(workspace)
            elif choice == "g":
                if not self.config.lipids and not self.config.solvate_only:
                    self.console.print(
                        "[red]Lipid composition is required. Configure it first "
                        "(option 2), or choose a solvate-only build there.[/red]"
                    )
                    continue
                return self._run_build(workspace)

    # ── Auto-configuration from workspace ────────────────────────────────

    def _auto_configure_from_workspace(self, workspace):
        """Pre-populate config from workspace state."""
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(workspace, self.console)
        pdb_path = selector.get_structure(silent=True)
        if pdb_path:
            self.config.protein_pdb = pdb_path

        # Protonation: always use tLEaP to add hydrogens, bypassing reduce.
        # tLEaP correctly handles both standard residues and custom forcefield
        # residues (e.g., HBO, FHO from redox site transformers), and respects
        # protonation state assignments (HID/HIE/HIP, ASH, GLH) from ProPrep's
        # Protonation State Analyzer.
        self.config.skip_protonation = True

        # Auto-compute charge delta from redox sites
        # TODO: This feature is currently non-functional. RedoxSite objects do not
        # carry net_formal_charge or formal_charge attributes, so the delta is always 0
        # and --charge_pdb_delta is never passed to packmol-memgen. The formal charge
        # per redox state IS available in metadata.json (e.g., reduced: 0, oxidized: +1
        # for b-type heme) but is not propagated to RedoxSite or transformer_info during
        # transformation. In practice this has no effect on the final result because
        # the tLEaP topology generation step performs an addionsrand neutralization pass
        # that corrects any charge imbalance from packmol-memgen's ion placement.
        # To fix: read formal_charge from forcefield metadata during transformation
        # and store it on the RedoxSite object or in transformer_info.
        redox_sites = workspace.get("detected_redox_sites")
        if redox_sites:
            self.config.auto_charge_delta = self._compute_redox_charge_delta(redox_sites)

    def _compute_redox_charge_delta(self, redox_sites) -> int:
        """
        Compute the charge contribution from non-standard residues that
        packmol-memgen's hardcoded dictionary won't recognize.

        NOTE: Currently returns 0 for all inputs because RedoxSite objects
        lack formal_charge attributes. See TODO above.
        """
        delta = 0
        for site in redox_sites:
            if hasattr(site, "net_formal_charge"):
                delta += site.net_formal_charge
            elif hasattr(site, "formal_charge"):
                delta += site.formal_charge
            elif isinstance(site, dict):
                delta += site.get("net_formal_charge", site.get("formal_charge", 0))
        return delta

    # ── Main menu display ────────────────────────────────────────────────

    def _show_main_menu(self, workspace):
        """Display the top-level membrane builder menu."""

        self.console.print(
            Panel(
                "[bold black]Build membrane-protein systems using packmol-memgen.[/bold black]\n"
                "This tool orients your protein in a lipid bilayer, packs "
                "lipids/water/ions around it, and produces a packed PDB. "
                "Parametrization with tLEaP and minimization are handled by "
                "ProPrep's downstream modules.",
                title="Membrane Builder",
                border_style="bright_blue",
                width=58,
                padding=(0, 1),
            )
        )

        # Dashboard
        protein_str = Path(self.config.protein_pdb).name if self.config.protein_pdb else "[grey50]none[/grey50]"
        lipid_str = self.config.lipids or "[dark_orange3]not set[/dark_orange3]"
        water_str = self.config.effective_water_model.upper()
        salt_str = f"{self.config.salt_concentration} M {self.config.cation}/{self.config.anion}" if self.config.salt else "none"
        orient_str = "Pre-oriented" if self.config.preoriented else f"Auto ({self.config.orientation_method})"

        # Rich's automatic highlighter (highlight=True by default) recolors bare
        # numbers (cyan), name=value tokens (yellow), and symbol patterns like
        # "Cl-" (magenta) on its own, which reads as arbitrary. Build the whole
        # dashboard as one block and print it with highlight=False so the ONLY
        # colors are the explicit markup below: bright-blue (#1f6feb) for the
        # option keys, dark green (#1a7f37) for the run action, dark_orange3 for
        # "not set"/"Required". All values stay default text. Chosen to stay
        # legible on both white (manuscript) and black backgrounds.
        BLUE = "#1f6feb"
        if self.config.solvate_only:
            lipid_indicator = "Solvate only (no lipids)"
        elif self.config.lipids:
            lipid_indicator = lipid_str
        else:
            lipid_indicator = "[dark_orange3]Required[/dark_orange3]"
        solute_count = len(self.config.solutes)
        solute_str = f"{solute_count} solute(s)" if solute_count else "none"

        if hasattr(self, '_selected_forcefields') and self._selected_forcefields:
            ff_parts = []
            for cat in ['protein', 'lipids', 'water']:
                sel = self._selected_forcefields.get(cat)
                if sel and isinstance(sel, dict):
                    ff_parts.append(sel['name'])
            extras = sum(1 for cat in ['modified_aa', 'dna', 'rna', 'carbohydrates', 'small_molecules', 'ions']
                        if self._selected_forcefields.get(cat) not in (None, []))
            ff_str = " / ".join(ff_parts)
            if extras:
                ff_str += f" (+{extras} more)"
        else:
            ff_str = f"{self.config.ffprot} / {self.config.fflip}"

        lines = [
            "  [bold black]Current Configuration[/bold black]",
            f"  Protein:       {protein_str}",
            f"  Lipids:        {lipid_str}",
            f"  Solvent:       {water_str}",
            f"  Salt:          {salt_str}",
            f"  Orientation:   {orient_str}",
            f"  Box padding:   {self.config.box_padding} Å",
        ]
        if self.config.auto_charge_delta != 0:
            lines.append(
                f"  Charge corr.:  {self.config.auto_charge_delta:+d} "
                f"(auto-detected from redox sites)"
            )
        lines += [
            "",
            "  [bold black]SYSTEM[/bold black]",
            f"    [{BLUE}]1[/{BLUE}]   {'Protein Selection':<25}{protein_str}",
            f"    [{BLUE}]2[/{BLUE}]   {'Lipid Composition':<25}{lipid_indicator}",
            f"    [{BLUE}]3[/{BLUE}]   {'Solvent & Ions':<25}{water_str}, {salt_str}",
            "",
            "  [bold black]GEOMETRY[/bold black]",
            f"    [{BLUE}]4[/{BLUE}]   {'Box & Membrane Dims':<25}{self.config.box_padding} Å pad, {self.config.water_layer} Å water",
            f"    [{BLUE}]5[/{BLUE}]   {'Protein Orientation':<25}{orient_str}",
            f"    [{BLUE}]6[/{BLUE}]   {'Specialized Geometry':<25}Curvature, CompEL...",
            "",
            "  [bold black]PACKING[/bold black]",
            f"    [{BLUE}]7[/{BLUE}]   {'PACKMOL Settings':<25}tol={self.config.tolerance}, nloop={self.config.nloop}/{self.config.nloop_all}",
            "",
            "  [bold black]PARAMETRIZATION[/bold black]",
            f"    [{BLUE}]8[/{BLUE}]   {'Force Fields':<25}{ff_str}",
            f"    [{BLUE}]9[/{BLUE}]   {'Custom Params & Solutes':<25}Solutes: {solute_str}",
            "",
            f"    [{BLUE}]r[/{BLUE}]   Review Full Configuration",
            "    [bold #1a7f37]g   Generate & Run[/bold #1a7f37]",
            f"    [{BLUE}]b[/{BLUE}]   Back",
        ]
        self.console.print("\n".join(lines), highlight=False)

    # ── Option 1: Protein Selection ──────────────────────────────────────

    def _menu_protein_selection(self, workspace):
        """Configure protein input."""
        self.console.print(
            Panel(
                "The protein structure will be embedded in the membrane bilayer.\n"
                "If your protein is already loaded in the ProPrep workspace, you can\n"
                "use it directly. You can also build a protein-free bilayer.",
                title="Protein Selection",
                border_style="bright_blue",
                expand=False,
            )
        )

        current = Path(self.config.protein_pdb).name if self.config.protein_pdb else "none"
        self.console.print(f"  Current: {current}\n")

        from proprep.utils.structure_selector import StructureSelector
        selector = StructureSelector(workspace, self.console)
        ws_structure = selector.get_structure(silent=True)

        choices = ["1", "2", "3", "b"]
        if ws_structure:
            self.console.print(f"    1   Use workspace structure       {Path(ws_structure).name}")
        else:
            self.console.print("    1   [grey50]Use workspace structure       (none loaded)[/grey50]")
        self.console.print("    2   Load a different PDB file")
        self.console.print("    3   No protein (empty bilayer)")
        self.console.print("    b   Back")

        choice = prompt_with_context(
            self.processor,
            "\nSelect",
            choices=choices,
            default="1" if ws_structure else "3",
            module=MODULE_NAME,
            description="Protein source selection",
            options_map={
                "1": "Use workspace structure",
                "2": "Load a different PDB file",
                "3": "No protein (empty bilayer)",
                "b": "Back",
            },
        )

        if choice == "1" and ws_structure:
            self.config.protein_pdb = ws_structure
            self.console.print(f"[green]Using workspace structure: {Path(ws_structure).name}[/green]")
            self._show_protein_in_viewer()
        elif choice == "2":
            from proprep.structure_prep.pdb_loader import display_pdb_file_menu

            start_dir = "."
            if self.config.protein_pdb:
                start_dir = str(Path(self.config.protein_pdb).parent)

            path = display_pdb_file_menu(
                directory=start_dir,
                console=self.console,
                processor=self.processor,
            )
            if path and Path(path).exists():
                self.config.protein_pdb = str(Path(path).resolve())
                self.console.print(f"[green]Loaded: {Path(path).name}[/green]")
                self._show_protein_in_viewer()
            elif path:
                self.console.print(f"[red]File not found: {path}[/red]")
            else:
                self.console.print("[grey50]Canceled.[/grey50]")
        elif choice == "3":
            self.config.protein_pdb = None
            self.console.print("[green]No protein — building empty bilayer.[/green]")
            # packmol-memgen sizes the membrane patch around the protein; with
            # no protein it requires --distxy_fix. Prompt for it now so the
            # user can't drift into Generate & Run with an undersized config.
            if self.config.fixed_xy is None:
                self.console.print(
                    "[grey50]Empty bilayers need an explicit x/y patch size "
                    "(packmol-memgen --distxy_fix).[/grey50]"
                )
                self.config.fixed_xy = prompt_float_with_retry(
                    self.processor,
                    "Membrane patch size (Å, x and y)",
                    default=70.0,
                    module=MODULE_NAME,
                    description="Empty-bilayer patch size",
                )

    def _show_protein_in_viewer(self) -> None:
        """Show ``self.config.protein_pdb`` in the viewer if set.

        Best-effort — the viewer is informational only here, so any
        failure is silently dropped. Used by Protein Selection,
        Review, and post-build hooks.
        """
        protein = getattr(self.config, "protein_pdb", None)
        if not protein:
            return
        try:
            from proprep.structure_prep.viewer_coordinator import (
                viewer as _viewer,
            )
            _viewer.show_structure(protein)
        except Exception as exc:
            logger.debug("membrane viewer hook silenced: %s", exc)

    # ── Option 2: Lipid Composition ──────────────────────────────────────

    def _menu_lipid_composition(self):
        """Configure lipid composition."""
        self.console.print(
            Panel(
                "Define the lipid composition of your membrane. You can use a single\n"
                "lipid type or a mixture. For asymmetric membranes, define upper and\n"
                "lower leaflets separately. Lipid names are packmol-memgen's; browse\n"
                "or search its database (option 2) to see what is available.",
                title="Lipid Composition",
                border_style="bright_blue",
                expand=False,
            )
        )

        if self.config.solvate_only:
            self.console.print("  Current: solvate only (no lipids, no membrane)\n")
        else:
            current = self.config.lipids or "not set"
            ratio = f"  ratio: {self.config.ratio}" if self.config.ratio else ""
            self.console.print(f"  Current: {current}{ratio}\n")

        self.console.print("    1   Build custom composition")
        self.console.print("    2   Search / browse lipid library")
        self.console.print("    3   Enter raw lipid string (packmol-memgen syntax)")
        self.console.print("    4   No lipids (solvate only, no membrane)")
        self.console.print("    b   Back")

        choice = prompt_with_context(
            self.processor,
            "\nSelect",
            choices=["1", "2", "3", "4", "b"],
            default="1",
            module=MODULE_NAME,
            description="Lipid composition method",
            options_map={
                "1": "Build custom composition",
                "2": "Search / browse lipid library",
                "3": "Enter raw lipid string (packmol-memgen syntax)",
                "4": "No lipids (solvate only, no membrane)",
                "b": "Back",
            },
        )

        if choice == "1":
            self._lipid_custom_builder()
        elif choice == "2":
            self._lipid_browser()
        elif choice == "3":
            self._lipid_raw_string()
        elif choice == "4":
            self._lipid_solvate_only()
            return

        # Choosing an actual lipid composition turns solvate-only back off: the
        # two are mutually exclusive, since a solvated box has no bilayer.
        if self.config.lipids:
            self.config.solvate_only = False

    def _lipid_solvate_only(self):
        """Configure a solvent-only build (packmol-memgen ``--solvate``, no lipids).

        This is the no-lipid case of a membrane build: packmol-memgen packs a
        water and ion box around the solute with no bilayer, using the same
        engine and box conventions as a membrane system. It is mutually
        exclusive with a lipid composition, so any lipids already set are
        cleared when it is enabled.
        """
        self.console.print(
            Panel(
                "Builds a solvated system with no lipid bilayer. packmol-memgen packs\n"
                "only water and ions around the solute, using the same engine and box\n"
                "conventions as a membrane build. Useful for preparing a matched\n"
                "non-membrane reference of the same system.",
                title="Solvate Only (no membrane)",
                border_style="bright_blue",
                expand=False,
            )
        )
        enable = confirm_with_context(
            self.processor, "Build a solvent-only system (no lipids)?",
            default=not self.config.solvate_only,
            module=MODULE_NAME, description="Solvate only (no membrane)",
        )
        self.config.solvate_only = enable
        if enable:
            self.config.lipids = None
            self.config.ratio = None
            self.config.cubic = confirm_with_context(
                self.processor, "Use a cubic box?", default=self.config.cubic,
                module=MODULE_NAME, description="Cubic box",
            )
            self.console.print("[green]Solvate-only mode enabled (no lipids).[/green]")
        else:
            self.console.print("[grey50]Solvate-only mode disabled.[/grey50]")

    def _lipid_custom_builder(self):
        """Interactively build a lipid composition."""
        self.console.print("\n[bold]Build custom lipid composition[/bold]\n")

        # Symmetry
        sym_choice = prompt_with_context(
            self.processor,
            "Membrane symmetry: [1] Symmetric  [2] Asymmetric",
            choices=["1", "2"],
            default="1",
            module=MODULE_NAME,
            description="Membrane symmetry",
            options_map={"1": "Symmetric", "2": "Asymmetric"},
        )
        symmetric = sym_choice == "1"

        if symmetric:
            lipids, ratio = self._build_leaflet("membrane")
            if lipids:
                self.config.lipids = lipids
                self.config.ratio = ratio
                self.config.symmetric = True
        else:
            self.console.print("\n[#0f7f99]Upper leaflet:[/#0f7f99]")
            upper_lipids, upper_ratio = self._build_leaflet("upper leaflet")
            self.console.print("\n[#0f7f99]Lower leaflet:[/#0f7f99]")
            lower_lipids, lower_ratio = self._build_leaflet("lower leaflet")

            if upper_lipids and lower_lipids:
                self.config.lipids = f"{upper_lipids}//{lower_lipids}"
                self.config.ratio = f"{upper_ratio}//{lower_ratio}" if upper_ratio or lower_ratio else ""
                self.config.symmetric = False

    def _build_leaflet(self, label: str):
        """Build lipid list for one leaflet. Returns (lipids_str, ratio_str)."""
        lipids = []
        ratios = []

        self.console.print(f"  Add lipids to the {label} (enter 'd' when done, 's' to search):\n")

        while True:
            if lipids:
                table = Table(show_header=True, header_style="bold")
                table.add_column("#", width=4)
                table.add_column("Lipid", min_width=8)
                table.add_column("Ratio", width=6)
                table.add_column("Charge", width=7)
                for i, (lip, rat) in enumerate(zip(lipids, ratios), 1):
                    entry = self.lipid_library.get(lip) if self.lipid_library.is_loaded else None
                    charge = entry.charge_str if entry else "?"
                    table.add_row(str(i), lip, str(rat), charge)
                self.console.print(table)

            name = prompt_with_context(
                self.processor,
                f"  Lipid name (or 's' to search, 'd' when done)",
                module=MODULE_NAME,
                description=f"Add lipid to {label}",
            ).strip().upper()

            if name == "D":
                break
            if name == "S":
                self._lipid_browser()
                continue

            # Validate
            if self.lipid_library.is_loaded or self.lipid_library.load():
                entry = self.lipid_library.get(name)
                if entry is None:
                    self.console.print(f"[red]Unknown lipid: {name}. Use 's' to search.[/red]")
                    continue
                self.console.print(f"  [grey50]{entry.full_name} (charge: {entry.charge_str})[/grey50]")

            ratio = prompt_int_with_retry(
                self.processor,
                f"  Molar ratio for {name}",
                default=1,
                min_value=1,
                module=MODULE_NAME,
                description=f"Ratio for {name}",
            )

            lipids.append(name)
            ratios.append(ratio)

        if not lipids:
            return "", ""

        return ":".join(lipids), ":".join(str(r) for r in ratios)

    def _lipid_browser(self):
        """Browse and search the lipid library."""
        if not self.lipid_library.load():
            self.console.print(
                f"[red]{self.lipid_library.load_error}[/red]"
            )
            return

        self.console.print(
            Panel(
                f"packmol-memgen includes {len(self.lipid_library.get_all())} lipid types.\n"
                "Search by name or keyword, or browse by category.",
                title="Lipid Library",
                border_style="bright_blue",
                expand=False,
            )
        )

        while True:
            self.console.print("\n    1   Search by name/keyword")
            self.console.print("    2   Browse by category")
            self.console.print("    b   Back")

            choice = prompt_with_context(
                self.processor,
                "\nSelect",
                choices=["1", "2", "b"],
                default="1",
                module=MODULE_NAME,
                description="Lipid library action",
            )

            if choice == "b":
                return
            elif choice == "1":
                self._lipid_search()
            elif choice == "2":
                self._lipid_browse_categories()

    def _lipid_search(self):
        """Search lipids by keyword."""
        query = prompt_with_context(
            self.processor,
            "Search lipids",
            module=MODULE_NAME,
            description="Lipid search query",
        )

        results = self.lipid_library.search(query)

        if not results:
            self.console.print(f"[dark_orange3]No lipids matching '{query}'[/dark_orange3]")
            return

        self._display_lipid_table(results, f"Search results for '{query}'")

    def _lipid_browse_categories(self):
        """Browse lipids by category."""
        categories = self.lipid_library.get_categories()

        self.console.print("\n[bold]Lipid categories:[/bold]\n")
        for i, (name, desc, count) in enumerate(categories, 1):
            self.console.print(f"    {i:>2}   {name:<40} {desc} ({count})")

        choices = [str(i) for i in range(1, len(categories) + 1)] + ["b"]
        choice = prompt_with_context(
            self.processor,
            "\nSelect category (or 'b')",
            choices=choices,
            default="1",
            module=MODULE_NAME,
            description="Lipid category selection",
        )

        if choice != "b":
            cat_name = categories[int(choice) - 1][0]
            lipids = self.lipid_library.get_by_category(cat_name)
            self._display_lipid_table(lipids, cat_name)

    def _display_lipid_table(self, lipids, title: str):
        """Display a table of lipids."""
        table = Table(title=title, show_header=True, header_style="bold")
        table.add_column("Name", style="cyan", width=8)
        table.add_column("Full name", min_width=30)
        table.add_column("APL (Å²)", justify="right", width=10)
        table.add_column("Charge", justify="right", width=7)

        for entry in lipids:
            table.add_row(
                entry.name,
                entry.full_name,
                entry.apl_display,
                entry.charge_str,
            )

        self.console.print(table)

    def _lipid_raw_string(self):
        """Enter a raw packmol-memgen lipid string."""
        self.console.print(
            "\n[bold]Raw lipid string[/bold]\n"
            "[grey50]Syntax: LIPID1:LIPID2 for mixed leaflet\n"
            "        UPPER//LOWER for asymmetric\n"
            "        UPPER///LOWER for double bilayer[/grey50]\n"
        )

        lipid_str = prompt_with_context(
            self.processor,
            "Lipid string",
            module=MODULE_NAME,
            description="Raw lipid composition string",
        )

        # Validate if library is available
        if self.lipid_library.is_loaded or self.lipid_library.load():
            valid, error = self.lipid_library.validate_lipid_string(lipid_str)
            if not valid:
                self.console.print(f"[red]{error}[/red]")
                return

        self.config.lipids = lipid_str
        self.config.symmetric = "//" not in lipid_str

        # Ratio
        ratio_str = prompt_with_context(
            self.processor,
            "Ratio string (or Enter to skip)",
            default="",
            module=MODULE_NAME,
            description="Lipid ratio string",
        )
        self.config.ratio = ratio_str

        self.console.print(f"[green]Set: {lipid_str} (ratio: {ratio_str or 'equal'})[/green]")

    # ── Option 3: Solvent & Ions ─────────────────────────────────────────

    def _menu_solvent_ions(self):
        """Configure solvent and ions."""
        while True:
            self.console.print(
                Panel(
                    "Configure the aqueous phase surrounding the membrane. The water model\n"
                    "is auto-selected to match the protein force field, but can be overridden.\n"
                    "Ions neutralize the system charge and optionally add physiological salt.",
                    title="Solvent & Ions",
                    border_style="bright_blue",
                    expand=False,
                )
            )

            auto_label = f" (auto: matches {self.config.ffprot})" if self.config.water_model == "auto" else ""
            self.console.print(f"  Water model:       {self.config.effective_water_model.upper()}{auto_label}")
            self.console.print(f"  Salt:              {'Yes' if self.config.salt else 'No'}")
            if self.config.salt:
                self.console.print(f"  Salt concentration:{self.config.salt_concentration} M")
                self.console.print(f"  Cation/Anion:      {self.config.cation} / {self.config.anion}")
            self.console.print(f"  Counterions:       {'Disabled' if self.config.no_counterions else 'Automatic'}")
            if self.config.total_charge_delta != 0:
                self.console.print(f"  Charge correction: {self.config.total_charge_delta:+d}")
            self.console.print()

            self.console.print("    1   Change water model")
            self.console.print("    2   Toggle salt addition")
            self.console.print("    3   Change salt concentration")
            self.console.print("    4   Change ion types")
            self.console.print("    5   Toggle counterion addition")
            self.console.print("    6   Manual charge adjustment")
            self.console.print("    b   Back")

            choice = prompt_with_context(
                self.processor,
                "\nSelect",
                choices=["1", "2", "3", "4", "5", "6", "b"],
                default="b",
                module=MODULE_NAME,
                description="Solvent & ions option",
            )

            if choice == "b":
                return
            elif choice == "1":
                models = ["auto", "tip3p", "opc", "opc3", "spce", "spceb", "tip4pew", "tip4pd", "fb3"]
                self.console.print("\n[bold]Water models:[/bold]")
                for i, m in enumerate(models, 1):
                    label = m.upper()
                    if m == "auto":
                        label = f"Auto (currently {self.config.effective_water_model.upper()} from {self.config.ffprot})"
                    self.console.print(f"    {i}   {label}")

                wm_choice = prompt_with_context(
                    self.processor,
                    "Select water model",
                    choices=[str(i) for i in range(1, len(models) + 1)],
                    default="1",
                    module=MODULE_NAME,
                    description="Water model selection",
                )
                self.config.water_model = models[int(wm_choice) - 1]

            elif choice == "2":
                self.config.salt = not self.config.salt
                self.console.print(f"[green]Salt: {'enabled' if self.config.salt else 'disabled'}[/green]")

            elif choice == "3":
                self.config.salt_concentration = prompt_float_with_retry(
                    self.processor,
                    "Salt concentration (M)",
                    default=self.config.salt_concentration,
                    min_value=0.0,
                    max_value=5.0,
                    module=MODULE_NAME,
                    description="Salt concentration",
                )

            elif choice == "4":
                self.config.cation = prompt_with_context(
                    self.processor,
                    "Cation",
                    default=self.config.cation,
                    module=MODULE_NAME,
                    description="Cation type",
                )
                self.config.anion = prompt_with_context(
                    self.processor,
                    "Anion",
                    default=self.config.anion,
                    module=MODULE_NAME,
                    description="Anion type",
                )

            elif choice == "5":
                self.config.no_counterions = not self.config.no_counterions
                state = "disabled" if self.config.no_counterions else "enabled"
                self.console.print(f"[green]Counterion addition: {state}[/green]")

            elif choice == "6":
                self.console.print(
                    "\n[grey50]Manual charge correction for residues not recognized by\n"
                    "packmol-memgen's charge dictionary. This is added ON TOP of\n"
                    f"the auto-detected redox correction ({self.config.auto_charge_delta:+d}).[/grey50]\n"
                )
                self.config.charge_delta = prompt_int_with_retry(
                    self.processor,
                    "Manual charge delta",
                    default=self.config.charge_delta,
                    module=MODULE_NAME,
                    description="Manual charge correction",
                )

    # ── Option 4: Box & Membrane Dimensions ──────────────────────────────

    def _menu_box_dimensions(self):
        """Configure box and membrane dimensions."""
        while True:
            self.console.print(
                Panel(
                    "Controls the size of the simulation box. By default, dimensions are\n"
                    "calculated automatically from the protein size and lipid packing\n"
                    "parameters (area per lipid). You can override with fixed dimensions.\n\n"
                    "The water layer thickness controls how much bulk water sits above\n"
                    "and below the membrane — 17.5 Å is typical.",
                    title="Box & Membrane Dimensions",
                    border_style="bright_blue",
                    expand=False,
                )
            )

            self.console.print(f"  Box padding (XY):  {self.config.box_padding} Å")
            self.console.print(f"  Water layer (Z):   {self.config.water_layer} Å")
            self.console.print(f"  Leaflet width:     {self.config.leaflet_width} Å")
            sizing = "Fixed" if self.config.fixed_dims else ("Fixed XY" if self.config.fixed_xy else "Automatic")
            self.console.print(f"  Box sizing:        {sizing}")
            self.console.print(f"  Lipid counting:    {'Volume-based' if self.config.lipid_count_method == 'volume' else 'APL-based'}")
            self.console.print()

            self.console.print("    1   Box padding distance")
            self.console.print("    2   Water layer thickness")
            self.console.print("    3   Leaflet width")
            self.console.print("    4   Fix XY dimensions")
            self.console.print("    5   Fix box dimensions (X Y Z)")
            self.console.print("    6   Lipid counting method")
            self.console.print("    7   Clear fixed dimensions (return to auto)")
            self.console.print("    b   Back")

            choice = prompt_with_context(
                self.processor,
                "\nSelect",
                choices=["1", "2", "3", "4", "5", "6", "7", "b"],
                default="b",
                module=MODULE_NAME,
                description="Box dimensions option",
            )

            if choice == "b":
                return
            elif choice == "1":
                self.config.box_padding = prompt_float_with_retry(
                    self.processor, "Box padding (Å)", default=self.config.box_padding,
                    min_value=0.0, module=MODULE_NAME, description="Box padding",
                )
            elif choice == "2":
                self.config.water_layer = prompt_float_with_retry(
                    self.processor, "Water layer thickness (Å)", default=self.config.water_layer,
                    min_value=0.0, module=MODULE_NAME, description="Water layer thickness",
                )
            elif choice == "3":
                self.config.leaflet_width = prompt_float_with_retry(
                    self.processor, "Leaflet width (Å)", default=self.config.leaflet_width,
                    min_value=5.0, module=MODULE_NAME, description="Leaflet width",
                )
            elif choice == "4":
                self.config.fixed_xy = prompt_float_with_retry(
                    self.processor, "Fixed XY dimension (Å)", default=self.config.fixed_xy or 80.0,
                    min_value=10.0, module=MODULE_NAME, description="Fixed XY dimension",
                )
                self.config.fixed_dims = None
            elif choice == "5":
                x = prompt_float_with_retry(self.processor, "X (Å)", default=80.0, min_value=10.0,
                                             module=MODULE_NAME, description="Box X dimension")
                y = prompt_float_with_retry(self.processor, "Y (Å)", default=80.0, min_value=10.0,
                                             module=MODULE_NAME, description="Box Y dimension")
                z = prompt_float_with_retry(self.processor, "Z (Å)", default=100.0, min_value=10.0,
                                             module=MODULE_NAME, description="Box Z dimension")
                self.config.fixed_dims = [x, y, z]
                self.config.fixed_xy = None
            elif choice == "6":
                method = prompt_with_context(
                    self.processor, "Lipid counting method: [1] APL-based  [2] Volume-based",
                    choices=["1", "2"], default="1",
                    module=MODULE_NAME, description="Lipid counting method",
                )
                self.config.lipid_count_method = "apl" if method == "1" else "volume"
            elif choice == "7":
                self.config.fixed_xy = None
                self.config.fixed_dims = None
                self.console.print("[green]Reverted to automatic box sizing.[/green]")

    # ── Option 5: Protein Orientation ────────────────────────────────────

    def _menu_protein_orientation(self):
        """Configure protein orientation method."""
        while True:
            self.console.print(
                Panel(
                    "Membrane proteins must be oriented with their transmembrane region\n"
                    "aligned to the bilayer plane (Z=0). packmol-memgen can do this\n"
                    "automatically using MEMEMBED or PPM3, or you can provide a pre-oriented\n"
                    "structure (e.g. from the OPM database).",
                    title="Protein Orientation",
                    border_style="bright_blue",
                    expand=False,
                )
            )

            self.console.print("    1   Automatic (MEMEMBED)")
            self.console.print("    2   Automatic (PPM3)")
            self.console.print("    3   Pre-oriented (skip alignment)")
            self.console.print("    4   MEMEMBED options")
            self.console.print("    5   N-terminus orientation")
            self.console.print("    b   Back")

            choice = prompt_with_context(
                self.processor,
                "\nSelect",
                choices=["1", "2", "3", "4", "5", "b"],
                default="1",
                module=MODULE_NAME,
                description="Orientation method",
            )

            if choice == "b":
                return
            elif choice == "1":
                self.config.preoriented = False
                self.config.orientation_method = "memembed"
                self.console.print("[green]Set: Automatic orientation via MEMEMBED[/green]")
            elif choice == "2":
                self.config.preoriented = False
                self.config.orientation_method = "ppm3"
                self.console.print("[green]Set: Automatic orientation via PPM3[/green]")
            elif choice == "3":
                self.config.preoriented = True
                self.console.print("[green]Set: Pre-oriented (alignment skipped)[/green]")
            elif choice == "4":
                self._memembed_options()
            elif choice == "5":
                self._nter_orientation()

    def _memembed_options(self):
        """Configure MEMEMBED-specific options."""
        self.console.print("\n[bold]MEMEMBED options:[/bold]")
        self.console.print("    1   Optimization algorithm        Currently: " +
                          ["GA", "Grid", "Direct", "GA×5"][self.config.memembed_algorithm])
        self.console.print(f"    2   Beta barrel mode              Currently: {'Yes' if self.config.barrel_mode else 'No'}")
        self.console.print(f"    3   Keep ligands after alignment  Currently: {'Yes' if self.config.keep_ligands else 'No'}")
        self.console.print("    b   Back")

        choice = prompt_with_context(
            self.processor, "\nSelect", choices=["1", "2", "3", "b"],
            default="b", module=MODULE_NAME, description="MEMEMBED option",
        )

        if choice == "1":
            self.console.print("\n  [1] GA  [2] Grid  [3] Direct  [4] GA×5")
            alg = prompt_with_context(
                self.processor, "Algorithm", choices=["1", "2", "3", "4"],
                default="1", module=MODULE_NAME, description="MEMEMBED algorithm",
            )
            self.config.memembed_algorithm = int(alg) - 1
        elif choice == "2":
            self.config.barrel_mode = not self.config.barrel_mode
        elif choice == "3":
            self.config.keep_ligands = not self.config.keep_ligands

    def _nter_orientation(self):
        """Configure N-terminus orientation per protein chain."""
        self.console.print(
            "\n[grey50]Specify whether each chain's N-terminus faces in (cytoplasmic)\n"
            "or out (extracellular). Leave empty for automatic detection.[/grey50]\n"
        )

        nter_str = prompt_with_context(
            self.processor,
            "N-terminus orientations (comma-separated: in,out,...) or Enter for auto",
            default="",
            module=MODULE_NAME,
            description="N-terminus orientation",
        )

        if nter_str.strip():
            self.config.n_ter_orientation = [x.strip() for x in nter_str.split(",")]
        else:
            self.config.n_ter_orientation = None

    # ── Option 6: Specialized Geometry ───────────────────────────────────

    def _menu_specialized_geometry(self):
        """Configure specialized geometry options."""
        while True:
            self.console.print(
                Panel(
                    "These options control membrane curvature, packing, and less common\n"
                    "build types. The defaults produce flat bilayers suitable for most\n"
                    "simulations, so these settings are only needed for specific systems.",
                    title="Specialized Geometry",
                    border_style="bright_blue",
                    expand=False,
                )
            )

            self.console.print("    1   Membrane curvature")
            self.console.print("    2   Gaussian deformation")
            self.console.print("    3   Channel plug")
            self.console.print("    4   Head/tail plane boundaries")
            self.console.print("    5   APL / lipid offset multipliers")
            self.console.print("    6   Self-assembly mode")
            self.console.print("    7   Double bilayer (CompEL)")
            self.console.print("    8   Periodic boundary conditions")
            self.console.print("    b   Back")

            choice = prompt_with_context(
                self.processor, "\nSelect",
                choices=["1", "2", "3", "4", "5", "6", "7", "8", "b"],
                default="b", module=MODULE_NAME, description="Specialized geometry option",
            )

            if choice == "b":
                return
            elif choice == "1":
                self.config.curvature = prompt_float_with_retry(
                    self.processor, "Curvature value (0 for flat, or Enter for none)",
                    default=0.0, module=MODULE_NAME, description="Membrane curvature",
                )
                if self.config.curvature == 0.0:
                    self.config.curvature = None
            elif choice == "2":
                self.console.print("[grey50]Gaussian parameters: C (amplitude), D (width), H (height offset)[/grey50]")
                c = prompt_float_with_retry(self.processor, "C", default=10.0, module=MODULE_NAME, description="Gaussian C")
                d = prompt_float_with_retry(self.processor, "D", default=20.0, module=MODULE_NAME, description="Gaussian D")
                h = prompt_float_with_retry(self.processor, "H", default=0.0, module=MODULE_NAME, description="Gaussian H")
                self.config.gaussian_params = [c, d, h]
            elif choice == "3":
                self.config.channel_plug = prompt_float_with_retry(
                    self.processor, "Channel plug radius (Å, 0 to disable)",
                    default=0.0, module=MODULE_NAME, description="Channel plug radius",
                )
                if self.config.channel_plug == 0.0:
                    self.config.channel_plug = None
            elif choice == "4":
                self.config.head_plane = prompt_float_with_retry(
                    self.processor, "Head plane boundary (Å)", default=18.0,
                    module=MODULE_NAME, description="Head plane boundary",
                )
                self.config.tail_plane = prompt_float_with_retry(
                    self.processor, "Tail plane boundary (Å)", default=4.0,
                    module=MODULE_NAME, description="Tail plane boundary",
                )
            elif choice == "5":
                self.config.apl_offset = prompt_float_with_retry(
                    self.processor, "APL offset multiplier (1.0 = no change)",
                    default=1.0, module=MODULE_NAME, description="APL offset",
                )
                if self.config.apl_offset == 1.0:
                    self.config.apl_offset = None
                self.config.lip_offset = prompt_float_with_retry(
                    self.processor, "Lipid segment size multiplier (1.0 = default)",
                    default=self.config.lip_offset, module=MODULE_NAME, description="Lipid offset",
                )
            elif choice == "6":
                self.config.self_assembly = not self.config.self_assembly
                state = "enabled" if self.config.self_assembly else "disabled"
                self.console.print(f"[green]Self-assembly mode: {state}[/green]")
            elif choice == "7":
                self._compel_options()
            elif choice == "8":
                self.config.pbc = not self.config.pbc
                state = "enabled" if self.config.pbc else "disabled"
                self.console.print(f"[green]PBC: {state}[/green]")

    def _compel_options(self):
        """Configure computational electrophysiology options."""
        self.console.print(
            Panel(
                "Creates a system with two stacked bilayers, used for computational\n"
                "electrophysiology (CompEL) simulations. The two compartments can have\n"
                "different ion concentrations to create a transmembrane potential.\n\n"
                "Reference: Kutzner et al., Biophys J. 2011, 101(4):809-817",
                title="Double Bilayer (CompEL)",
                border_style="bright_blue",
                expand=False,
            )
        )

        self.config.double_bilayer = confirm_with_context(
            self.processor, "Enable double bilayer?",
            default=self.config.double_bilayer,
            module=MODULE_NAME, description="Enable double bilayer",
        )

        if self.config.double_bilayer:
            self.config.charge_imbalance = prompt_int_with_retry(
                self.processor, "Charge imbalance (electron charge units)",
                default=self.config.charge_imbalance,
                module=MODULE_NAME, description="Charge imbalance",
            )
            ion_choice = prompt_with_context(
                self.processor, "Imbalance ion: [1] Cation  [2] Anion",
                choices=["1", "2"], default="1",
                module=MODULE_NAME, description="Imbalance ion type",
            )
            self.config.imbalance_ion = "cat" if ion_choice == "1" else "an"

    # ── Option 7: PACKMOL Settings ───────────────────────────────────────

    def _menu_packmol_settings(self):
        """Configure PACKMOL packing parameters."""
        while True:
            self.console.print(
                Panel(
                    "Controls how PACKMOL places molecules. The defaults work well for\n"
                    "most systems. Increase nloop values if packing fails to converge.",
                    title="PACKMOL Settings",
                    border_style="bright_blue",
                    expand=False,
                )
            )

            self.console.print(f"  Clash tolerance:      {self.config.tolerance} Å")
            self.console.print(f"  Protein radius:       {self.config.protein_radius} Å")
            self.console.print(f"  Loops (individual):   {self.config.nloop}")
            self.console.print(f"  Loops (all-together): {self.config.nloop_all}")
            self.console.print(f"  GENCAN iterations:    {self.config.gencan_iterations}")
            self.console.print(f"  Move fraction:        {self.config.move_fraction}")
            self.console.print()

            self.console.print("    1   Clash tolerance")
            self.console.print("    2   Protein atom radius")
            self.console.print("    3   Packing iterations (nloop)")
            self.console.print("    4   All-together iterations")
            self.console.print("    5   GENCAN max iterations")
            self.console.print("    6   Move fraction")
            self.console.print("    7   Troubleshooting options")
            self.console.print("    8   Output options (trajectories, plots)")
            self.console.print("    b   Back")

            choice = prompt_with_context(
                self.processor, "\nSelect",
                choices=["1", "2", "3", "4", "5", "6", "7", "8", "b"],
                default="b", module=MODULE_NAME, description="PACKMOL option",
            )

            if choice == "b":
                return
            elif choice == "1":
                self.config.tolerance = prompt_float_with_retry(
                    self.processor, "Clash tolerance (Å)", default=self.config.tolerance,
                    min_value=0.1, module=MODULE_NAME, description="Clash tolerance",
                )
            elif choice == "2":
                self.config.protein_radius = prompt_float_with_retry(
                    self.processor, "Protein atom radius (Å)", default=self.config.protein_radius,
                    min_value=0.1, module=MODULE_NAME, description="Protein radius",
                )
            elif choice == "3":
                self.config.nloop = prompt_int_with_retry(
                    self.processor, "Individual packing loops", default=self.config.nloop,
                    min_value=1, module=MODULE_NAME, description="nloop",
                )
            elif choice == "4":
                self.config.nloop_all = prompt_int_with_retry(
                    self.processor, "All-together packing loops", default=self.config.nloop_all,
                    min_value=1, module=MODULE_NAME, description="nloop_all",
                )
            elif choice == "5":
                self.config.gencan_iterations = prompt_int_with_retry(
                    self.processor, "GENCAN iterations per loop", default=self.config.gencan_iterations,
                    min_value=1, module=MODULE_NAME, description="GENCAN iterations",
                )
            elif choice == "6":
                self.config.move_fraction = prompt_float_with_retry(
                    self.processor, "Move fraction", default=self.config.move_fraction,
                    min_value=0.001, max_value=1.0, module=MODULE_NAME, description="Move fraction",
                )
            elif choice == "7":
                self.config.move_bad_random = confirm_with_context(
                    self.processor, "Randomize badly-placed molecules?",
                    default=self.config.move_bad_random,
                    module=MODULE_NAME, description="movebadrandom",
                )
                self.config.short_penalty = confirm_with_context(
                    self.processor, "Add short-range overlap penalty?",
                    default=self.config.short_penalty,
                    module=MODULE_NAME, description="short_penalty",
                )
                self.config.pack_all = confirm_with_context(
                    self.processor, "Skip individual packing steps (packall)?",
                    default=self.config.pack_all,
                    module=MODULE_NAME, description="packall",
                )
            elif choice == "8":
                self.config.save_trajectory = confirm_with_context(
                    self.processor, "Save intermediate PDB snapshots?",
                    default=self.config.save_trajectory,
                    module=MODULE_NAME, description="Save trajectory",
                )
                self.config.plot_optimization = confirm_with_context(
                    self.processor, "Create optimization function plot?",
                    default=self.config.plot_optimization,
                    module=MODULE_NAME, description="Plot optimization",
                )
                self.config.random_seed = confirm_with_context(
                    self.processor, "Use random seed?",
                    default=self.config.random_seed,
                    module=MODULE_NAME, description="Random seed",
                )

    # ── Option 8: Force Fields ───────────────────────────────────────────

    def _menu_force_fields(self, workspace=None):
        """
        Configure force field selections using tLEaP's comprehensive options.

        Presents the same forcefield categories as the Topology Generator
        (protein, modified AA, lipids, water, ions, etc.) so the user makes all
        FF decisions here. The tLEaP module then reads these settings directly
        from workspace without re-prompting.
        """
        from rich.table import Table
        from proprep.tleap_prep.tleap_input_generator import TLeapInputGenerator
        from proprep.forcefield_params.forcefield_catalog import (
            recommended_water_for_protein,
            recommended_ions_for_water,
        )

        ff_options = TLeapInputGenerator.FORCEFIELD_OPTIONS

        while True:
            self.console.print(
                Panel(
                    "Select ALL force fields for this system. These will be passed directly\n"
                    "to the Topology Generator — you won't be asked again downstream.\n\n"
                    "Only protein, lipid, and water FFs affect packmol-memgen's packing.\n"
                    "Additional FFs (constant pH/Eh, modified AAs, etc.) are for tLEaP only.",
                    title="Force Fields (Comprehensive)",
                    border_style="bright_blue",
                    expand=False,
                )
            )

            # Show current selections
            if hasattr(self, '_selected_forcefields') and self._selected_forcefields:
                self.console.print("[bold]Current selections:[/bold]")
                for cat, sel in self._selected_forcefields.items():
                    if sel is None:
                        self.console.print(f"  {cat}: [grey50]None[/grey50]")
                    elif isinstance(sel, list):
                        names = [s['name'] for s in sel] if sel else ["None"]
                        self.console.print(f"  {cat}: {', '.join(names)}")
                    else:
                        self.console.print(f"  {cat}: {sel['name']}")
                self.console.print()

            self.console.print("    1   Select all force fields (guided walkthrough)")
            self.console.print("    b   Back")

            choice = prompt_with_context(
                self.processor, "\nSelect",
                choices=["1", "b"],
                default="1", module=MODULE_NAME, description="Force field option",
            )

            if choice == "b":
                return

            # Guided walkthrough through all categories
            selected = {}
            category_order = ['protein', 'modified_aa', 'dna', 'rna', 'carbohydrates',
                             'lipids', 'small_molecules', 'water', 'ions']

            # Redox-site force-field prerequisites. When a detected and
            # transformed cofactor requires particular leaprcs (a bis-his b-type
            # heme, for instance, needs both constph and conste), show the
            # prerequisites panel and tag the standard-FF options that satisfy
            # them. The Topology Generator owns this logic and reads the prereqs
            # from transformer_info in the shared workspace, so reusing it here
            # makes the annotation identical to what the user sees downstream.
            topgen = self.processor.get_module_instance("Topology Generator")
            cofactor_prereq_groups = []
            if topgen is not None:
                try:
                    cofactor_prereq_groups = topgen._collect_cofactor_prereq_groups() or []
                except Exception:
                    cofactor_prereq_groups = []
                if cofactor_prereq_groups:
                    topgen._show_cofactor_ff_prerequisites_panel()

            for category in category_order:
                if category not in ff_options:
                    continue

                cat_info = ff_options[category]
                options = cat_info['options']

                # Tag options that FULLY satisfy a cofactor's leaprc prereqs
                # (their leaprc set hits at least one entry in every AND-group).
                # For the satisfied category the tagged option becomes the sole
                # recommendation, so the default lands on the FF the cofactor
                # needs rather than the generic catalog default.
                if cofactor_prereq_groups:
                    def _opt_full_triggers(opt):
                        opt_leaprcs = topgen._option_leaprc_set(opt)
                        if not opt_leaprcs:
                            return set()
                        return {
                            cof["residue_name"]
                            for cof in cofactor_prereq_groups
                            if topgen._option_fully_satisfies(opt_leaprcs, cof["groups"])
                        }
                    if any(_opt_full_triggers(opt) for opt in options):
                        options = [dict(opt) for opt in options]
                        for opt in options:
                            triggers = _opt_full_triggers(opt)
                            if triggers:
                                res_join = " + ".join(sorted(triggers))
                                opt['recommended'] = True
                                opt['recommendation_reason'] = f"satisfies your {res_join} selection"
                                opt['_cofactor_trigger_count'] = len(triggers)
                            else:
                                opt.pop('recommended', None)

                # Recommend the water model that matches the protein FF just
                # picked. Water is processed after protein, so selected['protein']
                # is populated here. Without this, the catalog's static default
                # (OPC, for ff19SB) stays marked even when the user chose, e.g.,
                # the constant-pH/redox FF, whose matching water is TIP3P.
                if category == 'water':
                    protein_sel = selected.get('protein')
                    rec_water = (
                        recommended_water_for_protein(protein_sel.get('leaprc'))
                        if protein_sel else None
                    )
                    if rec_water:
                        options = [dict(opt) for opt in options]
                        for opt in options:
                            if opt['name'] == rec_water:
                                opt['recommended'] = True
                                opt['recommendation_reason'] = f"matches {protein_sel['name']}"
                            else:
                                opt.pop('recommended', None)

                # Recommend the divalent+ ion set that matches the water model
                # just picked. Ions are processed after water, so selected['water']
                # is populated here. Without this, the catalog's static default
                # (12-6-4 OPC) stays marked even against, e.g., TIP3P water.
                if category == 'ions':
                    water_sel = selected.get('water')
                    water_name = water_sel.get('name') if water_sel else None
                    rec_ion = recommended_ions_for_water(water_name)
                    if rec_ion:
                        options = [dict(opt) for opt in options]
                        for opt in options:
                            if opt['name'] == rec_ion:
                                opt['recommended'] = True
                                if rec_ion == 'Default only' and water_name:
                                    opt['recommendation_reason'] = f"no Li/Merz 12-6-4 set for {water_name}"
                                elif water_name:
                                    opt['recommendation_reason'] = f"matches {water_name}"
                                else:
                                    opt['recommendation_reason'] = "no explicit water model"
                            else:
                                opt.pop('recommended', None)

                self.console.print(f"\n[bold #0f7f99]━━━ {cat_info['title']} ━━━[/bold #0f7f99]")
                if 'description' in cat_info:
                    self.console.print(f"[grey50]{cat_info['description']}[/grey50]")

                # "None" is always row 1 when allowed, matching the Topology
                # Generator so the user doesn't have to count to a different
                # final index per category.
                allow_none = cat_info.get('allow_none', False)
                none_offset = 1 if allow_none else 0

                table = Table(show_header=False, box=None, padding=(0, 2))
                table.add_column("Choice", style="bold")
                table.add_column("Name", style="green")
                table.add_column("Description")

                if allow_none:
                    table.add_row("1", "None", cat_info.get('none_text', 'Skip this category'))

                for i, opt in enumerate(options, start=1 + none_offset):
                    name = opt['name']
                    if opt.get('recommended'):
                        reason = opt.get('recommendation_reason')
                        if reason:
                            name += f" [dark_orange3](Recommended — {reason})[/dark_orange3]"
                        else:
                            name += " [dark_orange3](Recommended)[/dark_orange3]"
                    table.add_row(str(i), name, opt['description'])

                self.console.print(table)

                if cat_info.get('multi_select', False):
                    if allow_none:
                        self.console.print("[grey50]Enter 1 for none, or comma-separated choices (e.g., 2,3), or press Enter for none[/grey50]")
                    else:
                        self.console.print("[grey50]Enter choices separated by commas (e.g., 1,2)[/grey50]")
                    choice_str = prompt_with_context(
                        self.processor, f"Select {cat_info['title'].lower()}", default="",
                        module=MODULE_NAME,
                        description=f"Select {cat_info['title'].lower()} forcefield(s)",
                    )
                    if not choice_str.strip():
                        selected[category] = []
                    else:
                        choices = [c.strip() for c in choice_str.split(',')]
                        selected_items = []
                        for c in choices:
                            try:
                                idx = int(c) - 1 - none_offset
                                if 0 <= idx < len(options):
                                    selected_items.append(options[idx])
                            except (ValueError, IndexError):
                                pass
                        selected[category] = selected_items
                else:
                    # Default to the recommended option (shifted past the None
                    # row if present), or the first real option if nothing is
                    # marked. When cofactor prereqs are in play, prefer the
                    # option that satisfies the most cofactors, matching the
                    # Topology Generator's default-selection rule.
                    default = str(1 + none_offset)
                    best_count = -1
                    for i, opt in enumerate(options, start=1 + none_offset):
                        if not opt.get('recommended'):
                            continue
                        count = opt.get('_cofactor_trigger_count', 0)
                        if count > best_count:
                            best_count = count
                            default = str(i)

                    choice_str = prompt_with_context(
                        self.processor, f"Select {cat_info['title'].lower()}", default=default,
                        module=MODULE_NAME,
                        description=f"Select {cat_info['title'].lower()} forcefield",
                    )

                    try:
                        choice = int(choice_str)
                        if allow_none and choice == 1:
                            selected[category] = None
                        elif 1 + none_offset <= choice <= len(options) + none_offset:
                            selected[category] = options[choice - 1 - none_offset]
                        else:
                            selected[category] = options[0]
                    except ValueError:
                        selected[category] = options[0]

            # Store the full selection
            self._selected_forcefields = selected

            # Sync MembraneConfig fields for packmol-memgen CLI args
            prot_sel = selected.get('protein')
            if prot_sel and isinstance(prot_sel, dict):
                # Extract short name (e.g., "ff14SB" from leaprc.protein.ff14SB)
                leaprc = prot_sel.get('leaprc', '')
                if isinstance(leaprc, str) and 'leaprc.protein.' in leaprc:
                    self.config.ffprot = leaprc.replace('leaprc.protein.', '')
                elif isinstance(leaprc, list):
                    # constph/conste case — packmol-memgen doesn't understand these,
                    # keep ffprot as ff14SB (the underlying protein FF)
                    pass

            lip_sel = selected.get('lipids')
            if lip_sel and isinstance(lip_sel, dict):
                leaprc = lip_sel.get('leaprc', '')
                if isinstance(leaprc, str) and 'leaprc.' in leaprc:
                    self.config.fflip = leaprc.replace('leaprc.', '')

            wat_sel = selected.get('water')
            if wat_sel and isinstance(wat_sel, dict):
                leaprc = wat_sel.get('leaprc', '')
                if isinstance(leaprc, str) and 'leaprc.water.' in leaprc:
                    self.config.water_model = leaprc.replace('leaprc.water.', '')

            # Store in workspace for tLEaP to read directly
            if workspace is not None:
                self.update_workspace(workspace, "selected_standard_forcefields", selected)

            # Summary
            self.console.print("\n[bold green]Selected Forcefields:[/bold green]")
            for category in category_order:
                sel = selected.get(category)
                if sel is None:
                    self.console.print(f"  {category}: [grey50]None[/grey50]")
                elif isinstance(sel, list):
                    names = [s['name'] for s in sel] if sel else ["None"]
                    self.console.print(f"  {category}: {', '.join(names)}")
                elif isinstance(sel, dict):
                    self.console.print(f"  {category}: {sel['name']}")

            # Custom force fields + bond directives for any detected redox /
            # specialized residues. Runs the Topology Generator's own pickers so
            # the choice (e.g. RESP vs CM5 for a bis-his b-type heme) is made
            # once here and reused downstream.
            self._select_redox_forcefields_and_bonds(workspace)
            return

    def _select_redox_forcefields_and_bonds(self, workspace):
        """Prompt for custom redox-site force-field parameter sets and gather
        tLEaP bond directives, reusing the Topology Generator.

        The membrane builder's pre-tLEaP hydrogen pass and the downstream
        Topology Generator both read the results from the shared workspace, so
        the user picks the parameter set (and the metal-ligand bonds are
        derived) exactly once. Results land in the same workspace keys the
        Topology Generator writes:
          - single_state_ff_requirements
          - single_state_selected_forcefields
          - ff_resolved_atom_types
          - combined_tleap_commands
        Gated on detected redox sites (which imply a protein); empty bilayers
        and site-free systems fall through untouched.
        """
        if workspace is None:
            return
        if not workspace.get("detected_redox_sites"):
            return

        topgen = self.processor.get_module_instance("Topology Generator")
        if topgen is None:
            self.console.print(
                "[dark_orange3]Topology Generator unavailable — skipping redox force-field "
                "selection; the hydrogen-addition pass will fall back to default sets.[/dark_orange3]"
            )
            return

        # (a) Custom force-field parameter sets for redox / specialized residues.
        try:
            requirements = topgen._get_single_state_forcefield_requirements()
        except Exception as e:
            self.console.print(f"[dark_orange3]Redox force-field requirements unavailable: {e}[/dark_orange3]")
            requirements = None

        if requirements:
            self.console.print("\n[bold #0f7f99]━━━ REDOX SITE FORCE FIELDS ━━━[/bold #0f7f99]")
            self.console.print(
                "[grey50]Custom parameter sets for detected redox/specialized residues "
                "(e.g. RESP vs CM5 charges). Your choice loads into the tLEaP hydrogen-addition "
                "pass and carries through to the Topology Generator.[/grey50]"
            )
            try:
                selected_ff = topgen._select_forcefields_for_single_state(requirements)
                selected_ff, ff_atom_types = topgen._resolve_ff_collisions(selected_ff)
                self.update_workspace(workspace, "single_state_ff_requirements", requirements)
                self.update_workspace(workspace, "single_state_selected_forcefields", selected_ff)
                self.update_workspace(workspace, "ff_resolved_atom_types", ff_atom_types)
            except Exception as e:
                self.console.print(f"[dark_orange3]Redox force-field selection failed: {e}[/dark_orange3]")

        # (b) Bond directives (metal-ligand coordinate/covalent bonds). Derived
        # from detected_redox_sites and stored in combined_tleap_commands.
        try:
            topgen.gather_bond_definitions()
        except Exception as e:
            self.console.print(f"[dark_orange3]Could not gather bond definitions: {e}[/dark_orange3]")

    # ── Option 9: Custom Parameters & Solutes ────────────────────────────

    def _menu_custom_parameters(self):
        """Configure custom parameters and solute molecules."""
        while True:
            self.console.print(
                Panel(
                    "Add solute molecules (drugs, ligands, cofactors) to the simulation box.\n"
                    "These are placed separately from the protein.\n\n"
                    "Each solute needs: a PDB file, its formal charge, and optionally\n"
                    "custom force field parameters (frcmod + lib files) for tLEaP.",
                    title="Custom Parameters & Solutes",
                    border_style="bright_blue",
                    expand=False,
                )
            )

            # Show auto-injected redox info if present
            if self.config.auto_charge_delta != 0:
                self.console.print(
                    f"  [#0f7f99]Auto-detected from workspace:[/#0f7f99]\n"
                    f"    Charge correction: {self.config.auto_charge_delta:+d} (redox sites)\n"
                )

            # Show current solutes
            if self.config.solutes:
                table = Table(show_header=True, header_style="bold")
                table.add_column("#", width=4)
                table.add_column("File")
                table.add_column("Charge", justify="right")
                table.add_column("Count")
                table.add_column("Placement")
                table.add_column("Parameters")
                for i, s in enumerate(self.config.solutes, 1):
                    placement = "membrane" if s.in_membrane else "solution"
                    params = f"{Path(s.frcmod).name} + {Path(s.lib).name}" if s.frcmod and s.lib else "none"
                    table.add_row(str(i), Path(s.pdb_file).name, str(s.charge), s.count, placement, params)
                self.console.print(table)
            else:
                self.console.print("  [grey50]No solutes added.[/grey50]")

            self.console.print()
            self.console.print("    1   Add solute")
            self.console.print("    2   Remove solute")
            self.console.print("    b   Back")

            choice = prompt_with_context(
                self.processor, "\nSelect",
                choices=["1", "2", "b"],
                default="b", module=MODULE_NAME, description="Custom parameters option",
            )

            if choice == "b":
                return
            elif choice == "1":
                self._add_solute()
            elif choice == "2":
                if self.config.solutes:
                    idx = prompt_int_with_retry(
                        self.processor, "Remove solute #",
                        default=1, min_value=1, max_value=len(self.config.solutes),
                        module=MODULE_NAME, description="Remove solute index",
                    )
                    removed = self.config.solutes.pop(idx - 1)
                    self.console.print(f"[green]Removed: {Path(removed.pdb_file).name}[/green]")
                else:
                    self.console.print("[dark_orange3]No solutes to remove.[/dark_orange3]")

    def _add_solute(self):
        """Add a solute molecule."""
        pdb_file = prompt_with_context(
            self.processor, "Solute PDB file path",
            module=MODULE_NAME, description="Solute PDB file",
        )
        if not Path(pdb_file).exists():
            self.console.print(f"[red]File not found: {pdb_file}[/red]")
            return

        charge = prompt_int_with_retry(
            self.processor, "Formal charge of this molecule",
            default=0, module=MODULE_NAME, description="Solute charge",
        )

        count_str = prompt_with_context(
            self.processor, "Count (integer) or concentration (e.g. '0.05M')",
            default="1", module=MODULE_NAME, description="Solute count/concentration",
        )

        in_membrane = confirm_with_context(
            self.processor, "Place in membrane (vs. in solution)?",
            default=False, module=MODULE_NAME, description="Solute placement",
        )

        prot_distance = None
        if not in_membrane and self.config.protein_pdb:
            set_dist = confirm_with_context(
                self.processor, "Set minimum distance from protein?",
                default=False, module=MODULE_NAME, description="Solute distance restraint",
            )
            if set_dist:
                prot_distance = prompt_float_with_retry(
                    self.processor, "Minimum cylindrical distance from protein (Å)",
                    default=10.0, min_value=0.0,
                    module=MODULE_NAME, description="Solute-protein distance",
                )

        # Force field parameters
        frcmod = None
        lib = None
        self.console.print("\n  Force field parameters for tLEaP:")
        self.console.print("    1   Provide frcmod + lib files (non-standard molecule)")
        self.console.print("    2   Standard residue, covered by the loaded force fields")
        self.console.print(
            "  [grey50]A non-standard molecule (drug, ligand, cofactor) needs option 1, "
            "or it must be parameterized before topology generation. packmol-memgen packs "
            "it either way, but tLEaP cannot build the topology without parameters.[/grey50]"
        )

        param_choice = prompt_with_context(
            self.processor, "Select",
            choices=["1", "2"], default="2",
            module=MODULE_NAME, description="Solute parameter source",
        )

        use_gaff2 = False
        if param_choice == "1":
            frcmod = prompt_with_context(
                self.processor, "frcmod file path",
                module=MODULE_NAME, description="Solute frcmod path",
            )
            lib = prompt_with_context(
                self.processor, "lib file path",
                module=MODULE_NAME, description="Solute lib path",
            )
            if not Path(frcmod).exists():
                self.console.print(f"[dark_orange3]Warning: frcmod not found: {frcmod}[/dark_orange3]")
            if not Path(lib).exists():
                self.console.print(f"[dark_orange3]Warning: lib not found: {lib}[/dark_orange3]")

        solute = SoluteConfig(
            pdb_file=str(Path(pdb_file).resolve()),
            charge=charge,
            count=count_str,
            in_membrane=in_membrane,
            prot_distance=prot_distance,
            frcmod=str(Path(frcmod).resolve()) if frcmod else None,
            lib=str(Path(lib).resolve()) if lib else None,
            use_gaff2=use_gaff2,
        )

        self.config.solutes.append(solute)
        self.console.print(f"[green]Added solute: {Path(pdb_file).name} (charge {charge:+d}, count {count_str})[/green]")

    # ── Review ───────────────────────────────────────────────────────────

    def _show_review(self, workspace):
        """Show full configuration review."""
        # Defensively re-show the configured protein. Other modules
        # may have swapped the active structure between Protein
        # Selection and now; this brings the viewer back in sync with
        # what's about to be built.
        self._show_protein_in_viewer()

        lines = []

        lines.append("[bold]PROTEIN[/bold]")
        if self.config.protein_pdb:
            lines.append(f"  Structure:          {Path(self.config.protein_pdb).name}")
        else:
            lines.append("  Structure:          [grey50]none (empty bilayer)[/grey50]")
        lines.append("  Protonation:        tLEaP (hydrogen addition from AMBER library templates)")
        orient = "Pre-oriented" if self.config.preoriented else f"Auto ({self.config.orientation_method})"
        lines.append(f"  Orientation:        {orient}")
        solute_count = len(self.config.solutes)
        lines.append(f"  Solutes:            {solute_count if solute_count else 'None'}")

        lines.append("")
        if self.config.solvate_only:
            lines.append("[bold]MEMBRANE[/bold]")
            lines.append("  Solvate only:       yes (no lipids, water + ions only)")
        else:
            lines.append("[bold]MEMBRANE[/bold]")
            sym_label = "symmetric" if self.config.symmetric else "asymmetric"
            lines.append(f"  Lipids ({sym_label}): {self.config.lipids}")
            if self.config.ratio:
                lines.append(f"  Ratio:              {self.config.ratio}")
            lines.append(f"  Leaflet width:      {self.config.leaflet_width} Å")

        lines.append("")
        lines.append("[bold]BOX[/bold]")
        if self.config.fixed_dims:
            lines.append(f"  Dimensions:         {self.config.fixed_dims[0]} × {self.config.fixed_dims[1]} × {self.config.fixed_dims[2]} Å (fixed)")
        elif self.config.fixed_xy:
            lines.append(f"  XY dimensions:      {self.config.fixed_xy} × {self.config.fixed_xy} Å (fixed)")
            lines.append(f"  Z:                  Auto (water layer {self.config.water_layer} Å)")
        else:
            lines.append(f"  XY padding:         {self.config.box_padding} Å")
            lines.append(f"  Water layer:        {self.config.water_layer} Å")
            lines.append("  Sizing:             Automatic")

        lines.append("")
        lines.append("[bold]SOLVENT & IONS[/bold]")
        lines.append(f"  Water model:        {self.config.effective_water_model.upper()}")
        if self.config.salt:
            lines.append(f"  Salt:               {self.config.salt_concentration} M {self.config.cation}/{self.config.anion}")
        else:
            lines.append("  Salt:               None")
        lines.append(f"  Neutralization:     {'Disabled' if self.config.no_counterions else 'Automatic'}")
        if self.config.total_charge_delta != 0:
            parts = []
            if self.config.auto_charge_delta != 0:
                parts.append(f"{self.config.auto_charge_delta:+d} auto (redox)")
            if self.config.charge_delta != 0:
                parts.append(f"{self.config.charge_delta:+d} manual")
            lines.append(f"  Charge correction:  {self.config.total_charge_delta:+d} ({', '.join(parts)})")

        lines.append("")
        lines.append("[bold]FORCE FIELDS[/bold]")
        if hasattr(self, '_selected_forcefields') and self._selected_forcefields:
            for cat in ['protein', 'modified_aa', 'dna', 'rna', 'carbohydrates',
                       'lipids', 'small_molecules', 'water', 'ions']:
                sel = self._selected_forcefields.get(cat)
                if sel is None:
                    continue
                elif isinstance(sel, list) and sel:
                    names = [s['name'] for s in sel]
                    lines.append(f"  {cat:20s} {', '.join(names)}")
                elif isinstance(sel, dict):
                    lines.append(f"  {cat:20s} {sel['name']}")
        else:
            lines.append(f"  Protein:            {self.config.ffprot}")
            lines.append(f"  Lipid:              {self.config.fflip}")
        lines.append(f"  leaprc needed:      {', '.join(self.config.get_leaprc_requirements())}")

        lines.append("")
        lines.append("[bold]PACKMOL[/bold]")
        lines.append(f"  Tolerance: {self.config.tolerance} Å    "
                     f"Loops: {self.config.nloop}/{self.config.nloop_all}    "
                     f"Protein radius: {self.config.protein_radius} Å")

        # Equivalent command
        cli_args = self.config.to_cli_args()
        cli_cmd = "packmol-memgen " + " ".join(cli_args)

        lines.append("")
        lines.append("[bold]Equivalent command:[/bold]")
        lines.append(f"  [grey50]{cli_cmd}[/grey50]")

        self.console.print(
            Panel(
                "\n".join(lines),
                title="Membrane Builder — Full Configuration Review",
                border_style="green",
                expand=False,
            )
        )

    # ── Generate & Run ───────────────────────────────────────────────────

    @staticmethod
    def _fix_ter_records(pdb_path: str) -> int:
        """
        Remove spurious TER records from packmol-memgen output.

        packmol-memgen inserts TER records around every non-standard residue
        (HBO, FHO, RHO, PRN, etc.), breaking tLEaP's chain parsing. This
        method removes TER records that appear between protein/cofactor atoms,
        keeping only TERs at true boundaries (e.g., before lipid MEMB atoms,
        between lipid molecules, before water/ions).

        The heuristic: a TER is spurious if the ATOM/HETATM line immediately
        before it and the one immediately after it are both non-lipid/non-solvent
        (i.e., neither has 'MEMB' segment nor is WAT/ion).
        """
        solvent_ion_residues = {"WAT", "HOH", "Na+", "Cl-", "K+", "Cs+", "Mg2", "Ca2", "Zn2"}

        with open(pdb_path, 'r') as f:
            lines = f.readlines()

        def is_membrane_or_solvent(line):
            """Check if an ATOM/HETATM line is lipid (MEMB segment) or solvent/ion."""
            # MEMB segment ID can appear at varying columns depending on PDB formatting.
            # Search after the coordinate/B-factor fields (column 55+).
            if "MEMB" in line[54:]:
                return True
            resname = line[17:20].strip() if len(line) > 20 else ""
            return resname in solvent_ion_residues

        def get_chain_resid(line):
            """Extract (chain, resid) from an ATOM/HETATM line."""
            return (line[21:22], line[22:26].strip())

        def has_oxt_in_residue(all_lines, atom_line):
            """Check if any atom in the same residue as atom_line has atom name OXT."""
            chain_resid = get_chain_resid(atom_line)
            # Scan backwards through recent lines for same residue
            for j in range(len(all_lines) - 1, max(0, len(all_lines) - 50) - 1, -1):
                ln = all_lines[j]
                if ln.startswith(("ATOM", "HETATM")):
                    if get_chain_resid(ln) == chain_resid:
                        if ln[12:16].strip() == "OXT":
                            return True
                    elif get_chain_resid(ln) != chain_resid:
                        break  # Past this residue
            return False

        output_lines = []
        removed = 0

        for i, line in enumerate(lines):
            if not line.startswith("TER"):
                output_lines.append(line)
                continue

            # Find the previous ATOM/HETATM line
            prev_atom = None
            for j in range(len(output_lines) - 1, -1, -1):
                if output_lines[j].startswith(("ATOM", "HETATM")):
                    prev_atom = output_lines[j]
                    break

            # Find the next ATOM/HETATM line
            next_atom = None
            for j in range(i + 1, len(lines)):
                if lines[j].startswith(("ATOM", "HETATM")):
                    next_atom = lines[j]
                    break

            # Keep TER if:
            #   - either neighbor is membrane/solvent (true chain boundary)
            #   - no neighbor found (end of file)
            #   - previous residue contains OXT (protein C-terminus)
            if prev_atom is None or next_atom is None:
                output_lines.append(line)
            elif is_membrane_or_solvent(prev_atom) or is_membrane_or_solvent(next_atom):
                output_lines.append(line)
            elif has_oxt_in_residue(output_lines, prev_atom):
                output_lines.append(line)
            else:
                # Both neighbors are protein/cofactor, not at chain end — spurious TER
                removed += 1

        if removed > 0:
            with open(pdb_path, 'w') as f:
                f.writelines(output_lines)

        return removed

    def _run_pre_tleap_hydrogen_pass(self, workspace, work_dir: str) -> Optional[str]:
        """
        Run tLEaP to add hydrogen atoms using AMBER library templates.

        This replaces packmol-memgen's default use of reduce for protonation.
        tLEaP correctly handles both standard residues and custom forcefield
        residues (e.g., HBO, FHO, RHO from redox site transformers), and
        respects protonation state assignments (HID/HIE/HIP) from upstream
        modules. The hydrogen-added PDB is then passed to packmol-memgen
        with --notprotonate --nottrim.

        Returns:
            Path to the hydrogen-added PDB, or None if the pass fails.
        """
        import subprocess

        # Empty-bilayer case: nothing to protonate. Lipid templates already
        # carry hydrogens, so packmol-memgen's --notprotonate path is fine.
        if not self.config.protein_pdb:
            return None

        self.console.print("\n[bold dark_orange3]tLEaP hydrogen addition[/bold dark_orange3]")
        self.console.print(
            "[grey50]Adding hydrogens via tLEaP using AMBER library templates[/grey50]"
        )

        transformer_info = workspace.get("transformer_info") or []

        lib_files = []
        frcmod_files = []
        atom_type_lines = []

        def _extend_files(ff_set):
            fm = ff_set.get("frcmod")
            if fm:
                frcmod_files.extend(fm if isinstance(fm, list) else [fm])
            lb = ff_set.get("lib")
            if lb:
                lib_files.extend(lb if isinstance(lb, list) else [lb])

        # Prefer the parameter sets the user chose during FF selection
        # (_select_redox_forcefields_and_bonds stores them via the Topology
        # Generator). This is what loads the exact RESP/CM5 set the user picked.
        selected_redox_ffs = workspace.get("single_state_selected_forcefields") or {}
        if selected_redox_ffs:
            requirements = workspace.get("single_state_ff_requirements") or {}
            resolved_atom_types = workspace.get("ff_resolved_atom_types") or []
            for ff_set in selected_redox_ffs.values():
                if isinstance(ff_set, dict):
                    _extend_files(ff_set)
            # Atom types: each selected set's metadata types + resolver-renamed
            # types, matching the Topology Generator's atom-types section.
            for key, info in requirements.items():
                if key in selected_redox_ffs:
                    for entry in info.get("atom_types", []):
                        if entry not in atom_type_lines:
                            atom_type_lines.append(entry)
            for entry in resolved_atom_types:
                if entry not in atom_type_lines:
                    atom_type_lines.append(entry)
        else:
            # Fallback (no recorded choice): auto-pick the default set so
            # non-interactive / site-free flows don't regress.
            from proprep.forcefield_params.loader import discover_forcefield_files
            from proprep.redoxsite_prep.transformation.redox_transformer_framework import redox_transformer_registry

            sites_with_ff = [s for s in transformer_info if s.get("has_transformer") and s.get("transformer_type")]
            seen_keys = set()
            for site in sites_with_ff:
                ttype = site["transformer_type"]
                rstate = site.get("redox_state", "")
                sstate = site.get("spin_state", "")
                key = (ttype, rstate, sstate)
                if key in seen_keys:
                    continue
                seen_keys.add(key)

                transformer_class = redox_transformer_registry.get_transformer(ttype)
                if not transformer_class or not transformer_class.FORCEFIELD_PATH:
                    continue

                ff_sets = discover_forcefield_files(transformer_class.FORCEFIELD_PATH, rstate, sstate)
                if not ff_sets:
                    continue

                # Use the default set (first one marked is_default, or just first)
                ff_set = next((fs for fs in ff_sets if fs.get("is_default")), ff_sets[0])
                _extend_files(ff_set)

                # Atom types from transformer_info
                for at in site.get("atom_types", []):
                    atom_type_lines.append(at)

        # Also collect preprocessing lib/frcmod files (small molecules, etc.)
        preproc_libs = workspace.get("preprocessing_lib_files", []) or []
        preproc_frcmods = workspace.get("preprocessing_frcmod_files", []) or []
        preproc_atom_types = workspace.get("preprocessing_atom_types", []) or []

        lib_files.extend(preproc_libs)
        frcmod_files.extend(preproc_frcmods)
        atom_type_lines.extend(preproc_atom_types)

        # Build tLEaP script using the full FF selection from the membrane builder
        protein_pdb = os.path.abspath(self.config.protein_pdb)
        output_prefix = os.path.join(work_dir, "pretleap_hydrogen")

        script_lines = ["# Pre-tLEaP hydrogen addition for membrane builder"]

        # Use full FF selection if available, otherwise fall back to config
        selected_ffs = getattr(self, '_selected_forcefields', None) or {}
        category_order = ['protein', 'modified_aa', 'dna', 'rna', 'carbohydrates',
                         'lipids', 'small_molecules', 'water', 'ions']

        has_leaprc = False
        for category in category_order:
            sel = selected_ffs.get(category)
            if sel is None:
                continue
            if isinstance(sel, list):
                for item in sel:
                    if isinstance(item, dict) and 'leaprc' in item:
                        leaprc = item['leaprc']
                        if isinstance(leaprc, list):
                            for lr in leaprc:
                                script_lines.append(f"source {lr}")
                        else:
                            script_lines.append(f"source {leaprc}")
                        has_leaprc = True
            elif isinstance(sel, dict):
                if 'leaprc' in sel:
                    leaprc = sel['leaprc']
                    if isinstance(leaprc, list):
                        for lr in leaprc:
                            script_lines.append(f"source {lr}")
                    else:
                        script_lines.append(f"source {leaprc}")
                    has_leaprc = True
                if 'frcmod' in sel and sel['frcmod']:
                    frcmod = sel['frcmod']
                    if isinstance(frcmod, list):
                        for fm in frcmod:
                            script_lines.append(f'loadamberparams "{fm}"')
                    else:
                        script_lines.append(f'loadamberparams "{frcmod}"')

        # Fallback if no FF selection was made
        if not has_leaprc:
            script_lines.append(f"source leaprc.protein.{self.config.ffprot}")
            script_lines.append(f"source leaprc.water.{self.config.effective_water_model}")

        # Atom types
        if atom_type_lines:
            script_lines.append("")
            script_lines.append("# Custom atom types")
            type_entries = " ".join(atom_type_lines)
            script_lines.append(f"addAtomTypes {{ {type_entries} }}")

        # Load custom libs and frcmods
        for lib_file in lib_files:
            ext = os.path.splitext(lib_file)[1].lower()
            if ext == ".lib":
                script_lines.append(f'loadoff "{lib_file}"')
            elif ext == ".mol2":
                script_lines.append(f'loadmol2 "{lib_file}"')

        for frcmod_file in frcmod_files:
            script_lines.append(f'loadamberparams "{frcmod_file}"')

        # Load structure and save
        script_lines.extend([
            "",
            f"mol = loadpdb {protein_pdb}",
        ])

        # Add bond definitions from workspace if available. Mirror the Topology
        # Generator's category ordering and coverage so the hydrogen-addition
        # topology has the same metal-ligand, disulfide, and peptide-backbone
        # bonds the final topology will.
        bond_commands = workspace.get("combined_tleap_commands", {})
        if bond_commands:
            for bond_type in ["covalent", "disulfide", "coordinate",
                              "metal-metal", "peptide_backbone", "other"]:
                for cmd in bond_commands.get(bond_type, []):
                    script_lines.append(cmd)

        script_lines.extend([
            "",
            f"saveamberparm mol {output_prefix}.prmtop {output_prefix}.rst7",
            "quit",
        ])

        script_content = "\n".join(script_lines)
        script_path = os.path.join(work_dir, "pretleap_hydrogen.in")

        with open(script_path, "w") as f:
            f.write(script_content)

        self.console.print(f"[grey50]tLEaP script: {script_path}[/grey50]")

        # Run tLEaP
        tleap_exe = "tleap"
        try:
            proc = subprocess.run(
                [tleap_exe, "-f", script_path],
                capture_output=True, text=True, timeout=120,
                cwd=work_dir
            )
        except FileNotFoundError:
            self.console.print("[red]tleap not found on PATH — skipping pre-tLEaP hydrogen pass[/red]")
            return None
        except subprocess.TimeoutExpired:
            self.console.print("[red]tLEaP timed out — skipping pre-tLEaP hydrogen pass[/red]")
            return None

        prmtop_file = f"{output_prefix}.prmtop"
        rst7_file = f"{output_prefix}.rst7"

        if not os.path.exists(prmtop_file) or not os.path.exists(rst7_file):
            self.console.print("[red]tLEaP did not produce topology files — skipping[/red]")
            if proc.stderr:
                # Show last few lines of stderr for diagnosis
                err_lines = proc.stderr.strip().split("\n")[-5:]
                for line in err_lines:
                    self.console.print(f"[red]  {line}[/red]")
            return None

        # Convert back to PDB with ambpdb
        h_pdb = os.path.join(work_dir, "protein_with_h.pdb")
        try:
            proc = subprocess.run(
                ["ambpdb", "-p", prmtop_file, "-c", rst7_file],
                capture_output=True, text=True, timeout=60,
                cwd=work_dir
            )
            if proc.returncode == 0 and proc.stdout:
                with open(h_pdb, "w") as f:
                    f.write(proc.stdout)
                self.console.print(f"[green]Hydrogens added via tLEaP → {os.path.basename(h_pdb)}[/green]")
                return h_pdb
            else:
                self.console.print("[red]ambpdb failed — skipping pre-tLEaP hydrogen pass[/red]")
                return None
        except FileNotFoundError:
            self.console.print("[red]ambpdb not found on PATH — skipping pre-tLEaP hydrogen pass[/red]")
            return None

    def _run_build(self, workspace) -> bool:
        """Execute the membrane build."""
        from .packmol_runner import find_packmol_memgen, run_packmol_memgen

        # Show review first
        self._show_review(workspace)

        proceed = confirm_with_context(
            self.processor,
            "\nProceed with build?",
            default=True,
            module=MODULE_NAME,
            description="Confirm build",
        )

        if not proceed:
            return False

        # Empty-bilayer guard: packmol-memgen has no protein to scale the
        # membrane patch around, so --distxy_fix must be set.
        if not self.config.protein_pdb and self.config.fixed_xy is None:
            self.console.print(
                "[red]Empty bilayer requires an x/y patch size. "
                "Set it in option 4 (Box & Membrane Dimensions) or "
                "re-enter option 1 (Protein Selection → No protein).[/red]"
            )
            return False

        # Check executable
        exe = find_packmol_memgen()
        if exe is None:
            self.console.print(
                "[red]packmol-memgen not found. Ensure AmberTools is installed "
                "and $AMBERHOME/bin is on your PATH.[/red]"
            )
            return False

        self.console.print(f"\n[#0f7f99]Using: {exe}[/#0f7f99]")

        # Determine working directory
        if self.config.protein_pdb:
            work_dir = str(Path(self.config.protein_pdb).parent)
        else:
            work_dir = os.getcwd()

        # tLEaP hydrogen pass: add H atoms using AMBER library templates.
        # This replaces reduce (packmol-memgen's default), ensuring correct
        # hydrogen placement for both standard and custom residues.
        # Empty-bilayer builds skip this — there's no protein to protonate
        # and lipid templates already carry hydrogens.
        if self.config.protein_pdb:
            h_pdb = self._run_pre_tleap_hydrogen_pass(workspace, work_dir)
            if h_pdb:
                self.config.protein_pdb = h_pdb
            else:
                self.console.print("[dark_orange3]tLEaP hydrogen pass failed — falling back to reduce[/dark_orange3]")
                self.config.skip_protonation = False

        # Build CLI args
        args = self.config.to_cli_args()

        self.console.print(f"[#0f7f99]Working directory: {work_dir}[/#0f7f99]")
        self.console.print("[#0f7f99]Running packmol-memgen...[/#0f7f99]\n")

        # Run
        result = run_packmol_memgen(args, work_dir, self.console)

        if not result.success:
            self.console.print(f"\n[red]Build failed: {result.error_message}[/red]")
            if result.log_file:
                self.console.print(f"[grey50]Log file: {result.log_file}[/grey50]")
            return False

        # Show warnings summary if any
        if result.warnings:
            self.console.print(
                Panel(
                    "\n".join(w.strip() for w in result.warnings if w.strip()),
                    title="Warnings from packmol-memgen",
                    border_style="yellow",
                    expand=False,
                )
            )

        # Post-process: remove spurious TER records
        # packmol-memgen inserts TER records around non-standard residues
        # (e.g., HBO, FHO, RHO, PRN) that break tLEaP's chain parsing.
        # Remove TERs between protein/cofactor atoms (same chain, non-lipid).
        if result.output_pdb and os.path.exists(result.output_pdb):
            removed = self._fix_ter_records(result.output_pdb)
            if removed > 0:
                self.console.print(f"[grey50]Removed {removed} spurious TER record(s) from output PDB[/grey50]")

        # Success — update workspace
        self.console.print(f"\n[green]Build successful![/green]")
        self.console.print(f"  Output PDB: {result.output_pdb}")

        # Auto-fired post-build refresh of the assembled membrane-protein
        # system. No force= so it stays silent in CLI when no viewer is
        # open; pushes live to one that is. Web shell auto-launches via env
        # flag. Completing a build is a milestone, not an explicit "show me
        # the viewer" command, so an unbidden browser pop would be intrusive.
        if result.output_pdb and os.path.exists(result.output_pdb):
            try:
                from proprep.structure_prep.viewer_coordinator import (
                    viewer as _viewer,
                )
                _viewer.show_structure(result.output_pdb)
            except Exception as exc:
                logger.debug("post-build viewer hook silenced: %s", exc)

        if result.box_dimensions:
            self.console.print(
                f"  Box: {result.box_dimensions[0]:.1f} × "
                f"{result.box_dimensions[1]:.1f} × "
                f"{result.box_dimensions[2]:.1f} Å"
            )
        if result.water_count:
            self.console.print(f"  Water molecules: {result.water_count:,}")
        if result.ion_counts:
            ion_str = ", ".join(f"{count} {name}" for name, count in result.ion_counts.items())
            self.console.print(f"  Ions: {ion_str}")

        # Write workspace keys
        self.update_workspace(workspace, "membrane_packed_pdb", result.output_pdb)
        self.update_workspace(workspace, "membrane_config", self.config.to_dict())
        self.update_workspace(workspace, "membrane_leaprc_requirements", self.config.get_leaprc_requirements())
        self.update_workspace(workspace, "membrane_box_dimensions", result.box_dimensions)
        self.update_workspace(workspace, "membrane_ion_summary", {
            "ion_counts": result.ion_counts or {},
            "water_count": result.water_count or 0,
            "net_charge_delta": self.config.total_charge_delta,
        })
        self.update_workspace(workspace, "membrane_solutes",
                              [s.to_dict() for s in self.config.solutes])

        # Wire user-supplied solute force-field files through to the topology
        # stage. The Topology Generator loads preprocessing_lib_files /
        # preprocessing_frcmod_files via loadoff / loadamberparams, so append
        # each solute's lib/frcmod there (deduped). Without this the files the
        # user provided under "Custom Parameters & Solutes" are recorded but
        # never loaded, and a non-standard solute goes unparameterized at
        # tLEaP. Solutes left without files fall through to ProPrep's normal
        # non-standard-residue handling downstream.
        solute_libs = [s.lib for s in self.config.solutes if s.lib]
        solute_frcmods = [s.frcmod for s in self.config.solutes if s.frcmod]
        if solute_libs:
            libs = list(workspace.get("preprocessing_lib_files", []) or [])
            for f in solute_libs:
                if f not in libs:
                    libs.append(f)
            self.update_workspace(workspace, "preprocessing_lib_files", libs)
        if solute_frcmods:
            frcmods = list(workspace.get("preprocessing_frcmod_files", []) or [])
            for f in solute_frcmods:
                if f not in frcmods:
                    frcmods.append(f)
            self.update_workspace(workspace, "preprocessing_frcmod_files", frcmods)

        self.update_workspace(workspace, "is_membrane_system", True)

        self._last_result = result

        self.console.print(
            "\n[green]Workspace updated. Next steps:[/green]\n"
            "  1. Run [bold]Topology Generator[/bold] to parametrize the system\n"
            "  2. Run [bold]Molecular Dynamics Manager[/bold] for minimization/equilibration"
        )

        return True
