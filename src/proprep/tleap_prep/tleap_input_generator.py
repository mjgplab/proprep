"""
Topology Generator Module

Generates tLEaP input files for MD simulations with proper bond definitions
for disulfide bonds, metal coordination, and other special bonds.
"""

import logging
import os
import parmed as pmd
import subprocess
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, List, Tuple

from rich.panel import Panel
from proprep.utils.prompts import prompt_with_context, confirm_with_context
from proprep.utils.file_browser import (
    remap_recorded_index, annotate_selected_path,
    remap_recorded_index_by_key, annotate_recorded_key,
)
from rich.table import Table

from parmed.tools.checkvalidity import check_validity

from proprep.utils.module_registry import ProcessingModule, register_module
from proprep.utils.tleap_utils import tleap_safe_unit_var
from proprep.forcefield_params.forcefield_catalog import FORCEFIELD_OPTIONS as _FF_CATALOG
from proprep.forcefield_params.forcefield_catalog import (
    recommended_water_for_protein,
    recommended_ions_for_water,
)
from proprep.utils.prompts import prompt_with_context, confirm_with_context
from proprep.tleap_prep.tleap_commands import (
    GatherBondDefinitionsCommand,
    EditCombinedBondsCommand,
    GenerateSingleStateTLeapCommand,
    GenerateMicrostateInputsCommand,
    GenerateTopologyFilesCommand,
    GenerateCpinCommand,
    RunPBTitrateCommand,
)

# Setup logging
logger = logging.getLogger(__name__)


# Counter-ion options exposed in the Topology Generator menu.
# `unit` is the tLEaP unit name from atomic_ions.lib (e.g., LI, MG, IOD use the
# legacy uppercase names; Na+, K+, Cl- use the modern signed names).
# `charge` is the integer formal charge.
# Every water leaprc except TIP5P loads parameters for all of these via
# atomic_ions.lib + ionsjc/ionslm/ions234lm frcmods (see leaprc.water.*).
_COUNTER_ION_CATIONS = [
    {"label": "Na+",  "unit": "Na+", "charge": +1},
    {"label": "K+",   "unit": "K+",  "charge": +1},
    {"label": "Li+",  "unit": "LI",  "charge": +1},
    {"label": "Cs+",  "unit": "CS",  "charge": +1},
    {"label": "Rb+",  "unit": "RB",  "charge": +1},
    {"label": "Mg2+", "unit": "MG",  "charge": +2},
    {"label": "Ca2+", "unit": "CA",  "charge": +2},
    {"label": "Zn2+", "unit": "ZN",  "charge": +2},
]
_COUNTER_ION_ANIONS = [
    {"label": "Cl-", "unit": "Cl-", "charge": -1},
    {"label": "Br-", "unit": "BR",  "charge": -1},
    {"label": "F-",  "unit": "F",   "charge": -1},
    {"label": "I-",  "unit": "IOD", "charge": -1},
]


def _counter_ion_default(label: str, kind: str) -> dict:
    """Look up a counter-ion entry by label, with a sane fallback.

    `kind` is "cation" or "anion". Used when reconstructing ion identity from
    a saved session that may predate counter-ion selection.
    """
    pool = _COUNTER_ION_CATIONS if kind == "cation" else _COUNTER_ION_ANIONS
    for entry in pool:
        if entry["label"] == label:
            return dict(entry)
    return dict(pool[0])


def _salt_formula_units(cation_charge: int, anion_charge: int) -> tuple:
    """Return (n_cat_per_formula, n_an_per_formula) for a salt M_a X_b.

    For a 1:1 salt (e.g. NaCl) returns (1, 1). For MgCl2 returns (1, 2).
    """
    from math import gcd
    cp = abs(cation_charge)
    ca = abs(anion_charge)
    g = gcd(cp, ca) or 1
    return (ca // g, cp // g)


def _salt_label(cation: dict, anion: dict) -> str:
    """Build a compact chemical label like 'NaCl', 'MgCl2', 'CaBr2'."""
    a, b = _salt_formula_units(cation["charge"], anion["charge"])
    c_sym = cation["label"].rstrip("+").rstrip("0123456789")
    a_sym = anion["label"].rstrip("-").rstrip("0123456789")
    a_str = str(a) if a > 1 else ""
    b_str = str(b) if b > 1 else ""
    return f"{c_sym}{a_str}{a_sym}{b_str}"


def _normalize_salts(solvation_params: dict) -> tuple:
    """Return (salts_list, neutralize_index) from solvation_parameters.

    Accepts both the new multi-salt shape (`salts: [...]`) and the older
    single-salt shape (`cation`, `anion`, `target_molarity`). Old-style
    blobs from saved sessions are wrapped into a 1-element list with
    neutralize_index=0 so downstream code only deals with one shape.

    Each salt entry is `{cation: {label, unit, charge}, anion: {...},
    concentration: float (M)}`. Returns an empty list + index 0 when no
    salt info is present (caller decides whether that's an error).
    """
    if not solvation_params:
        return [], 0

    salts = solvation_params.get('salts')
    if salts:
        idx = solvation_params.get('neutralize_index', 0)
        if idx < 0 or idx >= len(salts):
            idx = 0
        return [dict(s) for s in salts], idx

    # Legacy single-salt shape
    cation = solvation_params.get('cation') or _counter_ion_default("Na+", "cation")
    anion = solvation_params.get('anion') or _counter_ion_default("Cl-", "anion")
    conc = solvation_params.get('target_molarity', 0.15)
    return [{
        'cation': dict(cation),
        'anion': dict(anion),
        'concentration': float(conc),
    }], 0


def _salts_summary(salts: list) -> str:
    """Compact one-line summary like '150 mM NaCl + 50 mM MgCl2'."""
    if not salts:
        return "no salt"
    parts = []
    for s in salts:
        conc_mM = s['concentration'] * 1000
        if conc_mM > 0:
            parts.append(f"{conc_mM:.0f} mM {_salt_label(s['cation'], s['anion'])}")
        else:
            parts.append(f"neutralize-only {_salt_label(s['cation'], s['anion'])}")
    return " + ".join(parts)


def _format_buffer_for_tleap(buffer, buffer_xyz=None, oct_diagonal=0.0,
                              use_octahedron=True, iso=False) -> str:
    """Format a buffer argument for solvateBox / solvateOct.

    - Uniform: returns "10.00" (optionally with " iso").
    - Per-axis (solvateBox): returns "{ 10 12 14 }".
    - Per-axis (solvateOct): returns "{ 10 12 14 0 }" (4 numbers; tLEaP uses
      the 4th as a diagonal clearance, with 0 meaning 'just report it').
    `iso` only applies when buffer is a single number; tLEaP requires `iso`
    after the buffer to force an isometric box.
    """
    if buffer_xyz:
        bx, by, bz = (float(x) for x in buffer_xyz)
        if use_octahedron:
            return f"{{ {bx:g} {by:g} {bz:g} {float(oct_diagonal or 0.0):g} }}"
        return f"{{ {bx:g} {by:g} {bz:g} }}"
    suffix = " iso" if iso else ""
    return f"{float(buffer):.2f}{suffix}"


def _buffer_describe(buffer, buffer_xyz=None, iso=False) -> str:
    """Human-readable buffer summary for status lines."""
    if buffer_xyz:
        bx, by, bz = buffer_xyz
        return f"{bx:g}/{by:g}/{bz:g} A (x/y/z)"
    if iso:
        return f"{float(buffer):.1f} A iso"
    return f"{float(buffer):.1f} A"


@register_module
class TLeapInputGenerator(ProcessingModule):
    """Module for generating tLEaP input files for MD simulations"""

    NAME = "Topology Generator"
    DESCRIPTION = "Generate topology+coordinate files for MD"
    CATEGORY = "preparation"
    VERSION = "1.0.0"

    def __init__(self):
        """Initialize the Topology Generator module"""
        self.processor = None
        self.combined_bonds = {}
        self.tleap_parameters = {}
        self.output_file = None

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Helper methods to centralize interactions with the workspace

    def get_workspace(self):
        """Get the current workspace object"""
        return self.processor.workspace

    def get_from_workspace(self, key, default=None):
        """Get values from the processor's workspace"""
        return self.processor.workspace.get(key, default)

    def update_workspace(self, key, value):
        """Update the processor's workspace"""
        self.processor.workspace.set(key, value)

    @staticmethod
    def _read_mol2_resname(mol2_file):
        """Return the residue name of the first atom in a mol2 file, or None.

        Used to name the tLEaP unit variable when a mol2 is loaded directly
        (`unit = loadmol2 ...`); the variable name becomes the template's
        match-name for loadpdb.
        """
        try:
            with open(mol2_file) as f:
                in_atoms = False
                for line in f:
                    s = line.strip()
                    if s.startswith('@<TRIPOS>ATOM'):
                        in_atoms = True
                        continue
                    if in_atoms:
                        if s.startswith('@<TRIPOS>'):
                            break
                        parts = s.split()
                        # mol2 ATOM record: id name x y z type subst_id subst_name ...
                        if len(parts) >= 8:
                            return parts[7]
        except Exception:
            pass
        return None

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Module Menu

    def get_menu_options(self) -> Dict[str, str]:
        """Get module menu options"""
        menu = {
            "edit_bonds": "Edit redox site bond definitions",
            "generate_single_state": "Generate tLEaP input for single state",
            "generate_microstate_inputs": "Generate tLEaP inputs for all redox microstates",
            "generate_topology": "Generate prmtop/rst7 files from tLEaP input files",
            "pb_titrate": "Refine protonation states via PB (PBSA)",
            "generate_cpin": "Generate cpin file for constant pH MD",
        }

        return menu

    def get_enhanced_menu_options(self, workspace):
        """
        Get menu options with enhanced status information.
        Menu adapts based on microstate generation choices from Redox Site Prep.

        Args:
            workspace: Current workspace

        Returns:
            List of MenuOption objects with status
        """
        from proprep.utils.enhanced_menu import MenuOption, OptionStatus
        import glob
        import os

        options = []

        # Check workspace state from Redox Site Preparation
        has_single_microstate = workspace.get("transformed_pdb_file") is not None
        has_batch_microstates = workspace.get("generated_microstate_pdbs") is not None

        # Check tLEaP workflow state (single state OR batch microstates)
        single_tleap_file = workspace.get("tleap_input_file")
        microstate_tleap_files = workspace.get("generated_microstate_tleap_files")
        has_tleap_input = (single_tleap_file is not None or microstate_tleap_files is not None)
        output_dir = workspace.get("output_dir", ".")
        prmtop_files = glob.glob(os.path.join(output_dir, "**/*.prmtop"), recursive=True)

        # Determine if topology generation is complete based on context
        if microstate_tleap_files:
            # For batch microstates, check if we have topologies for most/all of them
            num_tleap = len(microstate_tleap_files)
            num_prmtop = len(prmtop_files)
            # Consider complete if we have at least 80% of expected topologies
            has_topology = num_prmtop >= (num_tleap * 0.8)
        else:
            # For single state, any prmtop is fine
            has_topology = len(prmtop_files) > 0

        has_cpin = workspace.get("cpin_file") is not None

        # Option 1: Edit bonds - always available
        options.append(MenuOption(
            key="1",
            description="Edit redox site bond definitions",
            status=OptionStatus.AVAILABLE
        ))

        # Option 2: Generate tLEaP input for single state
        if has_batch_microstates:
            # User already generated batch - block single state option
            options.append(MenuOption(
                key="2",
                description="Generate tLEaP input for single state",
                status=OptionStatus.BLOCKED,
                dependency_text="[Batch microstates already generated] ○"
            ))
        else:
            # Needs a PDB-yielding structure (or packed membrane). NOT the
            # broader can_process, which also passes for a parm7/rst7-only
            # PB-Titrate resume workspace that this single-state path cannot use.
            from proprep.utils.structure_selector import StructureSelector
            _console = self.processor.console if self.processor else None
            has_pdb = (StructureSelector(workspace, _console)
                       .get_structure_status().get("has_any", False)
                       or bool(workspace.get("membrane_packed_pdb")))
            if has_pdb:
                options.append(MenuOption(
                    key="2",
                    description="Generate tLEaP input for single state",
                    status=OptionStatus.COMPLETED if has_tleap_input else OptionStatus.AVAILABLE
                ))
            else:
                options.append(MenuOption(
                    key="2",
                    description="Generate tLEaP input for single state",
                    status=OptionStatus.BLOCKED,
                    dependency_text="Load a structure first"
                ))

        # Option 3: Generate tLEaP inputs for all redox microstates
        if has_batch_microstates:
            options.append(MenuOption(
                key="3",
                description="Generate tLEaP inputs for all redox microstates",
                status=OptionStatus.COMPLETED if has_tleap_input else OptionStatus.AVAILABLE
            ))
        else:
            # User generated single or nothing - block batch option
            options.append(MenuOption(
                key="3",
                description="Generate tLEaP inputs for all redox microstates",
                status=OptionStatus.BLOCKED,
                dependency_text="[Need batch microstate metadata from Redox Site Prep] ○"
            ))

        # Option 4: Generate topology - requires tLEaP input
        if has_tleap_input:
            status = OptionStatus.COMPLETED if has_topology else OptionStatus.READY
            dep_text = ""
        else:
            status = OptionStatus.BLOCKED
            dep_text = "[Need to generate tLEaP input first] ○"

        options.append(MenuOption(
            key="4",
            description="Generate prmtop/rst7 files from tLEaP input files",
            status=status,
            dependency_text=dep_text
        ))

        # Option 5: Refine protonation states via PB (PBSA) - requires topology.
        # Optional refinement step. Recommendations feed into option 6 (cpin)
        # as a third "Use titrate recommendations" choice for initial states.
        has_titrate_recs = workspace.get("titrate_recommendations") is not None
        if has_topology:
            status = (OptionStatus.COMPLETED if has_titrate_recs
                      else OptionStatus.READY)
            dep_text = ""
        else:
            status = OptionStatus.BLOCKED
            dep_text = "[Need to generate topology files first] ○"

        options.append(MenuOption(
            key="5",
            description="Refine protonation states via PB (PBSA)",
            status=status,
            dependency_text=dep_text
        ))

        # Option 6: Generate cpin - requires topology
        if has_topology:
            status = OptionStatus.COMPLETED if has_cpin else OptionStatus.READY
            dep_text = ""
        else:
            status = OptionStatus.BLOCKED
            dep_text = "[Need to generate topology files first] ○"

        options.append(MenuOption(
            key="6",
            description="Generate cpin file for constant pH MD",
            status=status,
            dependency_text=dep_text
        ))

        return options

    def get_menu_suggestion(self, workspace):
        """
        Get a suggestion for the next recommended action.

        Args:
            workspace: Current workspace

        Returns:
            Suggestion text or None
        """
        import glob
        import os

        # Check workspace state
        has_single_microstate = workspace.get("transformed_pdb_file") is not None
        has_batch_microstates = workspace.get("generated_microstate_pdbs") is not None
        microstate_tleap_files = workspace.get("generated_microstate_tleap_files")
        has_tleap_input = (workspace.get("tleap_input_file") is not None or
                          microstate_tleap_files is not None)
        output_dir = workspace.get("output_dir", ".")
        prmtop_files = glob.glob(os.path.join(output_dir, "**/*.prmtop"), recursive=True)
        num_prmtop = len(prmtop_files)

        # Determine if topology generation is complete based on context
        if microstate_tleap_files:
            num_tleap = len(microstate_tleap_files)
            has_topology = num_prmtop >= (num_tleap * 0.8)
        else:
            has_topology = num_prmtop > 0

        has_cpin = workspace.get("cpin_file") is not None

        if not has_tleap_input:
            if has_batch_microstates:
                return "Edit bonds (option 1) if needed, then generate tLEaP inputs for all microstates (option 3)"
            else:
                return "Edit bonds (option 1) if needed, then generate tLEaP input for single state (option 2)"
        elif not has_topology:
            if microstate_tleap_files:
                return f"✓ tLEaP inputs generated ({len(microstate_tleap_files)} files). Generate topology files (option 4) to create prmtop/rst7 files"
            else:
                return "✓ tLEaP input generated. Generate topology files (option 4) to create prmtop/rst7 files"
        elif has_cpin:
            return f"✓ All steps complete ({num_prmtop} topologies, cpin generated). Press [m] to return to the main menu"
        else:
            return (f"✓ Topology files generated ({num_prmtop} prmtop files). "
                    f"Optionally refine protonation states via PB (option 5) or "
                    f"generate cpin for constant pH (option 6), or press [m] to "
                    f"return to the main menu")

    def handle_menu_option(self, option: str) -> bool:
        """Handle a menu option selection"""
        if option == "edit_bonds":
            command = EditCombinedBondsCommand(self.processor)
            return command.execute()
        elif option == "generate_single_state":
            command = GenerateSingleStateTLeapCommand(self.processor)
            return command.execute()
        elif option == "generate_microstate_inputs":
            command = GenerateMicrostateInputsCommand(self.processor)
            return command.execute()
        elif option == "generate_topology":
            command = GenerateTopologyFilesCommand(self.processor)
            return command.execute()
        elif option == "pb_titrate":
            command = RunPBTitrateCommand(self.processor)
            return command.execute()
        elif option == "generate_cpin":
            command = GenerateCpinCommand(self.processor)
            return command.execute()

        return False

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Bond definition gathering and management

    def _display_bond_summary(self):
        """Display a summary of all bond definitions with detailed formatting"""
        total_bonds = sum(len(bonds) for bonds in self.combined_bonds.values())

        if total_bonds == 0:
            self.processor.console.print("[yellow]No bond definitions found.[/yellow]")
            return

        self.processor.console.print("\n[bold]Bond Definitions Summary[/bold]")
        self.processor.console.print(
            f"[green]Successfully gathered {total_bonds} bond definitions:[/green]"
        )

        # Create a table for each category that has bonds
        for category, bonds in self.combined_bonds.items():
            if not bonds:
                continue

            # Create a descriptive title based on category
            if category == "covalent":
                title = "Covalent Bonds (Non-Metal to Non-Metal)"
            elif category == "disulfide":
                title = "Disulfide Bonds (CYS/CYX SG-SG)"
            elif category == "coordinate":
                title = "Coordination Bonds (Metal to Non-Metal)"
            elif category == "metal-metal":
                title = "Metal-Metal Bonds"
            elif category == "peptide_backbone":
                title = "Peptide Backbone Bonds"
            else:
                title = "Other Custom Bonds"

            # Show count in the category header
            self.processor.console.print(
                f"\n[bold blue]{title} ({len(bonds)} bonds)[/bold blue]"
            )

            table = Table(title=title)
            table.add_column("Index", style="blue", width=6)
            table.add_column("Bond Command", style="green")
            table.add_column("Comment", style="yellow")

            # Add each bond to the table
            for i, bond in enumerate(bonds, 1):
                # Split into command and comment if present
                parts = bond.split("#", 1)
                command = parts[0].strip()
                comment = parts[1].strip() if len(parts) > 1 else ""

                table.add_row(str(i), command, comment)

            # Print the table
            self.processor.console.print(table)

            # For large tables, show only the first few bonds
            max_display = 5
            if len(bonds) > max_display:
                self.processor.console.print(
                    f"  ... and {len(bonds) - max_display} more bonds"
                )

    def gather_bond_definitions(self) -> bool:
        """Gather bond definitions from all modules.

        Returns True if the caller should proceed (bonds gathered, or user
        confirmed proceeding without bonds, or membrane-only system where
        bond definitions don't apply). Returns False only if the user
        explicitly aborted.
        """
        workspace = self.get_workspace()

        # Initialize combined bonds dictionary - matches comprehensive_redox_detector categories
        self.combined_bonds = {
            "covalent": [],           # 0 metals (2 non-metals) - excludes disulfides
            "coordinate": [],         # 1 metal + 1 non-metal
            "metal-metal": [],        # 2 metals (matches detector naming)
            "disulfide": [],          # SG-SG bonds between CYS residues
            "peptide_backbone": [],   # Generated for non-standard residues
            "other": [],              # Legacy/manual bonds
        }

        # Collect bonds from comprehensive RedoxSite objects
        redox_sites = self.get_from_workspace("detected_redox_sites")
        if redox_sites:
            # Diagnostic: show what we found
            for _i, _site in enumerate(redox_sites):
                _n_bonds = len(_site.bonds) if hasattr(_site, 'bonds') else 0
                _n_coords = len(_site.coord_to_pdb) if hasattr(_site, 'coord_to_pdb') else 0
                _site_id = _site.site_id if hasattr(_site, 'site_id') else f'site_{_i}'
                self.processor.console.print(
                    f"[grey50]  {_site_id}: {_n_bonds} bonds, {_n_coords} coord_to_pdb entries[/grey50]"
                )
            self.processor.console.print(f"[grey50]Processing {len(redox_sites)} RedoxSite objects for tLEaP bond definitions...[/grey50]")
            redox_bond_commands = self._convert_redox_sites_to_tleap_commands(redox_sites)
            
            # Merge RedoxSite bonds into combined_bonds
            for category, commands in redox_bond_commands.items():
                if commands:  # Only add if there are commands
                    self.combined_bonds[category].extend(commands)

            # Display detailed bond summary
            total_bonds = sum(len(bonds) for bonds in self.combined_bonds.values())
            if total_bonds > 0:
                # Use the improved display method
                self._display_bond_summary()

            # Store in workspace
            self.update_workspace("combined_tleap_commands", self.combined_bonds)
            return True
        else:
            # Membrane-only systems can't have redox sites (no protein), so
            # there's nothing to detect and nothing to prompt about — proceed
            # silently with empty bonds.
            if self.get_from_workspace("is_membrane_system", False) and \
                    not self.get_from_workspace("membrane_config", {}).get("protein_pdb"):
                self.processor.console.print(
                    "[grey50]Empty bilayer — no bond definitions needed.[/grey50]"
                )
                return True

            self.processor.console.print(
                "[yellow]No bond definitions found. Run comprehensive RedoxSite detection first.[/yellow]"
            )

            # Ask if user wants to continue anyway
            if confirm_with_context(
                self.processor,
                "Do you want to continue without bond definitions?",
                default=False,
                module="Topology Generator",
                description="Continue without bond definitions",
            ):
                self.processor.console.print(
                    "[blue]Continuing with empty bond definitions.[/blue]"
                )
                return True
            else:
                self.processor.console.print(
                    "[yellow]Please run RedoxSite detection and try again.[/yellow]"
                )
                return False

    def edit_combined_bonds(self):
        """Edit the combined bond definitions using single-page dashboard UX"""
        # Load bonds if not already loaded
        if not self.combined_bonds:
            self.combined_bonds = self.get_from_workspace("combined_tleap_commands", {})

            if not self.combined_bonds:
                self.processor.console.print(
                    "[yellow]No bond definitions found. Gathering from modules...[/yellow]"
                )
                self.gather_bond_definitions()
                if not self.combined_bonds:
                    self.processor.console.print(
                        "[yellow]No bond definitions found in any module.[/yellow]"
                    )
                    return

        # Main dashboard loop
        filter_category = None  # None = show all, or specific category name
        
        while True:
            # Clear screen for dashboard effect
            self.processor.console.clear()
            
            # Display dashboard
            self._display_bond_dashboard(filter_category)
            
            # Get user command
            command = prompt_with_context(
                self.processor,
                "\n[bold blue]Command[/bold blue]",
                default="q",
                module="Topology Generator",
                description="Bond dashboard command",
                options_map={
                    "q": "Quit dashboard",
                    "c": "Filter: covalent",
                    "m": "Filter: coordinate",
                    "mm": "Filter: metal-metal",
                    "p": "Filter: peptide backbone",
                    "o": "Filter: other",
                    "all": "Clear filter",
                    "a": "Add bond interactively",
                    "s": "Save bonds to workspace",
                    "h": "Help",
                    "#": "Edit bond by number",
                    "x#": "Delete bond by number",
                },
            ).lower().strip()
            
            # Handle commands
            if command == 'q' or command == 'quit':
                break
            elif command == 'c':
                filter_category = 'covalent'
            elif command == 'd':
                filter_category = 'disulfide'
            elif command == 'm':
                filter_category = 'coordinate'
            elif command == 'mm':
                filter_category = 'metal-metal'
            elif command == 'p':
                filter_category = 'peptide_backbone'
            elif command == 'o':
                filter_category = 'other'
            elif command == 'all' or command == '':
                filter_category = None
            elif command == 'a' or command == 'add':
                self._add_bond_interactive()
            elif command == 's' or command == 'save':
                self._save_bonds()
                self.processor.console.print("[green]✓ Bonds saved to workspace[/green]")
                prompt_with_context(
                    self.processor,
                    "Press Enter to continue...",
                    default="",
                    module="Topology Generator",
                    description="Pause after dashboard action",
                )
            elif command.isdigit():
                # Edit specific bond by number
                bond_num = int(command)
                self._edit_bond_by_number(bond_num, filter_category)
            elif command.startswith('x') and len(command) > 1 and command[1:].isdigit():
                # Delete bond by number (x3 = delete bond 3)
                bond_num = int(command[1:])
                self._delete_bond_by_number(bond_num, filter_category)
            elif command in ['h', 'help']:
                self._show_dashboard_help()
            else:
                self.processor.console.print(f"[red]Unknown command: {command}[/red]")
                self.processor.console.print("Type 'h' for help or 'q' to quit")
                prompt_with_context(
                    self.processor,
                    "Press Enter to continue...",
                    default="",
                    module="Topology Generator",
                    description="Pause after dashboard action",
                )
        
        # Save bonds when exiting (ensures workspace is always updated)
        self._save_bonds()

    def _handle_global_actions(self):
        """Handle global actions that apply to all categories."""
        while True:
            self.processor.console.print("\n[bold]Global Bond Actions[/bold]")
            self.processor.console.print("1. Manage intra-residue bonds (view and delete across all categories)", highlight=False)
            self.processor.console.print("2. Back to category selection", highlight=False)

            choice = prompt_with_context(
                self.processor,
                "Enter your choice",
                choices=["1", "2"],
                default="1",
                module="Topology Generator",
                description="Global bond actions menu",
                options_map={
                    "1": "Manage intra-residue bonds (view and delete across all categories)",
                    "2": "Back to category selection",
                },
            )

            if choice == "1":
                self._manage_global_intra_residue_bonds()
            elif choice == "2":
                break

    def _manage_global_intra_residue_bonds(self):
        """Manage intra-residue bonds across all categories."""
        # Collect all intra-residue bonds from all categories
        intra_bonds = {}
        total_intra = 0
        
        for category, bonds in self.combined_bonds.items():
            category_intra = []
            for bond in bonds:
                if self._is_intra_residue_bond(bond):
                    category_intra.append(bond)
                    total_intra += 1
            if category_intra:
                intra_bonds[category] = category_intra

        if total_intra == 0:
            self.processor.console.print("[green]No intra-residue bonds found across all categories.[/green]")
            return

        # Display all intra-residue bonds organized by category
        self.processor.console.print(f"\n[bold]Found {total_intra} intra-residue bonds:[/bold]")
        
        for category, bonds in intra_bonds.items():
            # Create a descriptive title based on category
            if category == "covalent":
                title = "Covalent Bonds (Non-Metal to Non-Metal)"
            elif category == "disulfide":
                title = "Disulfide Bonds (CYS/CYX SG-SG)"
            elif category == "coordinate":
                title = "Coordination Bonds (Metal to Non-Metal)"
            elif category == "metal-metal":
                title = "Metal-Metal Bonds"
            elif category == "peptide_backbone":
                title = "Peptide Backbone Bonds"
            else:
                title = "Other Custom Bonds"
                
            self.processor.console.print(f"\n[blue]{title} ({len(bonds)} intra-residue bonds):[/blue]")
            for i, bond in enumerate(bonds, 1):
                self.processor.console.print(f"  {i}. {bond}")

        # Ask if user wants to delete all intra-residue bonds
        if confirm_with_context(
            self.processor,
            f"Delete all {total_intra} intra-residue bonds?",
            default=False,
            module="Topology Generator",
            description="Delete all intra-residue bonds across categories",
        ):
            deleted_count = 0
            for category in self.combined_bonds:
                original_count = len(self.combined_bonds[category])
                self.combined_bonds[category] = [
                    bond for bond in self.combined_bonds[category] 
                    if not self._is_intra_residue_bond(bond)
                ]
                deleted_from_category = original_count - len(self.combined_bonds[category])
                deleted_count += deleted_from_category

            self.processor.console.print(f"[green]Deleted {deleted_count} intra-residue bonds across all categories.[/green]")
            # Update workspace
            self.update_workspace("combined_tleap_commands", self.combined_bonds)
        else:
            self.processor.console.print("[yellow]No bonds were deleted.[/yellow]")

    def _is_intra_residue_bond(self, bond):
        """Check if a bond is an intra-residue bond."""
        import re
        # Pattern to extract residue IDs from bond commands
        # Matches: bond mol.RESID.ATOM mol.RESID.ATOM
        bond_pattern = re.compile(r'bond\s+mol\.(\d+)\.\w+\s+mol\.(\d+)\.\w+')
        
        match = bond_pattern.search(bond)
        if match:
            resid1, resid2 = match.groups()
            return resid1 == resid2
        return False

    def _display_bond_summary(self):
        """Display a summary of all bond definitions"""
        self.processor.console.print("\n[bold]Bond Definitions Summary[/bold]")

        total_bonds = sum(len(bonds) for bonds in self.combined_bonds.values())
        self.processor.console.print(f"Total bond definitions: {total_bonds}")

        for category, bonds in self.combined_bonds.items():
            if bonds:  # Only show categories with bonds
                self.processor.console.print(f"{category}: {len(bonds)} bonds")

    def _view_bonds(self, category):
        """View bonds in a specific category"""
        bonds = self.combined_bonds.get(category, [])

        if not bonds:
            self.processor.console.print(
                f"[yellow]No bonds in category: {category}[/yellow]"
            )
            return

        table = Table(title=f"{category} Bonds")
        table.add_column("Index", style="blue", width=6)
        table.add_column("Bond Command", style="green")
        table.add_column("Comment", style="yellow")

        for i, bond in enumerate(bonds, 1):
            # Split into command and comment if present
            parts = bond.split("#", 1)
            command = parts[0].strip()
            comment = parts[1].strip() if len(parts) > 1 else ""

            table.add_row(str(i), command, comment)

        self.processor.console.print(table)
        prompt_with_context(
            self.processor,
            "Press Enter to continue...",
            default="",
            module="Topology Generator",
            description="Pause after viewing bonds in category",
        )

    def _add_bond(self, category):
        """Add a new bond to a category"""
        self.processor.console.print(f"\n[bold]Add Bond to {category}[/bold]")

        # Provide example based on category
        if category == "covalent":
            self.processor.console.print("Example: bond mol.145.SG mol.167.SG")
        elif category == "disulfide":
            self.processor.console.print("Example: bond mol.6.SG mol.127.SG")
        elif category == "coordinate":
            self.processor.console.print("Example: bond mol.204.ZN mol.201.ND1")
        elif category == "metal-metal":
            self.processor.console.print("Example: bond mol.204.ZN mol.205.FE")
        elif category == "peptide_backbone":
            self.processor.console.print("Example: bond mol.156.C mol.157.N")
        else:
            self.processor.console.print("Example: bond mol.100.X mol.200.Y")

        # Get bond command
        bond_cmd = prompt_with_context(
            self.processor,
            "Enter bond command",
            module="Topology Generator",
            description=f"Add bond to {category} category",
        )

        # Ensure it starts with 'bond'
        if not bond_cmd.strip().lower().startswith("bond "):
            bond_cmd = "bond " + bond_cmd

        # Get optional comment
        comment = prompt_with_context(
            self.processor,
            "Enter comment (optional)",
            default="",
            module="Topology Generator",
            description=f"Comment for new {category} bond",
        )
        if comment:
            bond_cmd += f" # {comment}"

        # Add to the appropriate category
        self.combined_bonds[category].append(bond_cmd)
        self.processor.console.print(f"[green]Added bond to {category}[/green]")

    def _edit_bond(self, category):
        """Edit an existing bond"""
        bonds = self.combined_bonds.get(category, [])

        if not bonds:
            self.processor.console.print(
                f"[yellow]No bonds in category: {category}[/yellow]"
            )
            return

        # Show bonds for selection
        self._view_bonds(category)

        # Get bond index to edit
        max_idx = len(bonds)
        edit_options_map = {"0": "Cancel"}
        for i, bond in enumerate(bonds, 1):
            edit_options_map[str(i)] = bond
        idx = prompt_with_context(
            self.processor,
            "Enter index of bond to edit (or 0 to cancel)",
            choices=["0"] + [str(i) for i in range(1, max_idx + 1)],
            default="0",
            module="Topology Generator",
            description=f"Select bond to edit in {category} category",
            options_map=edit_options_map,
        )

        if idx == "0":
            return

        # Get the current bond
        bond_idx = int(idx) - 1
        current_bond = bonds[bond_idx]

        # Split into command and comment
        parts = current_bond.split("#", 1)
        command = parts[0].strip()
        comment = parts[1].strip() if len(parts) > 1 else ""

        # Get updated values
        self.processor.console.print(f"Current command: {command}")
        new_command = prompt_with_context(
            self.processor,
            "Enter new command",
            default=command,
            module="Topology Generator",
            description=f"New bond command for bond {idx} in {category}",
        )

        self.processor.console.print(f"Current comment: {comment}")
        new_comment = prompt_with_context(
            self.processor,
            "Enter new comment",
            default=comment,
            module="Topology Generator",
            description=f"New bond comment for bond {idx} in {category}",
        )

        # Update the bond
        if new_comment:
            updated_bond = f"{new_command} # {new_comment}"
        else:
            updated_bond = new_command

        self.combined_bonds[category][bond_idx] = updated_bond
        self.processor.console.print(f"[green]Updated bond {idx}[/green]")

    def _delete_bond(self, category):
        """Delete a bond"""
        bonds = self.combined_bonds.get(category, [])

        if not bonds:
            self.processor.console.print(
                f"[yellow]No bonds in category: {category}[/yellow]"
            )
            return

        # Show bonds for selection
        self._view_bonds(category)

        # Get bond index to delete
        max_idx = len(bonds)
        self.processor.console.print("\n[grey50]Tip: Enter -1 to remove all intra-residue bonds[/grey50]")
        delete_options_map = {"0": "Cancel", "-1": "Delete all intra-residue bonds"}
        for i, bond in enumerate(bonds, 1):
            delete_options_map[str(i)] = bond
        idx = prompt_with_context(
            self.processor,
            "Enter index of bond to delete (or 0 to cancel, -1 for intra-residue)",
            module="Topology Generator",
            description=f"Select bond to delete in {category} category",
            options_map=delete_options_map,
        )

        # Handle special case for intra-residue bonds
        if idx == "-1":
            self._delete_intra_residue_bonds(category)
            return
        
        # Validate input for normal deletion
        if idx not in ["0"] + [str(i) for i in range(1, max_idx + 1)]:
            self.processor.console.print("[red]Invalid index[/red]")
            return

        if idx == "0":
            return

        # Confirm deletion
        bond_idx = int(idx) - 1
        bond_to_delete = bonds[bond_idx]

        if confirm_with_context(
            self.processor,
            f"Are you sure you want to delete: {bond_to_delete}?",
            default=False,
            module="Topology Generator",
            description=f"Confirm delete bond from {category}",
        ):
            deleted = self.combined_bonds[category].pop(bond_idx)
            self.processor.console.print(f"[green]Deleted bond: {deleted}[/green]")
        else:
            self.processor.console.print("[yellow]Deletion cancelled[/yellow]")
    
    def _delete_intra_residue_bonds(self, category):
        """Delete all bonds between atoms of the same residue"""
        import re
        
        bonds = self.combined_bonds.get(category, [])
        if not bonds:
            self.processor.console.print(
                f"[yellow]No bonds in category: {category}[/yellow]"
            )
            return
        
        # Pattern to extract residue IDs from bond commands
        # Matches: bond mol.RESID.ATOM mol.RESID.ATOM
        bond_pattern = re.compile(r'bond\s+mol\.(\d+)\.\w+\s+mol\.(\d+)\.\w+')
        
        # Identify intra-residue bonds
        intra_residue_bonds = []
        for i, bond in enumerate(bonds):
            match = bond_pattern.search(bond)
            if match:
                resid1, resid2 = match.groups()
                if resid1 == resid2:
                    intra_residue_bonds.append((i, bond, resid1))
        
        if not intra_residue_bonds:
            self.processor.console.print(
                "[yellow]No intra-residue bonds found in this category[/yellow]"
            )
            return
        
        # Display bonds to be removed
        self.processor.console.print(f"\n[bold]Found {len(intra_residue_bonds)} intra-residue bond(s):[/bold]")
        table = Table(show_header=True, header_style="bold magenta")
        table.add_column("Index", style="blue")
        table.add_column("Residue ID", style="yellow")
        table.add_column("Bond Command", style="white")
        
        for idx, bond, resid in intra_residue_bonds:
            table.add_row(str(idx + 1), resid, bond.split('#')[0].strip())
        
        self.processor.console.print(table)
        
        # Confirm deletion
        if confirm_with_context(
            self.processor,
            f"\nRemove all {len(intra_residue_bonds)} intra-residue bonds?",
            default=False,
            module="Topology Generator",
            description=f"Remove all intra-residue bonds in {category}",
        ):
            # Remove bonds in reverse order to maintain indices
            for idx, bond, resid in reversed(intra_residue_bonds):
                self.combined_bonds[category].pop(idx)
            
            self.processor.console.print(
                f"[green]Removed {len(intra_residue_bonds)} intra-residue bond(s)[/green]"
            )
            
            # Show remaining bonds count
            remaining = len(self.combined_bonds[category])
            self.processor.console.print(
                f"[grey50]{remaining} bond(s) remaining in {category}[/grey50]"
            )
        else:
            self.processor.console.print("[yellow]Deletion cancelled[/yellow]")

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # tLEaP configuration and file generation

    def _find_transformed_pdb_files(self):
        """Find transformed PDB files matching *_*_transformed.pdb pattern"""
        import glob
        from pathlib import Path
        
        # Look for transformed PDB files
        pattern = "*_*_transformed.pdb"
        transformed_files = glob.glob(pattern)
        
        return transformed_files
    
    def _select_priority_pdb_file(self, silent=False):
        """Select the appropriate PDB file with priority: transformed > protonation-updated > repaired > filtered > interactive"""
        console = self.processor.console
        workspace = self.get_workspace()

        # Priority order: reordered -> preprocessing -> protonation-updated -> transformed -> repaired -> filtered -> interactive
        # reordered_pdb_file is highest because reordering is the last structure transformation
        # protonation_pdb_file > transformed_pdb_file because protonation runs after
        # transformation and may rename residues (e.g., ASP->AS4 for constant pH)
        priority_keys = [
            # Highest: a PB-titrate-renamed PDB written for a modern-FF rebuild
            # (standard residue names encoding the recommended protonation
            # states). This is the explicit "build production topology from
            # these states" hand-off, so it outranks all prep-step structures.
            ("pb_rename_pdb_file", "PB-titrate protonation-renamed"),
            ("membrane_packed_pdb", "membrane-packed"),  # full membrane system from membrane builder
            ("reordered_pdb_file", "reordered"),  # From structure preparation step
            ("prepared_pdb", "prepared (MCPB preprocessing)"),
            ("preprocessing_protein_input", "preprocessing protein"),
            ("protonation_pdb_file", "protonation-updated"),
            ("transformed_pdb_file", "transformed"),
            ("repaired_pdb_file", "repaired"),
            ("filtered_pdb_file", "filtered"),
            # Fallbacks: untouched user-loaded structures. Picked up when none of
            # the processing-step keys above are populated (e.g. user reloaded a
            # post-transformed PDB and is going straight to topology generation).
            ("local_pdb_file", "user-loaded local PDB"),
            ("rcsb_pdb_file", "user-loaded RCSB PDB"),
        ]

        for workspace_key, description in priority_keys:
            pdb_file = self.get_from_workspace(workspace_key)
            if pdb_file and os.path.exists(pdb_file):
                if not silent:
                    console.print(f"[green]Using {description} PDB file: {pdb_file}[/green]")
                return pdb_file

        # No processed structure found - ask user to select interactively
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(workspace, console, self.processor)
        pdb_file = selector.get_structure(interactive=True)

        if pdb_file:
            # StructureSelector already printed the selection message
            return pdb_file

        # No valid structure found
        console.print("[red]No valid PDB file found in workspace[/red]")
        console.print("[yellow]Please load a PDB file first[/yellow]")
        return None

    def _provisional_solvation_needs_prompt(self, solvation_params, is_membrane):
        """Whether a provisional (preprocessing-origin) implicit solvation should
        be cleared so the user is prompted for a real solvation choice.

        Returns False for membrane systems, while preprocessing runs its own
        metal-free tleap (`_preprocessing_tleap_active` set — that build is
        intentionally box-less), and when the solvation is a genuine user choice
        (no `provisional` marker). Returns True only for a user-initiated
        production build that inherited preprocessing's placeholder implicit.
        """
        if is_membrane:
            return False
        if self.get_from_workspace("_preprocessing_tleap_active", False):
            return False
        return isinstance(solvation_params, dict) and bool(solvation_params.get("provisional"))

    def _reconcile_template_loadpdb(self, template):
        """Rewrite a cached template's loadpdb path if the priority structure drifted.

        The tLEaP template caches its `mol = loadpdb <path>` line as literal text.
        When the template is reused (user answered "no" to Reconfigure), that path
        can be stale — e.g. it was resolved before MCPB Force Field Integration wrote
        the renamed ``prepared_structure.pdb``. This re-resolves the current priority
        PDB and, only if it genuinely differs from the baked-in path, rewrites the
        load line (and the saveamberparm output prefix that was derived from it),
        surfacing the change to the user rather than silently building from an
        outdated (metal-free / un-renamed) structure.

        Returns the (possibly rewritten) template; saves it back to the workspace
        when changed.
        """
        if not template:
            return template

        selected_pdb = self._select_priority_pdb_file(silent=True)
        if not selected_pdb:
            return template

        new_path = os.path.abspath(selected_pdb)
        new_prefix = os.path.splitext(os.path.basename(new_path))[0]

        import re
        m = re.search(r'^(\s*mol\s*=\s*loadpdb\s+)(\S+)\s*$', template, re.MULTILINE)
        if not m:
            return template

        old_path = m.group(2)
        if os.path.abspath(old_path) == new_path:
            return template  # already current — no-op

        console = self.processor.console
        console.print(
            "[yellow]Cached tLEaP template referenced a different structure than the "
            "current priority PDB:[/yellow]\n"
            f"  [grey50]was:[/grey50] {old_path}\n"
            f"  [green]now:[/green] {new_path}"
        )
        console.print(
            "[grey50]Updating the loadpdb path (and output prefix) so the reused "
            "template tracks the current structure.[/grey50]"
        )

        updated = template[:m.start()] + m.group(1) + new_path + template[m.end():]
        # Rewrite the saveamberparm output prefix (derived from the old basename)
        # to match the new structure, keeping output filenames consistent.
        updated = re.sub(
            r'(saveamberparm\s+mol\s+)(\S+)(\.prmtop\s+)(\S+)(\.rst7)',
            lambda mm: mm.group(1) + new_prefix + mm.group(3) + new_prefix + mm.group(5),
            updated,
        )
        self._save_tleap_template(updated)
        return updated

    def _configure_solvation_parameters(self) -> dict:
        """
        Interactive configuration of solvation parameters.

        Prompts the user for:
        1. Solvent model (explicit vs implicit)
        2. Box shape (truncated octahedron vs rectangular) - explicit only
        3. Buffer distance - explicit only
        4. Salt concentration - explicit only

        Returns:
            Dictionary with solvation parameters
        """
        console = self.processor.console

        from rich.panel import Panel

        console.print()
        console.print(Panel(
            "[bold]SOLVATION CONFIGURATION[/bold]\n\n"
            "Configure how your system will be solvated for MD simulation.\n"
            "These settings control the solvent model, simulation box geometry,\n"
            "buffer distance, and ionic strength.",
            title="Solvation Settings",
            border_style="blue",
            expand=False
        ))

        # If the user picked "None" for the water forcefield, that's an implicit
        # solvent signal — skip the redundant solvent-model prompt.
        selected_ffs = self.get_from_workspace("selected_standard_forcefields", None) or {}
        if "water" in selected_ffs and selected_ffs["water"] is None:
            console.print("[blue]No water forcefield selected → using implicit solvent (GB/PB continuum model)[/blue]")
            params = {'solvent_model': 'implicit'}
            self.update_workspace("solvation_parameters", params)
            return params

        # Step 1: Solvent model
        console.print("\n[bold]Solvent Model[/bold]")
        console.print("  1. Explicit solvent (solvation box with water molecules)")
        console.print("  2. Implicit solvent (GB/PB continuum model, no water box)")

        solvent_choice = prompt_with_context(
            processor=self.processor,
            prompt="Select solvent model",
            choices=["1", "2"],
            default="1",
            module="Topology Generator",
            description="Select solvent model"
        )

        if solvent_choice == "2":
            console.print("[blue]Using implicit solvent - no solvation box will be added[/blue]")
            params = {'solvent_model': 'implicit'}
            self.update_workspace("solvation_parameters", params)
            return params

        # Explicit solvent - continue with box configuration
        console.print("[blue]Using explicit solvent[/blue]")

        # Step 2: Box shape
        console.print("\n[bold]Box Shape[/bold]")
        console.print("  1. Truncated octahedron [yellow](Recommended)[/yellow] - ~29% fewer waters than rectangular")
        console.print("  2. Rectangular box")

        box_choice = prompt_with_context(
            processor=self.processor,
            prompt="Select box shape",
            choices=["1", "2"],
            default="1",
            module="Topology Generator",
            description="Select box shape"
        )
        use_octahedron = (box_choice == "1")

        # Step 3: Buffer distance
        buffer, buffer_xyz, oct_diagonal, iso = self._prompt_buffer(use_octahedron)

        # Step 4: Salt composition (one or more salts)
        salts, neutralize_index = self._prompt_salt_composition()

        # Summary
        box_name = "truncated octahedron" if use_octahedron else "rectangular"
        buf_desc = _buffer_describe(buffer, buffer_xyz, iso)
        salts_summary = _salts_summary(salts)
        console.print(f"\n[green]Solvation: {box_name}, {buf_desc} buffer, {salts_summary}[/green]")

        params = {
            'solvent_model': 'explicit',
            'use_octahedron': use_octahedron,
            'buffer': buffer,
            'buffer_xyz': buffer_xyz,
            'oct_diagonal': oct_diagonal,
            'iso': iso,
            'salts': salts,
            'neutralize_index': neutralize_index,
        }
        self.update_workspace("solvation_parameters", params)
        return params

    def _prompt_salt_composition(self) -> tuple:
        """Pick one or more salts (cation/anion/concentration each) and choose
        which one provides the protein-charge neutralizers. Returns
        (salts list, neutralize_index)."""
        console = self.processor.console

        console.print("\n[bold]Salt Composition[/bold]")
        console.print("[grey50]Add one or more salts. Each salt's cation, anion, and concentration "
                      "are independent.[/grey50]")

        salts = []
        while True:
            console.print(f"\n[bold blue]Salt #{len(salts) + 1}[/bold blue]")
            cation, anion = self._prompt_counter_ions()
            conc = self._prompt_concentration(_salt_label(cation, anion))
            salts.append({'cation': cation, 'anion': anion, 'concentration': conc})

            console.print("\n[blue]Configured salts:[/blue]")
            for i, s in enumerate(salts, 1):
                conc_mM = s['concentration'] * 1000
                lbl = _salt_label(s['cation'], s['anion'])
                if conc_mM > 0:
                    console.print(f"  {i}. {conc_mM:.0f} mM {lbl}")
                else:
                    console.print(f"  {i}. neutralize-only {lbl}")

            if not confirm_with_context(
                self.processor,
                "Add another salt?",
                default=False,
                module="Topology Generator",
                description="Add another salt to the mix",
            ):
                break

        # Pick the neutralizer salt when there's more than one
        neutralize_index = 0
        if len(salts) > 1:
            console.print("\n[bold]Charge Neutralization[/bold]")
            console.print("[grey50]Choose which salt's ions are used to balance the protein charge "
                          "via tLEaP's `addions` command.[/grey50]")
            for i, s in enumerate(salts, 1):
                lbl = _salt_label(s['cation'], s['anion'])
                console.print(f"  {i}. {lbl}")
            choice = prompt_with_context(
                processor=self.processor,
                prompt="Use which salt for neutralization",
                choices=[str(i) for i in range(1, len(salts) + 1)],
                default="1",
                module="Topology Generator",
                description="Select neutralizer salt"
            )
            neutralize_index = int(choice) - 1

        return salts, neutralize_index

    def _prompt_concentration(self, salt_label: str) -> float:
        """Prompt for the concentration of one salt. Returns molarity (mol/L)."""
        console = self.processor.console
        console.print(f"\n[bold]{salt_label} concentration[/bold]")
        console.print(f"  1. 150 mM {salt_label} [yellow](Recommended)[/yellow] - physiological ionic strength")
        console.print(f"  2. 100 mM {salt_label}")
        console.print(f"  3. 50 mM {salt_label}")
        console.print(f"  4. Neutralize only (0 mM)")
        console.print(f"  5. Custom concentration")

        choice = prompt_with_context(
            processor=self.processor,
            prompt="Select concentration",
            choices=["1", "2", "3", "4", "5"],
            default="1",
            module="Topology Generator",
            description=f"Select {salt_label} concentration"
        )
        molarity_map = {"1": 0.15, "2": 0.10, "3": 0.05, "4": 0.0}
        if choice != "5":
            return molarity_map[choice]

        while True:
            raw = prompt_with_context(
                processor=self.processor,
                prompt="Enter concentration in mM",
                default="150",
                module="Topology Generator",
                description=f"Enter custom {salt_label} concentration"
            )
            try:
                M = float(raw) / 1000.0
                if 0.0 <= M <= 1.0:
                    return M
                console.print("[yellow]Concentration should be between 0 and 1000 mM[/yellow]")
            except ValueError:
                console.print("[yellow]Please enter a number in mM[/yellow]")

    def _prompt_buffer(self, use_octahedron: bool) -> tuple:
        """Prompt for buffer geometry. Returns (buffer, buffer_xyz, oct_diagonal, iso).

        For rectangular boxes (solvateBox), the user can pick uniform or
        per-axis (separate x/y/z buffers). For truncated octahedral boxes,
        only a uniform buffer is offered: tLEaP's solvateOct rejects the
        list-form buffer documented in its `help` output (verified against
        AmberTools 25.2.0 — it always errors with "iso requires a single
        clearance value"). Iso mode is also not exposed here because
        solvateOct already orients onto principal axes internally.
        """
        console = self.processor.console

        console.print("\n[bold]Buffer Distance[/bold]")
        console.print("[grey50]Distance from solute to box edge. 8-10 A typical, "
                      "12-15 A for large conformational changes.[/grey50]")

        if use_octahedron:
            console.print("[grey50]solvateOct accepts only a single uniform buffer in this "
                          "AmberTools build. Switch to a rectangular box if you need "
                          "per-axis buffers.[/grey50]")
            buf = self._prompt_buffer_value("Buffer distance (Angstroms)", "10.0",
                                             "Enter buffer distance")
            return buf, None, 0.0, False

        console.print("\n  1. Uniform [yellow](Recommended)[/yellow] - same buffer in x, y, z")
        console.print("  2. Per-axis - different x, y, z buffers")

        mode = prompt_with_context(
            processor=self.processor,
            prompt="Select buffer mode",
            choices=["1", "2"],
            default="1",
            module="Topology Generator",
            description="Select buffer mode"
        )

        if mode == "2":
            bx = self._prompt_buffer_value("Buffer in X direction (Angstroms)", "10.0",
                                            "Enter X buffer")
            by = self._prompt_buffer_value("Buffer in Y direction (Angstroms)", "10.0",
                                            "Enter Y buffer")
            bz = self._prompt_buffer_value("Buffer in Z direction (Angstroms)", "10.0",
                                            "Enter Z buffer")
            return None, (bx, by, bz), 0.0, False

        buf = self._prompt_buffer_value("Buffer distance (Angstroms)", "10.0",
                                         "Enter buffer distance")
        return buf, None, 0.0, False

    def _prompt_buffer_value(self, prompt: str, default: str, description: str,
                              min_val: float = 5.0, max_val: float = 25.0) -> float:
        """Prompt for one buffer-axis value, validating against [min_val, max_val]."""
        console = self.processor.console
        while True:
            raw = prompt_with_context(
                processor=self.processor,
                prompt=prompt,
                default=default,
                module="Topology Generator",
                description=description
            )
            try:
                v = float(raw)
                if min_val <= v <= max_val:
                    return v
                console.print(f"[yellow]Value should be between {min_val} and {max_val} Angstroms[/yellow]")
            except ValueError:
                console.print("[yellow]Please enter a number[/yellow]")

    def _available_counter_ions(self) -> tuple:
        """Return (cations, anions, water_name) usable with the current water model.

        Every water leaprc except TIP5P auto-loads parameters for all common
        monovalent + divalent counter-ions. For TIP5P we still return the full
        list and warn at the prompt — the user has to supply ion params manually.
        """
        selected_ffs = self.get_from_workspace("selected_standard_forcefields", None) or {}
        water_sel = selected_ffs.get("water")
        water_name = water_sel.get("name") if isinstance(water_sel, dict) else None
        if not water_name:
            membrane = self.get_from_workspace("membrane_config", {}) or {}
            water_name = membrane.get("effective_water_model") or ""
        return list(_COUNTER_ION_CATIONS), list(_COUNTER_ION_ANIONS), water_name

    def _prompt_counter_ions(self) -> tuple:
        """Interactively pick a cation and anion. Returns two dicts with label/unit/charge."""
        console = self.processor.console
        cations, anions, water_name = self._available_counter_ions()

        if water_name and water_name.upper().replace("-", "") in {"TIP5P"}:
            console.print(
                "[yellow]Warning: TIP5P does not auto-load ion parameters. "
                "You may need to source an ion frcmod manually after solvation.[/yellow]"
            )

        console.print("\n[bold]Cation[/bold]")
        console.print("[grey50]Monovalent (Na+/K+/Li+/Cs+/Rb+) for typical salt; divalent (Mg2+/Ca2+/Zn2+) "
                      "yields a 1:2 salt with the chosen anion.[/grey50]")
        for i, c in enumerate(cations, 1):
            rec = " [yellow](Recommended)[/yellow]" if c["label"] == "Na+" else ""
            console.print(f"  {i}. {c['label']}{rec}")
        cation_choice = prompt_with_context(
            processor=self.processor,
            prompt="Select cation",
            choices=[str(i) for i in range(1, len(cations) + 1)],
            default="1",
            module="Topology Generator",
            description="Select cation"
        )
        cation = cations[int(cation_choice) - 1]

        console.print("\n[bold]Anion[/bold]")
        for i, a in enumerate(anions, 1):
            rec = " [yellow](Recommended)[/yellow]" if a["label"] == "Cl-" else ""
            console.print(f"  {i}. {a['label']}{rec}")
        anion_choice = prompt_with_context(
            processor=self.processor,
            prompt="Select anion",
            choices=[str(i) for i in range(1, len(anions) + 1)],
            default="1",
            module="Topology Generator",
            description="Select anion"
        )
        anion = anions[int(anion_choice) - 1]
        return cation, anion

    def configure_tleap_parameters(self):
        """Configure tLEaP template - show editable tLEaP script"""
        self.processor.console.print()

        # Check if this is a membrane system (already solvated by packmol-memgen)
        is_membrane = self.get_from_workspace("is_membrane_system", False)

        # Check if previous configuration exists
        selected_standard = self.get_from_workspace("selected_standard_forcefields", None)
        solvation_params = self.get_from_workspace("solvation_parameters", None)
        template = self.get_from_workspace("tleap_template", None)

        # MCPB preprocessing writes a PROVISIONAL implicit solvation (and a water
        # FF with box='none') into these shared keys so its own metal-free tleap
        # builds without a box. That is not a user solvation choice: the
        # production topology still needs one. Treat it as unconfigured so the
        # solvation prompt runs, and drop the cached preprocessing template (which
        # has no solvation section) so it is regenerated with the chosen box.
        # Without this, the `if not solvation_params` guard below is silently
        # satisfied and the production topology is built with no periodic box —
        # pmemd then aborts with "peek_ewald_inpcrd: Box info not found".
        #
        # EXCEPTION: while preprocessing runs its OWN metal-free tleap it sets
        # `_preprocessing_tleap_active`; there we must keep the provisional
        # implicit solvation (that build is intentionally box-less) and NOT
        # prompt — the prompt belongs to the user's production run only.
        if self._provisional_solvation_needs_prompt(solvation_params, is_membrane):
            solvation_params = None
            template = None
            self.update_workspace("solvation_parameters", None)
            self.update_workspace("tleap_template", None)

        if is_membrane:
            # Membrane system: auto-populate forcefields and skip solvation
            if not selected_standard:
                selected_standard = self._get_membrane_forcefields()
            else:
                # Membrane builder already stored FF selection — augment with lipid_ext
                selected_standard = self._get_membrane_forcefields()
            solvation_params = {"solvent_model": "membrane_pre_solvated"}
            self.update_workspace("solvation_parameters", solvation_params)
            template = None  # Force template regeneration with membrane settings

        elif selected_standard and solvation_params and template:
            # Summarize what was previously configured
            ff_summary = []
            for key in ["protein", "water", "ions"]:
                val = selected_standard.get(key)
                if val and val != "None":
                    ff_summary.append(f"{key}: {val}")
            # Describe the solvation honestly: an implicit model has no box, so
            # don't dress it up with default box/buffer/salt values (which used to
            # print "rectangular, 10.0 A buffer, ..." for a plain implicit config
            # and mask that no solvation box was set).
            if solvation_params.get("solvent_model") == "implicit":
                solvation_desc = "implicit (no solvation box)"
            else:
                box_shape = "truncated octahedron" if solvation_params.get("use_octahedron") else "rectangular"
                buffer = solvation_params.get("buffer", 10.0)
                buffer_xyz = solvation_params.get("buffer_xyz")
                iso = solvation_params.get("iso", False)
                buf_desc = _buffer_describe(buffer if buffer is not None else 10.0, buffer_xyz, iso)
                salts_list, _ = _normalize_salts(solvation_params)
                solvation_desc = f"{box_shape}, {buf_desc} buffer, {_salts_summary(salts_list)}"
            self.processor.console.print(
                f"[blue]Previous settings:[/blue] {', '.join(ff_summary)}\n"
                f"[blue]Solvation:[/blue] {solvation_desc}"
            )
            reconfigure = confirm_with_context(
                self.processor,
                "Reconfigure all settings (forcefields, solvation, template)?",
                default=False,
                module="Topology Generator",
                description="Reconfigure forcefield and solvation settings"
            )
            if reconfigure:
                selected_standard = None
                solvation_params = None
                template = None
                self.update_workspace("selected_standard_forcefields", None)
                self.update_workspace("solvation_parameters", None)
                self.update_workspace("tleap_template", None)

        # Step 1: Select standard forcefields interactively
        # This must happen BEFORE template generation so we know what to include
        if not selected_standard:
            # Layer 1 — show prereq panel up front (before FF picks) so users
            # know which AMBER leaprcs their cofactor selections will require.
            self._show_cofactor_ff_prerequisites_panel()
            selected_standard = self._select_standard_forcefields_interactive()

        # Step 2: Pick custom forcefields for any non-standard (e.g. redox) sites.
        # This MUST run before _generate_tleap_template so that the info pass
        # emitted inside template generation includes loadoff/loadamberparams
        # for the custom libs — otherwise `charge mol` reports the wrong net
        # charge (non-standard residues get charge 0) and ion counts are off.
        requirements = self._get_single_state_forcefield_requirements()
        if requirements:
            # If the Membrane Builder already ran this picker, its choices are in
            # the workspace. Reuse them by default so the user isn't asked twice,
            # but offer to revise.
            prior_selection = self.get_from_workspace("single_state_selected_forcefields")
            reuse_prior = False
            if prior_selection:
                self.processor.console.print(
                    "[grey50]Using redox-site force-field selections from the "
                    "Membrane Builder.[/grey50]"
                )
                reuse_prior = not confirm_with_context(
                    self.processor,
                    "Revise these force-field selections?",
                    default=False,
                    module="Topology Generator",
                    description="Revise membrane-builder force-field selections",
                )

            if reuse_prior:
                selected_forcefields = prior_selection
            else:
                selected_forcefields = self._select_forcefields_for_single_state(requirements)
                selected_forcefields, ff_atom_types = self._resolve_ff_collisions(
                    selected_forcefields
                )
                self.update_workspace("single_state_ff_requirements", requirements)
                self.update_workspace("single_state_selected_forcefields", selected_forcefields)
                self.update_workspace("ff_resolved_atom_types", ff_atom_types)

        # Step 3: Configure solvation parameters (box shape, buffer, salt)
        if not solvation_params:
            solvation_params = self._configure_solvation_parameters()

        # Step 4: Generate or retrieve template (now uses selected forcefields + solvation params)
        if not template:
            template = self._generate_tleap_template(solvation_params=solvation_params)
            self._save_tleap_template(template)
        else:
            # A cached template freezes its `loadpdb` path at generation time. The
            # priority structure can change AFTER that (most importantly, MCPB
            # Force Field Integration writes the renamed prepared_structure.pdb
            # only once mcpb-4 runs — a template cached before then still points at
            # the metal-free / un-renamed preprocessing PDB). Answering "no" to
            # Reconfigure would then silently build the topology from a stale
            # structure. Re-resolve the priority PDB and rewrite the load path if
            # it drifted, so the reused template always tracks the current best
            # structure.
            template = self._reconcile_template_loadpdb(template)

        # Show current template
        self.processor.console.print("\n[bold]Current tLEaP Template:[/bold]")
        self.processor.console.print("[grey50]" + "─" * 60 + "[/grey50]")

        from rich.syntax import Syntax
        syntax = Syntax(template, "bash", theme="monokai", line_numbers=True)
        self.processor.console.print(syntax)

        self.processor.console.print("[grey50]" + "─" * 60 + "[/grey50]")

        # Let user edit the template
        if confirm_with_context(self.processor, "Edit the tLEaP template?", default=True,
                               module="Topology Generator", description="Edit the tLEaP template"):
            self._edit_tleap_template(template)

        self.processor.console.print(
            "\n[grey50]The template will be used to generate tLEaP input files with your custom settings and ProPrep bond definitions.[/grey50]"
        )
        
        # Initialize basic parameters for file handling
        self.tleap_parameters = self.get_from_workspace("tleap_parameters", {})
        
        # For single state template, we need the transformed PDB file
        selected_pdb = self._select_priority_pdb_file()
        if not selected_pdb:
            return  # Error message already shown in _select_priority_pdb_file
        
        self.tleap_parameters["pdb_file"] = selected_pdb
        
        # Set output prefix based on selected PDB filename
        base_name = os.path.splitext(os.path.basename(selected_pdb))[0]
        output_prefix = base_name  # Don't add _tleap here, it gets added in the filename
        self.tleap_parameters["output_prefix"] = output_prefix
        
        # Store in workspace
        self.update_workspace("tleap_parameters", self.tleap_parameters)


    def write_tleap_input_file_from_template(self):
        """Write tLEaP input file using the template with placeholder substitution"""
        # Get template from workspace
        template = self.get_from_workspace("tleap_template", None)
        if not template:
            self.processor.console.print(
                "[yellow]No template found. Creating default template...[/yellow]"
            )
            template = self._generate_tleap_template()
            self._save_tleap_template(template)
        
        # Get bond definitions from RedoxSite objects
        # Check for actual bonds, not just an initialized dict with empty lists
        has_bonds = self.combined_bonds and any(self.combined_bonds.values())
        if not has_bonds:
            # gather_bond_definitions handles the "no bonds" prompt internally
            # and returns False only if the user explicitly aborted.
            if not self.gather_bond_definitions():
                return
        
        # Get parameters for file info
        self.tleap_parameters = self.get_from_workspace("tleap_parameters", {})
        output_prefix = self.tleap_parameters.get("output_prefix", "system")
        
        # Build atom types section
        atom_types_section = self._build_atom_types_section()
        
        # Build forcefield parameters section
        forcefield_section = self._build_forcefield_parameters_section()
        
        # Build bond definitions section
        bond_section = self._build_bond_definitions_section()
        
        # Substitute placeholders in template
        final_content = template.replace("# ATOM_TYPES_SECTION", atom_types_section)
        final_content = final_content.replace("# FORCEFIELD_PARAMETERS_SECTION", forcefield_section)
        final_content = final_content.replace("# BOND_DEFINITIONS_SECTION", bond_section)
        
        # Get output filename
        default_filename = f"{output_prefix}_tleap.in"
        output_file = prompt_with_context(
            self.processor,
            "Enter output filename for tLEaP input",
            default=default_filename,
            module="Topology Generator",
            description="Output filename for tLEaP input file",
        )

        try:
            with open(output_file, "w") as f:
                f.write(final_content)

            self.output_file = output_file
            self.processor.console.print(
                f"[green]Successfully wrote tLEaP input file to: {output_file}[/green]"
            )

            # Save output file to workspace
            self.update_workspace("tleap_input_file", output_file)

            self.processor.console.print(
                "[green]tLEaP template configuration complete![/green]"
            )

            # Ask if user wants to view the file
            if confirm_with_context(
                self.processor,
                "View the tLEaP input file?",
                default=True,
                module="Topology Generator",
                description="View the generated tLEaP input file",
            ):
                try:
                    with open(output_file, "r") as f:
                        content = f.read()
                    self.processor.console.print(
                        Panel(
                            content,
                            title=f"tLEaP Input File: {output_file}",
                            border_style="green",
                            expand=False
                        )
                    )
                except Exception as e:
                    self.processor.console.print(
                        f"[yellow]Error reading file: {str(e)}[/yellow]"
                    )
            
            return True
        except Exception as e:
            self.processor.console.print(
                f"[red]Error writing tLEaP input file: {str(e)}[/red]"
            )
            return False
    
    def generate_single_state_tleap(self):
        """Generate tLEaP input for a single state - configure and write in one step"""
        from rich.panel import Panel

        self.processor.console.print("\n[bold blue]Generate tLEaP Input for Single State[/bold blue]")
        self.processor.console.print()

        # Educational overview panel
        self.processor.console.print(Panel(
            "[bold blue]What is tLEaP?[/bold blue]\n"
            "tLEaP is AMBER's force-field setup tool. It takes your prepared PDB and\n"
            "writes the topology (prmtop) and coordinate (rst7) files an MD run needs.\n\n"
            "[bold blue]What you'll be asked to choose:[/bold blue]\n"
            "  1. [blue]Forcefields[/blue]  -- protein, water, ion, and any specialty force fields\n"
            "  2. [blue]Solvation[/blue]    -- solvent model, box shape, buffer, and salt concentration\n"
            "                    (ProPrep computes exact ion counts via a tLEaP info pass)\n\n"
            "[bold blue]What ProPrep does for you:[/bold blue]\n"
            "  3. [blue]Template[/blue]     -- writes an editable tLEaP input file from your choices\n"
            "  4. [blue]Bond defs[/blue]    -- fills in bond commands from your redox-site definitions\n"
            "  5. [blue]Parameters[/blue]   -- fills in custom force-field parameters for non-standard residues\n\n"
            "[grey50]The result: a complete, editable tLEaP input file, ready to build MD topology.[/grey50]",
            title="tLEaP Input Generation",
            border_style="blue",
            expand=False
        ))

        # ═══════════════════════════════════════════════════════════════════
        # STEP 1: Forcefields, Solvation & Template
        # ═══════════════════════════════════════════════════════════════════
        self.processor.console.print("\n" + "═" * 70)
        self.processor.console.print("[bold blue]STEP 1: Forcefields, Solvation & Template[/bold blue]")
        self.processor.console.print("═" * 70 + "\n")

        self.configure_tleap_parameters()

        # ═══════════════════════════════════════════════════════════════════
        # STEP 2: Structure Preparation (PDB Reordering & RedoxSite Sync)
        # ═══════════════════════════════════════════════════════════════════
        self.processor.console.print("\n" + "═" * 70)
        self.processor.console.print("[bold blue]STEP 2: Structure Preparation[/bold blue]")
        self.processor.console.print("═" * 70 + "\n")

        # Perform PDB reordering and RedoxSite synchronization BEFORE bond generation
        if not self._prepare_structure_for_tleap():
            self.processor.console.print("[yellow]Structure preparation skipped or failed[/yellow]")
            # Continue anyway - user may not need reordering

        # ═══════════════════════════════════════════════════════════════════
        # STEP 3: tLEaP Input File Generation
        # ═══════════════════════════════════════════════════════════════════
        self.processor.console.print("\n" + "═" * 70)
        self.processor.console.print("[bold blue]STEP 3: tLEaP Input File Generation[/bold blue]")
        self.processor.console.print("═" * 70 + "\n")

        return self.write_tleap_input_file_from_template()
    
    def _prepare_structure_for_tleap(self) -> bool:
        """
        Prepare PDB structure for tLEaP input generation.

        This performs PDB reordering and RedoxSite synchronization BEFORE
        bond definitions are gathered, ensuring bond commands reference
        correct residue numbers.

        Returns:
            bool: True if preparation succeeded or wasn't needed, False if cancelled
        """
        from proprep.tleap_prep.pdb_molecule_configurator import PDBMoleculeConfigurator
        from pathlib import Path

        console = self.processor.console

        # Membrane systems: packmol-memgen output is already a single packed PDB
        # with its own atom ordering. Re-running PDB reordering / TER-record
        # fixing on a packed bilayer would treat lipids as separate "molecules"
        # and create noise. Reordering belongs upstream of membrane building.
        #
        # REVISIT FOR PROTEIN-MEMBRANE SYSTEMS [2026-05-13]:
        # This skip is currently scoped to ALL membrane systems, including
        # protein-membrane builds. When we add the first redox-protein-in-bilayer
        # workflow, decide whether to:
        #   (a) keep skipping (treat packmol-memgen output as authoritative —
        #       any reordering/TER fixing must happen on the protein PDB BEFORE
        #       the membrane is built), or
        #   (b) narrow the skip to empty bilayers only
        #       (`is_membrane_system AND not membrane_config.protein_pdb`) and
        #       let STEP 2 run TER fixing on the protein portion of the
        #       packed PDB — would need the TER-fixer to be lipid-aware so it
        #       doesn't insert TERs between every lipid copy.
        # Option (a) is the safer default; option (b) is more user-friendly
        # if reordering wasn't done upstream.
        if self.get_from_workspace("is_membrane_system", False):
            console.print(
                "[grey50]Membrane system — packmol-memgen output used as-is "
                "(structure preparation not applicable).[/grey50]"
            )
            return True

        # Get PDB file using priority logic (same as template generation)
        # Also track which workspace key it came from so we can update the right one
        workspace = self.get_workspace()
        pdb_file = None
        pdb_workspace_key = None

        priority_keys = [
            ("reordered_pdb_file", "reordered"),
            ("protonation_pdb_file", "protonation-updated"),
            ("transformed_pdb_file", "transformed"),
            ("repaired_pdb_file", "repaired"),
            ("filtered_pdb_file", "filtered"),
            ("local_pdb_file", "user-loaded local PDB"),
            ("rcsb_pdb_file", "user-loaded RCSB PDB"),
        ]

        for workspace_key, description in priority_keys:
            candidate_pdb = self.get_from_workspace(workspace_key)
            if candidate_pdb and os.path.exists(candidate_pdb):
                pdb_file = candidate_pdb
                pdb_workspace_key = workspace_key
                break

        if not pdb_file:
            console.print("[grey50]No PDB file found in workspace - skipping structure preparation[/grey50]")
            return True

        pdb_path = Path(pdb_file)
        config_file = pdb_path.parent / f"{pdb_path.stem}_molecule_config.json"
        reordered_pdb = pdb_path.parent / f"{pdb_path.stem}_reordered.pdb"

        # Initialize configurator
        configurator = PDBMoleculeConfigurator(console=console, processor=self.processor)

        # Analyze PDB structure
        console.print(f"[blue]Analyzing PDB structure: {pdb_path.name}[/blue]\n")
        segments = configurator.analyze_pdb_segments(pdb_file)

        if not segments:
            console.print("[yellow]No chains found in PDB file[/yellow]")
            return True

        # Check if structure is complex (multiple chains)
        if len(segments) <= 1:
            # Display structure summary for simple structures
            self._display_structure_summary(segments, console)
            console.print("\n[green]✓ Simple structure - no molecule reordering needed[/green]")

            # Validate and auto-fix TER records (always needed for tleap)
            needs_fixing, issues = self._check_ter_records(pdb_file)
            if needs_fixing:
                console.print(f"\n[blue]Checking TER records in {Path(pdb_file).name}...[/blue]")
                console.print(f"[yellow]Found {len(issues)} TER record issue(s):[/yellow]")
                for issue in issues:
                    console.print(f"  • {issue}")
                console.print(f"\n[grey50]TER records separate molecules in PDB files. TLeaP requires them to "
                              f"correctly identify molecular boundaries (protein chain termini, "
                              f"ligands, cofactors, ions). Missing TER records cause TLeaP to treat "
                              f"separate molecules as one continuous chain, leading to incorrect bonds.[/grey50]")
                if self._fix_ter_records(pdb_file, issues):
                    console.print(f"[green]✓ Inserted {len(issues)} TER record(s) in {Path(pdb_file).name}[/green]")
                else:
                    console.print(f"[yellow]⚠ Could not fix TER records in {Path(pdb_file).name}[/yellow]")
            else:
                console.print(f"\n[green]✓ TER records valid in {Path(pdb_file).name}[/green]")

            return True

        # Complex structure detected
        console.print("[yellow]⚠ Multi-chain structure detected[/yellow]\n")

        # Ask whether to proceed with reordering
        if not confirm_with_context(
            self.processor,
            "Reorder PDB chains/molecules for tLEaP?",
            default=True,
            module="Topology Generator",
            description="Structure reordering",
        ):
            console.print("[grey50]Skipping structure preparation[/grey50]\n")
            self.update_workspace("reordering_skipped", True)
            return True

        # If reordered PDB already exists, offer use/redo/skip
        if reordered_pdb.exists():
            console.print(f"[blue]Found existing reordered PDB: {reordered_pdb.name}[/blue]")
            action = prompt_with_context(
                self.processor,
                "Use existing reordered PDB, redo reordering, or skip?",
                choices=["use", "redo", "skip"],
                default="use",
                module="Topology Generator",
                description="Reordered PDB action",
            )
            if action == "use":
                self.update_workspace("reordered_pdb_file", str(reordered_pdb))
                console.print(f"[green]✓ Using {reordered_pdb.name} for tLEaP input[/green]\n")
                return True
            elif action == "skip":
                console.print(f"[grey50]Skipping structure preparation (using {pdb_path.name} as-is)[/grey50]\n")
                return True

        # Check for existing configuration
        config = None
        if config_file.exists():
            console.print(f"[blue]Found existing configuration: {config_file.name}[/blue]")
            if confirm_with_context(
                self.processor,
                "Load existing configuration?",
                default=True,
                module="Topology Generator",
                description="Load existing tLEaP configuration file",
            ):
                config = configurator.load_configuration(str(config_file))

        # If no existing config, run interactive configuration
        # The configurator will handle all prompts and displays
        if config is None:
            config = configurator.configure_molecule_grouping(pdb_file)

            if config is None:
                console.print("[yellow]Configuration cancelled[/yellow]")
                return False

        # Reorder PDB based on configuration (with renumbering)
        console.print()
        reordered_file, mapper = configurator.reorder_pdb(
            pdb_file,
            config,
            output_file=str(reordered_pdb)
        )

        if not reordered_file or not mapper:
            console.print("[red]✗ Failed to reorder PDB[/red]")
            return False

        # Validate and fix TER records in reordered PDB
        console.print()
        console.print("[blue]Validating TER records in reordered PDB...[/blue]")
        if not self._validate_and_fix_single_pdb_ter_records(reordered_file):
            console.print("[yellow]⚠ TER record validation had issues (continuing anyway)[/yellow]")

        # Synchronize RedoxSites with reordered/renumbered structure
        console.print()
        redox_sites = workspace.get("detected_redox_sites", [])
        if redox_sites:
            sync_summary = configurator.synchronize_redox_sites(redox_sites, mapper)
            # Save updated redox sites back to workspace
            workspace.set("detected_redox_sites", redox_sites)
            console.print(
                f"[green]✓ RedoxSites synchronized with reordered structure[/green]"
            )
        else:
            console.print("[grey50]No RedoxSites to synchronize[/grey50]")

        # Note: Protonation results synchronization is handled in CPIN generation
        # using coordinate-based mapping from structure_with_prot_resnames to the
        # PDB corresponding to the topology (prmtop/rst7)

        # Store reordered PDB under its own dedicated workspace key
        self.update_workspace("reordered_pdb_file", reordered_file)

        # Update template to reference reordered PDB
        template = self.get_from_workspace("tleap_template", None)
        if template:
            old_basename = os.path.basename(pdb_file)
            new_basename = os.path.basename(reordered_file)
            updated_template = template.replace(old_basename, new_basename)
            self.update_workspace("tleap_template", updated_template)
            console.print(f"[green]✓ Updated template to use {new_basename}[/green]")

        console.print(f"[green]✓ Structure preparation complete[/green]\n")

        return True

    def write_tleap_input_file(self):
        """Write tLEaP input file - uses template-based approach"""
        self.processor.console.print(
            "[yellow]Using template-based tLEaP input generation...[/yellow]"
        )
        return self.write_tleap_input_file_from_template()

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Core module methods

    def get_workspace_requirements(self) -> List[str]:
        """Get workspace requirements"""
        # Requires any PDB structure to be available
        return []  # No hard requirements - will use priority selection

    def get_workspace_outputs(self) -> List[str]:
        """Get workspace outputs"""
        return [
            "combined_tleap_commands",
            "solvation_parameters",
            "single_state_ff_requirements",
            "single_state_selected_forcefields",
            "tleap_parameters",
            "tleap_input_file",
            "tleap_template",
            "microstate_tleap_template",
            "generated_microstate_tleap_files",
            "parm7_file",
            "rst7_file",
            "selected_standard_forcefields",
            "_active_tleap_input_file",
            "cpin_config",
            "cpin_file",
        ]

    def availability_note(self, workspace):
        """Menu note when unavailable (○). Mirrors can_process exactly."""
        return None if self.can_process(workspace) else \
            "Needs a structure, a packed membrane, or a parm7/rst7 pair"

    def can_process(self, workspace: Dict[str, Any]) -> bool:
        """Check if the module can process the current workspace.

        Available if a structure (PDB) is present — needed to generate a
        topology from scratch — OR a prmtop/rst7 pair is loaded. The latter is
        what the PB Titrate sub-workflow operates on, so a run resumed from a
        reloaded prmtop/rst7 with no PDB must still reach this module.
        """
        from proprep.utils.structure_selector import StructureSelector

        # Membrane systems use membrane_packed_pdb (not in the StructureSelector
        # registry — that registry is for protein-stage modules).
        membrane_pdb = workspace.get("membrane_packed_pdb")
        if membrane_pdb and os.path.exists(membrane_pdb):
            return True

        selector = StructureSelector(workspace, self.processor.console if self.processor else None)
        status = selector.get_structure_status()
        if status.get("has_any", False):
            return True

        # Fall back to a loaded topology/coordinate pair (parm7_file + rst7_file)
        # — the input PB Titrate and topology-loading paths consume.
        return bool(workspace.get("parm7_file") and workspace.get("rst7_file"))

    def process(self, workspace):
        """Process the workspace"""
        if self.can_process(workspace):
            # Check for RedoxSite objects
            redox_sites = workspace.get("detected_redox_sites")

            if redox_sites:
                # Create combined bonds dictionary - matches comprehensive_redox_detector categories
                combined_bonds = {
                    "covalent": [],
                    "coordinate": [],
                    "metal-metal": [],
                    "disulfide": [],
                    "peptide_backbone": [],
                    "other": [],
                }

                # Convert RedoxSite objects to bond commands
                redox_bond_commands = self._convert_redox_sites_to_tleap_commands(redox_sites)
                
                # Merge into combined bonds
                for bond_type, bonds in redox_bond_commands.items():
                    if bond_type in combined_bonds:
                        combined_bonds[bond_type].extend(bonds)

                # Store combined bonds in workspace
                workspace.set("combined_tleap_commands", combined_bonds)

            # Add basic tLEaP parameters if not present
            if not workspace.get("tleap_parameters"):
                # Set defaults
                pdb_file = workspace.get("pdb_file", "")
                output_prefix = "system"
                if pdb_file:
                    output_prefix = os.path.splitext(os.path.basename(pdb_file))[0]

                tleap_parameters = {
                    "forcefield": "leaprc.constph",
                    "water_model": "leaprc.water.tip3p",
                    "ions_model": "leaprc.water.tip3p",
                    "gaff": "leaprc.gaff2",
                    "box_type": "oct",
                    "box_distance": 10.0,
                    "neutralize": True,
                    "salt_conc": 0.15,
                    "positive_ion": "Na+",
                    "negative_ion": "Cl-",
                    "pdb_file": pdb_file,
                    "unit_name": "mol",
                    "complex_name": "complex",
                    "output_prefix": output_prefix,
                }

                workspace.set("tleap_parameters", tleap_parameters)

        return workspace
    
    def generate_microstate_inputs_template_based(self, metadata_file: str = None) -> bool:
        """
        Generate tLEaP input files for all redox microstates using template-based approach.

        Args:
            metadata_file: Path to microstate metadata JSON file. If None, will try to find it.

        Returns:
            True if successful, False otherwise
        """
        try:
            import json
            import os
            from pathlib import Path
            from proprep.forcefield_params import discover_forcefield_files, ForcefieldNotFoundError

            console = self.processor.console

            console.print("\n[bold blue]Generate tLEaP Inputs for All Redox Microstates[/bold blue]")
            console.print()

            # ═══════════════════════════════════════════════════════════════════
            # STEP 1: Template Configuration
            # ═══════════════════════════════════════════════════════════════════
            console.print("═" * 70)
            console.print("[bold blue]STEP 1: Template Configuration[/bold blue]")
            console.print("═" * 70 + "\n")

            console.print(
                "[yellow]This will generate tLEaP inputs for all microstates using an editable template.[/yellow]"
            )
            console.print(
                "[grey50]ProPrep will automatically fill in PDB files, atom types, forcefield files, and bond definitions.[/grey50]"
            )
            console.print()

            # Select standard forcefields interactively (must happen BEFORE template generation)
            selected_standard = self.get_from_workspace("selected_standard_forcefields", None)
            if not selected_standard:
                selected_standard = self._select_standard_forcefields_interactive()

            # Get or create microstate template (now uses selected forcefields)
            template = self.get_from_workspace("microstate_tleap_template", None)
            if not template:
                template = self._generate_microstate_tleap_template()
                self.update_workspace("microstate_tleap_template", template)

            # Show current template
            console.print("[bold]Current Microstate tLEaP Template:[/bold]")
            console.print("[grey50]" + "─" * 60 + "[/grey50]")

            from rich.syntax import Syntax
            syntax = Syntax(template, "bash", theme="monokai", line_numbers=True)
            console.print(syntax)

            console.print("[grey50]" + "─" * 60 + "[/grey50]")

            # Let user edit the template
            if confirm_with_context(
                processor=self.processor,
                prompt="Edit the microstate tLEaP template?",
                default=False,
                module="Topology Generator",
                description="Edit microstate template"
            ):
                self._edit_microstate_template(template)
                # Get updated template
                template = self.get_from_workspace("microstate_tleap_template", template)

            # ═══════════════════════════════════════════════════════════════════
            # STEP 2: Structure Preparation (applies to all microstates)
            # ═══════════════════════════════════════════════════════════════════
            console.print("\n" + "═" * 70)
            console.print("[bold blue]STEP 2: Structure Preparation[/bold blue]")
            console.print("[grey50](All microstates have identical structure - preparation done once)[/grey50]")
            console.print("═" * 70 + "\n")

            # Prepare microstate structures (reorder all if needed, sync RedoxSites once)
            if not self._prepare_microstate_structures(metadata_file):
                console.print("[yellow]Structure preparation failed or was cancelled[/yellow]")
                return False

            # ═══════════════════════════════════════════════════════════════════
            # STEP 3: Bond Definition Generation
            # ═══════════════════════════════════════════════════════════════════
            console.print("\n" + "═" * 70)
            console.print("[bold blue]STEP 3: Bond Definition Generation[/bold blue]")
            console.print("[grey50](Same bonds apply to all microstates - generated once)[/grey50]")
            console.print("═" * 70 + "\n")

            # Gather bond definitions if not already present
            combined_bonds = self.get_from_workspace("combined_tleap_commands", {})
            if not combined_bonds or not any(len(bonds) > 0 for bonds in combined_bonds.values()):
                console.print("[blue]Gathering bond definitions from RedoxSite objects...[/blue]")
                self.gather_bond_definitions()

            # ═══════════════════════════════════════════════════════════════════
            # STEP 4: tLEaP Input File Generation
            # ═══════════════════════════════════════════════════════════════════
            console.print("\n" + "═" * 70)
            console.print("[bold blue]STEP 4: tLEaP Input File Generation[/bold blue]")
            console.print("═" * 70 + "\n")

            # Continue with the actual processing
            return self._process_microstate_inputs_with_template(metadata_file, template)
            
        except Exception as e:
            logger.error(f"Error generating microstate tLEaP inputs: {e}")
            console.print(f"[red]Error: {str(e)}[/red]")
            return False
    
    def _prepare_microstate_structures(self, metadata_file: str = None) -> bool:
        """
        Prepare all microstate PDB structures for tLEaP input generation.

        Since all microstates have identical structure (only residue names differ),
        we analyze one microstate and apply the same reordering to all.

        Args:
            metadata_file: Path to microstate metadata JSON file

        Returns:
            bool: True if preparation succeeded or wasn't needed, False if cancelled
        """
        from proprep.tleap_prep.pdb_molecule_configurator import PDBMoleculeConfigurator
        from pathlib import Path
        import json

        console = self.processor.console

        # Find metadata file if not provided
        if not metadata_file:
            metadata_files = list(Path('.').glob('*microstate*metadata*.json'))
            if not metadata_files:
                console.print("[yellow]No microstate metadata file found - skipping structure preparation[/yellow]")
                return True
            elif len(metadata_files) == 1:
                metadata_file = str(metadata_files[0])
            else:
                console.print("Multiple metadata files found:")
                for i, file in enumerate(metadata_files, 1):
                    console.print(f"{i}. {file}")
                from proprep.utils.prompts import prompt_with_context
                choice = prompt_with_context(self.processor,
                    "Select metadata file",
                    choices=[str(i) for i in range(1, len(metadata_files) + 1)],
                    default="1"
                )
                choice = remap_recorded_index(self.processor, metadata_files, str(choice))
                metadata_file = str(metadata_files[int(choice) - 1])
                annotate_selected_path(self.processor, metadata_files[int(choice) - 1])

        # Load metadata
        if not os.path.exists(metadata_file):
            console.print(f"[yellow]Metadata file not found: {metadata_file}[/yellow]")
            return True

        with open(metadata_file, 'r') as f:
            metadata = json.load(f)

        microstates = metadata.get('microstates', [])
        if not microstates:
            console.print("[yellow]No microstates found in metadata[/yellow]")
            return True

        # Get first microstate PDB (all have identical structure)
        first_microstate = microstates[0]
        sample_pdb = first_microstate.get('filename')

        if not sample_pdb or not os.path.exists(sample_pdb):
            console.print(f"[yellow]Sample PDB not found: {sample_pdb}[/yellow]")
            return True

        pdb_path = Path(sample_pdb)
        config_file = pdb_path.parent / f"{pdb_path.stem.rsplit('_', 1)[0]}_molecule_config.json"

        console.print(f"[blue]Analyzing microstate structures (using {pdb_path.name} as representative)[/blue]\n")

        # Check if already reordered
        if sample_pdb.endswith("_reordered.pdb"):
            console.print("[green]✓ Microstate PDBs already reordered[/green]\n")
            self.update_workspace("reordered_pdb_file", sample_pdb)
            return True

        # Check for existing reordered files
        sample_reordered = pdb_path.parent / f"{pdb_path.stem}_reordered.pdb"
        if sample_reordered.exists():
            console.print(f"[blue]Found existing reordered PDB: {sample_reordered.name}[/blue]")
            action = prompt_with_context(
                processor=self.processor,
                prompt="Use existing reordered PDBs, redo reordering, or skip?",
                choices=["use", "redo", "skip"],
                default="use",
                module="Topology Generator",
                description="Reordered PDB action (microstates)",
            )
            if action == "use":
                # Update all microstate filenames in metadata
                for microstate_info in microstates:
                    orig_path = Path(microstate_info['filename'])
                    reordered_path = orig_path.parent / f"{orig_path.stem}_reordered.pdb"
                    if reordered_path.exists():
                        microstate_info['filename'] = str(reordered_path)
                        reordered_json = reordered_path.parent / f"{reordered_path.stem}_redox_sites_updated.json"
                        if reordered_json.exists():
                            microstate_info['redox_sites_json'] = reordered_json.name

                # Save updated metadata
                with open(metadata_file, 'w') as f:
                    json.dump(metadata, indent=2, fp=f)

                console.print(f"[green]✓ Updated metadata to use reordered PDBs[/green]\n")
                self.update_workspace("reordered_pdb_file", str(sample_reordered))
                return True
            elif action == "skip":
                console.print(f"[grey50]Skipping reordering (using microstates as-is)[/grey50]\n")
                self.update_workspace("transformed_pdb_file", sample_pdb)
                return True

        # Initialize configurator
        configurator = PDBMoleculeConfigurator(console=console, processor=self.processor)

        # Analyze PDB structure
        segments = configurator.analyze_pdb_segments(sample_pdb)

        if not segments:
            console.print("[yellow]No chains found in PDB file[/yellow]")
            return True

        # Check if structure is complex (multiple chains)
        if len(segments) <= 1:
            console.print("\n[green]✓ Simple structure - no molecule reordering needed[/green]\n")
            return True

        # Complex structure detected
        console.print("\n[yellow]⚠ Multi-chain structure detected[/yellow]\n")

        # Ask whether to proceed with reordering
        if not confirm_with_context(
            self.processor,
            "Reorder PDB chains/molecules for tLEaP?",
            default=True,
            module="Topology Generator",
            description="Structure reordering (microstates)",
        ):
            console.print("[grey50]Skipping structure preparation[/grey50]\n")
            self.update_workspace("reordering_skipped", True)
            # Set transformed_pdb_file to the representative microstate PDB
            # so bond generation can find correct post-transformation residue IDs
            self.update_workspace("transformed_pdb_file", sample_pdb)
            return True

        console.print("[grey50]Note: Configuration will apply to ALL microstates (they have identical structure - only oxidation states differ)[/grey50]\n")

        # Check for existing configuration
        config = None
        if config_file.exists():
            console.print(f"[blue]Found existing configuration: {config_file.name}[/blue]")
            if confirm_with_context(
                self.processor,
                "Load existing configuration?",
                default=True,
                module="Topology Generator",
                description="Load existing tLEaP configuration file",
            ):
                config = configurator.load_configuration(str(config_file))

        # If no existing config, run interactive configuration
        # The configurator will handle all prompts and displays
        if config is None:
            config = configurator.configure_molecule_grouping(sample_pdb)

            if config is None:
                console.print("[yellow]Configuration cancelled[/yellow]")
                return False

        # Reorder all microstate PDBs with the same configuration
        console.print(f"\n[blue]Reordering {len(microstates)} microstate PDBs...[/blue]\n")

        mapper = None
        reordered_count = 0
        reordered_files = []

        for i, microstate_info in enumerate(microstates, 1):
            pdb_file = microstate_info['filename']
            if not os.path.exists(pdb_file):
                console.print(f"[yellow]  {i}. Skipping {Path(pdb_file).name} (not found)[/yellow]")
                continue

            pdb_path = Path(pdb_file)
            reordered_pdb = pdb_path.parent / f"{pdb_path.stem}_reordered.pdb"

            console.print(f"  {i}/{len(microstates)}: {pdb_path.name} → {reordered_pdb.name}")

            # Reorder this microstate PDB
            reordered_file, current_mapper = configurator.reorder_pdb(
                pdb_file,
                config,
                output_file=str(reordered_pdb)
            )

            if reordered_file:
                microstate_info['filename'] = reordered_file
                reordered_files.append(reordered_file)
                reordered_count += 1

                # Save the mapper from the first successful reordering
                if mapper is None:
                    mapper = current_mapper
            else:
                console.print(f"[red]    ✗ Failed to reorder {pdb_path.name}[/red]")

        if reordered_count == 0:
            console.print("[red]✗ Failed to reorder any microstate PDBs[/red]")
            return False

        console.print(f"\n[green]✓ Reordered {reordered_count}/{len(microstates)} microstate PDBs[/green]")

        # Validate and fix TER records in all reordered PDBs
        console.print(f"\n[blue]Validating TER records in {len(reordered_files)} reordered PDBs...[/blue]")

        # First pass: check all files for TER issues
        files_with_issues = []
        for reordered_file in reordered_files:
            needs_fixing, issues = self._check_ter_records(reordered_file)
            if needs_fixing:
                files_with_issues.append((reordered_file, issues))

        if files_with_issues:
            console.print(f"\n[yellow]TER record issues found in {len(files_with_issues)}/{len(reordered_files)} files[/yellow]")
            # Show all issues from first file (issues are identical across all microstates)
            first_file, first_issues = files_with_issues[0]
            console.print(f"[grey50]Issues (identical across all microstates):[/grey50]")
            for issue in first_issues:
                console.print(f"  • {issue}")
            console.print(f"\n[grey50]TER records separate molecules in PDB files. TLeaP requires them to "
                          f"correctly identify molecular boundaries. Inserting missing TER records.[/grey50]")

            fixed_count = 0
            failed_count = 0
            for pdb_file, issues in files_with_issues:
                if self._fix_ter_records(pdb_file, issues):
                    fixed_count += 1
                else:
                    failed_count += 1
                    console.print(f"[red]✗ Failed to fix: {Path(pdb_file).name}[/red]")

            if failed_count == 0:
                console.print(f"[green]✓ Inserted TER records in all {fixed_count} files[/green]")
            else:
                console.print(f"[yellow]⚠ Fixed {fixed_count} files, {failed_count} failed[/yellow]")
        else:
            console.print(f"[green]✓ TER records valid in all {len(reordered_files)} reordered PDBs[/green]")

        # Synchronize RedoxSites once (same mapping applies to all microstates)
        console.print()
        workspace = self.processor.workspace
        redox_sites = workspace.get("detected_redox_sites", [])
        if redox_sites and mapper:
            sync_summary = configurator.synchronize_redox_sites(redox_sites, mapper)
            workspace.set("detected_redox_sites", redox_sites)
            console.print(f"[green]✓ RedoxSites synchronized (applies to all microstates)[/green]")

            # Re-export per-microstate redox JSONs so they match the reordered PDBs
            from proprep.redoxsite_prep.transformation.redox_transformation_manager import export_redox_sites_to_json
            exported = 0
            # Embed transformer assignments so each per-microstate JSON can be
            # read back after a restart to recover the site→forcefield mapping.
            ms_transformer_info = self.get_from_workspace("transformer_info", None)
            for microstate_info in metadata['microstates']:
                pdb_path = microstate_info.get('filename')
                if not pdb_path or not os.path.exists(pdb_path):
                    continue
                out_path = export_redox_sites_to_json(
                    redox_sites, pdb_path,
                    transformer_info=ms_transformer_info)
                microstate_info['redox_sites_json'] = out_path.name
                exported += 1
            if exported:
                console.print(f"[green]✓ Re-exported {exported} redox-site JSON file(s) for reordered structures[/green]")
        else:
            console.print("[grey50]No RedoxSites to synchronize[/grey50]")

        # Note: Protonation results synchronization is handled in CPIN generation
        # using coordinate-based mapping (same approach as single-state)

        # Save updated metadata with reordered filenames
        with open(metadata_file, 'w') as f:
            json.dump(metadata, indent=2, fp=f)

        console.print(f"[green]✓ Updated metadata with reordered PDB filenames[/green]")

        # Set reordered_pdb_file to the first reordered microstate PDB
        # so bond generation uses correct post-reorder residue IDs
        if reordered_files:
            self.update_workspace("reordered_pdb_file", reordered_files[0])

        console.print(f"[green]✓ Structure preparation complete[/green]\n")

        return True

    def _edit_microstate_template(self, current_template):
        """Allow user to edit the microstate tLEaP template"""
        import tempfile
        import subprocess
        import os
        
        # Create temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tleap', delete=False) as tmp:
            tmp.write(current_template)
            tmp_path = tmp.name
        
        try:
            # Get user's preferred editor
            editor = os.environ.get('EDITOR', 'nano')
            
            self.processor.console.print(f"[yellow]Opening microstate template in {editor}...[/yellow]")
            self.processor.console.print("[grey50]Save and exit the editor when done.[/grey50]")
            
            # Open editor
            subprocess.run([editor, tmp_path], check=True, stdin=sys.stdin, stdout=sys.stdout, stderr=sys.stderr)
            
            # Read back the edited content
            with open(tmp_path, 'r') as f:
                edited_template = f.read()
            
            # Save the edited template
            self.update_workspace("microstate_tleap_template", edited_template)
            
            self.processor.console.print("[green]Microstate template updated successfully![/green]")
            
            # Show the updated template
            if confirm_with_context(
                processor=self.processor,
                prompt="Show updated microstate template?",
                default=True,
                module="Topology Generator",
                description="Show updated template"
            ):
                self.processor.console.print(f"\n[bold blue]Updated Microstate tLEaP Template:[/bold blue]")
                self.processor.console.print("[grey50]" + "─" * 60 + "[/grey50]")
                
                from rich.syntax import Syntax
                syntax = Syntax(edited_template, "bash", theme="monokai", line_numbers=True)
                self.processor.console.print(syntax)
                
                self.processor.console.print("[grey50]" + "─" * 60 + "[/grey50]")
            
        except subprocess.CalledProcessError:
            self.processor.console.print("[red]Editor was cancelled or failed[/red]")
        
        except FileNotFoundError:
            self.processor.console.print(f"[red]Editor '{editor}' not found. Template not modified.[/red]")
        
        finally:
            # Clean up temp file
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
    
    def _process_microstate_inputs_with_template(self, metadata_file, template):
        """Process microstate inputs using the template approach"""
        import json
        import os
        from pathlib import Path
        
        console = self.processor.console
        
        # Find metadata file if not provided
        if not metadata_file:
            # Look for metadata files in current directory
            metadata_files = list(Path('.').glob('*microstate*metadata*.json'))
            
            if not metadata_files:
                console.print("[red]No microstate metadata file found. Please generate microstates first.[/red]")
                return False
            elif len(metadata_files) == 1:
                metadata_file = str(metadata_files[0])
                console.print(f"[blue]Found metadata file: {metadata_file}[/blue]")
            else:
                # Multiple files found, let user choose
                console.print("Multiple metadata files found:")
                for i, file in enumerate(metadata_files, 1):
                    console.print(f"{i}. {file}")
                
                choice = prompt_with_context(self.processor,
                    "Select metadata file",
                    choices=[str(i) for i in range(1, len(metadata_files) + 1)],
                    default="1"
                )
                choice = remap_recorded_index(self.processor, metadata_files, str(choice))
                metadata_file = str(metadata_files[int(choice) - 1])
                annotate_selected_path(self.processor, metadata_files[int(choice) - 1])
        
        # Load metadata
        if not os.path.exists(metadata_file):
            console.print(f"[red]Metadata file not found: {metadata_file}[/red]")
            return False
        
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)
        
        microstates = metadata.get('microstates', [])
        if not microstates:
            console.print("[red]No microstates found in metadata file.[/red]")
            return False
        
        console.print(f"[blue]Found {len(microstates)} microstates to process[/blue]")

        # Collect unique forcefield requirements
        unique_combinations = self._collect_unique_forcefield_requirements(microstates)

        # Select forcefield files for each unique combination
        selected_forcefields = self._select_forcefields_for_microstates(unique_combinations)
        if not selected_forcefields:
            return False

        # Configure solvation parameters (solvent model, box shape, buffer, salt)
        solvation_params = self.get_from_workspace("solvation_parameters", None)
        if not solvation_params:
            solvation_params = self._configure_solvation_parameters()

        is_explicit = solvation_params.get('solvent_model') == 'explicit'

        # Extract parameters for explicit solvent
        import re
        water_box = "TIP3PBOX"
        # Extract water box from template
        water_match = re.search(r'solvate(?:oct|Box)\s+mol\s+(\w+)', template)
        if water_match:
            water_box = water_match.group(1)

        buffer = solvation_params.get('buffer', 10.0)
        buffer_xyz = solvation_params.get('buffer_xyz')
        oct_diagonal = solvation_params.get('oct_diagonal', 0.0)
        iso = solvation_params.get('iso', False)
        use_octahedron = solvation_params.get('use_octahedron', True)
        salts, neutralize_index = _normalize_salts(solvation_params)

        # Generate tLEaP input for each microstate using template
        generated_inputs = []
        failed_inputs = []
        microstate_ion_counts = {}

        if is_explicit:
            console.print("\n[bold yellow]Running tLEaP info passes for accurate ion counts...[/bold yellow]")
            console.print("[grey50]This runs tLEaP once per microstate with early termination[/grey50]\n")

            from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, MofNCompleteColumn

            shown_educational_panel = False

            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                MofNCompleteColumn(),
                console=console
            ) as progress:
                task = progress.add_task("Info pass...", total=len(microstates))

                for ms_idx, microstate_info in enumerate(microstates, 1):
                    code = microstate_info['code']
                    progress.update(task, description=f"[blue]Microstate {ms_idx:03d}: Info pass[/blue]")

                    system_info = self._run_info_pass_for_microstate(
                        microstate_info, selected_forcefields,
                        water_box=water_box, buffer=buffer, use_octahedron=use_octahedron,
                        buffer_xyz=buffer_xyz, oct_diagonal=oct_diagonal, iso=iso
                    )

                    n_waters = system_info.get('n_waters')
                    net_charge = system_info.get('net_charge')

                    if n_waters is not None and net_charge is not None:
                        ion_counts = self._calculate_multi_salt_ions(
                            n_waters, int(round(net_charge)),
                            salts=salts, neutralize_index=neutralize_index
                        )
                        microstate_ion_counts[code] = ion_counts

                        # Show educational panel once (for the first microstate)
                        if not shown_educational_panel:
                            progress.stop()
                            console.print(f"\n[grey50]Showing ion calculation for first microstate (microstate {ms_idx:03d}):[/grey50]")
                            self._display_ion_calculation(ion_counts)
                            if len(microstates) > 1:
                                console.print("[grey50]The same per-salt math is applied to each microstate with its own charge and water count.[/grey50]\n")
                            shown_educational_panel = True
                            progress.start()
                    else:
                        microstate_ion_counts[code] = None

                    progress.advance(task)

            # Show summary of charge range across all microstates
            charges = [ic['net_charge'] for ic in microstate_ion_counts.values() if ic]
            if charges:
                min_charge = min(charges)
                max_charge = max(charges)
                console.print(f"\n[green]Info pass complete for {len(microstates)} microstates[/green]")
                console.print(f"  Charge range: {min_charge:+d} to {max_charge:+d}")
                if min_charge != max_charge:
                    console.print(f"  [grey50]Ion counts vary per microstate due to different charges[/grey50]")

        console.print("\n[bold]Generating tLEaP input files...[/bold]")

        for ms_idx, microstate_info in enumerate(microstates, 1):
            code = microstate_info['code']
            ms_label = f"microstate_{ms_idx:03d}"
            script_filename = f"{ms_label}_tleap.in"

            ion_counts = microstate_ion_counts.get(code) if is_explicit else None

            success = self._generate_microstate_tleap_from_template(
                microstate_info, selected_forcefields, script_filename, template,
                ion_counts=ion_counts, enable_solvation=is_explicit,
                water_box=water_box, buffer=buffer, use_octahedron=use_octahedron,
                buffer_xyz=buffer_xyz, oct_diagonal=oct_diagonal, iso=iso
            )

            if success:
                console.print(f"[green]✓ Generated: {script_filename}[/green]")
                generated_inputs.append(script_filename)
            else:
                console.print(f"[red]✗ Failed: {script_filename}[/red]")
                failed_inputs.append(script_filename)
        
        # Report results
        if generated_inputs:
            console.print(f"\n[bold green]Successfully generated {len(generated_inputs)} tLEaP input files![/bold green]")
            for filename in generated_inputs:
                console.print(f"  ✓ {filename}")

            # Save generated tLEaP files to workspace for menu state tracking
            self.update_workspace("generated_microstate_tleap_files", generated_inputs)

        if failed_inputs:
            console.print(f"\n[bold red]Failed to generate {len(failed_inputs)} files:[/bold red]")
            for filename in failed_inputs:
                console.print(f"  ✗ {filename}")

        return len(generated_inputs) > 0
    
    def _generate_microstate_tleap_from_template(self, microstate_info, selected_forcefields, script_filename, template,
                                                     ion_counts=None, enable_solvation=False,
                                                     water_box="TIP3PBOX", buffer=10.0, use_octahedron=True,
                                                     buffer_xyz=None, oct_diagonal=0.0, iso=False):
        """
        Generate tLEaP input script for a microstate using template with placeholder substitution.

        Counter-ion identities and counts are read from `ion_counts` (the
        multi-salt dict produced by _calculate_multi_salt_ions). When
        `ion_counts` is None — e.g. info pass failed — falls back to a
        minimal Na+/Cl- neutralizer.
        """
        try:
            from pathlib import Path

            # Build all the auto-filled sections
            atom_types_section = self._build_microstate_atom_types_section(microstate_info, selected_forcefields)
            forcefield_section = self._build_microstate_forcefield_section(microstate_info, selected_forcefields)
            pdb_section = self._build_microstate_pdb_section(microstate_info)
            bond_section = self._build_microstate_bond_section(microstate_info)

            # Substitute placeholders in template
            final_content = template.replace("# ATOM_TYPES_SECTION", atom_types_section)
            final_content = final_content.replace("# FORCEFIELD_PARAMETERS_SECTION", forcefield_section)
            final_content = final_content.replace("# PDB_FILE_SECTION", pdb_section)
            final_content = final_content.replace("# BOND_DEFINITIONS_SECTION", bond_section)

            # Update the microstate label in the template (used for prmtop/rst7 filenames)
            # Derive label from PDB filename (e.g., "transformed_microstate_001")
            # or fall back to script filename stem
            pdb_stem = Path(microstate_info['filename']).stem
            ms_label = pdb_stem if pdb_stem else Path(script_filename).stem.replace('_tleap', '')
            final_content = final_content.replace("MICROSTATE", ms_label)

            # Replace solvation placeholder based on solvent model choice
            if enable_solvation:
                solvate_cmd = "solvateoct" if use_octahedron else "solvateBox"
                ion_commands = self._build_ion_commands_for_template(ion_counts)
                buf_arg = _format_buffer_for_tleap(buffer, buffer_xyz, oct_diagonal,
                                                    use_octahedron, iso)
                solvation_block = f"# === SOLVATION ===\n{solvate_cmd} mol {water_box} {buf_arg}\n\n{ion_commands}"
                final_content = final_content.replace("# SOLVATION_SECTION", solvation_block)
                # Also handle legacy templates that have separate ION placeholder
                final_content = final_content.replace("# SOLVATION_ION_SECTION", ion_commands)
            else:
                final_content = final_content.replace("# SOLVATION_SECTION", "# Implicit solvent - no solvation")

            # Write the final script
            with open(script_filename, 'w') as f:
                f.write(final_content)

            return True

        except Exception as e:
            logger.error(f"Error generating microstate tLEaP script for {microstate_info['code']}: {e}")
            return False
    
    def generate_microstate_inputs(self, metadata_file: str = None) -> bool:
        """
        Generate tLEaP input files for all redox microstates - uses template-based approach.
        
        Args:
            metadata_file: Path to microstate metadata JSON file. If None, will try to find it.
            
        Returns:
            True if successful, False otherwise
        """
        self.processor.console.print(
            "[yellow]Using template-based microstate tLEaP input generation...[/yellow]"
        )
        return self.generate_microstate_inputs_template_based(metadata_file)
    
    def _collect_unique_forcefield_requirements(self, microstates):
        """Collect unique cofactor+redox+spin combinations across all microstates."""
        unique_combinations = {}
        
        for microstate in microstates:
            for site in microstate['sites']:
                key = (site['transformer_type'], site['redox_state'], site['spin_state'])
                if key not in unique_combinations:
                    unique_combinations[key] = {
                        'transformer_type': site['transformer_type'],
                        'redox_state': site['redox_state'],
                        'spin_state': site['spin_state'],
                        'residue_name': site['residue_name'],
                        'atom_types': site['forcefield_info'].get('atom_types', [])
                    }
        
        return unique_combinations
    
    def _select_forcefields_for_microstates(self, unique_combinations):
        """Let user select forcefield files for each unique combination."""
        self._ensure_user_transformers_registered()
        console = self.processor.console
        selected_forcefields = {}
        
        console.print("\n[bold]Select Forcefield Parameters[/bold]")
        
        for key, info in unique_combinations.items():
            transformer_type, redox_state, spin_state = key

            # Skip forcefield selection for no_transformation - it doesn't need forcefields
            if transformer_type == 'no_transformation':
                console.print(f"\n[blue]Skipping {info['residue_name']} (no_transformation - no forcefield files needed)[/blue]")
                continue

            console.print(f"\n[blue]Forcefield for {info['residue_name']} ({transformer_type}, {redox_state}, {spin_state}):[/blue]")

            # Find available forcefield files
            try:
                from proprep.redoxsite_prep.transformation.redox_transformer_framework import redox_transformer_registry
                from proprep.forcefield_params import discover_forcefield_files
                from pathlib import Path
                from proprep.utils.prompts import prompt_with_context, confirm_with_context

                transformer_class = redox_transformer_registry.get_transformer(transformer_type)

                if not transformer_class or not getattr(transformer_class, 'FORCEFIELD_PATH', None):
                    console.print(f"[yellow]No forcefield files needed for {transformer_type}[/yellow]")
                    continue

                options = discover_forcefield_files(
                    transformer_class.FORCEFIELD_PATH, redox_state, spin_state
                )

                # Honour the Stage-1 fixed-pH/constant-pH choice (see the
                # single-state picker): show only sets matching the treatment a
                # site recorded for this combo, leaving just the charge-model
                # axis. Gated on a real choice + non-empty match so legacy /
                # treatment-agnostic cofactors fall through unchanged.
                combo_ph_treatment = self._preferred_ph_treatment_for_combo(
                    transformer_type, redox_state, spin_state)
                if combo_ph_treatment and options:
                    ph_matched = [
                        o for o in options
                        if o.get('ph_treatment') == combo_ph_treatment
                    ]
                    if ph_matched:
                        options = ph_matched

                if not options:
                    console.print(f"[yellow]No forcefield files found for {transformer_type}/{redox_state}/{spin_state}[/yellow]")
                    console.print(f"[yellow]Please add .frcmod and .lib files to the appropriate directory[/yellow]")
                    continue
                
                if len(options) == 1:
                    # Only one option, confirm with user
                    option = options[0]
                    console.print(f"Found forcefield set: {self._ff_set_display_title(option)}")
                    console.print(f"  - frcmod: {self._format_ff_file_basenames(option['frcmod'])}")
                    console.print(f"  - lib: {self._format_ff_file_basenames(option.get('lib'))}")
                    ph_phrase = self._ph_treatment_phrase(option.get('ph_treatment'))
                    if ph_phrase:
                        console.print(f"  - Propionate treatment: {ph_phrase}")
                    if option.get('description'):
                        console.print(f"  - Description: {option['description']}")
                    if option.get('version'):
                        console.print(f"  - Version: {option['version']}")
                    if option.get('reference'):
                        console.print(f"  - Reference: {option['reference']}")
                    
                    if confirm_with_context(
                        processor=self.processor,
                        prompt="Use this forcefield set?",
                        default=True,
                        module="Topology Generator",
                        description="Confirm forcefield selection"
                    ):
                        selected_forcefields[key] = option
                    else:
                        console.print("[yellow]Skipping this combination[/yellow]")
                        continue
                else:
                    # Multiple options, let user choose. Pre-select the set the
                    # site's Stage-1 pH-treatment choice implies (gated on an
                    # actual fork; no-op otherwise → existing option-1 default).
                    preferred_name = self._preferred_ff_set_for_combo(
                        transformer_type, redox_state, spin_state)
                    preferred_pos = next(
                        (i for i, o in enumerate(options, 1)
                         if o.get('name') == preferred_name), None)
                    default_choice = str(preferred_pos) if preferred_pos else "1"
                    # Name the chosen treatment once when the list is filtered to
                    # one; fall back to a per-line chip if a mix survives. Mirrors
                    # the single-state picker.
                    ph_present = {o.get('ph_treatment') for o in options}
                    show_ph_chip = len([p for p in ph_present if p]) > 1
                    combo_ph_phrase = self._ph_treatment_phrase(combo_ph_treatment)
                    if combo_ph_phrase and not show_ph_chip:
                        console.print(
                            "Forcefield sets for your chosen propionate treatment "
                            f"({combo_ph_phrase}) — differing only in charge model:")
                    else:
                        console.print("Multiple forcefield sets available:")
                    for i, option in enumerate(options, 1):
                        title = option.get('display_name') or option['name']
                        default_marker = " [default]" if option.get('is_default') else ""
                        ph_chip = (
                            {'fixed_pH': ' [fixed-pH]',
                             'constant_pH': ' [constant-pH]'}.get(
                                option.get('ph_treatment'), '')
                            if show_ph_chip else '')
                        preferred_marker = (
                            " [selected — matches your pH-treatment choice]"
                            if preferred_pos == i else "")
                        console.print(f"{i}. {title}{ph_chip}{default_marker}{preferred_marker}")
                        console.print(f"   - frcmod: {self._format_ff_file_basenames(option['frcmod'])}")
                        console.print(f"   - lib: {self._format_ff_file_basenames(option.get('lib'))}")
                        if option.get('description'):
                            console.print(f"   - Description: {option['description']}")
                        if option.get('version'):
                            console.print(f"   - Version: {option['version']}")
                        if option.get('reference'):
                            console.print(f"   - Reference: {option['reference']}")
                        console.print()  # Empty line between options

                    choice = prompt_with_context(
                        processor=self.processor,
                        prompt="Select forcefield set",
                        choices=[str(i) for i in range(1, len(options)+1)],
                        default=default_choice,
                        module="Topology Generator",
                        description="Select forcefield set"
                    )
                    # Replay by FF-set name so a changed/reordered option list
                    # re-selects the same set (options are dicts, not file paths).
                    choice = remap_recorded_index_by_key(
                        self.processor, options, lambda o: o.get('name', ''), str(choice))
                    selected_forcefields[key] = options[int(choice)-1]
                    annotate_recorded_key(self.processor, options[int(choice)-1].get('name', ''))
                    
            except Exception as e:
                console.print(f"[red]Error finding forcefields for {transformer_type}: {e}[/red]")
                continue
        
        return selected_forcefields
    
    def _get_bond_commands_for_microstate(self, microstate_info):
        """Get bond commands applicable to this microstate."""
        workspace = self.get_workspace()
        
        # Get existing combined bond commands
        combined_bonds = workspace.get("combined_tleap_commands", {})
        
        bond_commands = []
        
        # Add all bond types
        for bond_type, bonds in combined_bonds.items():
            if bonds:
                bond_commands.extend([f"# {bond_type.replace('_', ' ').title()}"])
                bond_commands.extend(bonds)
                bond_commands.append("")
        
        return bond_commands

    def _convert_redox_sites_to_tleap_commands(self, redox_sites):
        """
        Convert RedoxSite objects to tLEaP bond commands with intelligent backbone detection.
        
        Args:
            redox_sites: List of RedoxSite objects from workspace
            
        Returns:
            Dictionary of categorized bond commands
        """
        from proprep.utils.prompts import confirm_with_context

        # Normalize dict-form sites to RedoxSite objects. detected_redox_sites
        # round-trips through JSON workspace state, so a resumed session (or one
        # that imported a redox-sites JSON) can hand us dicts; the code below
        # accesses .site_id / .bonds / .atoms as attributes. The converter is
        # idempotent (objects pass through unchanged).
        from proprep.structure_prep.comprehensive_redox_detector import (
            dict_to_redox_site)
        redox_sites = [dict_to_redox_site(s) for s in redox_sites]

        # Initialize bond commands dictionary
        bond_commands = {
            "covalent": [],
            "coordinate": [],
            "metal-metal": [],
            "disulfide": [],
            "peptide_backbone": [],
            "other": []
        }
        
        # Track user decisions for backbone bonds to avoid repeated prompts
        user_backbone_decisions = {}  # {(site_id, chain, resid, resname): bool}
        
        # Global tracking of peptide bonds to prevent duplicates across all sites
        global_peptide_bonds = set()  # Track all peptide bonds across sites
        
        self.processor.console.print(f"\n[bold blue]Processing {len(redox_sites)} RedoxSite objects for tLEaP bond definitions...[/bold blue]")

        # Track statistics for summary
        sites_with_bonds = 0
        sites_without_bonds = 0
        total_structural_bonds = 0
        total_peptide_bonds = 0
        sites_processed = 0

        for site in redox_sites:
            site_id = site.site_id
            sites_processed += 1
            
            # DEBUG: Show coordinate mapping info for this site
            # self.processor.console.print(f"[blue]DEBUG: Site {site_id} has {len(site.coord_to_pdb)} coordinate mappings[/blue]")
            # self.processor.console.print(f"[blue]DEBUG: Site {site_id} has {len(site.atoms)} atoms[/blue]")
            
            # DEBUG: Check if coord_to_pdb should be populated from atoms
            # if len(site.coord_to_pdb) == 0 and len(site.atoms) > 0:
            #     self.processor.console.print(f"[red]DEBUG: coord_to_pdb is empty but {len(site.atoms)} atoms exist! This suggests import/loading issue.[/red]")
            #     sample_atoms = site.atoms[:3]
            #     for atom in sample_atoms:
            #         self.processor.console.print(f"[red]  Sample atom: {atom.coords} -> {atom.resname} {atom.resid} {atom.atom_name}[/red]")
            
            # if len(site.coord_to_pdb) > 0:
            #     # Show first few mappings as examples
            #     sample_coords = list(site.coord_to_pdb.keys())[:3]
            #     for coord in sample_coords:
            #         info = site.coord_to_pdb[coord]
            #         self.processor.console.print(f"[blue]  Sample mapping: {coord} -> {info.get('resname', '?')} {info.get('resid', '?')} {info.get('atom_name', '?')}[/blue]")
            
            # Step 1: Convert all RedoxSiteBond objects to tLEaP commands
            site_bond_count = 0
            site_bond_types = {"coordinate": 0, "covalent": 0, "disulfide": 0, "metal-metal": 0, "other": 0}

            for bond_idx, bond in enumerate(site.bonds):
                # A restrained metal-ligand contact (e.g. a coordinated water
                # held by an MD distance restraint) is realized as a nonbonded
                # residue: emit NO explicit tleap `bond` command for it. Bonding
                # it would make tleap build a metal-ligand bonded term and then
                # demand OW/HW-metal parameters the frcmod deliberately omits
                # ("Could not find bond/angle parameter: OW - M2"). The MD
                # restraint holds it in place instead. This mirrors the gates in
                # structure_preprocessor / fingerprint_generator so all bond
                # emitters agree.
                if getattr(bond, 'treatment', 'bonded') == 'restrained':
                    self.processor.console.print(
                        f"[grey50]  Skipping bond for restrained ligand "
                        f"(held by MD restraint, not a tleap bond).[/grey50]"
                    )
                    continue

                # Get current PDB info from coordinate mapping
                atom1_info = site.coord_to_pdb.get(bond.atom1_coords)
                atom2_info = site.coord_to_pdb.get(bond.atom2_coords)

                if not atom1_info or not atom2_info:
                    self.processor.console.print(f"[yellow]  Warning: Missing coordinate mapping for bond in {site_id}[/yellow]")

                    # DEBUG: Show what's missing
                    if not atom1_info:
                        self.processor.console.print(f"[red]  Missing atom1: {bond.atom1_coords}[/red]")
                    if not atom2_info:
                        self.processor.console.print(f"[red]  Missing atom2: {bond.atom2_coords}[/red]")

                    # DEBUG: Try to find close coordinates (coordinate precision issue)
                    self._debug_find_close_coordinates(site, bond.atom1_coords, bond.atom2_coords)
                    continue

                # Extract tLEaP command components
                resid1 = atom1_info['resid']
                atom1 = atom1_info['atom_name']
                resid2 = atom2_info['resid']
                atom2 = atom2_info['atom_name']

                # Safety net: skip bonds that ended up intra-residue after the
                # transformer's atom migration. The standard path is that the
                # transformer specs bond requirements only for bonds that will
                # be inter-residue post-transformation, so user-defined bonds
                # naturally fall into the explicit `bond` cmd output. But if a
                # user (or an outdated transformer spec) defined a bond between
                # atoms that have since been migrated into the same residue,
                # the residue's lib template already declares the bond and
                # emitting an explicit `bond` directive triggers tleap's
                # "1-4: cannot add bond" fatal error.
                chain1 = atom1_info.get('chain')
                chain2 = atom2_info.get('chain')
                ic1 = atom1_info.get('insertion_code', '')
                ic2 = atom2_info.get('insertion_code', '')
                if chain1 == chain2 and resid1 == resid2 and ic1 == ic2:
                    self.processor.console.print(
                        f"[grey50]  Skipping intra-residue bond "
                        f"{chain1}:{atom1_info.get('resname','?')}{resid1}:{atom1} "
                        f"-> {atom2_info.get('resname','?')}{resid2}:{atom2} "
                        f"(both endpoints migrated to the same residue; "
                        f"lib template already declares this bond).[/grey50]"
                    )
                    continue

                # Generate tLEaP command
                bond_cmd = f"bond mol.{resid1}.{atom1} mol.{resid2}.{atom2}"

                # Classify bonds based on comprehensive_redox_detector chemical_type definitions
                if bond.chemical_type == "covalent":
                    bond_commands["covalent"].append(bond_cmd)
                    site_bond_types["covalent"] += 1
                elif bond.chemical_type == "coordinate":
                    bond_commands["coordinate"].append(bond_cmd)
                    site_bond_types["coordinate"] += 1
                elif bond.chemical_type == "metal-metal":
                    bond_commands["metal-metal"].append(bond_cmd)
                    site_bond_types["metal-metal"] += 1
                elif bond.chemical_type == "disulfide":
                    bond_commands["disulfide"].append(bond_cmd)
                    site_bond_types["disulfide"] += 1
                else:
                    # Legacy or unrecognized types
                    bond_commands["other"].append(bond_cmd)
                    site_bond_types["other"] += 1

                site_bond_count += 1
            
            # Step 2: Intelligent backbone detection and peptide bond generation
            backbone_bonds_added = self._generate_peptide_backbone_bonds_for_site(
                site, bond_commands, user_backbone_decisions, global_peptide_bonds
            )

            # Print concise per-site summary with DEBUG info
            total_bonds = site_bond_count + backbone_bonds_added
            if total_bonds > 0:
                # Build bond type breakdown
                bond_parts = []
                if site_bond_types["coordinate"] > 0:
                    bond_parts.append(f"{site_bond_types['coordinate']} coordinate")
                if site_bond_types["covalent"] > 0:
                    bond_parts.append(f"{site_bond_types['covalent']} covalent")
                if site_bond_types["disulfide"] > 0:
                    bond_parts.append(f"{site_bond_types['disulfide']} disulfide")
                if site_bond_types["metal-metal"] > 0:
                    bond_parts.append(f"{site_bond_types['metal-metal']} metal-metal")
                if site_bond_types["other"] > 0:
                    bond_parts.append(f"{site_bond_types['other']} other")
                if backbone_bonds_added > 0:
                    bond_parts.append(f"{backbone_bonds_added} peptide")

                bond_breakdown = ", ".join(bond_parts)
                self.processor.console.print(f"[green]Site {site_id}: {total_bonds} bonds ({bond_breakdown})[/green]")

                sites_with_bonds += 1
                total_structural_bonds += site_bond_count
                total_peptide_bonds += backbone_bonds_added
            else:
                self.processor.console.print(f"[blue]Site {site_id}: 0 bonds (standard residues only)[/blue]")
                sites_without_bonds += 1
        
        # Remove duplicates from all categories
        for category in bond_commands:
            bond_commands[category] = list(dict.fromkeys(bond_commands[category]))

        # Print final summary
        total_bonds_after_dedup = sum(len(bonds) for bonds in bond_commands.values())
        self.processor.console.print(f"\n[bold green]Bond definitions complete:[/bold green]")
        self.processor.console.print(f"  • {sites_with_bonds} sites with bonds (Total: {total_bonds_after_dedup} bonds after deduplication)")
        self.processor.console.print(f"  • {sites_without_bonds} sites with standard residues only")
        if user_backbone_decisions:
            self.processor.console.print(f"  • {len(user_backbone_decisions)} residue type(s) required user decisions: {', '.join(user_backbone_decisions.keys())}")

        return bond_commands

    def _generate_peptide_backbone_bonds_for_site(self, site, bond_commands, user_decisions, global_peptide_bonds):
        """
        Generate peptide backbone bonds for residues in a RedoxSite with intelligent detection.
        
        Args:
            site: RedoxSite object
            bond_commands: Dictionary to add peptide bonds to
            user_decisions: Cache of user decisions to avoid repeated prompts
            global_peptide_bonds: Set to track peptide bonds across all sites (prevents duplicates)
            
        Returns:
            int: Number of peptide bonds added
        """
        from proprep.utils.prompts import confirm_with_context
        from proprep.forcefield_params import residue_has_template_connectivity

        # Group atoms by residue
        residues_in_site = {}
        for atom in site.atoms:
            res_key = (atom.chain, atom.resid, atom.resname, atom.insertion_code)
            if res_key not in residues_in_site:
                residues_in_site[res_key] = []
            residues_in_site[res_key].append(atom.atom_name)

        peptide_bonds_added = 0

        # Define standard amino acids - these are handled by forcefield
        STANDARD_AMINO_ACIDS = {
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
        }

        # Track residue names already covered by their .lib templates so we
        # only print the "skipping, lib has head/tail" note once per name.
        announced_template_skips = set()

        # Process each residue for potential peptide backbone bonds
        for (chain, resid, resname, icode), atom_names in residues_in_site.items():
            # Check for backbone atoms
            backbone_atoms = {"N", "CA", "C", "O"}
            has_all_backbone = backbone_atoms.issubset(set(atom_names))

            if has_all_backbone:
                # Skip standard amino acids - forcefield handles their peptide bonds
                if resname in STANDARD_AMINO_ACIDS:
                    continue

                # Skip non-standard residues whose bundled .lib already declares
                # head (connect0) + tail (connect1) atoms. tleap will auto-form
                # the peptide bonds from the residue templates, and emitting an
                # explicit `bond mol.<i>.C mol.<i+1>.N` would duplicate them
                # (silent no-op in newer tleap, but a hard crash in some
                # versions and a noisy redundancy at minimum).
                if residue_has_template_connectivity(resname):
                    if resname not in announced_template_skips:
                        self.processor.console.print(
                            f"[grey50]Skipping explicit peptide bonds for {resname}: "
                            f"its .lib template defines head/tail atoms, so tleap "
                            f"will auto-form the peptide bonds.[/grey50]"
                        )
                        announced_template_skips.add(resname)
                    continue

                # Non-standard residue with backbone atoms - ask user once per residue name
                decision_key = resname  # Cache by residue name only, not per-residue instance

                if decision_key not in user_decisions:
                    # First time seeing this residue type - ask user
                    self.processor.console.print(f"\n[yellow]Non-standard residue with backbone atoms:[/yellow] {resname} {resid} in Site {site.site_id}")
                    self.processor.console.print(f"[grey50]Residue {resname} has peptide backbone atoms (N, CA, C, O)[/grey50]")
                    self.processor.console.print(
                        f"[grey50]Select [bold]y[/bold] only if the chosen force field residue definition for "
                        f"{resname} does NOT specify head/tail atoms (connect0/connect1) that bond to the\n"
                        f"preceding and following residues. In that case ProPrep will auto-generate the explicit "
                        f"C-N peptide bonds connecting this residue to the rest of the (poly)peptide chain.\n"
                        f"Select [bold]n[/bold] if the template already declares head/tail atoms — tleap forms "
                        f"those peptide bonds itself, and an explicit bond would duplicate them.[/grey50]"
                    )
                    should_add_backbone = confirm_with_context(
                        processor=self.processor,
                        prompt=f"Should peptide backbone bonds be defined for ALL {resname} residues?",
                        default=True,  # Default True since backbone is present
                        module="Topology Generator",
                        description=f"Define backbone bonds for {resname}"
                    )
                    # Cache decision for this residue name
                    user_decisions[decision_key] = should_add_backbone
                else:
                    should_add_backbone = user_decisions[decision_key]

                if should_add_backbone:
                    # Generate C-N peptide bonds to adjacent residues IN THE PROTEIN STRUCTURE
                    bonds_for_residue = self._generate_backbone_bonds_for_residue_in_protein(
                        chain, resid, resname, icode
                    )

                    # Add to bond commands, avoiding global duplicates across all sites
                    for bond_cmd in bonds_for_residue:
                        if bond_cmd not in global_peptide_bonds:
                            bond_commands["peptide_backbone"].append(bond_cmd)
                            global_peptide_bonds.add(bond_cmd)  # Track globally
                            peptide_bonds_added += 1
                        # Silently skip duplicates - they'll be counted in the final summary
        
        return peptide_bonds_added

    def _generate_backbone_bonds_for_residue(self, chain, resid, resname, icode, all_residues):
        """
        Generate C-N peptide bonds from a specific residue to adjacent residues.
        
        Args:
            chain: Chain identifier
            resid: Residue number
            resname: Residue name
            icode: Insertion code
            all_residues: Dict of all residues in site with their atoms
            
        Returns:
            List of tLEaP bond commands
        """
        bonds = []

        # Check for C-terminal bond (current residue C to next residue N)
        next_resid = resid + 1

        # Find next residue (resname may differ)
        next_residue_key = None
        for res_key in all_residues.keys():
            if res_key[0] == chain and res_key[1] == next_resid and res_key[3] == icode:
                next_residue_key = res_key
                break

        if next_residue_key:
            next_atoms = all_residues[next_residue_key]
            current_atoms = all_residues[(chain, resid, resname, icode)]
            # Both residues must have required atoms
            if "C" in current_atoms and "N" in next_atoms:
                bond_cmd = f"bond mol.{resid}.C mol.{next_resid}.N"
                bonds.append(bond_cmd)

        # Check for N-terminal bond (previous residue C to current residue N)
        prev_resid = resid - 1
        if prev_resid > 0:  # Ensure valid residue number
            # Find previous residue
            prev_residue_key = None
            for res_key in all_residues.keys():
                if res_key[0] == chain and res_key[1] == prev_resid and res_key[3] == icode:
                    prev_residue_key = res_key
                    break

            if prev_residue_key:
                prev_atoms = all_residues[prev_residue_key]
                current_atoms = all_residues[(chain, resid, resname, icode)]
                # Both residues must have required atoms
                if "C" in prev_atoms and "N" in current_atoms:
                    bond_cmd = f"bond mol.{prev_resid}.C mol.{resid}.N"
                    bonds.append(bond_cmd)

        return bonds

    def _generate_backbone_bonds_for_residue_in_protein(self, chain, resid, resname, icode):
        """
        Generate C-N peptide bonds from a specific residue to adjacent residues in the entire protein structure.
        
        Args:
            chain: Chain identifier
            resid: Residue number
            resname: Residue name
            icode: Insertion code
            
        Returns:
            List of tLEaP bond commands
        """
        bonds = []
        
        # Get the protein structure from workspace to check for adjacent residues
        # Priority: reordered > protonation-updated > transformed > repaired > filtered > local-loaded > original
        # Must check reordered_pdb_file first to match the structure that RedoxSites were synced to
        structure = (self.get_from_workspace("reordered_pdb_file") or
                    self.get_from_workspace("protonation_pdb_file") or
                    self.get_from_workspace("transformed_pdb_file") or
                    self.get_from_workspace("transformed_structure") or
                    self.get_from_workspace("structure_with_prot_resnames") or  # legacy key
                    self.get_from_workspace("repaired_pdb_file") or
                    self.get_from_workspace("filtered_pdb_file") or
                    self.get_from_workspace("local_pdb_file") or
                    self.get_from_workspace("rcsb_pdb_file") or
                    self.get_from_workspace("original_pdb_file"))
        
        if not structure:
            self.processor.console.print(f"[yellow]Warning: No structure available for peptide bond generation[/yellow]")
            return bonds
        
        try:
            # Access the structure (could be file path or BioPython object)
            if isinstance(structure, str):
                # It's a file path, parse it
                from Bio.PDB import PDBParser
                parser = PDBParser(QUIET=True)
                structure = parser.get_structure("protein", structure)
            
            # Get the first model and chain
            model = structure[0]
            # Try the given chain ID; fall back to first chain if empty/missing
            try:
                protein_chain = model[chain]
            except KeyError:
                # Chain ID may be empty (' ') after MCPB preprocessing
                chains = list(model.get_chains())
                if chains:
                    protein_chain = chains[0]
                else:
                    return bonds
            
            # Check for C-terminal bond (current residue C to next residue N)
            next_resid = resid + 1
            if next_resid in protein_chain:
                next_residue = protein_chain[next_resid]
                # Check if both residues have required atoms
                current_residue = protein_chain[resid]
                if 'C' in current_residue and 'N' in next_residue:
                    bond_cmd = f"bond mol.{resid}.C mol.{next_resid}.N"
                    bonds.append(bond_cmd)
                
            # Check for N-terminal bond (previous residue C to current residue N)  
            prev_resid = resid - 1
            if prev_resid > 0 and prev_resid in protein_chain:
                prev_residue = protein_chain[prev_resid]
                current_residue = protein_chain[resid]
                # Check if both residues have required atoms
                if 'C' in prev_residue and 'N' in current_residue:
                    bond_cmd = f"bond mol.{prev_resid}.C mol.{resid}.N"
                    bonds.append(bond_cmd)
        
        except (KeyError, TypeError) as e:
            self.processor.console.print(f"[yellow]Warning: Could not access protein structure for residue {resname} {resid}: {e}[/yellow]")
        
        return bonds

    def _debug_find_close_coordinates(self, site, coord1, coord2):
        """
        Debug helper to find coordinates that are close to the missing ones.
        This helps identify coordinate precision issues.
        """
        import numpy as np
        
        def find_closest_coord(target_coord, coord_dict, threshold=0.01):
            """Find coordinates within threshold distance."""
            target = np.array(target_coord)
            close_coords = []
            
            for coord, info in coord_dict.items():
                distance = np.linalg.norm(np.array(coord) - target)
                if distance <= threshold:
                    close_coords.append((coord, info, distance))
            
            return sorted(close_coords, key=lambda x: x[2])  # Sort by distance
        
        # Check for close coordinates
        close1 = find_closest_coord(coord1, site.coord_to_pdb)
        close2 = find_closest_coord(coord2, site.coord_to_pdb)
        
        if close1:
            self.processor.console.print(f"[yellow]  Found {len(close1)} coordinates close to atom1:[/yellow]")
            for coord, info, dist in close1[:3]:  # Show top 3
                self.processor.console.print(f"[yellow]    {coord} -> {info.get('resname', '?')} {info.get('resid', '?')} {info.get('atom_name', '?')} (dist: {dist:.6f})[/yellow]")
        
        if close2:
            self.processor.console.print(f"[yellow]  Found {len(close2)} coordinates close to atom2:[/yellow]")
            for coord, info, dist in close2[:3]:  # Show top 3
                self.processor.console.print(f"[yellow]    {coord} -> {info.get('resname', '?')} {info.get('resid', '?')} {info.get('atom_name', '?')} (dist: {dist:.6f})[/yellow]")
        
        if not close1 and not close2:
            self.processor.console.print(f"[red]  No close coordinates found within 0.01Å threshold[/red]")

    def generate_topology_files(self) -> bool:
        """
        Generate prmtop/rst7 files from tLEaP input files in the working directory.
        
        This method:
        1. Finds all tLEaP input files (.in, .inp, .leap files)
        2. Parses each file to find the PDB file it references
        3. Shows a table of tLEaP/PDB pairs and their status
        4. Lets user select which tLEaP files to run
        5. Executes tleap for selected files
        6. Validates generated prmtop/rst7 files
        """
        import os
        import glob
        import subprocess
        from pathlib import Path
        from rich.table import Table
        from proprep.utils.prompts import prompt_with_context, confirm_with_context
        
        console = self.processor.console
        console.print("\n[bold blue]tLEaP Topology File Generation[/bold blue]")
        
        # Step 1: Find all tLEaP input files
        tleap_patterns = ["*.in", "*.inp", "*.leap", "*.tleap"]
        tleap_files = []
        
        for pattern in tleap_patterns:
            tleap_files.extend(sorted(glob.glob(pattern)))
        
        if not tleap_files:
            console.print("[red]No tLEaP input files found in current directory.[/red]")
            console.print("Expected extensions: .in, .inp, .leap, .tleap")
            return False
        
        console.print(f"Found {len(tleap_files)} tLEaP input files")
        
        # Step 2: Parse each tLEaP file to find referenced PDB files and expected outputs
        tleap_info = []
        
        for tleap_file in tleap_files:
            file_info = self._parse_tleap_input_file(tleap_file)
            if file_info:
                tleap_info.append(file_info)
        
        if not tleap_info:
            console.print("[red]No valid tLEaP input files found.[/red]")
            return False
        
        # Step 3: Show table with file status
        self._display_tleap_status_table(tleap_info)
        
        # Step 4: Let user select which files to run
        selected_files = self._select_tleap_files_to_run(tleap_info)
        
        if not selected_files:
            console.print("[yellow]No files selected for execution.[/yellow]")
            return False
        
        # Step 5: Run tLEaP for selected files
        success = self._execute_tleap_files(selected_files)
        
        # Step 6: Validate generated files
        if success:
            self._validate_generated_topology_files(selected_files)
        
        return success

    def _parse_tleap_input_file(self, tleap_file: str) -> dict:
        """Parse a tLEaP input file to extract PDB files and expected outputs."""
        try:
            with open(tleap_file, 'r') as f:
                content = f.read()
            
            info = {
                'tleap_file': tleap_file,
                'pdb_files': [],
                'expected_prmtop': None,
                'expected_rst7': None,
                'pdb_exists': [],
                'outputs_exist': {'prmtop': False, 'rst7': False}
            }
            
            lines = content.split('\n')
            for line in lines:
                line = line.strip()

                # Skip empty lines and comments
                if not line or line.startswith('#'):
                    continue

                # Look for PDB loading commands
                if 'loadpdb' in line.lower() or 'loadPdb' in line:
                    # Extract PDB filename: loadpdb mol filename.pdb
                    parts = line.split()
                    if len(parts) >= 3:
                        pdb_file = parts[-1]  # Last part should be filename
                        info['pdb_files'].append(pdb_file)
                        info['pdb_exists'].append(os.path.exists(pdb_file))
                
                # Look for output commands
                elif 'saveamberparm' in line.lower():
                    # Extract output files: saveamberparm mol output.prmtop output.rst7
                    parts = line.split()
                    if len(parts) >= 4:
                        info['expected_prmtop'] = parts[-2]  # Second to last
                        info['expected_rst7'] = parts[-1]   # Last
                        info['outputs_exist']['prmtop'] = os.path.exists(parts[-2])
                        info['outputs_exist']['rst7'] = os.path.exists(parts[-1])
            
            return info if info['pdb_files'] else None
            
        except Exception as e:
            self.processor.console.print(f"[red]Error parsing {tleap_file}: {e}[/red]")
            return None

    def _display_tleap_status_table(self, tleap_info: list):
        """Display a table showing tLEaP file status."""
        table = Table(title="tLEaP Input Files Status")
        table.add_column("tLEaP File", style="blue")
        table.add_column("PDB Files", style="blue")
        table.add_column("PDB Status", style="green")
        table.add_column("Expected Output", style="yellow")
        table.add_column("Output Status", style="magenta")
        
        for info in tleap_info:
            # Format PDB files and status
            pdb_display = []
            for pdb_file, exists in zip(info['pdb_files'], info['pdb_exists']):
                status = "✓" if exists else "✗"
                color = "green" if exists else "red"
                pdb_display.append(f"[{color}]{pdb_file} {status}[/{color}]")
            
            pdb_files_str = "\n".join(pdb_display)
            pdb_status = "All exist" if all(info['pdb_exists']) else "Missing files"
            pdb_status_color = "green" if all(info['pdb_exists']) else "red"
            
            # Format expected outputs
            outputs = []
            if info['expected_prmtop']:
                outputs.append(info['expected_prmtop'])
            if info['expected_rst7']:
                outputs.append(info['expected_rst7'])
            expected_output = "\n".join(outputs)
            
            # Format output status
            prmtop_status = "✓" if info['outputs_exist']['prmtop'] else "✗"
            rst7_status = "✓" if info['outputs_exist']['rst7'] else "✗"
            output_status = f"prmtop {prmtop_status} rst7 {rst7_status}"
            
            table.add_row(
                info['tleap_file'],
                pdb_files_str,
                f"[{pdb_status_color}]{pdb_status}[/{pdb_status_color}]",
                expected_output,
                output_status
            )
        
        self.processor.console.print(table)

    def _select_tleap_files_to_run(self, tleap_info: list) -> list:
        """Let user select which tLEaP files to run."""
        console = self.processor.console

        # Filter to files that have all required PDBs
        runnable_files = [info for info in tleap_info if all(info['pdb_exists'])]

        if not runnable_files:
            console.print("[red]No tLEaP files can be run (missing PDB files).[/red]")
            return []

        # Check which files already have topology files generated
        files_with_topology = []
        files_without_topology = []
        for info in runnable_files:
            prmtop = info.get('expected_prmtop')
            rst7 = info.get('expected_rst7')
            # Consider complete if prmtop exists (rst7 is usually generated together)
            if prmtop and os.path.exists(prmtop):
                files_with_topology.append(info)
            else:
                files_without_topology.append(info)

        num_total = len(runnable_files)
        num_done = len(files_with_topology)
        num_missing = len(files_without_topology)

        console.print(f"\n[bold]Found {num_total} tLEaP files ready to run.[/bold]")
        if num_done > 0:
            console.print(f"  [green]• {num_done} already have topology files[/green]")
        if num_missing > 0:
            console.print(f"  [yellow]• {num_missing} need generation[/yellow]")

        # Single file case
        if num_total == 1:
            if num_done == 1:
                if confirm_with_context(
                    processor=self.processor,
                    prompt=f"Topology already exists. Regenerate {os.path.basename(runnable_files[0]['tleap_file'])}?",
                    default=False,
                    module="Topology Generator",
                    description="Regenerate existing topology"
                ):
                    return runnable_files
                return []
            else:
                if confirm_with_context(
                    processor=self.processor,
                    prompt=f"Run {os.path.basename(runnable_files[0]['tleap_file'])}?",
                    default=True,
                    module="Topology Generator",
                    description="Generate topology"
                ):
                    return runnable_files
                return []

        # Multiple files - show options based on what's available
        console.print("\n[bold]Run options:[/bold]")
        options = []
        choices = []

        if num_missing > 0 and num_done > 0:
            # Both missing and existing - offer choice
            options.append(("missing", f"Generate missing only ({num_missing} files)"))
            options.append(("all", f"Regenerate all ({num_total} files)"))
            options.append(("select", "Select specific files"))
            options.append(("skip", "Skip - continue without regenerating"))
        elif num_missing > 0:
            # All need generation
            options.append(("all", f"Generate all ({num_total} files)"))
            options.append(("select", "Select specific files"))
            options.append(("skip", "Skip - continue without generating"))
        else:
            # All already exist
            console.print("[green]All topology files already exist.[/green]")
            options.append(("all", f"Regenerate all ({num_total} files)"))
            options.append(("select", "Select specific files"))
            options.append(("skip", "Skip - all topologies already exist"))

        for i, (key, label) in enumerate(options, 1):
            console.print(f"  {i}. {label}")
            choices.append(str(i))

        choice = prompt_with_context(
            processor=self.processor,
            prompt="Select option",
            choices=choices,
            default=str(len(options)),  # Default to skip
            module="Topology Generator",
            description="Topology generation options"
        )

        try:
            choice_num = int(choice)
            if 1 <= choice_num <= len(options):
                selected_option = options[choice_num - 1][0]

                if selected_option == "skip":
                    return []
                elif selected_option == "missing":
                    return files_without_topology
                elif selected_option == "all":
                    return runnable_files
                elif selected_option == "select":
                    return self._select_specific_tleap_files(runnable_files, files_with_topology)
        except ValueError:
            pass

        return []

    def _select_specific_tleap_files(self, runnable_files: list, files_with_topology: list) -> list:
        """Let user select specific tLEaP files to run."""
        console = self.processor.console
        files_with_topology_set = set(id(f) for f in files_with_topology)

        console.print("\nSelect files to run (comma-separated numbers, 'all', or 'skip'):")
        for i, info in enumerate(runnable_files, 1):
            base_name = os.path.basename(info['tleap_file'])
            status = "[green]✓[/green]" if id(info) in files_with_topology_set else "[yellow]○[/yellow]"
            console.print(f"  {i:2d}. {status} {base_name}")

        console.print("\n[grey50]✓ = topology exists, ○ = needs generation[/grey50]")

        choice = prompt_with_context(
            processor=self.processor,
            prompt="Enter file numbers (e.g., 1,3,5 or 'all' or 'skip')",
            default="skip",
            module="Topology Generator",
            description="Select specific topology files"
        )

        if choice.lower() in ['skip', 's', 'cancel', 'c']:
            return []
        elif choice.lower() in ['all', 'a']:
            return runnable_files
        else:
            try:
                indices = [int(x.strip()) for x in choice.split(',')]
                selected = []
                valid = True
                for idx in indices:
                    if 1 <= idx <= len(runnable_files):
                        selected.append(runnable_files[idx - 1])
                    else:
                        console.print(f"[red]Invalid number: {idx}[/red]")
                        valid = False
                        break
                if valid and selected:
                    return selected
            except ValueError:
                console.print("[red]Please enter comma-separated numbers.[/red]")

        return []

    def _prompt_parallel_workers(self, num_files: int) -> int:
        """Prompt user for number of parallel tLEaP workers."""
        console = self.processor.console
        cpu_count = os.cpu_count() or 1

        console.print(f"\n[bold]Parallel Execution[/bold]")
        console.print(f"[grey50]{num_files} tLEaP files to process. Available CPUs: {cpu_count}[/grey50]")

        while True:
            workers_str = prompt_with_context(
                processor=self.processor,
                prompt="Number of parallel tLEaP executions",
                default="1",
                module="Topology Generator",
                description="Enter number of parallel workers"
            )
            try:
                workers = int(workers_str)
                if 1 <= workers <= num_files:
                    return workers
                console.print(f"[yellow]Please enter a number between 1 and {num_files}[/yellow]")
            except ValueError:
                console.print("[yellow]Please enter a number[/yellow]")

    def _swap_md_pair_to_cpin_prmtop(self, original_prmtop: str,
                                        cpin_prmtop: str) -> bool:
        """After cpinutil writes a `_cpin.prmtop` for explicit-solvent CpHMD,
        update the matching md_structure_pairs entry to point at it.

        The raw constph prmtop produced by tleap is the *input* to cpinutil,
        not a production-ready topology. _register_topology_outputs stored
        the raw prmtop in md_structure_pairs when tleap finished. Once
        cpinutil has run with -op, the cpin-modified prmtop is the correct
        target for any production CpHMD launch (it carries the corrected
        carboxylate radii needed by sander/pmemd's PB-flavor electrostatics),
        so swap the entry's `prmtop` field in place. The rst7 is unchanged
        (cpinutil does not modify coordinates).

        Matches by resolved absolute path. Returns True if a swap occurred.
        """
        if not (original_prmtop and cpin_prmtop):
            return False
        if not os.path.exists(cpin_prmtop):
            return False
        original_abs = os.path.abspath(original_prmtop)
        cpin_abs = os.path.abspath(cpin_prmtop)
        pairs = self.get_from_workspace("md_structure_pairs", []) or []
        swapped = False
        for entry in pairs:
            try:
                entry_prmtop_abs = os.path.abspath(entry.get("prmtop", ""))
            except (TypeError, ValueError):
                continue
            if entry_prmtop_abs != original_abs:
                continue
            entry["prmtop"] = cpin_abs
            entry["name"] = os.path.splitext(os.path.basename(cpin_abs))[0]
            swapped = True
        if swapped:
            self.update_workspace("md_structure_pairs", pairs)
        return swapped

    def _register_topology_outputs(self, info: dict) -> bool:
        """Record a freshly generated prmtop/rst7 pair in the workspace.

        Updates `parm7_file`/`rst7_file` (overwritten so the most recent
        success wins) and appends a deduplicated entry to
        `md_structure_pairs`, which the MD Manager's Step 0 picker consumes
        to enumerate all available microstate topologies.
        """
        prmtop = info.get('expected_prmtop')
        if not prmtop or not os.path.exists(prmtop):
            return False
        prmtop_abs = os.path.abspath(prmtop)
        rst7 = info.get('expected_rst7')
        rst7_abs = os.path.abspath(rst7) if rst7 and os.path.exists(rst7) else None

        self.update_workspace("parm7_file", prmtop_abs)
        if rst7_abs:
            self.update_workspace("rst7_file", rst7_abs)

            existing_pairs = self.get_from_workspace("md_structure_pairs", []) or []
            new_entry = {
                "name": os.path.splitext(os.path.basename(prmtop_abs))[0],
                "prmtop": prmtop_abs,
                "rst7": rst7_abs,
            }
            if not any(
                p.get("prmtop") == new_entry["prmtop"]
                and p.get("rst7") == new_entry["rst7"]
                for p in existing_pairs
            ):
                existing_pairs.append(new_entry)
                self.update_workspace("md_structure_pairs", existing_pairs)
        return True

    def _execute_tleap_files(self, selected_files: list) -> bool:
        """
        Execute tLEaP for selected files (single-pass - templates already have accurate ions).

        Since template generation now runs an info pass to calculate accurate ion counts,
        topology generation only needs a single tLEaP run.
        """
        console = self.processor.console
        success_count = 0

        # Ask for parallel workers when multiple files
        max_workers = 1
        if len(selected_files) > 1:
            max_workers = self._prompt_parallel_workers(len(selected_files))

        # Use batch mode for multiple files (more than 3) or parallel execution
        batch_mode = len(selected_files) > 3 or max_workers > 1

        if batch_mode:
            return self._execute_tleap_files_batch(selected_files, max_workers)

        for info in selected_files:
            tleap_file = info['tleap_file']
            console.print(f"\n[blue]Running tLEaP for {tleap_file}...[/blue]")

            # Check and configure molecule grouping (if needed)
            molecule_config_success = self._check_and_configure_molecules(info)
            if not molecule_config_success:
                console.print(f"[red]✗ Molecule configuration cancelled for {tleap_file}[/red]")
                continue

            # Validate and fix TER records in PDB files before running tLEaP
            ter_validation_success = self._validate_and_fix_ter_records(info)
            if not ter_validation_success:
                console.print(f"[red]✗ TER record validation failed for {tleap_file} - incorrect bonds will be generated in tLEaP[/red]")
                continue

            try:
                # Single-pass execution (templates already have accurate ion counts)
                console.print(f"[blue]Running: tleap -s -f {tleap_file}[/blue]")
                result = subprocess.run(
                    ["tleap", "-s", "-f", tleap_file],
                    capture_output=True,
                    text=True
                )

                if result.returncode == 0:
                    console.print(f"[green]✓ tLEaP completed successfully[/green]")

                    # Parse and display leap.log messages
                    self._display_leap_log_messages(tleap_file)

                    # Run ParmEd validation
                    if info.get('expected_prmtop') and os.path.exists(info['expected_prmtop']):
                        rst7_file = info.get('expected_rst7') if info.get('expected_rst7') and os.path.exists(info.get('expected_rst7')) else None
                        self._run_parmed_validation(info['expected_prmtop'], rst7_file)

                        # Hand off topology paths to downstream modules (MD Manager, QM/MM).
                        # Appends to md_structure_pairs so batches of microstates all survive.
                        self._register_topology_outputs(info)

                    success_count += 1
                else:
                    console.print(f"[red]✗ tLEaP failed (exit code: {result.returncode})[/red]")
                    self._display_leap_log_messages(tleap_file)
                    if result.stderr:
                        console.print(f"[red]{result.stderr}[/red]")

            except FileNotFoundError:
                console.print(f"[red]✗ tLEaP command not found. Please ensure tLEaP is installed and in PATH.[/red]")
                return False
            except Exception as e:
                console.print(f"[red]✗ Error running tLEaP for {tleap_file}: {e}[/red]")

        console.print(f"\n[bold]Results: {success_count}/{len(selected_files)} files processed successfully[/bold]")
        return success_count > 0

    def _execute_tleap_two_pass(self, info: dict, tleap_file: str) -> bool:
        """
        Execute tLEaP in two passes with early termination on pass 1.

        Pass 1 (Info Pass): Run tLEaP to get water count and net charge,
                           then KILL immediately once info is captured
        Pass 2 (Production): Full run with correct ion counts

        Args:
            info: File info dictionary with template details
            tleap_file: Path to original tLEaP input file

        Returns:
            True if successful, False otherwise
        """
        console = self.processor.console

        console.print("[bold yellow]Two-pass execution: Gathering system info...[/bold yellow]")
        console.print("[grey50]Pass 1: Info gathering (with early termination)[/grey50]")

        # Pass 1: Get system info with early termination
        system_info = self._execute_tleap_info_pass(tleap_file)

        n_waters = system_info.get('n_waters')
        net_charge = system_info.get('net_charge')

        if n_waters is None:
            console.print("[yellow]⚠ Could not determine water count, using original template[/yellow]")
            final_template = tleap_file
        else:
            # Use charge from tLEaP if available, else fall back to workspace
            if net_charge is None:
                net_charge = self.get_from_workspace("net_charge", 0) or 0
                console.print(f"[grey50]Using workspace charge: {net_charge:+.0f}[/grey50]")

            # Look up the user's salt selection (back-compat for legacy single-salt blobs)
            sp = self.get_from_workspace("solvation_parameters", {}) or {}
            salts, neutralize_index = _normalize_salts(sp)

            # Calculate accurate ion counts
            ion_counts = self._calculate_multi_salt_ions(
                n_waters, int(round(net_charge)),
                salts=salts, neutralize_index=neutralize_index
            )

            console.print(f"[blue]Calculated ions ({_salts_summary(salts)}):[/blue]")
            console.print(f"  Waters: {n_waters}")
            console.print(f"  Net charge: {ion_counts['net_charge']:+d}")
            for s in ion_counts.get('salts', []):
                lbl = _salt_label(s['cation'], s['anion'])
                conc_mM = s['concentration'] * 1000
                console.print(f"  {lbl} ({conc_mM:.0f} mM): {s['n_pairs']} formula units "
                              f"-> {s['n_bulk_cations']} {s['cation']['label']} + "
                              f"{s['n_bulk_anions']} {s['anion']['label']}")
            n_neut_cat = ion_counts.get('n_neutralize_cation', 0)
            n_neut_an = ion_counts.get('n_neutralize_anion', 0)
            if n_neut_cat or n_neut_an:
                neut_salt = ion_counts['salts'][ion_counts['neutralize_index']]
                console.print(f"  Neutralizer (via {_salt_label(neut_salt['cation'], neut_salt['anion'])}): "
                              f"{n_neut_cat} {neut_salt['cation']['label']} + {n_neut_an} {neut_salt['anion']['label']}")

            # Update template with accurate ion counts
            final_template = self._update_template_ion_counts(tleap_file, ion_counts)
            if not final_template or final_template == tleap_file:
                console.print("[yellow]⚠ Could not update ion counts, using original template[/yellow]")
                final_template = tleap_file

        # Pass 2: Full run with accurate ions
        console.print("\n[grey50]Pass 2: Full topology generation[/grey50]")
        cmd = ["tleap", "-s", "-f", final_template]
        console.print(f"[blue]Running: {' '.join(cmd)}[/blue]")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            console.print(f"[green]✓ tLEaP completed successfully[/green]")

            # Parse and display leap.log messages
            self._display_leap_log_messages(final_template)

            # Run ParmEd validation
            if info.get('expected_prmtop') and os.path.exists(info['expected_prmtop']):
                rst7_file = info.get('expected_rst7') if info.get('expected_rst7') and os.path.exists(info.get('expected_rst7')) else None
                self._run_parmed_validation(info['expected_prmtop'], rst7_file)

            return True
        else:
            console.print(f"[red]✗ Pass 2 failed (exit code: {result.returncode})[/red]")
            self._display_leap_log_messages(final_template)
            if result.stderr:
                console.print(f"[red]{result.stderr}[/red]")
            return False

    def _execute_tleap_info_pass(self, tleap_file: str) -> dict:
        """
        Run tLEaP info pass with early termination.

        Monitors tLEaP output in real-time and kills it immediately once we have:
        - Water count (from solvate command output)
        - Net charge (from 'charge' command output)

        Args:
            tleap_file: Path to original tLEaP input file

        Returns:
            Dictionary with 'n_waters', 'net_charge', 'volume' (values may be None)
        """
        import re

        console = self.processor.console

        # Create info-pass template with charge command added
        info_template = self._create_info_pass_template(tleap_file)
        if not info_template:
            console.print("[red]✗ Failed to create info pass template[/red]")
            return {'n_waters': None, 'net_charge': None, 'volume': None}

        info = {'n_waters': None, 'net_charge': None, 'volume': None}
        required_fields = {'n_waters', 'net_charge'}  # Minimum needed for ion calculation

        # Patterns to match in tLEaP output
        patterns = {
            'water_count': re.compile(r'Added\s+(\d+)\s+residues\.'),
            'charge': re.compile(r'Total unperturbed charge:\s*([-\d.]+)'),
            'volume': re.compile(r'Volume:\s*([\d.]+)\s+A\^3'),
        }

        proc = None
        try:
            # Start tLEaP with line-buffered output
            proc = subprocess.Popen(
                ['tleap', '-f', info_template],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1  # Line buffered for real-time reading
            )

            # Monitor output line by line
            for line in iter(proc.stdout.readline, ''):
                if not line:  # Empty string means EOF
                    break

                # Parse water count from solvate output
                if match := patterns['water_count'].search(line):
                    info['n_waters'] = int(match.group(1))
                    console.print(f"[grey50]  → Water count: {info['n_waters']}[/grey50]")

                # Parse charge from 'charge mol' output
                if match := patterns['charge'].search(line):
                    info['net_charge'] = float(match.group(1))
                    console.print(f"[grey50]  → Net charge: {info['net_charge']:+.1f}[/grey50]")

                # Parse volume (optional, for information)
                if match := patterns['volume'].search(line):
                    info['volume'] = float(match.group(1))

                # Check if we have all required info - TERMINATE EARLY
                if all(info.get(field) is not None for field in required_fields):
                    console.print("[green]  ✓ Got required info, terminating tLEaP early[/green]")
                    proc.terminate()
                    try:
                        proc.wait(timeout=2)  # Give it 2 seconds to clean up
                    except subprocess.TimeoutExpired:
                        proc.kill()  # Force kill if it doesn't respond
                        proc.wait()
                    break

            # If we exited the loop but process is still running, clean up
            if proc.poll() is None:
                proc.terminate()
                try:
                    proc.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()

        except FileNotFoundError:
            console.print("[red]tLEaP not found in PATH[/red]")
            return info
        except Exception as e:
            console.print(f"[red]Error in info pass: {e}[/red]")
            if proc is not None and proc.poll() is None:
                proc.kill()
                proc.wait()
        finally:
            # Clean up temp file
            if info_template and os.path.exists(info_template):
                try:
                    os.remove(info_template)
                except Exception:
                    pass

        return info

    def _create_info_pass_template(self, original_template: str) -> str:
        """
        Create a tLEaP template for info gathering with charge command.

        This template:
        1. Loads forcefields and structure (same as original)
        2. Adds 'charge mol' BEFORE solvation (structure charge)
        3. Runs solvation command (to get water count)
        4. Adds 'charge mol' AFTER solvation (total system charge)
        5. Skips ion commands and file saving (we'll terminate early anyway)

        Args:
            original_template: Path to original tLEaP input file

        Returns:
            Path to info pass template, or None on failure
        """
        try:
            with open(original_template, 'r') as f:
                lines = f.readlines()

            # Create info template path
            base = os.path.splitext(original_template)[0]
            info_file = f"{base}_info_pass.in"

            with open(info_file, 'w') as f:
                f.write("# INFO PASS - will be terminated early after gathering water count and charge\n")

                found_solvate = False
                mol_name = "mol"  # Default molecule name, will try to detect

                for line in lines:
                    stripped = line.strip().lower()

                    # Try to detect molecule name from loadpdb/loadmol2
                    if 'loadpdb' in stripped or 'loadmol2' in stripped:
                        # Parse: mol = loadpdb filename.pdb
                        parts = line.split('=')
                        if len(parts) >= 1:
                            mol_name = parts[0].strip()
                        f.write(line)
                        continue

                    # Insert charge command before solvate
                    if ('solvate' in stripped) and not found_solvate:
                        found_solvate = True
                        f.write(f"\n# Get charge before solvation\n")
                        f.write(f"charge {mol_name}\n\n")
                        f.write(line)  # Write solvate command
                        f.write(f"\n# Get charge after solvation (includes water)\n")
                        f.write(f"charge {mol_name}\n")
                        continue

                    # Skip ion commands - we don't need them for info pass
                    if stripped.startswith('addions') or stripped.startswith('addionsrand'):
                        f.write(f"# [SKIPPED FOR INFO PASS] {line}")
                        continue

                    # Skip file saving commands
                    if stripped.startswith('saveamberparm') or stripped.startswith('savepdb'):
                        f.write(f"# [SKIPPED FOR INFO PASS] {line}")
                        continue

                    # Write all other lines as-is
                    f.write(line)

                # Ensure quit at end (though we'll likely kill it before this)
                f.write("\nquit\n")

            logger.info(f"Created info pass template: {info_file}")
            return info_file

        except Exception as e:
            logger.error(f"Error creating info pass template: {e}")
            return None

    def _execute_tleap_info_pass_quiet(self, tleap_file: str) -> dict:
        """
        Run tLEaP info pass with early termination (quiet version for batch mode).

        Same as _execute_tleap_info_pass but without console output.

        Args:
            tleap_file: Path to original tLEaP input file

        Returns:
            Dictionary with 'n_waters', 'net_charge', 'volume' (values may be None)
        """
        import re

        # Create info-pass template with charge command added
        info_template = self._create_info_pass_template(tleap_file)
        if not info_template:
            return {'n_waters': None, 'net_charge': None, 'volume': None}

        info = {'n_waters': None, 'net_charge': None, 'volume': None}
        required_fields = {'n_waters', 'net_charge'}

        # Patterns to match in tLEaP output
        patterns = {
            'water_count': re.compile(r'Added\s+(\d+)\s+residues\.'),
            'charge': re.compile(r'Total unperturbed charge:\s*([-\d.]+)'),
            'volume': re.compile(r'Volume:\s*([\d.]+)\s+A\^3'),
        }

        proc = None
        try:
            proc = subprocess.Popen(
                ['tleap', '-f', info_template],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1
            )

            for line in iter(proc.stdout.readline, ''):
                if not line:
                    break

                if match := patterns['water_count'].search(line):
                    info['n_waters'] = int(match.group(1))

                if match := patterns['charge'].search(line):
                    info['net_charge'] = float(match.group(1))

                if match := patterns['volume'].search(line):
                    info['volume'] = float(match.group(1))

                # Terminate early once we have required info
                if all(info.get(field) is not None for field in required_fields):
                    proc.terminate()
                    try:
                        proc.wait(timeout=2)
                    except subprocess.TimeoutExpired:
                        proc.kill()
                        proc.wait()
                    break

            if proc.poll() is None:
                proc.terminate()
                try:
                    proc.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()

        except Exception:
            if proc is not None and proc.poll() is None:
                proc.kill()
                proc.wait()
        finally:
            if info_template and os.path.exists(info_template):
                try:
                    os.remove(info_template)
                except Exception:
                    pass

        return info

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Info Pass for Template Generation
    # Run tLEaP during template generation to get accurate water count and charge
    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#

    def _build_info_pass_script_for_template(self, pdb_path: str, water_box: str = "TIP3PBOX",
                                              buffer: float = 10.0, use_octahedron: bool = True,
                                              buffer_xyz=None, oct_diagonal: float = 0.0,
                                              iso: bool = False) -> str:
        """
        Build a minimal tLEaP script for info gathering during template generation.

        This script:
        1. Loads standard forcefields (selected by user)
        2. Loads custom forcefields (from redox site parameterization)
        3. Loads the PDB structure
        4. Applies bond definitions
        5. Runs solvation (to get water count)
        6. Runs 'charge mol' (to get net charge)
        7. Skips ion commands and file saving (we'll terminate early anyway)

        Args:
            pdb_path: Path to the PDB file
            water_box: Water box type (e.g., TIP3PBOX, OPCBOX)
            buffer: Buffer distance in Angstroms
            use_octahedron: If True, use solvateoct; else use solvateBox

        Returns:
            Path to the info pass script file
        """
        import tempfile

        # Build forcefield section from user selections
        forcefield_section, _ = self._build_standard_forcefield_section()

        # Build custom forcefield section
        custom_ff_section = self._build_forcefield_parameters_section()

        # Build atom types section
        atom_types_section = self._build_atom_types_section()

        # Build bond definitions section
        bond_section = self._build_bond_definitions_section()

        # Solvation command
        solvate_cmd = "solvateoct" if use_octahedron else "solvateBox"
        buf_arg = _format_buffer_for_tleap(buffer, buffer_xyz, oct_diagonal,
                                            use_octahedron, iso)

        script = f"""# INFO PASS - gathering water count and charge for template generation
# This script will be terminated early once we have the required information

# === STANDARD FORCEFIELDS ===
{forcefield_section}

# === CUSTOM ATOM TYPES ===
{atom_types_section}

# === CUSTOM FORCEFIELD PARAMETERS ===
{custom_ff_section}

# === LOAD STRUCTURE ===
mol = loadpdb {pdb_path}

# === BOND DEFINITIONS ===
{bond_section}

# === GET CHARGE BEFORE SOLVATION ===
charge mol

# === SOLVATION ===
{solvate_cmd} mol {water_box} {buf_arg}

# === GET CHARGE AFTER SOLVATION ===
charge mol

quit
"""

        # Write to temp file
        fd, temp_path = tempfile.mkstemp(suffix='_info_pass.in', prefix='tleap_')
        os.close(fd)
        with open(temp_path, 'w') as f:
            f.write(script)

        return temp_path

    def _run_info_pass_for_template(self, pdb_path: str, water_box: str = "TIP3PBOX",
                                     buffer: float = 10.0, use_octahedron: bool = True,
                                     quiet: bool = False, buffer_xyz=None,
                                     oct_diagonal: float = 0.0, iso: bool = False) -> dict:
        """
        Run tLEaP info pass during template generation to get accurate water count and charge.

        This is used during template generation (Options 2/3) to get real values from tLEaP
        instead of the vdW estimation. The info pass terminates early once we have the
        required information.

        Args:
            pdb_path: Path to the PDB file
            water_box: Water box type (e.g., TIP3PBOX, OPCBOX)
            buffer: Buffer distance in Angstroms
            use_octahedron: If True, use solvateoct; else use solvateBox
            quiet: If True, suppress console output

        Returns:
            Dictionary with 'n_waters', 'net_charge', 'volume' (values may be None)
        """
        import re

        console = self.processor.console if not quiet else None

        # Build info pass script
        info_script = self._build_info_pass_script_for_template(
            pdb_path, water_box, buffer, use_octahedron,
            buffer_xyz=buffer_xyz, oct_diagonal=oct_diagonal, iso=iso
        )
        if not info_script:
            if console:
                console.print("[red]✗ Failed to create info pass script[/red]")
            return {'n_waters': None, 'net_charge': None, 'volume': None}

        # Log the script path for debugging
        logger.debug("Info pass script: %s", info_script)

        info = {'n_waters': None, 'net_charge': None, 'volume': None}
        required_fields = {'n_waters', 'net_charge'}

        # Patterns to match in tLEaP output
        patterns = {
            'water_count': re.compile(r'Added\s+(\d+)\s+residues\.'),
            'charge': re.compile(r'Total unperturbed charge:\s*([-\d.]+)'),
            'volume': re.compile(r'Volume:\s*([\d.]+)\s+A\^3'),
        }

        proc = None
        output_lines = []  # Capture all output for diagnostics
        try:
            # Start tLEaP with line-buffered output
            # Use DEVNULL for stdin to prevent tLEaP from reading terminal on error
            proc = subprocess.Popen(
                ['tleap', '-f', info_script],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                stdin=subprocess.DEVNULL,
                text=True,
                bufsize=1
            )

            # Monitor output line by line
            import time
            start_time = time.time()
            timeout_seconds = 120  # 2 minute overall timeout

            for line in iter(proc.stdout.readline, ''):
                if not line:
                    break

                output_lines.append(line.rstrip())

                # Check overall timeout
                if time.time() - start_time > timeout_seconds:
                    if console:
                        console.print("[red]✗ tLEaP info pass timed out[/red]")
                    proc.kill()
                    proc.wait()
                    break

                # Parse water count from solvate output
                if match := patterns['water_count'].search(line):
                    info['n_waters'] = int(match.group(1))
                    if console:
                        console.print(f"[grey50]  → Water count: {info['n_waters']}[/grey50]")

                # Parse charge from 'charge mol' output
                if match := patterns['charge'].search(line):
                    info['net_charge'] = float(match.group(1))
                    if console:
                        console.print(f"[grey50]  → Net charge: {info['net_charge']:+.1f}[/grey50]")

                # Parse volume (optional)
                if match := patterns['volume'].search(line):
                    info['volume'] = float(match.group(1))

                # Terminate early once we have required info
                if all(info.get(field) is not None for field in required_fields):
                    if console:
                        console.print("[green]  ✓ Got required info, terminating tLEaP early[/green]")
                    proc.terminate()
                    try:
                        proc.wait(timeout=2)
                    except subprocess.TimeoutExpired:
                        proc.kill()
                        proc.wait()
                    break

            # Clean up if still running
            if proc.poll() is None:
                proc.terminate()
                try:
                    proc.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()

            # If we didn't get required info, show diagnostic output
            if any(info.get(field) is None for field in required_fields):
                if console:
                    # Show tLEaP errors/warnings for debugging
                    error_lines = [l for l in output_lines
                                   if any(kw in l.upper() for kw in ['ERROR', 'FATAL', 'FAIL', 'WARNING', 'COULD NOT'])]
                    if error_lines:
                        console.print("[red]  tLEaP errors/warnings:[/red]")
                        for el in error_lines[:10]:
                            console.print(f"[red]    {el}[/red]")
                    else:
                        # No obvious errors — show last 15 lines for context
                        console.print("[yellow]  tLEaP output (last 15 lines):[/yellow]")
                        for el in output_lines[-15:]:
                            console.print(f"[grey50]    {el}[/grey50]")
                logger.warning("tLEaP info pass failed to extract required fields. "
                             "Captured %d output lines.", len(output_lines))

        except FileNotFoundError:
            if console:
                console.print("[red]tLEaP not found in PATH[/red]")
            return info
        except Exception as e:
            if console:
                from rich.markup import escape
                console.print(f"[red]Error in info pass: {escape(str(e))}[/red]")
            if proc is not None and proc.poll() is None:
                proc.kill()
                proc.wait()
        finally:
            # Clean up temp file — keep on failure for debugging
            if info_script and os.path.exists(info_script):
                if all(info.get(field) is not None for field in required_fields):
                    try:
                        os.remove(info_script)
                    except Exception:
                        pass
                else:
                    if console:
                        console.print(f"[grey50]  Info pass script saved for inspection: {info_script}[/grey50]")

        return info

    def _build_microstate_info_pass_script(self, microstate_info: dict, selected_forcefields: dict,
                                            water_box: str = "TIP3PBOX", buffer: float = 10.0,
                                            use_octahedron: bool = True, buffer_xyz=None,
                                            oct_diagonal: float = 0.0, iso: bool = False) -> str:
        """
        Build a minimal tLEaP script for info gathering for a specific microstate.

        Args:
            microstate_info: Microstate metadata dictionary
            selected_forcefields: Dictionary of selected forcefield files
            water_box: Water box type
            buffer: Buffer distance in Angstroms
            use_octahedron: If True, use solvateoct; else use solvateBox

        Returns:
            Path to the info pass script file
        """
        import tempfile

        # Build forcefield section from user selections
        forcefield_section, _ = self._build_standard_forcefield_section()

        # Build microstate-specific sections
        atom_types_section = self._build_microstate_atom_types_section(microstate_info, selected_forcefields)
        custom_ff_section = self._build_microstate_forcefield_section(microstate_info, selected_forcefields)
        pdb_section = self._build_microstate_pdb_section(microstate_info)
        bond_section = self._build_microstate_bond_section(microstate_info)

        # Solvation command
        solvate_cmd = "solvateoct" if use_octahedron else "solvateBox"
        buf_arg = _format_buffer_for_tleap(buffer, buffer_xyz, oct_diagonal,
                                            use_octahedron, iso)

        from pathlib import Path
        ms_pdb_stem = Path(microstate_info['filename']).stem
        script = f"""# INFO PASS for {ms_pdb_stem} - gathering water count and charge
# This script will be terminated early once we have the required information

# === STANDARD FORCEFIELDS ===
{forcefield_section}

# === CUSTOM ATOM TYPES ===
{atom_types_section}

# === CUSTOM FORCEFIELD PARAMETERS ===
{custom_ff_section}

# === LOAD STRUCTURE ===
{pdb_section}

# === BOND DEFINITIONS ===
{bond_section}

# === GET CHARGE BEFORE SOLVATION ===
charge mol

# === SOLVATION ===
{solvate_cmd} mol {water_box} {buf_arg}

# === GET CHARGE AFTER SOLVATION ===
charge mol

quit
"""

        # Write to temp file
        fd, temp_path = tempfile.mkstemp(suffix='_info_pass.in', prefix=f'tleap_{ms_pdb_stem}_')
        os.close(fd)
        with open(temp_path, 'w') as f:
            f.write(script)

        return temp_path

    def _run_info_pass_for_microstate(self, microstate_info: dict, selected_forcefields: dict,
                                       water_box: str = "TIP3PBOX", buffer: float = 10.0,
                                       use_octahedron: bool = True, buffer_xyz=None,
                                       oct_diagonal: float = 0.0, iso: bool = False) -> dict:
        """
        Run tLEaP info pass for a specific microstate to get water count and charge.

        Args:
            microstate_info: Microstate metadata dictionary
            selected_forcefields: Dictionary of selected forcefield files
            water_box: Water box type
            buffer: Buffer distance in Angstroms
            use_octahedron: If True, use solvateoct; else use solvateBox

        Returns:
            Dictionary with 'n_waters', 'net_charge' (values may be None)
        """
        import re

        # Build info pass script
        info_script = self._build_microstate_info_pass_script(
            microstate_info, selected_forcefields, water_box, buffer, use_octahedron,
            buffer_xyz=buffer_xyz, oct_diagonal=oct_diagonal, iso=iso
        )
        if not info_script:
            return {'n_waters': None, 'net_charge': None}

        info = {'n_waters': None, 'net_charge': None}
        required_fields = {'n_waters', 'net_charge'}

        # Patterns to match
        patterns = {
            'water_count': re.compile(r'Added\s+(\d+)\s+residues\.'),
            'charge': re.compile(r'Total unperturbed charge:\s*([-\d.]+)'),
        }

        proc = None
        try:
            proc = subprocess.Popen(
                ['tleap', '-f', info_script],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1
            )

            for line in iter(proc.stdout.readline, ''):
                if not line:
                    break

                if match := patterns['water_count'].search(line):
                    info['n_waters'] = int(match.group(1))

                if match := patterns['charge'].search(line):
                    info['net_charge'] = float(match.group(1))

                # Terminate early
                if all(info.get(field) is not None for field in required_fields):
                    proc.terminate()
                    try:
                        proc.wait(timeout=2)
                    except subprocess.TimeoutExpired:
                        proc.kill()
                        proc.wait()
                    break

            if proc.poll() is None:
                proc.terminate()
                try:
                    proc.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()

        except Exception:
            if proc is not None and proc.poll() is None:
                proc.kill()
                proc.wait()
        finally:
            if info_script and os.path.exists(info_script):
                try:
                    os.remove(info_script)
                except Exception:
                    pass

        return info

    def _execute_tleap_files_batch(self, selected_files: list, max_workers: int = 1) -> bool:
        """
        Execute tLEaP for multiple files in batch mode with compact output.

        Single-pass execution since templates now have accurate ion counts from info pass.
        When max_workers > 1, delegates to parallel execution with isolated directories.
        """
        if max_workers > 1:
            return self._execute_tleap_files_parallel(selected_files, max_workers)

        from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, MofNCompleteColumn
        import subprocess
        import re

        console = self.processor.console
        success_count = 0
        failed_files = []
        all_unique_warnings = set()
        all_unique_errors = set()
        parmed_reordered_count = 0

        console.print(f"\n[bold blue]Processing {len(selected_files)} tLEaP files in batch mode...[/bold blue]")
        console.print("[grey50]Templates already have accurate ion counts from template generation[/grey50]\n")

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            MofNCompleteColumn(),
            console=console
        ) as progress:
            task = progress.add_task("Generating topologies...", total=len(selected_files))

            for info in selected_files:
                tleap_file = info['tleap_file']
                # Extract microstate name from filename (e.g., "tleap_OOOR.in" -> "OOOR")
                base_name = os.path.basename(tleap_file)
                microstate_id = base_name.replace("tleap_", "").replace(".in", "").replace("_tleap", "")
                # Truncate long redox state strings (e.g., "OLOLOLOLOL..." -> "OLOL...OLOL")
                if len(microstate_id) > 12:
                    microstate_id = microstate_id[:5] + "..." + microstate_id[-4:]

                def update_step(step_name):
                    progress.update(task, description=f"[blue]{microstate_id}: {step_name}[/blue]")

                # Step 1: Validate PDB structure
                update_step("Validating PDB")
                molecule_config_success = self._check_and_configure_molecules(info, batch_mode=True)
                if not molecule_config_success:
                    failed_files.append((tleap_file, "Molecule configuration failed"))
                    progress.advance(task)
                    continue

                # Step 2: Check TER records
                update_step("Checking TER records")
                ter_ok = self._validate_and_fix_ter_records(info, quiet=True)
                if not ter_ok:
                    failed_files.append((tleap_file, "TER validation failed"))
                    progress.advance(task)
                    continue

                try:
                    # Single-pass tLEaP execution (templates have accurate ions)
                    update_step("Running tLEaP")
                    cmd = ["tleap", "-s", "-f", tleap_file]
                    result = subprocess.run(
                        cmd,
                        capture_output=True,
                        text=True
                    )

                    if result.returncode == 0:
                        # Parse warnings/errors from leap.log silently
                        warnings, errors, notes, summary = self._parse_leap_log("leap.log")

                        # Collect unique warning types (just the message text, not line numbers)
                        for warning_block in warnings:
                            if warning_block:
                                # Extract the key part of the warning
                                warning_text = ' '.join(warning_block).strip()
                                # tleap prefixes every message with the full
                                # "/path/to/teLeap: " program path (~50 chars);
                                # strip it so the actual message (e.g. the atom
                                # pair in a close-contact warning) isn't pushed
                                # past the truncation cap.
                                warning_text = re.sub(
                                    r'^(?:\S*/)?te?Leap:\s*', '', warning_text,
                                    flags=re.IGNORECASE)
                                # Simplify to get unique warning types
                                if "One sided connection" in warning_text:
                                    # Extract residue type
                                    match = re.search(r"Residue \((\w+)\)", warning_text)
                                    if match:
                                        all_unique_warnings.add(f"One sided connection: Residue ({match.group(1)})")
                                elif "addIons" in warning_text:
                                    all_unique_warnings.add("addIons: charges of the same sign")
                                else:
                                    # Truncate long warnings (cap raised so the
                                    # message survives after prefix stripping).
                                    short_warn = warning_text[:200]
                                    all_unique_warnings.add(short_warn)

                        for error_block in errors:
                            if error_block:
                                error_text = re.sub(
                                    r'^(?:\S*/)?te?Leap:\s*', '',
                                    ' '.join(error_block).strip(),
                                    flags=re.IGNORECASE)[:200]
                                all_unique_errors.add(error_text)

                        # Step 4: Validate topology with ParmEd
                        if info.get('expected_prmtop') and os.path.exists(info['expected_prmtop']):
                            update_step("Validating topology (ParmEd)")
                            rst7_file = info.get('expected_rst7') if info.get('expected_rst7') and os.path.exists(info.get('expected_rst7')) else None
                            reordered = self._run_parmed_validation_quiet(info['expected_prmtop'], rst7_file)
                            if reordered:
                                parmed_reordered_count += 1

                            self._register_topology_outputs(info)

                        success_count += 1
                    else:
                        failed_files.append((tleap_file, f"Exit code {result.returncode}"))

                except FileNotFoundError:
                    console.print(f"\n[red]✗ tLEaP command not found. Please ensure tLEaP is installed and in PATH.[/red]")
                    return False
                except Exception as e:
                    failed_files.append((tleap_file, str(e)))

                progress.advance(task)

        # Print summary
        console.print(f"\n[bold]{'═' * 60}[/bold]")
        console.print(f"[bold blue]Batch tLEaP Summary[/bold blue]")
        console.print(f"[bold]{'═' * 60}[/bold]")

        console.print(f"\n[green]✓ Successful:[/green] {success_count}/{len(selected_files)} topologies generated")

        if parmed_reordered_count > 0:
            console.print(f"[yellow]⚠ ParmEd reordered atoms in {parmed_reordered_count} files[/yellow] [grey50](molecules were not contiguous)[/grey50]")

        if failed_files:
            console.print(f"\n[red]✗ Failed ({len(failed_files)}):[/red]")
            for fname, reason in failed_files[:5]:
                console.print(f"  [red]• {os.path.basename(fname)}: {reason}[/red]")
            if len(failed_files) > 5:
                console.print(f"  [grey50]... and {len(failed_files) - 5} more[/grey50]")

        if all_unique_errors:
            console.print(f"\n[bold red]Unique Errors ({len(all_unique_errors)}):[/bold red]")
            for error in sorted(all_unique_errors):
                console.print(f"  [red]• {error}[/red]")

        if all_unique_warnings:
            console.print(f"\n[bold yellow]Unique Warning Types ({len(all_unique_warnings)}):[/bold yellow]")
            for warning in sorted(all_unique_warnings):
                console.print(f"  [yellow]• {warning}[/yellow]")
            console.print(f"[grey50]  (These warnings appeared across all microstates - typically expected for non-standard residues)[/grey50]")

        console.print()
        return success_count > 0

    def _execute_tleap_files_parallel(self, selected_files: list, max_workers: int) -> bool:
        """
        Execute tLEaP files in parallel using isolated working directories.

        Each tLEaP process runs in its own temp directory (via symlinks to the
        source directory) so that leap.log files don't conflict between workers.
        """
        from concurrent.futures import ThreadPoolExecutor, as_completed
        from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, MofNCompleteColumn
        import subprocess
        import re
        import tempfile
        import shutil

        console = self.processor.console
        source_dir = os.getcwd()
        success_count = 0
        failed_files = []
        all_unique_warnings = set()
        all_unique_errors = set()
        parmed_reordered_count = 0

        # Phase 1: Pre-validate all files sequentially
        console.print(f"\n[bold blue]Validating {len(selected_files)} tLEaP files...[/bold blue]")
        validated = []
        for info in selected_files:
            tleap_file = info['tleap_file']
            if not self._check_and_configure_molecules(info, batch_mode=True):
                failed_files.append((tleap_file, "Molecule configuration failed"))
                continue
            if not self._validate_and_fix_ter_records(info, quiet=True):
                failed_files.append((tleap_file, "TER validation failed"))
                continue
            validated.append(info)

        if not validated:
            console.print("[red]No files passed validation.[/red]")
            return False

        if failed_files:
            console.print(f"[yellow]⚠ {len(failed_files)} files skipped during validation[/yellow]")

        # Phase 2: Run tLEaP in parallel with isolated directories
        console.print(f"\n[bold blue]Running {len(validated)} tLEaP jobs ({max_workers} parallel workers)...[/bold blue]")
        console.print("[grey50]Templates already have accurate ion counts from template generation[/grey50]\n")

        def run_tleap_isolated(info):
            """Run a single tLEaP in an isolated temp directory."""
            tleap_file = info['tleap_file']
            temp_dir = tempfile.mkdtemp(prefix="tleap_parallel_")
            try:
                # Symlink all files from source directory into temp dir
                for item in os.listdir(source_dir):
                    src = os.path.join(source_dir, item)
                    dst = os.path.join(temp_dir, item)
                    try:
                        os.symlink(src, dst)
                    except OSError:
                        pass

                # Run tLEaP in isolated directory
                cmd = ["tleap", "-s", "-f", os.path.basename(tleap_file)]
                result = subprocess.run(cmd, capture_output=True, text=True, cwd=temp_dir)

                warnings_set = set()
                errors_set = set()

                if result.returncode == 0:
                    # Parse leap.log from temp directory
                    leap_log = os.path.join(temp_dir, "leap.log")
                    if os.path.exists(leap_log):
                        warnings, errors, notes, summary = self._parse_leap_log(leap_log)
                        for warning_block in warnings:
                            if warning_block:
                                warning_text = ' '.join(warning_block).strip()
                                if "One sided connection" in warning_text:
                                    match = re.search(r"Residue \((\w+)\)", warning_text)
                                    if match:
                                        warnings_set.add(f"One sided connection: Residue ({match.group(1)})")
                                elif "addIons" in warning_text:
                                    warnings_set.add("addIons: charges of the same sign")
                                else:
                                    warnings_set.add(warning_text[:100])
                        for error_block in errors:
                            if error_block:
                                errors_set.add(' '.join(error_block).strip()[:100])

                    # Copy output files back to source directory
                    for key in ['expected_prmtop', 'expected_rst7']:
                        if info.get(key):
                            basename = os.path.basename(info[key])
                            temp_out = os.path.join(temp_dir, basename)
                            dest = os.path.join(source_dir, basename)
                            if os.path.exists(temp_out):
                                shutil.copy2(temp_out, dest)

                    return (True, warnings_set, errors_set)
                else:
                    return (False, set(), {f"Exit code {result.returncode}"})

            except FileNotFoundError:
                return (False, set(), {"tLEaP command not found"})
            except Exception as e:
                return (False, set(), {str(e)})
            finally:
                shutil.rmtree(temp_dir, ignore_errors=True)

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            MofNCompleteColumn(),
            console=console
        ) as progress:
            task = progress.add_task("Generating topologies...", total=len(validated))

            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = {executor.submit(run_tleap_isolated, info): info for info in validated}
                for future in as_completed(futures):
                    info = futures[future]
                    try:
                        success, warns, errs = future.result()
                        if success:
                            success_count += 1
                            all_unique_warnings.update(warns)
                        else:
                            failed_files.append((info['tleap_file'], next(iter(errs), "Unknown error")))
                            all_unique_errors.update(errs)
                    except Exception as e:
                        failed_files.append((info['tleap_file'], str(e)))
                    progress.advance(task)

        # Phase 3: ParmEd validation (sequential - fast, accesses shared workspace)
        console.print(f"[grey50]Running ParmEd validation...[/grey50]")
        for info in validated:
            if info.get('expected_prmtop'):
                prmtop_path = os.path.join(source_dir, os.path.basename(info['expected_prmtop']))
                rst7_path = None
                if info.get('expected_rst7'):
                    rst7_path = os.path.join(source_dir, os.path.basename(info['expected_rst7']))
                    if not os.path.exists(rst7_path):
                        rst7_path = None
                if os.path.exists(prmtop_path):
                    reordered = self._run_parmed_validation_quiet(prmtop_path, rst7_path)
                    if reordered:
                        parmed_reordered_count += 1
                    self._register_topology_outputs(info)

        # Phase 4: Summary
        console.print(f"\n[bold]{'═' * 60}[/bold]")
        console.print(f"[bold blue]Parallel tLEaP Summary[/bold blue]")
        console.print(f"[bold]{'═' * 60}[/bold]")

        console.print(f"\n[green]✓ Successful:[/green] {success_count}/{len(selected_files)} topologies generated")

        if parmed_reordered_count > 0:
            console.print(f"[yellow]⚠ ParmEd reordered atoms in {parmed_reordered_count} files[/yellow] [grey50](molecules were not contiguous)[/grey50]")

        if failed_files:
            console.print(f"\n[red]✗ Failed ({len(failed_files)}):[/red]")
            for fname, reason in failed_files[:5]:
                console.print(f"  [red]• {os.path.basename(fname)}: {reason}[/red]")
            if len(failed_files) > 5:
                console.print(f"  [grey50]... and {len(failed_files) - 5} more[/grey50]")

        if all_unique_errors:
            console.print(f"\n[bold red]Unique Errors ({len(all_unique_errors)}):[/bold red]")
            for error in sorted(all_unique_errors):
                console.print(f"  [red]• {error}[/red]")

        if all_unique_warnings:
            console.print(f"\n[bold yellow]Unique Warning Types ({len(all_unique_warnings)}):[/bold yellow]")
            for warning in sorted(all_unique_warnings):
                console.print(f"  [yellow]• {warning}[/yellow]")
            console.print(f"[grey50]  (These warnings appeared across all microstates - typically expected for non-standard residues)[/grey50]")

        console.print()
        return success_count > 0

    def _get_12_6_4_params(self) -> dict:
        """
        Check if the user selected a 12-6-4 ion model requiring ParmEd add12_6_4.

        The 12-6-4 ion parameters (Li & Merz) were developed with ff14SB, but the
        correction transfers to ff19SB: the AmberTools polarizability file
        (lj_1264_pol.dat) ships the ff19SB-renamed types (e.g. XC) alongside the
        ff14SB names (CX) with identical polarizabilities, and ProPrep's
        _build_augmented_polfile infers any still-missing type from element and
        bond count as a backstop. Verified empirically: params1264 with OPC water
        applies cleanly to an ff19SB topology.

        ff15ipq remains gated, but NOT for a technical reason: params1264 runs to
        completion once _build_augmented_polfile supplies the novel IPQ atom types
        (TA/TH/TJ/TM/TP), which it infers correctly from their shared LJ slots.
        The gate is a modeling decision. The Li/Merz 12-6-4 sets were calibrated
        against ion interactions built on ff14SB-family point charges, whereas
        ff15ipq uses the implicitly-polarized IPQ charge model with systematically
        different electrostatics. Layering an ff14SB-calibrated ion correction on
        top of IPQ charges is not a validated combination, so ProPrep declines it
        by default rather than silently produce an unvalidated model. A user who
        has validated the pairing can apply the correction manually.

        Returns:
            dict with params if 12-6-4 should be applied, or None if not needed.
            If the forcefield is incompatible, returns a dict with 'incompatible': True
            instead of 'needed': True.
        """
        selected = self.get_from_workspace("selected_standard_forcefields", {})
        ion_sel = selected.get('ions')
        water_sel = selected.get('water')
        protein_sel = selected.get('protein')

        if not ion_sel or not isinstance(ion_sel, dict):
            return None

        # Check if frcmod contains '1264' (the 12-6-4 indicator)
        frcmod = ion_sel.get('frcmod', '') or ''
        # frcmod can be a string or list of strings
        frcmod_str = ' '.join(frcmod) if isinstance(frcmod, list) else frcmod
        if '1264' not in frcmod_str:
            return None

        # Check protein forcefield compatibility. Only ff15ipq is gated, and by
        # modeling decision rather than technical limitation: the Li/Merz 12-6-4
        # sets were calibrated on ff14SB-family point charges, while ff15ipq uses
        # the implicitly-polarized IPQ charge model, so the combination is not
        # validated (see docstring). ff19SB is NOT gated — it stays in the ff14SB
        # electrostatic family and params1264 applies cleanly (confirmed).
        incompatible_ffs = {'ff15ipq', 'ff15IPQ'}
        protein_name = protein_sel.get('name', '') if isinstance(protein_sel, dict) else ''
        if protein_name in incompatible_ffs:
            return {
                'incompatible': True,
                'protein_ff': protein_name,
                'ion_selection_name': ion_sel.get('name', '12-6-4'),
            }

        # Map ProPrep water model names to ParmEd's expected names
        water_model_map = {
            'OPC': 'OPC',
            'OPC3': 'OPC3',
            'OPC3-pol': 'OPC3',
            'TIP3P': 'TIP3P',
            'TIP4P-Ew': 'TIP4PEW',
            'SPC/E': 'SPCE',
            'SPC/Eb': 'SPCE',
            'TIP3P-FB': 'FB3',
            'TIP4P-FB': 'FB4',
            'TIP5P': 'TIP3P',  # Fallback; TIP5P not supported by add12_6_4
        }

        water_name = water_sel.get('name', 'OPC') if isinstance(water_sel, dict) else 'OPC'
        parmed_water = water_model_map.get(water_name, 'OPC')

        return {
            'needed': True,
            'water_model': parmed_water,
            'ion_selection_name': ion_sel.get('name', '12-6-4'),
        }

    def _detect_ion_mask(self, parm, console=None) -> str:
        """
        Dynamically build an AMBER mask covering FREE ions in the topology.

        Scans the prmtop for single-atom residues with atomic number > 2
        (excludes H and He). Skips bonded-model metals: any metal atom with
        explicit bond partners is part of a parameterized metal site
        (MCPB-style Zn(Cys)4, heme Fe-coordination shell, etc.) and already
        has custom LJ values from its frcmod. Such atoms must not be passed
        to ParmEd's add12_6_4, which is calibrated for free divalent ions in
        their +2 state and crashes on fractional metal charges (e.g., the
        Zn(Cys)4 Zn carries +0.96 e after S→Zn donation, so add1264 tries
        to look up 'Zn0' in its 12-6-4 parameter table and KeyErrors).

        Args:
            parm: ParmEd AmberParm object
            console: Optional Rich console for output

        Returns:
            AMBER mask string (e.g., ':Ca2+,Na+,Cl-') or None if no ions found
        """
        ion_names = set()
        for residue in parm.residues:
            if len(residue.atoms) != 1:
                continue
            atom = residue.atoms[0]
            if atom.atomic_number <= 2:
                continue
            # Bonded-model metals have explicit bonds to their ligand atoms;
            # free solvent ions (Na+/K+/Cl-/free Mg2+/Ca2+/Zn2+) have none.
            # This catches Zn(Cys)4 (ZM type, bonded to 4 SG), heme Fe (FO/FR,
            # bonded to 4 pyrrole Ns + axial S), MCPB metals (M1/M2/...,
            # bonded to whatever the site's ligands are), and any future
            # bonded-model metal added to the bundle.
            if len(atom.bond_partners) > 0:
                if console:
                    n = len(atom.bond_partners)
                    console.print(
                        f"[grey50]  Skipping bonded-model metal {residue.name} "
                        f"(type {atom.type.strip()}, {n} explicit bond(s)) — "
                        f"custom LJ already from frcmod; 12-6-4 not applicable[/grey50]"
                    )
                continue
            ion_names.add(residue.name)

        if not ion_names:
            if console:
                console.print("[grey50]  No free ions detected — skipping add12_6_4[/grey50]")
            return None

        mask = ':' + ','.join(sorted(ion_names))
        if console:
            console.print(f"[grey50]  Detected free ions: {', '.join(sorted(ion_names))}[/grey50]")
        return mask

    def _build_augmented_polfile(self, parm, console=None) -> str:
        """
        Build an augmented polarizability file for ParmEd's add12_6_4 action.

        The default lj_1264_pol.dat only covers standard AMBER atom types.
        Non-standard types from modified amino acids, custom parameterizations,
        etc. will cause add12_6_4 to crash. This method:
          1. Reads the default pol file
          2. Identifies atom types in the prmtop that are missing
          3. Infers polarizabilities from element and hybridization (bond count)
          4. Writes an augmented file to the working directory

        Polarizabilities are from Miller, JACS 112, 8533 (1990), the same source
        used by the original pol file. See also:
        http://archive.ambermd.org/202303/0095.html

        Args:
            parm: ParmEd AmberParm object
            console: Optional Rich console for output

        Returns:
            Path to the polarizability file (augmented or default)
        """
        import sys

        # Locate default pol file
        default_polfile = None
        amberhome = os.environ.get('AMBERHOME', '')
        if amberhome:
            candidate = os.path.join(amberhome, 'dat', 'leap', 'parm', 'lj_1264_pol.dat')
            if os.path.exists(candidate):
                default_polfile = candidate
        if not default_polfile:
            candidate = os.path.join(sys.prefix, 'dat', 'leap', 'parm', 'lj_1264_pol.dat')
            if os.path.exists(candidate):
                default_polfile = candidate
        if not default_polfile:
            raise FileNotFoundError(
                "Could not find lj_1264_pol.dat. Set $AMBERHOME or ensure "
                "AmberTools is installed in the active conda environment."
            )

        # Parse default pol file: {atom_type: polarizability}
        known_types = {}
        default_lines = []
        with open(default_polfile, 'r') as f:
            for line in f:
                default_lines.append(line)
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        known_types[parts[0]] = float(parts[1])
                    except ValueError:
                        pass

        # Find all unique atom types in the prmtop
        prmtop_types = set(parm.parm_data['AMBER_ATOM_TYPE'])
        typs = parm.parm_data['AMBER_ATOM_TYPE']
        typinds = parm.parm_data['ATOM_TYPE_INDEX']

        # ─── Polfile/prmtop reconciliation pass ────────────────────────
        # ParmEd's params1264 indexes the C4 matrix by ATOM_TYPE_INDEX (the
        # prmtop's LJ slot), not by AMBER_ATOM_TYPE name. It requires every
        # type name in a given LJ slot to report the same α. Bundles that
        # add a custom atom-type name without shipping a separate NONBON
        # entry (e.g. heme's Cp/VB inheriting CA's LJ via equivalence)
        # collapse multiple names to one slot — and if the default polfile
        # has entries for BOTH names with different α, parmed crashes with
        # "Polarizability parameter of X is not the same as Y, but their
        # VDW parameters are the same".
        #
        # Detect that case: for each LJ slot, if >1 known type maps to it
        # AND those types have different α in the default polfile, pick the
        # most-populated type's α as canonical (because that name is what
        # the slot actually represents in this prmtop) and override the
        # losing types' α to match. Stale legacy polfile entries (e.g.
        # ancient "Cp 1.061" for some unrelated pre-fork meaning) get
        # rewritten to the in-use slot's α.
        import collections
        slot_to_known_pols = {}
        for at, ti in zip(typs, typinds):
            if at in known_types:
                slot_to_known_pols.setdefault(ti, {})[at] = known_types[at]
        type_counts = collections.Counter(typs)
        polfile_overrides = {}  # type_name -> (new_pol, canonical_type, slot)
        for ti, pols in slot_to_known_pols.items():
            if len(set(pols.values())) <= 1:
                continue
            canonical_type = max(pols, key=lambda t: type_counts.get(t, 0))
            canonical_pol = pols[canonical_type]
            for t, p in pols.items():
                if t == canonical_type or p == canonical_pol:
                    continue
                polfile_overrides[t] = (canonical_pol, canonical_type, ti)
                known_types[t] = canonical_pol
        if polfile_overrides and console:
            for t, (newp, ref, ti) in sorted(polfile_overrides.items()):
                console.print(
                    f"[yellow]  ⚠ Polfile override: {t} α={newp:.3f} "
                    f"(was inconsistent with LJ-slot-{ti} partner {ref}; "
                    f"matched to canonical α)[/yellow]"
                )
        # ────────────────────────────────────────────────────────────────

        # Find missing types (those still absent from known_types after
        # the reconciliation pass)
        missing_types = prmtop_types - set(known_types.keys())

        if not missing_types and not polfile_overrides:
            # Default polfile is fine as-is
            return default_polfile

        # Build index: atom_type -> representative atom (for element/bond info)
        type_to_atom = {}
        for atom in parm.atoms:
            atype = parm.parm_data['AMBER_ATOM_TYPE'][atom.idx]
            if atype in missing_types and atype not in type_to_atom:
                type_to_atom[atype] = atom
            if len(type_to_atom) == len(missing_types):
                break

        # ParmEd's params1264 requires all AMBER_ATOM_TYPEs sharing an
        # ATOM_TYPE_INDEX to have the same polarizability. Build a map from
        # type index → known polarizability so missing types that share a
        # VDW type with a known type get the correct value.
        typind_to_pol = {}
        for at, ti in zip(typs, typinds):
            if at in known_types and ti not in typind_to_pol:
                typind_to_pol[ti] = (known_types[at], at)
        missing_typind = {}
        for at, ti in zip(typs, typinds):
            if at in missing_types and at not in missing_typind:
                missing_typind[at] = ti

        # Infer polarizabilities using element + hybridization
        # Values from Miller, JACS 112, 8533 (1990)
        new_entries = []
        for atype in sorted(missing_types):
            ti = missing_typind.get(atype)
            if ti is not None and ti in typind_to_pol:
                # Must match the known type sharing this VDW type index
                pol, ref = typind_to_pol[ti]
                comment = f"matched to {ref} (shared ATOM_TYPE_INDEX {ti})"
            else:
                atom = type_to_atom.get(atype)
                if atom is None:
                    pol = 0.000
                    comment = "Unknown (no representative atom found)"
                else:
                    z = atom.atomic_number
                    n_bonds = len(atom.bond_partners)
                    pol, comment = self._infer_polarizability(z, n_bonds, atom)
                    comment = f"element={atom.element_name}, bonds={n_bonds}, {comment}"
                if ti is not None:
                    typind_to_pol[ti] = (pol, atype)

            new_entries.append((atype, pol, comment))

            if console:
                console.print(
                    f"[yellow]  ⚠ Inferred polarizability for atom type "
                    f"{atype}: {pol:.3f} ({comment})[/yellow]"
                )
            logger.debug("Inferred polarizability for atom type %s: %.3f (%s)",
                        atype, pol, comment)

        # Write augmented file to working directory
        augmented_path = os.path.join(os.getcwd(), 'lj_1264_pol_augmented.dat')
        with open(augmented_path, 'w') as f:
            # Copy default file contents, rewriting α on any line whose type
            # was overridden by the reconciliation pass above.
            for line in default_lines:
                parts = line.split()
                if len(parts) >= 2 and parts[0] in polfile_overrides:
                    new_pol, ref, ti = polfile_overrides[parts[0]]
                    # Preserve the original tail comment if any (everything
                    # after the 2nd whitespace-separated token).
                    after_pol = line.split(None, 2)
                    tail_orig = after_pol[2].rstrip("\n") if len(after_pol) >= 3 else ""
                    override_note = (
                        f"ProPrep override (was {parts[1]}, matched to {ref} via LJ slot {ti})"
                    )
                    if tail_orig:
                        tail = f"    {tail_orig}    {override_note}"
                    else:
                        tail = f"    {override_note}"
                    f.write(f"{parts[0]:<6s}{new_pol:.3f}{tail}\n")
                else:
                    f.write(line)
            # Ensure trailing newline before appended entries
            if default_lines and not default_lines[-1].endswith('\n'):
                f.write('\n')
            # Append inferred entries (no comment-only lines — ParmEd's parser
            # doesn't handle them; comments go after the two data columns)
            for atype, pol, comment in new_entries:
                f.write(f"{atype:<6s}{pol:.3f}    Inferred by ProPrep ({comment})\n")

        if console:
            console.print(f"[grey50]  Augmented polarizability file: {os.path.basename(augmented_path)}[/grey50]")

        return augmented_path

    @staticmethod
    def _infer_polarizability(atomic_number: int, n_bonds: int, atom=None) -> tuple:
        """
        Infer atomic polarizability from element and bond count (hybridization).

        Values from Miller, JACS 112, 8533 (1990) — the same source used by
        AMBER's lj_1264_pol.dat.

        Args:
            atomic_number: Atomic number of the element
            n_bonds: Number of bonded partners
            atom: Optional ParmEd atom for additional context

        Returns:
            Tuple of (polarizability: float, comment: str)
        """
        # Carbon: hybridization-dependent
        if atomic_number == 6:
            if n_bonds >= 4:
                return 1.061, "sp3 carbon"
            elif n_bonds == 3:
                return 1.352, "sp2 carbon"
            elif n_bonds == 2:
                return 1.283, "sp carbon"
            else:
                return 1.352, "carbon (assumed sp2)"

        # Hydrogen
        if atomic_number == 1:
            # Check if bonded to oxygen (HO convention: pol=0.000)
            if atom is not None and n_bonds >= 1:
                for partner in atom.bond_partners:
                    if partner.atomic_number == 8:
                        return 0.000, "hydrogen on oxygen"
            return 0.387, "hydrogen"

        # Nitrogen
        if atomic_number == 7:
            return 1.090, "nitrogen"

        # Oxygen
        if atomic_number == 8:
            if n_bonds <= 1:
                return 0.569, "carbonyl/terminal oxygen"
            else:
                return 0.637, "bridging/hydroxyl oxygen"

        # Sulfur
        if atomic_number == 16:
            return 3.000, "sulfur"

        # Phosphorus
        if atomic_number == 15:
            return 1.538, "phosphorus"

        # Halogens
        halogen_map = {9: (0.32, "fluorine"), 17: (1.91, "chlorine"),
                       35: (2.88, "bromine"), 53: (4.69, "iodine")}
        if atomic_number in halogen_map:
            return halogen_map[atomic_number]

        # Metals and other elements: 0.000 is safe
        # (ions already have entries in the default pol file)
        return 0.000, "default (no Miller data for this element)"

    def _write_parmed_replay_script(
        self,
        prmtop_file: str,
        rst7_file: str = None,
        ion_mask: str = None,
        water_model: str = None,
        polfile: str = None,
        console=None,
    ) -> str:
        """Write a Python script that replays the parmed validation pass.

        ProPrep runs ParmEd's check_validity / rediscover_molecules /
        add12_6_4 in-process via the Python API. Users can't re-run that
        without re-running ProPrep. This emits a standalone Python script
        next to the prmtop that uses the parmed Python API to reproduce
        the same operations — easier to edit (tweak mask, water model,
        polfile, output filename) than a CLI .in script, and avoids
        parmed CLI's argument-parsing quirks (rediscover_molecules is
        Python-API-only; add12_6_4 CLI uses ``key value`` not ``-flag``).

        Args:
            prmtop_file: Path to the topology the script targets.
            rst7_file: Optional coords; loaded alongside the prmtop.
            ion_mask: AMBER mask for add12_6_4. None → skip 12-6-4.
            water_model: Water model name for add12_6_4 (e.g., 'TIP3P').
            polfile: Path to the augmented polarizability file. None →
                let parmed use its built-in default.
            console: Optional Rich console.

        Returns:
            Path to the written script.
        """
        from pathlib import Path
        p = Path(prmtop_file)
        script_path = p.with_name(p.stem + "_parmed_replay.py")
        rst7_name = Path(rst7_file).name if rst7_file else None
        polfile_name = Path(polfile).name if polfile else None
        amberhome = os.environ.get("AMBERHOME", "/path/to/amber")

        # Build the script body. Python triple-quoted with explicit Python
        # source — emit raw so the user can edit any of the configuration
        # constants up top.
        body = f'''#!/usr/bin/env python3
"""Auto-generated by ProPrep — replays the parmed validation + 12-6-4 pass
ProPrep ran in-process on this topology, using parmed's Python API.

Run:
    python {script_path.name}

Or with conda activated and AMBERHOME exported, just:
    ./{script_path.name}

Edit the configuration block below to tweak any knob — the mask, water
model, polfile, output filenames, etc. The default is to overwrite the
input prmtop (and rst7) in place; change PRMTOP_OUT / RST7_OUT to write
to new filenames if you want to keep the originals.
"""
import os
import sys
import warnings as _w

# ─── Configuration (edit anything below as needed) ──────────────────────
os.environ.setdefault("AMBERHOME", "{amberhome}")

PRMTOP_IN   = "{p.name}"
RST7_IN     = {repr(rst7_name)}      # None if no coordinates
PRMTOP_OUT  = PRMTOP_IN               # overwrite in place; change to keep original
RST7_OUT    = RST7_IN                 # ditto

ION_MASK    = {repr(ion_mask)}     # AmberMask covering free divalent+ ions; None → skip add12_6_4
WATER_MODEL = {repr(water_model)}     # one of TIP3P / TIP4PEW / SPCE / OPC3 / OPC / FB3 / FB4
POLFILE     = {repr(polfile_name)}    # augmented polfile (None → parmed default)
# ────────────────────────────────────────────────────────────────────────

import parmed as pmd
from parmed.tools.checkvalidity import check_validity
from parmed.exceptions import ParmedWarning

# 1. Load
print(f"Loading {{PRMTOP_IN}}...")
if RST7_IN and os.path.exists(RST7_IN):
    parm = pmd.load_file(PRMTOP_IN, xyz=RST7_IN)
else:
    parm = pmd.load_file(PRMTOP_IN)
print(f"  {{len(parm.atoms)}} atoms, {{len(parm.residues)}} residues")

# 2. Validity check
print("Running check_validity()...")
with _w.catch_warnings(record=True) as wlist:
    _w.simplefilter("ignore")
    _w.simplefilter("always", category=ParmedWarning)
    check_validity(parm, _w)
    if wlist:
        print(f"  {{len(wlist)}} warning(s):")
        for warning in wlist[:5]:
            print(f"    • {{warning.message}}")
        if len(wlist) > 5:
            print(f"    ... and {{len(wlist) - 5}} more")
    else:
        print("  ✓ no warnings")

# 3. Rediscover molecules (fixes ATOMS_PER_MOLECULE when tleap left atoms
#    of the same molecule non-contiguous — common after explicit `bond`
#    directives across residues). Python-only; not available via parmed CLI.
print("Rediscovering molecules...")
atom_reorder = parm.rediscover_molecules(solute_ions=True, fix_broken=True)
if atom_reorder is not None:
    print("  ⚠ Atoms were reordered; rebuilding parm_data arrays")
    parm._xfer_atom_info()
else:
    print("  ✓ molecule definitions already correct")

# 4. 12-6-4 LJ correction for divalent+ free ions (Li/Merz parameters).
if ION_MASK:
    from parmed.tools.add1264 import params1264
    from parmed.amber import AmberMask

    print(f"Applying 12-6-4 (mask={{ION_MASK}}, water={{WATER_MODEL}})...")
    mask = AmberMask(parm, ION_MASK)
    if "LENNARD_JONES_CCOEF" not in parm.flag_list:
        parm.add_flag(
            "LENNARD_JONES_CCOEF", "5E16.8",
            num_items=len(parm.parm_data["LENNARD_JONES_ACOEF"]),
            comments=["For 12-6-4 potential used for ions"],
        )
    c4_terms = params1264(parm, mask, None, WATER_MODEL, POLFILE, 1.0)
    for i, param in enumerate(c4_terms):
        parm.parm_data["LENNARD_JONES_CCOEF"][i] = param
    print(f"  ✓ added LENNARD_JONES_CCOEF (C4 terms)")
else:
    print("Skipping add12_6_4 (no free divalent+ ions in topology)")

# 5. Save
print(f"Saving to {{PRMTOP_OUT}}...")
parm.save(PRMTOP_OUT, overwrite=True)
if RST7_OUT and parm.coordinates is not None:
    parm.save(RST7_OUT, overwrite=True)
    print(f"  ✓ wrote prmtop + rst7")
else:
    print(f"  ✓ wrote prmtop")
'''
        script_path.write_text(body)
        # Make it executable
        try:
            import stat
            current = script_path.stat().st_mode
            script_path.chmod(current | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        except OSError:
            pass
        if console:
            console.print(
                f"[grey50]  Wrote parmed replay script: {script_path.name}[/grey50]"
            )
        return str(script_path)

    def _apply_12_6_4(self, parm, params: dict, console=None,
                       prmtop_file: str = None, rst7_file: str = None) -> bool:
        """
        Apply ParmEd add12_6_4 action to a loaded topology.

        Dynamically detects ions in the system and builds an augmented
        polarizability file if needed for non-standard atom types.

        Args:
            parm: ParmEd AmberParm object
            params: dict from _get_12_6_4_params() (must have 'water_model')
            console: Optional Rich console for output
            prmtop_file: Path to the topology on disk — used to emit the
                parmed replay script alongside it.
            rst7_file: Optional coords path for the replay script.

        Returns:
            bool: True if successfully applied
        """
        try:
            from parmed.tools.add1264 import params1264
            from parmed.amber import AmberMask

            water_model = params['water_model']

            if console:
                console.print(f"[blue]  Applying 12-6-4 LJ potential (water model: {water_model})...[/blue]")

            # Detect ions dynamically from the topology
            ion_mask_str = self._detect_ion_mask(parm, console)
            if ion_mask_str is None:
                if console:
                    console.print("[yellow]  ⚠ No ions found in topology — skipping add12_6_4[/yellow]")
                # Still write a replay script (just check+rediscover+outparm)
                if prmtop_file:
                    self._write_parmed_replay_script(
                        prmtop_file, rst7_file, None, None, None, console
                    )
                return True

            # Build augmented polfile if needed for non-standard atom types
            polfile = self._build_augmented_polfile(parm, console)

            # Emit the parmed CLI replay script BEFORE calling params1264 so
            # users have the script even if the in-process call crashes
            # (which is exactly the situation that motivated adding this).
            if prmtop_file:
                self._write_parmed_replay_script(
                    prmtop_file, rst7_file, ion_mask_str, water_model, polfile, console
                )

            # Call params1264 directly to avoid ParmEd's arg parser issues
            # (the arg parser can't reliably separate mask from polfile path)
            mask = AmberMask(parm, ion_mask_str)

            # Add LENNARD_JONES_CCOEF flag if not present
            if 'LENNARD_JONES_CCOEF' not in parm.flag_list:
                parm.add_flag(
                    'LENNARD_JONES_CCOEF', '5E16.8',
                    num_items=len(parm.parm_data['LENNARD_JONES_ACOEF']),
                    comments=['For 12-6-4 potential used for ions']
                )

            c4_terms = params1264(parm, mask, None, water_model, polfile, 1.0)
            for i, param in enumerate(c4_terms):
                parm.parm_data['LENNARD_JONES_CCOEF'][i] = param

            if console:
                console.print(f"[green]  ✓ Added LENNARD_JONES_CCOEF (C4 terms) for 12-6-4 potential[/green]")

            return True

        except Exception as e:
            msg = f"add12_6_4 failed: {e}"
            if console:
                from rich.markup import escape
                console.print(f"[red]  ✗ {escape(msg)}[/red]")
                console.print(f"[yellow]    You may need to run ParmEd manually with a custom polfile[/yellow]")
            logger.warning(msg)
            return False

    def _run_parmed_validation_quiet(self, prmtop_file: str, rst7_file: str = None) -> bool:
        """
        Run ParmEd validation silently, returning True if atoms were reordered.

        Args:
            prmtop_file: Path to the AMBER topology file (.prmtop)
            rst7_file: Optional path to the coordinate file (.rst7)

        Returns:
            bool: True if atoms were reordered, False otherwise
        """
        import warnings as py_warnings

        if not os.path.exists(prmtop_file):
            return False

        try:
            import parmed as pmd

            if rst7_file and os.path.exists(rst7_file):
                parm = pmd.load_file(prmtop_file, xyz=rst7_file)
            else:
                parm = pmd.load_file(prmtop_file)

            # Run rediscover_molecules to fix any ordering issues
            atoms_reordered = False
            atom_reorder = parm.rediscover_molecules(solute_ions=True, fix_broken=True)
            if atom_reorder is not None:
                atoms_reordered = True
                # rediscover_molecules reorders parm.atoms but does NOT rebuild
                # parm.parm_data arrays (CHARGE, AMBER_ATOM_TYPE, etc.).
                # Rebuild so downstream code (e.g., params1264) sees correct values.
                parm._xfer_atom_info()

            # Apply 12-6-4 C4 terms if user selected 12-6-4 ion parameters
            lj_params = self._get_12_6_4_params()
            if lj_params and lj_params.get('needed'):
                self._apply_12_6_4(parm, lj_params)
            elif lj_params and lj_params.get('incompatible'):
                logger.warning("Skipping 12-6-4: %s is incompatible with %s",
                             lj_params.get('ion_selection_name'),
                             lj_params.get('protein_ff'))

            # Save corrected topology
            parm.save(prmtop_file, overwrite=True)
            if rst7_file and parm.coordinates is not None:
                parm.save(rst7_file, overwrite=True)

            return atoms_reordered

        except Exception:
            return False

    def _validate_generated_topology_files(self, selected_files: list):
        """Validate that expected prmtop and rst7 files were generated."""
        console = self.processor.console
        console.print("\n[blue]Validating generated topology files...[/blue]")
        
        validation_table = Table(title="Generated File Validation")
        validation_table.add_column("tLEaP File", style="blue")
        validation_table.add_column("Expected prmtop", style="yellow")
        validation_table.add_column("prmtop Status", style="green")
        validation_table.add_column("Expected rst7", style="yellow")
        validation_table.add_column("rst7 Status", style="green")
        
        all_success = True
        
        for info in selected_files:
            prmtop_status = "✓ Found" if os.path.exists(info['expected_prmtop']) else "✗ Missing"
            rst7_status = "✓ Found" if os.path.exists(info['expected_rst7']) else "✗ Missing"
            
            prmtop_color = "green" if "Found" in prmtop_status else "red"
            rst7_color = "green" if "Found" in rst7_status else "red"
            
            if "Missing" in prmtop_status or "Missing" in rst7_status:
                all_success = False
            
            validation_table.add_row(
                info['tleap_file'],
                info['expected_prmtop'] or "N/A",
                f"[{prmtop_color}]{prmtop_status}[/{prmtop_color}]",
                info['expected_rst7'] or "N/A", 
                f"[{rst7_color}]{rst7_status}[/{rst7_color}]"
            )
        
        console.print(validation_table)
        
        if all_success:
            console.print("[bold green]All expected topology files were generated successfully![/bold green]")
        else:
            console.print("[bold yellow]Some expected files are missing. Check tLEaP output for errors.[/bold yellow]")

        return all_success

    def _print_parmed_info_panel(self, prmtop_file: str, console) -> None:
        """Show a short orientation panel before the ParmEd validation pass.

        The pass produces a lot of output — especially the polarizability
        inference cascade when bundled cofactor force fields introduce many
        forked atom types — and new users can mistake the informational
        messages for errors. This panel frames what's about to happen so the
        warnings that follow read as expected bookkeeping, not failures.
        """
        lj_params = self._get_12_6_4_params() or {}
        will_apply_1264 = bool(lj_params.get('needed'))
        incompatible_1264 = bool(lj_params.get('incompatible'))

        # Build the body as Rich-renderable lines. Each bullet is one logical
        # paragraph — Rich handles the wrapping. Blank strings produce blank
        # lines in the panel.
        body = [
            "[bold]ParmEd is Amber's topology editor.[/bold] ProPrep runs it here as a final pass over"
            " the prmtop tLEaP just wrote. It does not change the force field — only fixes topology"
            " bookkeeping and (optionally) adds the Li/Merz 12-6-4 ion correction.",
            "",
            "[bold blue]1. Validity check[/bold blue] — scans the prmtop for known structural issues."
            " Informational warnings (e.g. close Cys–Cys sulfur pairs that may want a disulfide bond)"
            " are flagged but do not stop the run.",
            "",
            "[bold blue]2. Rediscover molecules[/bold blue] — recomputes the ATOMS_PER_MOLECULE list"
            " from the actual bond graph. Required when tLEaP's grouping is stale (e.g. after custom"
            " bonds were added). If atoms are reordered, the topology arrays are rebuilt accordingly.",
            "",
        ]

        if will_apply_1264:
            body += [
                "[bold blue]3. Apply 12-6-4 LJ correction[/bold blue] — installs the Li/Merz C4 ion"
                " terms for free metal ions. Bonded-model metals are skipped because their LJ already"
                " comes from the bundled frcmod.",
                "",
                "ParmEd indexes the C4 matrix by LJ slot, so every atom-type name sharing a slot must"
                " report the same polarizability (α). ProPrep builds an augmented polfile for this"
                " prmtop with two kinds of entries you'll see scroll past:",
                "",
                "  • [bold]Polfile override[/bold] — resolves α conflicts where a forked type"
                " (e.g. heme [italic]Cp[/italic]) and its parent (e.g. [italic]CA[/italic])"
                " collapsed into one LJ slot in the prmtop.",
                "  • [bold]Inferred polarizability[/bold] — supplies α for atom-type forks the default"
                " polfile doesn't know about. Derived from the parent type sharing the LJ slot, or"
                " from element + bond count when no parent is available.",
                "",
                "[grey50]Both are informational, not errors.[/grey50]",
            ]
        elif incompatible_1264:
            body += [
                "[bold blue]3. 12-6-4 LJ correction[/bold blue] — [yellow]skipped[/yellow]:"
                " the selected ion parameter set is incompatible with the chosen protein force field."
                " Standard 12-6 LJ is retained.",
            ]
        else:
            body += [
                "[bold blue]3. 12-6-4 LJ correction[/bold blue] — not requested for this system."
                " Standard 12-6 LJ is retained.",
            ]

        body += [
            "",
            "When the pass finishes, ProPrep overwrites [italic]{0}[/italic] with the corrected"
            " topology and writes a [italic]<stem>_parmed_replay.py[/italic] script next to it so"
            " the same pass can be rerun by hand.".format(os.path.basename(prmtop_file)),
        ]

        console.print(Panel(
            "\n".join(body),
            title="[bold blue]About to run: ParmEd validation[/bold blue]",
            border_style="blue",
            padding=(1, 2),
            width=min(100, console.width),
        ))

    def _run_parmed_validation(self, prmtop_file: str, rst7_file: str = None) -> bool:
        """
        Run ParmEd validation on generated topology files.

        This method:
        1. Loads the topology file with ParmEd
        2. Runs check_validity() to identify issues
        3. Runs rediscover_molecules() to fix ATOMS_PER_MOLECULE section
        4. Saves the corrected topology file

        Args:
            prmtop_file: Path to the AMBER topology file (.prmtop)
            rst7_file: Optional path to the coordinate file (.rst7)

        Returns:
            bool: True if validation succeeded and file was saved
        """
        console = self.processor.console

        if not os.path.exists(prmtop_file):
            console.print(f"[red]Cannot validate topology - file not found: {prmtop_file}[/red]")
            return False

        # Show an info panel explaining what this step does. The output below
        # can look intimidating to a new user (especially the polarizability
        # inference cascade when many forked atom types are present), so a
        # short orientation up front helps frame what's expected vs. alarming.
        self._print_parmed_info_panel(prmtop_file, console)

        console.print(f"[blue]Running ParmEd validation on {os.path.basename(prmtop_file)}...[/blue]")

        try:
            # Load topology file
            if rst7_file and os.path.exists(rst7_file):
                parm = pmd.load_file(prmtop_file, xyz=rst7_file)
            else:
                parm = pmd.load_file(prmtop_file)

            # Capture warnings during validity check
            # Filter to ParmEd warnings only — Python DeprecationWarnings from
            # ParmEd internals (e.g., "'count' is passed as positional argument")
            # would otherwise flood the output with 100k+ spurious warnings.
            from parmed.exceptions import ParmedWarning
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("ignore")
                warnings.simplefilter("always", category=ParmedWarning)

                # Run validity check
                console.print("[blue]  Running validity check...[/blue]")
                check_validity(parm, warnings)

                # Display any warnings found
                if w:
                    console.print(f"[yellow]  ⚠ Found {len(w)} validation warning(s):[/yellow]")
                    for warning in w[:5]:  # Show first 5 warnings
                        console.print(f"[yellow]    • {warning.message}[/yellow]")
                    if len(w) > 5:
                        console.print(f"[yellow]    ... and {len(w) - 5} more[/yellow]")
                else:
                    console.print("[green]  ✓ Initial validity check passed[/green]")

            # Run rediscover_molecules to fix ATOMS_PER_MOLECULE section
            console.print("[blue]  Running rediscover_molecules() to fix molecule definitions...[/blue]")

            atom_reorder = parm.rediscover_molecules(solute_ions=True, fix_broken=True)

            if atom_reorder is not None:
                console.print("[yellow]  ⚠ Molecule atoms were not contiguous - atoms have been reordered[/yellow]")
                # rediscover_molecules reorders parm.atoms but does NOT rebuild
                # parm.parm_data arrays (CHARGE, ATOM_TYPE_INDEX, etc.).
                # Rebuild them so downstream code (e.g., params1264) sees correct values.
                parm._xfer_atom_info()
                if parm.coordinates is None:
                    console.print("[yellow]  ⚠ No coordinates loaded - coordinate file may need regeneration[/yellow]")
            else:
                console.print("[green]  ✓ Molecule definitions are correct[/green]")

            # Apply 12-6-4 C4 terms if user selected 12-6-4 ion parameters
            lj_params = self._get_12_6_4_params()
            if lj_params and lj_params.get('needed'):
                self._apply_12_6_4(parm, lj_params, console=console,
                                   prmtop_file=prmtop_file, rst7_file=rst7_file)
            else:
                if lj_params and lj_params.get('incompatible'):
                    console.print(f"[yellow]  ⚠ Skipping 12-6-4 LJ potential: {lj_params.get('ion_selection_name')} "
                                f"is not a validated pairing with {lj_params.get('protein_ff')}[/yellow]")
                    console.print("[grey50]    The Li/Merz 12-6-4 ion parameters were calibrated on ff14SB-family point "
                                "charges; ff15ipq uses the implicitly-polarized IPQ charge model, so the combination "
                                "is unvalidated. ParmEd can apply it technically — do so manually if you have validated it.[/grey50]")
                # Write a parmed replay script anyway (no add12_6_4 line; just
                # checkValidity + rediscover + outparm) so the user can rerun
                # the validation pass manually.
                self._write_parmed_replay_script(
                    prmtop_file, rst7_file, None, None, None, console
                )

            # Save the corrected topology
            console.print(f"[blue]  Saving corrected topology to {os.path.basename(prmtop_file)}...[/blue]")
            parm.save(prmtop_file, overwrite=True)

            # Also save coordinate file if it was loaded
            if rst7_file and os.path.exists(rst7_file) and parm.coordinates is not None:
                parm.save(rst7_file, overwrite=True)
                console.print(f"[green]  ✓ Saved corrected topology and coordinates[/green]")
            else:
                console.print(f"[green]  ✓ Saved corrected topology[/green]")

            console.print("[bold green]ParmEd validation completed successfully[/bold green]")
            return True

        except Exception as e:
            from rich.markup import escape
            console.print(f"[red]✗ ParmEd validation failed: {escape(str(e))}[/red]")
            import traceback
            logger.error(f"ParmEd validation error: {traceback.format_exc()}")
            return False

    def _check_and_configure_molecules(self, tleap_info: dict, batch_mode: bool = False) -> bool:
        """
        Check if PDB structure needs molecule grouping configuration.

        NOTE: This method is now primarily used for the TER validation workflow (option 4).
        Structure preparation for option 2 (generate_single_state_tleap) happens in
        _prepare_structure_for_tleap() instead.

        This runs BEFORE TER record validation to ensure atoms are ordered
        correctly by molecule (bonded components must be contiguous).

        Args:
            tleap_info: Dictionary containing tLEaP file info including PDB files
            batch_mode: If True, skip interactive prompts and use existing configuration

        Returns:
            bool: True if configuration succeeded or wasn't needed, False if cancelled
        """
        from proprep.tleap_prep.pdb_molecule_configurator import PDBMoleculeConfigurator
        from pathlib import Path

        console = self.processor.console

        # Get PDB files from tleap_info
        pdb_files = tleap_info.get('pdb_files', [])
        if not pdb_files:
            # No PDB files to process
            return True

        # Process each PDB file
        for pdb_file in pdb_files:
            if not os.path.exists(pdb_file):
                if not batch_mode:
                    console.print(f"[yellow]Warning: PDB file not found: {pdb_file}[/yellow]")
                continue

            # Check if configuration already exists
            pdb_path = Path(pdb_file)
            config_file = pdb_path.parent / f"{pdb_path.stem}_molecule_config.json"
            reordered_pdb = pdb_path.parent / f"{pdb_path.stem}_reordered.pdb"

            # If PDB is already reordered (ends with _reordered.pdb), skip
            if pdb_file.endswith("_reordered.pdb"):
                if not batch_mode:
                    console.print(f"[green]✓ PDB already reordered: {pdb_path.name}[/green]")
                    console.print("[grey50]  (Structure preparation was done in option 2)[/grey50]")
                continue

            # If user explicitly skipped reordering in option 2, don't ask again
            if self.get_from_workspace("reordering_skipped"):
                if not batch_mode:
                    console.print(f"[grey50]Reordering was skipped in structure preparation (using {pdb_path.name} as-is)[/grey50]")
                continue

            # If reordered PDB already exists, offer use/redo/skip
            if reordered_pdb.exists():
                if batch_mode:
                    # In batch mode, automatically use the reordered PDB
                    self._update_tleap_file_pdb_reference(tleap_info['tleap_file'], pdb_file, str(reordered_pdb))
                    self.update_workspace("reordered_pdb_file", str(reordered_pdb))
                    continue
                console.print(f"[blue]Found existing reordered PDB: {reordered_pdb.name}[/blue]")
                action = prompt_with_context(
                    self.processor,
                    "Use existing reordered PDB, redo reordering, or skip?",
                    choices=["use", "redo", "skip"],
                    default="use",
                    module="Topology Generator",
                    description="Reordered PDB action",
                )
                if action == "use":
                    self._update_tleap_file_pdb_reference(tleap_info['tleap_file'], pdb_file, str(reordered_pdb))
                    self.update_workspace("reordered_pdb_file", str(reordered_pdb))
                    console.print(f"[green]✓ Using {reordered_pdb.name} in tLEaP file[/green]")
                    continue
                elif action == "skip":
                    console.print(f"[grey50]Skipping reordering (using {pdb_path.name} as-is)[/grey50]")
                    continue

            # Initialize configurator
            configurator = PDBMoleculeConfigurator(console=console)

            # Analyze PDB structure
            if not batch_mode:
                console.print(f"\n[blue]Analyzing PDB structure: {pdb_path.name}[/blue]\n")
            segments = configurator.analyze_pdb_segments(pdb_file)

            if not segments:
                if not batch_mode:
                    console.print("[yellow]No chains found in PDB file[/yellow]")
                continue

            # Display structure summary
            if not batch_mode:
                self._display_structure_summary(segments, console)

            # Check if structure is complex (multiple chains)
            if len(segments) <= 1:
                if not batch_mode:
                    console.print("\n[green]✓ Simple structure - no molecule reordering needed[/green]")
                    console.print("  Proceeding with TER record validation...\n")
                continue

            # Complex structure detected
            if not batch_mode:
                console.print("\n[yellow]⚠ Multi-chain structure detected[/yellow]\n")

            # Check for existing configuration
            config = None
            if config_file.exists():
                if batch_mode:
                    # In batch mode, automatically load existing configuration
                    config = configurator.load_configuration(str(config_file))
                else:
                    console.print(f"[blue]Found existing configuration: {config_file.name}[/blue]")
                    if confirm_with_context(
                self.processor,
                "Load existing configuration?",
                default=True,
                module="Topology Generator",
                description="Load existing tLEaP configuration file",
            ):
                        config = configurator.load_configuration(str(config_file))

            # In batch mode, skip interactive configuration - ParmEd will handle reordering
            if batch_mode and config is None:
                # No config exists, skip for batch mode (ParmEd will fix contiguity issues later)
                continue

            # Ask if user wants to configure molecule grouping
            if config is None:
                if not confirm_with_context(
                    self.processor,
                    "\n[bold]Configure molecule grouping for AMBER?[/bold]\n"
                    "(Reorders PDB atoms to group chains that bond together)\n\n"
                    "[grey50]This is needed when chains form covalent bonds across chain boundaries[/grey50]\n"
                    "[grey50](e.g., disulfides, metal coordination, cofactor bridges)[/grey50]\n\n"
                    "Configure molecule grouping?",
                    default=False,
                    module="Topology Generator",
                    description="Configure molecule grouping for AMBER",
                ):
                    console.print("\n[yellow]Skipping molecule configuration[/yellow]")
                    console.print("[grey50]Note: If AMBER fails with 'intermolecular PRFs' error, you'll need molecule grouping or ParmEd post-processing[/grey50]\n")
                    continue

                # Show info panel
                console.print(Panel(
                    "[bold blue]Molecule Grouping for AMBER[/bold blue]\n\n"
                    "[bold]What it does:[/bold]\n"
                    "Reorders PDB atoms so bonded chains have contiguous atoms.\n"
                    "AMBER requires this for any covalent bonds crossing chain boundaries.\n\n"
                    "[bold]How to think about it:[/bold]\n\n"
                    "[yellow]Example 1:[/yellow] Chains A and B do NOT bond to each other\n"
                    "  → Molecule 1: (protein A, cofactors A)\n"
                    "  → Molecule 2: (protein B, cofactors B)\n"
                    "  → Output: protein_A, cofactors_A, TER, protein_B, cofactors_B, TER\n\n"
                    "[yellow]Example 2:[/yellow] Chains A and B DO bond together (e.g., via cofactor bridge)\n"
                    "  → Molecule 1: (protein A, protein B, cofactors A, cofactors B)\n"
                    "  → Output: protein_A, protein_B, cofactors_A, cofactors_B, TER\n\n"
                    "[bold]Common scenarios:[/bold]\n"
                    "  • Disulfide bonds between protein chains\n"
                    "  • Metal centers coordinated by multiple chains\n"
                    "  • Cofactors that bridge protein chains\n"
                    "  • Any covalent crosslinks across chains\n\n"
                    "[grey50]You'll specify which chains bond together based on your structure's chemistry.[/grey50]",
                    border_style="blue",
                    expand=False
                ))

                # Run interactive configuration
                config = configurator.configure_molecule_grouping(pdb_file)

                if config is None:
                    console.print("[yellow]Configuration cancelled[/yellow]")
                    return False

            # Reorder PDB based on configuration (with renumbering)
            if not batch_mode:
                console.print()
            reordered_file, mapper = configurator.reorder_pdb(
                pdb_file,
                config,
                output_file=str(reordered_pdb)
            )

            if not reordered_file or not mapper:
                if not batch_mode:
                    console.print("[red]✗ Failed to reorder PDB[/red]")
                return False

            # Synchronize RedoxSites with reordered/renumbered structure
            workspace = self.processor.workspace
            redox_sites = workspace.get("detected_redox_sites", [])
            if redox_sites:
                sync_summary = configurator.synchronize_redox_sites(redox_sites, mapper)
                # Save updated redox sites back to workspace
                workspace.set("detected_redox_sites", redox_sites)
                if not batch_mode:
                    console.print(
                        f"[green]✓ RedoxSites synchronized with reordered structure[/green]"
                    )
            elif not batch_mode:
                console.print("[grey50]No RedoxSites to synchronize[/grey50]")

            # Store reordered PDB under its own dedicated workspace key
            workspace.set("reordered_pdb_file", reordered_file)

            # Update tLEaP file to use reordered PDB
            self._update_tleap_file_pdb_reference(tleap_info['tleap_file'], pdb_file, reordered_file)
            if not batch_mode:
                console.print(f"[green]✓ Updated tLEaP file to use {Path(reordered_file).name}[/green]\n")

        return True

    def _display_structure_summary(self, segments, console):
        """Display a summary table of PDB structure."""
        from rich.table import Table

        table = Table(title="Structure Summary")
        table.add_column("Chain", style="blue")
        table.add_column("Type", style="yellow")
        table.add_column("Protein Range", style="green")
        table.add_column("Protein Count", justify="right", style="green")
        table.add_column("Cofactor Residues", style="magenta")

        for seg in segments:
            seg_type = "Protein" if seg.is_protein else "Hetero"
            protein_range = seg.get_protein_range_str()
            protein_count = str(len(seg.protein_residues)) if seg.protein_residues else "-"
            cofactor_str = seg.get_cofactor_str()

            table.add_row(
                seg.chain_id,
                seg_type,
                protein_range,
                protein_count,
                cofactor_str
            )

        console.print(table)

    def _update_tleap_file_pdb_reference(self, tleap_file: str, old_pdb: str, new_pdb: str):
        """Update tLEaP file to reference reordered PDB instead of original."""
        try:
            with open(tleap_file, 'r') as f:
                content = f.read()

            # Replace PDB filename in loadpdb commands
            old_filename = os.path.basename(old_pdb)
            new_filename = os.path.basename(new_pdb)

            content = content.replace(old_filename, new_filename)

            with open(tleap_file, 'w') as f:
                f.write(content)

        except Exception as e:
            logger.error(f"Error updating tLEaP file: {e}")

    def _validate_and_fix_single_pdb_ter_records(self, pdb_file: str) -> bool:
        """
        Validate and fix TER records in a single PDB file.

        Args:
            pdb_file: Path to PDB file

        Returns:
            bool: True if validation passed or fixes were applied successfully
        """
        console = self.processor.console

        if not os.path.exists(pdb_file):
            console.print(f"[red]PDB file not found: {pdb_file}[/red]")
            return False

        # Check if TER record validation is needed
        needs_fixing, issues = self._check_ter_records(pdb_file)

        if needs_fixing:
            console.print(f"[yellow]TER record issues found in {os.path.basename(pdb_file)}:[/yellow]")
            for issue in issues:
                console.print(f"  • {issue}")
            console.print(f"\n[grey50]TER records separate molecules in PDB files. TLeaP requires them to "
                          f"correctly identify molecular boundaries (protein chain termini, "
                          f"ligands, cofactors, ions). Missing TER records cause TLeaP to treat "
                          f"separate molecules as one continuous chain, leading to incorrect bonds.[/grey50]")

            if self._fix_ter_records(pdb_file, issues):
                console.print(f"[green]✓ Inserted {len(issues)} TER record(s) in {os.path.basename(pdb_file)}[/green]")
                return True
            else:
                console.print(f"[red]✗ Failed to fix TER records in {os.path.basename(pdb_file)}[/red]")
                return False
        else:
            console.print(f"[green]✓ TER records valid in {os.path.basename(pdb_file)}[/green]")
            return True

    def _validate_and_fix_ter_records(self, tleap_info: dict, quiet: bool = False) -> bool:
        """
        Validate and fix TER records in PDB files to prevent tLEaP bonding issues.

        Args:
            tleap_info: Dictionary containing tLEaP file info including PDB files
            quiet: If True, auto-fix issues without prompts or verbose output

        Returns:
            bool: True if validation passed or fixes were applied successfully
        """
        console = self.processor.console
        all_files_valid = True

        for pdb_file in tleap_info['pdb_files']:
            if not os.path.exists(pdb_file):
                if not quiet:
                    console.print(f"[red]PDB file not found: {pdb_file}[/red]")
                all_files_valid = False
                continue

            # Check if TER record validation is needed
            needs_fixing, issues = self._check_ter_records(pdb_file)

            if needs_fixing:
                if not quiet:
                    console.print(f"[yellow]TER record issues found in {os.path.basename(pdb_file)}:[/yellow]")
                    for issue in issues:
                        console.print(f"  • {issue}")

                if self._fix_ter_records(pdb_file, issues):
                    if not quiet:
                        console.print(f"[green]✓ Inserted {len(issues)} TER record(s) in {os.path.basename(pdb_file)}[/green]")
                else:
                    if not quiet:
                        console.print(f"[red]✗ Failed to fix TER records in {os.path.basename(pdb_file)}[/red]")
                    all_files_valid = False
            elif not quiet:
                console.print(f"[green]✓ TER records valid in {pdb_file}[/green]")

        return all_files_valid

    def _check_ter_records(self, pdb_file: str) -> tuple:
        """
        Check if TER records are properly placed in PDB file.
        
        Returns:
            tuple: (needs_fixing: bool, issues: list)
        """
        issues = []
        
        try:
            with open(pdb_file, 'r') as f:
                lines = f.readlines()
            
            # Track chains, residue info, and TER positions
            prev_chain = None
            prev_record_type = None
            prev_residue_id = None
            prev_residue_name = None
            prev_atom_name = None
            line_num = 0
            
            for line in lines:
                line_num += 1
                
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    chain = line[21:22].strip() if len(line) > 21 else ''
                    residue_name = line[17:20].strip() if len(line) > 19 else ''
                    residue_id = line[22:26].strip() if len(line) > 25 else ''
                    atom_name = line[12:16].strip() if len(line) > 15 else ''
                    record_type = line[:6].strip()

                    # Check for chain changes without TER (after ATOM records)
                    if prev_chain is not None and chain != prev_chain and prev_record_type == 'ATOM':
                        issues.append(f"Missing TER record between chain {prev_chain} and {chain} at line {line_num}")

                    # Check for ATOM -> HETATM transition without TER
                    if prev_record_type == 'ATOM' and record_type == 'HETATM' and chain == prev_chain:
                        issues.append(f"Missing TER record between ATOM and HETATM in chain {chain} at line {line_num}")

                    # Check for new residue after OXT (C-terminus → next molecule)
                    if prev_atom_name == 'OXT' and residue_id != prev_residue_id:
                        issues.append(f"Missing TER record after OXT (C-terminus of {prev_residue_name} {prev_residue_id}) at line {line_num}")

                    # Check for HETATM residue changes without TER (different residue name OR different residue ID)
                    if (prev_record_type == 'HETATM' and record_type == 'HETATM' and
                        chain == prev_chain and prev_residue_id and
                        (residue_name != prev_residue_name or residue_id != prev_residue_id)):
                        issues.append(f"Missing TER record between HETATM {prev_residue_name} {prev_residue_id} and {residue_name} {residue_id} in chain {chain} at line {line_num}")

                    prev_chain = chain
                    prev_record_type = record_type
                    prev_residue_id = residue_id
                    prev_residue_name = residue_name
                    prev_atom_name = atom_name
                    
                elif line.startswith('TER'):
                    # Reset tracking after TER. prev_atom_name must reset too —
                    # otherwise a TER following an OXT atom leaves prev_atom_name='OXT'
                    # and the next ATOM line trips the "Missing TER after OXT" check
                    # (printing 'None None' for the already-reset residue identity).
                    prev_record_type = 'TER'
                    prev_residue_id = None
                    prev_residue_name = None
                    prev_atom_name = None
                    prev_chain = None
            
        except Exception as e:
            issues.append(f"Error reading PDB file: {e}")
        
        return len(issues) > 0, issues

    def _fix_ter_records(self, pdb_file: str, issues: list) -> bool:
        """
        Fix TER record issues in PDB file.
        
        Args:
            pdb_file: Path to PDB file
            issues: List of issues to fix
            
        Returns:
            bool: True if fixes were applied successfully
        """
        try:
            with open(pdb_file, 'r') as f:
                lines = f.readlines()
            
            # Create backup
            backup_file = f"{pdb_file}.ter_backup"
            with open(backup_file, 'w') as f:
                f.writelines(lines)
            
            # Apply fixes
            fixed_lines = []
            prev_chain = None
            prev_record_type = None
            prev_residue_id = None
            prev_residue_name = None
            prev_atom_name = None

            for line in lines:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    chain = line[21:22].strip() if len(line) > 21 else ''
                    residue_name = line[17:20].strip() if len(line) > 19 else ''
                    residue_id = line[22:26].strip() if len(line) > 25 else ''
                    atom_name = line[12:16].strip() if len(line) > 15 else ''
                    record_type = line[:6].strip()

                    # Insert TER before chain changes (after ATOM records)
                    if (prev_chain is not None and chain != prev_chain and
                        prev_record_type == 'ATOM'):
                        fixed_lines.append("TER\n")

                    # Insert TER before ATOM -> HETATM transition
                    elif (prev_record_type == 'ATOM' and record_type == 'HETATM' and
                          chain == prev_chain):
                        fixed_lines.append("TER\n")

                    # Insert TER after OXT (C-terminus → next molecule)
                    elif prev_atom_name == 'OXT' and residue_id != prev_residue_id:
                        fixed_lines.append("TER\n")

                    # Insert TER between different HETATM residues (different name OR different ID)
                    elif (prev_record_type == 'HETATM' and record_type == 'HETATM' and
                          chain == prev_chain and prev_residue_id and
                          (residue_name != prev_residue_name or residue_id != prev_residue_id)):
                        fixed_lines.append("TER\n")

                    prev_chain = chain
                    prev_record_type = record_type
                    prev_residue_id = residue_id
                    prev_residue_name = residue_name
                    prev_atom_name = atom_name
                    
                elif line.startswith('TER'):
                    prev_record_type = 'TER'
                    prev_residue_id = None
                    prev_residue_name = None
                
                fixed_lines.append(line)
            
            # Write fixed file
            with open(pdb_file, 'w') as f:
                f.writelines(fixed_lines)
            
            return True
            
        except Exception as e:
            self.processor.console.print(f"[red]Error fixing TER records: {e}[/red]")
            return False

    def _display_leap_log_messages(self, tleap_file: str):
        """
        Parse and display warnings, errors, and notes from leap.log file.
        
        Args:
            tleap_file: Path to the tLEaP input file (used to find corresponding leap.log)
        """
        console = self.processor.console
        
        # Look for leap.log in current directory
        leap_log_file = "leap.log"
        
        if not os.path.exists(leap_log_file):
            console.print(f"[yellow]No leap.log file found for {tleap_file}[/yellow]")
            return
        
        try:
            # Parse leap.log file
            warnings, errors, notes, summary = self._parse_leap_log(leap_log_file)
            
            # Display summary first if available
            if summary:
                console.print()  # Blank line before tLEaP Summary
                console.print(f"[bold blue]tLEaP Summary: {summary}[/bold blue]")
            
            # Display errors (most important)
            if errors:
                console.print(f"\n[bold red]Errors ({len(errors)}):[/bold red]")
                for i, error in enumerate(errors, 1):
                    console.print(f"[red]Error {i}:[/red]")
                    for line in error:
                        console.print(f"  {line}")
                    console.print()  # Blank line between errors
            
            # Display warnings
            if warnings:
                console.print(f"\n[bold yellow]Warnings ({len(warnings)}):[/bold yellow]")
                for i, warning in enumerate(warnings, 1):
                    console.print(f"[yellow]Warning {i}:[/yellow]")
                    for line in warning:
                        console.print(f"  {line}")
                    console.print()  # Blank line between warnings
            
            # Display notes
            if notes:
                console.print(f"\n[bold blue]Notes ({len(notes)}):[/bold blue]")
                for i, note in enumerate(notes, 1):
                    console.print(f"[blue]Note {i}:[/blue]")
                    for line in note:
                        console.print(f"  {line}")
                    console.print()  # Blank line between notes
            
            # If no messages found but log exists
            if not warnings and not errors and not notes and not summary:
                console.print(f"[green]No warnings, errors, or notes found in leap.log[/green]")
                
        except Exception as e:
            console.print(f"[red]Error reading leap.log: {e}[/red]")

    def _parse_leap_log(self, leap_log_file: str) -> tuple:
        """
        Parse leap.log file to extract warning, error, and note blocks.
        
        Returns:
            tuple: (warnings, errors, notes, summary)
                - warnings: list of warning text blocks (each block is list of lines)
                - errors: list of error text blocks (each block is list of lines)  
                - notes: list of note text blocks (each block is list of lines)
                - summary: final summary line or None
        """
        warnings = []
        errors = []
        notes = []
        summary = None
        
        try:
            with open(leap_log_file, 'r') as f:
                lines = f.readlines()
            
            current_block = []
            current_type = None
            
            for line in lines:
                line = line.rstrip()  # Remove trailing whitespace
                
                # Check for start of new message block
                if ": Warning!" in line:
                    # Save previous block if exists
                    if current_block and current_type:
                        if current_type == "warning":
                            warnings.append(current_block)
                        elif current_type == "error":
                            errors.append(current_block)
                        elif current_type == "note":
                            notes.append(current_block)
                    
                    # Start new warning block
                    current_block = [line]
                    current_type = "warning"
                    
                elif ": Error!" in line:
                    # Save previous block if exists
                    if current_block and current_type:
                        if current_type == "warning":
                            warnings.append(current_block)
                        elif current_type == "error":
                            errors.append(current_block)
                        elif current_type == "note":
                            notes.append(current_block)
                    
                    # Start new error block
                    current_block = [line]
                    current_type = "error"
                    
                elif ": Note." in line:
                    # Save previous block if exists
                    if current_block and current_type:
                        if current_type == "warning":
                            warnings.append(current_block)
                        elif current_type == "error":
                            errors.append(current_block)
                        elif current_type == "note":
                            notes.append(current_block)
                    
                    # Start new note block
                    current_block = [line]
                    current_type = "note"
                    
                elif line.startswith("Exiting LEaP:") and ("Errors =" in line and "Warnings =" in line):
                    # Final summary line
                    summary = line
                    
                elif current_block and current_type:
                    # Check if this line ends the current block
                    stripped_line = line.strip()
                    if not stripped_line or stripped_line == ">" or stripped_line == ">>":
                        # Empty line or tLEaP prompt - end of current block
                        if current_type == "warning":
                            warnings.append(current_block)
                        elif current_type == "error":
                            errors.append(current_block)
                        elif current_type == "note":
                            notes.append(current_block)
                        
                        current_block = []
                        current_type = None
                    else:
                        # Continue current block with non-empty, non-prompt line
                        # Only add one line after the path-containing line to avoid irregular output
                        if len(current_block) == 1:  # Only the initial line exists
                            current_block.append(line)
                        # Skip additional lines to keep output clean
            
            # Handle last block if file doesn't end with blank line
            if current_block and current_type:
                if current_type == "warning":
                    warnings.append(current_block)
                elif current_type == "error":
                    errors.append(current_block)
                elif current_type == "note":
                    notes.append(current_block)
            
        except Exception as e:
            # Return empty results on error
            pass
        
        return warnings, errors, notes, summary

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
    # Two-pass tLEaP execution for accurate ion calculation
    # We always run tLEaP twice:
    # Pass 1: Get actual water count and charge (with early termination)
    # Pass 2: Full run with correct ion counts based on actual values

    def _calculate_ions_from_water_count(self, n_waters: int, net_charge: int,
                                          target_molarity: float = 0.15,
                                          cation_charge: int = +1,
                                          anion_charge: int = -1) -> dict:
        """Single-salt entry point (kept for back-compat). See
        _calculate_multi_salt_ions for the general case."""
        cation = {'label': 'cat', 'unit': 'cat', 'charge': int(cation_charge)}
        anion = {'label': 'an', 'unit': 'an', 'charge': int(anion_charge)}
        salts = [{'cation': cation, 'anion': anion, 'concentration': float(target_molarity)}]
        return self._calculate_multi_salt_ions(n_waters, net_charge, salts, neutralize_index=0)

    def _calculate_multi_salt_ions(self, n_waters: int, net_charge: int,
                                    salts: list, neutralize_index: int = 0) -> dict:
        """Compute per-salt bulk ion counts plus a single neutralizer pair.

        For a single salt this preserves the SPLIT method (Machado & Pantano,
        JCTC 2020) so existing single-salt outputs are unchanged. For
        multi-salt mixes each salt contributes pure bulk at its formula
        ratio (a*n_pair cations + b*n_pair anions); the protein charge is
        balanced once using the chosen salt's ions. The error vs. true SPLIT
        is O(|Q|/N_w) — a few mM at most, which is well below the
        equilibration-density inflation users already accept.

        Args:
            n_waters: Water count from the tLEaP info pass
            net_charge: System net charge (signed)
            salts: list of {'cation', 'anion', 'concentration'} dicts
            neutralize_index: which salt's ions are used for neutralization

        Returns:
            ion_counts dict with:
                n_waters, net_charge, neutralize_index,
                n_neutralize_cation, n_neutralize_anion,
                salts: [{cation, anion, concentration, n_pairs,
                         n_bulk_cations, n_bulk_anions}, ...]
        """
        water_molarity = 55.5
        Q = int(round(net_charge))
        salts = salts or []
        if not salts:
            # No salt configured — emit empty structure for caller to handle
            return {
                'n_waters': n_waters,
                'net_charge': Q,
                'neutralize_index': 0,
                'n_neutralize_cation': max(0, -Q),
                'n_neutralize_anion': max(0, Q),
                'salts': [],
            }

        idx = neutralize_index if 0 <= neutralize_index < len(salts) else 0
        neut_salt = salts[idx]
        cp = abs(int(neut_salt['cation']['charge']))
        ca = abs(int(neut_salt['anion']['charge']))

        # ---- Single-salt path: preserve SPLIT exactly as before ----
        if len(salts) == 1:
            s = salts[0]
            target_molarity = float(s['concentration'])
            symmetric = (cp == 1 and ca == 1)
            a, b = _salt_formula_units(cp, ca)

            if target_molarity > 0:
                n_salt_exact = n_waters * target_molarity / water_molarity
                n_salt = int(round(n_salt_exact))
            else:
                n_salt_exact = 0.0
                n_salt = 0

            if symmetric:
                if target_molarity > 0:
                    n_cations_exact = n_salt_exact - Q / 2.0
                    n_cations = max(0, int(round(n_cations_exact)))
                    n_anions = max(0, n_cations + Q)
                    if n_cations + Q < 0:
                        n_anions = 0
                        n_cations = abs(Q)
                else:
                    n_cations = max(0, -Q)
                    n_anions = max(0, Q)
                n_neut_cat = abs(Q) if Q < 0 else 0
                n_neut_an = abs(Q) if Q > 0 else 0
                n_bulk_cat = max(0, n_cations - n_neut_cat)
                n_bulk_an = max(0, n_anions - n_neut_an)
            else:
                n_bulk_cat = a * n_salt
                n_bulk_an = b * n_salt
                if Q < 0:
                    n_neut_cat = -(-(-Q) // cp)
                    excess = cp * n_neut_cat + Q
                    n_neut_an = excess // ca
                elif Q > 0:
                    n_neut_an = -(-Q // ca)
                    excess = ca * n_neut_an - Q
                    n_neut_cat = excess // cp
                else:
                    n_neut_cat = 0
                    n_neut_an = 0

            return {
                'n_waters': n_waters,
                'net_charge': Q,
                'neutralize_index': 0,
                'n_neutralize_cation': n_neut_cat,
                'n_neutralize_anion': n_neut_an,
                'salts': [{
                    'cation': dict(s['cation']),
                    'anion': dict(s['anion']),
                    'concentration': target_molarity,
                    'n_pairs': n_salt,
                    'n_bulk_cations': n_bulk_cat,
                    'n_bulk_anions': n_bulk_an,
                }],
            }

        # ---- Multi-salt path: bulk = a:b ratio per salt, neutralize via chosen salt ----
        salt_results = []
        for s in salts:
            conc = float(s['concentration'])
            cp_i = abs(int(s['cation']['charge']))
            ca_i = abs(int(s['anion']['charge']))
            a_i, b_i = _salt_formula_units(cp_i, ca_i)
            n_pair = int(round(n_waters * conc / water_molarity)) if conc > 0 else 0
            salt_results.append({
                'cation': dict(s['cation']),
                'anion': dict(s['anion']),
                'concentration': conc,
                'n_pairs': n_pair,
                'n_bulk_cations': a_i * n_pair,
                'n_bulk_anions': b_i * n_pair,
            })

        # Neutralize via salts[idx]'s ions; same asymmetric handling we use
        # for single asymmetric salts.
        if Q < 0:
            n_neut_cat = -(-(-Q) // cp)
            excess = cp * n_neut_cat + Q
            n_neut_an = excess // ca
        elif Q > 0:
            n_neut_an = -(-Q // ca)
            excess = ca * n_neut_an - Q
            n_neut_cat = excess // cp
        else:
            n_neut_cat = 0
            n_neut_an = 0

        return {
            'n_waters': n_waters,
            'net_charge': Q,
            'neutralize_index': idx,
            'n_neutralize_cation': n_neut_cat,
            'n_neutralize_anion': n_neut_an,
            'salts': salt_results,
        }

    def _update_template_ion_counts(self, original_template: str, ion_counts: dict,
                                     cation: dict = None, anion: dict = None) -> str:
        """
        Update a tLEaP template with correct ion counts.

        Replaces the existing addions/addionsrand commands with updated counts.
        Reuses the multi-salt builder so the rewritten template matches what
        we'd emit from scratch.

        Args:
            original_template: Path to original tLEaP input file
            ion_counts: Dictionary from _calculate_multi_salt_ions()
            cation, anion: Legacy single-salt overrides (only used if ion_counts
                lacks a salts list, e.g. from very old test fixtures).

        Returns:
            Path to updated template, or None on failure
        """
        import re

        try:
            with open(original_template, 'r') as f:
                content = f.read()

            n_waters = ion_counts.get('n_waters')
            header = f"# Ion calculation updated with actual water count ({n_waters} waters)"
            ion_block = self._build_ion_commands_for_template(
                ion_counts, cation=cation, anion=anion
            ).rstrip()
            new_ion_commands = f"{header}\n{ion_block}"

            # Replace existing ion commands. The block may span:
            #   - leading "# ..." comments,
            #   - one addions line,
            #   - one or more addionsrand lines (multi-salt emits one per salt),
            #     each optionally preceded by a "# ..." comment.
            pattern = (r'(#[^\n]*ion[^\n]*\n)*'
                       r'addions[^\n]+'
                       r'(\n(#[^\n]*\n)?addionsrand[^\n]+)*')
            content, n_subs = re.subn(pattern, new_ion_commands, content,
                                       count=1, flags=re.IGNORECASE)

            if n_subs == 0:
                # Fallback: two consecutive addions lines (legacy minimal seed)
                pattern = r'addions[^\n]+\naddions[^\n]+'
                content, n_subs = re.subn(pattern, new_ion_commands, content,
                                           count=1, flags=re.IGNORECASE)

            if n_subs == 0:
                logger.warning("Could not find ion commands to replace in template")
                return original_template  # Return original if can't update

            # Write updated template
            base = os.path.splitext(original_template)[0]
            updated_file = f"{base}_accurate_ions.in"

            with open(updated_file, 'w') as f:
                f.write(content)

            logger.info(f"Created template with accurate ion counts: {updated_file}")
            return updated_file

        except Exception as e:
            logger.error(f"Error updating template ion counts: {e}")
            return None

    def _calculate_protein_dimensions(self, pdb_path: str) -> dict:
        """
        Calculate protein bounding box dimensions from PDB file.

        Returns:
            Dictionary with dimensions, volume estimates, and ion calculations
        """
        try:
            # Load structure with parmed
            structure = pmd.load_file(pdb_path)

            # Get all atom coordinates
            coords = structure.coordinates
            if coords is None or len(coords) == 0:
                return None

            # Calculate bounding box
            x_coords = coords[:, 0]
            y_coords = coords[:, 1]
            z_coords = coords[:, 2]

            x_min, x_max = x_coords.min(), x_coords.max()
            y_min, y_max = y_coords.min(), y_coords.max()
            z_min, z_max = z_coords.min(), z_coords.max()

            # Protein dimensions
            lx = x_max - x_min
            ly = y_max - y_min
            lz = z_max - z_min

            # Get net charge - priority order:
            # 1. From workspace (Protonation State Analyzer result)
            # 2. Estimate from residue formal charges at pH 7
            net_charge = None
            charge_source = None
            charge_warnings = []

            # Try workspace first (from Protonation State Analyzer)
            workspace_charge = self.get_from_workspace("net_charge")
            if workspace_charge is not None:
                net_charge = round(workspace_charge)
                charge_source = "protonation_analyzer"
            else:
                # Estimate from residue formal charges at physiological pH (~7)
                # Standard amino acid charges: Asp/Glu = -1, Lys/Arg = +1, His = 0 (mostly)
                # N-terminus = +1, C-terminus = -1
                net_charge, charge_warnings = self._estimate_formal_charge(structure)
                charge_source = "residue_estimate"

            return {
                'lx': lx,
                'ly': ly,
                'lz': lz,
                'x_range': (x_min, x_max),
                'y_range': (y_min, y_max),
                'z_range': (z_min, z_max),
                'net_charge': net_charge,
                'charge_source': charge_source,
                'charge_warnings': charge_warnings
            }

        except Exception as e:
            logger.warning(f"Could not calculate protein dimensions: {e}")
            return None

    def _estimate_formal_charge(self, structure) -> tuple:
        """
        Estimate net formal charge from residue names at physiological pH (~7).

        This is a rough estimate based on standard amino acid pKa values.
        For accurate charges, run the Protonation State Analyzer first.

        Returns:
            Tuple of (net_charge: int, warnings: list of strings)
        """
        # Standard amino acids (3-letter codes)
        standard_amino_acids = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
            # Common variants
            'HID', 'HIE', 'HIP',  # Histidine protonation states
            'CYX', 'CYM',  # Cysteine variants
            'ASH', 'GLH',  # Protonated Asp/Glu
        }

        # Formal charges at pH 7 for standard amino acids
        charged_residues = {
            'ASP': -1, 'GLU': -1,  # Acidic (deprotonated at pH 7)
            'LYS': +1, 'ARG': +1,  # Basic (protonated at pH 7)
            'HIS': 0, 'HID': 0, 'HIE': 0,  # Neutral histidine
            'HIP': +1,  # Protonated histidine
            'CYM': -1,  # Deprotonated cysteine
            'ASH': 0, 'GLH': 0,  # Protonated Asp/Glu (neutral)
        }

        # Common solvent/ion residues to ignore
        solvent_ions = {'WAT', 'HOH', 'TIP', 'TIP3', 'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'SOL'}

        # Build set of residues to exclude (from detected redox sites)
        redox_residues = set()
        detected_sites = self.get_from_workspace("detected_redox_sites")
        if detected_sites:
            for site in detected_sites:
                # RedoxSite objects have residue_groups: Dict[(chain, resid, ins_code), coords]
                if hasattr(site, 'residue_groups'):
                    for (chain, resid, ins_code) in site.residue_groups.keys():
                        redox_residues.add((chain, resid))
                # Also exclude center residues
                if hasattr(site, 'centers'):
                    for center in site.centers:
                        redox_residues.add((center.chain, center.resid))

        net_charge = 0
        chain_termini = set()
        warnings = []
        unknown_hetero_residues = []
        excluded_redox_residues = []

        for residue in structure.residues:
            resname = residue.name.upper()

            # Skip solvent and common ions
            if resname in solvent_ions:
                continue

            # Get residue identifiers
            chain = residue.chain if hasattr(residue, 'chain') else ''
            resid = residue.number if hasattr(residue, 'number') else 0

            # Check if this residue is part of a detected redox site
            if (chain, resid) in redox_residues:
                excluded_redox_residues.append(f"{resname} {chain}:{resid}")
                continue

            # Count charged standard residues
            if resname in charged_residues:
                net_charge += charged_residues[resname]

            # Check for non-standard residues (hetero groups)
            elif resname not in standard_amino_acids:
                # This is a non-standard residue we can't estimate charge for
                unknown_hetero_residues.append(f"{resname} {chain}:{resid}")

            # Track chains for terminal charges
            if chain:
                chain_termini.add(chain)

        # Build warnings
        if excluded_redox_residues:
            warnings.append(f"Excluded {len(excluded_redox_residues)} redox site residue(s) from charge calculation")

        if unknown_hetero_residues:
            # Limit display to first few
            display_limit = 5
            hetero_list = unknown_hetero_residues[:display_limit]
            if len(unknown_hetero_residues) > display_limit:
                hetero_list.append(f"...and {len(unknown_hetero_residues) - display_limit} more")
            warnings.append(
                f"Cannot estimate charge for {len(unknown_hetero_residues)} non-standard residue(s): "
                f"{', '.join(hetero_list)}"
            )

        # N-terminus (+1) and C-terminus (-1) per chain cancel out

        return net_charge, warnings

    def _calculate_vdw_dimensions_with_forcefield(self, pdb_path: str,
                                                     custom_forcefields: dict = None) -> tuple:
        """
        Calculate vdW bounding box dimensions using forcefield R* values.

        Uses ForceFieldLoader to get accurate R* values from selected forcefields
        and any custom frcmod files (e.g., from redox site parameterization).

        Args:
            pdb_path: Path to PDB file
            custom_forcefields: Optional dict of custom forcefields (for multi-microstate).
                               If None, loads from workspace (for single-state).

        Returns:
            Tuple of (vdw_lx, vdw_ly, vdw_lz) or None if forcefields not available
        """
        try:
            # Lazy import to avoid circular dependency
            from proprep.forcefield_prep import ForceFieldLoader

            # Check if forcefields have been selected
            selected_ffs = self.get_from_workspace("selected_standard_forcefields")
            if not selected_ffs:
                logger.debug("No forcefields selected, cannot use accurate R* values")
                return None

            # Initialize ForceFieldLoader
            loader = ForceFieldLoader()

            # Load standard forcefields
            loader.load_from_selected_forcefields(selected_ffs)

            # Load custom forcefields
            if custom_forcefields:
                # Multi-microstate: custom forcefields passed as parameter
                loader.load_from_custom_forcefield_dict(custom_forcefields)
            else:
                # Single-state: custom forcefields from workspace
                single_state_ffs = self.get_from_workspace("single_state_selected_forcefields")
                if single_state_ffs:
                    loader.load_from_custom_forcefield_dict(single_state_ffs)

            # Check if we have enough parameters
            stats = loader.get_statistics()
            if stats['n_nonbonded_params'] == 0:
                logger.warning("No R* parameters loaded from forcefields")
                return None

            logger.info(f"Loaded {stats['n_nonbonded_params']} R* parameters for vdW calculation")

            # Load structure and calculate vdW bounding box
            structure = pmd.load_file(pdb_path)
            coords = structure.coordinates
            if coords is None or len(coords) == 0:
                return None

            # Calculate vdW bounding box: for each atom, add R* to position
            x_min_vdw = float('inf')
            x_max_vdw = float('-inf')
            y_min_vdw = float('inf')
            y_max_vdw = float('-inf')
            z_min_vdw = float('inf')
            z_max_vdw = float('-inf')

            for i, atom in enumerate(structure.atoms):
                x, y, z = coords[i]

                # Get R* value for this atom
                resname = atom.residue.name if atom.residue else 'UNK'
                atom_name = atom.name
                element = atom.element_name if hasattr(atom, 'element_name') else ''

                r_star = loader.get_r_star_for_atom(resname, atom_name, element)

                # Update vdW bounding box
                x_min_vdw = min(x_min_vdw, x - r_star)
                x_max_vdw = max(x_max_vdw, x + r_star)
                y_min_vdw = min(y_min_vdw, y - r_star)
                y_max_vdw = max(y_max_vdw, y + r_star)
                z_min_vdw = min(z_min_vdw, z - r_star)
                z_max_vdw = max(z_max_vdw, z + r_star)

            vdw_lx = x_max_vdw - x_min_vdw
            vdw_ly = y_max_vdw - y_min_vdw
            vdw_lz = z_max_vdw - z_min_vdw

            # Log any warnings about missing atom types
            ff_warnings = loader.get_warnings()
            for warning in ff_warnings:
                logger.warning(warning)

            return (vdw_lx, vdw_ly, vdw_lz)

        except Exception as e:
            logger.warning(f"Could not calculate vdW dimensions with forcefield: {e}")
            return None

    def _calculate_diagonal_cut_factor(self, pdb_path: str, buffer: float) -> float:
        """
        Calculate diagonal cut scaling factor for truncated octahedron.

        Implements tLEaP's ToolOctBoxCheck() algorithm. For a truncated octahedron,
        the buffer distance must be maintained not just along axes but also along
        the diagonal faces.

        Args:
            pdb_path: Path to PDB file
            buffer: Solvation buffer distance in Angstroms

        Returns:
            Scaling factor (typically 1.0-1.3, default 1.2525 if calculation fails)
        """
        import math

        DEFAULT_FACTOR = 1.2525  # Empirically derived fallback

        try:
            structure = pmd.load_file(pdb_path)
            coords = structure.coordinates
            if coords is None or len(coords) == 0:
                return DEFAULT_FACTOR

            # Get atom centers bounding box
            x_coords = coords[:, 0]
            y_coords = coords[:, 1]
            z_coords = coords[:, 2]

            lx = x_coords.max() - x_coords.min()
            ly = y_coords.max() - y_coords.min()
            lz = z_coords.max() - z_coords.min()

            # Half-box dimensions with buffer
            hx = 0.5 * lx + buffer
            hy = 0.5 * ly + buffer
            hz = 0.5 * lz + buffer

            # Unit vector along diagonal direction
            # Components proportional to INVERSE of half-box dims (from tLEaP source)
            dx = hy * hz
            dy = hx * hz
            dz = hx * hy
            norm = math.sqrt(dx**2 + dy**2 + dz**2)
            if norm < 1e-10:
                return DEFAULT_FACTOR
            ux, uy, uz = dx/norm, dy/norm, dz/norm

            # Find center of atom positions (tLEaP centers the molecule)
            cx = (x_coords.max() + x_coords.min()) / 2
            cy = (y_coords.max() + y_coords.min()) / 2
            cz = (z_coords.max() + z_coords.min()) / 2

            # Find maximum atom distance from center along diagonal direction
            max_diag_extent = 0.0
            for i in range(len(coords)):
                # Position relative to center
                rx = x_coords[i] - cx
                ry = y_coords[i] - cy
                rz = z_coords[i] - cz
                # Project onto diagonal unit vector (using absolute values as tLEaP does)
                diag_proj = abs(rx) * ux + abs(ry) * uy + abs(rz) * uz
                max_diag_extent = max(max_diag_extent, diag_proj)

            # Distance from origin to diagonal face
            face_dist = 0.5 * math.sqrt(hx**2 + hy**2 + hz**2)

            # Required distance along diagonal to maintain buffer
            required_dist = max_diag_extent + buffer

            if required_dist <= face_dist:
                return 1.0  # No scaling needed

            # Scale factor to meet diagonal criterion
            scale_factor = required_dist / face_dist

            logger.debug(f"Diagonal cut factor: {scale_factor:.4f} "
                        f"(max_diag={max_diag_extent:.1f}, face_dist={face_dist:.1f})")

            return scale_factor

        except Exception as e:
            logger.debug(f"Could not calculate diagonal cut factor: {e}, using default")
            return DEFAULT_FACTOR

    def _estimate_solvation_parameters(self, dimensions: dict, buffer: float = 10.0,
                                        use_octahedron: bool = True,
                                        pdb_path: str = None,
                                        custom_forcefields: dict = None) -> dict:
        """
        Estimate solvation box parameters and ion counts using the SPLIT method.

        The SPLIT method (Machado et al.) correctly accounts for ion depletion near
        charged solutes, avoiding overestimation of salt concentration.

        tLEaP's solvateoct algorithm (derived from actual tLEaP output):
        1. Add hydrogens and complete incomplete residues
        2. Calculate vdW bounding box (atom centers + vdW radii, ~2.7 Å/side)
        3. Make box cubic (use max dimension)
        4. Apply diagonal cut scaling: cubic_dim = vdw_max + 2 × buffer × 1.2525
        5. Oct volume = cubic³ × 0.519
        6. Waters ≈ volume / 39.7 Å³

        Args:
            dimensions: Dictionary from _calculate_protein_dimensions()
            buffer: Solvation buffer distance in Angstroms
            use_octahedron: If True, use truncated octahedron; otherwise rectangular box
            pdb_path: Path to PDB file for accurate vdW calculation (optional)
            custom_forcefields: Custom forcefield dict for multi-microstate workflow (optional)

        Returns:
            Dictionary with volume, water count, and ion calculations
        """
        if dimensions is None:
            return None

        lx = dimensions['lx']
        ly = dimensions['ly']
        lz = dimensions['lz']
        Q = dimensions['net_charge']

        # Try to use accurate vdW dimensions from forcefield R* values
        vdw_source = "empirical"
        accurate_vdw = None
        if pdb_path:
            accurate_vdw = self._calculate_vdw_dimensions_with_forcefield(
                pdb_path, custom_forcefields=custom_forcefields
            )

        if accurate_vdw:
            # Use accurate dimensions from forcefield
            vdw_lx, vdw_ly, vdw_lz = accurate_vdw
            h_vdw_correction = None  # Not applicable when using accurate values
            vdw_source = "forcefield"
            logger.info(f"Using accurate vdW dimensions from forcefield R* values")
        else:
            # PDB bounding box is atom centers only. tLEaP adds:
            # 1. Hydrogens (extend ~1.0-1.1 Å beyond heavy atoms)
            # 2. vdW radii (~1.5-1.7 Å for heavy atoms, 0-1.5 Å for H depending on type)
            # Combined effect: ~2.7 Å per side, or 5.4 Å per dimension
            # Note: LJ radii for H atoms range from 0.0 to 1.487 Å depending on atom type
            h_vdw_correction = 5.4  # Å per dimension (empirically derived)
            vdw_lx = lx + h_vdw_correction
            vdw_ly = ly + h_vdw_correction
            vdw_lz = lz + h_vdw_correction

        # Rectangular box volume with buffer (for reference)
        box_x = vdw_lx + 2 * buffer
        box_y = vdw_ly + 2 * buffer
        box_z = vdw_lz + 2 * buffer
        rect_volume = box_x * box_y * box_z

        if use_octahedron:
            # tLEaP makes the box cubic using max vdW dimension
            vdw_max = max(vdw_lx, vdw_ly, vdw_lz)

            # Calculate diagonal cut scaling factor
            # tLEaP scales the BUFFER (not the whole box) to meet the diagonal
            # cut criterion for truncated octahedron geometry
            # Formula: cubic_dim = vdw_max + 2 × buffer × diagonal_cut_factor
            if pdb_path:
                # Calculate actual factor from atom positions
                diagonal_cut_factor = self._calculate_diagonal_cut_factor(pdb_path, buffer)
            else:
                # Use empirical default (derived from typical proteins)
                diagonal_cut_factor = 1.2525

            cubic_dim = vdw_max + 2 * buffer * diagonal_cut_factor

            cubic_volume = cubic_dim ** 3

            # Truncated octahedron is ~51.9% of the enclosing cubic box volume
            # (Empirically derived from tLEaP output: 260734 / 79.49³ = 0.519)
            oct_factor = 0.519
            box_volume = cubic_volume * oct_factor
            box_type = "truncated octahedron"
        else:
            vdw_max = None
            cubic_dim = None
            cubic_volume = None
            diagonal_cut_factor = None
            box_volume = rect_volume
            oct_factor = 1.0
            box_type = "rectangular"

        # Apply conservative scaling factor (0.90) for ion calculation
        # Reasons: (1) H/vdW estimates are approximate, (2) incomplete residues
        # may be built by tLEaP, (3) better to underestimate ions since box
        # shrinks ~12% during equilibration, increasing concentration
        conservative_factor = 0.90
        box_volume_for_ions = box_volume * conservative_factor

        # Estimate number of water molecules
        # From tLEaP output: volume/n_waters ≈ 39.7 Å³ (at ~0.88 g/cc initial density)
        water_volume = 39.7  # Å³ per water molecule
        n_waters = int(box_volume_for_ions / water_volume)

        # SPLIT method for 150 mM NaCl (Machado et al.)
        # No = Nw * M / 56
        # Where 56 ≈ molarity of pure water (55.5 M)
        # This assumes water density of 1 g/mL at 25°C
        molarity = 0.15  # 150 mM NaCl
        water_molarity = 55.5  # More precise value

        No = n_waters * molarity / water_molarity

        # SPLIT the solute charge:
        # N+ = No - Q/2  (cations)
        # N- = No + Q/2  (anions)
        # Note: Q is the NET charge, so if protein is negative (Q < 0),
        # we need MORE cations to neutralize
        n_cations = No - Q / 2
        n_anions = No + Q / 2

        # Round to integers (round up for safety)
        import math
        n_cations_int = math.ceil(n_cations)
        n_anions_int = math.ceil(n_anions)

        # Ensure non-negative
        n_cations_int = max(0, n_cations_int)
        n_anions_int = max(0, n_anions_int)

        # Check SPLIT validity: No >= |Q|
        split_valid = No >= abs(Q) if Q != 0 else True
        split_accurate = No >= 2 * abs(Q) if Q != 0 else True  # <1% error

        # Build result dictionary
        result = {
            'buffer': buffer,
            'box_type': box_type,
            'atom_dimensions': (lx, ly, lz),
            'vdw_dimensions': (vdw_lx, vdw_ly, vdw_lz),
            'vdw_source': vdw_source,  # 'forcefield' or 'empirical'
            'h_vdw_correction': h_vdw_correction,  # None if using forcefield
            'box_dimensions': (box_x, box_y, box_z),
            'rect_volume': rect_volume,
            'box_volume': box_volume,
            'conservative_factor': conservative_factor,
            'box_volume_for_ions': box_volume_for_ions,
            'oct_factor': oct_factor,
            'n_waters': n_waters,
            'water_volume': water_volume,
            'molarity': molarity,
            'No': No,
            'net_charge': Q,
            'n_cations': n_cations_int,
            'n_anions': n_anions_int,
            'n_cations_exact': n_cations,
            'n_anions_exact': n_anions,
            'split_valid': split_valid,
            'split_accurate': split_accurate
        }

        # Add octahedron-specific info if applicable
        if use_octahedron:
            result['vdw_max'] = vdw_max
            result['diagonal_cut_factor'] = diagonal_cut_factor
            result['cubic_dim'] = cubic_dim
            result['cubic_volume'] = cubic_volume

        return result

    def _format_ion_calculation_section(self, dimensions: dict, solvation: dict) -> str:
        """
        Format the ion calculation section for the tLEaP template.

        Shows all calculations transparently so users can understand and verify.
        """
        if dimensions is None or solvation is None:
            return """# ============================================================================
# ION CALCULATION (could not calculate - PDB not available)
# ============================================================================
# To calculate ions for a specific concentration, you need:
#   1. Number of water molecules (Nw) - from tLEaP output after solvation
#   2. Target molarity (M) - e.g., 0.15 M for physiological salt
#   3. Solute net charge (Q)
#
# SPLIT method (Machado et al.):
#   No = Nw * M / 55.5
#   Total Na+ = No - Q/2
#   Total Cl- = No + Q/2
#
# TWO-STAGE ION PLACEMENT (recommended):
#   Stage A: Neutralize with 'addions' (electrostatic placement)
#     - If Q < 0: addions mol Na+ 0  (adds |Q| Na+ near negative charges)
#     - If Q > 0: addions mol Cl- 0  (adds |Q| Cl- near positive charges)
#   Stage B: Add bulk salt with 'addionsrand' (random placement)
#     - addionsrand mol Na+ X Cl- Y  (where X,Y = remaining ions for target conc.)
#
# Reference: http://archive.ambermd.org/202002/0194.html
# ============================================================================
"""

        lx, ly, lz = dimensions['lx'], dimensions['ly'], dimensions['lz']
        Q = solvation['net_charge']
        buffer = solvation['buffer']
        charge_source = dimensions.get('charge_source', 'unknown')
        charge_warnings = dimensions.get('charge_warnings', [])

        # Get vdW-corrected dimensions
        vdw_lx, vdw_ly, vdw_lz = solvation.get('vdw_dimensions', (lx + 5.4, ly + 5.4, lz + 5.4))
        h_vdw_correction = solvation.get('h_vdw_correction', 5.4)
        vdw_source = solvation.get('vdw_source', 'empirical')
        conservative_factor = solvation.get('conservative_factor', 0.90)

        # Format charge source for display
        if charge_source == 'protonation_analyzer':
            charge_note = "(from Protonation State Analyzer)"
        elif charge_source == 'residue_estimate':
            charge_note = "(estimated from residues at pH 7 - run Protonation Analyzer for accuracy)"
        else:
            charge_note = "(source unknown)"

        # Format charge warnings
        warnings_text = ""
        if charge_warnings:
            warnings_text = "\n#\n# ⚠ CHARGE ESTIMATION WARNINGS:\n"
            for warning in charge_warnings:
                warnings_text += f"#   • {warning}\n"

        # Format vdW dimension source
        if vdw_source == 'forcefield':
            vdw_note = "vdW box (from forcefield R* values)"
            vdw_line = f"#   {vdw_note}: {vdw_lx:.1f} × {vdw_ly:.1f} × {vdw_lz:.1f} Å"
        else:
            vdw_note = f"+H/vdW correction (+{h_vdw_correction:.1f} Å empirical)"
            vdw_line = f"#   {vdw_note}: {vdw_lx:.1f} × {vdw_ly:.1f} × {vdw_lz:.1f} Å"

        section = f"""# ============================================================================
# ION CALCULATION FOR 150 mM NaCl (SPLIT method, Machado et al.)
# ============================================================================
# This is an ESTIMATE. tLEaP adds hydrogens and may complete incomplete residues,
# which can change the bounding box. Check tLEaP output and recalculate if needed.
#
# STEP 1: Protein dimensions from PDB
#   Atom centers: {lx:.1f} × {ly:.1f} × {lz:.1f} Å
{vdw_line}
#   Net charge Q = {Q:+d} {charge_note}{warnings_text}
#
# STEP 2: Estimated box volume (tLEaP solvateoct with {buffer:.1f} Å buffer)
"""

        if solvation['box_type'] == "truncated octahedron":
            # Show the tLEaP-accurate octahedron calculation
            vdw_max = solvation.get('vdw_max', max(vdw_lx, vdw_ly, vdw_lz))
            diagonal_cut_factor = solvation.get('diagonal_cut_factor', 1.2525)
            cubic_dim = solvation.get('cubic_dim', vdw_max + 2 * buffer * diagonal_cut_factor)
            cubic_volume = solvation.get('cubic_volume', cubic_dim ** 3)
            section += f"""#   tLEaP formula: cubic_dim = vdw_max + 2 × buffer × {diagonal_cut_factor:.4f}
#                 cubic_dim = {vdw_max:.1f} + 2 × {buffer:.1f} × {diagonal_cut_factor:.4f} = {cubic_dim:.1f} Å
#   Cubic volume: {cubic_dim:.1f}³ = {cubic_volume:.0f} Å³
#   Octahedron factor: ×{solvation['oct_factor']:.3f} (truncated octahedron geometry)
#   Estimated oct volume: {solvation['box_volume']:.0f} Å³
#   Conservative factor: ×{conservative_factor:.2f} (for ion calculation uncertainty)
#   Volume for ion calc: {solvation['box_volume_for_ions']:.0f} Å³
"""
        else:
            box_x, box_y, box_z = solvation['box_dimensions']
            section += f"""#   Box: {box_x:.1f} × {box_y:.1f} × {box_z:.1f} Å
#   Rectangular volume: {solvation['rect_volume']:.0f} Å³
#   Conservative factor: ×{conservative_factor:.2f} (for ion calculation uncertainty)
#   Volume for ion calc: {solvation['box_volume_for_ions']:.0f} Å³
"""

        section += f"""#
# STEP 3: Estimate water molecules (at tLEaP initial density ~0.88 g/cc)
#   Volume per water ≈ {solvation['water_volume']:.1f} Å³
#   Nw ≈ {solvation['box_volume_for_ions']:.0f} / {solvation['water_volume']:.1f} = {solvation['n_waters']} waters
#
# STEP 4: SPLIT method for {solvation['molarity']*1000:.0f} mM NaCl
#   No = Nw × M / 55.5 = {solvation['n_waters']} × {solvation['molarity']} / 55.5 = {solvation['No']:.1f}
#   Total Na+ = No - Q/2 = {solvation['No']:.1f} - ({Q}/2) = {solvation['n_cations_exact']:.1f} → {solvation['n_cations']}
#   Total Cl- = No + Q/2 = {solvation['No']:.1f} + ({Q}/2) = {solvation['n_anions_exact']:.1f} → {solvation['n_anions']}
#
# STEP 5: Two-stage ion placement (recommended procedure)
#   Stage A: Neutralize with 'addions' (electrostatic placement near charged groups)
#   Stage B: Add bulk salt with 'addionsrand' (random placement for ionic strength)
#
"""

        # Add the two-stage breakdown
        if Q < 0:
            n_neutralize = abs(Q)
            n_bulk_cations = max(0, solvation['n_cations'] - n_neutralize)
            n_bulk_anions = solvation['n_anions']
            section += f"""#   Your system (Q = {Q:+d}):
#     A. addions mol Na+ 0      → adds {n_neutralize} Na+ near negative charges
#     B. addionsrand mol Na+ {n_bulk_cations} Cl- {n_bulk_anions}  → adds bulk salt randomly
#
"""
        elif Q > 0:
            n_neutralize = abs(Q)
            n_bulk_cations = solvation['n_cations']
            n_bulk_anions = max(0, solvation['n_anions'] - n_neutralize)
            section += f"""#   Your system (Q = {Q:+d}):
#     A. addions mol Cl- 0      → adds {n_neutralize} Cl- near positive charges
#     B. addionsrand mol Na+ {n_bulk_cations} Cl- {n_bulk_anions}  → adds bulk salt randomly
#
"""
        else:
            section += f"""#   Your system (Q = 0, neutral):
#     No neutralization needed - just add bulk salt randomly:
#     addionsrand mol Na+ {solvation['n_cations']} Cl- {solvation['n_anions']}
#
"""

        # Add validity warnings
        if not solvation['split_valid']:
            section += f"""# ⚠ WARNING: SPLIT may be inaccurate (No={solvation['No']:.1f} < |Q|={abs(Q)})
#   Consider using SLTCAP server for precise calculation.
#
"""
        elif not solvation['split_accurate']:
            section += f"""# NOTE: SPLIT approximation (<5% error, No < 2×|Q|)
#
"""

        section += """# APPROXIMATIONS IN THIS ESTIMATE:
#   • H/vdW correction: +5.4 Å/dimension (H atom LJ radii range 0-1.49 Å by type)
#   • Diagonal cut factor: 1.2525 (tLEaP scales buffer for octahedron geometry)
#   • Octahedron factor: 0.519 (empirical ratio of oct to cubic volume)
#   • Conservative factor: 0.90 (reduces estimate to avoid over-counting ions)
#   • tLEaP may build atoms for incomplete residues, changing the bounding box
#
# DENSITY AND EQUILIBRATION:
#   • tLEaP solvates on a grid less dense than liquid water
#   • NPT equilibration contracts the box ~5-10% to reach ~1.0 g/cc
#   • This INCREASES concentration: a 150 mM target rises ~5-10% after equilibration
#   • The conservative factor partially compensates for this
#
# FOR PRECISE WORK:
#   • Check tLEaP output for actual water count ("Added X residues")
#   • Recalculate: No = Nw × M / 55.5, Na+ = No - Q/2, Cl- = No + Q/2
#
# Reference: Machado et al., AMBER mailing list (2020)
#   http://archive.ambermd.org/202002/0194.html
# ============================================================================
"""
        return section

    def _format_ion_calculation_section_from_info_pass(self, ion_counts: dict) -> str:
        """
        Format the ion calculation section using actual values from tLEaP info pass.

        This is more accurate than the vdW estimation method because it uses
        the real water count and charge from tLEaP itself.

        Args:
            ion_counts: Dictionary from _calculate_ions_from_water_count()

        Returns:
            Formatted ion calculation section for the template
        """
        n_waters = ion_counts['n_waters']
        Q = ion_counts['net_charge']
        n_salt = ion_counts['n_salt_pairs']
        n_cations = ion_counts['n_cations']
        n_anions = ion_counts['n_anions']
        n_neutralize = ion_counts['n_neutralize']
        n_bulk_cations = ion_counts['n_bulk_cations']
        n_bulk_anions = ion_counts['n_bulk_anions']
        molarity = ion_counts['target_molarity']

        section = f"""# ============================================================================
# ION CALCULATION FOR {molarity*1000:.0f} mM NaCl (SPLIT method, Machado et al.)
# ============================================================================
# Values calculated from tLEaP info pass (ACCURATE - using actual tLEaP values)
#
# STEP 1: System info from tLEaP
#   Water molecules (Nw) = {n_waters}
#   Net charge (Q) = {Q:+d}
#
# STEP 2: SPLIT method calculation
#   No = Nw × M / 55.5 = {n_waters} × {molarity} / 55.5 = {n_salt}
#   Total Na+ = No - Q/2 = {n_salt} - ({Q}/2) = {n_cations}
#   Total Cl- = No + Q/2 = {n_salt} + ({Q}/2) = {n_anions}
#
# STEP 3: Two-stage ion placement
"""

        if Q < 0:
            section += f"""#   Stage A: Neutralize with 'addions'
#     addions mol Na+ 0  → adds {n_neutralize} Na+ near negative charges
#   Stage B: Add bulk salt with 'addionsrand'
#     addionsrand mol Na+ {n_bulk_cations} Cl- {n_bulk_anions}  → bulk salt
#
"""
        elif Q > 0:
            section += f"""#   Stage A: Neutralize with 'addions'
#     addions mol Cl- 0  → adds {n_neutralize} Cl- near positive charges
#   Stage B: Add bulk salt with 'addionsrand'
#     addionsrand mol Na+ {n_bulk_cations} Cl- {n_bulk_anions}  → bulk salt
#
"""
        else:
            section += f"""#   System is neutral - just add bulk salt:
#     addionsrand mol Na+ {n_cations} Cl- {n_anions}
#
"""

        section += """# NOTE: tLEaP solvates on a grid less dense than liquid water.
#       NPT equilibration contracts the box ~5-10% to reach ~1.0 g/cc.
#       This INCREASES concentration: a 150 mM target rises ~5-10% after equilibration.
#
# Reference: Machado et al., AMBER mailing list (2020)
#   http://archive.ambermd.org/202002/0194.html
# ============================================================================
"""
        return section

    def _build_standard_forcefield_section(self) -> tuple:
        """
        Build the standard forcefield source commands from user selections.

        Returns:
            Tuple of (forcefield_lines: str, water_box_type: str)
        """
        selected = self.get_from_workspace("selected_standard_forcefields", {})

        # Fallback if no selections made (shouldn't happen in normal workflow)
        if not selected:
            logger.warning("No standard forcefields selected, using defaults")
            return (
                "# Default forcefields (no selection made)\n"
                "source leaprc.protein.ff19SB\n"
                "source leaprc.gaff2\n"
                "source leaprc.water.opc",
                "OPCBOX"
            )

        lines = []
        water_box = "OPCBOX"  # Default

        # Order matters for AMBER: protein -> modified AA -> nucleic -> carbs -> gaff -> water
        # AMBER load order: protein -> modified AA -> nucleic acid -> carbs -> gaff -> water -> ions
        category_order = ['protein', 'modified_aa', 'dna', 'rna', 'carbohydrates',
                         'lipids', 'small_molecules', 'water', 'ions']

        for category in category_order:
            sel = selected.get(category)
            if sel is None:
                continue

            if isinstance(sel, list):
                # Multi-select category (e.g., 'modified_aa')
                for item in sel:
                    if 'leaprc' in item:
                        lines.append(f"source {item['leaprc']}")
            elif isinstance(sel, dict):
                if 'leaprc' in sel:
                    leaprc_val = sel['leaprc']
                    if isinstance(leaprc_val, list):
                        for leaprc in leaprc_val:
                            lines.append(f"source {leaprc}")
                    else:
                        lines.append(f"source {leaprc_val}")
                if 'frcmod' in sel and sel['frcmod']:
                    frcmod_val = sel['frcmod']
                    if isinstance(frcmod_val, list):
                        for frcmod_file in frcmod_val:
                            lines.append(f"loadamberparams {frcmod_file}")
                    else:
                        lines.append(f"loadamberparams {frcmod_val}")
                # Track water box type
                if category == 'water':
                    box_val = sel.get('box')
                    # A water FF selected during MCPB preprocessing carries
                    # box='none' (the metal-free preprocessing build adds no
                    # box). If we are emitting an explicit solvation box, that
                    # placeholder must resolve to the model's real solvent unit
                    # (e.g. tip3p -> TIP3PBOX), otherwise tleap gets
                    # `solvateBox mol none ...` and fails.
                    if box_val and str(box_val).lower() != 'none':
                        water_box = box_val
                    else:
                        water_box = self._default_water_box(sel)

        return '\n'.join(lines), water_box

    @staticmethod
    def _default_water_box(water_sel: dict) -> str:
        """Map a water FF selection to its AMBER solvent box unit name.

        Used when the stored box is missing or 'none' (e.g. a water model chosen
        during MCPB preprocessing, which deliberately records no box).
        """
        name = (water_sel.get('name') or '') if isinstance(water_sel, dict) else ''
        leaprc = (water_sel.get('leaprc') or '') if isinstance(water_sel, dict) else ''
        token = (str(name) + ' ' + str(leaprc)).lower()
        # Order matters: check longer/more-specific names before their prefixes
        # (tip4pew before tip4p, opc3 before opc).
        mapping = [
            ('tip4pew', 'TIP4PEWBOX'),
            ('tip4p', 'TIP4PBOX'),
            ('tip5p', 'TIP5PBOX'),
            ('opc3', 'OPC3BOX'),
            ('opc', 'OPCBOX'),
            ('spce', 'SPCBOX'),
            ('spc', 'SPCBOX'),
            ('fb3', 'TIP3PFBOX'),
            ('fb4', 'TIP4PFBOX'),
            ('tip3p', 'TIP3PBOX'),
        ]
        for key, box in mapping:
            if key in token:
                return box
        return 'TIP3PBOX'  # safe, widely-compatible default

    def _generate_tleap_template(self, solvation_params: dict = None):
        """
        Generate the tLEaP template using selected forcefields and solvation parameters.

        For explicit solvent: runs info pass to get accurate ion counts.
        For implicit solvent: no solvation section in template.

        Args:
            solvation_params: Dict from _configure_solvation_parameters()
        """
        console = self.processor.console

        if solvation_params is None:
            solvation_params = self.get_from_workspace("solvation_parameters", {})

        # Safety net: refuse to emit a tleap script if the final set of frcmod /
        # lib files-to-load still contains atom-type parameter collisions. The
        # picker-side resolver runs at FF selection time and should already have
        # rewritten any colliding sets, but this catches paths that bypass the
        # picker (e.g., MCPB-derived FFs from preprocessing that collide with a
        # bundled cofactor FF) or out-of-band workspace edits.
        self._verify_no_ff_collisions_or_abort()

        # Use the same priority logic as _select_priority_pdb_file (silent to avoid duplicate messages)
        selected_pdb = self._select_priority_pdb_file(silent=True)
        if selected_pdb:
            pdb_path = os.path.abspath(selected_pdb)
            pdb_basename = os.path.basename(selected_pdb)
            output_prefix = os.path.splitext(pdb_basename)[0]
        else:
            pdb_path = "STRUCTURE_FILE.pdb"
            output_prefix = "system"

        # Build forcefield section from user selections
        forcefield_section, water_box = self._build_standard_forcefield_section()

        # Build solvation + ion section based on solvent model
        solvation_section = ""
        is_membrane = solvation_params.get('solvent_model') == 'membrane_pre_solvated'
        is_explicit = solvation_params.get('solvent_model') == 'explicit'

        if is_membrane:
            # Membrane system: already solvated by packmol-memgen.
            # Following packmol-memgen's own tLEaP convention:
            #   1. Set box from vdW surface (packmol-memgen's default behavior;
            #      explicit dimensions are only used with --dims/tight-box mode)
            #   2. Neutralization pass: packmol-memgen adds ions for target salt
            #      concentration, but this doesn't guarantee exact charge neutrality.
            #      addionsrand ... 0 adds the minimum extra ions to reach net zero.
            membrane_config = self.get_from_workspace("membrane_config", {})

            solvation_section = "# === PERIODIC BOX & NEUTRALIZATION (membrane system, pre-solvated) ===\n"
            solvation_section += "setBox mol vdw\n"

            # Neutralization pass — same as packmol-memgen's tLEaP script
            cation = membrane_config.get("cation", "K+")
            anion = membrane_config.get("anion", "Cl-")
            solvation_section += f"addionsrand mol {cation} 0\n"
            solvation_section += f"addionsrand mol {anion} 0\n"

            console.print("[grey50]Membrane system: vdW box + neutralization pass[/grey50]")
        elif is_explicit:
            use_octahedron = solvation_params.get('use_octahedron', True)
            buffer = solvation_params.get('buffer', 10.0)
            buffer_xyz = solvation_params.get('buffer_xyz')
            oct_diagonal = solvation_params.get('oct_diagonal', 0.0)
            iso = solvation_params.get('iso', False)
            salts, neutralize_index = _normalize_salts(solvation_params)
            solvate_cmd = "solvateoct" if use_octahedron else "solvateBox"
            buf_arg = _format_buffer_for_tleap(buffer, buffer_xyz, oct_diagonal,
                                                use_octahedron, iso)

            # Solvation command
            solvation_section = f"# === SOLVATION ===\n{solvate_cmd} mol {water_box} {buf_arg}\n"

            # Run info pass for accurate ion counts
            ion_counts = None
            if selected_pdb and os.path.exists(pdb_path):
                console.print("\n[bold yellow]Running tLEaP info pass for accurate ion counts...[/bold yellow]")
                console.print("[grey50]This runs tLEaP with early termination to get real water count and charge[/grey50]")

                system_info = self._run_info_pass_for_template(
                    pdb_path, water_box=water_box, buffer=buffer,
                    use_octahedron=use_octahedron, quiet=False,
                    buffer_xyz=buffer_xyz, oct_diagonal=oct_diagonal, iso=iso
                )

                n_waters = system_info.get('n_waters')
                net_charge = system_info.get('net_charge')

                if n_waters is not None and net_charge is not None:
                    ion_counts = self._calculate_multi_salt_ions(
                        n_waters, int(round(net_charge)),
                        salts=salts, neutralize_index=neutralize_index
                    )
                    self._display_ion_calculation(ion_counts)
                else:
                    console.print("[yellow]Could not get water count/charge from tLEaP info pass.[/yellow]")

            # Build ion commands
            solvation_section += "\n" + self._build_ion_commands_for_template(ion_counts)
        else:
            console.print("[grey50]Implicit solvent selected - no solvation section in template[/grey50]")

        template = f"""# ProPrep-generated tLEaP Input File

# === STANDARD FORCEFIELDS ===
{forcefield_section}

# === CUSTOM ATOM TYPES (auto-filled by ProPrep) ===
# ATOM_TYPES_SECTION

# === CUSTOM FORCEFIELD PARAMETERS (auto-filled by ProPrep) ===
# FORCEFIELD_PARAMETERS_SECTION

# === LOAD STRUCTURE ===
mol = loadpdb {pdb_path}

# === BOND DEFINITIONS (auto-filled by ProPrep) ===
# BOND_DEFINITIONS_SECTION

{solvation_section}
# === VALIDATION ===
check mol

# === OUTPUT ===
saveamberparm mol {output_prefix}.prmtop {output_prefix}.rst7

quit"""

        return template

    def _display_ion_calculation(self, ion_counts: dict, target_molarity: float = None,
                                   cation: dict = None, anion: dict = None):
        """Display the ion calculation with educational explanation in a rich panel.

        Handles both single-salt (preserves the original SPLIT-method
        narrative) and multi-salt (a per-salt breakdown).
        """
        from rich.panel import Panel

        console = self.processor.console
        salts = ion_counts.get('salts') or []

        lines = []
        if not salts:
            self._build_neutralize_only_panel(lines, ion_counts.get('net_charge', 0),
                                                ion_counts.get('n_neutralize_cation', 0)
                                                + ion_counts.get('n_neutralize_anion', 0),
                                                cation, anion)
            title = "[bold blue]Ion Calculation -- Neutralize Only[/bold blue]"
        elif len(salts) == 1 and salts[0]['concentration'] > 0:
            self._build_split_method_panel(lines, ion_counts, salts[0])
            conc_mM = salts[0]['concentration'] * 1000
            label = _salt_label(salts[0]['cation'], salts[0]['anion'])
            title = f"[bold blue]Ion Calculation -- {conc_mM:.0f} mM {label}, SPLIT Method[/bold blue]"
        elif len(salts) == 1 and salts[0]['concentration'] == 0:
            n_neut = ion_counts.get('n_neutralize_cation', 0) + ion_counts.get('n_neutralize_anion', 0)
            self._build_neutralize_only_panel(lines, ion_counts['net_charge'], n_neut,
                                                salts[0]['cation'], salts[0]['anion'])
            title = "[bold blue]Ion Calculation -- Neutralize Only[/bold blue]"
        else:
            self._build_multi_salt_panel(lines, ion_counts)
            title = "[bold blue]Ion Calculation -- Multi-Salt[/bold blue]"

        body = "\n".join(lines)
        panel = Panel(body, title=title, border_style="blue",
                      expand=True, padding=(1, 2))
        console.print(panel)

    def _build_split_method_panel(self, lines, ion_counts, salt):
        """Educational SPLIT method panel for the single-salt case."""

        n_waters = ion_counts['n_waters']
        Q = ion_counts['net_charge']
        n_neut_cat = ion_counts.get('n_neutralize_cation', 0)
        n_neut_an = ion_counts.get('n_neutralize_anion', 0)
        n_neutralize = max(n_neut_cat, n_neut_an)
        cation, anion = salt['cation'], salt['anion']
        cl, al = cation['label'], anion['label']
        cu, au = cation['unit'], anion['unit']
        n_salt = salt.get('n_pairs', 0)
        n_bulk_cations = salt.get('n_bulk_cations', 0)
        n_bulk_anions = salt.get('n_bulk_anions', 0)
        n_cations = n_bulk_cations + n_neut_cat
        n_anions = n_bulk_anions + n_neut_an
        target_molarity = salt['concentration']
        conc_mM = target_molarity * 1000
        salt_label = _salt_label(cation, anion)
        symmetric = (abs(cation['charge']) == 1 and abs(anion['charge']) == 1)

        # --- What and why ---
        lines.append("[bold blue]Why add ions?[/bold blue]")
        lines.append("")
        lines.append(
            f"Your system has net charge [bold dark_orange3]{Q:+d}[/bold dark_orange3]; PME with periodic "
            "boundaries requires a net-neutral box, so counterions are added. Adding "
            "salt also sets the bulk ionic strength, which affects electrostatic "
            "screening, intermolecular interactions, and conformational stability.")
        lines.append("")

        # --- The method ---
        lines.append("[bold blue]How: the SPLIT method[/bold blue]")
        lines.append(
            "Machado, M. R. & Pantano, S. \"Split the charge difference in two! "
            "A rule of thumb for adding proper amounts of ions in MD simulations.\" "
            "J. Chem. Theory Comput. 16(3), 1367-1372 (2020)")
        lines.append("")
        lines.append(
            "Neutralizing first and then adding salt overcounts ions and gives the "
            "wrong bulk concentration. SPLIT distributes the system charge across "
            "both ion types in one pass, so the counts neutralize AND hit the target "
            "at once:")
        lines.append("")

        n_salt_exact = n_waters * target_molarity / 55.5
        lines.append(f"  [bold blue]Step 1[/bold blue]  water count (tLEaP info pass):  Nw = [bold dark_orange3]{n_waters}[/bold dark_orange3]")
        lines.append(f"  [bold blue]Step 2[/bold blue]  formula units:  No = Nw · M / 55.5 = {n_waters} · {target_molarity} / 55.5 = {n_salt_exact:.1f} -> [bold dark_orange3]{n_salt}[/bold dark_orange3]")
        lines.append("          (M = your target molarity; 55.5 mol/L = the molarity of pure water)")
        if symmetric:
            n_cat_exact = n_salt_exact - Q / 2.0
            lines.append(f"  [bold blue]Step 3[/bold blue]  split charge Q = {Q:+d}:  {cl} = No - (Q / 2) = {n_salt_exact:.1f} - ({Q} / 2) = {n_cat_exact:.1f} -> [bold dark_orange3]{n_cations}[/bold dark_orange3]")
            lines.append(f"          {al} = {cl} + Q = {n_cations} + ({Q:+d}) = [bold dark_orange3]{n_anions}[/bold dark_orange3]   (guarantees exact neutrality)")
        else:
            a, b = _salt_formula_units(cation['charge'], anion['charge'])
            lines.append(f"  [bold blue]Step 3[/bold blue]  bulk at {a}:{b} ratio + neutralize:  {a} · {n_salt} = [bold dark_orange3]{n_bulk_cations}[/bold dark_orange3] {cl}, {b} · {n_salt} = [bold dark_orange3]{n_bulk_anions}[/bold dark_orange3] {al}")
            lines.append(f"          neutralizers added on the deficit side to satisfy Q = {Q:+d}   (non-1:1 salt)")
        lines.append("")

        # --- Two-stage placement ---
        lines.append("[bold blue]Ion placement[/bold blue]  (two stages in tLEaP)")
        lines.append("")
        lines.append("  [bold blue]A. Neutralize[/bold blue] with addions -- counterions at lowest-potential sites")
        if Q < 0:
            lines.append(f"     addions mol {cu} 0   --> {n_neutralize} {cl} near the negative charges")
        elif Q > 0:
            lines.append(f"     addions mol {au} 0   --> {n_neutralize} {al} near the positive charges")
        else:
            lines.append("     system is neutral -- no neutralization needed")
        lines.append("  [bold blue]B. Bulk salt[/bold blue] with addionsrand -- remaining ions at random solvent sites")
        lines.append(f"     addionsrand mol {cu} {n_bulk_cations} {au} {n_bulk_anions}")
        lines.append("")

        # --- Result summary ---
        ion_charge = n_cations * cation['charge'] + n_anions * anion['charge']
        lines.append(
            f"[bold blue]Result:[/bold blue] total charge = system "
            f"([bold dark_orange3]{Q:+d}[/bold dark_orange3]) + ions "
            f"([bold dark_orange3]{ion_charge:+d}[/bold dark_orange3]) = 0  ->  "
            f"[bold dark_orange3]{n_cations}[/bold dark_orange3] {cl} + "
            f"[bold dark_orange3]{n_anions}[/bold dark_orange3] {al} "
            f"at ~{conc_mM:.0f} mM {salt_label}")
        lines.append("")

        # --- Caveats ---
        lines.append("[bold blue]Caveat: Concentration After Equilibration[/bold blue]")
        lines.append("")
        lines.append(
            "tLEaP solvates on a grid that is less dense than liquid water. During "
            "NPT equilibration the box contracts ~5-10% to reach liquid density "
            "(~1.0 g/cc), which raises the effective salt concentration: a "
            f"{conc_mM:.0f} mM target becomes roughly {conc_mM * 1.05:.0f}-"
            f"{conc_mM * 1.10:.0f} mM. This is standard practice and generally "
            "acceptable for most applications.")

    def _build_multi_salt_panel(self, lines, ion_counts):
        """Per-salt breakdown panel for multi-salt mixes."""

        salts = ion_counts['salts']
        Q = ion_counts['net_charge']
        n_waters = ion_counts['n_waters']
        n_neut_cat = ion_counts.get('n_neutralize_cation', 0)
        n_neut_an = ion_counts.get('n_neutralize_anion', 0)
        idx = ion_counts.get('neutralize_index', 0)
        neut_salt = salts[idx]
        nc, na = neut_salt['cation'], neut_salt['anion']

        lines.append("[bold blue]Why a multi-salt mix?[/bold blue]")
        lines.append("")
        lines.append(f"Your system has net charge [bold dark_orange3]{Q:+d}[/bold dark_orange3] and {n_waters} waters. You're adding "
                     f"{len(salts)} salts simultaneously; each contributes its own bulk ions at "
                     f"its target concentration.")
        lines.append("")

        lines.append("[bold blue]Per-salt bulk (a:b formula ratio at target M)[/bold blue]")
        lines.append("")
        for i, s in enumerate(salts):
            cl, al = s['cation']['label'], s['anion']['label']
            conc_mM = s['concentration'] * 1000
            lbl = _salt_label(s['cation'], s['anion'])
            n_pair = s.get('n_pairs', 0)
            n_bc = s.get('n_bulk_cations', 0)
            n_ba = s.get('n_bulk_anions', 0)
            tag = "  [bold dark_orange3](neutralizer)[/bold dark_orange3]" if i == idx else ""
            lines.append(f"  {i+1}. {lbl} @ {conc_mM:.0f} mM -> {n_pair} formula units"
                         f" ([bold dark_orange3]{n_bc}[/bold dark_orange3] {cl} + [bold dark_orange3]{n_ba}[/bold dark_orange3] {al}){tag}")
        lines.append("")

        lines.append(f"[bold blue]Charge neutralization (via {_salt_label(nc, na)})[/bold blue]")
        lines.append("")
        if Q == 0:
            lines.append("  System is already neutral; no extra counterions needed.")
        elif n_neut_cat and n_neut_an:
            lines.append(f"  addions mol {nc['unit']} {n_neut_cat} {na['unit']} {n_neut_an}")
            lines.append(f"  --> {n_neut_cat} {nc['label']} + {n_neut_an} {na['label']} "
                         f"to absorb Q={Q:+d}")
        elif n_neut_cat:
            unit = nc['unit']
            lbl = nc['label']
            lines.append(f"  addions mol {unit} 0")
            lines.append(f"  --> places {n_neut_cat} {lbl} near negative charges")
        elif n_neut_an:
            unit = na['unit']
            lbl = na['label']
            lines.append(f"  addions mol {unit} 0")
            lines.append(f"  --> places {n_neut_an} {lbl} near positive charges")
        lines.append("")

        # Aggregated totals across all salts + neutralizer
        agg = {}
        def _add(key_unit, key_label, charge, count):
            entry = agg.setdefault(key_unit, {'label': key_label, 'charge': charge, 'count': 0})
            entry['count'] += count
        for s in salts:
            _add(s['cation']['unit'], s['cation']['label'], s['cation']['charge'],
                 s.get('n_bulk_cations', 0))
            _add(s['anion']['unit'], s['anion']['label'], s['anion']['charge'],
                 s.get('n_bulk_anions', 0))
        _add(nc['unit'], nc['label'], nc['charge'], n_neut_cat)
        _add(na['unit'], na['label'], na['charge'], n_neut_an)

        lines.append("[bold blue]Total ions in box[/bold blue]")
        lines.append("")
        check = Q
        for unit, info in sorted(agg.items(), key=lambda kv: -kv[1]['charge']):
            if info['count'] == 0:
                continue
            lines.append(f"  {info['label']:>5}: [bold dark_orange3]{info['count']}[/bold dark_orange3]  (z = {info['charge']:+d})")
            check += info['count'] * info['charge']
        if check != 0:
            lines.append(f"  [red]Charge balance check failed: residual {check:+d}[/red]")
        else:
            lines.append("")
            lines.append("  Charge balance verified: total = 0")

    def _build_neutralize_only_panel(self, lines, Q, n_neutralize, cation=None, anion=None):
        """Build the explanation panel for neutralize-only (0 mM salt)."""

        cation = cation or _counter_ion_default("Na+", "cation")
        anion = anion or _counter_ion_default("Cl-", "anion")
        cl, al = cation['label'], anion['label']
        cu, au = cation['unit'], anion['unit']
        salt_label = _salt_label(cation, anion)

        lines.append("[bold blue]Why add counterions?[/bold blue]")
        lines.append("")
        lines.append(
            f"Your system has a net charge of [bold dark_orange3]{Q:+d}[/bold dark_orange3]. MD simulations with "
            "periodic boundary conditions and PME electrostatics require a net-neutral "
            "box. Counterions are added to exactly cancel the system charge.")
        lines.append("")
        lines.append(
            "No additional salt is being added (0 mM). If you want physiological ionic "
            f"strength, re-run with a non-zero salt concentration (e.g. 150 mM {salt_label}).")
        lines.append("")

        if Q < 0:
            lines.append("[bold blue]Ion Placement:[/bold blue]")
            lines.append("")
            lines.append(f"  addions mol {cu} 0")
            lines.append(f"  --> places {n_neutralize} {cl} near negative charges (electrostatic placement)")
        elif Q > 0:
            lines.append("[bold blue]Ion Placement:[/bold blue]")
            lines.append("")
            lines.append(f"  addions mol {au} 0")
            lines.append(f"  --> places {n_neutralize} {al} near positive charges (electrostatic placement)")
        else:
            lines.append("System is neutral -- no counterions needed.")

    def _build_ion_commands_for_template(self, ion_counts: dict, target_molarity: float = None,
                                          cation: dict = None, anion: dict = None) -> str:
        """Build tLEaP ion commands from a multi-salt ion_counts dict.

        `target_molarity`, `cation`, `anion` are accepted for back-compat but
        ignored when ion_counts already contains the salts list. Emits one
        addions-style neutralizer command (using the chosen salt's ions)
        followed by one addionsrand per salt with non-zero bulk content.
        """
        if ion_counts is None:
            cation = cation or _counter_ion_default("Na+", "cation")
            anion = anion or _counter_ion_default("Cl-", "anion")
            return (f"# === ADD IONS ===\n"
                    f"addions mol {cation['unit']} 0\n"
                    f"addions mol {anion['unit']} 0\n")

        salts = ion_counts.get('salts') or []
        Q = ion_counts['net_charge']
        n_waters = ion_counts['n_waters']
        n_neut_cat = ion_counts.get('n_neutralize_cation', abs(Q) if Q < 0 else 0)
        n_neut_an = ion_counts.get('n_neutralize_anion', abs(Q) if Q > 0 else 0)
        idx = ion_counts.get('neutralize_index', 0)

        if salts and 0 <= idx < len(salts):
            neut_salt = salts[idx]
        else:
            neut_salt = {
                'cation': cation or _counter_ion_default("Na+", "cation"),
                'anion': anion or _counter_ion_default("Cl-", "anion"),
            }
        nc, na = neut_salt['cation'], neut_salt['anion']
        nc_u, na_u = nc['unit'], na['unit']
        nc_l, na_l = nc['label'], na['label']
        symmetric_neut = (abs(nc['charge']) == 1 and abs(na['charge']) == 1)

        lines = []
        summary = _salts_summary(salts) if salts else "neutralize only"
        lines.append(f"# === ADD IONS ({summary}, {n_waters} waters, Q={Q:+d}) ===")

        # --- Step 1: neutralize protein charge ---
        if Q != 0:
            if symmetric_neut:
                neut_unit = nc_u if Q < 0 else na_u
                neut_label = nc_l if Q < 0 else na_l
                count = n_neut_cat if Q < 0 else n_neut_an
                lines.append(f"# Neutralize {Q:+d} charge with {count} {neut_label} (from {_salt_label(nc, na)})")
                lines.append(f"addions mol {neut_unit} 0")
            else:
                lines.append(f"# Neutralize {Q:+d} charge with {n_neut_cat} {nc_l} + {n_neut_an} {na_l} (from {_salt_label(nc, na)})")
                if n_neut_cat and n_neut_an:
                    lines.append(f"addions mol {nc_u} {n_neut_cat} {na_u} {n_neut_an}")
                elif n_neut_cat:
                    lines.append(f"addions mol {nc_u} {n_neut_cat}")
                elif n_neut_an:
                    lines.append(f"addions mol {na_u} {n_neut_an}")

        # --- Step 2: per-salt bulk (one addionsrand per salt) ---
        any_bulk = False
        for s in salts:
            n_bc = s.get('n_bulk_cations', 0)
            n_ba = s.get('n_bulk_anions', 0)
            if not (n_bc or n_ba):
                continue
            any_bulk = True
            cu, au = s['cation']['unit'], s['anion']['unit']
            cl, al = s['cation']['label'], s['anion']['label']
            conc_mM = s.get('concentration', 0.0) * 1000
            lbl = _salt_label(s['cation'], s['anion'])
            lines.append(f"# Bulk {lbl} ({conc_mM:.0f} mM): {n_bc} {cl} + {n_ba} {al}")
            lines.append(f"addionsrand mol {cu} {n_bc} {au} {n_ba}")

        if Q == 0 and not any_bulk:
            lines.append("# System is neutral, no ions needed")

        return "\n".join(lines) + "\n"

    def _save_tleap_template(self, template):
        """Save the tLEaP template to workspace"""
        self.update_workspace("tleap_template", template)
        self.processor.console.print("[green]tLEaP template saved[/green]")

    def _edit_tleap_template(self, current_template):
        """Allow user to edit the tLEaP template"""
        import tempfile
        import subprocess
        import os
        
        # Create temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tleap', delete=False) as tmp:
            tmp.write(current_template)
            tmp_path = tmp.name
        
        try:
            # Get user's preferred editor
            editor = os.environ.get('EDITOR', 'nano')
            
            self.processor.console.print(f"[yellow]Opening template in {editor}...[/yellow]")
            self.processor.console.print("[grey50]Save and exit the editor when done.[/grey50]")
            
            # Open editor
            subprocess.run([editor, tmp_path], check=True, stdin=sys.stdin, stdout=sys.stdout, stderr=sys.stderr)
            
            # Read back the edited content
            with open(tmp_path, 'r') as f:
                edited_template = f.read()
            
            # Save the edited template
            self._save_tleap_template(edited_template)

            self.processor.console.print("[green]Template updated successfully![/green]")
        
        except subprocess.CalledProcessError:
            self.processor.console.print("[red]Editor was cancelled or failed[/red]")
            # Fall back to the original template
            self._save_tleap_template(current_template)
        
        except FileNotFoundError:
            self.processor.console.print(f"[red]Editor '{editor}' not found. Using simple text input instead.[/red]")
            # Fall back to simple text input
            self._simple_text_edit(current_template)
        
        finally:
            # Clean up temp file
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)

    def _simple_text_edit(self, current_template):
        """Simple line-by-line text editing fallback"""
        self.processor.console.print("[yellow]Simple template editor[/yellow]")
        self.processor.console.print("[grey50]Enter 'done' on a line by itself to finish[/grey50]")
        
        lines = current_template.split('\n')
        self.processor.console.print(f"[grey50]Current template has {len(lines)} lines[/grey50]")
        
        if confirm_with_context(
            self.processor,
            "Edit template line by line?",
            default=False,
            module="Topology Generator",
            description="Enter line-by-line template editor",
        ):
            new_lines = []
            for i, line in enumerate(lines, 1):
                self.processor.console.print(f"[grey50]Line {i}: {line}[/grey50]")
                new_line = prompt_with_context(
                    self.processor,
                    f"Edit line {i} (or press Enter to keep as-is)",
                    default=line,
                    module="Topology Generator",
                    description=f"Edit template line {i}",
                )
                new_lines.append(new_line)
            
            edited_template = '\n'.join(new_lines)
            self._save_tleap_template(edited_template)
        else:
            # Just save the current template
            self._save_tleap_template(current_template)

    # =========================================================================
    # Standard Forcefield Selection UX
    # =========================================================================

    # Forcefield options organized by category (mirrors the comprehensive template guide)
    # Based on AMBER Manual Chapter 3 and standard leaprc files
    # Curated FF catalog now lives in a shared module so the Force Field
    # Explorer can reuse it for display enrichment (see
    # proprep.forcefield_params.forcefield_catalog). Kept as a class
    # attribute for backward compatibility: structure_preprocessor and
    # membrane_builder read TLeapInputGenerator.FORCEFIELD_OPTIONS.
    FORCEFIELD_OPTIONS = _FF_CATALOG

    def _get_membrane_forcefields(self) -> dict:
        """
        Get forcefields for a membrane system.

        First checks if the membrane builder already stored a comprehensive
        selection in workspace (via its own FF menu). Falls back to
        reconstructing a basic selection from membrane_config if not.
        """
        console = self.processor.console

        # Check if membrane builder already stored full FF selection
        existing = self.get_from_workspace("selected_standard_forcefields", None)
        if existing:
            # Augment lipid leaprc with membrane-specific requirements (lipid_ext, etc.)
            leaprc_reqs = self.get_from_workspace("membrane_leaprc_requirements", [])
            lip_sel = existing.get("lipids")
            if lip_sel and isinstance(lip_sel, dict) and leaprc_reqs:
                lip_leaprc = lip_sel.get("leaprc", "")
                if isinstance(lip_leaprc, str):
                    lip_leaprc = [lip_leaprc]
                for req in leaprc_reqs:
                    if req not in lip_leaprc:
                        lip_leaprc.append(req)
                existing["lipids"] = dict(lip_sel, leaprc=lip_leaprc)
                self.update_workspace("selected_standard_forcefields", existing)

            console.print("[bold green]Membrane system — using forcefields selected in membrane builder:[/bold green]")
            for cat, sel in existing.items():
                if sel is None:
                    continue
                elif isinstance(sel, list):
                    if sel:
                        names = [s['name'] for s in sel]
                        console.print(f"  {cat}: {', '.join(names)}")
                elif isinstance(sel, dict):
                    console.print(f"  {cat}: {sel['name']}")
            console.print("[grey50]Solvation skipped — system already solvated by membrane builder[/grey50]")
            return existing

        # Fallback: reconstruct from membrane_config (old-style ffprot/fflip only)
        membrane_config = self.get_from_workspace("membrane_config", {})
        leaprc_reqs = self.get_from_workspace("membrane_leaprc_requirements", [])

        ffprot = membrane_config.get("ffprot", "ff14SB")
        fflip = membrane_config.get("fflip", "lipid21")
        water_model = membrane_config.get("effective_water_model", "tip3p")

        prot_map = {
            "ff19SB": {"name": "ff19SB", "leaprc": "leaprc.protein.ff19SB"},
            "ff14SB": {"name": "ff14SB", "leaprc": "leaprc.protein.ff14SB"},
            "ff15ipq": {"name": "ff15ipq", "leaprc": "leaprc.protein.ff15ipq"},
            "fb15": {"name": "fb15", "leaprc": "leaprc.protein.fb15"},
        }
        protein_ff = prot_map.get(ffprot, {"name": ffprot, "leaprc": f"leaprc.protein.{ffprot}"})

        lipid_leaprc = [f"leaprc.{fflip}"]
        for req in leaprc_reqs:
            if req not in lipid_leaprc:
                lipid_leaprc.append(req)
        lipid_ff = {"name": fflip, "leaprc": lipid_leaprc}

        water_map = {
            "opc": {"name": "OPC", "leaprc": "leaprc.water.opc", "box": "OPCBOX"},
            "tip3p": {"name": "TIP3P", "leaprc": "leaprc.water.tip3p", "box": "TIP3PBOX"},
            "spceb": {"name": "SPC/Eb", "leaprc": "leaprc.water.spceb", "box": "SPCBOX"},
            "opc3": {"name": "OPC3", "leaprc": "leaprc.water.opc3", "box": "OPC3BOX"},
            "spce": {"name": "SPC/E", "leaprc": "leaprc.water.spce", "box": "SPCBOX"},
            "tip4pew": {"name": "TIP4P-Ew", "leaprc": "leaprc.water.tip4pew", "box": "TIP4PEWBOX"},
        }
        water_ff = water_map.get(water_model, {"name": water_model, "leaprc": f"leaprc.water.{water_model}", "box": "TIP3PBOX"})

        ion_map = {
            "OPC": {"name": "12-6-4 OPC", "frcmod": "frcmod.ionslm_1264_opc"},
            "TIP3P": {"name": "12-6-4 TIP3P", "frcmod": ["frcmod.ions1lm_1264_tip3p", "frcmod.ions234lm_1264_tip3p"]},
            "SPC/E": {"name": "12-6-4 SPC/E", "frcmod": ["frcmod.ions1lm_1264_spce", "frcmod.ions234lm_1264_spce"]},
            "SPC/Eb": {"name": "12-6-4 SPC/E", "frcmod": ["frcmod.ions1lm_1264_spce", "frcmod.ions234lm_1264_spce"]},
            "TIP4P-Ew": {"name": "12-6-4 TIP4P-Ew", "frcmod": ["frcmod.ions1lm_1264_tip4pew", "frcmod.ions234lm_1264_tip4pew"]},
        }
        ion_ff = ion_map.get(water_ff["name"], {"name": "Default only", "frcmod": None})

        selected = {
            "protein": protein_ff,
            "lipids": lipid_ff,
            "water": water_ff,
            "ions": ion_ff,
        }

        self.update_workspace("selected_standard_forcefields", selected)

        # Only list leaprc additions that aren't the base lipid leaprc itself
        # (e.g. leaprc.lipid_ext, addPath for extended lipids, leaprc.extra_solvents).
        base_lipid_leaprc = f"leaprc.{fflip}"
        extra_leaprc = [r for r in leaprc_reqs if r != base_lipid_leaprc]
        if extra_leaprc:
            lipids_line = f"{fflip} (+ {', '.join(extra_leaprc)})"
        else:
            lipids_line = fflip

        console.print("[bold green]Membrane system — forcefields auto-configured from membrane builder:[/bold green]")
        console.print(f"  Protein: {protein_ff['name']}")
        console.print(f"  Lipids:  {lipids_line}")
        console.print(f"  Water:   {water_ff['name']}")
        console.print(f"  Ions:    {ion_ff['name']}")
        console.print("[grey50]Solvation skipped — system already solvated by membrane builder[/grey50]")

        return selected

    def _select_standard_forcefields_interactive(self) -> dict:
        """
        Interactive UX for selecting standard AMBER forcefields.

        Presents each category sequentially (matching the template guide structure).
        User selects one option per category, or "None" if not applicable.

        Returns:
            Dictionary with selected forcefields, stored in workspace as
            'selected_standard_forcefields'
        """
        from rich.panel import Panel
        from rich.table import Table
        from proprep.forcefield_params.forcefield_menu import render_forcefield_category
        console = self.processor.console

        console.print(Panel(
            "[bold]FORCEFIELD SELECTION[/bold]\n\n"
            "ProPrep will guide you through forcefield selection. Your choices will be\n"
            "used to generate the tLEaP input file and estimate solvation parameters.",
            title="Standard Forcefields",
            border_style="blue",
            expand=False
        ))

        selected = {}

        # Process each category in order
        # AMBER load order: protein -> modified AA -> nucleic acid -> carbs -> gaff -> water -> ions
        category_order = ['protein', 'modified_aa', 'dna', 'rna', 'carbohydrates',
                         'lipids', 'small_molecules', 'water', 'ions']

        # Check workspace for constant pH data to adjust protein FF recommendation
        has_constant_ph = bool(
            self.get_from_workspace("constant_ph_residues", None)
            or self.get_from_workspace("constant_ph_data", None)
        )

        # Layer 2 — gather the per-cofactor AND-group structure for the
        # selected redox-site transformers. A picker option is tagged
        # "satisfies your <residues> selections" iff it FULLY SATISFIES the
        # cofactor's prereqs (its leaprc set hits ≥1 entry in every AND-group
        # of at least one cofactor). The fully-satisfying rule fixes the
        # multi-protein-FF-REQUIRED bug: for an OR-group like Zn(Cys)4's "any
        # parm10-compatible protein FF", every compatible row is tagged
        # individually so the user sees "any of these works for Zn"; for an
        # AND-of-ANDs like bis-his hemes (needs both constph AND conste in
        # one pick), only the combined "Constant pH + Redox" row gets tagged
        # because options 8/9 alone leave one group unsatisfied.
        cofactor_prereq_groups = self._collect_cofactor_prereq_groups()
        cofactor_prereq_requesters = self._collect_cofactor_prereq_requesters()
        cofactor_prereq_leaprcs = set(cofactor_prereq_requesters.keys())
        cofactor_any_loaded = bool(
            self.get_from_workspace("transformer_info", []) or []
        ) and any(
            si.get("has_transformer") and si.get("cofactor_path")
            and si.get("transformer_type") != "no_transformation"
            for si in (self.get_from_workspace("transformer_info", []) or [])
        )

        for category in category_order:
            if category not in self.FORCEFIELD_OPTIONS:
                continue

            cat_info = self.FORCEFIELD_OPTIONS[category]
            options = cat_info['options']

            # If constant pH data found, recommend constph instead of ff19SB
            if category == 'protein' and has_constant_ph:
                console.print("\n[yellow]Note: Constant pH residues detected in workspace. "
                              "Recommending Constant pH (ff10) forcefield.[/yellow]")
                options = [dict(opt) for opt in options]  # shallow copy to avoid mutating class data
                for opt in options:
                    if opt.get('leaprc') == 'leaprc.constph':
                        opt['recommended'] = True
                    else:
                        opt.pop('recommended', None)

            # Layer 2 — mark a picker option REQUIRED iff it FULLY SATISFIES
            # at least one cofactor (its leaprc set covers every AND-group of
            # that cofactor). The trigger label is the set of cofactor residue
            # names this option fully satisfies. Some options carry a
            # list-valued `leaprc` (e.g. the combined "Constant pH + Redox"
            # protein FF); the set-intersection logic in
            # `_option_fully_satisfies` handles that uniformly.
            def _opt_full_triggers(opt):
                """Set of residue names this option fully satisfies."""
                opt_leaprcs = self._option_leaprc_set(opt)
                if not opt_leaprcs:
                    return set()
                return {
                    cof["residue_name"]
                    for cof in cofactor_prereq_groups
                    if self._option_fully_satisfies(opt_leaprcs, cof["groups"])
                }
            if cofactor_prereq_groups and any(
                _opt_full_triggers(opt) for opt in options
            ):
                options = [dict(opt) for opt in options]
                for opt in options:
                    triggers = _opt_full_triggers(opt)
                    if triggers:
                        res_join = " + ".join(sorted(triggers))
                        # "satisfies" not "REQUIRED": for an OR-group (e.g.
                        # Zn(Cys)4 accepts any parm10 protein FF) several rows
                        # each satisfy the cofactor, so no single one is
                        # required — saying "REQUIRED" misleads (e.g. constph
                        # reads as "must use constant-pH" when it's just one
                        # acceptable parm10 source).
                        reason = f"satisfies your {res_join} selections"
                        opt['recommended'] = True
                        opt['recommendation_reason'] = reason
                        opt['cofactor_satisfier'] = True
                        # Stash the trigger count so the default-selection
                        # loop below can prefer the option that satisfies the
                        # most cofactors (avoids defaulting to a row that
                        # leaves another cofactor's prereqs unmet).
                        opt['_cofactor_trigger_count'] = len(triggers)
                    else:
                        # Don't clobber existing 'recommended' marks (e.g. constph above)
                        if not opt.get('recommendation_reason'):
                            opt.pop('recommended', None)

            # Layer 2b — when cofactors are loaded, the protein FF is
            # effectively required (its parm10 base covers the cofactors'
            # standard atom types). Mark "None" with a visible incompatibility
            # warning so the user is reminded as they scan the menu.
            if category == 'protein' and cofactor_any_loaded:
                cat_info = dict(cat_info)
                cat_info['none_text'] = (
                    "← Not compatible with your cofactor selections "
                    "(cofactors need parm10 types from a protein FF)"
                )

            # Recommend the water model that matches the protein FF the user
            # just picked. Water is processed after protein in category_order,
            # so selected['protein'] is already populated here. Pairings live in
            # the catalog (recommended_water_for_protein): ff19SB→OPC,
            # ff14SB→TIP3P, constant-pH / constant-redox / both→TIP3P, etc. When
            # no protein FF was chosen or its pairing is unknown, we leave the
            # catalog's own default (OPC) marked.
            if category == 'water':
                protein_sel = selected.get('protein')
                rec_water = (
                    recommended_water_for_protein(protein_sel.get('leaprc'))
                    if protein_sel else None
                )
                if rec_water:
                    options = [dict(opt) for opt in options]  # don't mutate class data
                    for opt in options:
                        if opt['name'] == rec_water:
                            opt['recommended'] = True
                            opt['recommendation_reason'] = (
                                f"matches {protein_sel['name']}"
                            )
                        else:
                            opt.pop('recommended', None)

            # Recommend the divalent+ ion set that matches the water model the
            # user just picked. Ions are processed after water in
            # category_order, so selected['water'] is populated here. Without
            # this, the catalog's static default (12-6-4 OPC) stays marked even
            # when the user chose, e.g., TIP3P — recommending an OPC-calibrated
            # ion set against TIP3P water. The matching set is resolved from
            # each option's `for_water` tag (recommended_ions_for_water),
            # preferring the most accurate 12-6-4 variant and falling back to
            # "Default only" for water models with no dedicated Li/Merz set.
            if category == 'ions':
                water_sel = selected.get('water')
                water_name = water_sel.get('name') if water_sel else None
                rec_ion = recommended_ions_for_water(water_name)
                options = [dict(opt) for opt in options]  # don't mutate class data
                for opt in options:
                    if opt['name'] == rec_ion:
                        opt['recommended'] = True
                        if rec_ion == 'Default only' and water_name:
                            opt['recommendation_reason'] = (
                                f"no Li/Merz 12-6-4 set for {water_name}"
                            )
                        elif water_name:
                            opt['recommendation_reason'] = f"matches {water_name}"
                        else:
                            opt['recommendation_reason'] = "no explicit water model"
                    else:
                        opt.pop('recommended', None)

            # "None" is always row 1 when allowed, so the user doesn't have to
            # count to a different final index per category.
            allow_none = cat_info.get('allow_none', False)

            # Shared renderer (also used by the component parameterizers) so the
            # curated-catalog picker looks identical everywhere it appears.
            none_offset = render_forcefield_category(
                console,
                cat_info['title'],
                cat_info.get('description'),
                options,
                allow_none=allow_none,
                none_text=cat_info.get('none_text', 'Skip this category'),
            )

            # Handle multi-select categories (like 'special')
            if cat_info.get('multi_select', False):
                if allow_none:
                    console.print("[grey50]Enter 1 for none, or comma-separated choices (e.g., 2,3), or press Enter for none[/grey50]")
                else:
                    console.print("[grey50]Enter choices separated by commas (e.g., 1,2)[/grey50]")
                choice_str = prompt_with_context(
                    self.processor, f"Select {cat_info['title'].lower()}", default="",
                    module="Topology Generator",
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
                # Single-select category.
                # Default to the recommended option with the HIGHEST cofactor
                # trigger count (fully satisfies the most cofactors); falls
                # back to the first recommended option, then to the first
                # real option if nothing is marked recommended. The trigger
                # count lives on `_cofactor_trigger_count` and is set by
                # Layer 2 above; options not touched by Layer 2 default to 0.
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
                    module="Topology Generator",
                    description=f"Select {cat_info['title'].lower()} forcefield",
                )

                try:
                    choice = int(choice_str)
                    if allow_none and choice == 1:
                        # User selected "None" (row 1 when allowed)
                        selected[category] = None
                    elif 1 + none_offset <= choice <= len(options) + none_offset:
                        selected[category] = options[choice - 1 - none_offset]
                    else:
                        # Invalid choice, use default (first real option)
                        selected[category] = options[0]
                except ValueError:
                    # Invalid input, use default
                    selected[category] = options[0]

        # Layer 3 — verify every cofactor-prereq AND-group is satisfied by at
        # least one sourced leaprc. With OR-groups, a flat leaprc-by-leaprc
        # check is wrong (e.g. Zn's "any of 6 protein FFs" group is satisfied
        # by ANY one pick, not all). The group-aware resolver prompts only
        # for groups that remain unsatisfied after the picker.
        if cofactor_prereq_groups:
            selected = self._resolve_cofactor_prereq_mismatches(
                selected, cofactor_prereq_groups, category_order
            )

        # Track which leaprcs the standard-FF block sourced — used by the
        # cofactor params section to know what's already loaded.
        sourced_leaprcs = self._extract_sourced_leaprcs(selected)
        self.update_workspace("standard_ff_leaprcs_sourced", sourced_leaprcs)

        # Store in workspace
        self.update_workspace("selected_standard_forcefields", selected)

        # Display summary
        console.print("\n[bold green]Selected Forcefields:[/bold green]")
        for category in category_order:
            if category not in selected:
                continue
            sel = selected[category]
            if sel is None:
                console.print(f"  {category}: [grey50]None[/grey50]")
            elif isinstance(sel, list):
                if sel:
                    names = [s['name'] for s in sel]
                    console.print(f"  {category}: {', '.join(names)}")
                else:
                    console.print(f"  {category}: [grey50]None[/grey50]")
            else:
                console.print(f"  {category}: {sel['name']}")

        return selected

    # ========================================================================
    # Cofactor prerequisite-leaprc plumbing (Layers 1/2/3 + emission helper)
    # ========================================================================

    def _collect_cofactor_prereq_leaprcs(self) -> set:
        """Gather the union of `prerequisites.leaprcs` from every selected
        redox-site transformer's metadata.

        Reads `transformer_info` from workspace (populated by
        :py:meth:`redox_transformation_manager._store_transformer_info_in_workspace`).
        Each entry's `cofactor_path` (also threaded in by the same method)
        identifies which metadata.json to read.

        Returns:
            Set of leaprc names like {"leaprc.gaff2", "leaprc.RNA.OL3"}, or
            an empty set if no transformers declare prereqs.
        """
        return set(self._collect_cofactor_prereq_requesters().keys())

    def _collect_cofactor_prereq_requesters(self) -> dict:
        """Like :py:meth:`_collect_cofactor_prereq_leaprcs` but keyed by leaprc
        with the value being the deduplicated set of residue names that
        triggered the prereq. Uses the flat union of all group satisfiers, so
        it's only meaningful as a "which leaprcs are mentioned anywhere in
        these cofactors' prereqs" lookup — not as a per-leaprc requirement
        marker (use :py:meth:`_collect_cofactor_prereq_groups` for that).

        Returns:
            ``{leaprc: {residue_name, ...}}`` mapping. Empty when no
            transformer declares a prereq leaprc.
        """
        try:
            from proprep.forcefield_params import get_prerequisite_leaprcs
        except ImportError:
            return {}

        transformer_info = self.get_from_workspace("transformer_info", []) or []
        requesters: dict = {}
        for site_info in transformer_info:
            cofactor_path = site_info.get("cofactor_path")
            if not cofactor_path:
                continue
            try:
                leaprcs = get_prerequisite_leaprcs(cofactor_path)
            except Exception as e:
                logger.warning(
                    f"Could not load prereqs for {cofactor_path}: {e}"
                )
                continue
            resname = site_info.get("residue_name") or "?"
            for leaprc in leaprcs:
                requesters.setdefault(leaprc, set()).add(resname)
        return requesters

    def _collect_cofactor_prereq_groups(self) -> list:
        """Return the per-cofactor AND-group structure used by the picker's
        fully-satisfying rule.

        Each entry corresponds to one selected redox-site transformer that
        declares prerequisites, with its AND-groups (each group is an OR list
        of satisfying leaprcs). Cofactors with no prerequisites are skipped.

        Returns:
            List of dicts, each with::

                {
                    "residue_name": "ZNC",
                    "cofactor_path": "zinc/cys4",
                    "groups": [["leaprc.protein.ff14SB", "leaprc.protein.ff19SB", ...]],
                }

            Empty list when no transformer declares any prereq groups.
        """
        try:
            from proprep.forcefield_params import get_prerequisite_leaprc_groups
        except ImportError:
            return []

        transformer_info = self.get_from_workspace("transformer_info", []) or []
        cofactors: list = []
        for site_info in transformer_info:
            cofactor_path = site_info.get("cofactor_path")
            if not cofactor_path:
                continue
            # Per-set prerequisites: a fixed-pH set requires ff14SB/ff19SB (not
            # leaprc.constph), so resolve which set this site's parameters imply
            # and read that set's prereqs. Falls back to the cofactor-global block
            # when the set has none / can't be resolved.
            set_name = self._resolve_site_forcefield_set(site_info)
            try:
                groups = get_prerequisite_leaprc_groups(cofactor_path, set_name)
            except Exception as e:
                logger.warning(
                    f"Could not load prereq groups for {cofactor_path}"
                    f"{f' (set {set_name})' if set_name else ''}: {e}"
                )
                continue
            if not groups:
                continue
            chain = site_info.get("chain")
            resid = site_info.get("residue_id")
            cofactors.append({
                "residue_name": site_info.get("residue_name") or "?",
                "location": f"{chain}.{resid}" if chain and resid is not None else "",
                "cofactor_path": cofactor_path,
                "groups": groups,
            })
        return cofactors

    def _resolve_site_forcefield_set(self, site_info: dict):
        """Best-effort name of the forcefield set a site's transformer parameters
        imply, so per-set leaprc prerequisites apply (e.g. fixed-pH → ff14SB,
        not constph). Prefers a value already recorded in transformer_info;
        otherwise asks the transformer. Returns None on any failure → caller
        falls back to the cofactor-global prerequisites.
        """
        recorded = site_info.get("forcefield_set")
        if recorded:
            return recorded
        try:
            from proprep.redoxsite_prep.transformation.redox_transformer_framework import (
                redox_transformer_registry,
            )
            tclass = redox_transformer_registry.get_transformer(site_info.get("transformer_type"))
            if tclass and hasattr(tclass, "select_forcefield_set_name"):
                return tclass.select_forcefield_set_name(site_info.get("parameters", {}) or {})
        except Exception as e:
            logger.debug("Could not resolve forcefield set for prereqs: %s", e)
        return None

    def _preferred_ff_set_for_combo(self, transformer_type, redox_state, spin_state):
        """Forcefield set implied by a site's Stage-1 pH-treatment choice for this
        exact (transformer, redox, spin) combo, read from workspace
        transformer_info. Used to pre-select the matching set in the FF-set
        picker. Returns None when no matching site recorded a pH-treatment choice
        (→ no pre-selection, existing default behavior).
        """
        for si in (self.get_from_workspace("transformer_info", []) or []):
            if (si.get("transformer_type") == transformer_type
                    and si.get("redox_state") == redox_state
                    and si.get("spin_state") == spin_state
                    and si.get("ph_treatment")):
                return si.get("forcefield_set")
        return None

    def _preferred_ph_treatment_for_combo(self, transformer_type, redox_state, spin_state):
        """The fixed-pH/constant-pH choice a site recorded in the redox site
        preparer (Stage 1) for this exact (transformer, redox, spin) combo, read
        from workspace transformer_info. Used to filter the multi-microstate
        FF-set picker to the chosen treatment. Returns None when no matching site
        recorded a choice (→ no filtering, existing behavior). Mirrors
        _preferred_ff_set_for_combo, which returns the implied set name.
        """
        for si in (self.get_from_workspace("transformer_info", []) or []):
            if (si.get("transformer_type") == transformer_type
                    and si.get("redox_state") == redox_state
                    and si.get("spin_state") == spin_state
                    and si.get("ph_treatment")):
                return si.get("ph_treatment")
        return None

    @staticmethod
    def _ph_treatment_phrase(ph_treatment):
        """User-facing label for a set's propionate pH treatment. The FF-set
        metadata describes only the chemistry (charge model + provenance); the
        picker is the single source of truth for the fixed-pH/constant-pH label,
        so it stays consistent across every cofactor and can't drift. Returns
        None for legacy/treatment-agnostic sets (ph_treatment=None)."""
        return {
            'fixed_pH': 'fixed-pH (static PRP/PRD)',
            'constant_pH': 'constant-pH (titratable PRN)',
        }.get(ph_treatment)

    @classmethod
    def _ff_set_display_title(cls, option: dict) -> str:
        """Display title for a forcefield set with a pH-treatment tag appended.

        A cofactor's fixed-pH and constant-pH sets often share an identical
        metadata ``name`` (the shared parameter provenance, e.g. "... (conste)"),
        so the bare name reads as ambiguous — a user can't tell which treatment
        they're being offered. The pH treatment is the axis that distinguishes
        them, so surface it in the title. No-op for treatment-agnostic sets.
        """
        title = option.get('display_name') or option.get('name') or '?'
        tag = {'fixed_pH': ' [fixed-pH]',
               'constant_pH': ' [constant-pH]'}.get(option.get('ph_treatment'), '')
        return f"{title}{tag}"

    @staticmethod
    def _option_leaprc_set(opt: dict) -> set:
        """Normalize a picker option's `leaprc` field to a set, handling
        scalar string, list, or missing values."""
        lr = opt.get('leaprc')
        if isinstance(lr, list):
            return set(lr)
        if lr:
            return {lr}
        return set()

    @staticmethod
    def _option_fully_satisfies(opt_leaprcs: set, groups: list) -> bool:
        """True iff the option's leaprc set hits at least one entry in EVERY
        AND-group (i.e. it fully satisfies the cofactor on its own)."""
        if not groups:
            return False
        return all(any(lr in opt_leaprcs for lr in g) for g in groups)

    # Compact, user-facing labels for the cofactors carried by the v1.1
    # fragment-typed library. Used by the cofactor-requirements panel below.
    def _show_cofactor_ff_prerequisites_panel(self) -> None:
        """One-shot summary shown before the standard-FF picker: which leaprcs
        the selected cofactors' parameter sets require.

        This used to be an essay. It asserted that every cofactor needs a
        protein force field (only zinc/cys4 declares one), pattern-matched
        residue NAMES to decide GAFF2 was needed (while flavin/fad declares
        exactly that in its metadata, unread), and explained the requirement in
        terms of a ribitol tail and which bond would fail -- none of which is
        known about an arbitrary parameter set, and all of which was wrong for
        an imported one.

        What is knowable is what the sets declare, so that is what is shown.
        Requirements are enforced after selection by
        _resolve_cofactor_prereq_mismatches against the same declarations; this
        exists only so the requirement is visible BEFORE picking.
        """
        from rich.table import Table

        console = self.processor.console
        cofactors = self._collect_cofactor_prereq_groups()
        if not cofactors:
            return

        table = Table(title="Cofactor force-field prerequisites", expand=False)
        table.add_column("Cofactor", style="cyan")
        table.add_column("Requires", style="green")

        for entry in cofactors:
            residue = entry.get("residue_name", "?")
            location = entry.get("location") or ""
            label = f"{residue} {location}".strip()
            for i, group in enumerate(entry["groups"]):
                # One row per AND-group; members are alternatives.
                requirement = " or ".join(group)
                table.add_row(label if i == 0 else "", requirement)

        console.print()
        console.print(table)
        console.print(
            "[grey50]Every line must be satisfied by the force fields selected "
            "below. Anything unmet is flagged after selection.[/grey50]")

    def _resolve_cofactor_prereq_mismatches(
        self, selected: dict, cofactor_prereq_groups: list, category_order: list
    ) -> dict:
        """Layer 3 — after the user finishes the picker, verify every
        cofactor's AND-groups are satisfied. For each unsatisfied group,
        prompt for an explicit resolution (add a satisfier / reconsider the
        category pick / force-skip).

        Group-aware semantics: a group is satisfied iff the sourced leaprc
        set intersects the group's `satisfied_by` (OR list). All groups must
        be satisfied (AND across groups) for a cofactor to be cleanly loaded.

        Updates `selected` in-place and also returns it.
        """
        from rich.panel import Panel

        console = self.processor.console
        sourced = self._extract_sourced_leaprcs(selected)

        # Build the de-duplicated list of unsatisfied groups, each tagged
        # with the cofactors that need it. A group satisfied by ANY source
        # leaprc is dropped. Groups with the same satisfied_by signature from
        # different cofactors are merged so the user is asked once.
        unsatisfied: list = []  # [{"satisfied_by": frozenset, "requesters": [str, ...]}]
        sig_to_index: dict = {}
        for cof in cofactor_prereq_groups:
            resname = cof.get("residue_name", "?")
            for group in cof.get("groups", []):
                if any(lr in sourced for lr in group):
                    continue  # already satisfied
                sig = frozenset(group)
                requester = f"{cof.get('cofactor_path', '?')} ({resname})"
                if sig in sig_to_index:
                    unsatisfied[sig_to_index[sig]]["requesters"].append(requester)
                else:
                    sig_to_index[sig] = len(unsatisfied)
                    unsatisfied.append({
                        "satisfied_by": list(group),
                        "requesters": [requester],
                    })

        if not unsatisfied:
            return selected

        # For each unsatisfied group, find the picker option that's the
        # natural "add this to fix it" injection target. Prefer the first
        # satisfier (in OR-list order) that appears as a real picker option.
        # Track which category it lives in so we can offer "reconsider this
        # category" as option 2.
        def _find_injection_target(satisfied_by):
            """Return (category, option_dict, satisfier_leaprc) or (None, None, None)."""
            for satisfier in satisfied_by:
                for category in category_order:
                    cat_info = self.FORCEFIELD_OPTIONS.get(category, {})
                    for opt in cat_info.get('options', []) or []:
                        lr = opt.get('leaprc')
                        lrs = lr if isinstance(lr, list) else ([lr] if lr else [])
                        if satisfier in lrs:
                            return category, opt, satisfier
            return None, None, None

        confirmed_extras: list = []
        forced_skips: list = []
        for group_info in unsatisfied:
            satisfied_by = group_info["satisfied_by"]
            requesters = sorted(set(group_info["requesters"]))
            category, injected_opt_template, primary_satisfier = (
                _find_injection_target(satisfied_by)
            )
            category_label = category or "unknown"
            current_pick = selected.get(category) if category else None
            current_pick_name = (
                "None" if current_pick is None
                else current_pick.get('name', '?') if isinstance(current_pick, dict)
                else str(current_pick)
            )

            # Build the OR-list display for the panel
            satisfiers_str = ", ".join(f"[bold]{lr}[/bold]" for lr in satisfied_by)
            is_or_group = len(satisfied_by) > 1
            primary_label = primary_satisfier or satisfied_by[0]

            if is_or_group:
                requirement_line = (
                    f"Your standard-FF picks did not include any of: {satisfiers_str}."
                )
            else:
                requirement_line = (
                    f"Your standard-FF picks did not include [bold]{primary_label}[/bold]."
                )
            panel_body = [
                requirement_line,
                "",
                f"  Category: {category_label}",
                f"  Your pick: [yellow]{current_pick_name}[/yellow]",
                "",
                "But the following selected redox-site cofactor(s) require it:",
            ]
            for req in requesters:
                panel_body.append(f"  • {req}")
            panel_body.extend([
                "",
                "Without it, tleap will fail at saveAmberParm with errors",
                "like \"atom type X for atom Y in residue Z was not found\".",
                "",
                "Options:",
                f"  [bold]1[/bold]. Add [bold]{primary_label}[/bold] to my forcefield set    [Recommended]",
                f"  [bold]2[/bold]. Go back and reconsider my pick for this category",
                "  [bold]3[/bold]. Proceed without it (tleap will fail; useful only for debugging",
                "      or if you plan to edit the tleap.in manually)",
            ])
            panel_title = (
                "[bold yellow]Forcefield prerequisite mismatch — "
                f"{primary_label}[/bold yellow]"
            )

            console.print()
            console.print(Panel(
                "\n".join(panel_body),
                title=panel_title,
                border_style="yellow",
                expand=False,
            ))

            choice = prompt_with_context(
                self.processor,
                f"Resolution for {primary_label}",
                default="1",
                module="Topology Generator",
                description=f"Resolve missing cofactor prereq {primary_label}",
            ).strip()

            if choice == "2":
                # Re-prompt the category. Recursive: re-run the FULL picker.
                console.print(
                    f"[blue]Returning to forcefield picker — reconsider the "
                    f"{category_label} category.[/blue]"
                )
                return self._select_standard_forcefields_interactive()
            elif choice == "3":
                forced_skips.append(primary_label)
                console.print(
                    f"[red]WARNING: proceeding without {primary_label}. "
                    f"tleap will fail. This pick is preserved in the generated tleap.in "
                    f"as a comment.[/red]"
                )
            else:
                # Option 1 (default): inject the chosen satisfier's
                # FORCEFIELD_OPTIONS entry into selected[category] so the
                # leaprc lands in the standard-FF block. For non-multi-select
                # categories, this overrides whatever was previously picked
                # (typically "None" or a row that didn't satisfy this group).
                if injected_opt_template is not None and category is not None:
                    injected_opt = dict(injected_opt_template)
                    injected_opt['recommendation_reason'] = (
                        "confirmed at FF-picker as a cofactor prerequisite"
                    )
                    cat_info = self.FORCEFIELD_OPTIONS.get(category, {})
                    if cat_info.get('multi_select'):
                        cur = selected.get(category)
                        if not isinstance(cur, list):
                            cur = []
                        cur.append(injected_opt)
                        selected[category] = cur
                    else:
                        selected[category] = injected_opt
                    # Update sourced set so subsequent groups see this pick
                    # (a later group might be auto-satisfied by it).
                    sourced.update(self._option_leaprc_set(injected_opt))
                    confirmed_extras.append(primary_label)
                else:
                    # No matching category option (shouldn't happen for the
                    # bundled metadata, but stay defensive): fall back to the
                    # sidecar list so the cofactor params section emits it.
                    confirmed_extras.append(primary_label)
                    selected.setdefault("cofactor_prereqs", []).append(primary_label)

        if forced_skips:
            selected.setdefault("cofactor_prereqs_skipped", []).extend(forced_skips)

        return selected

    @staticmethod
    def _option_has_methodology_content(option: dict) -> bool:
        """True if there's anything for `_show_ff_methodology_panel` to render
        — used to gate the picker's methodology peek prompt."""
        return bool(
            option.get('methodology')
            or option.get('cofactor_methodology')
            or option.get('cofactor_set_comparison_guidance')
        )

    def _show_ff_methodology_panel(self, display_name: str, option: dict) -> None:
        """Render a methodology info panel for one force-field set.

        Composes three sections (in order):
          1. Cofactor-level "Parameterization protocol" — shared across all
             bundled variants of this cofactor. Populated when the cofactor
             metadata declares a top-level `methodology` field (typical for
             cofactors whose variants share most of the parameterization
             pipeline and differ only in detail, like Fe4S4 BS guesses).
          2. Per-set "This variant" — what's specific to this particular set
             (e.g., the spin-coupling pattern, or the difference between AMBER
             library bonded vs QM Hessian bonded variants for the cofactor
             library).
          3. Cofactor-level "Choosing between bundled variants" — the chooser
             guidance text users need to pick informedly.

        Sections are omitted when the corresponding metadata field is empty.
        """
        from rich.panel import Panel
        from rich.console import Group
        from rich.text import Text
        from rich.table import Table
        from rich.padding import Padding
        import re

        console = self.processor.console
        cofactor_methodology = option.get('cofactor_methodology', '') or ''
        per_set_methodology = option.get('methodology', '') or ''
        chooser_guidance = option.get('cofactor_set_comparison_guidance', '') or ''

        def _fmt_block(text: str) -> list:
            """Render a prose blob as scannable renderables: paragraphs wrap as
            Text; a run of dash bullets becomes a hanging-indented list whose
            continuation lines align under the bullet text. Each bullet's short
            lead-in (its topic, up to the first ':' or '(') is bold-blue, as are
            the 'Note:'/'References:' markers. Bold blue (not bold-default) keeps
            the emphasis legible on a white background."""
            out = []
            for i, para in enumerate(re.split(r'\n\s*\n', text.strip())):
                if i:
                    out.append(Text(""))  # blank line between paragraphs
                lines = [l for l in para.split('\n') if l.strip()]
                if lines and all(l.lstrip().startswith('- ') for l in lines):
                    grid = Table.grid(padding=(0, 1, 0, 0), expand=True)
                    grid.add_column(no_wrap=True, vertical="top")  # bullet glyph
                    grid.add_column(overflow="fold", ratio=1)       # wrapped text
                    for j, l in enumerate(lines):
                        if j:
                            grid.add_row("", "")  # blank line between bullets
                        rest = l.lstrip()[2:]
                        m = re.match(r'(.{1,48}?)\s*([:(])', rest)
                        if m:
                            lead = m.group(1).strip()
                            rest = f"[bold blue]{lead}[/bold blue]" + rest[len(m.group(1)):]
                        grid.add_row("•", Text.from_markup(rest))
                    out.append(Padding(grid, (0, 0, 0, 2)))
                else:
                    t = para
                    if t.startswith('Note:'):
                        t = '[bold blue]Note:[/bold blue]' + t[5:]
                    elif t.startswith('References:'):
                        t = '[bold blue]References:[/bold blue]' + t[11:]
                    out.append(Text.from_markup(t))
            return out

        sections = [
            ("Parameterization protocol", cofactor_methodology),
            ("This variant", per_set_methodology),
            ("Choosing between bundled variants", chooser_guidance),
        ]
        parts = []
        for header, content in sections:
            if not content:
                continue
            if parts:
                parts.append(Text(""))  # spacer between sections
            parts.append(Text.from_markup(f"[bold blue]{header}[/bold blue]"))
            parts.append(Text(""))
            parts.extend(_fmt_block(content))
        if not parts:
            parts.append(Text.from_markup(
                "[grey50]No methodology information available for this set.[/grey50]"))

        console.print()
        console.print(Panel(
            Group(*parts),
            title=f"[bold blue]{display_name} -- parameterization details[/bold blue]",
            border_style="blue",
            width=min(console.width, 92),
            expand=False,
        ))
        console.print()

    def _extract_sourced_leaprcs(self, selected: dict) -> set:
        """Return the set of leaprc names that the user's picks will source.

        Walks `selected` (the output of the standard-FF picker) and pulls
        the `leaprc` field from each picked option. Used by Layer 3 to detect
        missing prereqs, and stored to workspace for the cofactor-params
        section's emission logic.
        """
        sourced: set = set()
        for category, pick in selected.items():
            if pick is None:
                continue
            if isinstance(pick, list):
                for item in pick:
                    if isinstance(item, dict) and item.get('leaprc'):
                        lr = item['leaprc']
                        if isinstance(lr, list):
                            sourced.update(lr)
                        else:
                            sourced.add(lr)
            elif isinstance(pick, dict) and pick.get('leaprc'):
                lr = pick['leaprc']
                if isinstance(lr, list):
                    sourced.update(lr)
                else:
                    sourced.add(lr)
            elif isinstance(pick, str):
                # Some entries store the raw leaprc string
                if pick.startswith("leaprc."):
                    sourced.add(pick)
        # Plus any extras the user explicitly confirmed adding via Layer 3
        for extra in selected.get("cofactor_prereqs", []) or []:
            sourced.add(extra)
        return sourced

    # Fixed-pH propionate residue names: present in a structure only when the
    # heme propionates have been switched off constant-pH (PRN) to static forms
    # — e.g. by the PB-Titrate modern-FF rebuild (pb_titrate/pdb_rename.py).
    _FIXED_PH_PROPIONATE_RESNAMES = {"PRP", "PRD"}

    def _build_uses_fixed_ph_propionate(self) -> bool:
        """True if the structure about to be built contains PRP/PRD propionates.

        These names appear only after the titratable PRN has been baked into a
        static protonation (the PB-Titrate handoff renames PRN → PRP/PRD). Their
        presence is the ground-truth signal that any heme site still recorded as
        ``constant_pH`` must be rebuilt with its fixed-pH companion set instead.
        Scans by residue name, so it is independent of residue renumbering (the
        PB topology is often tleap-renumbered relative to the ids recorded in
        transformer_info).

        Single-state only: PB-Titrate emits one renamed PDB (workspace key
        ``pb_rename_pdb_file``), consumed by single-state generation. The
        microstate path is not (yet) wired to receive PB-Titrate output — when
        a future per-microstate rebuild lands, it must run the same
        reconciliation over each microstate's structure.
        """
        from pathlib import Path

        seen = set()
        for p in (self.get_from_workspace("pb_rename_pdb_file", None),
                  self._select_priority_pdb_file(silent=True)):
            if not p:
                continue
            try:
                rp = str(Path(p).resolve())
            except Exception:
                continue
            if rp in seen or not Path(p).exists():
                seen.add(rp)
                continue
            seen.add(rp)
            try:
                structure = pmd.load_file(str(p))
            except Exception as e:
                logger.debug(f"Could not scan {p} for fixed-pH propionates: {e}")
                continue
            if any(r.name and r.name.strip() in self._FIXED_PH_PROPIONATE_RESNAMES
                   for r in structure.residues):
                return True
        return False

    def _reconcile_ph_treatment_with_structure(self, transformer_info):
        """Switch stale constant-pH heme sites to fixed-pH when the structure's
        propionates are actually static (PRP/PRD).

        The redox-site preparer records the Stage-1 pH-treatment choice on each
        site. When PB-Titrate later bakes recommended protonation states into a
        modern-FF PDB, it renames the titratable PRN propionates to static
        PRP/PRD — but the recorded ``ph_treatment``/``forcefield_set`` still say
        ``constant_pH`` and the constant-pH set. Left unreconciled, the Topology
        Generator would load the constant-pH lib (which defines PRN, not
        PRP/PRD) under leaprc.constph/conste, so tleap fails on the missing
        templates and uses the wrong (ff10) backbone.

        This mutates ``transformer_info`` in place (and re-persists the
        workspace copy so the prerequisite-emission path sees it too): for each
        heme site recorded ``constant_pH`` whose cofactor has a fixed-pH
        companion, flips ``ph_treatment`` → ``fixed_pH`` and ``forcefield_set``
        → the companion (same parameterization/charge model). The companion lib
        shares identical center/ligand units and merely adds PRP/PRD, so the
        already-transformed center and ligands are unaffected. No-op unless the
        structure actually contains PRP/PRD.
        """
        from proprep.forcefield_params import find_companion_set
        console = self.processor.console

        if not transformer_info:
            return transformer_info
        if not self._build_uses_fixed_ph_propionate():
            return transformer_info  # constant-pH path: leave everything alone

        changed = False
        for site in transformer_info:
            if not site.get("has_transformer"):
                continue
            if site.get("ph_treatment") != "constant_pH":
                continue
            cofactor_path = site.get("cofactor_path")
            redox = site.get("redox_state")
            spin = site.get("spin_state")
            if not (cofactor_path and redox and spin):
                continue
            companion = find_companion_set(
                cofactor_path, redox, spin,
                site.get("forcefield_set"), target_treatment="fixed_pH")
            if not companion:
                console.print(
                    f"[yellow]⚠ Site {site.get('site_id')} "
                    f"({site.get('residue_name')}) has static PRP/PRD "
                    f"propionates in the structure, but cofactor "
                    f"'{cofactor_path}' defines no fixed-pH forcefield set for "
                    f"{redox}/{spin}. tleap will fail on PRP/PRD — supply a "
                    f"fixed-pH parameter set.[/yellow]")
                continue
            old_set = site.get("forcefield_set")
            site["ph_treatment"] = "fixed_pH"
            site["forcefield_set"] = companion
            changed = True
            console.print(
                f"[blue]ⓘ Site {site.get('site_id')} "
                f"({site.get('residue_name')}): structure has static PRP/PRD "
                f"propionates — switching from constant-pH "
                f"({old_set or 'default'}) to fixed-pH set '{companion}' "
                f"(loads ff14SB/ff19SB, not leaprc.constph/conste).[/blue]")

        if changed:
            try:
                self.update_workspace("transformer_info", transformer_info)
            except Exception as e:
                logger.warning(f"Could not persist reconciled transformer_info: {e}")
        return transformer_info

    @staticmethod
    def _format_ff_file_basenames(value) -> str:
        """Human-readable basename(s) for a FF set's ``lib`` or ``frcmod`` field.

        ``discover_forcefield_files`` returns both as EITHER a single path string
        or a LIST (a metal site has one .lib per renamed residue, and a bonded
        frcmod plus each ligand's GAFF frcmod). ``Path(list)`` raises, so
        normalize both shapes to a comma-joined basename list.
        """
        if not value:
            return "(none)"
        items = value if isinstance(value, (list, tuple)) else [value]
        return ", ".join(Path(p).name for p in items)

    def _ensure_user_transformers_registered(self):
        """Register user transformers from ~/.proprep/transformers into the global
        registry.

        The transformation manager loads these in its ``__init__``, but topology
        generation can be reached without ever constructing that manager (e.g. a
        fresh ProPrep run that loads an already-transformed structure, or any
        restart between the redox-site preparer and here). Without this, an
        auto-emitted reuse transformer isn't in the registry, so
        ``get_transformer(...)`` returns None and the FF picker silently skips
        the site's deposited parameters. Idempotent and best-effort.
        """
        try:
            from proprep.redoxsite_prep.transformation.auto_rename import (
                load_user_transformers,
            )
            load_user_transformers()
        except Exception as e:
            logger.debug("load_user_transformers (topology) silenced: %s", e)

    def _get_single_state_forcefield_requirements(self):
        """Get forcefield requirements for the current transformed state"""
        console = self.processor.console
        self._ensure_user_transformers_registered()

        # Get stored transformer information from workspace
        transformer_info = self.get_from_workspace("transformer_info", [])
        requirements = {}

        if not transformer_info:
            # First fall-through: the redox-site preparer embeds the transformer
            # assignments in its exported JSON (final*…redox_sites_updated.json).
            # The workspace copy is session-only and lost on a ProPrep restart,
            # so recover it from that sidecar before resorting to residue-name
            # inference. This preserves cofactor assignments (NADPH/FAD/FMN/Ca²⁺)
            # whose residue names don't uniquely encode the transformer.
            console.print("[grey50]No transformer info in workspace — checking the "
                          "redox-sites JSON on disk...[/grey50]")
            transformer_info = self._load_transformer_info_from_sidecar()

        if not transformer_info:
            # Second fall-through: auto-detect from residue names. Handles the
            # "loaded a post-transformed structure, skipped transformation, no
            # sidecar" case; residue names that encode redox/spin state (HCO ≠
            # HCR) map unambiguously.
            console.print("[grey50]No transformer JSON on disk — attempting "
                          "auto-detect from residue names...[/grey50]")
            requirements = self._auto_detect_redox_requirements()
            if requirements:
                return requirements
            console.print("[yellow]No transformer information found in workspace or on disk.[/yellow]")
            console.print("[grey50]This is normal if no transformations were performed.[/grey50]")
            return requirements

        console.print(f"[green]Found transformer information for {len(transformer_info)} sites[/green]")

        # Reconcile the recorded pH-treatment against the actual structure: if a
        # heme's propionates are now static PRP/PRD (e.g. after the PB-Titrate
        # rebuild) but the site is still recorded constant_pH, switch it to its
        # fixed-pH companion set so the right lib + ff14SB/ff19SB prereqs load.
        transformer_info = self._reconcile_ph_treatment_with_structure(transformer_info)

        # Process sites that have transformers
        for site_info in transformer_info:
            if site_info.get('has_transformer', False):
                transformer_type = site_info.get('transformer_type')
                redox_state = site_info.get('redox_state', 'unknown')
                spin_state = site_info.get('spin_state', 'unknown')
                atom_types = site_info.get('atom_types', [])
                residue_name = site_info.get('residue_name', 'UNK')

                # Display N/A for no_transformation since it has no redox/spin state parameters
                if transformer_type == 'no_transformation':
                    display_redox = 'N/A'
                    display_spin = 'N/A'
                else:
                    display_redox = redox_state
                    display_spin = spin_state

                console.print(f"[blue]Site {site_info.get('site_id')}: {transformer_type}, {display_redox}, {display_spin}[/blue]")
                
                key = (transformer_type, redox_state, spin_state)
                if key not in requirements:
                    requirements[key] = {
                        'transformer_type': transformer_type,
                        'redox_state': redox_state,
                        'spin_state': spin_state,
                        'residue_name': residue_name,
                        'atom_types': atom_types,
                        # Pre-selection hint: the set the site's Stage-1 parameters
                        # imply (None when no pH-treatment fork). First site wins
                        # for a given (transformer, redox, spin) key.
                        'forcefield_set': site_info.get('forcefield_set'),
                        'ph_treatment': site_info.get('ph_treatment'),
                    }

        console.print(f"[green]Built requirements for {len(requirements)} unique transformer combinations[/green]")
        return requirements

    def _load_transformer_info_from_sidecar(self):
        """Recover the per-site transformer_info list from the redox-site
        preparer's JSON on disk, for when the workspace copy was lost (e.g. a
        ProPrep restart between the transformation step and topology generation).

        The preparer writes ``final*…redox_sites_updated.json`` (and, since the
        round-trip fix, embeds a top-level ``transformer_info`` key). This finds
        that sidecar in the working directory — preferring the one whose stem
        matches the structure topology is about to build, else the most recently
        modified — reads the embedded list, rehydrates the workspace key so the
        rest of the session sees it too, and returns the list (or None).
        """
        import glob as _glob
        import json as _json
        from pathlib import Path

        console = self.processor.console

        # Working directory: alongside the structure we're about to build.
        pdb = self._select_priority_pdb_file(silent=True)
        if not pdb:
            return None
        work_dir = Path(pdb).parent
        loaded_stem = Path(pdb).stem

        # The preparer's output is signed: starts with "final", ends with
        # "redox_sites_updated.json".
        candidates = [Path(p) for p in
                      _glob.glob(str(work_dir / "final*redox_sites_updated.json"))]
        # Be permissive if the strict signature finds nothing: any
        # *_redox_sites_updated.json carrying the embedded key is still usable.
        if not candidates:
            candidates = [Path(p) for p in
                          _glob.glob(str(work_dir / "*_redox_sites_updated.json"))]
        if not candidates:
            return None

        # Prefer a sidecar whose stem is a prefix of the loaded structure's stem
        # (the renamed/transformed PDB derives from it); else newest by mtime.
        def _rank(p: Path):
            stem = p.stem.replace("_redox_sites_updated", "")
            matches = loaded_stem.startswith(stem) or stem.startswith(loaded_stem)
            return (1 if matches else 0, p.stat().st_mtime)
        candidates.sort(key=_rank, reverse=True)

        for cand in candidates:
            try:
                with open(cand) as f:
                    data = _json.load(f)
            except Exception as e:
                console.print(f"[grey50]Could not read {cand.name}: {e}[/grey50]")
                continue
            tinfo = data.get("transformer_info")
            if tinfo:
                console.print(
                    f"[green]Recovered transformer assignments for "
                    f"{len(tinfo)} site(s) from {cand.name}.[/green]")
                # Rehydrate the workspace so later steps don't re-trigger this.
                try:
                    self.update_workspace("transformer_info", tinfo)
                except Exception as e:
                    logger.warning(f"Could not rehydrate transformer_info: {e}")
                return tinfo

        console.print(
            "[grey50]Found redox-sites JSON on disk but it carries no embedded "
            "transformer assignments (older format).[/grey50]")
        return None

    def _auto_detect_redox_requirements(self) -> dict:
        """Build a single-state requirements dict by matching residue names
        in the workspace against the residue names defined in each registered
        transformer's lib files.

        Used when the user has loaded a previously transformed structure
        (e.g. a redox-sites JSON) and skipped the transformation step itself,
        leaving `transformer_info` empty. The residue names in the structure
        already encode the redox/spin state (HCO ≠ HCR), so the mapping is
        unambiguous in the typical case.
        """
        console = self.processor.console

        # 1. Collect non-standard residue names from workspace state
        resnames = self._collect_non_standard_residue_names()
        if not resnames:
            return {}

        # 2. Build the residue → transformer index
        res_index = self._build_redox_residue_index()
        if not res_index:
            return {}

        # 3. Match resnames against the index
        matched_combos = {}     # (transformer, redox, spin) -> set of resnames
        ambiguous = []          # [(resname, [combos])]
        for rn in sorted(resnames):
            entries = res_index.get(rn)
            if not entries:
                continue
            combos = {(e['transformer_type'], e['redox_state'], e['spin_state']) for e in entries}
            if len(combos) > 1:
                ambiguous.append((rn, sorted(combos)))
                continue
            combo = combos.pop()
            matched_combos.setdefault(combo, set()).add(rn)

        if ambiguous:
            console.print("[yellow]Some residue names map to multiple transformer/redox/spin combos:[/yellow]")
            for rn, combos in ambiguous:
                console.print(f"  {rn}: {combos}")
            console.print("[yellow]Skipping ambiguous residues — please run the redox transformation step manually.[/yellow]")

        if not matched_combos:
            return {}

        # 4. Pull atom_types from each transformer's metadata.json
        atom_types_by_combo = self._lookup_atom_types_for_combos(matched_combos.keys())

        # 5. Build requirements dict in the format _select_forcefields_for_single_state expects
        requirements = {}
        for combo, resname_set in matched_combos.items():
            transformer_type, redox_state, spin_state = combo
            requirements[combo] = {
                'transformer_type': transformer_type,
                'redox_state': redox_state,
                'spin_state': spin_state,
                'residue_name': sorted(resname_set)[0],
                'atom_types': atom_types_by_combo.get(combo, []),
            }

        # 6. Show what we matched + persist a synthesized transformer_info so any
        # other consumer of that workspace key (e.g. membrane_builder) sees it too.
        console.print(f"\n[green]Auto-detected {len(requirements)} redox-site forcefield requirement(s) from residue names:[/green]")
        for combo, resname_set in matched_combos.items():
            t, r, s = combo
            console.print(f"  [blue]{t}[/blue] ({r}, {s}) — residues: {', '.join(sorted(resname_set))}")

        try:
            workspace = self.processor._get_workspace()
            transformer_info = []
            for i, (combo, info) in enumerate(requirements.items()):
                transformer_info.append({
                    'site_id': f'auto_{i}',
                    'has_transformer': True,
                    'transformer_type': info['transformer_type'],
                    'redox_state': info['redox_state'],
                    'spin_state': info['spin_state'],
                    'residue_name': info['residue_name'],
                    'atom_types': info['atom_types'],
                    'parameters': {},
                })
            workspace.update({'transformer_info': transformer_info})
        except Exception as e:
            logger.warning(f"Could not persist auto-detected transformer_info: {e}")

        return requirements

    def _collect_non_standard_residue_names(self) -> set:
        """Return the set of residue names in the workspace that aren't
        standard amino acids / nucleic acids / water / ions.

        First source: `detected_redox_sites` (which the redox-sites JSON
        populates via its loader). Each site exposes `coord_to_pdb` (dict
        on RedoxSite objects, dict in JSON) plus a `centers` list.
        Fallback: scan the priority PDB if those are empty.
        """
        from pathlib import Path

        STANDARD = {
            'ALA','ARG','ASN','ASP','CYS','CYX','CYM','GLN','GLU','GLY','HIS','HID','HIE','HIP',
            'ILE','LEU','LYS','LYN','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','HYP',
            'ACE','NME','NHE','ASH','GLH',
            'WAT','HOH','TIP3','TIP4','SOL','TP3','TP4','TP5','OPC','OPC3','SPC','SPCE',
            'Na+','K+','Cl-','LI','CS','RB','MG','CA','ZN','BR','IOD','F','FE','FE2','MN',
            'A','C','G','U','T','DA','DC','DG','DT','RA','RC','RG','RU',
        }

        resnames = set()

        sites = self.get_from_workspace('detected_redox_sites', []) or []
        for site in sites:
            # RedoxSite objects expose attributes; JSON-loaded sites are dicts.
            coord_to_pdb = (site.get('coord_to_pdb') if isinstance(site, dict)
                            else getattr(site, 'coord_to_pdb', None)) or {}
            for atom_info in coord_to_pdb.values():
                rn = atom_info.get('resname') if isinstance(atom_info, dict) else None
                if rn:
                    resnames.add(rn.strip())
            centers = (site.get('centers') if isinstance(site, dict)
                       else getattr(site, 'centers', None)) or []
            for c in centers:
                rn = c.get('resname') if isinstance(c, dict) else getattr(c, 'resname', None)
                if rn:
                    resnames.add(rn.strip())

        if not resnames:
            # Fallback: scan the priority PDB for residue names.
            pdb_path = None
            for k in ('priority_pdb_file', 'transformed_pdb_file', 'final_pdb_file', 'pdb_file'):
                pdb_path = self.get_from_workspace(k, None)
                if pdb_path:
                    break
            if pdb_path and Path(pdb_path).exists():
                try:
                    structure = pmd.load_file(pdb_path)
                    for res in structure.residues:
                        if res.name:
                            resnames.add(res.name.strip())
                except Exception as e:
                    logger.debug(f"Could not scan PDB for residues: {e}")

        return {r for r in resnames if r and r not in STANDARD}

    def _build_redox_residue_index(self) -> dict:
        """Walk every registered transformer's FORCEFIELD_PATH/<redox>/<spin>
        and parse the .lib files, building a map from residue name to the
        list of (transformer, redox, spin, lib, frcmod) entries that define it.
        """
        import re as _re
        from pathlib import Path
        from proprep.redoxsite_prep.transformation.redox_transformer_framework import redox_transformer_registry
        from proprep.forcefield_params.loader import get_forcefield_base_path

        index = {}
        try:
            base = Path(get_forcefield_base_path())
        except Exception as e:
            logger.warning(f"Could not resolve forcefield base path: {e}")
            return index

        entry_pattern = _re.compile(r'^!entry\.([A-Za-z0-9_-]+)\.unit\.atoms\s+table')

        for tname, cls in redox_transformer_registry.get_all_transformers().items():
            ff_rel = getattr(cls, 'FORCEFIELD_PATH', None)
            if not ff_rel:
                continue
            ff_root = base / ff_rel
            if not ff_root.is_dir():
                continue
            for redox_dir in sorted(p for p in ff_root.iterdir() if p.is_dir()):
                for spin_dir in sorted(p for p in redox_dir.iterdir() if p.is_dir()):
                    for lib_path in sorted(spin_dir.glob('*.lib')):
                        frcmod = lib_path.with_suffix('.frcmod')
                        try:
                            with open(lib_path) as f:
                                for line in f:
                                    m = entry_pattern.match(line)
                                    if not m:
                                        continue
                                    rn = m.group(1)
                                    index.setdefault(rn, []).append({
                                        'transformer_type': tname,
                                        'redox_state': redox_dir.name,
                                        'spin_state': spin_dir.name,
                                        'lib': str(lib_path),
                                        'frcmod': str(frcmod) if frcmod.exists() else None,
                                    })
                        except Exception as e:
                            logger.debug(f"Skipping lib {lib_path}: {e}")

        return index

    def _lookup_atom_types_for_combos(self, combos) -> dict:
        """For each (transformer_type, redox_state, spin_state) combo, return
        the atom_types list from the transformer's metadata.json, or [] on miss."""
        import json as _json
        from pathlib import Path
        from proprep.redoxsite_prep.transformation.redox_transformer_framework import redox_transformer_registry
        from proprep.forcefield_params.loader import get_forcefield_base_path

        try:
            base = Path(get_forcefield_base_path())
        except Exception:
            return {c: [] for c in combos}

        out = {}
        for combo in combos:
            tname, redox_state, spin_state = combo
            cls = redox_transformer_registry.get_transformer(tname)
            ff_rel = getattr(cls, 'FORCEFIELD_PATH', None) if cls else None
            if not ff_rel:
                out[combo] = []
                continue
            metadata_path = base / ff_rel / 'metadata.json'
            if not metadata_path.exists():
                out[combo] = []
                continue
            try:
                with open(metadata_path) as f:
                    metadata = _json.load(f)
                redox_data = metadata.get('redox_states', {}).get(redox_state, {})
                spin_data = redox_data.get('spin_states', {}).get(spin_state, {})
                out[combo] = spin_data.get('atom_types', []) or []
            except Exception as e:
                logger.debug(f"metadata read failed for {combo}: {e}")
                out[combo] = []
        return out

    def _select_forcefields_for_single_state(self, requirements):
        """Select forcefield files for single state (similar to microstate logic)"""
        console = self.processor.console
        selected_forcefields = {}
        
        if not requirements:
            return selected_forcefields
            
        console.print("\n[bold]Select Forcefield Parameters for Single State[/bold]")
        
        for key, info in requirements.items():
            transformer_type, redox_state, spin_state = key

            # Skip forcefield selection for no_transformation - it doesn't need forcefields
            if transformer_type == 'no_transformation':
                console.print(f"\n[blue]Skipping {info['residue_name']} (no_transformation - no forcefield files needed)[/blue]")
                continue

            console.print(f"\n[blue]Forcefield for {info['residue_name']} ({transformer_type}, {redox_state}, {spin_state}):[/blue]")

            # Find available forcefield files (reuse microstate logic)
            try:
                from proprep.redoxsite_prep.transformation.redox_transformer_framework import redox_transformer_registry
                from proprep.forcefield_params import discover_forcefield_files
                from pathlib import Path
                from proprep.utils.prompts import prompt_with_context, confirm_with_context

                transformer_class = redox_transformer_registry.get_transformer(transformer_type)

                if not transformer_class or not hasattr(transformer_class, 'FORCEFIELD_PATH'):
                    console.print(f"[yellow]No forcefield path defined for {transformer_type}[/yellow]")
                    continue

                # Check if FORCEFIELD_PATH is None before using it
                if transformer_class.FORCEFIELD_PATH is None:
                    console.print(f"[yellow]No forcefield files needed for {transformer_type}[/yellow]")
                    continue

                options = discover_forcefield_files(
                    transformer_class.FORCEFIELD_PATH, redox_state, spin_state
                )

                # The fixed-pH vs constant-pH propionate axis is already decided
                # in the redox site preparer (Stage 1, persisted on the site's
                # `ph_treatment`). Honour that here: show only the sets matching
                # that choice so the picker offers the remaining axis (charge
                # model) and never re-surfaces the pH variants the user ruled
                # out. Gated on an actual choice AND a non-empty match, so legacy
                # / treatment-agnostic cofactors (ph_treatment=None) and any
                # cofactor whose sets predate the fork fall through unchanged.
                chosen_ph_treatment = info.get('ph_treatment')
                if chosen_ph_treatment and options:
                    ph_matched = [
                        o for o in options
                        if o.get('ph_treatment') == chosen_ph_treatment
                    ]
                    if ph_matched:
                        options = ph_matched

                if not options:
                    console.print(f"[yellow]No forcefield files found for {transformer_type}/{redox_state}/{spin_state}[/yellow]")
                    continue

                if len(options) == 1:
                    # Only one option, confirm with user
                    option = options[0]
                    display_name = option.get('display_name') or option['name']
                    console.print(f"Found forcefield set: {display_name}")
                    console.print(f"  - frcmod: {self._format_ff_file_basenames(option['frcmod'])}")
                    console.print(f"  - lib: {self._format_ff_file_basenames(option.get('lib'))}")
                    # When the Stage-1 pH filter narrows a cofactor to one set,
                    # this branch is what the user sees — surface the propionate
                    # treatment here too (the metadata text no longer carries it).
                    ph_phrase = self._ph_treatment_phrase(option.get('ph_treatment'))
                    if ph_phrase:
                        console.print(f"  - Propionate treatment: {ph_phrase}")
                    if option.get('description'):
                        console.print(f"  - Description: {option['description']}")
                    if self._option_has_methodology_content(option):
                        # Single option — offer methodology peek without an
                        # explicit selection prompt
                        if confirm_with_context(
                            processor=self.processor,
                            prompt="Show parameterization methodology?",
                            default=False,
                            module="Topology Generator",
                            description="View FF parameterization methodology",
                        ):
                            self._show_ff_methodology_panel(display_name, option)

                    if confirm_with_context(
                        processor=self.processor,
                        prompt="Use this forcefield set?",
                        default=True,
                        module="Topology Generator",
                        description="Confirm forcefield selection"
                    ):
                        # Thread cofactor_path so prereq-leaprc emission knows
                        # which metadata to consult for this entry.
                        option = dict(option)
                        option['cofactor_path'] = transformer_class.FORCEFIELD_PATH
                        selected_forcefields[key] = option
                else:
                    # Multiple options — loop until the user makes a selection.
                    # Accept either a numeric choice (1..N) or `i<N>` to view
                    # the methodology for option N before deciding.
                    has_methodology = any(
                        self._option_has_methodology_content(opt) for opt in options
                    )
                    # Pre-selection: the set the site's Stage-1 pH-treatment choice
                    # implies (e.g. the fixed-pH set). Shown highlighted and used as
                    # the prompt default; ALL options remain selectable. Gated on an
                    # actual ph_treatment choice so cofactors WITHOUT a fork keep
                    # their existing default (option 1, no marker).
                    preferred_name = info.get('forcefield_set') if info.get('ph_treatment') else None
                    preferred_pos = next(
                        (i for i, opt in enumerate(options, 1)
                         if opt.get('name') == preferred_name), None)
                    default_choice = str(preferred_pos) if preferred_pos else "1"
                    from rich.padding import Padding
                    # pH treatment is normally pinned by the Stage-1 filter above,
                    # so the surviving sets share one treatment — name it once in
                    # the header and skip the per-line tag. If a mix survives (no
                    # Stage-1 choice, or a treatment-agnostic cofactor), fall back
                    # to a per-line chip so every set stays labelled.
                    ph_present = {o.get('ph_treatment') for o in options}
                    show_ph_chip = len([p for p in ph_present if p]) > 1
                    chosen_ph_phrase = self._ph_treatment_phrase(chosen_ph_treatment)
                    if chosen_ph_phrase and not show_ph_chip:
                        list_header = (
                            "Forcefield sets for your chosen propionate treatment "
                            f"({chosen_ph_phrase}) — differing only in charge model:"
                        )
                    else:
                        list_header = "Multiple forcefield sets available:"
                    while True:
                        console.print(list_header)
                        for i, option in enumerate(options, 1):
                            # Prefer the short, user-facing copy from the
                            # cofactor metadata's forcefield_set_summary block
                            # when present. Falls back to the longer-form
                            # methodology description otherwise.
                            summary_label = option.get('summary_label')
                            user_summary = option.get('user_summary')
                            header = summary_label or (
                                option.get('display_name') or option['name']
                            )
                            default_marker = " [default]" if option.get('is_default') else ""
                            ph_chip = (
                                {'fixed_pH': ' [fixed-pH]',
                                 'constant_pH': ' [constant-pH]'}.get(
                                    option.get('ph_treatment'), '')
                                if show_ph_chip else '')
                            preferred_marker = (
                                " [selected — matches your pH-treatment choice]"
                                if preferred_pos == i else "")
                            console.print(f"[bold]{i}. {header}[/bold]{ph_chip}{default_marker}{preferred_marker}")
                            if user_summary:
                                # `Padding(..., left=3)` re-indents EVERY wrapped
                                # continuation line to column 3, matching the
                                # first-character position. Plain `f"   {text}"`
                                # only indents line 1.
                                console.print(Padding(user_summary, (0, 0, 0, 3)))
                            elif option.get('description'):
                                console.print(Padding(
                                    f"- Description: {option['description']}",
                                    (0, 0, 0, 3),
                                ))
                            console.print()

                        prompt_text = "Select forcefield set"
                        if has_methodology:
                            prompt_text += " (or 'i<N>' for methodology details, e.g. 'i1')"
                        choice_raw = prompt_with_context(
                            processor=self.processor,
                            prompt=prompt_text,
                            default=default_choice,
                            module="Topology Generator",
                            description="Select forcefield set",
                        )
                        choice_str = choice_raw.strip().lower()

                        # Handle `i<N>` — show methodology for option N, then re-prompt
                        if choice_str.startswith('i') and len(choice_str) > 1:
                            try:
                                methodology_idx = int(choice_str[1:]) - 1
                                if 0 <= methodology_idx < len(options):
                                    opt = options[methodology_idx]
                                    self._show_ff_methodology_panel(
                                        opt.get('display_name') or opt['name'], opt,
                                    )
                                    continue
                                else:
                                    console.print(f"[yellow]No option {choice_str[1:]} — pick 1..{len(options)}[/yellow]")
                                    continue
                            except ValueError:
                                console.print(f"[yellow]Could not parse '{choice_raw}' — use 1..{len(options)} or i<N>[/yellow]")
                                continue

                        # Otherwise treat as a numeric selection
                        try:
                            choice_idx = int(choice_str) - 1
                            if 0 <= choice_idx < len(options):
                                chosen_option = dict(options[choice_idx])
                                chosen_option['cofactor_path'] = transformer_class.FORCEFIELD_PATH
                                selected_forcefields[key] = chosen_option
                                break
                            else:
                                console.print(f"[yellow]Pick 1..{len(options)}[/yellow]")
                        except ValueError:
                            console.print(f"[yellow]Could not parse '{choice_raw}' — use 1..{len(options)} or i<N>[/yellow]")

            except Exception as e:
                console.print(f"[red]Error finding forcefields for {transformer_type}: {e}[/red]")
                continue

        return selected_forcefields

    def _build_single_state_atom_types_section(self, requirements, selected_forcefields):
        """Build atom types section for single state template.

        Combines two sources:
          1. Each picked FF set's metadata.json ``atom_types`` (unchanged types
             — for renamed sets, these don't include the renamed codes since
             the renamed codes come from upstream parameterizations like GAFF2
             and aren't enumerated in metadata).
          2. The resolver's emitted addAtomTypes lines for renamed codes
             (workspace key ``ff_resolved_atom_types``).
        """
        all_atom_types = []
        seen = set()

        for key, info in requirements.items():
            if key in selected_forcefields:
                for entry in info.get('atom_types', []):
                    if entry not in seen:
                        all_atom_types.append(entry)
                        seen.add(entry)

        for entry in self.get_from_workspace("ff_resolved_atom_types", []) or []:
            if entry not in seen:
                all_atom_types.append(entry)
                seen.add(entry)

        if not all_atom_types:
            return "# No custom atom types needed"

        section_lines = [
            "# Custom atom types (auto-generated by ProPrep)",
            "addAtomTypes {"
        ]
        for atom_type in sorted(all_atom_types):
            section_lines.append(f"    {atom_type}")
        section_lines.append("}")

        return '\n'.join(section_lines)

    def _verify_no_ff_collisions_or_abort(self):
        """Last-line-of-defense check: parse every frcmod/lib that will be
        loaded by the upcoming tleap script and refuse to proceed if any pair
        conflicts on shared atom-type / bonded-term keys with different values.

        Raises ``RuntimeError`` with a human-readable diagnostic if any
        collision is found. Silent on success.
        """
        from proprep.ff_compat.matrix import compare_signatures
        from proprep.ff_compat.parser import parse_set

        # Collect every (label, frcmod, libs[]) tuple the generator will emit.
        units: List[tuple] = []

        # Preprocessing-derived bundles (e.g., MCPB).
        pre_libs = self.get_from_workspace("preprocessing_lib_files", []) or []
        pre_frcmods = self.get_from_workspace("preprocessing_frcmod_files", []) or []
        # The preprocessing arrays aren't keyed by source; treat each frcmod as
        # its own unit, pairing it with no specific lib (the lib files are
        # parsed separately below). This over-counts a little but only matters
        # when the lib introduces new type usage that the frcmod doesn't list.
        for i, frcmod in enumerate(pre_frcmods):
            try:
                # Pair each frcmod with all preprocessing lib files since we
                # don't have a 1:1 mapping; parse_set tolerates extra libs.
                libs = [lf for lf in pre_libs if not str(lf).endswith(".mol2")]
                units.append((f"preprocessing:{Path(frcmod).name}", frcmod, libs))
            except Exception:
                continue

        # Picker-selected FF sets.
        selected = self.get_from_workspace("single_state_selected_forcefields", {}) or {}
        for key, entry in selected.items():
            libs = entry.get('lib', [])
            if isinstance(libs, str):
                libs = [libs]
            label = entry.get('name') or f"selected:{key}"
            units.append((label, entry.get('frcmod'), libs))

        # Parse each unit's signature. Skip silently on parse error — the
        # safety net should never crash a real workflow over an unparseable
        # frcmod; that just means we can't check it.
        signatures = {}
        for label, frcmod, libs in units:
            if not frcmod:
                continue
            try:
                signatures[label] = parse_set(label, frcmod, libs or [])
            except Exception:
                continue

        labels = sorted(signatures)
        collisions = []
        for i, a in enumerate(labels):
            for b in labels[i + 1:]:
                report = compare_signatures(signatures[a], signatures[b])
                if report:
                    collisions.append((a, b, report))

        if not collisions:
            return

        console = self.processor.console
        console.print()
        console.print(
            "[red]✗ Refusing to write the tleap script: residual atom-type "
            "collisions detected across the files-to-load.[/red]"
        )
        console.print(
            "[grey50]Either the FF picker's resolver didn't run on these inputs, "
            "or a non-bundled FF (preprocessing-derived) shares parameter keys "
            "with a bundled set. Re-run the FF picker so the resolver can "
            "rewrite the colliding sets, or build the resolver bundle manually "
            "with `python -m proprep.ff_compat.resolve`.[/grey50]"
        )
        for a, b, report in collisions:
            from collections import Counter
            by_section = Counter(e["section"] for e in report["entries"])
            shared = ", ".join(report["shared_types"][:6])
            if len(report["shared_types"]) > 6:
                shared += f", … (+{len(report['shared_types']) - 6} more)"
            console.print(f"  • {a}  vs  {b}")
            console.print(f"      shared types: {shared}")
            console.print(
                f"      differing entries: {len(report['entries'])} "
                f"({', '.join(f'{s}={n}' for s, n in sorted(by_section.items()))})"
            )

        raise RuntimeError(
            f"FF collision safety net tripped: {len(collisions)} pair(s) of "
            f"to-be-loaded FF sets still conflict after the picker resolver."
        )

    def _resolve_ff_collisions(self, selected_forcefields):
        """Detect + resolve atom-type collisions among the picked FF sets.

        For each pair of selections, parse signatures and find shared atom-type
        / bonded-term keys with different values. If any pair conflicts, prompt
        the user to designate which side's atom types should be renamed and run
        the resolver. Returns ``(rewritten_selected_forcefields, atom_type_lines)``
        where the dict's frcmod/lib paths point at the rewritten files (when a
        side was renamed) and ``atom_type_lines`` is the list of resolver-emitted
        addAtomTypes entries to inject into the topology's atom-types block.

        Raises if a conflict can't be resolved (e.g., no rename target chosen
        in batch mode).
        """
        from proprep.ff_compat.matrix import compare_signatures
        from proprep.ff_compat.parser import parse_set
        from proprep.ff_compat.resolver import resolve
        from proprep.utils.prompts import prompt_with_context

        if not selected_forcefields:
            return selected_forcefields, []

        console = self.processor.console

        # Set 2 single-cofactor guard. The fragment-typed cofactor library
        # ships two FF sets per residue: Set 1 (production minimalist; multi-
        # cofactor safe by construction) and Set 2 (full QM Hessian; single-
        # cofactor only because residue-specific bond force constants on shared
        # GAFF2 types would silently clobber under tleap's last-loaded-wins).
        # Block any Set-2 pick when more than one cofactor residue is in the
        # topology — the user must downgrade to Set 1 or reduce the count.
        set2_keys = [
            key for key, entry in selected_forcefields.items()
            if entry.get('name') == 'Set2'
        ]
        if set2_keys and len(selected_forcefields) > 1:
            console.print()
            console.print(
                "[red]✗ Set 2 (full QM Hessian) is single-cofactor only.[/red]"
            )
            console.print(
                "[grey50]Set 2 derives every bond/angle from the whole-residue QM Hessian, "
                "including shared GAFF2 types (c3-c3, OS-P, etc.). Two Set-2 cofactors in "
                "one topology would silently let one residue's parameter values clobber "
                "the other's via tleap's last-loaded-wins rule. Pick Set 1 for at least "
                f"{len(selected_forcefields) - 1} of these, or reduce the cofactor count to 1.[/grey50]"
            )
            console.print()
            console.print("Selections that picked Set 2:")
            for key in set2_keys:
                console.print(f"  • {key} → {selected_forcefields[key].get('cofactor_path', '?')}")
            raise RuntimeError(
                f"Set 2 selected for {len(set2_keys)} cofactor(s) in a "
                f"{len(selected_forcefields)}-cofactor topology — Set 2 is single-cofactor only."
            )

        if len(selected_forcefields) < 2:
            return selected_forcefields, []

        # Build set_id per picker entry. Picker keys are
        # (transformer_type, redox_state, spin_state); each value carries
        # cofactor_path + name + frcmod + lib.
        sid_for_key = {}
        descriptors_inline = {}
        for key, entry in selected_forcefields.items():
            cofactor_path = entry.get('cofactor_path')
            set_name = entry.get('name')
            if not cofactor_path or not set_name:
                # Without a cofactor path we can't form a set_id; skip — the
                # picker entry must be for something outside the bundled FF
                # hierarchy (e.g. MCPB-derived) and the resolver doesn't apply.
                continue
            transformer_type, redox_state, spin_state = key
            sid = f"{cofactor_path}:{redox_state}:{spin_state}:{set_name}"
            sid_for_key[key] = sid
            libs = entry['lib'] if isinstance(entry['lib'], list) else [entry['lib']]
            descriptors_inline[sid] = {
                "frcmod": entry['frcmod'],
                "libs": libs,
            }

        if len(sid_for_key) < 2:
            return selected_forcefields, []

        signatures = {
            sid: parse_set(sid, descriptors_inline[sid]["frcmod"], descriptors_inline[sid]["libs"])
            for sid in sid_for_key.values()
        }

        sorted_sids = sorted(set(sid_for_key.values()))
        conflicts: List[tuple] = []
        for i, a in enumerate(sorted_sids):
            for b in sorted_sids[i + 1:]:
                report = compare_signatures(signatures[a], signatures[b])
                if report:
                    conflicts.append((a, b, report))

        if not conflicts:
            return selected_forcefields, []

        console.print()
        console.print(
            f"[yellow]⚠ FF parameter collisions across {len(conflicts)} pair(s) of "
            f"selected FF sets[/yellow]"
        )
        console.print(
            "[grey50]Same atom-type code (or bonded-term key) declared with different "
            "values in two selected sets. tleap silently lets the last-loaded "
            "definition win, which corrupts the resulting prmtop. ProPrep will "
            "rename the loser-side codes so both sets can coexist.[/grey50]"
        )
        console.print()

        rename_choices: Dict[str, str] = {}
        # Track per-set rename choices already made so we don't re-prompt for
        # the same set in subsequent pairs (auto-propagate the user's intent).
        renamed_so_far: set = set()
        for a, b, report in conflicts:
            shared = ", ".join(report["shared_types"][:8])
            if len(report["shared_types"]) > 8:
                shared += f", … (+{len(report['shared_types']) - 8} more)"
            from collections import Counter
            by_section = Counter(e["section"] for e in report["entries"])
            section_summary = ", ".join(
                f"{sec}={n}" for sec, n in sorted(by_section.items())
            )
            console.print(f"[bold]Conflict:[/bold]")
            console.print(f"  (a) {a}")
            console.print(f"  (b) {b}")
            console.print(f"  Shared types: {shared or '(name space disjoint; conflict via shared keys)'}")
            console.print(f"  Differing entries: {len(report['entries'])} ({section_summary})")

            # Auto-propagate prior choices: if exactly one side has already
            # been chosen for renaming in another pair, default to that side.
            if a in renamed_so_far and b not in renamed_so_far:
                chosen = a
                console.print(f"  [grey50](a) was chosen to be renamed in an earlier pair — applying same here.[/grey50]")
            elif b in renamed_so_far and a not in renamed_so_far:
                chosen = b
                console.print(f"  [grey50](b) was chosen to be renamed in an earlier pair — applying same here.[/grey50]")
            else:
                while True:
                    choice = prompt_with_context(
                        processor=self.processor,
                        prompt="Which set's atom types should be renamed?",
                        default="a",
                        module="Topology Generator",
                        description="Resolve FF parameter collision",
                        options_map={"a": f"rename {a}", "b": f"rename {b}"},
                    ).strip().lower()
                    if choice in ("a", "1"):
                        chosen = a
                        break
                    if choice in ("b", "2"):
                        chosen = b
                        break
                    console.print("[yellow]Please type 'a' or 'b'.[/yellow]")
            renamed_so_far.add(chosen)
            rename_choices[f"{a}|{b}"] = chosen
            console.print()

        # Determine where to write rewritten files: the session output_dir.
        workspace = self.processor._get_workspace()
        output_dir = Path(workspace.get("output_dir", ".")).expanduser().resolve()

        # The bundled resolver's discover_all_sets lookup expects set_ids
        # discoverable from the bundled/user FF hierarchy. Our inline-built
        # set_ids match that scheme as long as the picker entries came from
        # the bundled tree — which is the case for redox-site picks. Use the
        # standard resolve() entry point.
        result = resolve(
            set_ids=sorted_sids,
            rename_choices=rename_choices,
            workspace=output_dir,
            prompt_for_pair=None,
        )

        console.print(f"[green]✓ Resolved into {result.workspace_dir}[/green]")
        for sid, so in result.set_outputs.items():
            if so.was_renamed:
                renames = ", ".join(f"{k}→{v}" for k, v in sorted(so.rename_map.items()))
                console.print(f"  [renamed] {sid}: {renames}")

        # Rewrite the picker entries to point at the resolver's outputs.
        sid_to_output = result.set_outputs
        new_selected = dict(selected_forcefields)
        for key, sid in sid_for_key.items():
            so = sid_to_output.get(sid)
            if so is None or not so.was_renamed:
                continue
            entry = dict(selected_forcefields[key])
            entry["frcmod"] = str(so.frcmod_path)
            new_lib = [str(p) for p in so.lib_paths]
            entry["lib"] = new_lib if len(new_lib) > 1 else new_lib[0]
            new_selected[key] = entry

        # Collect resolver-emitted addAtomTypes lines for the atom-types block.
        atom_type_lines: List[str] = []
        for so in sid_to_output.values():
            atom_type_lines.extend(so.add_atom_types_entries)

        return new_selected, atom_type_lines
    
    def _build_forcefield_parameters_section(self):
        """Build the forcefield parameters section for template substitution.

        Includes parameters from two sources:
        1. Preprocessing: lib/frcmod files for small molecules (preprocessing_lib_files, preprocessing_frcmod_files)
        2. Redox site parameterization: lib/frcmod for redox states (single_state_selected_forcefields)
        """
        loaded_files = set()
        section_lines = []

        # 1. Load preprocessing lib/frcmod files (from structure_preprocessor)
        preprocessing_lib = self.get_from_workspace("preprocessing_lib_files", [])
        preprocessing_frcmod = self.get_from_workspace("preprocessing_frcmod_files", [])

        if preprocessing_lib or preprocessing_frcmod:
            section_lines.append("# Small molecule parameters (from preprocessing)")
            # Source the ligand base force field so tLEaP has the atom-type
            # DEFINITIONS (element, hybridization, GB radii) for the GAFF types
            # stored in the ligand lib. Without this, tLEaP emits non-fatal
            # "UNKNOWN ATOM TYPE" warnings and then SILENTLY mis-assigns the
            # element (e.g. Cl -> carbon in ATOMIC_NUMBER) and defaults every GB
            # radius/screen value. ProPrep's small-molecule parameterizer uses
            # GAFF2, so source it here.
            if "leaprc.gaff2" not in loaded_files:
                section_lines.append("source leaprc.gaff2")
                loaded_files.add("leaprc.gaff2")
            for lib_file in preprocessing_lib:
                if lib_file not in loaded_files:
                    # mol2 files use loadmol2, lib files use loadoff.
                    # Filenames are quoted: a ligand/cofactor lib named after a
                    # digit-leading residue code (e.g. 9E2.mol2) is otherwise
                    # lexed by tLEaP as a scientific-notation number.
                    if lib_file.endswith('.mol2'):
                        # A bare `loadmol2 <file>` registers NO template — the
                        # unit must be assigned to a variable. That variable name
                        # becomes the template's match-name for loadpdb, so derive
                        # it from the mol2's own residue name (kept tLEaP-safe).
                        unit = tleap_safe_unit_var(self._read_mol2_resname(lib_file) or "MOL")
                        section_lines.append(f'{unit} = loadmol2 "{lib_file}"')
                    else:
                        section_lines.append(f'loadoff "{lib_file}"')
                    loaded_files.add(lib_file)
            for frcmod_file in preprocessing_frcmod:
                if frcmod_file not in loaded_files:
                    section_lines.append(f'loadamberparams "{frcmod_file}"')
                    loaded_files.add(frcmod_file)

        # 2. Load redox site parameterization forcefields (from configure step)
        selected_forcefields = self.get_from_workspace("single_state_selected_forcefields", {})

        if selected_forcefields:
            if section_lines:
                section_lines.append("")  # Blank line separator
            section_lines.append("# Redox site parameters (from parameterization)")

            # 2a. Emit any cofactor-prereq leaprcs the user explicitly confirmed
            # at Layer 3 of the FF picker, with a comment trail explaining why.
            selected_standard = self.get_from_workspace("selected_standard_forcefields", {}) or {}
            confirmed_extras = selected_standard.get("cofactor_prereqs", []) or []
            forced_skips = selected_standard.get("cofactor_prereqs_skipped", []) or []
            for leaprc in confirmed_extras:
                if leaprc not in loaded_files:
                    section_lines.append(
                        f"# {leaprc} — confirmed at FF-picker as a cofactor prerequisite"
                    )
                    section_lines.append(f"source {leaprc}")
                    loaded_files.add(leaprc)
            for leaprc in forced_skips:
                section_lines.append(
                    f"# WARNING: {leaprc} required by a selected redox-site transformer "
                    f"but the user explicitly chose to skip it. tleap will fail."
                )

            for key, ff_info in selected_forcefields.items():
                # Load .lib files first (may be a list for MCPB multi-mol2 sites)
                lib_files = ff_info['lib']
                if isinstance(lib_files, list):
                    for lib_file in lib_files:
                        if lib_file not in loaded_files:
                            if lib_file.endswith('.mol2'):
                                section_lines.append(f'loadmol2 "{lib_file}"')
                            else:
                                section_lines.append(f'loadoff "{lib_file}"')
                            loaded_files.add(lib_file)
                else:
                    if lib_files not in loaded_files:
                        section_lines.append(f'loadoff "{lib_files}"')
                        loaded_files.add(lib_files)

                # Then load .frcmod files (may be a list: a metal site has a
                # bonded frcmod plus each organic ligand's own GAFF frcmod)
                frcmod_files = ff_info['frcmod']
                if isinstance(frcmod_files, list):
                    for frcmod_file in frcmod_files:
                        if frcmod_file not in loaded_files:
                            section_lines.append(f'loadamberparams "{frcmod_file}"')
                            loaded_files.add(frcmod_file)
                else:
                    if frcmod_files not in loaded_files:
                        section_lines.append(f'loadamberparams "{frcmod_files}"')
                        loaded_files.add(frcmod_files)

        if not section_lines:
            return "# No custom forcefield parameters needed"

        # Add header
        section_lines.insert(0, "# Custom forcefield parameters (auto-generated by ProPrep)")

        return '\n'.join(section_lines)
    
    def _build_atom_types_section(self):
        """Build the atom types section for single state template substitution"""
        # Get stored requirements and selected forcefields
        requirements = self.get_from_workspace("single_state_ff_requirements", {})
        selected_forcefields = self.get_from_workspace("single_state_selected_forcefields", {})

        result = self._build_single_state_atom_types_section(requirements, selected_forcefields)

        # Fallback: use atom types from MCPB preprocessing
        if result == "# No custom atom types needed":
            preprocessing_types = self.get_from_workspace("preprocessing_atom_types", [])
            if preprocessing_types:
                lines = [
                    "# Custom atom types (from MCPB preprocessing)",
                    "addAtomTypes {"
                ]
                for entry in preprocessing_types:
                    lines.append(f"    {entry}")
                lines.append("}")
                return '\n'.join(lines)

        return result
    
    def _build_bond_definitions_section(self):
        """Build the bond definitions section for template substitution"""
        has_bonds = self.combined_bonds and any(len(bonds) > 0 for bonds in self.combined_bonds.values())

        # Fallback: use pre-converted bond commands from preprocessing
        if not has_bonds:
            preprocessing_bonds = self.get_from_workspace("preprocessing_bond_commands", [])
            if preprocessing_bonds:
                section_lines = ["# Bond definitions (from preprocessing)"]
                section_lines.extend(preprocessing_bonds)
                return '\n'.join(section_lines)
            return "# No bond definitions found"

        section_lines = []
        
        section_lines.append("# Bond definitions (auto-generated by ProPrep)")
        section_lines.append("")
        
        # Write covalent bonds
        if self.combined_bonds.get("covalent"):
            section_lines.append("# Covalent bonds (non-metal to non-metal)")
            for bond in self.combined_bonds["covalent"]:
                section_lines.append(bond)
            section_lines.append("")

        # Write disulfide bonds (SG-SG between CYS/CYX). Previously omitted
        # here, so disulfide-categorized bonds were silently dropped from the
        # single-state tLEaP script even though they were detected and stored.
        if self.combined_bonds.get("disulfide"):
            section_lines.append("# Disulfide bonds (CYS/CYX SG-SG)")
            for bond in self.combined_bonds["disulfide"]:
                section_lines.append(bond)
            section_lines.append("")

        # Write coordination bonds
        if self.combined_bonds.get("coordinate"):
            section_lines.append("# Coordination bonds (metal to non-metal)")
            for bond in self.combined_bonds["coordinate"]:
                section_lines.append(bond)
            section_lines.append("")

        # Write metal-metal bonds. Category key is "metal-metal" (hyphen) to
        # match the detector and combined_bonds initialization; the previous
        # "metal_metal" (underscore) lookup never matched, dropping these bonds.
        if self.combined_bonds.get("metal-metal"):
            section_lines.append("# Metal-metal bonds")
            for bond in self.combined_bonds["metal-metal"]:
                section_lines.append(bond)
            section_lines.append("")
        
        # Write peptide backbone bonds
        if self.combined_bonds.get("peptide_backbone"):
            section_lines.append("# Peptide backbone bonds")
            for bond in self.combined_bonds["peptide_backbone"]:
                section_lines.append(bond)
            section_lines.append("")
        
        
        # Write other bonds
        if self.combined_bonds.get("other"):
            section_lines.append("# Other custom bonds")
            for bond in self.combined_bonds["other"]:
                section_lines.append(bond)
            section_lines.append("")
        
        return '\n'.join(section_lines)
    
    def _generate_microstate_tleap_template(self, microstate_code="MICROSTATE"):
        """Generate tLEaP template specifically for microstate processing"""
        # Build forcefield section from user selections (shared with single state)
        forcefield_section, water_box = self._build_standard_forcefield_section()

        template = f"""# ProPrep-generated tLEaP Input for Microstate: {microstate_code}

# === STANDARD FORCEFIELDS ===
{forcefield_section}

# === CUSTOM ATOM TYPES (auto-filled by ProPrep) ===
# ATOM_TYPES_SECTION

# === CUSTOM PARAMETERS (auto-filled by ProPrep) ===
# FORCEFIELD_PARAMETERS_SECTION

# === LOAD STRUCTURE (auto-filled by ProPrep) ===
# PDB_FILE_SECTION

# === BOND DEFINITIONS (auto-filled by ProPrep) ===
# BOND_DEFINITIONS_SECTION

# === SOLVATION (auto-filled by ProPrep) ===
# SOLVATION_SECTION

# === VALIDATION ===
check mol

# === OUTPUT ===
saveamberparm mol {microstate_code}.prmtop {microstate_code}.rst7

quit"""
        return template
    
    def _build_microstate_atom_types_section(self, microstate_info, selected_forcefields):
        """Build atom types section for microstate template substitution"""
        # Collect and deduplicate atom types
        all_atom_types = set()
        for site in microstate_info['sites']:
            key = (site['transformer_type'], site['redox_state'], site['spin_state'])
            if key in selected_forcefields:
                atom_types = site['forcefield_info'].get('atom_types', [])
                all_atom_types.update(atom_types)
        
        if not all_atom_types:
            return "# No custom atom types needed"
        
        section_lines = [
            "# Custom atom types (auto-generated by ProPrep)",
            "addAtomTypes {"
        ]
        for atom_type in sorted(all_atom_types):
            section_lines.append(f"    {atom_type}")
        section_lines.append("}")
        
        return '\n'.join(section_lines)
    
    def _build_microstate_forcefield_section(self, microstate_info, selected_forcefields):
        """Build forcefield parameters section for microstate template substitution.

        Includes parameters from two sources:
        1. Preprocessing: lib/frcmod files for small molecules (preprocessing_lib_files, preprocessing_frcmod_files)
        2. Microstate-specific: lib/frcmod for redox states (selected_forcefields)
        """
        loaded_files = set()
        section_lines = []

        # 1. Load preprocessing lib/frcmod files (from structure_preprocessor)
        preprocessing_lib = self.get_from_workspace("preprocessing_lib_files", [])
        preprocessing_frcmod = self.get_from_workspace("preprocessing_frcmod_files", [])

        if preprocessing_lib or preprocessing_frcmod:
            section_lines.append("# Small molecule parameters (from preprocessing)")
            # Source the ligand base force field so tLEaP has the atom-type
            # DEFINITIONS (element, hybridization, GB radii) for the GAFF types
            # stored in the ligand lib. Without this, tLEaP emits non-fatal
            # "UNKNOWN ATOM TYPE" warnings and then SILENTLY mis-assigns the
            # element (e.g. Cl -> carbon in ATOMIC_NUMBER) and defaults every GB
            # radius/screen value. ProPrep's small-molecule parameterizer uses
            # GAFF2, so source it here.
            if "leaprc.gaff2" not in loaded_files:
                section_lines.append("source leaprc.gaff2")
                loaded_files.add("leaprc.gaff2")
            for lib_file in preprocessing_lib:
                if lib_file not in loaded_files:
                    # mol2 files use loadmol2, lib files use loadoff.
                    # Filenames are quoted: a ligand/cofactor lib named after a
                    # digit-leading residue code (e.g. 9E2.mol2) is otherwise
                    # lexed by tLEaP as a scientific-notation number.
                    if lib_file.endswith('.mol2'):
                        # A bare `loadmol2 <file>` registers NO template — the
                        # unit must be assigned to a variable. That variable name
                        # becomes the template's match-name for loadpdb, so derive
                        # it from the mol2's own residue name (kept tLEaP-safe).
                        unit = tleap_safe_unit_var(self._read_mol2_resname(lib_file) or "MOL")
                        section_lines.append(f'{unit} = loadmol2 "{lib_file}"')
                    else:
                        section_lines.append(f'loadoff "{lib_file}"')
                    loaded_files.add(lib_file)
            for frcmod_file in preprocessing_frcmod:
                if frcmod_file not in loaded_files:
                    section_lines.append(f'loadamberparams "{frcmod_file}"')
                    loaded_files.add(frcmod_file)

        # 2. Load microstate-specific forcefields
        has_microstate_ffs = False
        for site in microstate_info['sites']:
            key = (site['transformer_type'], site['redox_state'], site['spin_state'])
            if key in selected_forcefields:
                if not has_microstate_ffs:
                    if section_lines:
                        section_lines.append("")  # Blank line separator
                    section_lines.append("# Redox site parameters (from parameterization)")
                    has_microstate_ffs = True

                ff_info = selected_forcefields[key]

                # Load .lib files first (may be a list for MCPB multi-mol2 sites)
                lib_files = ff_info['lib']
                if isinstance(lib_files, list):
                    for lib_file in lib_files:
                        if lib_file not in loaded_files:
                            if lib_file.endswith('.mol2'):
                                section_lines.append(f'loadmol2 "{lib_file}"')
                            else:
                                section_lines.append(f'loadoff "{lib_file}"')
                            loaded_files.add(lib_file)
                else:
                    if lib_files not in loaded_files:
                        section_lines.append(f'loadoff "{lib_files}"')
                        loaded_files.add(lib_files)

                # Then load .frcmod files (may be a list: a metal site has a
                # bonded frcmod plus each organic ligand's own GAFF frcmod)
                frcmod_files = ff_info['frcmod']
                if isinstance(frcmod_files, list):
                    for frcmod_file in frcmod_files:
                        if frcmod_file not in loaded_files:
                            section_lines.append(f'loadamberparams "{frcmod_file}"')
                            loaded_files.add(frcmod_file)
                else:
                    if frcmod_files not in loaded_files:
                        section_lines.append(f'loadamberparams "{frcmod_files}"')
                        loaded_files.add(frcmod_files)

        if not section_lines:
            return "# No custom forcefield parameters needed"

        # Add header
        section_lines.insert(0, "# Custom forcefield parameters (auto-generated by ProPrep)")

        return '\n'.join(section_lines)
    
    def _build_microstate_pdb_section(self, microstate_info):
        """Build PDB loading section for microstate template substitution"""
        return f"""# Load structure (auto-generated by ProPrep)
mol = loadpdb {microstate_info['filename']}"""
    
    def _build_microstate_bond_section(self, microstate_info):
        """Build bond definitions section for microstate template substitution"""
        bond_commands = self._get_bond_commands_for_microstate(microstate_info)

        if not bond_commands:
            # Fallback: use pre-converted bond commands from preprocessing
            preprocessing_bonds = self.get_from_workspace("preprocessing_bond_commands", [])
            if preprocessing_bonds:
                section_lines = ["# Bond definitions (from preprocessing)"]
                section_lines.extend(preprocessing_bonds)
                return '\n'.join(section_lines)
            return "# No bond definitions needed"

        section_lines = ["# Bond definitions (auto-generated by ProPrep)"]
        section_lines.extend(bond_commands)

        return '\n'.join(section_lines)
    
    def _display_bond_dashboard(self, filter_category=None):
        """Display the bond editing dashboard"""
        console = self.processor.console
        
        # Header
        console.print("\n[bold blue]═══ tLEaP Bond Editor Dashboard ═══[/bold blue]")
        
        # Current filter status
        if filter_category:
            category_names = {
                "covalent": "Covalent Bonds",
                "disulfide": "Disulfide Bonds",
                "coordinate": "Coordination Bonds",
                "metal-metal": "Metal-Metal Bonds",
                "peptide_backbone": "Peptide Backbone",
                "other": "Other Bonds"
            }
            filter_name = category_names.get(filter_category, filter_category)
            console.print(f"[yellow]Showing: {filter_name}[/yellow]")
        else:
            console.print("[yellow]Showing: All Bond Categories[/yellow]")
        
        # Bond display
        bond_counter = 1
        categories_to_show = [filter_category] if filter_category else self.combined_bonds.keys()
        
        for category in categories_to_show:
            if category not in self.combined_bonds:
                continue
                
            bonds = self.combined_bonds[category]
            if not bonds:
                continue
                
            # Category header
            category_names = {
                "covalent": "🔗 Covalent Bonds",
                "disulfide": "🔗 Disulfide Bonds (SG-SG)",
                "coordinate": "⚙️  Coordination Bonds",
                "metal-metal": "🔩 Metal-Metal Bonds",
                "peptide_backbone": "🧬 Peptide Backbone",
                "other": "🔧 Other Bonds"
            }
            category_display = category_names.get(category, category.replace("_", " ").title())
            
            if not filter_category:  # Only show category headers when showing all
                console.print(f"\n[bold]{category_display}[/bold]")
                
            # List bonds
            for bond in bonds:
                console.print(f"  [bright_white]{bond_counter:2d}.[/bright_white] {bond}")
                bond_counter += 1
        
        if bond_counter == 1:
            console.print("\n[grey50]No bonds to display[/grey50]")
        
        # Commands help
        console.print("\n[bold blue]Commands:[/bold blue]")
        console.print("[grey50]Filters:[/grey50] \\[d]isulfide \\[m]etal \\[p]eptide \\[o]ther \\[all] categories")
        console.print("[grey50]Actions:[/grey50] \\[a]dd bond, \\[#] edit bond, \\[x#] delete bond, \\[s]ave, \\[q]uit, \\[h]elp")

    def _add_bond_interactive(self):
        """Interactively add a new bond"""
        console = self.processor.console
        
        console.print("\n[bold blue]Add New Bond[/bold blue]")
        
        # Get bond command
        bond_cmd = prompt_with_context(
            self.processor,
            "Enter tLEaP bond command (e.g. 'bond mol.123.SG mol.456.SG')",
            module="Topology Generator",
            description="New tLEaP bond command",
        ).strip()
        
        if not bond_cmd:
            console.print("[yellow]Cancelled[/yellow]")
            return
            
        # Validate basic format
        if not bond_cmd.startswith("bond "):
            console.print("[red]Bond command must start with 'bond '[/red]")
            return
            
        # Get category
        console.print("\nSelect category:")
        categories = {
            "1": ("covalent", "Covalent Bonds"),
            "2": ("disulfide", "Disulfide Bonds (SG-SG)"),
            "3": ("coordinate", "Coordination Bonds"),
            "4": ("metal-metal", "Metal-Metal Bonds"),
            "5": ("peptide_backbone", "Peptide Backbone"),
            "6": ("other", "Other Bonds")
        }
        
        for key, (_, name) in categories.items():
            console.print(f"  {key}. {name}")
            
        category_options_map = {k: name for k, (_, name) in categories.items()}
        choice = prompt_with_context(
            self.processor,
            "Category",
            choices=list(categories.keys()),
            default="6",
            module="Topology Generator",
            description="Bond category for new bond",
            options_map=category_options_map,
        )
        category_key, category_name = categories[choice]
        
        # Check for duplicates
        if category_key not in self.combined_bonds:
            self.combined_bonds[category_key] = []
            
        if bond_cmd in self.combined_bonds[category_key]:
            console.print(f"[yellow]Bond already exists in {category_name}[/yellow]")
            return
            
        # Add bond
        self.combined_bonds[category_key].append(bond_cmd)
        console.print(f"[green]Added bond to {category_name}[/green]")
        prompt_with_context(
            self.processor,
            "Press Enter to continue...",
            default="",
            module="Topology Generator",
            description="Pause after adding bond",
        )

    def _edit_bond_by_number(self, bond_num, filter_category):
        """Edit a bond by its display number"""
        console = self.processor.console
        
        # Find the bond by number
        bond_counter = 1
        target_bond = None
        target_category = None
        
        categories_to_search = [filter_category] if filter_category else self.combined_bonds.keys()
        
        for category in categories_to_search:
            if category not in self.combined_bonds:
                continue
                
            for bond in self.combined_bonds[category]:
                if bond_counter == bond_num:
                    target_bond = bond
                    target_category = category
                    break
                bond_counter += 1
            
            if target_bond:
                break
        
        if not target_bond:
            console.print(f"[red]Bond #{bond_num} not found[/red]")
            prompt_with_context(
                self.processor,
                "Press Enter to continue...",
                default="",
                module="Topology Generator",
                description="Pause after bond-not-found error",
            )
            return

        console.print(f"\n[bold blue]Edit Bond #{bond_num}[/bold blue]")
        console.print(f"Current: [white]{target_bond}[/white]")

        # Get new bond command
        new_bond = prompt_with_context(
            self.processor,
            "New bond command",
            default=target_bond,
            module="Topology Generator",
            description=f"New command text for bond #{bond_num}",
        ).strip()

        if new_bond and new_bond != target_bond:
            # Replace the bond
            bond_index = self.combined_bonds[target_category].index(target_bond)
            self.combined_bonds[target_category][bond_index] = new_bond
            console.print("[green]Bond updated[/green]")
        else:
            console.print("[yellow]No changes made[/yellow]")

        prompt_with_context(
            self.processor,
            "Press Enter to continue...",
            default="",
            module="Topology Generator",
            description="Pause after bond edit",
        )

    def _delete_bond_by_number(self, bond_num, filter_category):
        """Delete a bond by its display number"""
        console = self.processor.console
        
        # Find the bond by number
        bond_counter = 1
        target_bond = None
        target_category = None
        
        categories_to_search = [filter_category] if filter_category else self.combined_bonds.keys()
        
        for category in categories_to_search:
            if category not in self.combined_bonds:
                continue
                
            for bond in self.combined_bonds[category]:
                if bond_counter == bond_num:
                    target_bond = bond
                    target_category = category
                    break
                bond_counter += 1
            
            if target_bond:
                break
        
        if not target_bond:
            console.print(f"[red]Bond #{bond_num} not found[/red]")
            prompt_with_context(
                self.processor,
                "Press Enter to continue...",
                default="",
                module="Topology Generator",
                description="Pause after bond-not-found error",
            )
            return

        console.print(f"\n[bold red]Delete Bond #{bond_num}[/bold red]")
        console.print(f"Bond: [white]{target_bond}[/white]")

        if confirm_with_context(
            self.processor,
            "Confirm deletion?",
            default=False,
            module="Topology Generator",
            description=f"Confirm delete bond #{bond_num}",
        ):
            self.combined_bonds[target_category].remove(target_bond)
            console.print("[green]Bond deleted[/green]")
        else:
            console.print("[yellow]Deletion cancelled[/yellow]")

        prompt_with_context(
            self.processor,
            "Press Enter to continue...",
            default="",
            module="Topology Generator",
            description="Pause after bond delete",
        )

    def _show_dashboard_help(self):
        """Show detailed help for the dashboard"""
        console = self.processor.console
        
        console.clear()
        console.print("\n[bold blue]═══ tLEaP Bond Editor Help ═══[/bold blue]")
        console.print("\n[bold]Navigation & Filtering:[/bold]")
        console.print("  [blue]d[/blue]     - Show only disulfide bonds")
        console.print("  [blue]m[/blue]     - Show only metal coordination bonds")
        console.print("  [blue]p[/blue]     - Show only peptide backbone bonds")
        console.print("  [blue]o[/blue]     - Show only other bonds")
        console.print("  [blue]all[/blue]   - Show all bond categories")
        
        console.print("\n[bold]Bond Operations:[/bold]")
        console.print("  [blue]a[/blue]     - Add a new bond interactively")
        console.print("  [blue]#[/blue]     - Edit bond number # (e.g. '5' edits bond #5)")
        console.print("  [blue]x#[/blue]    - Delete bond number # (e.g. 'x5' deletes bond #5)")
        
        console.print("\n[bold]File Operations:[/bold]")
        console.print("  [blue]s[/blue]     - Save changes to workspace")
        console.print("  [blue]q[/blue]     - Quit dashboard")
        
        console.print("\n[bold]Other:[/bold]")
        console.print("  [blue]h[/blue]     - Show this help message")
        console.print("  [blue]help[/blue]  - Show this help message")
        
        console.print("\n[bold]Bond Format:[/bold]")
        console.print("  tLEaP bonds use format: [white]bond mol.ResidueID.AtomName mol.ResidueID.AtomName[/white]")
        console.print("  Example: [grey50]bond mol.123.SG mol.456.SG[/grey50]")

        prompt_with_context(
            self.processor,
            "\nPress Enter to continue...",
            default="",
            module="Topology Generator",
            description="Pause after viewing dashboard help",
        )

    def _save_bonds(self):
        """Save the combined bonds to workspace"""
        console = self.processor.console
        
        # Update workspace
        self.update_workspace("combined_tleap_commands", self.combined_bonds)
        
        # Count total bonds
        total_bonds = sum(len(bonds) for bonds in self.combined_bonds.values())
        
        console.print(f"\n[green]✓ Saved {total_bonds} bond definitions to workspace[/green]")
        
        # Show summary
        for category, bonds in self.combined_bonds.items():
            if bonds:
                category_names = {
                    "covalent": "Covalent",
                    "disulfide": "Disulfide",
                    "coordinate": "Coordinate",
                    "metal-metal": "Metal-Metal",
                    "peptide_backbone": "Peptide",
                    "other": "Other"
                }
                name = category_names.get(category, category)
                console.print(f"  {name}: {len(bonds)} bonds")


    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # PB Titrate refinement (Poisson-Boltzmann pKa via PBSA)

    def run_pb_titrate(self):
        """Refine protonation states via Poisson-Boltzmann calculations.

        Delegates to the WorkflowChecklist defined in
        `proprep.pb_titrate.workflow`. Reads the prmtop+rst7 produced by
        the Topology Generator and the `detected_redox_sites` workspace
        key (so multi-residue MCPB groups like hemes are titrated as
        integer-charge envelopes). Writes:

          - `prmtop_titrated` / `rst7_titrated`: prmtop with recommended
            charges baked in for non-cpinutil residues.
          - `titrate_recommendations`: per-site (state_id, state_name,
            pKa_eff, populations) — read by the cpin step's "Use titrate
            recommendations" option (option 6).
          - `titrate_report_csv`: human-readable per-site summary.

        Returns True on success, False otherwise.
        """
        console = self.processor.console
        try:
            from proprep.pb_titrate.workflow import run_pb_titrate_workflow
        except ImportError as exc:
            console.print(f"[yellow]PB titrate workflow not yet available: "
                          f"{exc}[/yellow]")
            console.print("[grey50]This option is under active development. "
                          "It will run a Poisson-Boltzmann pKa refinement "
                          "across all titratable residues using AmberTools "
                          "pbsa, account for site-site coupling via "
                          "mean-field / Monte Carlo / exact enumeration, "
                          "and update the prmtop charges in place.[/grey50]")
            return False
        return run_pb_titrate_workflow(self.processor)


    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # CPIN File Generation for Constant pH MD

    def generate_cpin_file(self):
        """Generate cpin file for constant pH MD simulations"""
        console = self.processor.console
        workspace = self.get_workspace()

        console.print("\n[bold blue]═══ Generate CPIN File for Constant pH MD ═══[/bold blue]")

        # Helper to truncate long filenames for display
        def truncate_filename(fname, max_len=50):
            if len(fname) <= max_len:
                return fname
            name, ext = os.path.splitext(fname)
            available = max_len - len(ext) - 3  # 3 for "..."
            return name[:available//2] + "..." + name[-(available//2):] + ext

        # Use dynamic step counter
        step = 1

        # Check for multi-microstate case
        microstate_tleap_files = workspace.get("generated_microstate_tleap_files")
        is_multi_microstate = microstate_tleap_files and len(microstate_tleap_files) > 1

        if is_multi_microstate:
            console.print(f"\n[blue]Detected {len(microstate_tleap_files)} microstate topologies[/blue]")
            console.print("[grey50]All microstates share the same titratable residues.[/grey50]")

        # Step 1: Select simulation type FIRST (before topology selection)
        sim_type = self._select_simulation_type(step)
        if not sim_type:
            return False
        step += 1

        # Now handle topology selection based on simulation type
        if is_multi_microstate:
            # Find all valid topologies
            tleap_input_file = None
            tleap_info = None
            all_microstate_info = []

            for tleap_file in microstate_tleap_files:
                if os.path.exists(tleap_file):
                    info = self._parse_tleap_input_file(tleap_file)
                    if info and info.get('expected_prmtop') and os.path.exists(info['expected_prmtop']):
                        all_microstate_info.append({'tleap_file': tleap_file, 'info': info})
                        if tleap_input_file is None:
                            tleap_input_file = tleap_file
                            tleap_info = info

            if not tleap_input_file:
                console.print("[red]ERROR: No valid topology files found. Run 'Generate Topology Files' first.[/red]")
                return False

            prmtop_display = truncate_filename(os.path.basename(tleap_info['expected_prmtop']))

            console.print("[grey50]One CPIN file will be generated (shared by all microstates).[/grey50]")
            if sim_type == 'implicit':
                console.print(f"[green]Using representative topology:[/green] {prmtop_display}")
            else:
                console.print(f"[yellow]Explicit solvent: Modified topologies with custom radii will be generated for all {len(all_microstate_info)} microstates.[/yellow]")
                console.print(f"[grey50]Representative topology for CPIN:[/grey50] {prmtop_display}")

            # Store the active tLEaP file for coordinate mapping
            workspace.set("_active_tleap_input_file", tleap_input_file)
        else:
            all_microstate_info = None  # Not multi-microstate

            # Get tLEaP input file (use workspace if available)
            tleap_input_file = workspace.get("tleap_input_file")

            if tleap_input_file and os.path.exists(tleap_input_file):
                console.print(f"\n[green]Using tLEaP file from workspace:[/green] {os.path.basename(tleap_input_file)}")
                tleap_info = self._parse_tleap_input_file(tleap_input_file)
                if not tleap_info:
                    console.print("[red]ERROR: Failed to parse tLEaP input file from workspace.[/red]")
                    return False
            else:
                # No workspace data - prompt user
                console.print(f"\n[bold]Step {step}: Select tLEaP Input File[/bold]")
                console.print("[grey50]This file contains the topology information and residue mapping.[/grey50]")
                tleap_input_file, tleap_info = self._select_tleap_input_file_manual(step)
                if not tleap_input_file:
                    return False
                step += 1

            # Store for coordinate mapping
            workspace.set("_active_tleap_input_file", tleap_input_file)

        # Extract prmtop and rst7 from tLEaP file (no user prompt needed)
        prmtop_file = tleap_info['expected_prmtop']
        rst7_file = tleap_info['expected_rst7']

        console.print(f"[grey50]  Topology file: {truncate_filename(os.path.basename(prmtop_file))}[/grey50]")
        console.print(f"[grey50]  Coordinate file: {truncate_filename(os.path.basename(rst7_file))}[/grey50]")

        # Verify files exist
        if not os.path.exists(prmtop_file):
            console.print(f"[red]ERROR: Topology file not found: {prmtop_file}[/red]")
            console.print("[yellow]Run 'Generate Topology Files' first.[/yellow]")
            return False

        if not os.path.exists(rst7_file):
            console.print(f"[red]ERROR: Coordinate file not found: {rst7_file}[/red]")
            console.print("[yellow]Run 'Generate Topology Files' first.[/yellow]")
            return False

        # Configure GB model (for implicit solvent only)
        if sim_type == 'implicit':
            igb, intdiel = self._configure_gb_parameters(step)
            step += 2  # This method internally uses 2 steps (GB model + dielectric)
        else:
            igb, intdiel = None, 1.0

        # Step N: Select titratable residues (scan topology directly)
        selected_residues = self._select_titratable_residues(step, prmtop_file, rst7_file)
        if not selected_residues:
            console.print("[yellow]No residues selected. Aborting.[/yellow]")
            return False
        step += 1

        # Step N+1: Set initial protonation states
        initial_states = self._set_initial_protonation_states(selected_residues, step)
        step += 1

        # Step N+2: Generate CPIN file (and modified prmtop for first microstate if explicit)
        success = self._run_cpinutil(
            prmtop_file=prmtop_file,
            tleap_input_file=tleap_input_file,
            tleap_info=tleap_info,
            sim_type=sim_type,
            igb=igb,
            intdiel=intdiel,
            selected_residues=selected_residues,
            initial_states=initial_states,
            step=step
        )

        if not success:
            return False

        # For multi-microstate explicit solvent, generate modified prmtops for remaining microstates
        if is_multi_microstate and sim_type == 'explicit' and all_microstate_info and len(all_microstate_info) > 1:
            console.print(f"\n[bold]Generating modified prmtops for remaining {len(all_microstate_info) - 1} microstates...[/bold]")
            console.print("[grey50]Each microstate needs its own prmtop with corrected radii for explicit solvent CpHMD.[/grey50]\n")

            # Skip the first one (already processed above)
            for i, ms_info in enumerate(all_microstate_info[1:], 2):
                ms_prmtop = ms_info['info']['expected_prmtop']
                console.print(f"  [{i}/{len(all_microstate_info)}] Processing {os.path.basename(ms_prmtop)}...")

                # Generate modified prmtop using cpinutil (CPIN output will be discarded/overwritten)
                modified_prmtop = os.path.splitext(ms_prmtop)[0] + "_cpin.prmtop"

                # Build minimal cpinutil command just for radii modification
                cmd = [
                    "cpinutil.py",
                    "-p", ms_prmtop,
                    "-op", modified_prmtop,
                    "-o", os.devnull  # Discard CPIN output (we already have the shared one)
                ]

                # Add residue selection
                if selected_residues:
                    resnums = [str(r['resnum']) for r in selected_residues]
                    cmd.extend(["-resnums"] + resnums)

                try:
                    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
                    if result.returncode == 0:
                        console.print(f"    [green]✓ Generated: {os.path.basename(modified_prmtop)}[/green]")
                        # Swap the matching md_structure_pairs entry for this
                        # microstate so the MD manager sees the cpin-ready
                        # prmtop instead of the pre-cpinutil constph one.
                        self._swap_md_pair_to_cpin_prmtop(ms_prmtop,
                                                            modified_prmtop)
                    else:
                        console.print(f"    [red]✗ Failed: {result.stderr[:100]}[/red]")
                except Exception as e:
                    console.print(f"    [red]✗ Error: {e}[/red]")

            console.print(f"\n[green]✓ Modified prmtops generated for all {len(all_microstate_info)} microstates[/green]")

        return success

    def _select_tleap_input_file_manual(self, step):
        """Manually select tLEaP input file (fallback when not in workspace)"""
        console = self.processor.console
        workspace = self.get_workspace()

        # Search output directory for tLEaP files
        import glob
        output_dir = workspace.get("output_dir", ".")
        tleap_files = []
        for pattern in ["*.in", "*.inp", "*.leap", "*.tleap"]:
            tleap_files.extend(sorted(glob.glob(os.path.join(output_dir, pattern))))

        # Parse all tLEaP files
        parsed_files = []
        for tleap_file in tleap_files:
            info = self._parse_tleap_input_file(tleap_file)
            if info:
                parsed_files.append((tleap_file, info))

        if parsed_files:
            if len(parsed_files) == 1:
                tleap_file, tleap_info = parsed_files[0]
                console.print(f"\n[green]Found tLEaP file:[/green] {os.path.basename(tleap_file)}")
                if confirm_with_context(
                    processor=self.processor,
                    prompt="Use this tLEaP input file?",
                    default=True,
                    module="Topology Generator - CPIN",
                    description="Confirm tLEaP input file"
                ):
                    return tleap_file, tleap_info
            else:
                # Multiple files - let user choose
                from rich.table import Table
                console.print(f"\n[green]Found {len(parsed_files)} tLEaP files:[/green]")
                table = Table()
                table.add_column("#", style="blue")
                table.add_column("tLEaP File", style="green")
                table.add_column("Generates", style="yellow")

                for idx, (tleap_file, info) in enumerate(parsed_files, 1):
                    prmtop_name = os.path.basename(info['expected_prmtop']) if info['expected_prmtop'] else "N/A"
                    table.add_row(str(idx), os.path.basename(tleap_file), prmtop_name)

                console.print(table)
                choice = prompt_with_context(
                    processor=self.processor,
                    prompt="\nSelect tLEaP input file",
                    choices=[str(i) for i in range(1, len(parsed_files) + 1)],
                    default="1",
                    module="Topology Generator - CPIN",
                    description="Select tLEaP input file"
                )
                return parsed_files[int(choice) - 1]

        # No files found - prompt for path
        console.print("\n[yellow]No tLEaP input files found in output directory.[/yellow]")
        tleap_path = prompt_with_context(
            processor=self.processor,
            prompt="Enter path to tLEaP input file",
            default="",
            module="Topology Generator - CPIN",
            description="tLEaP input file path"
        )

        if not tleap_path or not os.path.exists(tleap_path):
            console.print("[red]ERROR: tLEaP input file is required for CPIN generation.[/red]")
            return None, None

        tleap_info = self._parse_tleap_input_file(tleap_path)
        if not tleap_info:
            console.print("[red]ERROR: Failed to parse tLEaP input file.[/red]")
            return None, None

        return tleap_path, tleap_info

    def _select_simulation_type(self, step):
        """Select simulation type (implicit vs explicit)"""
        console = self.processor.console

        console.print(f"\n[bold]Step {step}: Simulation Type[/bold]")
        console.print("  1. Implicit solvent (Generalized Born)")
        console.print("  2. Explicit solvent (requires custom radii in prmtop)")

        choice = prompt_with_context(
            processor=self.processor,
            prompt="\nSelect simulation type",
            choices=["1", "2"],
            default="1",
            module="Topology Generator - CPIN",
            description="Simulation type",
            options_map={"1": "Implicit solvent", "2": "Explicit solvent"}
        )

        if choice == "2":
            console.print("\n[yellow]⚠  Explicit Solvent Constant pH Setup[/yellow]")
            console.print("For explicit solvent simulations, custom radii are needed for")
            console.print("accurate carboxylate pKas (AS4/GL4 residues).")
            console.print("\ncpinutil will generate a modified prmtop with these radii.")
            console.print("[grey50]Note: Use this modified prmtop for constant pH simulations only![/grey50]")

        return 'implicit' if choice == "1" else 'explicit'

    def _configure_gb_parameters(self, step):
        """Configure GB model and internal dielectric"""
        console = self.processor.console

        console.print(f"\n[bold]Step {step}: Generalized Born Model[/bold]")
        console.print("Which GB model will you use for dynamics?")
        console.print("  1. igb=1 (Hawkins, Cramer, Truhlar pairwise GB)")
        console.print("  2. igb=2 (Modified GB model - default for constant pH)")
        console.print("  3. igb=5 (Modified GB model II)")
        console.print("  4. igb=7 (GBn model)")
        console.print("  5. igb=8 (GBn2 model)")

        choice = prompt_with_context(
            processor=self.processor,
            prompt="\nSelect GB model",
            choices=["1", "2", "3", "4", "5"],
            default="2",
            module="Topology Generator - CPIN",
            description="GB model",
            options_map={"1": "igb=1", "2": "igb=2 (recommended)", "3": "igb=5", "4": "igb=7", "5": "igb=8"}
        )

        igb_map = {"1": 1, "2": 2, "3": 5, "4": 7, "5": 8}
        igb = igb_map[choice]

        console.print(f"\n[bold]Step {step+1}: Internal Dielectric Constant[/bold]")
        intdiel_str = prompt_with_context(
            processor=self.processor,
            prompt="Internal dielectric for GB evaluation",
            default="1.0",
            module="Topology Generator - CPIN",
            description="Internal dielectric constant"
        )

        intdiel = float(intdiel_str)

        return igb, intdiel

    def _select_titratable_residues(self, step, prmtop_file, rst7_file):
        """
        Select titratable residues for cpin file.

        Uses direct topology scanning - no complex residue mapping needed!
        The residue numbers from the topology PDB are exactly what cpinutil expects.
        """
        console = self.processor.console
        workspace = self.get_workspace()

        console.print(f"\n[bold]Step {step}: Select Titratable Residues[/bold]")

        # DIRECT APPROACH: Scan the topology for titratable residues
        # This gives us the exact residue numbers that cpinutil expects
        titratable_residues = self._scan_topology_for_titratable_residues(prmtop_file, rst7_file)

        if not titratable_residues:
            console.print("[red]ERROR: No titratable residues found in topology.[/red]")
            console.print("[yellow]Make sure your topology contains titratable residues like AS4, GL4, HIS, LYS, etc.[/yellow]")
            return None

        # Display the found residues
        self._display_titratable_residues_table(titratable_residues)

        # Selection options
        console.print("\n[bold]How would you like to select residues?[/bold]")
        console.print("  1. Include all titratable residues (default)")
        console.print("  2. Include by residue type (e.g., ASP, GLU, HIS)")
        console.print("  3. Include by residue number")
        console.print("  4. Custom selection (exclude specific residues)")

        choices = ["1", "2", "3", "4"]

        choice = prompt_with_context(
            processor=self.processor,
            prompt="\nSelect option",
            choices=choices,
            default="1",
            module="Topology Generator - CPIN",
            description="Residue selection method"
        )

        if choice == "1":
            # All titratable residues
            selected = titratable_residues
        elif choice == "2":
            # By residue type
            selected = self._select_by_residue_type(titratable_residues)
        elif choice == "3":
            # By residue number
            selected = self._select_by_residue_number(titratable_residues)
        elif choice == "4":
            # Custom exclusion
            selected = self._select_custom(titratable_residues)
        else:
            selected = titratable_residues

        if not selected:
            console.print("[yellow]No residues selected.[/yellow]")
            return None

        # Confirm selection
        console.print(f"\n[green]Selected {len(selected)} titratable residues[/green]")
        
        if not confirm_with_context(
            processor=self.processor,
            prompt="Confirm selection?",
            default=True,
            module="Topology Generator - CPIN",
            description="Confirm residue selection"
        ):
            console.print("[yellow]Selection cancelled.[/yellow]")
            return None

        return selected

    def _extract_titratable_residues_from_analysis(self, protonation_data):
        """Extract titratable residues from protonation analysis"""
        titratable = []

        # Titratable residue types supported by AMBER constant pH MD
        # See: cpinutil.py --describe
        # AS4 (ASP), GL4 (GLU), HIS (HIP), CYS, TYR, LYS, PRN (heme propionate)
        # Note: ARG, N_TERM, C_TERM are NOT supported by AMBER CpHMD
        titratable_types = {
            'ASP', 'AS4',  # Aspartate (standard, constant pH)
            'GLU', 'GL4',  # Glutamate (standard, constant pH)
            'HIS', 'HID', 'HIE', 'HIP',  # Histidine variants
            'LYS',  # Lysine
            'CYS',  # Cysteine
            'TYR',  # Tyrosine
            'PRN'   # Propionate (for heme groups)
            # ARG is NOT included - not supported by AMBER CpHMD
        }

        for res_key, res_data in protonation_data.items():
            if isinstance(res_data, dict):
                # The actual keys are 'type' (not 'resname') and 'number' (not 'resnum')
                resname = res_data.get('type', '').upper()
                if resname in titratable_types:
                    # Determine protonation state from 'protonated' boolean
                    protonated = res_data.get('protonated', False)
                    state_at_ph = 'protonated' if protonated else 'deprotonated'

                    titratable.append({
                        'resnum': res_data.get('number'),
                        'chain': res_data.get('chain', 'A'),
                        'resname': resname,
                        'pka': res_data.get('pKa'),
                        'state_at_ph': state_at_ph,
                        'charge': res_data.get('charge', 0.0)
                    })

        return titratable

    def _display_titratable_residues_table(self, residues):
        """Display table of titratable residues"""
        console = self.processor.console

        console.print(f"\nDetected {len(residues)} titratable residues:")

        # Add a "Source" column only if any residue carries a pka_source tag
        # (set by the topology scanner when pb_titrate results are present).
        has_source = any(r.get('pka_source') for r in residues)

        table = Table()
        table.add_column("Res#", style="blue")
        table.add_column("Chain", style="blue")
        table.add_column("Residue", style="green")
        table.add_column("pKa", style="yellow")
        if has_source:
            table.add_column("pKa source", style="magenta")

        for res in residues:
            pka_str = f"{res['pka']:.2f}" if res.get('pka') is not None else "N/A"
            row = [
                str(res['resnum']),
                res['chain'],
                res['resname'],
                pka_str,
            ]
            if has_source:
                row.append(res.get('pka_source') or '—')
            table.add_row(*row)

        console.print(table)

    def _select_by_pka_range(self, residues):
        """Filter residues by pKa range"""
        console = self.processor.console
        
        console.print("\n[bold]Filter by pKa Range:[/bold]")
        
        min_pka_str = prompt_with_context(
            processor=self.processor,
            prompt="Minimum pKa to include (leave blank for no minimum)",
            default="",
            module="Topology Generator - CPIN",
            description="Minimum pKa"
        )
        
        max_pka_str = prompt_with_context(
            processor=self.processor,
            prompt="Maximum pKa to include (leave blank for no maximum)",
            default="",
            module="Topology Generator - CPIN",
            description="Maximum pKa"
        )
        
        min_pka = float(min_pka_str) if min_pka_str else None
        max_pka = float(max_pka_str) if max_pka_str else None
        
        selected = []
        for res in residues:
            pka = res.get('pka')
            if pka is None:
                continue
            if min_pka is not None and pka < min_pka:
                continue
            if max_pka is not None and pka > max_pka:
                continue
            selected.append(res)
        
        range_str = f"{min_pka if min_pka else 'any'} - {max_pka if max_pka else 'any'}"
        console.print(f"\n[green]Filtered selection: {len(selected)} residues (pKa {range_str})[/green]")
        
        return selected

    def _select_by_residue_type(self, residues):
        """Select residues by type"""
        console = self.processor.console
        
        console.print("\n[bold]Include by Residue Type[/bold]")
        console.print("Enter residue types to include (comma-separated)")
        console.print("Examples: ASP,GLU,HIS  or  HIS")
        
        types_str = prompt_with_context(
            processor=self.processor,
            prompt="Residue types to include",
            default="ASP,GLU,HIS,LYS,CYS,TYR",
            module="Topology Generator - CPIN",
            description="Residue types"
        )
        
        selected_types = {t.strip().upper() for t in types_str.split(',')}
        selected = [r for r in residues if r['resname'] in selected_types]
        
        console.print(f"\n[green]Selected {len(selected)} residues of types: {', '.join(selected_types)}[/green]")
        return selected

    def _select_by_residue_number(self, residues):
        """Select residues by number"""
        console = self.processor.console
        
        console.print("\n[bold]Include by Residue Number[/bold]")
        console.print("Enter residue numbers (comma-separated or ranges)")
        console.print("Examples: 25,42,56  or  25-50,78-100")
        
        nums_str = prompt_with_context(
            processor=self.processor,
            prompt="Residue numbers to include",
            default="",
            module="Topology Generator - CPIN",
            description="Residue numbers"
        )
        
        # Parse ranges and numbers
        selected_nums = set()
        for part in nums_str.split(','):
            part = part.strip()
            if '-' in part:
                start, end = part.split('-')
                selected_nums.update(range(int(start), int(end) + 1))
            else:
                selected_nums.add(int(part))
        
        selected = [r for r in residues if r['resnum'] in selected_nums]
        
        console.print(f"\n[green]Selected {len(selected)} residues[/green]")
        return selected

    def _select_custom(self, residues):
        """Custom selection with exclusions"""
        console = self.processor.console
        
        console.print("\n[bold]Custom Selection[/bold]")
        console.print("Enter residue numbers to EXCLUDE (comma-separated)")
        
        exclude_str = prompt_with_context(
            processor=self.processor,
            prompt="Residue numbers to exclude",
            default="",
            module="Topology Generator - CPIN",
            description="Exclude residues"
        )
        
        if not exclude_str:
            return residues
        
        exclude_nums = {int(n.strip()) for n in exclude_str.split(',')}
        selected = [r for r in residues if r['resnum'] not in exclude_nums]
        
        console.print(f"\n[green]Selected {len(selected)} residues (excluded {len(exclude_nums)})[/green]")
        return selected

    def _set_initial_protonation_states(self, selected_residues, step):
        """Set initial protonation states for constant pH MD."""
        console = self.processor.console
        workspace = self.get_workspace()

        console.print(f"\n[bold]Step {step}: Set Initial Protonation States[/bold]")

        # Standard pKa values from cpinutil
        pka_values = {
            'AS4': 4.0, 'GL4': 4.4, 'HIP': 6.6,
            'LYS': 10.4, 'CYS': 8.5, 'TYR': 9.6, 'PRN': 4.8
        }

        # Option 3 (Use titrate recommendations) is only offered when the
        # PB Titrate workflow has run and persisted a recommendations dict.
        titrate_recs = workspace.get("titrate_recommendations") or {}
        has_titrate = bool(titrate_recs)

        console.print("\nHow would you like to set initial protonation states?")
        console.print("  1. Use cpinutil defaults (reference states)")
        console.print("  2. Set states based on target pH (using standard pKa values)")
        if has_titrate:
            console.print(f"  3. Use PB Titrate recommendations "
                            f"({len(titrate_recs)} sites refined)")

        choices = ["1", "2", "3"] if has_titrate else ["1", "2"]
        options_map = {"1": "cpinutil defaults", "2": "pH-based"}
        if has_titrate:
            options_map["3"] = "PB titrate recommendations"
            default_choice = "3"
        else:
            default_choice = "1"

        choice = prompt_with_context(
            processor=self.processor,
            prompt="\nSelect option",
            choices=choices,
            default=default_choice,
            module="Topology Generator - CPIN",
            description="Initial state method",
            options_map=options_map,
        )

        if choice == "1":
            # Use cpinutil defaults - no -states flag
            console.print("\n[grey50]Using cpinutil default reference states (no -states flag).[/grey50]")
            console.print("")
            console.print("cpinutil default states:")
            console.print("  AS4/GL4/PRN: State 0 (deprotonated, COO⁻)")
            console.print("  LYS/CYS/TYR: State 0 (protonated)")
            console.print("  HIP: State 0 (doubly protonated, HIP⁺)")
            return None

        if choice == "3":
            # Use PB Titrate recommendations from the workspace.
            # titrate_recs schema: {(resname, resnum) -> {state_id, state_name,
            #                       prot_count, pka_corr, net_charge}}.
            # Selected residues that have no titrate recommendation fall back
            # to cpinutil's default state for their type (state 0 for all
            # supported resnames per cpinutil --describe).
            console.print(
                "\n[blue]Setting initial states from PB Titrate "
                "recommendations[/blue]")
            initial_states = {}
            n_recommended = 0
            n_default = 0
            state_summary = {}
            for res in selected_residues:
                resname = res['resname']
                resnum = res['resnum']
                rec = titrate_recs.get((resname, resnum))
                if rec is None:
                    # No recommendation — fall back to cpinutil default
                    state = 0
                    n_default += 1
                else:
                    state = rec["state_id"]
                    n_recommended += 1
                initial_states[resnum] = state
                state_summary.setdefault(resname, {}).setdefault(state, 0)
                state_summary[resname][state] += 1
            console.print(f"  PB-recommended: {n_recommended}  |  "
                            f"cpinutil-default fallback: {n_default}")
            console.print("")
            for resname in sorted(state_summary.keys()):
                for state, count in sorted(state_summary[resname].items()):
                    description = self._get_cpin_state_description(resname, state)
                    console.print(f"  {resname} state {state} = {description} "
                                    f"({count} residues)")
            return initial_states

        # Option 2: Set states based on target pH
        ph_str = prompt_with_context(
            processor=self.processor,
            prompt="Target pH",
            default="7.0",
            module="Topology Generator - CPIN",
            description="Target pH for initial states"
        )
        ph = float(ph_str)

        console.print(f"\n[blue]Setting initial states for pH {ph} based on standard pKa values:[/blue]")

        # Calculate states based on pH vs pKa
        initial_states = {}
        state_summary = {}  # Track states by residue type

        for res in selected_residues:
            resname = res['resname']
            pka = pka_values.get(resname, 7.0)

            # Determine state based on pH vs pKa
            if resname in ['AS4', 'GL4', 'PRN']:
                # Carboxylic acids: deprotonated (state 0) if pH > pKa
                if ph > pka:
                    state = 0  # Deprotonated (COO⁻)
                else:
                    state = 1  # Protonated (COOH)
            elif resname == 'HIP':
                # Histidine: doubly protonated (state 0) if pH < pKa, else singly protonated
                if ph < pka:
                    state = 0  # Doubly protonated (HIP⁺)
                else:
                    state = 2  # Singly protonated (HIE, neutral)
            elif resname == 'LYS':
                # Lysine: protonated (state 0) if pH < pKa
                if ph < pka:
                    state = 0  # Protonated (NH₃⁺)
                else:
                    state = 1  # Deprotonated (NH₂)
            elif resname in ['CYS', 'TYR']:
                # Cysteine/Tyrosine: protonated (state 0) if pH < pKa
                if ph < pka:
                    state = 0  # Protonated
                else:
                    state = 1  # Deprotonated
            else:
                state = 0

            initial_states[res['resnum']] = state

            # Track for summary
            if resname not in state_summary:
                state_summary[resname] = {}
            if state not in state_summary[resname]:
                state_summary[resname][state] = 0
            state_summary[resname][state] += 1

        # Display state summary
        console.print("")
        for resname in sorted(state_summary.keys()):
            pka = pka_values.get(resname, "?")
            for state, count in sorted(state_summary[resname].items()):
                description = self._get_cpin_state_description(resname, state)
                relation = ">" if ph > pka else "<" if ph < pka else "="
                console.print(f"  {resname} (pKa {pka}): pH {ph} {relation} pKa → State {state} = {description} ({count} residues)")

        return initial_states

    def _get_cpin_state_description(self, resname, state):
        """Get human-readable description of a cpin state"""
        state_descriptions = {
            'ASP': {
                0: "COO⁻ (deprotonated, -1)",
                1: "COOH (protonated, neutral)",
            },
            'AS4': {
                0: "COO⁻ (deprotonated, -1)",
                1: "COOH (protonated, neutral)",
                2: "COOH (protonated, neutral)",
                3: "COOH (protonated, neutral)",
                4: "COOH (protonated, neutral)",
            },
            'GLU': {
                0: "COO⁻ (deprotonated, -1)",
                1: "COOH (protonated, neutral)",
            },
            'GL4': {
                0: "COO⁻ (deprotonated, -1)",
                1: "COOH (protonated, neutral)",
                2: "COOH (protonated, neutral)",
                3: "COOH (protonated, neutral)",
                4: "COOH (protonated, neutral)",
            },
            'HIS': {
                0: "HIP (doubly protonated, +1)",
                1: "HID (proton on ND1, neutral)",
                2: "HIE (proton on NE2, neutral)",
            },
            'HIP': {
                0: "HIP (doubly protonated, +1)",
                1: "HID (proton on ND1, neutral)",
                2: "HIE (proton on NE2, neutral)",
            },
            'LYS': {
                0: "NH₃⁺ (protonated, +1)",
                1: "NH₂ (deprotonated, neutral)",
            },
            'CYS': {
                0: "SH (protonated, neutral)",
                1: "S⁻ (deprotonated, -1)",
            },
            'TYR': {
                0: "OH (protonated, neutral)",
                1: "O⁻ (deprotonated, -1)",
            },
            'PRN': {
                0: "COO⁻ (deprotonated, -1)",
                1: "COOH (protonated, neutral)",
                2: "COOH (protonated, neutral)",
                3: "COOH (protonated, neutral)",
                4: "COOH (protonated, neutral)",
            },
        }

        return state_descriptions.get(resname, {}).get(state, f"State {state}")

    def _map_protonation_state_to_cpin(self, residue):
        """
        Map protonation state to cpin state index.

        Based on cpinutil state definitions (from cpinutil.py --describe):
        - AS4: pKa=4.0, State 0 = deprotonated (-1), States 1-4 = protonated (neutral)
        - GL4: pKa=4.4, State 0 = deprotonated (-1), States 1-4 = protonated (neutral)
        - HIP: pKa=6.6, State 0 = doubly protonated (+1), State 1 = HID, State 2 = HIE
        - LYS: pKa=10.4, State 0 = protonated (+1), State 1 = deprotonated (neutral)
        - CYS: pKa=8.5, State 0 = protonated (neutral), State 1 = deprotonated (-1)
        - TYR: pKa=9.6, State 0 = protonated (neutral), State 1 = deprotonated (-1)
        - PRN: pKa=4.8, State 0 = deprotonated (-1), States 1-4 = protonated (neutral)
        """
        resname = residue['resname']
        state_at_ph = residue.get('state_at_ph', 'unknown').lower()
        charge = residue.get('charge', 0.0)

        # ASP and GLU: State 0 = deprotonated, State 1 = protonated
        if resname in ['ASP', 'AS4']:
            if 'deprotonated' in state_at_ph or charge < -0.5:
                return 0  # Deprotonated (COO-)
            else:
                return 1  # Protonated (COOH) - default to first tautomer

        elif resname in ['GLU', 'GL4']:
            if 'deprotonated' in state_at_ph or charge < -0.5:
                return 0  # Deprotonated (COO-)
            else:
                return 1  # Protonated (COOH) - default to first tautomer

        # HIS: State 0 = HIP (+1), State 1 = HID (ND1-H), State 2 = HIE (NE2-H)
        elif resname in ['HIS', 'HIP', 'HID', 'HIE']:
            if 'doubly protonated' in state_at_ph or charge > 0.5:
                return 0  # HIP (both nitrogens protonated, +1)
            elif 'hid' in state_at_ph.lower() or 'nd1' in state_at_ph.lower():
                return 1  # HID (proton on ND1)
            elif 'hie' in state_at_ph.lower() or 'ne2' in state_at_ph.lower():
                return 2  # HIE (proton on NE2)
            else:
                # Default to HIE if neutral
                return 2

        # LYS: State 0 = protonated (+1), State 1 = deprotonated (neutral)
        elif resname == 'LYS':
            if 'protonated' in state_at_ph or charge > 0.5:
                return 0  # Protonated (NH3+)
            else:
                return 1  # Deprotonated (NH2)

        # CYS: State 0 = protonated (SH), State 1 = deprotonated (S-)
        elif resname == 'CYS':
            if 'deprotonated' in state_at_ph or charge < -0.5:
                return 1  # Deprotonated (S-)
            else:
                return 0  # Protonated (SH)

        # TYR: State 0 = protonated (OH), State 1 = deprotonated (O-)
        elif resname == 'TYR':
            if 'deprotonated' in state_at_ph or charge < -0.5:
                return 1  # Deprotonated (O-)
            else:
                return 0  # Protonated (OH)

        # PRN (propionate): State 0 = deprotonated (COO-), States 1-4 = protonated (COOH)
        elif resname == 'PRN':
            if 'deprotonated' in state_at_ph or charge < -0.5:
                return 0  # Deprotonated (COO-)
            else:
                return 1  # Protonated (COOH) - default to first tautomer

        # Default fallback
        else:
            return 0

    def _parse_pdb_residues(self, pdb_file):
        """
        Parse PDB file to extract residue information directly from ATOM/HETATM records.

        This avoids BioPython which can lose atoms or mishandle non-standard residues
        like AS4, GL4, HIP, etc.

        Args:
            pdb_file: Path to PDB file

        Returns:
            List of dicts with 'chain', 'resnum', 'resname', 'ca_coords' for each unique residue
        """
        residues = []
        seen_residues = set()
        residue_coords = {}  # Store CA coordinates for each residue
        residue_atom_names = {}  # key -> set of atom names

        try:
            with open(pdb_file, 'r') as f:
                for line in f:
                    # Only process ATOM and HETATM records
                    if not (line.startswith('ATOM') or line.startswith('HETATM')):
                        continue

                    # PDB format column positions (1-indexed in spec, 0-indexed in Python)
                    # ATOM/HETATM record format:
                    # Columns:  1-6   7-11  13-16  17    18-20  22    23-26    27    31-38  39-46  47-54
                    # Fields:   ATOM  serial name   altLoc resName chain resSeq  iCode x      y      z

                    try:
                        atom_name = line[12:16].strip()  # Atom name (columns 13-16)
                        resname = line[17:20].strip()  # Residue name (columns 18-20)
                        chain = line[21:22].strip()     # Chain ID (column 22)
                        resnum_str = line[22:26].strip()  # Residue number (columns 23-26)
                        icode = line[26:27].strip() if len(line) > 26 else ''  # Insertion code (column 27)

                        # Handle chain - if empty, use space as default
                        if not chain:
                            chain = ' '

                        # Parse residue number (may have insertion code)
                        # Remove any non-numeric characters for the base number
                        resnum_clean = ''.join(c for c in resnum_str if c.isdigit() or c == '-')
                        if not resnum_clean:
                            continue  # Skip if no numeric part

                        resnum = int(resnum_clean)

                        # Include insertion code in the residue number if present
                        # Store as string to preserve insertion codes like "100A"
                        if icode:
                            resnum_with_icode = f"{resnum}{icode}"
                        else:
                            resnum_with_icode = resnum

                        # Create unique key
                        key = (chain, resnum_with_icode, resname)

                        # Extract coordinates for CA atoms (or first heavy atom for non-protein)
                        if len(line) >= 54:
                            try:
                                x = float(line[30:38].strip())
                                y = float(line[38:46].strip())
                                z = float(line[46:54].strip())
                                coords = (x, y, z)

                                # Prefer CA for proteins, or first atom for others
                                if atom_name == 'CA' or key not in residue_coords:
                                    residue_coords[key] = coords
                            except ValueError:
                                pass

                        # Only add if not seen before (avoid duplicates from multiple atoms)
                        if key not in seen_residues:
                            residues.append({
                                'chain': chain,
                                'resnum': resnum_with_icode,  # Store with insertion code
                                'resname': resname,
                                'ca_coords': None,  # Will be filled in below
                                'atom_names': None,  # Will be filled in below
                            })
                            seen_residues.add(key)
                            residue_atom_names[key] = set()
                        residue_atom_names[key].add(atom_name)

                    except (ValueError, IndexError):
                        # Skip malformed lines
                        continue

            # Add coordinates and atom-name sets to residue records
            for res in residues:
                key = (res['chain'], res['resnum'], res['resname'])
                res['ca_coords'] = residue_coords.get(key)
                res['atom_names'] = residue_atom_names.get(key, set())

        except Exception as e:
            logger.error(f"Error parsing PDB file {pdb_file}: {e}")
            return []

        return residues

    def _scan_topology_for_titratable_residues(self, prmtop_file, rst7_file):
        """
        Directly scan topology for titratable residues - no mapping needed!

        This is the simplest and most reliable approach: generate a PDB from the
        final topology and scan it directly for titratable residue types. The residue
        numbers in the PDB are exactly what cpinutil expects.

        Args:
            prmtop_file: Path to topology file (.prmtop)
            rst7_file: Path to coordinate file (.rst7)

        Returns:
            List of dicts with 'resnum', 'resname', 'chain' for titratable residues,
            or None if failed
        """
        console = self.processor.console
        workspace = self.get_workspace()

        # Titratable residue types supported by AMBER constant pH MD
        # IMPORTANT: Only specific residue names are titratable in CpHMD!
        # - HIP (not HIE/HID) - doubly protonated histidine
        # - AS4 (not ASP) - titratable aspartate
        # - GL4 (not GLU) - titratable glutamate
        # - LYS, CYS, TYR, PRN are titratable with their standard names
        # Standard pKa values from cpinutil.py --describe
        titratable_info = {
            'AS4': 4.0,   # Titratable aspartate (NOT ASP!)
            'GL4': 4.4,   # Titratable glutamate (NOT GLU!)
            'HIP': 6.6,   # Titratable histidine (NOT HIE/HID!)
            'LYS': 10.4,  # Lysine
            'CYS': 8.5,   # Cysteine (free thiol)
            'TYR': 9.6,   # Tyrosine
            'PRN': 4.8,   # Propionate (for heme groups)
        }
        titratable_types = set(titratable_info.keys())

        console.print("\n[grey50]Scanning topology for titratable residues...[/grey50]")

        # Generate PDB from topology using cpptraj
        import tempfile
        output_dir = workspace.get("output_dir", ".")

        topo_pdb = tempfile.NamedTemporaryFile(
            mode='w', suffix='_topology.pdb', delete=False, dir=output_dir
        )
        topo_pdb_path = topo_pdb.name
        topo_pdb.close()

        try:
            # Create cpptraj input
            cpptraj_input = f"""parm {prmtop_file}
trajin {rst7_file}
trajout {topo_pdb_path} pdb
go
quit
"""
            cpptraj_input_file = tempfile.NamedTemporaryFile(
                mode='w', suffix='.cpptraj', delete=False, dir=output_dir
            )
            cpptraj_input_file.write(cpptraj_input)
            cpptraj_input_file.close()

            # Run cpptraj
            result = subprocess.run(
                ['cpptraj', '-i', cpptraj_input_file.name],
                capture_output=True,
                text=True,
                timeout=60
            )

            os.unlink(cpptraj_input_file.name)

            if result.returncode != 0:
                console.print(f"[red]ERROR: cpptraj failed: {result.stderr}[/red]")
                return None

            if not os.path.exists(topo_pdb_path):
                console.print("[red]ERROR: cpptraj did not generate PDB file.[/red]")
                return None

            # Parse the topology PDB
            residues = self._parse_pdb_residues(topo_pdb_path)

            # Filter for titratable residues only (exclude water/ions)
            water_ion_names = {'WAT', 'HOH', 'Na+', 'Cl-', 'K+', 'Mg2', 'NA', 'CL', 'MG', 'CA', 'ZN'}

            # Terminal-residue detection. Amber's constph framework has no
            # multi-state libraries for N/C-terminal forms (constph.lib only
            # ships internal AS2/AS4/GL4; aminoct10.lib only has single-state
            # CLYS/CASP), so cpinutil cannot titrate them — including them
            # silently corrupts terminal-specific backbone charges. We detect
            # terminals here via the same H1/H2/H3 (N-term) / OXT (C-term)
            # atom markers used by pb_titrate.sites._is_terminal. PRN
            # (heme propionate) has no chain termini and is exempt.
            NTERM_MARKERS = {'H1', 'H2', 'H3'}
            CTERM_MARKERS = {'OXT'}

            # Optional cross-check: pb_titrate may have already flagged these
            # in its step 2 (workspace key 'pb_titrate_terminal_excluded').
            pbt_excluded = workspace.get("pb_titrate_terminal_excluded") or []
            pbt_excluded_keys = {
                (e['resname'], int(e['resnum'])) for e in pbt_excluded
                if isinstance(e, dict)
            }

            # pb_titrate per-site pKa map (from step 3). When present, those
            # values override the cpinutil-header defaults in the displayed
            # table so users see the actual computed pKa for their structure.
            pbt_pka_map = workspace.get("pb_titrate_pka") or {}
            # Workspaces sometimes round-trip dict keys via JSON, which kills
            # tuple keys. Normalize: accept either tuple or "RESNAME:RESNUM".
            normalized_pka: Dict[Tuple[str, int], float] = {}
            for k, v in (pbt_pka_map.items() if hasattr(pbt_pka_map, 'items') else []):
                if v is None:
                    continue
                if isinstance(k, tuple) and len(k) == 2:
                    normalized_pka[(k[0], int(k[1]))] = float(v)
                elif isinstance(k, str) and ':' in k:
                    rn, rnum = k.split(':', 1)
                    try:
                        normalized_pka[(rn.strip(), int(rnum.strip()))] = float(v)
                    except ValueError:
                        continue
            pb_source = "PB Titrate" if normalized_pka else "cpinutil --describe"

            titratable_residues = []
            terminal_excluded = []
            for res in residues:
                if res['resname'] in water_ion_names:
                    continue
                if res['resname'] not in titratable_types:
                    continue
                # Skip terminals (PRN exempt — no backbone termini)
                atom_names = res.get('atom_names') or set()
                terminal_kind = None
                if res['resname'] != 'PRN':
                    if atom_names & NTERM_MARKERS:
                        terminal_kind = 'N-terminal'
                    elif atom_names & CTERM_MARKERS:
                        terminal_kind = 'C-terminal'
                # Also honor pb_titrate's flagging
                try:
                    rkey = (res['resname'], int(res['resnum']))
                except (TypeError, ValueError):
                    rkey = None
                if rkey is not None and rkey in pbt_excluded_keys:
                    terminal_kind = terminal_kind or 'terminal (per pb_titrate)'
                if terminal_kind is not None:
                    terminal_excluded.append({
                        'resnum': res['resnum'],
                        'resname': res['resname'],
                        'chain': res['chain'],
                        'kind': terminal_kind,
                        'markers': sorted(atom_names
                                          & (NTERM_MARKERS | CTERM_MARKERS)),
                    })
                    continue
                # pKa source: PB Titrate per-site if available, else cpinutil default
                pka = (normalized_pka.get(rkey)
                       if rkey is not None
                       else None)
                if pka is None:
                    pka = titratable_info[res['resname']]
                    src = "cpinutil --describe"
                else:
                    src = "PB Titrate"
                titratable_residues.append({
                    'resnum': res['resnum'],
                    'chain': res['chain'],
                    'resname': res['resname'],
                    'pka': pka,
                    'pka_source': src,
                })

            # Clean up
            os.unlink(topo_pdb_path)

            console.print(f"[green]✓ Found {len(titratable_residues)} titratable residues in topology[/green]")

            # Show summary by type
            from collections import Counter
            type_counts = Counter(r['resname'] for r in titratable_residues)
            type_summary = ", ".join(f"{count} {restype}" for restype, count in sorted(type_counts.items()))
            if type_summary:
                console.print(f"[grey50]  Types: {type_summary}[/grey50]")

            # Report terminal exclusions
            if terminal_excluded:
                console.print(
                    f"\n[yellow]Excluded {len(terminal_excluded)} terminal "
                    f"residue(s) — Amber's constph has no multi-state "
                    f"libraries for chain termini:[/yellow]")
                for e in terminal_excluded:
                    markers = ', '.join(e['markers']) or '(no markers)'
                    console.print(
                        f"  [yellow]{e['resname']}-{e['resnum']}[/yellow] "
                        f"({e['kind']}; detected via {markers})")

            # Report pKa source so users know what they're looking at
            n_pbt = sum(1 for r in titratable_residues
                         if r.get('pka_source') == "PB Titrate")
            if n_pbt > 0:
                console.print(
                    f"[grey50]  pKa values: {n_pbt} from PB Titrate (computed "
                    f"for this structure), "
                    f"{len(titratable_residues) - n_pbt} from cpinutil "
                    f"--describe (residue-type defaults)[/grey50]")
            else:
                console.print(
                    f"[grey50]  pKa values: all from cpinutil --describe "
                    f"(residue-type defaults; run PB Titrate to get "
                    f"per-site values for this structure)[/grey50]")

            return titratable_residues

        except subprocess.TimeoutExpired:
            console.print("[red]ERROR: cpptraj timed out[/red]")
            return None
        except FileNotFoundError:
            console.print("[red]ERROR: cpptraj not found in PATH[/red]")
            return None
        except Exception as e:
            console.print(f"[red]ERROR scanning topology: {e}[/red]")
            logger.exception("Topology scan error:")
            return None

    def _create_residue_number_mapping(self, tleap_info):
        """
        Create mapping from original PDB residue numbers to tLEaP-renumbered residues.

        Uses cpptraj to generate a PDB from prmtop/rst7, then compares residue
        numbering between original and tLEaP structures.

        Args:
            tleap_info: Parsed tLEaP file info containing input PDB and output files

        Returns:
            Dictionary mapping {(chain, original_resnum): tleap_resnum} or None if failed
        """
        console = self.processor.console
        workspace = self.get_workspace()

        console.print("\n[grey50]Creating residue number mapping...[/grey50]")

        # Get original PDB file from tLEaP input
        if not tleap_info['pdb_files']:
            console.print("[red]ERROR: No input PDB found in tLEaP file.[/red]")
            return None

        original_pdb = tleap_info['pdb_files'][0]  # Use first PDB
        if not os.path.exists(original_pdb):
            console.print(f"[red]ERROR: Original PDB file not found: {original_pdb}[/red]")
            return None

        # Get prmtop and rst7 from tLEaP info
        prmtop_file = tleap_info['expected_prmtop']
        rst7_file = tleap_info['expected_rst7']

        if not prmtop_file or not os.path.exists(prmtop_file):
            console.print(f"[red]ERROR: Topology file not found: {prmtop_file}[/red]")
            return None

        if not rst7_file or not os.path.exists(rst7_file):
            console.print(f"[red]ERROR: Coordinate file not found: {rst7_file}[/red]")
            return None

        # Generate PDB from prmtop/rst7 using cpptraj
        import tempfile
        tleap_pdb = tempfile.NamedTemporaryFile(mode='w', suffix='_tleap.pdb', delete=False, dir=workspace.get("output_dir", "."))
        tleap_pdb_path = tleap_pdb.name
        tleap_pdb.close()

        try:
            # Create cpptraj input
            cpptraj_input = f"""parm {prmtop_file}
trajin {rst7_file}
trajout {tleap_pdb_path} pdb
go
quit
"""

            cpptraj_input_file = tempfile.NamedTemporaryFile(mode='w', suffix='.cpptraj', delete=False, dir=workspace.get("output_dir", "."))
            cpptraj_input_file.write(cpptraj_input)
            cpptraj_input_file.close()

            # Run cpptraj
            console.print(f"[grey50]Running cpptraj to generate PDB from topology...[/grey50]")
            result = subprocess.run(
                ['cpptraj', '-i', cpptraj_input_file.name],
                capture_output=True,
                text=True,
                timeout=60
            )

            # Clean up cpptraj input file
            os.unlink(cpptraj_input_file.name)

            if result.returncode != 0:
                console.print(f"[red]ERROR: cpptraj failed:[/red]")
                console.print(result.stderr)
                return None

            if not os.path.exists(tleap_pdb_path):
                console.print(f"[red]ERROR: cpptraj did not generate PDB file.[/red]")
                return None

            # Parse both PDB files directly (avoid BioPython to handle non-standard residues)
            console.print(f"[grey50]Parsing original PDB: {os.path.basename(original_pdb)}[/grey50]")
            original_residues = self._parse_pdb_residues(original_pdb)

            console.print(f"[grey50]Parsing tLEaP-generated PDB: {os.path.basename(tleap_pdb_path)}[/grey50]")
            tleap_residues = self._parse_pdb_residues(tleap_pdb_path)

            if not original_residues:
                console.print(f"[red]ERROR: No residues found in original PDB[/red]")
                return None

            if not tleap_residues:
                console.print(f"[red]ERROR: No residues found in tLEaP PDB[/red]")
                return None

            # Create mapping using coordinate-based matching
            # tLEaP/ParmEd may reorder residues, so we match by CA atom coordinates
            num_original = len(original_residues)
            num_tleap = len(tleap_residues)

            # Only consider non-water tLEaP residues for mapping
            water_ion_names = {'WAT', 'HOH', 'Na+', 'Cl-', 'K+', 'Mg2', 'NA', 'CL', 'MG', 'CA', 'ZN'}
            tleap_protein_residues = [r for r in tleap_residues if r['resname'] not in water_ion_names]

            console.print(f"[grey50]Original PDB: {num_original} residues[/grey50]")
            console.print(f"[grey50]tLEaP PDB: {len(tleap_protein_residues)} non-water residues (+ {num_tleap - len(tleap_protein_residues)} water/ions)[/grey50]")

            # Build coordinate map for tLEaP residues (using CA or first atom coords)
            tleap_coord_map = {}
            for res in tleap_protein_residues:
                if res.get('ca_coords'):
                    # Round coordinates to handle floating point differences
                    coords = tuple(round(c, 2) for c in res['ca_coords'])
                    tleap_coord_map[coords] = res['resnum']

            console.print(f"[grey50]Built coordinate map for {len(tleap_coord_map)} tLEaP residues[/grey50]")

            # Match original residues to tLEaP residues by coordinates
            mapping = {}
            matched = 0
            unmatched = 0

            for orig_res in original_residues:
                if orig_res.get('ca_coords'):
                    coords = tuple(round(c, 2) for c in orig_res['ca_coords'])
                    if coords in tleap_coord_map:
                        key = (orig_res['chain'], orig_res['resnum'])
                        mapping[key] = tleap_coord_map[coords]
                        matched += 1
                    else:
                        # Try with less precision
                        coords_1 = tuple(round(c, 1) for c in orig_res['ca_coords'])
                        found = False
                        for tleap_coords, tleap_resnum in tleap_coord_map.items():
                            tleap_coords_1 = tuple(round(c, 1) for c in tleap_coords)
                            if coords_1 == tleap_coords_1:
                                key = (orig_res['chain'], orig_res['resnum'])
                                mapping[key] = tleap_resnum
                                matched += 1
                                found = True
                                break
                        if not found:
                            unmatched += 1
                else:
                    unmatched += 1

            console.print(f"[green]✓ Mapped {matched} residues by coordinates[/green]")
            if unmatched > 0:
                console.print(f"[yellow]Warning: {unmatched} residues could not be mapped[/yellow]")

            # Show sample mappings
            if mapping:
                from collections import defaultdict
                by_chain = defaultdict(list)
                for (chain, orig_num), tleap_num in mapping.items():
                    by_chain[chain].append((orig_num, tleap_num))

                console.print("[grey50]Sample mappings:[/grey50]")
                for chain in sorted(by_chain.keys())[:3]:  # Show first 3 chains
                    chain_maps = by_chain[chain]
                    if chain_maps:
                        first = chain_maps[0]
                        console.print(f"[grey50]  {chain}:{first[0]} → {first[1]}[/grey50]")

            # Clean up temporary PDB file
            os.unlink(tleap_pdb_path)

            return mapping

        except subprocess.TimeoutExpired:
            console.print("[red]ERROR: cpptraj timed out[/red]")
            return None
        except FileNotFoundError:
            console.print("[red]ERROR: cpptraj not found in PATH[/red]")
            console.print("Make sure AmberTools is installed and in your PATH")
            return None
        except Exception as e:
            console.print(f"[red]ERROR creating residue mapping: {e}[/red]")
            logger.exception("Residue mapping error:")
            return None

    def _run_cpinutil(self, prmtop_file, tleap_input_file, tleap_info, sim_type, igb, intdiel, selected_residues, initial_states, step):
        """Run cpinutil to generate cpin file"""
        console = self.processor.console
        workspace = self.get_workspace()

        console.print(f"\n[bold]Step {step}: Generate CPIN File[/bold]")
        
        # Get output filename
        default_cpin = os.path.splitext(os.path.basename(prmtop_file))[0] + ".cpin"
        output_cpin = prompt_with_context(
            processor=self.processor,
            prompt="Output filename",
            default=default_cpin,
            module="Topology Generator - CPIN",
            description="Output cpin filename"
        )
        
        # Get system name
        system_name = prompt_with_context(
            processor=self.processor,
            prompt="System name (for documentation)",
            default="system",
            module="Topology Generator - CPIN",
            description="System name"
        )
        
        # Make output path absolute
        output_dir = workspace.get("output_dir", ".")
        output_cpin_path = os.path.join(output_dir, output_cpin)

        # NOTE: No mapping needed! Residues were scanned directly from the topology,
        # so their residue numbers are exactly what cpinutil expects.
        console.print(f"[grey50]Using {len(selected_residues)} residues with topology numbering[/grey50]")

        # Build cpinutil command
        cmd = ["cpinutil.py", "-p", prmtop_file]
        
        if sim_type == 'implicit':
            cmd.extend(["-igb", str(igb)])
            cmd.extend(["-intdiel", str(intdiel)])
        else:
            # Explicit solvent - generate modified prmtop
            modified_prmtop = os.path.splitext(prmtop_file)[0] + "_cpin.prmtop"
            cmd.extend(["-op", modified_prmtop])
        
        cmd.extend(["-o", output_cpin_path])
        cmd.extend(["-system", system_name])
        
        # Add residue selection if not all
        if selected_residues:
            resnums = [str(r['resnum']) for r in selected_residues]
            cmd.extend(["-resnums"] + resnums)
        
        # Add initial states if provided
        if initial_states:
            states = [str(initial_states.get(r['resnum'], 0)) for r in selected_residues]
            cmd.extend(["-states"] + states)
        
        # Display command
        console.print("\n[bold]Generating cpin file...[/bold]")
        console.print(f"[grey50]Running: {' '.join(cmd)}[/grey50]\n")
        
        # Run cpinutil
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            
            if result.returncode != 0:
                console.print(f"[red]Error running cpinutil:[/red]")
                console.print(result.stderr)
                return False
            
            console.print("[green]✓ Generated:[/green] " + output_cpin)

            if sim_type == 'explicit':
                console.print("[green]✓ Generated modified prmtop:[/green] " + os.path.basename(modified_prmtop))

            # Save command to a script for reproducibility
            script_name = os.path.splitext(output_cpin)[0] + "_generate.sh"
            script_path = os.path.join(output_dir, script_name)
            try:
                with open(script_path, 'w') as f:
                    f.write("#!/bin/bash\n")
                    f.write("# Auto-generated script to regenerate cpin file\n")
                    f.write(f"# Generated by ProPrep on {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                    f.write("#\n")
                    f.write("# Usage: bash " + script_name + "\n")
                    f.write("# Edit residue numbers or states below as needed\n\n")

                    # Write the command with line continuation for readability
                    # Format: one flag per line for long lists (resnums, states)
                    formatted_cmd = []
                    i = 0
                    while i < len(cmd):
                        part = cmd[i]
                        # Handle flags that take multiple arguments (resnums, states)
                        if part in ['-resnums', '-states'] and i + 1 < len(cmd):
                            # Start the flag on its own line
                            formatted_cmd.append(f"  {part}")
                            i += 1
                            # Collect all numeric arguments following this flag
                            args = []
                            while i < len(cmd) and not cmd[i].startswith('-'):
                                args.append(cmd[i])
                                i += 1
                            # Format args: 10 per line for readability
                            for j in range(0, len(args), 10):
                                chunk = ' '.join(args[j:j+10])
                                if j > 0:
                                    formatted_cmd.append(f"    {chunk}")
                                else:
                                    formatted_cmd.append(f" {chunk}")
                        else:
                            # Regular flag or value
                            formatted_cmd.append(f"  {part}")
                            i += 1

                    # Write with backslash line continuations
                    f.write("cpinutil.py \\\n")
                    f.write(" \\\n".join(formatted_cmd[1:]))  # Skip 'cpinutil.py' since we wrote it above
                    f.write('\n')

                # Make script executable
                import stat
                os.chmod(script_path, os.stat(script_path).st_mode | stat.S_IEXEC)

                console.print(f"[green]✓ Saved command script:[/green] {script_name}")
            except Exception as e:
                console.print(f"[yellow]⚠ Could not save script file: {e}[/yellow]")

            # Display summary
            console.print(f"\n[bold]Summary:[/bold]")
            console.print(f"  - {len(selected_residues)} titratable residues included")
            if initial_states:
                console.print(f"  - Initial states set from protonation analysis")
            if sim_type == 'implicit':
                console.print(f"  - Compatible with implicit solvent (igb={igb})")
            else:
                console.print(f"  - Compatible with explicit solvent")
            
            console.print(f"\n[bold blue]Next steps:[/bold blue]")
            console.print(
                f"  Run MD via [bold]MD Manager → Setup and configure simulations[/bold], "
                f"then pick a predefined protocol.")
            console.print(
                f"  ProPrep detects this CPIN and offers to enable constant pH MD on the "
                f"production")
            console.print(
                f"  steps for you — it sets icnstph/solvph and the "
                f"[bold]-cpin -cpout -cprestrt[/bold] flags automatically"
                + (" (and uses the modified prmtop)." if sim_type == 'explicit' else "."))
            console.print(
                f"  [grey50]Manual route: run {output_cpin} with sander/pmemd, set icnstph=1 "
                f"and solvph=X.X in your mdin"
                + (", and use the modified prmtop." if sim_type == 'explicit' else ".")
                + "[/grey50]")
            
            # Save to workspace automatically
            cpin_config = {
                'cpin_file': output_cpin_path,
                'prmtop_file': prmtop_file,
                'simulation_type': sim_type,
                'igb': igb if sim_type == 'implicit' else None,
                'intdiel': intdiel,
                'num_residues': len(selected_residues),
                'system_name': system_name,
                'selected_residues': selected_residues,
                'initial_states': initial_states
            }

            if sim_type == 'explicit':
                cpin_config['modified_prmtop'] = modified_prmtop
                # Swap the raw constph prmtop in md_structure_pairs to the
                # cpin-modified one — that's the production CpHMD topology
                # the MD manager should hand to pmemd/sander, not the
                # pre-cpinutil input.
                if self._swap_md_pair_to_cpin_prmtop(prmtop_file, modified_prmtop):
                    console.print(
                        f"[grey50]  Updated md_structure_pairs: "
                        f"{os.path.basename(prmtop_file)} → "
                        f"{os.path.basename(modified_prmtop)}[/grey50]")

            self.update_workspace("cpin_config", cpin_config)
            self.update_workspace("cpin_file", output_cpin_path)  # For menu status tracking
            console.print("[green]✓ Saved cpin settings to workspace[/green]")
            
            return True
            
        except subprocess.TimeoutExpired:
            console.print("[red]Error: cpinutil timed out after 5 minutes[/red]")
            return False
        except FileNotFoundError:
            console.print("[red]Error: cpinutil.py not found in PATH[/red]")
            console.print("Make sure AmberTools is installed and in your PATH")
            return False
        except Exception as e:
            console.print(f"[red]Error running cpinutil: {e}[/red]")
            return False

    def _get_all_titratable_residues_from_prmtop(self):
        """Fallback: Get all titratable residues when no protonation data available"""
        console = self.processor.console

        console.print("\n[yellow]No titratable residue data found in workspace.[/yellow]")
        console.print("\n[bold]To use this option, you need to:[/bold]")
        console.print("  1. Run the [blue]Protonation State Analyzer[/blue] module")
        console.print("  2. Select option 1: [blue]Analyze protonation states[/blue]")
        console.print("  3. Select option 4: [blue]Set titratable residues[/blue]")
        console.print("\nThis will populate the workspace with titratable residue data,")
        console.print("including residue types, pKa values, and predicted protonation states.")
        console.print("\n[yellow]Returning to main menu...[/yellow]\n")

        return None
