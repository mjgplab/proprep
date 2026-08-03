"""
Force Field Data Collector

Collects atom type definitions from AMBER force field libraries.
Parses leaprc files and extracts atom data from lib/off files.

Based on the approach in show_leaprc_files.py
"""

import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass

from rich.console import Console


@dataclass
class AtomTypeDefinition:
    """Atom type definition from force field library."""
    residue_name: str
    atom_name: str
    atom_type: str
    atom_charge: float

    def __repr__(self):
        return f"{self.residue_name}:{self.atom_name} -> {self.atom_type} ({self.atom_charge:.6f})"


@dataclass
class BondParameter:
    """Bond parameter from frcmod file."""
    atom_type1: str
    atom_type2: str
    force_constant: float  # kcal/mol/Å²
    eq_length: float       # Å
    source: str = ""


@dataclass
class AngleParameter:
    """Angle parameter from frcmod file."""
    atom_type1: str
    atom_type2: str
    atom_type3: str
    force_constant: float  # kcal/mol/rad²
    eq_angle: float        # degrees
    source: str = ""


class ForcefieldDataCollector:
    """
    Collects atom type definitions from AMBER force field libraries.

    Parses leaprc files to find lib/off files, then extracts atom type
    definitions for all standard residues.
    """

    def __init__(self, force_field: str = 'ff19SB', console: Optional[Console] = None, processor=None):
        """
        Initialize force field data collector.

        Args:
            force_field: Force field name (e.g., 'ff19SB', 'ff14SB') - used as default only
            console: Rich console for output
            processor: ProcessingModule for workspace access (optional)
        """
        self.force_field = force_field  # Will be updated when user selects leaprc files
        self.selected_leaprcs = []  # List of selected leaprc file names
        self.selected_frcmod_files = []  # List of .frcmod file paths from leaprc files
        self.selected_dat_files = []  # List of base .dat file paths (e.g., parm10.dat)
        self.console = console or Console()
        self.processor = processor  # For workspace access
        self.amberhome = self._get_amberhome()
        self.atom_definitions: Dict[Tuple[str, str], AtomTypeDefinition] = {}

        # Bond and angle parameters from frcmod files (unified source)
        self.bond_parameters: Dict[str, BondParameter] = {}  # "type1-type2" -> BondParameter
        self.angle_parameters: Dict[str, AngleParameter] = {}  # "type1-type2-type3" -> AngleParameter

        # Residue aliases for mapping (e.g., HEM -> HEH)
        self.residue_aliases: Dict[str, str] = {}

    def _get_amberhome(self) -> Path:
        """Get AMBERHOME directory from environment."""
        amberhome = os.environ.get('AMBERHOME')
        if not amberhome:
            raise RuntimeError("AMBERHOME environment variable not set")

        amberhome_path = Path(amberhome)
        if not amberhome_path.exists():
            raise RuntimeError(f"AMBERHOME directory does not exist: {amberhome}")

        return amberhome_path

    def _get_leaprc_path(self) -> Path:
        """Get path to leaprc file for the specified force field."""
        # Try different naming conventions
        leap_cmd_dir = self.amberhome / 'dat' / 'leap' / 'cmd'

        possible_names = [
            f'leaprc.protein.{self.force_field}',
            f'leaprc.{self.force_field}',
            f'oldff/leaprc.{self.force_field}',
        ]

        for name in possible_names:
            leaprc_path = leap_cmd_dir / name
            if leaprc_path.exists():
                return leaprc_path

        raise FileNotFoundError(
            f"Could not find leaprc file for force field: {self.force_field}\n"
            f"Searched in: {leap_cmd_dir}\n"
            f"Tried names: {', '.join(possible_names)}"
        )

    def _extract_loaded_files(self, leaprc_path: Path) -> List[str]:
        """
        Extract files loaded by a leaprc file.

        Args:
            leaprc_path: Path to leaprc file

        Returns:
            List of loaded file names
        """
        try:
            with open(leaprc_path, 'r') as f:
                content = f.read()

            loaded_files = []

            # Patterns for different load commands
            patterns = [
                r'loadOff\s+([^\s#]+)',
                r'loadAmberParams\s+([^\s#]+)',
                r'loadamberparams\s*=\s*loadamberparams\s+([^\s#]+)',
                r'loadoff\s+([^\s#]+)',
                r'source\s+([^\s#]+)',
                r'addPdbResMap\s*\{\s*\{[^}]+\}\s+[^}]+\}\s+([^\s#]+)',
                r'loadAmberPrep\s+([^\s#]+)',
                r'loadMol2\s+([^\s#]+)',
                r'loadPdb\s+([^\s#]+)',
            ]

            for pattern in patterns:
                matches = re.finditer(pattern, content, re.IGNORECASE | re.MULTILINE)
                for match in matches:
                    filename = match.group(1).strip()
                    # Filter out invalid entries
                    if filename and not filename.startswith('#') and filename not in ['{', '}', '']:
                        loaded_files.append(filename)

            return loaded_files

        except Exception as e:
            self.console.print(f"[yellow]Warning: Could not read leaprc file: {e}[/yellow]")
            return []

    def _find_lib_file(self, filename: str, leaprc_path: Path) -> Optional[Path]:
        """
        Find the full path to a lib/off file.

        Args:
            filename: Name of lib file
            leaprc_path: Path to leaprc file

        Returns:
            Full path to lib file or None
        """
        # Only look for .lib and .off files
        if not (filename.endswith('.lib') or filename.endswith('.off')):
            return None

        # Get the base AMBER directory structure
        # leaprc files can be in cmd/ or cmd/oldff/
        # lib files can be in lib/ or lib/oldff/
        # We need to find the 'leap' directory which is parent of 'cmd'
        leap_cmd_dir = leaprc_path.parent

        # Go up to the 'leap' directory
        if leap_cmd_dir.name == 'oldff':
            # leaprc is in cmd/oldff/, so go up twice to get to leap/
            leap_dir = leap_cmd_dir.parent.parent
        elif leap_cmd_dir.name == 'cmd':
            # leaprc is in cmd/, so go up once to get to leap/
            leap_dir = leap_cmd_dir.parent
        else:
            # Unexpected structure, fall back to old behavior
            leap_dir = leap_cmd_dir.parent

        lib_dir = leap_dir / 'lib'
        parm_dir = leap_dir / 'parm'

        # Search in common locations
        search_paths = [
            leap_cmd_dir / filename,        # Same directory as leaprc
            lib_dir / filename,              # leap/lib/filename
            lib_dir / 'oldff' / filename,    # leap/lib/oldff/filename (for oldff files)
            parm_dir / filename,             # leap/parm/filename
        ]

        for path in search_paths:
            if path.exists():
                return path

        return None

    def _find_frcmod_file(self, filename: str, leaprc_path: Path) -> Optional[Path]:
        """
        Find the full path to a .frcmod file.

        Args:
            filename: Name of frcmod file
            leaprc_path: Path to leaprc file

        Returns:
            Full path to frcmod file or None
        """
        # Look for .frcmod, .dat files, or files starting with "frcmod."
        # Note: Some AMBER frcmod files are named "frcmod.XXX" without an extension
        is_param_file = (filename.endswith('.frcmod') or
                        filename.endswith('.dat') or
                        filename.startswith('frcmod.'))

        if not is_param_file:
            return None

        leap_cmd_dir = leaprc_path.parent

        # Go up to the 'leap' directory
        if leap_cmd_dir.name == 'oldff':
            leap_dir = leap_cmd_dir.parent.parent
        elif leap_cmd_dir.name == 'cmd':
            leap_dir = leap_cmd_dir.parent
        else:
            leap_dir = leap_cmd_dir.parent

        parm_dir = leap_dir / 'parm'
        lib_dir = leap_dir / 'lib'

        # Search in common locations
        search_paths = [
            leap_cmd_dir / filename,        # Same directory as leaprc
            parm_dir / filename,             # leap/parm/filename (most common)
            parm_dir / 'oldff' / filename,   # leap/parm/oldff/filename
            lib_dir / filename,              # leap/lib/filename
        ]

        for path in search_paths:
            if path.exists():
                return path

        return None

    def _parse_lib_file(self, lib_path: Path) -> List[AtomTypeDefinition]:
        """
        Parse an AMBER lib (OFF format) file and extract atom information.

        Args:
            lib_path: Path to lib file

        Returns:
            List of AtomTypeDefinition objects
        """
        atoms_data = []

        try:
            with open(lib_path, 'r') as f:
                content = f.read()

            # Find all unit entries (residues)
            # Pattern to find unit.atoms tables
            atom_table_pattern = r'!entry\.([^.]+)\.unit\.atoms table[^\n]*\n((?:^[^!].*\n)*)'

            matches = re.finditer(atom_table_pattern, content, re.MULTILINE)

            for match in matches:
                residue_name = match.group(1)
                atoms_section = match.group(2)

                # Parse each atom line
                # Format: "atom_name" "atom_type" typex resx flags seq elmnt charge
                atom_line_pattern = r'"([^"]+)"\s+"([^"]+)"\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+([-\d.]+)'

                atom_lines = re.finditer(atom_line_pattern, atoms_section)

                for atom_match in atom_lines:
                    atom_name = atom_match.group(1)
                    atom_type = atom_match.group(2)
                    atom_charge = float(atom_match.group(3))

                    atoms_data.append(AtomTypeDefinition(
                        residue_name=residue_name,
                        atom_name=atom_name,
                        atom_type=atom_type,
                        atom_charge=atom_charge
                    ))

            return atoms_data

        except FileNotFoundError:
            return []
        except Exception as e:
            self.console.print(f"[yellow]Warning: Error parsing {lib_path}: {e}[/yellow]")
            return []

    def _parse_prep_file(self, prep_path: Path) -> List[AtomTypeDefinition]:
        """
        Parse an AMBER .prep (internal coordinate) file and extract atom information.

        Args:
            prep_path: Path to prep file

        Returns:
            List of AtomTypeDefinition objects
        """
        try:
            from proprep.oniom_prep.amber_lib_parser import AmberPrepParser

            parser = AmberPrepParser(str(prep_path))
            parser.parse()

            atoms_data = []
            for resname, res_data in parser.residues.items():
                for atom in res_data.get("atoms", []):
                    atoms_data.append(AtomTypeDefinition(
                        residue_name=resname,
                        atom_name=atom["name"],
                        atom_type=atom["type"],
                        atom_charge=atom["charge"]
                    ))

            return atoms_data

        except FileNotFoundError:
            return []
        except Exception as e:
            self.console.print(f"[yellow]Warning: Error parsing {prep_path}: {e}[/yellow]")
            return []

    def _get_leaprc_category(self, filename: str) -> str:
        """
        Get category for a leaprc file based on its name.

        Args:
            filename: Name of leaprc file

        Returns:
            Category string
        """
        name_lower = filename.lower()

        if 'dna' in name_lower:
            return "DNA"
        elif 'glycam' in name_lower:
            return "Carbohydrate"
        elif 'rna' in name_lower or 'modrna' in name_lower:
            return "RNA"
        elif 'amberdyes' in name_lower:
            return "Dye"
        elif 'constph' in name_lower or 'conste' in name_lower:
            return "Constant pH/Redox"
        elif 'gaff' in name_lower:
            return "GAFF (General)"
        elif 'lipid' in name_lower:
            return "Lipid"
        elif 'protein' in name_lower or 'phos' in name_lower or 'mimetic' in name_lower or 'fluorine' in name_lower:
            return "Protein"
        elif 'water' in name_lower:
            return "Water"
        else:
            return "Other"

    def select_and_collect_force_field_data(self) -> Dict[Tuple[str, str], AtomTypeDefinition]:
        """
        Three-step force field selection:
        1. User selects standard force field categories (Protein, Water, GAFF, etc.)
        2. User selects specific leaprc files from those categories
        3. User optionally loads custom parameter files (lib/prep/frcmod)

        Returns:
            Dictionary mapping (residue_name, atom_name) to AtomTypeDefinition
        """
        from rich.table import Table
        from proprep.utils.prompts import prompt_with_context, confirm_with_context

        leap_cmd_dir = self.amberhome / 'dat' / 'leap' / 'cmd'

        # Get all leaprc files
        leaprc_files = sorted([f for f in leap_cmd_dir.iterdir() if f.name.startswith('leaprc')])

        if not leaprc_files:
            raise FileNotFoundError(f"No leaprc files found in {leap_cmd_dir}")

        # Categorize all files
        categorized = {}
        for idx, leaprc_file in enumerate(leaprc_files, 1):
            category = self._get_leaprc_category(leaprc_file.name)
            if category not in categorized:
                categorized[category] = []
            categorized[category].append((idx, leaprc_file))

        # ============================================================
        # STEP 1: Select Standard Force Field Categories
        # ============================================================
        category_order = ["Protein", "DNA", "RNA", "Water", "GAFF (General)",
                         "Lipid", "Carbohydrate", "Dye", "Constant pH/Redox", "Other"]

        # Filter to only categories that exist
        available_categories = [cat for cat in category_order if cat in categorized]

        # Show category selection table
        cat_table = Table(
            title="Select Standard Force Field Categories",
            show_header=True,
            header_style="bold cyan"
        )
        cat_table.add_column("Index", style="yellow", no_wrap=True)
        cat_table.add_column("Category", style="cyan")
        cat_table.add_column("Available Files", style="green", justify="right")

        for cat_idx, category in enumerate(available_categories, 1):
            file_count = len(categorized[category])
            cat_table.add_row(str(cat_idx), category, str(file_count))

        self.console.print(cat_table)
        self.console.print("\n[grey50]Common combinations:[/grey50]")
        self.console.print("[grey50]  • Protein simulations: 1,4 (Protein + Water)[/grey50]")
        self.console.print("[grey50]  • Protein with ligands: 1,4,5 (Protein + Water + GAFF)[/grey50]")
        self.console.print("[grey50]  • DNA simulations: 2,4 (DNA + Water)[/grey50]\n")

        while True:
            cat_selection = prompt_with_context(
                self.processor,
                "Select categories (comma-separated)",
                default="1,4",
                module="Force Field Selection",
                description="Select force field categories",
            )

            # Parse category selections
            selected_category_indices = []
            try:
                for s in cat_selection.split(','):
                    idx = int(s.strip())
                    if idx < 1 or idx > len(available_categories):
                        raise ValueError(
                            f"Invalid category: {idx}. Please enter 1-{len(available_categories)}"
                        )
                    selected_category_indices.append(idx)
                break  # Valid input, exit loop
            except ValueError as e:
                self.console.print(f"[red]{e}. Please try again.[/red]\n")

        selected_categories = [available_categories[i - 1] for i in selected_category_indices]

        # ============================================================
        # STEP 2: Select Specific Force Fields from Categories
        # ============================================================
        # Collect all files from selected categories
        filtered_files = []
        for category in selected_categories:
            filtered_files.extend(categorized[category])

        # Sort by original index to maintain order
        filtered_files.sort(key=lambda x: x[0])

        # Create mapping from new index to original index
        display_idx_to_original = {}
        for display_idx, (original_idx, leaprc_file) in enumerate(filtered_files, 1):
            display_idx_to_original[display_idx] = original_idx

        # Show filtered force field table
        ff_table = Table(
            title=f"Select Force Fields from {', '.join(selected_categories)}",
            show_header=True,
            header_style="bold magenta",
            show_lines=True
        )
        ff_table.add_column("Index", style="yellow", no_wrap=True, width=6)
        ff_table.add_column("Category", style="magenta", no_wrap=True, width=18)
        ff_table.add_column("Leaprc File", style="cyan", no_wrap=True)
        ff_table.add_column("Loaded Files", style="green")

        for display_idx, (original_idx, leaprc_file) in enumerate(filtered_files, 1):
            category = self._get_leaprc_category(leaprc_file.name)
            loaded_files = self._extract_loaded_files(leaprc_file)

            if loaded_files:
                # Show only first 3 files to keep table compact
                files_str = "\n".join(loaded_files[:3])
                if len(loaded_files) > 3:
                    files_str += f"\n[grey50]... +{len(loaded_files)-3} more[/grey50]"
            else:
                files_str = "[grey50](none)[/grey50]"

            ff_table.add_row(str(display_idx), category, leaprc_file.name, files_str)

        self.console.print("\n")
        self.console.print(ff_table)
        self.console.print(f"\n[bold]Showing {len(filtered_files)} force field(s)[/bold]\n")

        # Prompt for force field selections
        self.console.print("[cyan]Select specific force fields (comma-separated indices)[/cyan]")
        self.console.print("[grey50]Tip: You can typically select one from each category[/grey50]\n")

        # Generate smart default based on selected categories
        smart_default = self._get_smart_default(selected_categories, filtered_files)

        ff_selection = prompt_with_context(
            self.processor,
            "Select force field(s)",
            default=smart_default,
            module="Force Field Selection",
            description="Select specific leaprc files",
        )

        # Parse force field selections
        ff_selections = []
        try:
            for s in ff_selection.split(','):
                display_idx = int(s.strip())
                if display_idx < 1 or display_idx > len(filtered_files):
                    raise ValueError(f"Invalid selection: {display_idx}. Please enter 1-{len(filtered_files)}")
                # Map display index back to original index
                original_idx = display_idx_to_original[display_idx]
                ff_selections.append(original_idx)
        except ValueError as e:
            raise ValueError(f"Invalid force field selection: {e}")

        # Get selected leaprc files using original indices
        selected_leaprcs = [leaprc_files[idx - 1] for idx in ff_selections]

        # Store selected leaprc filenames for display in error messages
        self.selected_leaprcs = [f.name for f in selected_leaprcs]

        self.console.print(f"\n[green]✓ Standard force fields: {', '.join(self.selected_leaprcs)}[/green]")

        # Collect data from all selected leaprc files
        for leaprc in selected_leaprcs:
            self._collect_from_leaprc(leaprc)

        # ============================================================
        # STEP 3: Custom Parameter Files (always offered)
        # ============================================================
        # Check if there are discoverable custom files before prompting
        discovered_files = self._discover_custom_files()
        has_discoverable = bool(discovered_files['lib'] or discovered_files['frcmod'])

        if has_discoverable:
            n_discovered = len(discovered_files['lib']) + len(discovered_files['frcmod'])
            self.console.print(f"\n[bold cyan]Custom Parameters[/bold cyan]")
            self.console.print(f"[grey50]Found {n_discovered} custom parameter file(s) in workspace/directory[/grey50]")
            self._load_custom_files(discovered_files=discovered_files)
        else:
            load_custom = confirm_with_context(
                self.processor,
                "Load custom parameter files (lib/prep/frcmod)?",
                default=False,
                module="Force Field Selection",
                description="Load custom parameter files",
            )
            if load_custom:
                self.console.print(f"\n[bold cyan]Custom Parameters[/bold cyan]")
                self._load_custom_files(discovered_files=discovered_files)

        return self.atom_definitions

    def _get_smart_default(self, selected_categories: List[str], filtered_files: List) -> str:
        """
        Generate smart default selections based on common patterns.

        Args:
            selected_categories: List of selected category names
            filtered_files: List of (original_idx, leaprc_file) tuples

        Returns:
            Comma-separated string of display indices
        """
        defaults = []

        # Look for common force field patterns
        for display_idx, (original_idx, leaprc_file) in enumerate(filtered_files, 1):
            filename = leaprc_file.name.lower()

            # Protein: prefer ff14SB
            if "Protein" in selected_categories and "ff14sb" in filename and "protein" in filename:
                defaults.append(display_idx)
            # Water: prefer tip3p
            elif "Water" in selected_categories and "tip3p" in filename:
                defaults.append(display_idx)
            # GAFF: prefer gaff2
            elif "GAFF (General)" in selected_categories and filename == "leaprc.gaff2":
                defaults.append(display_idx)
            # DNA: prefer OL15
            elif "DNA" in selected_categories and "ol15" in filename:
                defaults.append(display_idx)
            # RNA: prefer OL3
            elif "RNA" in selected_categories and "ol3" in filename:
                defaults.append(display_idx)

        # If we found smart defaults, return them
        if defaults:
            return ",".join(map(str, defaults))

        # Otherwise, default to first item from each category
        seen_categories = set()
        for display_idx, (original_idx, leaprc_file) in enumerate(filtered_files, 1):
            category = self._get_leaprc_category(leaprc_file.name)
            if category not in seen_categories:
                defaults.append(display_idx)
                seen_categories.add(category)

        return ",".join(map(str, defaults))

    def _load_custom_files(self, discovered_files=None):
        """
        Prompt user for custom lib/frcmod files and load them.

        Supports loading user-generated parameter files (e.g., from small molecule
        parameterizer or external tools like antechamber).

        Features:
        - Checks workspace for files from previous parameterizations
        - Scans current directory for lib/frcmod files
        - Allows manual path entry

        Args:
            discovered_files: Pre-computed discovery results to avoid redundant scanning.
                              If None, will run discovery automatically.
        """
        from proprep.utils.prompts import prompt_with_context, confirm_with_context, Confirm, IntPrompt
        from rich.table import Table

        self.console.print("[grey50]Load custom lib (.lib/.off/.prep) and/or frcmod files[/grey50]\n")

        # Track loaded files
        self.custom_lib_files = []
        self.custom_frcmod_files = []

        # Discover available files (use pre-computed if provided)
        if discovered_files is None:
            discovered_files = self._discover_custom_files()

        if discovered_files['lib'] or discovered_files['frcmod']:
            # Show discovered files in a table
            table = Table(title="Discovered Parameter Files", show_header=True)
            table.add_column("#", style="yellow", width=4)
            table.add_column("Type", style="cyan", width=8)
            table.add_column("File", style="green")
            table.add_column("Source", style="grey50")

            file_list = []
            idx = 1
            for lib_file, source in discovered_files['lib']:
                table.add_row(str(idx), "lib", Path(lib_file).name, source)
                file_list.append(('lib', lib_file))
                idx += 1
            for frcmod_file, source in discovered_files['frcmod']:
                table.add_row(str(idx), "frcmod", Path(frcmod_file).name, source)
                file_list.append(('frcmod', frcmod_file))
                idx += 1

            self.console.print(table)
            self.console.print("\n[grey50]Enter file numbers to load (comma-separated), 'all', 'none', or 'manual'[/grey50]")

            selection = prompt_with_context(
                self.processor,
                "Select files",
                default="all",
                module="Custom Parameters",
                description="Select custom parameter files",
            )

            if selection.lower() == 'all':
                selected_indices = list(range(len(file_list)))
            elif selection.lower() == 'none':
                selected_indices = []
            elif selection.lower() == 'manual':
                selected_indices = []  # Will prompt for manual entry below
            else:
                try:
                    selected_indices = [int(s.strip()) - 1 for s in selection.split(',')]
                except ValueError:
                    self.console.print("[yellow]Invalid selection, skipping discovered files[/yellow]")
                    selected_indices = []

            # Load selected discovered files
            for idx in selected_indices:
                if 0 <= idx < len(file_list):
                    file_type, file_path = file_list[idx]
                    self._load_single_file(file_type, file_path)

            # Ask for manual entry if user selected 'manual' or wants to add more
            if selection.lower() == 'manual' or (selected_indices and confirm_with_context(
                    self.processor,
                    "Add more files manually?",
                    default=False,
                    module="Custom Parameters",
                    description="Add more files manually",
                )):
                self._prompt_manual_file_entry()
        else:
            # No discovered files - go straight to manual entry
            self.console.print("[grey50]No lib/frcmod files found in workspace or current directory[/grey50]\n")
            self._prompt_manual_file_entry()

        # Summary
        total_custom = len(self.custom_lib_files) + len(self.custom_frcmod_files)
        if total_custom > 0:
            self.console.print(f"\n[green]✅ Loaded {total_custom} custom file(s)[/green]")
        else:
            self.console.print("\n[yellow]No custom files loaded[/yellow]")

    def _discover_custom_files(self) -> dict:
        """
        Discover lib/frcmod files from workspace and current directory.

        Returns:
            Dict with 'lib' and 'frcmod' lists of (path, source) tuples
        """
        discovered = {'lib': [], 'frcmod': []}

        # 1. Check workspace for files from previous parameterizations
        if self.processor:
            try:
                param_residues = self.processor.get_from_workspace("parameterized_residues", {})
                for res_name, res_data in param_residues.items():
                    output_files = res_data.get("output_files", {})
                    output_dir = res_data.get("output_dir", "")

                    # Look for lib files
                    for key in ['lib_file', 'prep_file', 'off_file']:
                        if key in output_files:
                            lib_path = Path(output_files[key])
                            if not lib_path.is_absolute() and output_dir:
                                lib_path = Path(output_dir) / lib_path
                            if lib_path.exists():
                                discovered['lib'].append((str(lib_path), f"workspace:{res_name}"))

                    # Look for frcmod files
                    for key in ['frcmod_file', 'frcmod']:
                        if key in output_files:
                            frcmod_path = Path(output_files[key])
                            if not frcmod_path.is_absolute() and output_dir:
                                frcmod_path = Path(output_dir) / frcmod_path
                            if frcmod_path.exists():
                                discovered['frcmod'].append((str(frcmod_path), f"workspace:{res_name}"))

                # Also check small_molecules in workspace
                small_mols = self.processor.get_from_workspace("small_molecules", [])
                for mol in small_mols:
                    if isinstance(mol, dict):
                        for key in ['lib_file', 'off_file']:
                            if key in mol and Path(mol[key]).exists():
                                discovered['lib'].append((mol[key], "workspace:small_molecule"))
                        if 'frcmod_file' in mol and Path(mol['frcmod_file']).exists():
                            discovered['frcmod'].append((mol['frcmod_file'], "workspace:small_molecule"))

            except Exception:
                pass  # Silently ignore workspace access errors

        # 2. Scan current directory and subdirectories (1 level deep)
        cwd = Path.cwd()
        search_dirs = [cwd] + [d for d in cwd.iterdir() if d.is_dir() and not d.name.startswith('.')]

        for search_dir in search_dirs:
            # Find lib/off/prep files
            for lib_file in search_dir.glob("*.lib"):
                path_str = str(lib_file)
                if not any(p == path_str for p, _ in discovered['lib']):
                    rel_path = lib_file.relative_to(cwd) if lib_file.is_relative_to(cwd) else lib_file
                    discovered['lib'].append((path_str, f"directory:{rel_path.parent}"))

            for off_file in search_dir.glob("*.off"):
                path_str = str(off_file)
                if not any(p == path_str for p, _ in discovered['lib']):
                    rel_path = off_file.relative_to(cwd) if off_file.is_relative_to(cwd) else off_file
                    discovered['lib'].append((path_str, f"directory:{rel_path.parent}"))

            for prep_file in search_dir.glob("*.prep"):
                path_str = str(prep_file)
                if not any(p == path_str for p, _ in discovered['lib']):
                    rel_path = prep_file.relative_to(cwd) if prep_file.is_relative_to(cwd) else prep_file
                    discovered['lib'].append((path_str, f"directory:{rel_path.parent}"))

            # Find frcmod files (both *.frcmod and frcmod.* naming conventions)
            for pattern in ["*.frcmod", "frcmod.*"]:
                for frcmod_file in search_dir.glob(pattern):
                    path_str = str(frcmod_file)
                    if not any(p == path_str for p, _ in discovered['frcmod']):
                        rel_path = frcmod_file.relative_to(cwd) if frcmod_file.is_relative_to(cwd) else frcmod_file
                        discovered['frcmod'].append((path_str, f"directory:{rel_path.parent}"))

        return discovered

    def _load_single_file(self, file_type: str, file_path: str):
        """Load a single lib or frcmod file."""
        path = Path(file_path)
        if not path.exists():
            self.console.print(f"[red]File not found: {file_path}[/red]")
            return

        if file_type == 'lib':
            self.console.print(f"[green]Loading lib: {path.name}[/green]")
            if path.suffix == '.prep':
                atoms = self._parse_prep_file(path)
            else:
                atoms = self._parse_lib_file(path)
            for atom_def in atoms:
                key = (atom_def.residue_name, atom_def.atom_name)
                self.atom_definitions[key] = atom_def
            self.custom_lib_files.append(str(path))
            self.console.print(f"  → Loaded {len(atoms)} atom definitions")
        elif file_type == 'frcmod':
            self.console.print(f"[green]Loading frcmod: {path.name}[/green]")
            self._load_frcmod_file(path)
            self.custom_frcmod_files.append(str(path))

    def _prompt_manual_file_entry(self):
        """Prompt user for manual file path entry."""
        from proprep.utils.prompts import prompt_with_context, confirm_with_context, Confirm

        # Ask for lib files
        if confirm_with_context(
            self.processor,
            "Enter lib/off file path(s) manually?",
            default=False,
            module="Custom Parameters",
            description="Enter lib file paths manually",
        ):
            lib_paths = prompt_with_context(
                self.processor,
                "Enter lib file path(s) (comma-separated)",
                default="",
                module="Custom Parameters",
                description="Lib file paths",
            )
            if lib_paths.strip():
                for lib_path_str in lib_paths.split(','):
                    lib_path = Path(lib_path_str.strip()).expanduser().resolve()
                    self._load_single_file('lib', str(lib_path))

        # Ask for frcmod files
        if confirm_with_context(
            self.processor,
            "Enter frcmod file path(s) manually?",
            default=False,
            module="Custom Parameters",
            description="Enter frcmod file paths manually",
        ):
            frcmod_paths = prompt_with_context(
                self.processor,
                "Enter frcmod file path(s) (comma-separated)",
                default="",
                module="Custom Parameters",
                description="Frcmod file paths",
            )
            if frcmod_paths.strip():
                for frcmod_path_str in frcmod_paths.split(','):
                    frcmod_path = Path(frcmod_path_str.strip()).expanduser().resolve()
                    self._load_single_file('frcmod', str(frcmod_path))

    def _load_frcmod_file(self, frcmod_path: Path, source: str = None):
        """
        Load parameters from an frcmod file.

        Stores bonded parameters in the unified bond_parameters and angle_parameters
        dictionaries for use in pre-frcmod generation and parameter lookup.

        Args:
            frcmod_path: Path to frcmod file
            source: Source identifier (defaults to filename)
        """
        if source is None:
            source = frcmod_path.name

        try:
            with open(frcmod_path, 'r') as f:
                lines = f.readlines()

            current_section = None
            bond_count = 0
            angle_count = 0

            for line in lines:
                line = line.rstrip()
                if not line or line.startswith('REMARK'):
                    continue

                # Detect section headers
                if line.startswith('MASS'):
                    current_section = 'masses'
                    continue
                elif line.startswith('BOND'):
                    current_section = 'bonds'
                    continue
                elif line.startswith('ANGL'):
                    current_section = 'angles'
                    continue
                elif line.startswith('DIHE'):
                    current_section = 'dihedrals'
                    continue
                elif line.startswith('IMPR'):
                    current_section = 'impropers'
                    continue
                elif line.startswith('NONB'):
                    current_section = 'nonbonded'
                    continue

                # Parse based on section
                if current_section == 'bonds' and len(line) >= 5 and '-' in line[:5]:
                    # Format: AA-BB  force_constant  eq_length
                    # Fixed-width: chars 0-1 type1, char 2 dash, chars 3-4 type2
                    try:
                        type1 = line[0:2].strip()
                        type2 = line[3:5].strip()
                        values = line[5:].split()
                        if len(values) >= 2 and type1 and type2:
                            force_constant = float(values[0])
                            eq_length = float(values[1])

                            # Store with both orderings
                            key1 = f"{type1}-{type2}"
                            key2 = f"{type2}-{type1}"
                            param = BondParameter(type1, type2, force_constant, eq_length, source)
                            self.bond_parameters[key1] = param
                            self.bond_parameters[key2] = param
                            bond_count += 1
                    except (ValueError, IndexError):
                        pass

                elif current_section == 'angles' and len(line) >= 8 and '-' in line[:8]:
                    # Format: AA-BB-CC  force_constant  eq_angle
                    # Fixed-width: chars 0-1, 3-4, 6-7 are atom types
                    try:
                        type1 = line[0:2].strip()
                        type2 = line[3:5].strip()
                        type3 = line[6:8].strip()
                        values = line[8:].split()
                        if len(values) >= 2 and type1 and type2 and type3:
                            force_constant = float(values[0])
                            eq_angle = float(values[1])

                            # Store with both orderings (reversible)
                            key1 = f"{type1}-{type2}-{type3}"
                            key2 = f"{type3}-{type2}-{type1}"
                            param = AngleParameter(type1, type2, type3, force_constant, eq_angle, source)
                            self.angle_parameters[key1] = param
                            self.angle_parameters[key2] = param
                            angle_count += 1
                    except (ValueError, IndexError):
                        pass

            self.console.print(f"  → Loaded {bond_count} bond, {angle_count} angle parameters from {source}")

        except Exception as e:
            self.console.print(f"[yellow]Warning: Error parsing {frcmod_path}: {e}[/yellow]")

    def _collect_from_leaprc(self, leaprc_path: Path, _processed_leaprcs: Optional[set] = None, _lib_count: Optional[list] = None) -> Dict[Tuple[str, str], AtomTypeDefinition]:
        """
        Collect atom definitions from a specific leaprc file.
        Recursively follows 'source' directives to load chained force fields.

        Args:
            leaprc_path: Path to leaprc file
            _processed_leaprcs: Set of already processed leaprc paths (for recursion tracking)
            _lib_count: List containing single int for tracking total lib files across recursion

        Returns:
            Dictionary of atom definitions
        """
        # Initialize on first call
        is_first_call = False
        if _processed_leaprcs is None:
            _processed_leaprcs = set()
            _lib_count = [0]  # Use list so it's mutable across recursive calls
            is_first_call = True
            self.console.print(f"[cyan]Collecting force field data from {leaprc_path.name}...[/cyan]")

        # Skip if already processed (prevents infinite loops)
        if leaprc_path in _processed_leaprcs:
            return self.atom_definitions

        _processed_leaprcs.add(leaprc_path)

        # Extract loaded files
        loaded_files = self._extract_loaded_files(leaprc_path)

        # Process files: separate lib files, frcmod files, and sourced leaprc files
        for filename in loaded_files:
            # Check if this is a sourced leaprc file
            # Examples: 'leaprc.constph', 'oldff/leaprc.ff10'
            is_leaprc = (filename.startswith('leaprc') or
                        ('/' in filename and 'leaprc' in filename))

            if is_leaprc:
                # This is a sourced leaprc - find and recursively process it
                sourced_leaprc = self._find_sourced_leaprc(filename, leaprc_path)
                if sourced_leaprc:
                    self.console.print(f"  [grey50]→ Following source: {sourced_leaprc.name}[/grey50]")
                    self._collect_from_leaprc(sourced_leaprc, _processed_leaprcs, _lib_count)
            elif (filename.endswith('.dat') and filename.startswith('parm')):
                # Base parameter file (e.g., parm10.dat) — different format from frcmod
                # Collect separately for ONIOM writer to load VDW and bonded parameters
                dat_path = self._find_frcmod_file(filename, leaprc_path)
                if dat_path:
                    if str(dat_path) not in self.selected_dat_files:
                        self.selected_dat_files.append(str(dat_path))
            elif (filename.endswith('.frcmod') or filename.startswith('frcmod.') or
                  filename.endswith('.dat')):
                # This is a parameter file (frcmod or .dat like amberdyes.dat)
                # Note: Some frcmod files are named "frcmod.XXX" without .frcmod extension
                frcmod_path = self._find_frcmod_file(filename, leaprc_path)
                if frcmod_path:
                    if str(frcmod_path) not in self.selected_frcmod_files:
                        self.selected_frcmod_files.append(str(frcmod_path))
                        # Actually parse the frcmod file to load bond/angle parameters
                        self._load_frcmod_file(frcmod_path, source=frcmod_path.name)
                else:
                    self.console.print(f"  [yellow]⚠ Could not find parameter file: {filename}[/yellow]")
            else:
                # This is a lib file
                lib_path = self._find_lib_file(filename, leaprc_path)
                if lib_path:
                    _lib_count[0] += 1
                    atoms = self._parse_lib_file(lib_path)

                    for atom_def in atoms:
                        key = (atom_def.residue_name, atom_def.atom_name)
                        self.atom_definitions[key] = atom_def

        # Only print summary for the top-level call
        if is_first_call:
            total_leaprcs = len(_processed_leaprcs)
            self.console.print(f"[green]✓ Loaded {len(self.atom_definitions)} atom definitions from {_lib_count[0]} library files[/green]")
            if total_leaprcs > 1:
                self.console.print(f"[grey50]  (followed {total_leaprcs - 1} sourced leaprc file(s))[/grey50]")

            # Show parameter summary
            if self.bond_parameters or self.angle_parameters:
                # Divide by 2 since we store both orderings
                n_bonds = len(self.bond_parameters) // 2
                n_angles = len(self.angle_parameters) // 2
                self.console.print(f"[green]✓ Loaded {n_bonds} bond, {n_angles} angle parameters from {len(self.selected_frcmod_files)} frcmod file(s)[/green]")

        return self.atom_definitions

    def _find_sourced_leaprc(self, filename: str, current_leaprc_path: Path) -> Optional[Path]:
        """
        Find a sourced leaprc file referenced in a source directive.

        Args:
            filename: Referenced file (e.g., 'oldff/leaprc.ff10')
            current_leaprc_path: Path to current leaprc file

        Returns:
            Path to sourced leaprc file or None
        """
        leap_cmd_dir = current_leaprc_path.parent

        # Try relative to current leaprc directory
        search_paths = [
            leap_cmd_dir / filename,
            # Also try without directory prefix if it has one
            leap_cmd_dir / Path(filename).name if '/' in filename else None
        ]

        for path in search_paths:
            if path and path.exists():
                return path

        return None

    def get_atom_type(self, residue_name: str, atom_name: str) -> Optional[AtomTypeDefinition]:
        """
        Get atom type definition for a specific residue and atom.

        Args:
            residue_name: Residue name (e.g., 'ALA', 'CYS')
            atom_name: Atom name (e.g., 'CA', 'CB')

        Returns:
            AtomTypeDefinition or None if not found
        """
        return self.atom_definitions.get((residue_name, atom_name))

    def get_residue_atoms(self, residue_name: str) -> List[AtomTypeDefinition]:
        """
        Get all atom type definitions for a residue.

        Args:
            residue_name: Residue name

        Returns:
            List of AtomTypeDefinition objects for the residue
        """
        return [
            atom_def for (res_name, _), atom_def in self.atom_definitions.items()
            if res_name == residue_name
        ]

    def has_residue(self, residue_name: str) -> bool:
        """
        Check if force field has definitions for a residue.

        Also checks residue aliases (e.g., HEM mapped to HEH).

        Args:
            residue_name: Residue name

        Returns:
            True if residue is defined in force field
        """
        # Check direct match
        if any(res_name == residue_name for res_name, _ in self.atom_definitions.keys()):
            return True

        # Check if there's an alias
        alias = self.residue_aliases.get(residue_name)
        if alias:
            return any(res_name == alias for res_name, _ in self.atom_definitions.keys())

        return False

    def get_bond_parameter(self, type1: str, type2: str) -> Optional[BondParameter]:
        """
        Get bond parameter for atom type pair from loaded frcmod files.

        Args:
            type1: First atom type
            type2: Second atom type

        Returns:
            BondParameter or None if not found
        """
        type1 = type1.strip()
        type2 = type2.strip()
        key = f"{type1}-{type2}"
        return self.bond_parameters.get(key)

    def get_angle_parameter(self, type1: str, type2: str, type3: str) -> Optional[AngleParameter]:
        """
        Get angle parameter for atom type triple from loaded frcmod files.

        Args:
            type1: First atom type
            type2: Second atom type (central)
            type3: Third atom type

        Returns:
            AngleParameter or None if not found
        """
        type1 = type1.strip()
        type2 = type2.strip()
        type3 = type3.strip()
        key = f"{type1}-{type2}-{type3}"
        return self.angle_parameters.get(key)

    def get_available_residues(self) -> List[str]:
        """
        Get list of all available residue names from loaded force field.

        Returns:
            Sorted list of unique residue names
        """
        residues = set()
        for res_name, _ in self.atom_definitions.keys():
            residues.add(res_name)
        return sorted(residues)

    def prompt_residue_mapping(self, unknown_residues: List[str], structure_atoms: Dict[str, List[str]] = None) -> Dict[str, str]:
        """
        Interactive interface for mapping unknown residues to available ones.

        Displays all available residues and lets user pick which ones to map
        unknown residues to (e.g., HEM -> HEH for heme).

        Args:
            unknown_residues: List of residue names not found in force field
            structure_atoms: Optional dict mapping residue names to list of atom names
                            in the actual structure (for checking which atoms are missing)

        Returns:
            Dict mapping unknown residue names to selected force field residues
        """
        from rich.table import Table
        from proprep.utils.prompts import prompt_with_context
        from rich.panel import Panel

        if not unknown_residues:
            return {}

        available = self.get_available_residues()
        if not available:
            self.console.print("[yellow]No residues loaded from force field[/yellow]")
            return {}

        # Show available residues in a compact table
        self.console.print(Panel(
            "[bold]Residue Mapping Required[/bold]\n"
            "Some residues in your structure are not in the loaded force field.\n"
            "You can map them to available residues (e.g., HEM → HEH for heme).",
            title="Residue Mapping"
        ))

        # Display available residues in columns
        table = Table(title="Available Residues in Force Field", show_header=True)
        table.add_column("#", style="yellow", width=4)
        table.add_column("Residue", style="cyan", width=8)
        table.add_column("#", style="yellow", width=4)
        table.add_column("Residue", style="cyan", width=8)
        table.add_column("#", style="yellow", width=4)
        table.add_column("Residue", style="cyan", width=8)
        table.add_column("#", style="yellow", width=4)
        table.add_column("Residue", style="cyan", width=8)

        # Fill table row by row (4 columns)
        for i in range(0, len(available), 4):
            row = []
            for j in range(4):
                if i + j < len(available):
                    row.extend([str(i + j + 1), available[i + j]])
                else:
                    row.extend(["", ""])
            table.add_row(*row)

        self.console.print(table)

        # Map each unknown residue
        mappings = {}
        for res_name in unknown_residues:
            self.console.print(f"\n[yellow]Unknown residue: {res_name}[/yellow]")
            self.console.print("[grey50]Tip: Enter multiple numbers separated by commas (e.g., '101,169' for HEH+PRN)[/grey50]")
            choice = prompt_with_context(
                self.processor,
                f"Enter number(s) to map '{res_name}' to, or 'skip' to leave unmapped",
                default="skip",
                module="Force Field Selection",
                description=f"Map unknown residue {res_name}",
            )

            if choice.lower() != 'skip':
                # Parse comma-separated indices
                try:
                    indices = [int(x.strip()) - 1 for x in choice.split(',')]
                    target_residues = []
                    valid = True

                    for idx in indices:
                        if 0 <= idx < len(available):
                            target_residues.append(available[idx])
                        else:
                            self.console.print(f"  [red]Invalid index: {idx + 1}[/red]")
                            valid = False
                            break

                    if valid and target_residues:
                        # Store the first as primary alias (for backward compatibility)
                        mappings[res_name] = target_residues[0]
                        self.residue_aliases[res_name] = target_residues[0]

                        # Copy atom definitions from ALL target residues
                        total_copied = 0
                        for target_res in target_residues:
                            copied = self._copy_residue_definitions(target_res, res_name)
                            total_copied += copied
                            self.console.print(f"  [green]✓ Copied {copied} atoms from {target_res}[/green]")

                        self.console.print(f"  [green]✓ Mapped {res_name} → {', '.join(target_residues)} ({total_copied} total atoms)[/green]")

                        # Check for missing atoms and offer atom-level mapping if needed
                        if structure_atoms and res_name in structure_atoms:
                            self._prompt_atom_level_mapping(res_name, target_residues, structure_atoms[res_name])
                except ValueError:
                    self.console.print(f"  [yellow]Invalid input, skipping[/yellow]")

        return mappings

    def _prompt_atom_level_mapping(self, target_res: str, source_residues: List[str], structure_atom_names: List[str]) -> None:
        """
        Prompt for atom-level mapping when atom names differ between structure and force field.

        For example, structure has HEM with atoms CAA, CBA, CGA but force field PRN
        defines CA, CB, CG. User can map: CAA=CA, CBA=CB, etc.

        Args:
            target_res: Structure residue name (e.g., 'HEM')
            source_residues: List of FF residue names mapped to (e.g., ['HEH', 'PRN'])
            structure_atom_names: List of atom names in the actual structure
        """
        from proprep.utils.prompts import prompt_with_context
        from rich.table import Table

        # Find which structure atoms still have no FF definition
        missing_atoms = []
        for atom_name in structure_atom_names:
            if (target_res, atom_name) not in self.atom_definitions:
                missing_atoms.append(atom_name)

        if not missing_atoms:
            self.console.print(f"  [green]✓ All {len(structure_atom_names)} structure atoms have FF definitions[/green]")
            return

        # Find which structure atoms already have definitions (matched by name)
        defined_struct_atoms = set()
        for atom_name in structure_atom_names:
            if (target_res, atom_name) in self.atom_definitions:
                defined_struct_atoms.add(atom_name)

        # Collect available FF atoms - only those NOT already used for name-matched definitions
        # This enforces 1:1 mapping: each FF atom can only define one structure atom
        ff_atoms_by_residue = {}
        all_ff_atoms = {}
        for source_res in source_residues:
            ff_atoms_by_residue[source_res] = []
            for (res_name, atom_name), atom_def in self.atom_definitions.items():
                if res_name == source_res:
                    # Only include if this FF atom name wasn't already used
                    if atom_name not in defined_struct_atoms:
                        ff_atoms_by_residue[source_res].append((atom_name, atom_def.atom_type))
                        all_ff_atoms[atom_name] = (source_res, atom_def)

        # Display: Structure atoms without FF definitions
        self.console.print(f"\n[yellow]⚠ {len(missing_atoms)} structure atoms have no FF definition:[/yellow]")
        self.console.print(f"  [yellow]{', '.join(sorted(missing_atoms))}[/yellow]")

        # Display: Available FF atoms in table format
        self.console.print(f"\n[cyan]Available FF atoms:[/cyan]")
        for source_res in source_residues:
            atoms = ff_atoms_by_residue.get(source_res, [])
            if atoms:
                self.console.print(f"  [cyan]{source_res}:[/cyan]")
                # Sort atoms and display in columns (4 per row)
                atoms_sorted = sorted(atoms, key=lambda x: x[0])
                row = []
                for name, atype in atoms_sorted:
                    row.append(f"{name:8s} ({atype})")
                    if len(row) == 4:
                        self.console.print(f"    {row[0]:20s} {row[1]:20s} {row[2]:20s} {row[3]:20s}")
                        row = []
                if row:  # Print remaining
                    self.console.print(f"    {''.join(f'{r:20s}' for r in row)}")

        self.console.print(f"\n[grey50]Unmapped atoms will need Seminario-derived bonded parameters.[/grey50]")
        self.console.print(f"[grey50]Enter mappings one per line (e.g., 'CAA=CA'). Empty line to finish.[/grey50]")

        # Collect mappings interactively, one per line
        atom_mappings = {}
        used_ff_atoms = set()

        while True:
            remaining = [a for a in missing_atoms if a not in atom_mappings]
            if not remaining:
                self.console.print(f"  [green]All structure atoms mapped![/green]")
                break

            mapping = prompt_with_context(
                self.processor,
                f"[grey50]Remaining: {', '.join(sorted(remaining)[:5])}{'...' if len(remaining) > 5 else ''}[/grey50]\nMapping",
                default="",
                module="Force Field Selection",
                description="Atom type mapping",
            )

            if not mapping.strip():
                break

            # Parse mapping (format: "CAA=CA")
            if '=' not in mapping:
                self.console.print(f"  [yellow]Invalid format. Use STRUCT_ATOM=FF_ATOM (e.g., CAA=CA)[/yellow]")
                continue

            struct_atom, ff_atom = mapping.split('=', 1)
            struct_atom = struct_atom.strip()
            ff_atom = ff_atom.strip()

            # Validate
            if struct_atom not in missing_atoms:
                self.console.print(f"  [yellow]'{struct_atom}' is not a missing structure atom[/yellow]")
                continue
            if struct_atom in atom_mappings:
                self.console.print(f"  [yellow]'{struct_atom}' already mapped[/yellow]")
                continue
            if ff_atom not in all_ff_atoms:
                self.console.print(f"  [yellow]'{ff_atom}' not found in available FF atoms[/yellow]")
                continue
            if ff_atom in used_ff_atoms:
                self.console.print(f"  [yellow]'{ff_atom}' already used for another mapping[/yellow]")
                continue

            # Apply mapping
            source_res, atom_def = all_ff_atoms[ff_atom]
            new_def = AtomTypeDefinition(
                residue_name=target_res,
                atom_name=struct_atom,
                atom_type=atom_def.atom_type,
                atom_charge=atom_def.atom_charge
            )
            self.atom_definitions[(target_res, struct_atom)] = new_def
            atom_mappings[struct_atom] = ff_atom
            used_ff_atoms.add(ff_atom)
            self.console.print(f"  [green]✓ {struct_atom} → {ff_atom} (type={atom_def.atom_type})[/green]")

        # Summary
        if atom_mappings:
            self.console.print(f"\n[green]✓ Mapped {len(atom_mappings)} atoms[/green]")

        still_missing = [a for a in missing_atoms if a not in atom_mappings]
        if still_missing:
            self.console.print(f"[grey50]Unmapped ({len(still_missing)}): {', '.join(sorted(still_missing))}[/grey50]")
            self.console.print(f"[grey50]Bonded parameters for these will come from Seminario.[/grey50]")

    def _copy_residue_definitions(self, source_res: str, target_res: str) -> int:
        """
        Copy atom definitions from source residue to target residue name.

        Args:
            source_res: Source residue name (e.g., 'HEH')
            target_res: Target residue name (e.g., 'HEM')

        Returns:
            Number of atom definitions copied
        """
        atoms_copied = 0
        for (res_name, atom_name), atom_def in list(self.atom_definitions.items()):
            if res_name == source_res:
                # Only copy if not already defined (avoid overwriting)
                if (target_res, atom_name) not in self.atom_definitions:
                    new_def = AtomTypeDefinition(
                        residue_name=target_res,
                        atom_name=atom_name,
                        atom_type=atom_def.atom_type,
                        atom_charge=atom_def.atom_charge
                    )
                    self.atom_definitions[(target_res, atom_name)] = new_def
                    atoms_copied += 1

        return atoms_copied

    def validate_redox_site_coverage(self, redox_sites: list, metal_elements: set = None) -> dict:
        """
        Validate that non-metal residues in RedoxSites have force field parameters.

        Checks each unique residue in the RedoxSites and reports which ones
        are missing from the loaded force field. Metal-containing residues
        are expected to get parameters from MCPB/Seminario, so they're checked
        but flagged differently.

        Args:
            redox_sites: List of RedoxSite objects to validate
            metal_elements: Set of metal element symbols (default: common metals)

        Returns:
            Dict with validation results:
                - valid: True if all non-metal residues have parameters
                - missing_residues: List of residues without parameters
                - metal_residues: List of metal-containing residues (expected to be parameterized)
                - covered_residues: List of residues with parameters
        """
        if metal_elements is None:
            metal_elements = {
                'FE', 'ZN', 'CU', 'MN', 'MG', 'CA', 'CO', 'NI', 'MO', 'W',
                'V', 'CR', 'CD', 'HG', 'PB', 'NA', 'K', 'LI', 'RB', 'CS',
                'BA', 'SR', 'AL', 'GA', 'IN', 'SN', 'TL', 'BI', 'PD', 'PT',
                'AG', 'AU', 'RU', 'RH', 'IR', 'OS', 'RE', 'TC'
            }

        # Collect unique residues and identify which contain metals
        residue_info = {}  # resname -> {'has_metal': bool, 'atoms': set}

        for site in redox_sites:
            for atom in site.atoms:
                resname = atom.resname
                if resname not in residue_info:
                    residue_info[resname] = {'has_metal': False, 'atoms': set()}

                residue_info[resname]['atoms'].add(atom.atom_name)

                # Check if this atom is a metal
                element = getattr(atom, 'element', None) or ''
                if element.upper() in metal_elements:
                    residue_info[resname]['has_metal'] = True

        # Categorize residues
        missing_residues = []
        metal_residues_with_params = []
        metal_residues_without_params = []
        covered_residues = []

        for resname, info in residue_info.items():
            has_params = self.has_residue(resname)

            if info['has_metal']:
                # Metal-containing residues - track whether they have FF params
                if has_params:
                    metal_residues_with_params.append(resname)
                else:
                    metal_residues_without_params.append(resname)
            elif has_params:
                covered_residues.append(resname)
            else:
                missing_residues.append(resname)

        # All metal residues (for backward compatibility)
        metal_residues = metal_residues_with_params + metal_residues_without_params

        # Include metal residues without params in missing list
        # User will be prompted to map OR proceed with Seminario+RESP
        all_missing = missing_residues + metal_residues_without_params

        # Determine overall validity
        # Valid only if ALL residues have parameters (including metal sites)
        valid = len(all_missing) == 0

        result = {
            'valid': valid,
            'missing_residues': all_missing,
            'metal_residues': metal_residues,
            'metal_residues_with_params': metal_residues_with_params,
            'metal_residues_without_params': metal_residues_without_params,
            'covered_residues': covered_residues,
            'total_residues': len(residue_info)
        }

        # Print results
        if valid:
            self.console.print(f"\n[green]✅ Force field coverage validated[/green]")
            self.console.print(f"[grey50]  Covered: {len(covered_residues)} residue(s)[/grey50]")
            if metal_residues_with_params:
                self.console.print(f"[grey50]  Metal sites (with FF params): {len(metal_residues_with_params)} residue(s)[/grey50]")
        else:
            self.console.print(f"\n[yellow]⚠️  Residues without force field parameters[/yellow]")
            if missing_residues:
                self.console.print(f"[yellow]  Non-metal missing: {', '.join(missing_residues)}[/yellow]")
            if metal_residues_without_params:
                self.console.print(f"[yellow]  Metal sites missing: {', '.join(metal_residues_without_params)}[/yellow]")
            self.console.print(f"[grey50]  Covered: {len(covered_residues)} residue(s)[/grey50]")
            if metal_residues_with_params:
                self.console.print(f"[grey50]  Metal sites (with FF params): {len(metal_residues_with_params)} residue(s)[/grey50]")

        return result

    def add_mol2_file(self, mol2_path: Path, residue_name: str) -> int:
        """
        Parse a mol2 file and add atom type definitions.

        Args:
            mol2_path: Path to mol2 file
            residue_name: Residue name to assign to atoms

        Returns:
            Number of atoms added
        """
        try:
            with open(mol2_path, 'r') as f:
                content = f.read()

            # Find @<TRIPOS>ATOM section
            atom_section_match = re.search(
                r'@<TRIPOS>ATOM\s*\n(.*?)(?:@<TRIPOS>|\Z)',
                content,
                re.DOTALL
            )

            if not atom_section_match:
                self.console.print(f"[yellow]Warning: No @<TRIPOS>ATOM section in {mol2_path}[/yellow]")
                return 0

            atom_section = atom_section_match.group(1)
            atom_count = 0

            # Parse atom lines
            # Format: atom_id atom_name x y z atom_type [subst_id [subst_name [charge [status_bit]]]]
            for line in atom_section.strip().split('\n'):
                parts = line.split()
                if len(parts) < 6:
                    continue

                atom_name = parts[1]
                atom_type = parts[5]
                atom_charge = float(parts[8]) if len(parts) > 8 else 0.0

                atom_def = AtomTypeDefinition(
                    residue_name=residue_name,
                    atom_name=atom_name,
                    atom_type=atom_type,
                    atom_charge=atom_charge
                )

                key = (residue_name, atom_name)
                self.atom_definitions[key] = atom_def
                atom_count += 1

            self.console.print(f"[green]✓ Added {atom_count} atoms from {mol2_path.name} for residue {residue_name}[/green]")
            return atom_count

        except FileNotFoundError:
            self.console.print(f"[red]Error: mol2 file not found: {mol2_path}[/red]")
            return 0
        except Exception as e:
            self.console.print(f"[red]Error parsing mol2 file {mol2_path}: {e}[/red]")
            return 0

    def add_manual_atom(self, residue_name: str, atom_name: str, atom_type: str, atom_charge: float):
        """
        Manually add an atom type definition.

        Args:
            residue_name: Residue name
            atom_name: Atom name
            atom_type: Atom type
            atom_charge: Partial charge
        """
        atom_def = AtomTypeDefinition(
            residue_name=residue_name,
            atom_name=atom_name,
            atom_type=atom_type,
            atom_charge=atom_charge
        )

        key = (residue_name, atom_name)
        self.atom_definitions[key] = atom_def

    def get_atom_types_for_element(self, element: str, include_gaff2: bool = True) -> Dict[str, List[Tuple[str, str, str, float, str]]]:
        """
        Get all atom types for a given element from the loaded force field.
        Optionally includes GAFF2 results as well.

        Useful for finding appropriate atom types for inorganic ligands.

        Args:
            element: Element symbol (e.g., 'S', 'N', 'O')
            include_gaff2: If True, also load and return GAFF2 atom types

        Returns:
            Dictionary with keys 'force_field' and optionally 'gaff2':
                'force_field': List of tuples (residue_name, atom_name, atom_type, charge, description)
                'gaff2': List of tuples (residue_name, atom_name, atom_type, charge, description)
                Note: description is empty string for force_field entries
        """
        import re

        element_upper = element.upper()

        def extract_for_element(atom_defs_dict):
            """Helper to extract atom types for element from a definitions dict."""
            results = []
            for (res_name, atom_name), atom_def in atom_defs_dict.items():
                # Extract element from atom name (usually first 1-2 characters)
                atom_element = re.match(r'([A-Z][a-z]?)', atom_name)
                if atom_element:
                    atom_elem_str = atom_element.group(1).upper()
                    if atom_elem_str == element_upper:
                        # Add empty description for force field entries
                        results.append((res_name, atom_name, atom_def.atom_type, atom_def.atom_charge, ""))
            results.sort(key=lambda x: (x[0], x[1]))
            return results

        # Get results from currently loaded force field
        ff_results = extract_for_element(self.atom_definitions)

        output = {'force_field': ff_results}

        # Optionally load GAFF2 results
        if include_gaff2:
            gaff2_results = self._get_gaff2_atom_types_for_element(element_upper)
            if gaff2_results:
                output['gaff2'] = gaff2_results

        return output

    def _get_gaff2_atom_types_for_element(self, element: str) -> List[Tuple[str, str, str, float, str]]:
        """
        Parse gaff2.dat file and extract atom types for an element.
        Does NOT modify self.atom_definitions.

        Args:
            element: Element symbol (uppercase)

        Returns:
            List of tuples: (residue_name, atom_name, atom_type, charge, description)
        """
        # GAFF2 mass-to-element mapping (all unique masses in gaff2.dat)
        MASS_TO_ELEMENT = {
            1.008: 'H',
            12.01: 'C',
            14.01: 'N',
            16.00: 'O',
            19.00: 'F',
            30.97: 'P',
            32.06: 'S',
            35.45: 'Cl',
            79.90: 'Br',
            126.9: 'I',
        }

        try:
            # Path to gaff2.dat file
            gaff2_dat = self.amberhome / 'dat' / 'leap' / 'parm' / 'gaff2.dat'

            if not gaff2_dat.exists():
                return []

            results = []

            with open(gaff2_dat, 'r') as f:
                for line in f:
                    line = line.strip()

                    # Skip header and empty lines
                    if not line or line.startswith('AMBER') or line.startswith('#'):
                        continue

                    # Parse line: atom_type  mass  radius  description
                    parts = line.split(None, 3)  # Split on whitespace, max 4 parts
                    if len(parts) < 4:
                        continue

                    atom_type = parts[0]
                    try:
                        mass = float(parts[1])
                    except ValueError:
                        continue

                    description = parts[3].strip()

                    # Map mass to element
                    elem = MASS_TO_ELEMENT.get(mass)
                    if elem and elem == element:
                        # For GAFF2, use generic "LIG" as residue name
                        results.append(("LIG", atom_type, atom_type, 0.0, description))

            return results

        except Exception as e:
            # Silently fail if gaff2.dat can't be read
            return []
