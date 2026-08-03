"""
Force Field Explorer — interactive browser for AMBER force field parameters.

Allows users to browse atom types, residues, and search for bond/angle/dihedral/VDW
parameters without needing a loaded structure. Can also type atoms from a loaded
structure against the selected force field and produce a diagnostic report.
"""

import os
from pathlib import Path
from typing import Dict, List, Optional, Callable, Set, Tuple, TYPE_CHECKING
from dataclasses import dataclass

from Bio.PDB import PDBParser
from rich.table import Table
from rich.panel import Panel

from proprep.utils.prompts import prompt_with_context, confirm_with_context
from proprep.forcefield_prep.forcefield_data import ForceFieldData

from .processor_command import MenuCommand

if TYPE_CHECKING:
    from .processor import Processor


# Standard amino acid residue names (for terminal detection)
STANDARD_AMINO_ACIDS = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'HID', 'HIE', 'HIP', 'CYX', 'CYM',
}

# PDB atom name aliases for N-terminal residues
NTERM_ATOM_ALIASES = {
    "H": "H1",  # PDB uses H for first N-H; AMBER uses H1
}


@dataclass
class TypedAtom:
    """Result of typing a single atom against the force field."""
    chain_id: str
    resid: int
    insertion_code: str
    resname: str         # original PDB residue name
    resolved_resname: str  # FF lookup name (e.g., NALA, HIE)
    atom_name: str
    x: float
    y: float
    z: float
    atom_type: Optional[str]   # AMBER atom type, or None if untyped
    charge: Optional[float]    # partial charge, or None if untyped
    source: str                # "force_field", "terminal_variant", "untyped"


PAGE_SIZE = 20


class ForceFieldExplorerCommand(MenuCommand):
    """Interactive force field parameter browser."""

    def __init__(self, processor: 'Processor'):
        super().__init__(processor, "Force Field Explorer", "Browse force field parameters")
        self.ff_data: Optional[ForceFieldData] = None

    def execute(self) -> None:
        self.enter_menu()
        try:
            if not self._load_forcefield():
                return
            self._explorer_loop()
        finally:
            self.exit_menu()

    # =========================================================================
    # Force field loading
    # =========================================================================

    def _load_forcefield(self) -> bool:
        """Load a force field. Returns True on success."""
        self.ff_data = ForceFieldData(processor=self.processor)
        try:
            n_atoms = self.ff_data.select_and_load(self.console, processor=self.processor)
        except RuntimeError as e:
            self.console.print(f"[red]Error: {e}[/red]")
            self.console.print(
                "[yellow]Please set AMBERHOME to your AmberTools installation directory.[/yellow]"
            )
            return False
        except FileNotFoundError as e:
            self.console.print(f"[red]Error: {e}[/red]")
            return False

        if n_atoms == 0:
            self.console.print("[yellow]No atom definitions loaded.[/yellow]")
            return False

        stats = self.ff_data.get_statistics()
        self.console.print(
            f"\n[green]Loaded: {stats['n_atom_definitions']} atoms, "
            f"{stats['n_residues']} residues, "
            f"{stats['n_bond_params']} bonds, "
            f"{stats['n_angle_params']} angles, "
            f"{stats['n_dihedral_params']} dihedrals, "
            f"{stats['n_nonbonded_params']} VDW params[/green]"
        )
        return True

    # =========================================================================
    # Main explorer loop
    # =========================================================================

    def _explorer_loop(self) -> None:
        while True:
            self.show_breadcrumbs()
            self._show_explorer_menu()

            choice = prompt_with_context(
                processor=self.processor,
                prompt="\nEnter your choice",
                choices=list(self.options_map.keys()),
                default="1",
                module="Force Field Explorer",
                description="Select explorer action",
                options_map=self.options_map,
            )

            if choice not in self.options_map:
                self.console.print(f"[red]Invalid choice: {choice}[/red]")
                continue

            action = self.options_map[choice]
            if action == "back":
                break
            elif action == "overview":
                self._show_overview()
            elif action == "atom_types":
                self._browse_atom_types()
            elif action == "residues":
                self._browse_residues()
            elif action == "search":
                self._search_parameters()
            elif action == "reverse":
                self._reverse_lookup()
            elif action == "reload":
                if not self._load_forcefield():
                    break
            elif action == "type_structure":
                self._type_structure_atoms()

    def _show_explorer_menu(self) -> None:
        leaprc_names = ", ".join(p.name for p in self.ff_data.loaded_leaprcs)
        self.console.print(f"\n[bold]===== Force Field Explorer =====[/bold]")
        self.console.print(f"[grey50]Loaded: {leaprc_names}[/grey50]\n")

        self.options_map = {
            "1": "overview",
            "2": "atom_types",
            "3": "residues",
            "4": "search",
            "5": "reverse",
            "6": "reload",
            "7": "type_structure",
            "b": "back",
        }

        self.console.print("  1. Overview (statistics & loaded files)")
        self.console.print("  2. Browse Atom Types")
        self.console.print("  3. Browse Residues")
        self.console.print("  4. Search Parameters (bond/angle/dihedral/VDW/mass)")
        self.console.print("  5. Reverse Lookup (all params for an atom type)")
        self.console.print("  6. Load Different Force Field")
        self.console.print("  7. Type Structure Atoms (atom typing report)")
        self.console.print("  b. Back")

    # =========================================================================
    # 1. Overview
    # =========================================================================

    def _show_overview(self) -> None:
        stats = self.ff_data.get_statistics()

        table = Table(title="Force Field Statistics", show_header=True, header_style="bold cyan")
        table.add_column("Category", style="cyan")
        table.add_column("Count", style="green", justify="right")

        table.add_row("Atom definitions", str(stats['n_atom_definitions']))
        table.add_row("Residues", str(stats['n_residues']))
        table.add_row("Residue aliases", str(stats['n_aliases']))
        table.add_row("Mass parameters", str(stats['n_mass_params']))
        table.add_row("Bond parameters", str(stats['n_bond_params']))
        table.add_row("Angle parameters", str(stats['n_angle_params']))
        table.add_row("Dihedral parameters", str(stats['n_dihedral_params']))
        table.add_row("VDW parameters", str(stats['n_nonbonded_params']))
        table.add_row("Improper parameters", str(len(self.ff_data.improper_parameters)))

        self.console.print(table)

        # Loaded files
        self.console.print("\n[bold]Loaded Files:[/bold]")
        for label, paths in [
            ("Leaprc", self.ff_data.loaded_leaprcs),
            ("Lib/Off", self.ff_data.loaded_libs),
            ("Parm.dat", self.ff_data.loaded_parms),
            ("Frcmod", self.ff_data.loaded_frcmods),
            ("Mol2", self.ff_data.loaded_mol2s),
        ]:
            if paths:
                self.console.print(f"  [cyan]{label}:[/cyan]")
                for p in paths:
                    self.console.print(f"    {p.name}")

        self.console.print()

    # =========================================================================
    # 2. Browse Atom Types
    # =========================================================================

    def _browse_atom_types(self) -> None:
        while True:
            filter_str = prompt_with_context(
                processor=self.processor,
                prompt="Enter atom type filter (e.g., 'CT'), Enter for all, 'b' to go back",
                default="",
                module="Force Field Explorer",
                description="Filter atom types",
            )

            if filter_str.lower() == "b":
                return

            # Collect all unique atom types from mass + nonbonded
            all_types = sorted(
                set(self.ff_data.mass_parameters.keys())
                | set(self.ff_data.nonbonded_parameters.keys())
            )

            if filter_str:
                filt = filter_str.upper()
                all_types = [t for t in all_types if filt in t.upper()]

            if not all_types:
                self.console.print("[yellow]No atom types match the filter.[/yellow]")
                continue

            def build_row(atom_type: str) -> List[str]:
                mass_p = self.ff_data.mass_parameters.get(atom_type)
                vdw_p = self.ff_data.nonbonded_parameters.get(atom_type)
                return [
                    atom_type,
                    f"{mass_p.mass:.3f}" if mass_p else "-",
                    f"{vdw_p.radius:.4f}" if vdw_p else "-",
                    f"{vdw_p.well_depth:.4f}" if vdw_p else "-",
                    (mass_p.comment if mass_p and mass_p.comment else
                     vdw_p.comment if vdw_p and vdw_p.comment else ""),
                ]

            columns = [
                ("Type", "cyan"),
                ("Mass", "green"),
                ("VDW Radius", "green"),
                ("VDW Depth", "green"),
                ("Description", "dim"),
            ]

            self._display_paginated(
                all_types, columns, build_row,
                f"Atom Types ({len(all_types)} total)"
            )

    # =========================================================================
    # 3. Browse Residues
    # =========================================================================

    def _browse_residues(self) -> None:
        while True:
            filter_str = prompt_with_context(
                processor=self.processor,
                prompt="Enter residue filter (e.g., 'ALA'), Enter for all, 'b' to go back",
                default="",
                module="Force Field Explorer",
                description="Filter residues",
            )

            if filter_str.lower() == "b":
                return

            residues = self.ff_data.get_available_residues()

            if filter_str:
                filt = filter_str.upper()
                residues = [r for r in residues if filt in r.upper()]

            if not residues:
                self.console.print("[yellow]No residues match the filter.[/yellow]")
                continue

            def build_row(resname: str) -> List[str]:
                atoms = self.ff_data.get_residue_atoms(resname)
                source = self.ff_data.get_source_file(resname)
                return [resname, str(len(atoms)), source]

            columns = [
                ("Residue", "cyan"),
                ("Atoms", "green"),
                ("Source", "dim"),
            ]

            selected = self._display_paginated(
                residues, columns, build_row,
                f"Residues ({len(residues)} total)",
                selectable=True,
            )

            if selected is not None and 0 <= selected < len(residues):
                self._show_residue_detail(residues[selected])

    def _show_residue_detail(self, resname: str) -> None:
        atoms = self.ff_data.get_residue_atoms(resname)
        if not atoms:
            self.console.print(f"[yellow]No atoms found for residue {resname}.[/yellow]")
            return

        table = Table(
            title=f"Residue: {resname} ({len(atoms)} atoms)",
            show_header=True,
            header_style="bold cyan",
        )
        table.add_column("Atom Name", style="cyan")
        table.add_column("Type", style="green")
        table.add_column("Charge", style="yellow", justify="right")
        table.add_column("Source", style="grey50")

        for atom_def in sorted(atoms, key=lambda a: a.atom_name):
            table.add_row(
                atom_def.atom_name,
                atom_def.atom_type,
                f"{atom_def.charge:.4f}" if atom_def.charge is not None else "-",
                atom_def.source_file,
            )

        self.console.print(table)
        self.console.print()

    # =========================================================================
    # 4. Search Parameters
    # =========================================================================

    def _search_parameters(self) -> None:
        while True:
            self.console.print("\n[bold]Search Parameters[/bold]")
            self.console.print("  a. Bond (2 atom types)")
            self.console.print("  b. Angle (3 atom types)")
            self.console.print("  c. Dihedral (4 atom types)")
            self.console.print("  d. VDW / Nonbonded (1 atom type)")
            self.console.print("  e. Mass (1 atom type)")
            self.console.print("  q. Back")

            search_map = {"a": "bond", "b": "angle", "c": "dihedral", "d": "vdw", "e": "mass", "q": "back"}

            choice = prompt_with_context(
                processor=self.processor,
                prompt="Select parameter type",
                choices=list(search_map.keys()),
                default="a",
                module="Force Field Explorer",
                description="Select parameter type to search",
                options_map=search_map,
            )

            action = search_map.get(choice, "back")
            if action == "back":
                return

            if action == "bond":
                self._search_bond()
            elif action == "angle":
                self._search_angle()
            elif action == "dihedral":
                self._search_dihedral()
            elif action == "vdw":
                self._search_vdw()
            elif action == "mass":
                self._search_mass()

    def _search_bond(self) -> None:
        types_str = prompt_with_context(
            processor=self.processor,
            prompt="Enter 2 atom types (e.g., 'CT N' or 'CT-N')",
            default="",
            module="Force Field Explorer",
            description="Bond parameter search",
        )
        types = self._parse_atom_types(types_str, 2)
        if not types:
            return

        param = self.ff_data.get_bond_parameter(types[0], types[1])
        if param:
            table = Table(title=f"Bond: {types[0]}-{types[1]}", header_style="bold cyan")
            table.add_column("Type 1", style="cyan")
            table.add_column("Type 2", style="cyan")
            table.add_column("Force Const.", style="green", justify="right")
            table.add_column("Eq Length", style="green", justify="right")
            table.add_column("Source", style="grey50")

            table.add_row(
                param.type1, param.type2,
                f"{param.force_constant:.1f}", f"{param.eq_length:.4f}",
                param.source,
            )
            self.console.print(table)
        else:
            self.console.print(f"[yellow]No bond parameter found for {types[0]}-{types[1]}[/yellow]")

    def _search_angle(self) -> None:
        types_str = prompt_with_context(
            processor=self.processor,
            prompt="Enter 3 atom types (e.g., 'CT CT N')",
            default="",
            module="Force Field Explorer",
            description="Angle parameter search",
        )
        types = self._parse_atom_types(types_str, 3)
        if not types:
            return

        param = self.ff_data.get_angle_parameter(types[0], types[1], types[2])
        if param:
            table = Table(title=f"Angle: {types[0]}-{types[1]}-{types[2]}", header_style="bold cyan")
            table.add_column("Type 1", style="cyan")
            table.add_column("Type 2 (center)", style="cyan")
            table.add_column("Type 3", style="cyan")
            table.add_column("Force Const.", style="green", justify="right")
            table.add_column("Eq Angle", style="green", justify="right")
            table.add_column("Source", style="grey50")

            table.add_row(
                param.type1, param.type2, param.type3,
                f"{param.force_constant:.1f}", f"{param.eq_angle:.2f}",
                param.source,
            )
            self.console.print(table)
        else:
            self.console.print(
                f"[yellow]No angle parameter found for {types[0]}-{types[1]}-{types[2]}[/yellow]"
            )

    def _search_dihedral(self) -> None:
        types_str = prompt_with_context(
            processor=self.processor,
            prompt="Enter 4 atom types (e.g., 'X CT CT X')",
            default="",
            module="Force Field Explorer",
            description="Dihedral parameter search",
        )
        types = self._parse_atom_types(types_str, 4)
        if not types:
            return

        params = self.ff_data.get_dihedral_parameters(types[0], types[1], types[2], types[3])
        if params:
            table = Table(
                title=f"Dihedral: {'-'.join(types)} ({len(params)} term{'s' if len(params) != 1 else ''})",
                header_style="bold cyan",
            )
            table.add_column("Types", style="cyan")
            table.add_column("IDIVF", style="green", justify="right")
            table.add_column("Barrier", style="green", justify="right")
            table.add_column("Phase", style="green", justify="right")
            table.add_column("Period.", style="green", justify="right")
            table.add_column("Source", style="grey50")

            for p in params:
                table.add_row(
                    f"{p.type1}-{p.type2}-{p.type3}-{p.type4}",
                    str(p.divider),
                    f"{p.barrier:.3f}",
                    f"{p.phase:.1f}",
                    f"{abs(p.periodicity):.1f}",
                    p.source,
                )
            self.console.print(table)
        else:
            self.console.print(
                f"[yellow]No dihedral parameter found for {'-'.join(types)}[/yellow]"
            )

    def _search_vdw(self) -> None:
        atom_type = prompt_with_context(
            processor=self.processor,
            prompt="Enter atom type (e.g., 'CT')",
            default="",
            module="Force Field Explorer",
            description="VDW parameter search",
        ).strip()

        if not atom_type:
            return

        param = self.ff_data.get_nonbonded_parameter(atom_type)
        if param:
            table = Table(title=f"VDW: {atom_type}", header_style="bold cyan")
            table.add_column("Type", style="cyan")
            table.add_column("Radius (Rmin/2)", style="green", justify="right")
            table.add_column("Well Depth", style="green", justify="right")
            table.add_column("Source", style="grey50")

            table.add_row(
                param.atom_type,
                f"{param.radius:.4f}",
                f"{param.well_depth:.4f}",
                param.source,
            )
            self.console.print(table)
        else:
            self.console.print(f"[yellow]No VDW parameter found for '{atom_type}'[/yellow]")

    def _search_mass(self) -> None:
        atom_type = prompt_with_context(
            processor=self.processor,
            prompt="Enter atom type (e.g., 'CT')",
            default="",
            module="Force Field Explorer",
            description="Mass parameter search",
        ).strip()

        if not atom_type:
            return

        param = self.ff_data.get_mass_parameter(atom_type)
        if param:
            table = Table(title=f"Mass: {atom_type}", header_style="bold cyan")
            table.add_column("Type", style="cyan")
            table.add_column("Mass", style="green", justify="right")
            table.add_column("Polarizability", style="green", justify="right")
            table.add_column("Description", style="grey50")
            table.add_column("Source", style="grey50")

            table.add_row(
                param.atom_type,
                f"{param.mass:.3f}",
                f"{param.polarizability:.4f}",
                param.comment,
                param.source,
            )
            self.console.print(table)
        else:
            self.console.print(f"[yellow]No mass parameter found for '{atom_type}'[/yellow]")

    # =========================================================================
    # 5. Reverse Lookup
    # =========================================================================

    def _reverse_lookup(self) -> None:
        atom_type = prompt_with_context(
            processor=self.processor,
            prompt="Enter atom type to search for (e.g., 'CT')",
            default="",
            module="Force Field Explorer",
            description="Reverse parameter lookup",
        ).strip()

        if not atom_type:
            return

        at = atom_type

        # Gather all results up front
        mass_p = self.ff_data.mass_parameters.get(at)
        vdw_p = self.ff_data.nonbonded_parameters.get(at)

        bonds = sorted(
            [p for p in self.ff_data.bond_parameters.values()
             if at in (p.type1.strip(), p.type2.strip())],
            key=lambda b: f"{b.type1}-{b.type2}",
        )

        angles = [p for p in self.ff_data.angle_parameters.values()
                  if at in (p.type1.strip(), p.type2.strip(), p.type3.strip())]

        dihedrals = []
        for params_list in self.ff_data.dihedral_parameters.values():
            for p in params_list:
                if at in (p.type1.strip(), p.type2.strip(), p.type3.strip(), p.type4.strip()):
                    dihedrals.append(p)
                    break

        impropers = [p for p in self.ff_data.improper_parameters
                     if at in (p.type1.strip(), p.type2.strip(), p.type3.strip(), p.type4.strip())]

        # Check if anything was found
        total = ((1 if mass_p else 0) + (1 if vdw_p else 0)
                 + len(bonds) + len(angles) + len(dihedrals) + len(impropers))
        if total == 0:
            self.console.print(f"[yellow]No parameters found for atom type '{at}'[/yellow]")
            return

        # Show summary with mass/VDW inline, then let user pick categories
        self.console.print(f"\n[bold]Parameters for atom type: {at}[/bold]")
        if mass_p:
            self.console.print(
                f"  Mass: {mass_p.mass:.3f}  "
                f"[grey50]({mass_p.comment}, {mass_p.source})[/grey50]"
            )
        if vdw_p:
            self.console.print(
                f"  VDW:  radius={vdw_p.radius:.4f}  "
                f"depth={vdw_p.well_depth:.4f}  [grey50]({vdw_p.source})[/grey50]"
            )

        # Build category menu
        categories = {}
        cat_options = {}
        opt = 1
        if bonds:
            categories[str(opt)] = ("bonds", bonds)
            cat_options[str(opt)] = f"Bonds ({len(bonds)})"
            opt += 1
        if angles:
            categories[str(opt)] = ("angles", angles)
            cat_options[str(opt)] = f"Angles ({len(angles)})"
            opt += 1
        if dihedrals:
            categories[str(opt)] = ("dihedrals", dihedrals)
            cat_options[str(opt)] = f"Dihedrals ({len(dihedrals)})"
            opt += 1
        if impropers:
            categories[str(opt)] = ("impropers", impropers)
            cat_options[str(opt)] = f"Impropers ({len(impropers)})"
            opt += 1

        if not categories:
            return

        cat_options["b"] = "Back"

        while True:
            self.console.print()
            for key, label in cat_options.items():
                if key != "b":
                    self.console.print(f"  {key}. {label}")
            self.console.print(f"  b. Back")

            choice = prompt_with_context(
                processor=self.processor,
                prompt=f"View category for {at}",
                choices=list(cat_options.keys()),
                default="b",
                module="Force Field Explorer",
                description="Reverse lookup category selection",
                options_map=cat_options,
            )

            if choice == "b" or choice not in categories:
                return

            cat_name, cat_items = categories[choice]
            self._show_reverse_category(at, cat_name, cat_items)

    def _show_reverse_category(self, at: str, category: str, items: list) -> None:
        """Display a single reverse-lookup category."""
        if category == "bonds":
            self._display_paginated(
                items,
                [("Types", "cyan"), ("Force Const.", "green"),
                 ("Eq Length", "green"), ("Source", "dim")],
                lambda p: [f"{p.type1}-{p.type2}", f"{p.force_constant:.1f}",
                           f"{p.eq_length:.4f}", p.source],
                f"Bonds involving {at} ({len(items)} found)",
            )
        elif category == "angles":
            self._display_paginated(
                items,
                [("Types", "cyan"), ("Force Const.", "green"),
                 ("Eq Angle", "green"), ("Source", "dim")],
                lambda p: [f"{p.type1}-{p.type2}-{p.type3}", f"{p.force_constant:.1f}",
                           f"{p.eq_angle:.2f}", p.source],
                f"Angles involving {at} ({len(items)} found)",
            )
        elif category == "dihedrals":
            self._display_paginated(
                items,
                [("Types", "cyan"), ("IDIVF", "green"), ("Barrier", "green"),
                 ("Phase", "green"), ("Period.", "green"), ("Source", "dim")],
                lambda p: [f"{p.type1}-{p.type2}-{p.type3}-{p.type4}", str(p.divider),
                           f"{p.barrier:.3f}", f"{p.phase:.1f}",
                           f"{abs(p.periodicity):.1f}", p.source],
                f"Dihedrals involving {at} ({len(items)} found)",
            )
        elif category == "impropers":
            self._display_paginated(
                items,
                [("Types", "cyan"), ("Barrier", "green"), ("Phase", "green"),
                 ("Period.", "green"), ("Source", "dim")],
                lambda p: [f"{p.type1}-{p.type2}-{p.type3}-{p.type4}", f"{p.barrier:.3f}",
                           f"{p.phase:.1f}", str(p.periodicity), p.source],
                f"Impropers involving {at} ({len(items)} found)",
            )

    # =========================================================================
    # 7. Type Structure Atoms
    # =========================================================================

    def _type_structure_atoms(self) -> None:
        """Type atoms from a loaded structure against the current force field."""

        # --- Get PDB file ---
        pdb_file = self._get_structure_pdb()
        if not pdb_file:
            return

        # --- Parse structure with BioPython ---
        parser = PDBParser(QUIET=True)
        try:
            structure = parser.get_structure('structure', pdb_file)
        except Exception as e:
            self.console.print(f"[red]Failed to parse PDB file: {e}[/red]")
            return

        # --- Detect terminals ---
        terminal_map = self._detect_terminals(structure)

        # --- Gather per-residue atom lists for HIS/CYS resolution ---
        residue_atoms: Dict[Tuple[str, int, str, str], list] = {}
        for model in structure:
            for chain in model:
                for residue in chain:
                    het, resseq, icode = residue.get_id()
                    if het.strip():
                        continue  # skip HETATM water/ions unless they have residue names
                    resname = residue.get_resname().strip()
                    chain_id = chain.get_id().strip()
                    icode_s = icode.strip()
                    key = (chain_id, resseq, icode_s, resname)
                    if key not in residue_atoms:
                        residue_atoms[key] = []
                    for atom in residue:
                        residue_atoms[key].append(atom)

        # --- Pre-resolve ambiguous residues (HIS, CYS) ---
        resname_cache: Dict[Tuple[str, int, str], str] = {}
        for (chain_id, resseq, icode_s, resname), atoms in residue_atoms.items():
            res_key = (chain_id, resseq, icode_s)
            resolved = self._resolve_structure_resname(resname, atoms)
            resname_cache[res_key] = resolved

        # --- Type all atoms ---
        typed_atoms: List[TypedAtom] = []

        for model in structure:
            for chain in model:
                for residue in chain:
                    het, resseq, icode = residue.get_id()
                    resname = residue.get_resname().strip()
                    chain_id = chain.get_id().strip()
                    icode_s = icode.strip()
                    res_key = (chain_id, resseq, icode_s)

                    resolved = resname_cache.get(res_key, resname)

                    # Apply terminal prefix
                    terminal_type = terminal_map.get(res_key)
                    if terminal_type == "NTERM":
                        lookup_resname = f"N{resolved}"
                    elif terminal_type == "CTERM":
                        lookup_resname = f"C{resolved}"
                    else:
                        lookup_resname = resolved

                    for atom in residue:
                        coords = atom.get_coord()
                        atom_name = atom.get_name()

                        atom_def = self.ff_data.get_atom_definition(lookup_resname, atom_name)

                        # N-terminal H → H1 fallback
                        if not atom_def and terminal_type == "NTERM":
                            alias = NTERM_ATOM_ALIASES.get(atom_name)
                            if alias:
                                atom_def = self.ff_data.get_atom_definition(lookup_resname, alias)

                        # Fallback: try without terminal prefix
                        if not atom_def and lookup_resname != resolved:
                            atom_def = self.ff_data.get_atom_definition(resolved, atom_name)

                        # Fallback: try original PDB resname
                        if not atom_def and resolved != resname:
                            atom_def = self.ff_data.get_atom_definition(resname, atom_name)

                        if atom_def:
                            typed_atoms.append(TypedAtom(
                                chain_id=chain_id, resid=resseq, insertion_code=icode_s,
                                resname=resname, resolved_resname=lookup_resname,
                                atom_name=atom_name,
                                x=float(coords[0]), y=float(coords[1]), z=float(coords[2]),
                                atom_type=atom_def.atom_type,
                                charge=atom_def.charge,
                                source="force_field",
                            ))
                        else:
                            typed_atoms.append(TypedAtom(
                                chain_id=chain_id, resid=resseq, insertion_code=icode_s,
                                resname=resname, resolved_resname=lookup_resname,
                                atom_name=atom_name,
                                x=float(coords[0]), y=float(coords[1]), z=float(coords[2]),
                                atom_type=None, charge=None, source="untyped",
                            ))

        # --- Display results ---
        self._display_typing_report(typed_atoms, pdb_file)

    def _get_structure_pdb(self) -> Optional[str]:
        """Get PDB file path from workspace or user input."""
        pdb_file = None

        # Try workspace first
        try:
            ws = self.workspace
            if ws:
                from proprep.utils.structure_selector import StructureSelector
                selector = StructureSelector(ws, self.console, processor=self.processor)
                pdb_file = selector.get_structure(silent=True)
        except Exception:
            pass

        if pdb_file and os.path.isfile(pdb_file):
            self.console.print(f"\n[green]Using workspace structure:[/green] {Path(pdb_file).name}")
            use_it = confirm_with_context(
                processor=self.processor,
                prompt="Use this structure?",
                default=True,
                module="Force Field Explorer",
                description="Confirm workspace structure for atom typing",
            )
            if use_it:
                return pdb_file

        # File browser fallback
        from proprep.structure_prep.pdb_loader import display_pdb_file_menu
        return display_pdb_file_menu(
            directory=".", console=self.console, processor=self.processor
        )

    def _detect_terminals(self, structure) -> Dict[Tuple[str, int, str], str]:
        """
        Lightweight terminal detection from PDB structure.

        Returns dict mapping (chain_id, resid, insertion_code) -> "NTERM" | "CTERM".
        Only standard amino acids are considered.
        """
        terminal_map: Dict[Tuple[str, int, str], str] = {}

        for model in structure:
            for chain in model:
                # Collect protein residues in sequence order
                protein_residues = []
                for residue in chain:
                    het, resseq, icode = residue.get_id()
                    resname = residue.get_resname().strip()
                    if resname in STANDARD_AMINO_ACIDS:
                        protein_residues.append((resseq, icode.strip(), resname, residue))

                if len(protein_residues) < 2:
                    continue

                # First and last protein residues are terminals
                first = protein_residues[0]
                last = protein_residues[-1]
                chain_id = chain.get_id().strip()

                # Confirm N-terminal: check for H1/H2/H3 or single H
                first_atoms = {a.get_name() for a in first[3]}
                if "H1" in first_atoms or "H2" in first_atoms or "H3" in first_atoms or "H" in first_atoms:
                    terminal_map[(chain_id, first[0], first[1])] = "NTERM"

                # Confirm C-terminal: check for OXT
                last_atoms = {a.get_name() for a in last[3]}
                if "OXT" in last_atoms:
                    terminal_map[(chain_id, last[0], last[1])] = "CTERM"

        return terminal_map

    def _resolve_structure_resname(self, resname: str, atoms) -> str:
        """
        Resolve ambiguous residue names (HIS, CYS) based on atom content.

        Args:
            resname: PDB residue name
            atoms: list of BioPython Atom objects for this residue
        """
        atom_names = {a.get_name() for a in atoms}

        if resname == "HIS":
            has_hd1 = "HD1" in atom_names
            has_he2 = "HE2" in atom_names
            if has_hd1 and has_he2:
                return "HIP"
            elif has_hd1:
                return "HID"
            elif has_he2:
                return "HIE"
            # No diagnostic hydrogens — default to HIE (most common)
            return "HIE"

        if resname == "CYS":
            # No HG (thiol H) means disulfide → CYX
            if "HG" not in atom_names and "SG" in atom_names:
                return "CYX"

        return resname

    def _display_typing_report(self, typed_atoms: List[TypedAtom], pdb_file: str) -> None:
        """Display the atom typing report with summary and paginated detail."""
        if not typed_atoms:
            self.console.print("[yellow]No atoms found in structure.[/yellow]")
            return

        typed = [a for a in typed_atoms if a.atom_type is not None]
        untyped = [a for a in typed_atoms if a.atom_type is None]
        unique_types = {a.atom_type for a in typed}

        # Halo untyped atoms in the viewer so the diagnostic question
        # the report is trying to answer ("where are my untyped atoms?")
        # has a visible answer in the structure pane. The halo persists
        # through the interactive detail view and is cleared when the
        # user picks 'done'.
        self._highlight_untyped_atoms(untyped, pdb_file)

        # --- Summary ---
        self.console.print(f"\n[bold cyan]{'=' * 65}[/bold cyan]")
        self.console.print(f"[bold cyan]  Atom Typing Report: {Path(pdb_file).name}[/bold cyan]")
        self.console.print(f"[bold cyan]{'=' * 65}[/bold cyan]")

        summary = Table(show_header=False, box=None, padding=(0, 2))
        summary.add_column("Label", style="bold")
        summary.add_column("Value", style="green")
        summary.add_row("Total atoms", str(len(typed_atoms)))
        summary.add_row("Typed", f"{len(typed)} ({100 * len(typed) / len(typed_atoms):.1f}%)")
        summary.add_row("Untyped", f"{len(untyped)} ({100 * len(untyped) / len(typed_atoms):.1f}%)")
        summary.add_row("Unique atom types", str(len(unique_types)))

        if typed:
            total_charge = sum(a.charge for a in typed)
            summary.add_row("Total charge (typed)", f"{total_charge:+.4f}")

        self.console.print(summary)
        self.console.print()

        # --- Untyped residue summary ---
        if untyped:
            untyped_residues: Dict[str, int] = {}
            for a in untyped:
                key = f"{a.chain_id}:{a.resname}{a.resid}"
                untyped_residues[key] = untyped_residues.get(key, 0) + 1

            self.console.print("[bold yellow]Untyped residues:[/bold yellow]")
            for res_label, count in sorted(untyped_residues.items()):
                # highlight=False: the label is data, not markup. Rich's repr
                # highlighter reads chain:resname labels like "A:FE1" / "B:FES501"
                # as IPv6 addresses (hex letters + colon) and paints them bold
                # bright-green — and untyped residues are exactly the metals /
                # Fe-S clusters / cofactors that trip it. Disable auto-highlight.
                self.console.print(f"  {res_label} ({count} atoms)", highlight=False)
            self.console.print()

        # --- Per-residue charge summary ---
        self._display_residue_charge_summary(typed_atoms)

        # --- Interactive detail view ---
        while True:
            detail_choice = prompt_with_context(
                processor=self.processor,
                prompt="View \\[a]ll atoms, \\[t]yped only, \\[u]ntyped only, \\[e]xport report, \\[d]one",
                choices=["a", "t", "u", "e", "d"],
                default="d",
                module="Force Field Explorer",
                description="Atom typing report view",
            )

            if detail_choice == "d":
                # Clean up the untyped-atoms halo so the next module to
                # visit the viewer doesn't inherit a stale FF-Explorer
                # overlay.
                from proprep.structure_prep.viewer_coordinator import viewer as _viewer
                _viewer.unhighlight("ff_explorer_untyped")
                return
            elif detail_choice == "e":
                self._export_typing_report(typed_atoms, pdb_file)
            else:
                if detail_choice == "t":
                    display_atoms = typed
                    label = "Typed Atoms"
                elif detail_choice == "u":
                    display_atoms = untyped
                    label = "Untyped Atoms"
                else:
                    display_atoms = typed_atoms
                    label = "All Atoms"

                if not display_atoms:
                    self.console.print("[yellow]No atoms to display.[/yellow]")
                    continue

                self._display_atom_table(display_atoms, label)

    def _highlight_untyped_atoms(self, untyped: List[TypedAtom], pdb_file: str) -> None:
        """Halo the untyped atoms in the viewer for the typing-report diagnostic.

        Builds a per-residue NGL selection grouping atom names so the
        emitted clause stays compact even for cofactors with dozens of
        untyped atoms (e.g. an unparameterized HEM with ~70 atoms).
        Uses the same halo style as the redox detector's bonded-atom
        rep — translucent yellow spacefill — so visually the user
        recognises "this is a question that needs attention".
        """
        if not untyped:
            return

        from proprep.structure_prep.viewer_coordinator import viewer as _viewer
        _viewer.show_structure(pdb_file)

        by_residue: Dict[Tuple[str, int], set] = {}
        for atom in untyped:
            key = (atom.chain_id, atom.resid)
            by_residue.setdefault(key, set()).add(atom.atom_name)

        clauses = []
        for (chain, resid), names in by_residue.items():
            sorted_names = sorted(n for n in names if n)
            if not sorted_names:
                continue
            names_clause = " or ".join(f".{n}" for n in sorted_names)
            clauses.append(f"(:{chain} and {resid} and ({names_clause}))")

        if not clauses:
            return

        _viewer.highlight(
            " or ".join(clauses),
            style="halo",
            color="#ffff00",
            label="ff_explorer_untyped",
        )

    def _display_residue_charge_summary(self, typed_atoms: List[TypedAtom]) -> None:
        """Show per-residue charge summary.

        "Sum of charges" is the literal sum of partial charges over the atoms
        that are present and typed in this residue. It equals AMBER's formal
        integer residue charge only when the residue is chemically complete
        (all atoms, including hydrogens, present and typed); for a heavy-atom-
        only structure the sums are non-integer because the H charges are absent.
        """
        # label -> [charge, n_typed, n_total, resname, chain_id, resid]
        residue_charges: Dict[str, list] = {}

        for a in typed_atoms:
            label = f"{a.chain_id}:{a.resname}{a.resid}"
            if label not in residue_charges:
                residue_charges[label] = [0.0, 0, 0, a.resname, a.chain_id, a.resid]
            entry = residue_charges[label]
            entry[2] += 1  # total
            if a.charge is not None:
                entry[0] += a.charge
                entry[1] += 1

        # Sort by (chain, residue number) so e.g. ALA10 < ALA11 < ALA107,
        # not the lexicographic ALA10 < ALA107 < ALA11.
        items = sorted(
            residue_charges.items(),
            key=lambda kv: (kv[1][4], kv[1][5], kv[1][3]),
        )

        def build_row(item):
            label, (charge, n_typed, n_total, resname, _chain, _resid) = item
            coverage = f"{n_typed}/{n_total}"
            charge_str = f"{charge:+.4f}" if n_typed > 0 else "-"
            status = "complete" if n_typed == n_total else "partial" if n_typed > 0 else "none"
            style = "" if status == "complete" else "[yellow]" if status == "partial" else "[red]"
            end = "" if status == "complete" else "[/yellow]" if status == "partial" else "[/red]"
            return [f"{style}{label}{end}", coverage, charge_str]

        columns = [
            ("Residue", "cyan"),
            ("Coverage (typed/present)", "green"),
            ("Sum of charges", "yellow"),
        ]

        # Context legend: total charge over everything typed, plus a note on what
        # the per-residue number means and why it may not be a whole number.
        total_charge = sum(e[0] for e in residue_charges.values())
        total_typed = sum(e[1] for e in residue_charges.values())
        total_atoms = sum(e[2] for e in residue_charges.values())
        self.console.print(
            f"[grey50]System total: {total_charge:+.4f} e over {total_typed}/{total_atoms} "
            f"atoms typed.[/grey50]"
        )
        self.console.print(
            "[grey50]\"Sum of charges\" = sum of partial charges over the atoms typed in "
            "each residue.\n"
            "It matches AMBER's formal integer residue charge only when the residue is "
            "chemically complete\n(all atoms incl. hydrogens present); heavy-atom-only "
            "structures give non-integer sums.[/grey50]"
        )

        self._display_paginated(
            items, columns, build_row,
            f"Per-Residue Charge Summary ({len(items)} residues)",
        )

    def _display_atom_table(self, atoms: List[TypedAtom], title: str) -> None:
        """Display paginated atom detail table."""

        def build_row(a: TypedAtom) -> List[str]:
            res_label = f"{a.chain_id}:{a.resname}{a.resid}"
            resolved = a.resolved_resname if a.resolved_resname != a.resname else ""
            atype = a.atom_type if a.atom_type else "[red]???[/red]"
            charge = f"{a.charge:+.4f}" if a.charge is not None else "[red]-[/red]"
            coords = f"{a.x:8.3f} {a.y:8.3f} {a.z:8.3f}"
            return [res_label, resolved, a.atom_name, atype, charge, coords]

        columns = [
            ("Residue", "cyan"),
            ("FF Name", "dim"),
            ("Atom", "cyan"),
            ("Type", "green"),
            ("Charge", "yellow"),
            ("X        Y        Z", "dim"),
        ]

        self._display_paginated(atoms, columns, build_row, title)

    def _export_typing_report(self, typed_atoms: List[TypedAtom], pdb_file: str) -> None:
        """Export atom typing report as a fixed-width aligned text file."""
        default_name = Path(pdb_file).stem + "_atom_typing.txt"

        out_path = prompt_with_context(
            processor=self.processor,
            prompt="Export path",
            default=default_name,
            module="Force Field Explorer",
            description="Export atom typing report file path",
        ).strip()

        if not out_path:
            return

        out_path = os.path.expanduser(out_path)

        # Column widths — generous enough for real data
        # Chain ResID ResName FF_Name AtomName Type   Charge       X        Y        Z
        header = (
            f"{'Chain':<6}{'ResID':>6}  {'ResName':<8}{'FF_Name':<8}"
            f"{'Atom':<6}{'Type':<6}{'Charge':>8}"
            f"{'X':>9}{'Y':>9}{'Z':>9}"
        )
        sep = "-" * len(header)

        try:
            with open(out_path, "w") as f:
                f.write(f"# Atom Typing Report: {Path(pdb_file).name}\n")
                f.write(f"# Force fields: "
                        f"{', '.join(p.name for p in self.ff_data.loaded_leaprcs)}\n")

                typed = [a for a in typed_atoms if a.atom_type is not None]
                untyped = [a for a in typed_atoms if a.atom_type is None]
                f.write(f"# Total: {len(typed_atoms)}  "
                        f"Typed: {len(typed)}  Untyped: {len(untyped)}\n")
                if typed:
                    f.write(f"# Total charge (typed): "
                            f"{sum(a.charge for a in typed):+.4f}\n")
                f.write(f"#\n")
                f.write(f"{header}\n")
                f.write(f"{sep}\n")

                for a in typed_atoms:
                    atype = a.atom_type if a.atom_type else "???"
                    charge = f"{a.charge:+.4f}" if a.charge is not None else "???"
                    f.write(
                        f"{a.chain_id:<6}{a.resid:>6}  {a.resname:<8}"
                        f"{a.resolved_resname:<8}{a.atom_name:<6}{atype:<6}"
                        f"{charge:>8}{a.x:>9.3f}{a.y:>9.3f}{a.z:>9.3f}\n"
                    )

            self.console.print(f"[green]Report exported to: {out_path}[/green]")
        except OSError as e:
            self.console.print(f"[red]Export failed: {e}[/red]")

    # =========================================================================
    # Helpers
    # =========================================================================

    def _parse_atom_types(self, text: str, expected: int) -> Optional[List[str]]:
        """Parse space- or dash-separated atom types."""
        if not text.strip():
            return None

        # Try dash-separated first, then space-separated
        if "-" in text:
            parts = [t.strip() for t in text.split("-") if t.strip()]
        else:
            parts = text.split()

        if len(parts) != expected:
            self.console.print(
                f"[red]Expected {expected} atom types, got {len(parts)}[/red]"
            )
            return None

        return parts

    def _display_paginated(
        self,
        items: list,
        columns: List[Tuple[str, str]],
        build_row: Callable,
        title: str,
        selectable: bool = False,
    ) -> Optional[int]:
        """
        Display items in a paginated Rich Table.

        Args:
            items: List of items to display
            columns: List of (name, style) tuples for table columns
            build_row: Callable that takes an item and returns a list of cell strings
            title: Table title
            selectable: If True, user can select an item by number

        Returns:
            Selected item index if selectable, else None
        """
        page = 0
        total_pages = max(1, (len(items) + PAGE_SIZE - 1) // PAGE_SIZE)

        while True:
            start = page * PAGE_SIZE
            end = min(start + PAGE_SIZE, len(items))
            page_items = items[start:end]

            table = Table(
                title=f"{title}" if total_pages == 1 else f"{title} — Page {page + 1}/{total_pages}",
                show_header=True,
                header_style="bold cyan",
            )

            if selectable:
                table.add_column("#", style="yellow", justify="right", width=4)

            for col_name, col_style in columns:
                table.add_column(col_name, style=col_style)

            for i, item in enumerate(page_items):
                row = build_row(item)
                if selectable:
                    table.add_row(str(start + i + 1), *row)
                else:
                    table.add_row(*row)

            self.console.print(table)

            # Single page with no selection: just show the table and return
            if total_pages == 1 and not selectable:
                return None

            # Navigation prompt — escape brackets so Rich doesn't treat them as markup
            nav_parts = []
            if page > 0:
                nav_parts.append("\\[b]ack")
            if page < total_pages - 1:
                nav_parts.append("\\[n]ext")
            if selectable:
                nav_parts.append("\\[#] select")
            if total_pages > 1 or selectable:
                nav_parts.append("\\[q]uit")

            nav_hint = "  ".join(nav_parts)
            self.console.print(f"[grey50]{nav_hint}[/grey50]")

            cmd = prompt_with_context(
                processor=self.processor,
                prompt="Navigate",
                default="n" if page < total_pages - 1 else "q",
                module="Force Field Explorer",
                description="Page navigation",
            ).strip().lower()

            if cmd == "n" and page < total_pages - 1:
                page += 1
            elif cmd == "b" and page > 0:
                page -= 1
            elif cmd == "q":
                return None
            elif selectable and cmd.isdigit():
                idx = int(cmd) - 1
                if 0 <= idx < len(items):
                    return idx
                else:
                    self.console.print(f"[red]Invalid selection: {cmd}[/red]")
            else:
                return None
