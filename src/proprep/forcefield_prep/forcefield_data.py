"""
Unified Force Field Data Manager

This module provides the SINGLE source of truth for all AMBER force field data.
Replaces the three overlapping systems:
- ForcefieldDataCollector (leaprc, lib/off files)
- ForceFieldLibrary (mol2 files)
- FFParameterReader (parmXX.dat, frcmod files)

All code should use ForceFieldData instead of the old classes.

Author: ProPrep Development Team
"""

import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Set
from dataclasses import dataclass, field

from rich.console import Console
from rich.table import Table
from rich.panel import Panel

from proprep.utils.prompts import prompt_with_context, confirm_with_context

from .ff_parsers import (
    parse_leaprc,
    parse_lib_file,
    parse_mol2_file,
    parse_parm_dat,
    parse_frcmod,
    find_amber_file,
    resolve_leaprc_path,
    AtomDefinition,
    LeaprcContents,
    ParmContents,
    FrcmodContents,
)
from .ff_types import (
    MassParameter,
    BondParameter,
    AngleParameter,
    DihedralParameter,
    ImproperParameter,
    NonbondedParameter,
)


@dataclass
class ForceFieldData:
    """
    Single source of truth for all AMBER force field data.

    Follows AMBER's leaprc loading convention:
    - Parses leaprc to find what files to load
    - Loads lib/off files for atom definitions + charges
    - Loads parmXX.dat for base bond/angle/dihedral parameters
    - Loads frcmod files for parameter modifications/additions

    Usage:
        ff_data = ForceFieldData()
        ff_data.select_and_load(console)  # Interactive selection

        # or non-interactive:
        ff_data.load_leaprc('leaprc.protein.ff14SB')
        ff_data.load_leaprc('leaprc.water.tip3p')

        # Lookups
        atom_def = ff_data.get_atom_definition('ALA', 'CA')
        bond = ff_data.get_bond_parameter('CT', 'N')
    """

    # Atom definitions from lib/off/mol2 files
    # Key: (residue_name, atom_name)
    atom_definitions: Dict[Tuple[str, str], AtomDefinition] = field(default_factory=dict)

    # Residue aliases for mapping (e.g., HEM → HEH)
    residue_aliases: Dict[str, str] = field(default_factory=dict)

    # Parameters from parmXX.dat + frcmod files
    mass_parameters: Dict[str, MassParameter] = field(default_factory=dict)
    bond_parameters: Dict[str, BondParameter] = field(default_factory=dict)
    angle_parameters: Dict[str, AngleParameter] = field(default_factory=dict)
    dihedral_parameters: Dict[str, List[DihedralParameter]] = field(default_factory=dict)
    improper_parameters: List[ImproperParameter] = field(default_factory=list)
    nonbonded_parameters: Dict[str, NonbondedParameter] = field(default_factory=dict)

    # Tracking what was loaded
    loaded_leaprcs: List[Path] = field(default_factory=list)
    loaded_libs: List[Path] = field(default_factory=list)
    loaded_frcmods: List[Path] = field(default_factory=list)
    loaded_parms: List[Path] = field(default_factory=list)
    loaded_mol2s: List[Path] = field(default_factory=list)

    # Environment
    amberhome: Optional[Path] = None
    console: Optional[Console] = None
    processor: Optional[Any] = None  # ProPrep processor for session recording context

    def __post_init__(self):
        """Initialize after dataclass creation."""
        if self.amberhome is None:
            amberhome_env = os.environ.get('AMBERHOME')
            if amberhome_env:
                self.amberhome = Path(amberhome_env)

        if self.console is None:
            self.console = Console()

    # =========================================================================
    # Loading Methods
    # =========================================================================

    def load_leaprc(self, leaprc_name: str, follow_sources: bool = True) -> int:
        """
        Load force field data from a leaprc file.

        Parses the leaprc file to find lib, frcmod, and parm.dat files to load.
        Optionally follows 'source' directives to load chained leaprc files.

        Args:
            leaprc_name: Name of leaprc file (e.g., 'leaprc.protein.ff14SB')
            follow_sources: Whether to follow 'source' directives

        Returns:
            Number of atom definitions loaded
        """
        if self.amberhome is None:
            raise RuntimeError("AMBERHOME not set")

        leaprc_path = resolve_leaprc_path(leaprc_name, self.amberhome)
        if leaprc_path is None:
            raise FileNotFoundError(f"Leaprc file not found: {leaprc_name}")

        return self._load_leaprc_recursive(leaprc_path, follow_sources, set())

    def _load_leaprc_recursive(
        self,
        leaprc_path: Path,
        follow_sources: bool,
        processed: Set[Path]
    ) -> int:
        """Recursively load leaprc file and its sources."""
        if leaprc_path in processed:
            return 0

        processed.add(leaprc_path)
        self.loaded_leaprcs.append(leaprc_path)

        contents = parse_leaprc(leaprc_path)
        atoms_loaded = 0

        # Follow source directives first (they define base force fields)
        if follow_sources:
            for source_file in contents.source_files:
                source_path = self._find_sourced_leaprc(source_file, leaprc_path)
                if source_path:
                    atoms_loaded += self._load_leaprc_recursive(
                        source_path, follow_sources, processed
                    )

        # Load lib files
        for lib_file in contents.lib_files:
            lib_path = self._find_lib_file(lib_file, leaprc_path)
            if lib_path:
                atoms_loaded += self.load_lib_file(lib_path)

        # Load parm.dat files
        for parm_file in contents.parm_files:
            parm_path = find_amber_file(parm_file, self.amberhome, 'parm')
            if parm_path:
                self.load_parm_file(parm_path)

        # Load frcmod files (these override parm.dat)
        for frcmod_file in contents.frcmod_files:
            frcmod_path = self._find_frcmod_file(frcmod_file, leaprc_path)
            if frcmod_path:
                self.load_frcmod_file(frcmod_path)

        # Load mol2 files
        for mol2_file in contents.mol2_files:
            mol2_path = find_amber_file(mol2_file, self.amberhome, 'mol2')
            if mol2_path:
                atoms_loaded += self.load_mol2_file(mol2_path)

        # Store PDB residue mappings
        for pdb_name, ff_name in contents.pdb_res_map.items():
            self.residue_aliases[pdb_name] = ff_name

        # Register bare LEaP unit aliases (e.g. `HOH = TP3`). Libs are already
        # loaded above, so the target unit's atoms are present now. Gate on the
        # target actually being a loaded unit so a stray assignment in a custom
        # leaprc can't create a dangling alias, and route through
        # add_residue_alias so the alias name also copies the target's atom
        # definitions -- this makes e.g. HOH both type AND appear in the
        # residue browser, matching how tleap treats the aliased name.
        for alias_name, unit_name in contents.unit_aliases.items():
            if alias_name in self.residue_aliases:
                continue  # an explicit addPdbResMap mapping takes precedence
            if any(key[0] == unit_name for key in self.atom_definitions):
                self.add_residue_alias(alias_name, unit_name)

        return atoms_loaded

    def load_lib_file(self, path: Path) -> int:
        """
        Load atom definitions from a lib/off file.

        Args:
            path: Path to lib file

        Returns:
            Number of atom definitions loaded
        """
        if path in self.loaded_libs:
            return 0

        self.loaded_libs.append(path)

        atoms = parse_lib_file(path, path.name)
        for atom_def in atoms:
            self.atom_definitions[atom_def.key] = atom_def

        return len(atoms)

    def load_mol2_file(self, path: Path, residue_name: str = None) -> int:
        """
        Load atom definitions from a mol2 file.

        Args:
            path: Path to mol2 file
            residue_name: Optional override for residue name

        Returns:
            Number of atom definitions loaded
        """
        if path in self.loaded_mol2s:
            return 0

        self.loaded_mol2s.append(path)

        atoms = parse_mol2_file(path, path.name)

        for atom_def in atoms:
            if residue_name:
                # Override residue name if specified
                atom_def = AtomDefinition(
                    residue_name=residue_name,
                    atom_name=atom_def.atom_name,
                    atom_type=atom_def.atom_type,
                    charge=atom_def.charge,
                    source_file=atom_def.source_file
                )
            self.atom_definitions[atom_def.key] = atom_def

        return len(atoms)

    def add_mol2_file(self, mol2_path: Path, residue_name: str = None) -> int:
        """
        Alias for load_mol2_file for backwards compatibility.

        Args:
            mol2_path: Path to mol2 file
            residue_name: Optional residue name override

        Returns:
            Number of atom definitions loaded
        """
        return self.load_mol2_file(Path(mol2_path), residue_name)

    def add_manual_atom(
        self,
        residue_name: str,
        atom_name: str,
        atom_type: str,
        atom_charge: float
    ) -> None:
        """
        Manually add an atom type definition.

        Args:
            residue_name: Residue name
            atom_name: Atom name
            atom_type: Atom type
            atom_charge: Partial charge
        """
        atom_def = AtomDefinition(
            residue_name=residue_name,
            atom_name=atom_name,
            atom_type=atom_type,
            charge=atom_charge,
            source_file="manual"
        )
        self.atom_definitions[atom_def.key] = atom_def

    def load_parm_file(self, path: Path) -> None:
        """
        Load parameters from a parmXX.dat file.

        Args:
            path: Path to parm.dat file
        """
        if path in self.loaded_parms:
            return

        self.loaded_parms.append(path)

        contents = parse_parm_dat(path, path.name)

        # Merge parameters (don't overwrite existing)
        for key, param in contents.mass_parameters.items():
            if key not in self.mass_parameters:
                self.mass_parameters[key] = param

        for key, param in contents.bond_parameters.items():
            if key not in self.bond_parameters:
                self.bond_parameters[key] = param

        for key, param in contents.angle_parameters.items():
            if key not in self.angle_parameters:
                self.angle_parameters[key] = param

        for key, params in contents.dihedral_parameters.items():
            if key not in self.dihedral_parameters:
                self.dihedral_parameters[key] = params

        for param in contents.improper_parameters:
            self.improper_parameters.append(param)

        for key, param in contents.nonbonded_parameters.items():
            if key not in self.nonbonded_parameters:
                self.nonbonded_parameters[key] = param

    def load_frcmod_file(self, path: Path) -> None:
        """
        Load parameters from a frcmod file.

        Frcmod parameters OVERRIDE existing parameters.

        Args:
            path: Path to frcmod file
        """
        if path in self.loaded_frcmods:
            return

        self.loaded_frcmods.append(path)

        contents = parse_frcmod(path, path.name)

        # Merge parameters (frcmod overrides existing)
        self.mass_parameters.update(contents.mass_parameters)
        self.bond_parameters.update(contents.bond_parameters)
        self.angle_parameters.update(contents.angle_parameters)

        for key, params in contents.dihedral_parameters.items():
            self.dihedral_parameters[key] = params  # Override

        for param in contents.improper_parameters:
            self.improper_parameters.append(param)

        self.nonbonded_parameters.update(contents.nonbonded_parameters)

    # =========================================================================
    # File Finding Helpers
    # =========================================================================

    def _find_sourced_leaprc(self, filename: str, current_leaprc: Path) -> Optional[Path]:
        """Find a sourced leaprc file."""
        leap_cmd_dir = current_leaprc.parent

        # Try relative to current directory
        path = leap_cmd_dir / filename
        if path.exists():
            return path

        # Try just the basename if it has a directory prefix
        if '/' in filename:
            path = leap_cmd_dir / Path(filename).name
            if path.exists():
                return path

        return None

    def _find_lib_file(self, filename: str, current_leaprc: Path) -> Optional[Path]:
        """Find a lib file."""
        # Common lib directories
        lib_dirs = [
            self.amberhome / 'dat' / 'leap' / 'lib',
            self.amberhome / 'dat' / 'leap' / 'parm',
            current_leaprc.parent.parent / 'lib',
        ]

        for lib_dir in lib_dirs:
            path = lib_dir / filename
            if path.exists():
                return path

        return None

    def _find_frcmod_file(self, filename: str, current_leaprc: Path) -> Optional[Path]:
        """Find a frcmod file."""
        parm_dir = self.amberhome / 'dat' / 'leap' / 'parm'

        # Try exact filename
        path = parm_dir / filename
        if path.exists():
            return path

        # Try with frcmod prefix if not present
        if not filename.startswith('frcmod'):
            path = parm_dir / f'frcmod.{filename}'
            if path.exists():
                return path

        return None

    # =========================================================================
    # Lookup Methods
    # =========================================================================

    def get_atom_definition(
        self,
        residue_name: str,
        atom_name: str
    ) -> Optional[AtomDefinition]:
        """
        Get atom definition for a residue/atom pair.

        Checks residue aliases if direct lookup fails.

        Args:
            residue_name: Residue name (e.g., 'ALA', 'HEM')
            atom_name: Atom name (e.g., 'CA', 'FE')

        Returns:
            AtomDefinition or None
        """
        # Direct lookup
        key = (residue_name, atom_name)
        if key in self.atom_definitions:
            return self.atom_definitions[key]

        # Try alias
        aliased_res = self.residue_aliases.get(residue_name)
        if aliased_res:
            key = (aliased_res, atom_name)
            if key in self.atom_definitions:
                return self.atom_definitions[key]

        return None

    def get_atom_type(self, residue_name: str, atom_name: str) -> Optional[str]:
        """Get atom type for a residue/atom pair."""
        atom_def = self.get_atom_definition(residue_name, atom_name)
        return atom_def.atom_type if atom_def else None

    def get_atom_charge(self, residue_name: str, atom_name: str) -> Optional[float]:
        """Get atom charge for a residue/atom pair."""
        atom_def = self.get_atom_definition(residue_name, atom_name)
        return atom_def.charge if atom_def else None

    def has_residue(self, residue_name: str) -> bool:
        """Check if residue is defined in force field."""
        for key in self.atom_definitions.keys():
            if key[0] == residue_name:
                return True

        # Check aliases
        if residue_name in self.residue_aliases:
            aliased = self.residue_aliases[residue_name]
            for key in self.atom_definitions.keys():
                if key[0] == aliased:
                    return True

        return False

    def get_residue_atoms(self, residue_name: str) -> List[AtomDefinition]:
        """Get all atom definitions for a residue."""
        atoms = []

        # Direct lookup
        for key, atom_def in self.atom_definitions.items():
            if key[0] == residue_name:
                atoms.append(atom_def)

        # If empty, try alias
        if not atoms and residue_name in self.residue_aliases:
            aliased = self.residue_aliases[residue_name]
            for key, atom_def in self.atom_definitions.items():
                if key[0] == aliased:
                    atoms.append(atom_def)

        return atoms

    def get_available_residues(self) -> List[str]:
        """Get list of all available residue names."""
        residues = set()
        for key in self.atom_definitions.keys():
            residues.add(key[0])
        return sorted(residues)

    def get_source_file(self, residue_name: str) -> str:
        """
        Get the source file for a residue's atom definitions.

        Args:
            residue_name: Residue name

        Returns:
            Source file name or 'unknown' if not found
        """
        # Check direct residue first
        for key, atom_def in self.atom_definitions.items():
            if key[0] == residue_name:
                return atom_def.source_file

        # Check alias
        aliased = self.residue_aliases.get(residue_name)
        if aliased:
            for key, atom_def in self.atom_definitions.items():
                if key[0] == aliased:
                    return atom_def.source_file

        return 'unknown'

    def get_mass_parameter(self, atom_type: str) -> Optional[MassParameter]:
        """Get mass parameter for an atom type."""
        return self.mass_parameters.get(atom_type.strip())

    def get_bond_parameter(self, type1: str, type2: str) -> Optional[BondParameter]:
        """
        Get bond parameter for an atom type pair.

        Handles both orderings (A-B and B-A).
        """
        type1 = type1.strip()
        type2 = type2.strip()

        # Try both orderings
        for t1, t2 in [(type1, type2), (type2, type1)]:
            key = f"{t1}-{t2}"
            if key in self.bond_parameters:
                return self.bond_parameters[key]

        # Try canonical key
        types = sorted([type1, type2])
        key = f"{types[0]}-{types[1]}"
        return self.bond_parameters.get(key)

    def get_angle_parameter(
        self,
        type1: str,
        type2: str,
        type3: str
    ) -> Optional[AngleParameter]:
        """
        Get angle parameter for an atom type triple.

        Handles both orderings (A-B-C and C-B-A).
        """
        type1 = type1.strip()
        type2 = type2.strip()
        type3 = type3.strip()

        # Try both orderings
        for t1, t3 in [(type1, type3), (type3, type1)]:
            key = f"{t1}-{type2}-{t3}"
            if key in self.angle_parameters:
                return self.angle_parameters[key]

        return None

    def get_dihedral_parameters(
        self,
        type1: str,
        type2: str,
        type3: str,
        type4: str
    ) -> List[DihedralParameter]:
        """
        Get dihedral parameters for an atom type quadruple.

        Dihedrals can have multiple terms (Fourier series).
        Handles wildcards (X) and both orderings.
        """
        type1 = type1.strip()
        type2 = type2.strip()
        type3 = type3.strip()
        type4 = type4.strip()

        # Try exact match first
        for t1, t2, t3, t4 in [
            (type1, type2, type3, type4),
            (type4, type3, type2, type1)
        ]:
            key = f"{t1}-{t2}-{t3}-{t4}"
            if key in self.dihedral_parameters:
                return self.dihedral_parameters[key]

        # Try wildcard matches (X-A-B-X pattern)
        for t1, t2, t3, t4 in [
            ('X', type2, type3, 'X'),
            ('X', type3, type2, 'X')
        ]:
            key = f"{t1}-{t2}-{t3}-{t4}"
            if key in self.dihedral_parameters:
                return self.dihedral_parameters[key]

        return []

    def get_nonbonded_parameter(self, atom_type: str) -> Optional[NonbondedParameter]:
        """Get nonbonded (VDW) parameter for an atom type."""
        return self.nonbonded_parameters.get(atom_type.strip())

    def get_atom_types_for_element(
        self,
        element: str,
        include_gaff2: bool = True
    ) -> Dict[str, List[Tuple[str, str, str, float, str]]]:
        """
        Get all atom types for a given element.

        Args:
            element: Element symbol (e.g., 'N', 'Fe', 'Zn')
            include_gaff2: Whether to include GAFF2 types (not yet implemented)

        Returns:
            Dict with 'force_field' and optionally 'gaff2' keys, each containing
            list of (residue, atom_name, atom_type, charge, description) tuples
        """
        result = {'force_field': [], 'gaff2': []}
        element_upper = element.upper()

        # Search through atom definitions for matching element
        for (resname, atom_name), atom_def in self.atom_definitions.items():
            # Infer element from atom name (first 1-2 characters)
            inferred_element = atom_name[0].upper()
            if len(atom_name) > 1 and atom_name[1].islower():
                inferred_element = atom_name[:2].upper()

            if inferred_element == element_upper:
                description = f"From {atom_def.source_file}" if atom_def.source_file else ""
                result['force_field'].append((
                    resname,
                    atom_name,
                    atom_def.atom_type,
                    atom_def.charge or 0.0,
                    description
                ))

        # Sort by residue name, then atom name
        result['force_field'].sort(key=lambda x: (x[0], x[1]))

        return result

    # =========================================================================
    # Residue Mapping
    # =========================================================================

    def add_residue_alias(self, from_name: str, to_name: str) -> int:
        """
        Add a residue alias and copy atom definitions.

        Args:
            from_name: Source residue name in structure (e.g., 'HEM')
            to_name: Target residue name in force field (e.g., 'HEH')

        Returns:
            Number of atom definitions copied
        """
        self.residue_aliases[from_name] = to_name

        # Copy atom definitions from target to source name
        copied = 0
        for key, atom_def in list(self.atom_definitions.items()):
            if key[0] == to_name:
                new_key = (from_name, key[1])
                if new_key not in self.atom_definitions:
                    self.atom_definitions[new_key] = AtomDefinition(
                        residue_name=from_name,
                        atom_name=atom_def.atom_name,
                        atom_type=atom_def.atom_type,
                        charge=atom_def.charge,
                        source_file=f"alias:{to_name}"
                    )
                    copied += 1

        return copied

    def prompt_residue_mapping(
        self,
        unknown_residues: List[str],
        structure_atoms: Dict[str, List[str]] = None
    ) -> Dict[str, str]:
        """
        Interactive interface for mapping unknown residues.

        Args:
            unknown_residues: List of residue names not in force field
            structure_atoms: Optional dict mapping residue names to atom names

        Returns:
            Dict mapping unknown residue names to force field residues
        """
        if not unknown_residues:
            return {}

        available = self.get_available_residues()
        if not available:
            self.console.print("[yellow]No residues loaded from force field[/yellow]")
            return {}

        # Show panel explaining what's happening
        self.console.print(Panel(
            "[bold]Residue Mapping Required[/bold]\n\n"
            "Some residues in your structure are not in the loaded force field.\n"
            "You can map them to available residues (e.g., HEM → HEH for heme).",
            title="Residue Mapping",
            expand=False,
        ))

        # Display available residues in columns
        table = Table(title="Available Residues in Force Field", show_header=True)
        for _ in range(4):
            table.add_column("#", style="yellow", width=4)
            table.add_column("Residue", style="cyan", width=8)

        # Fill table row by row
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

            if structure_atoms and res_name in structure_atoms:
                atom_list = structure_atoms[res_name][:10]
                atoms_str = ", ".join(atom_list)
                if len(structure_atoms[res_name]) > 10:
                    atoms_str += f" ... (+{len(structure_atoms[res_name]) - 10} more)"
                self.console.print(f"[grey50]Atoms: {atoms_str}[/grey50]")

            self.console.print(
                "[grey50]Enter number(s) to map to, comma-separated for multiple "
                "(e.g., '101,169' for HEH+PRN)[/grey50]"
            )

            choice = prompt_with_context(
                self.processor,
                f"Map '{res_name}' to",
                default="skip",
                module="Force Field Data",
                description=f"Map residue '{res_name}' to force-field residue",
            )

            if choice.lower() == 'skip':
                continue

            try:
                indices = [int(x.strip()) - 1 for x in choice.split(',')]
                target_residues = []

                for idx in indices:
                    if 0 <= idx < len(available):
                        target_residues.append(available[idx])
                    else:
                        self.console.print(f"[red]Invalid index: {idx + 1}[/red]")

                if target_residues:
                    # Add aliases for all target residues
                    total_copied = 0
                    for target_res in target_residues:
                        copied = self.add_residue_alias(res_name, target_res)
                        total_copied += copied
                        self.console.print(
                            f"  [green]✓ Copied {copied} atoms from {target_res}[/green]"
                        )

                    mappings[res_name] = target_residues[0]  # Primary mapping
                    self.console.print(
                        f"  [green]✓ Mapped {res_name} → {', '.join(target_residues)} "
                        f"({total_copied} total atoms)[/green]"
                    )

                    # Check for missing atoms
                    if structure_atoms and res_name in structure_atoms:
                        self._prompt_atom_level_mapping(
                            res_name, target_residues, structure_atoms[res_name]
                        )

            except ValueError:
                self.console.print("[yellow]Invalid input, skipping[/yellow]")

        return mappings

    def _prompt_atom_level_mapping(
        self,
        res_name: str,
        source_residues: List[str],
        structure_atom_names: List[str]
    ) -> None:
        """Prompt for atom-level mapping when atom names differ."""
        # Find which structure atoms still have no FF definition
        missing_atoms = []
        for atom_name in structure_atom_names:
            if (res_name, atom_name) not in self.atom_definitions:
                missing_atoms.append(atom_name)

        if not missing_atoms:
            self.console.print(
                f"  [green]✓ All {len(structure_atom_names)} structure atoms "
                "have FF definitions[/green]"
            )
            return

        self.console.print(
            f"\n  [yellow]⚠ {len(missing_atoms)} atoms need mapping:[/yellow] "
            f"{', '.join(missing_atoms[:10])}"
            f"{'...' if len(missing_atoms) > 10 else ''}"
        )

        # Collect available FF atoms from source residues
        ff_atoms = {}
        for source_res in source_residues:
            for key, atom_def in self.atom_definitions.items():
                if key[0] == source_res:
                    ff_atoms[atom_def.atom_name] = atom_def

        if not ff_atoms:
            return

        # Ask if user wants to do atom-level mapping
        if not confirm_with_context(
            self.processor,
            "  Do you want to map individual atoms?",
            default=False,
            module="Force Field Data",
            description="Map individual atoms for residue",
        ):
            return

        # Show available FF atoms
        ff_atom_list = sorted(ff_atoms.keys())
        self.console.print(f"  [cyan]Available FF atoms: {', '.join(ff_atom_list)}[/cyan]")

        # Map each missing atom
        for missing in missing_atoms[:20]:  # Limit to first 20
            choice = prompt_with_context(
                self.processor,
                f"  Map '{missing}' to (or 'skip')",
                default="skip",
                module="Force Field Data",
                description=f"Map atom '{missing}' to force-field atom",
            )

            if choice.lower() == 'skip':
                continue

            if choice in ff_atoms:
                ff_atom = ff_atoms[choice]
                self.atom_definitions[(res_name, missing)] = AtomDefinition(
                    residue_name=res_name,
                    atom_name=missing,
                    atom_type=ff_atom.atom_type,
                    charge=ff_atom.charge,
                    source_file=f"mapped:{choice}"
                )
                self.console.print(
                    f"    [green]✓ {missing} → {choice} "
                    f"(type={ff_atom.atom_type}, charge={ff_atom.charge:.4f})[/green]"
                )
            else:
                self.console.print(f"    [yellow]'{choice}' not found[/yellow]")

    # =========================================================================
    # Interactive Selection
    # =========================================================================

    def select_and_load(self, console: Console = None, processor=None) -> int:
        """
        Interactive force field selection.

        Two-step selection:
        1. Select categories (Protein, Water, GAFF, etc.)
        2. Select specific force fields from those categories

        Args:
            console: Rich console for output

        Returns:
            Number of atom definitions loaded
        """
        if console:
            self.console = console

        if self.amberhome is None:
            raise RuntimeError("AMBERHOME not set")

        leap_cmd_dir = self.amberhome / 'dat' / 'leap' / 'cmd'

        # Get all leaprc files
        leaprc_files = sorted([
            f for f in leap_cmd_dir.iterdir()
            if f.name.startswith('leaprc')
        ])

        if not leaprc_files:
            raise FileNotFoundError(f"No leaprc files found in {leap_cmd_dir}")

        # Categorize files
        categorized = {}
        for leaprc_file in leaprc_files:
            category = self._get_leaprc_category(leaprc_file.name)
            if category not in categorized:
                categorized[category] = []
            categorized[category].append(leaprc_file)

        # =====================================================================
        # Step 1: Select Categories
        # =====================================================================
        category_order = [
            "Protein", "Modified Amino Acids", "DNA", "RNA", "Water", "GAFF (General)",
            "Lipid", "Carbohydrate", "Dye", "Constant pH/Redox", "Other"
        ]

        available_categories = [cat for cat in category_order if cat in categorized]

        bundled_count = self._count_bundled_redox_sets()
        if bundled_count > 0:
            available_categories.append("Redox Cofactors (built-in)")
        available_categories.append("Custom/User-Generated")

        cat_table = Table(
            title="Step 1: Select Force Field Categories",
            show_header=True,
            header_style="bold cyan"
        )
        cat_table.add_column("Index", style="yellow", no_wrap=True)
        cat_table.add_column("Category", style="cyan")
        cat_table.add_column("Available Files", style="green", justify="right")

        for idx, category in enumerate(available_categories, 1):
            if category == "Custom/User-Generated":
                cat_table.add_row(str(idx), category, "[grey50]User files[/grey50]")
            elif category == "Redox Cofactors (built-in)":
                cat_table.add_row(str(idx), category, str(bundled_count))
            else:
                file_count = len(categorized[category])
                cat_table.add_row(str(idx), category, str(file_count))

        self.console.print(cat_table)

        # Build the hint/default indices dynamically from the category list so
        # they stay correct regardless of which categories are present or where
        # a new category (e.g. "Modified Amino Acids") sits in the order.
        def _ci(name):
            return available_categories.index(name) + 1 if name in available_categories else None
        prot_i, wat_i = _ci("Protein"), _ci("Water")
        gaff_i, dna_i = _ci("GAFF (General)"), _ci("DNA")

        self.console.print("\n[grey50]Common combinations:[/grey50]")
        if prot_i and wat_i:
            self.console.print(f"[grey50]  • Protein simulations: {prot_i},{wat_i} (Protein + Water)[/grey50]")
        if prot_i and wat_i and gaff_i:
            self.console.print(f"[grey50]  • Protein with ligands: {prot_i},{wat_i},{gaff_i} (Protein + Water + GAFF)[/grey50]")
        if dna_i and wat_i:
            self.console.print(f"[grey50]  • DNA simulations: {dna_i},{wat_i} (DNA + Water)[/grey50]")
        self.console.print("")

        default_cat = ",".join(str(i) for i in (prot_i, wat_i) if i) or "1"
        cat_selection = prompt_with_context(processor,
            "Select categories (comma-separated)",
            default=default_cat,
            module="Force Field Explorer",
            description="Select force field categories",
        )

        # Parse category selections
        selected_category_indices = []
        for s in cat_selection.split(','):
            try:
                idx = int(s.strip())
                if 1 <= idx <= len(available_categories):
                    selected_category_indices.append(idx)
            except ValueError:
                pass

        selected_categories = [
            available_categories[i - 1] for i in selected_category_indices
        ]

        # Handle special categories (not driven by `categorized`)
        has_custom = "Custom/User-Generated" in selected_categories
        if has_custom:
            selected_categories.remove("Custom/User-Generated")

        has_bundled = "Redox Cofactors (built-in)" in selected_categories
        if has_bundled:
            selected_categories.remove("Redox Cofactors (built-in)")

        if not selected_categories and not has_custom and not has_bundled:
            self.console.print("[yellow]No categories selected[/yellow]")
            return 0

        # =====================================================================
        # Step 2: Select Specific Force Fields
        # =====================================================================
        total_atoms = 0

        if selected_categories:
            filtered_files = []
            for category in selected_categories:
                filtered_files.extend(categorized[category])

            # Enrich the scanned files with the curated catalog shared with the
            # Topology Generator: friendly names, descriptions, recommended (★)
            # marks, and the add-on (⚠) flag for modified-AA FFs that need a
            # base protein FF. Files not in the catalog (the long tail the
            # Explorer still exposes) fall back to a derived label.
            from proprep.forcefield_params.forcefield_catalog import build_leaprc_index
            catalog = build_leaprc_index()

            self.console.print(
                f"\n[bold magenta]Step 2: Select Force Fields[/bold magenta] "
                f"[grey50]({', '.join(selected_categories)})[/grey50]\n"
            )

            # Walk categories in the same order filtered_files was built so the
            # printed index matches filtered_files[idx-1] used for selection.
            idx = 0
            for category in selected_categories:
                files = categorized.get(category, [])
                if not files:
                    continue
                header = f"[bold cyan]{category}[/bold cyan]"
                if category == "Modified Amino Acids":
                    header += "  [yellow](add-ons — combine with a base protein FF)[/yellow]"
                self.console.print(header)
                for leaprc_file in files:
                    idx += 1
                    meta = catalog.get(leaprc_file.name)
                    if meta:
                        disp, desc = meta['name'], meta['description']
                        recommended, addon = meta['recommended'], meta['is_addon']
                    else:
                        disp = leaprc_file.name.replace('leaprc.', '')
                        desc = ""
                        recommended = False
                        addon = (category == "Modified Amino Acids")
                    if len(desc) > 50:
                        desc = desc[:47] + "..."
                    if recommended:
                        marker = "[green]★[/green]"
                    elif addon:
                        marker = "[yellow]⚠[/yellow]"
                    else:
                        marker = " "
                    self.console.print(
                        f"  [yellow]{idx:>3}[/yellow] {marker} "
                        f"[cyan]{disp:<26}[/cyan] [grey50]{desc}[/grey50]"
                    )
                self.console.print("")

            self.console.print(
                f"[grey50]Showing {len(filtered_files)} force field(s).  "
                f"[green]★[/green] recommended · [yellow]⚠[/yellow] add-on "
                f"(needs a base protein FF)[/grey50]\n"
            )

            # Generate smart default
            default = self._get_smart_default(selected_categories, filtered_files)

            ff_selection = prompt_with_context(processor,
                "Select force field(s) (comma-separated)",
                default=default,
                module="Force Field Explorer",
                description="Select force field files",
            )

            # Parse selections
            selected_leaprcs = []
            for s in ff_selection.split(','):
                try:
                    idx = int(s.strip())
                    if 1 <= idx <= len(filtered_files):
                        selected_leaprcs.append(filtered_files[idx - 1])
                except ValueError:
                    pass

            if not selected_leaprcs and not has_custom:
                self.console.print("[yellow]No force fields selected[/yellow]")
                return 0

            if selected_leaprcs:
                # Load selected force fields
                self.console.print(f"\n[cyan]Loading: {', '.join(f.name for f in selected_leaprcs)}[/cyan]\n")

                for leaprc in selected_leaprcs:
                    atoms = self._load_leaprc_recursive(leaprc, True, set())
                    total_atoms += atoms

        # =====================================================================
        # Step 3: Load Bundled Redox Cofactor Parameters (if selected)
        # =====================================================================
        if has_bundled:
            bundled_atoms = self._load_bundled_redox_centers(self.console, processor=processor)
            total_atoms += bundled_atoms

        # =====================================================================
        # Step 4: Load Custom Force Field Files (if selected)
        # =====================================================================
        if has_custom:
            custom_atoms = self._load_custom_files(self.console, processor=processor)
            total_atoms += custom_atoms

        if total_atoms == 0:
            self.console.print("[yellow]No atom definitions loaded.[/yellow]")
            return 0

        self.console.print(
            f"\n[green]✓ Loaded {len(self.atom_definitions)} atom definitions "
            f"from {len(self.loaded_libs)} library files[/green]"
        )

        if self.bond_parameters:
            n_bonds = len(set(p.key for p in self.bond_parameters.values()))
            n_angles = len(set(p.key for p in self.angle_parameters.values()))
            self.console.print(
                f"[green]✓ Loaded {n_bonds} bond, {n_angles} angle parameters "
                f"from {len(self.loaded_frcmods)} frcmod file(s)[/green]"
            )

        return total_atoms

    # Supported custom force field file extensions
    CUSTOM_FF_EXTENSIONS = {'.lib', '.off', '.frcmod', '.mol2', '.prep'}

    def load_prep_file(self, path: Path) -> int:
        """
        Load atom definitions from an AMBER .prep file.

        Args:
            path: Path to prep file

        Returns:
            Number of atom definitions loaded
        """
        if path in self.loaded_libs:
            return 0

        from proprep.oniom_prep.amber_lib_parser import AmberPrepParser

        parser = AmberPrepParser(str(path))
        residues = parser.parse()

        count = 0
        for resname, res_data in residues.items():
            for atom in res_data.get("atoms", []):
                atom_def = AtomDefinition(
                    residue_name=resname,
                    atom_name=atom["name"],
                    atom_type=atom["type"],
                    charge=atom["charge"],
                    source_file=path.name,
                )
                self.atom_definitions[atom_def.key] = atom_def
                count += 1

        self.loaded_libs.append(path)
        return count

    def _load_custom_files(self, console: Console, processor=None) -> int:
        """
        Interactive loader for custom force field files using a file browser.

        Uses a terminal file browser (matching ProPrep's UX elsewhere) to let
        the user navigate to and select .lib, .off, .frcmod, .mol2, and .prep files.

        Returns:
            Number of atom definitions loaded from custom files
        """
        console.print(Panel(
            "[bold]Custom Force Field File Loader[/bold]\n\n"
            "Browse and select your custom force field files.\n"
            "Supported formats: [cyan].lib[/cyan], [cyan].off[/cyan], "
            "[cyan].frcmod[/cyan], [cyan].mol2[/cyan], [cyan].prep[/cyan]\n\n"
            "Select one or more files (e.g. [cyan]2,3[/cyan] or [cyan]2-4[/cyan]). "
            "The browser remembers the last directory across visits.",
            title="Custom/User-Generated",
            border_style="magenta",
            expand=False,
        ))

        total_atoms = 0
        loaded_files = []
        last_dir = "."

        while True:
            result = self._browse_ff_files(console, directory=last_dir, processor=processor)

            if not result:
                break

            picked, last_dir = result

            for raw_path in picked:
                path = Path(raw_path)
                count = self._load_single_custom_file(path, console)
                if count is not None:
                    total_atoms += count
                    loaded_files.append(path)

            if not confirm_with_context(
                processor,
                "Load another file?",
                default=True,
                module="Force Field Explorer",
                description="Load another custom force field file",
            ):
                break

        # Summary
        if loaded_files:
            console.print(f"\n[green]Loaded {len(loaded_files)} custom file(s) "
                          f"({total_atoms} atom definitions)[/green]")
            for p in loaded_files:
                console.print(f"  [grey50]{p.name}[/grey50]")
        else:
            console.print("[grey50]No custom files loaded.[/grey50]")

        return total_atoms

    def _enumerate_bundled_redox_sets(self) -> List[Dict[str, Any]]:
        """Discover all bundled redox-cofactor forcefield sets.

        Walks ``forcefield_params/specialized_residues`` (and the user override at
        ``~/.proprep/forcefield_params/specialized_residues``) via the loader API and
        returns a flat list of selectable rows. Each row carries enough metadata
        to display in a table and locate the actual .frcmod / .lib file paths.
        """
        try:
            from proprep.forcefield_params.loader import (
                discover_forcefield_files,
                get_available_cofactor_types,
            )
        except ImportError:
            return []

        rows: List[Dict[str, Any]] = []
        try:
            cofactor_types = get_available_cofactor_types()
        except Exception:
            return []

        for cofactor_path, info in cofactor_types.items():
            if not info.get("valid", False):
                continue
            for redox_state, rdata in info.get("redox_states", {}).items():
                for spin_state in rdata.get("spin_states", {}):
                    try:
                        sets = discover_forcefield_files(
                            cofactor_path, redox_state, spin_state
                        )
                    except Exception:
                        continue
                    for ff_set in sets:
                        rows.append({
                            "cofactor": cofactor_path,
                            "redox": redox_state,
                            "spin": spin_state,
                            "name": ff_set["name"],
                            "description": ff_set.get("description", ""),
                            "frcmod": ff_set["frcmod"],
                            "lib": ff_set["lib"],
                            "is_default": ff_set.get("is_default", False),
                        })
        rows.sort(key=lambda r: (r["cofactor"], r["redox"], r["spin"], r["name"]))
        return rows

    def _count_bundled_redox_sets(self) -> int:
        """Cheap count for the Step 1 category table; safe if discovery fails."""
        try:
            return len(self._enumerate_bundled_redox_sets())
        except Exception:
            return 0

    def _load_bundled_redox_centers(self, console: Console, processor=None) -> int:
        """Interactive picker for built-in redox-cofactor forcefield sets."""
        rows = self._enumerate_bundled_redox_sets()
        if not rows:
            console.print("[yellow]No bundled redox cofactor parameters found.[/yellow]")
            return 0

        console.print(Panel(
            "[bold]Built-in Redox Cofactor Parameters[/bold]\n\n"
            "Bundled .frcmod / .lib sets for hemes, iron-sulfur clusters, and other "
            "redox-active cofactors. Each row is one parameter set; selecting it "
            "loads the matching .frcmod and .lib file(s).",
            title="Redox Cofactors (built-in)",
            border_style="magenta",
            expand=False,
        ))

        table = Table(
            title="Available Redox Cofactor Parameter Sets",
            show_header=True,
            header_style="bold cyan",
            show_lines=False,
        )
        table.add_column("Index", style="yellow", no_wrap=True, width=6)
        table.add_column("Cofactor", style="cyan", no_wrap=True)
        table.add_column("Redox", style="magenta", no_wrap=True)
        table.add_column("Spin / Variant", style="magenta", no_wrap=True)
        table.add_column("Set", style="green", no_wrap=True)
        table.add_column("Notes", style="grey50")

        for idx, row in enumerate(rows, 1):
            note_bits = []
            if row["is_default"]:
                note_bits.append("default")
            if row["description"]:
                note_bits.append(row["description"])
            note = " — ".join(note_bits)
            table.add_row(
                str(idx),
                row["cofactor"],
                row["redox"],
                row["spin"],
                row["name"],
                note,
            )
        console.print(table)
        console.print(
            f"\n[bold]Showing {len(rows)} parameter set(s)[/bold]  "
            "[grey50](use commas or ranges, e.g. 1,3 or 2-4)[/grey50]\n"
        )

        selection = prompt_with_context(
            processor,
            "Select set(s) (comma/range or 'none')",
            default="none",
            module="Force Field Explorer",
            description="Select bundled redox cofactor parameter sets",
        ).strip()

        if not selection or selection.lower() in ("none", "skip", "exit", "0"):
            console.print("[grey50]No bundled sets selected.[/grey50]")
            return 0

        indices = self._parse_select_spec(selection, len(rows), console)
        if not indices:
            return 0

        total_atoms = 0
        loaded_files: List[Path] = []
        for n in indices:
            row = rows[n - 1]
            console.print(
                f"\n[cyan]Loading {row['cofactor']} / {row['redox']} / "
                f"{row['spin']} / {row['name']}[/cyan]"
            )
            for raw_path in self._resolve_set_files(row):
                p = Path(raw_path)
                count = self._load_single_custom_file(p, console)
                if count is not None:
                    total_atoms += count
                    loaded_files.append(p)

        if loaded_files:
            console.print(
                f"\n[green]Loaded {len(loaded_files)} bundled file(s) "
                f"({total_atoms} atom definitions)[/green]"
            )
        return total_atoms

    @staticmethod
    def _resolve_set_files(row: Dict[str, Any]) -> List[str]:
        """Flatten a forcefield set's frcmod + lib reference(s) into a path list.

        ``lib`` may be a single path string or a list of paths (the multi-mol2
        MCPB case). frcmod is always a single path.
        """
        files: List[str] = []
        frcmod = row.get("frcmod")
        if frcmod:
            files.append(frcmod)
        lib_ref = row.get("lib")
        if isinstance(lib_ref, list):
            files.extend(lib_ref)
        elif lib_ref:
            files.append(lib_ref)
        return files

    def _load_single_custom_file(self, path: Path, console: Console) -> Optional[int]:
        """Load a single custom FF file. Returns atom count or None on error."""
        suffix = path.suffix.lower()

        try:
            if suffix in ('.lib', '.off'):
                atoms = self.load_lib_file(path)
                console.print(f"[green]✓ Loaded {atoms} atom definitions from {path.name}[/green]")
                return atoms
            elif suffix == '.prep':
                atoms = self.load_prep_file(path)
                console.print(f"[green]✓ Loaded {atoms} atom definitions from {path.name}[/green]")
                return atoms
            elif suffix == '.frcmod':
                self.load_frcmod_file(path)
                console.print(f"[green]✓ Loaded parameters from {path.name}[/green]")
                return 0
            elif suffix == '.mol2':
                atoms = self.load_mol2_file(path)
                console.print(f"[green]✓ Loaded {atoms} atom definitions from {path.name}[/green]")
                return atoms
            else:
                console.print(
                    f"[yellow]Unsupported file type: {suffix}  "
                    f"(expected .lib, .off, .frcmod, .mol2, or .prep)[/yellow]"
                )
                return None
        except Exception as e:
            console.print(f"[red]Error loading {path.name}: {e}[/red]")
            return None

    def _browse_ff_files(
        self, console: Console, directory: str = ".", processor=None
    ) -> Optional[Tuple[List[str], str]]:
        """
        Terminal file browser for force field files (.lib, .off, .frcmod, .mol2,
        .prep). Thin wrapper over the shared file browser: unified bare-N / q UX,
        comma+range multi-select, and filename-based session replay.

        Args:
            console: Rich console for display
            directory: Starting directory
            processor: Processor for session recording

        Returns:
            (selected_paths, exit_directory) tuple, or None if canceled.
            exit_directory is the directory of the last selection so callers can
            resume the browser from where the user left off; it is reconstructed
            from the chosen files (multi-select picks share one directory).
        """
        from proprep.utils.file_browser import file_browser, default_size_detail

        selected = file_browser(
            directory=directory,
            extensions=list(self.CUSTOM_FF_EXTENSIONS),
            console=console,
            processor=processor,
            multi=True,
            label="FF file",
            entry_detail=default_size_detail,
            module="Force Field Explorer",
        )
        if not selected:
            return None
        return selected, os.path.dirname(selected[0])



    def _parse_select_spec(
        self, spec: str, max_index: int, console: Console
    ) -> Optional[List[int]]:
        """Parse a 'select' argument supporting 'N', 'N,M', and 'N-M' forms.

        Returns a deduplicated, ordered list of 1-based indices, or None if
        the spec was malformed or out of range (after printing an error).
        """
        if not spec:
            console.print("[red]Usage: select N, select N,M, or select N-M[/red]")
            return None

        seen = set()
        result: List[int] = []
        for token in spec.split(','):
            token = token.strip()
            if not token:
                continue
            if '-' in token:
                lo_str, hi_str = token.split('-', 1)
                try:
                    lo, hi = int(lo_str), int(hi_str)
                except ValueError:
                    console.print(f"[red]Invalid range: {token!r}[/red]")
                    return None
                if lo > hi:
                    lo, hi = hi, lo
                rng = range(lo, hi + 1)
            else:
                try:
                    rng = [int(token)]
                except ValueError:
                    console.print(f"[red]Invalid number: {token!r}[/red]")
                    return None
            for n in rng:
                if not (1 <= n <= max_index):
                    console.print(
                        f"[red]Out of range: {n} (choose 1-{max_index})[/red]"
                    )
                    return None
                if n not in seen:
                    seen.add(n)
                    result.append(n)
        if not result:
            console.print("[red]No items selected[/red]")
            return None
        return result

    def _get_leaprc_category(self, filename: str) -> str:
        """Categorize a leaprc file by its name."""
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
        elif ('phos' in name_lower or 'modaa' in name_lower
              or 'mimetic' in name_lower or 'fluorine' in name_lower):
            # Modified / non-standard amino-acid ADD-ONS (phosphorylated, modAA,
            # mimetic, fluorinated). These require a base protein FF and must
            # NOT be grouped with standalone protein FFs — selecting one alone
            # leaves the protein untyped. Checked before 'protein' because
            # several of these filenames also contain 'protein'.
            return "Modified Amino Acids"
        elif 'protein' in name_lower:
            return "Protein"
        elif 'water' in name_lower:
            return "Water"
        else:
            return "Other"

    def _get_smart_default(
        self,
        selected_categories: List[str],
        filtered_files: List[Path]
    ) -> str:
        """Generate smart default selections based on common patterns."""
        defaults = []
        picked_protein = picked_water = picked_gaff = False

        for idx, leaprc_file in enumerate(filtered_files, 1):
            filename = leaprc_file.name.lower()

            # Prefer ff14SB for protein. Match the exact stem so we don't also
            # grab ff14SBonlysc / ff14SB_modAA (substring "ff14sb"), which would
            # default to loading two conflicting protein FFs at once.
            if "Protein" in selected_categories and not picked_protein:
                if filename == "leaprc.protein.ff14sb":
                    defaults.append(idx)
                    picked_protein = True
                    continue

            # Prefer TIP3P for water
            if "Water" in selected_categories and not picked_water:
                if filename == "leaprc.water.tip3p":
                    defaults.append(idx)
                    picked_water = True
                    continue

            # Prefer GAFF2 for general
            if "GAFF (General)" in selected_categories and not picked_gaff:
                if "gaff2" in filename:
                    defaults.append(idx)
                    picked_gaff = True
                    continue

        # If no smart defaults, use first file from each category
        if not defaults:
            seen_categories = set()
            for idx, leaprc_file in enumerate(filtered_files, 1):
                category = self._get_leaprc_category(leaprc_file.name)
                if category not in seen_categories:
                    defaults.append(idx)
                    seen_categories.add(category)

        return ",".join(str(d) for d in defaults)

    # =========================================================================
    # Statistics
    # =========================================================================

    def get_statistics(self) -> Dict:
        """Get statistics about loaded data."""
        return {
            'n_atom_definitions': len(self.atom_definitions),
            'n_residues': len(self.get_available_residues()),
            'n_aliases': len(self.residue_aliases),
            'n_mass_params': len(self.mass_parameters),
            'n_bond_params': len(set(p.key for p in self.bond_parameters.values())),
            'n_angle_params': len(set(p.key for p in self.angle_parameters.values())),
            'n_dihedral_params': len(self.dihedral_parameters),
            'n_nonbonded_params': len(self.nonbonded_parameters),
            'n_leaprcs_loaded': len(self.loaded_leaprcs),
            'n_libs_loaded': len(self.loaded_libs),
            'n_frcmods_loaded': len(self.loaded_frcmods),
        }

    def __repr__(self) -> str:
        stats = self.get_statistics()
        return (
            f"ForceFieldData("
            f"atoms={stats['n_atom_definitions']}, "
            f"residues={stats['n_residues']}, "
            f"bonds={stats['n_bond_params']}, "
            f"angles={stats['n_angle_params']})"
        )
