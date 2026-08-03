"""
Biological Assembly Generator

Generates biological assemblies from crystal asymmetric units by applying
REMARK 350 transformation matrices or user-provided matrices.
"""

import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

import numpy as np

from proprep.utils.module_registry import ProcessingModule, register_module
from proprep.utils.prompts import prompt_with_context, confirm_with_context, int_prompt_with_context, float_prompt_with_context

logger = logging.getLogger(__name__)

# Chain ID sequence: A-Z, a-z, 0-9 (62 total)
CHAIN_IDS = (
    list('ABCDEFGHIJKLMNOPQRSTUVWXYZ') +
    list('abcdefghijklmnopqrstuvwxyz') +
    list('0123456789')
)


def parse_remark_350(remark_lines: List[str]) -> Dict[int, Dict]:
    """
    Parse REMARK 350 records into structured assembly data.

    Args:
        remark_lines: List of REMARK 350 lines from PDB file

    Returns:
        Dictionary mapping assembly number to assembly data:
        {
            1: {
                'description': 'DIMERIC',
                'author_determined': True,
                'software_determined': False,
                'chains': ['A', 'B'],
                'transforms': [
                    {'id': 1, 'matrix': np.array(3x3), 'translation': np.array(3)},
                    ...
                ]
            },
            ...
        }
    """
    assemblies = {}
    current_assembly = None
    current_transform_id = None
    current_matrix_rows = []
    current_translation = []

    for line in remark_lines:
        line = line.strip()

        # New biomolecule
        match = re.match(r'BIOMOLECULE:\s*(\d+)', line)
        if match:
            # Save previous assembly if exists
            if current_assembly is not None and current_transform_id is not None:
                _finalize_transform(assemblies, current_assembly, current_transform_id,
                                   current_matrix_rows, current_translation)

            current_assembly = int(match.group(1))
            assemblies[current_assembly] = {
                'description': '',
                'author_determined': False,
                'software_determined': False,
                'chains': [],
                'transforms': []
            }
            current_transform_id = None
            current_matrix_rows = []
            current_translation = []
            continue

        if current_assembly is None:
            continue

        # Author determined biological unit
        if 'AUTHOR DETERMINED BIOLOGICAL UNIT' in line:
            assemblies[current_assembly]['author_determined'] = True
            # Extract description (e.g., "DIMERIC")
            match = re.search(r'UNIT:\s*(\w+)', line)
            if match:
                assemblies[current_assembly]['description'] = match.group(1)
            continue

        # Software determined
        if 'SOFTWARE DETERMINED' in line:
            assemblies[current_assembly]['software_determined'] = True
            match = re.search(r'QUATERNARY STRUCTURE:\s*(\w+)', line)
            if match:
                if not assemblies[current_assembly]['description']:
                    assemblies[current_assembly]['description'] = match.group(1)
            continue

        # Chains to apply transformation to
        match = re.match(r'APPLY THE FOLLOWING TO CHAINS:\s*(.+)', line)
        if match:
            chains_str = match.group(1)
            chains = [c.strip() for c in chains_str.split(',')]
            assemblies[current_assembly]['chains'].extend(chains)
            continue

        # Additional chains (continuation)
        match = re.match(r'AND CHAINS:\s*(.+)', line)
        if match:
            chains_str = match.group(1)
            chains = [c.strip() for c in chains_str.split(',')]
            assemblies[current_assembly]['chains'].extend(chains)
            continue

        # BIOMT records (transformation matrices)
        # Format: BIOMT1   1  1.000000  0.000000  0.000000        0.00000
        match = re.match(r'BIOMT(\d)\s+(\d+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', line)
        if match:
            row_num = int(match.group(1))  # 1, 2, or 3
            transform_id = int(match.group(2))

            # If this is a new transform, save the previous one
            if transform_id != current_transform_id:
                if current_transform_id is not None:
                    _finalize_transform(assemblies, current_assembly, current_transform_id,
                                       current_matrix_rows, current_translation)
                current_transform_id = transform_id
                current_matrix_rows = []
                current_translation = []

            # Extract matrix row and translation component
            m1 = float(match.group(3))
            m2 = float(match.group(4))
            m3 = float(match.group(5))
            t = float(match.group(6))

            current_matrix_rows.append([m1, m2, m3])
            current_translation.append(t)

    # Finalize last transform
    if current_assembly is not None and current_transform_id is not None:
        _finalize_transform(assemblies, current_assembly, current_transform_id,
                           current_matrix_rows, current_translation)

    # Remove duplicates from chains lists
    for assembly in assemblies.values():
        assembly['chains'] = list(dict.fromkeys(assembly['chains']))

    return assemblies


def _finalize_transform(assemblies: Dict, assembly_num: int, transform_id: int,
                        matrix_rows: List, translation: List):
    """Helper to finalize a transformation matrix."""
    if len(matrix_rows) == 3 and len(translation) == 3:
        assemblies[assembly_num]['transforms'].append({
            'id': transform_id,
            'matrix': np.array(matrix_rows),
            'translation': np.array(translation)
        })


def extract_remark_350_from_pdb(pdb_path: str) -> List[str]:
    """
    Extract REMARK 350 lines from a PDB file.

    Args:
        pdb_path: Path to PDB file

    Returns:
        List of REMARK 350 content lines (without the "REMARK 350" prefix)
    """
    remark_lines = []

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('REMARK 350'):
                # Remove "REMARK 350" prefix and keep the rest
                content = line[10:].strip()
                if content:
                    remark_lines.append(content)

    return remark_lines


def apply_assembly_transforms(
    pdb_path: str,
    assembly_data: Dict,
    output_path: str,
    include_hetatm: bool = True,
    include_water: bool = False
) -> Dict[str, Any]:
    """
    Apply transformation matrices to generate biological assembly.

    Args:
        pdb_path: Path to input PDB file (asymmetric unit)
        assembly_data: Assembly data from parse_remark_350 for one assembly
        output_path: Path for output PDB file
        include_hetatm: Include HETATM records (ligands, cofactors)
        include_water: Include water molecules

    Returns:
        Dictionary with statistics:
        {
            'original_chains': ['A', 'B'],
            'original_atoms': 4521,
            'new_chains': ['A', 'B', 'C', 'D'],
            'new_atoms': 9042,
            'transforms_applied': 2
        }
    """
    chains_to_transform = set(assembly_data['chains'])
    transforms = assembly_data['transforms']

    # Read atoms from applicable chains
    atoms_by_chain = {}  # chain_id -> list of atom lines
    other_records = []   # HEADER, TITLE, etc.

    with open(pdb_path, 'r') as f:
        for line in f:
            record_type = line[:6].strip()

            if record_type in ('ATOM', 'HETATM'):
                chain_id = line[21] if len(line) > 21 else ' '

                # Skip waters if not requested
                if not include_water:
                    res_name = line[17:20].strip()
                    if res_name in ('HOH', 'WAT', 'TIP', 'TIP3', 'SOL'):
                        continue

                # Skip HETATM if not requested
                if record_type == 'HETATM' and not include_hetatm:
                    continue

                if chain_id in chains_to_transform:
                    if chain_id not in atoms_by_chain:
                        atoms_by_chain[chain_id] = []
                    atoms_by_chain[chain_id].append(line)

            elif record_type in ('HEADER', 'TITLE', 'COMPND', 'SOURCE', 'EXPDTA',
                                'AUTHOR', 'REVDAT', 'REMARK', 'CRYST1'):
                other_records.append(line)

    # Calculate total chains needed
    original_chains = sorted(atoms_by_chain.keys())
    num_original_chains = len(original_chains)
    num_transforms = len(transforms)
    total_chains_needed = num_original_chains * num_transforms

    if total_chains_needed > len(CHAIN_IDS):
        logger.warning(
            f"Assembly requires {total_chains_needed} chains but only {len(CHAIN_IDS)} "
            f"chain IDs available. Some chains may have duplicate IDs."
        )

    # Generate new chain ID mapping
    # For each transform, assign new chain IDs to the original chains
    chain_mapping = {}  # (original_chain, transform_id) -> new_chain
    chain_idx = 0

    for transform in transforms:
        transform_id = transform['id']
        for orig_chain in original_chains:
            if chain_idx < len(CHAIN_IDS):
                new_chain = CHAIN_IDS[chain_idx]
            else:
                # Overflow - reuse chain IDs with warning
                new_chain = CHAIN_IDS[chain_idx % len(CHAIN_IDS)]
            chain_mapping[(orig_chain, transform_id)] = new_chain
            chain_idx += 1

    # Apply transformations and write output
    output_lines = []

    # Add header records
    output_lines.extend(other_records)

    # Add remark about assembly generation
    output_lines.append(f"REMARK   1 BIOLOGICAL ASSEMBLY GENERATED BY PROPREP\n")
    output_lines.append(f"REMARK   1 TRANSFORMS APPLIED: {num_transforms}\n")
    output_lines.append(f"REMARK   1 ORIGINAL CHAINS: {', '.join(original_chains)}\n")

    atom_serial = 1
    new_chains_list = []

    for transform in transforms:
        transform_id = transform['id']
        R = transform['matrix']
        T = transform['translation']

        for orig_chain in original_chains:
            new_chain = chain_mapping[(orig_chain, transform_id)]
            if new_chain not in new_chains_list:
                new_chains_list.append(new_chain)

            for atom_line in atoms_by_chain[orig_chain]:
                # Parse coordinates
                x = float(atom_line[30:38])
                y = float(atom_line[38:46])
                z = float(atom_line[46:54])

                # Apply transformation
                coords = np.array([x, y, z])
                new_coords = R @ coords + T

                # Build new atom line
                new_line = (
                    atom_line[:6] +                           # Record type
                    f"{atom_serial:5d}" +                     # Atom serial
                    atom_line[11:21] +                        # Atom name, altLoc, resName
                    new_chain +                               # Chain ID
                    atom_line[22:30] +                        # resSeq, iCode
                    f"{new_coords[0]:8.3f}" +                 # X
                    f"{new_coords[1]:8.3f}" +                 # Y
                    f"{new_coords[2]:8.3f}" +                 # Z
                    atom_line[54:]                            # occupancy, tempFactor, element, charge
                )
                output_lines.append(new_line)
                atom_serial += 1

            # Add TER record after each chain
            output_lines.append(f"TER   {atom_serial:5d}\n")
            atom_serial += 1

    output_lines.append("END\n")

    # Write output file
    with open(output_path, 'w') as f:
        f.writelines(output_lines)

    # Calculate statistics
    original_atom_count = sum(len(atoms) for atoms in atoms_by_chain.values())

    return {
        'original_chains': original_chains,
        'original_atoms': original_atom_count,
        'new_chains': new_chains_list,
        'new_atoms': original_atom_count * num_transforms,
        'transforms_applied': num_transforms
    }


def create_manual_transform() -> Dict:
    """
    Interactively create a transformation matrix from user input.

    Returns:
        Transform dictionary with 'id', 'matrix', and 'translation' keys
    """
    print("\nEnter transformation matrix (3x3 rotation matrix):")
    matrix_rows = []
    for i in range(3):
        while True:
            row_str = prompt_with_context(
                self.processor, f"Row {i+1} (3 values separated by spaces)",
                module="Biological Assembly",
                description="Rotation matrix row",
            ).strip()
            try:
                values = [float(x) for x in row_str.split()]
                if len(values) == 3:
                    matrix_rows.append(values)
                    break
                else:
                    print("    Please enter exactly 3 values.")
            except ValueError:
                print("    Invalid input. Please enter 3 numeric values.")

    print("\nEnter translation vector (3 values):")
    while True:
        trans_str = prompt_with_context(
            self.processor, "Translation (x y z)",
            module="Biological Assembly",
            description="Translation vector",
        ).strip()
        try:
            translation = [float(x) for x in trans_str.split()]
            if len(translation) == 3:
                break
            else:
                print("    Please enter exactly 3 values.")
        except ValueError:
            print("    Invalid input. Please enter 3 numeric values.")

    return {
        'id': 1,
        'matrix': np.array(matrix_rows),
        'translation': np.array(translation)
    }


@register_module
class BiologicalAssemblyGenerator(ProcessingModule):
    """Module for generating biological assemblies from asymmetric units."""

    NAME = "Biological Assembly Generator"
    DESCRIPTION = "Generate biological assembly from crystal asymmetric unit"
    VERSION = "1.0.0"
    CATEGORY = "preparation"
    PRIORITY = 2  # After PDB Loader

    @property
    def console(self):
        """Get console from processor or create new one"""
        if self.processor and hasattr(self.processor, 'console'):
            return self.processor.console
        else:
            from rich.console import Console
            return Console()

    def get_menu_options(self) -> Dict[str, str]:
        """Get module menu options."""
        return {
            "generate": "Generate biological assembly from asymmetric unit",
        }

    def get_enhanced_menu_options(self, workspace):
        """Get menu options with enhanced status information."""
        from proprep.utils.enhanced_menu import MenuOption, OptionStatus
        from proprep.utils.structure_selector import StructureSelector

        options = []

        # Check if a structure is loaded
        selector = StructureSelector(workspace, self.console)
        has_structure = selector.get_structure(silent=True) is not None

        # Check if assembly already generated
        has_assembly = workspace.get("biological_assembly_pdb_file") is not None

        if has_assembly:
            status = OptionStatus.COMPLETED
            desc = "Generate biological assembly from asymmetric unit"
        elif has_structure:
            status = OptionStatus.AVAILABLE
            desc = "Generate biological assembly from asymmetric unit"
        else:
            status = OptionStatus.BLOCKED
            desc = "Generate biological assembly from asymmetric unit"

        options.append(MenuOption(
            key="1",
            description=desc,
            status=status,
            dependency_text="[Load a structure first] ○" if not has_structure else ""
        ))

        return options

    def handle_menu_option(self, option: str) -> bool:
        """Handle a menu option selection."""
        if option == "generate":
            return self.process(self.processor.workspace)
        return False

    def get_workspace_outputs(self) -> list:
        """Get workspace outputs"""
        return ["biological_assembly_pdb_file"]

    def availability_note(self, workspace):
        """Menu note when unavailable (○). Mirrors can_process exactly."""
        return None if self.can_process(workspace) else "Needs a loaded structure"

    def can_process(self, workspace) -> bool:
        """Check if a structure is available for assembly generation."""
        from proprep.utils.structure_selector import StructureSelector
        selector = StructureSelector(workspace, self.console)
        return selector.get_structure(silent=True) is not None

    def process(self, workspace) -> bool:
        """Generate biological assembly from loaded structure."""
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(workspace, self.console)

        # Get input structure
        result = selector.get_structure(silent=False, return_key=True)
        if not result:
            self.console.print("[red]No structure available in workspace.[/red]")
            return False

        pdb_path, source_key = result
        self.console.print(f"\n[cyan]Input structure:[/cyan] {Path(pdb_path).name}")
        self.console.print(f"[grey50]Source: {source_key}[/grey50]")

        # Extract REMARK 350 from PDB
        remark_lines = extract_remark_350_from_pdb(pdb_path)

        if remark_lines:
            assemblies = parse_remark_350(remark_lines)
            self.console.print(f"\n[green]Found {len(assemblies)} biological assembly definition(s) in REMARK 350[/green]")
        else:
            assemblies = {}
            self.console.print("\n[yellow]No REMARK 350 records found in PDB file.[/yellow]")

        # Display available assemblies
        if assemblies:
            self.console.print("\n[bold]Available biological assemblies:[/bold]")
            for num, data in sorted(assemblies.items()):
                desc = data['description'] or 'Unknown'
                author = " (author determined)" if data['author_determined'] else ""
                software = " (software determined)" if data['software_determined'] else ""
                chains = ', '.join(data['chains']) if data['chains'] else 'All'
                num_transforms = len(data['transforms'])

                self.console.print(
                    f"  [bright_blue]{num}.[/bright_blue] {desc}{author}{software} - "
                    f"{num_transforms} transform(s), chains: {chains}"
                )

        # Get user choice
        self.console.print("\n[bold]Options:[/bold]")
        if assemblies:
            self.console.print("  [cyan]1-N[/cyan]  Select assembly from REMARK 350")
        self.console.print("  [cyan]m[/cyan]    Enter manual transformation matrix")
        self.console.print("  [cyan]c[/cyan]    Cancel")

        valid_choices = ['m', 'c']
        if assemblies:
            valid_choices.extend([str(n) for n in assemblies.keys()])

        choice = prompt_with_context(
            self.processor,
            "\nSelect option",
            default='1' if assemblies else 'm',
            module=self.NAME,
            description="Select biological assembly"
        )

        if choice.lower() == 'c':
            self.console.print("[yellow]Cancelled.[/yellow]")
            return False

        # Determine assembly data to use
        if choice.lower() == 'm':
            # Manual matrix entry
            self.console.print("\n[bold]Manual transformation matrix entry[/bold]")

            # Get chains to transform
            chains_str = prompt_with_context(
                self.processor,
                "Chains to transform (comma-separated, or 'all')",
                default='all',
                module=self.NAME,
                description="Chains for transformation"
            )

            if chains_str.lower() == 'all':
                # Get all chains from PDB
                chains = self._get_all_chains(pdb_path)
            else:
                chains = [c.strip() for c in chains_str.split(',')]

            # Get transforms
            transforms = []
            while True:
                transform = create_manual_transform()
                transform['id'] = len(transforms) + 1
                transforms.append(transform)

                if not confirm_with_context(
                    self.processor,
                    "Add another transformation?",
                    default=False,
                    module=self.NAME,
                    description="Add more transforms"
                ):
                    break

            assembly_data = {
                'description': 'Manual',
                'author_determined': False,
                'software_determined': False,
                'chains': chains,
                'transforms': transforms
            }

        else:
            # Use REMARK 350 assembly
            try:
                assembly_num = int(choice)
                if assembly_num not in assemblies:
                    self.console.print(f"[red]Invalid assembly number: {assembly_num}[/red]")
                    return False
                assembly_data = assemblies[assembly_num]
            except ValueError:
                self.console.print(f"[red]Invalid choice: {choice}[/red]")
                return False

        # Options for assembly generation
        include_hetatm = confirm_with_context(
            self.processor,
            "Include ligands/cofactors (HETATM)?",
            default=True,
            module=self.NAME,
            description="Include HETATM records"
        )

        include_water = confirm_with_context(
            self.processor,
            "Include water molecules?",
            default=False,
            module=self.NAME,
            description="Include water molecules"
        )

        # Warn about large assemblies
        num_transforms = len(assembly_data['transforms'])
        estimated_chains = len(assembly_data['chains']) * num_transforms

        if estimated_chains > 62:
            self.console.print(
                f"\n[yellow]Warning: Assembly will have {estimated_chains} chains, "
                f"but PDB format only supports 62 unique chain IDs.[/yellow]"
            )
            if not confirm_with_context(
                self.processor,
                "Continue anyway (some chain IDs will be duplicated)?",
                default=False,
                module=self.NAME,
                description="Confirm large assembly"
            ):
                return False

        # Generate output filename
        input_name = Path(pdb_path).stem
        output_name = f"{input_name}_biological_assembly.pdb"
        output_path = Path(pdb_path).parent / output_name

        # Apply transformations
        self.console.print(f"\n[cyan]Generating biological assembly...[/cyan]")

        try:
            stats = apply_assembly_transforms(
                pdb_path=pdb_path,
                assembly_data=assembly_data,
                output_path=str(output_path),
                include_hetatm=include_hetatm,
                include_water=include_water
            )

            self.console.print(f"\n[green]Assembly generated successfully![/green]")
            self.console.print(f"  Original: {len(stats['original_chains'])} chains, {stats['original_atoms']:,} atoms")
            self.console.print(f"  Assembly: {len(stats['new_chains'])} chains ({', '.join(stats['new_chains'])}), {stats['new_atoms']:,} atoms")
            self.console.print(f"  Transforms applied: {stats['transforms_applied']}")
            self.console.print(f"\n[green]Saved:[/green] {output_path}")

            # Save to workspace via update_workspace so the coordinator's
            # structure-key hook fires and the viewer auto-swaps from
            # the input PDB to the freshly-generated assembly. A raw
            # workspace.set would bypass that notification.
            workspace = self.update_workspace(
                workspace, "biological_assembly_pdb_file", str(output_path),
            )
            self.console.print("[green]Registered in workspace as 'biological_assembly_pdb_file'[/green]")

            # Also cache the parsed Structure object so downstream consumers
            # (e.g. PDB Filter) analyze the full multi-chain assembly rather
            # than re-resolving to a different cached structure. Mirrors
            # structure_loader's biological_assembly_structure caching.
            # Best-effort: the file is already saved and its path registered,
            # so a parse failure here only means the structure is re-parsed
            # on demand (PDB Filter handles a missing cached object).
            try:
                from Bio.PDB import PDBParser
                structure = PDBParser(QUIET=True).get_structure(
                    f"{input_name}_assembly", str(output_path)
                )
                workspace = self.update_workspace(
                    workspace, "biological_assembly_structure", structure,
                )
            except Exception as e:
                logger.warning(f"Could not cache biological_assembly_structure: {e}")

            return True

        except Exception as e:
            self.console.print(f"[red]Error generating assembly: {e}[/red]")
            logger.exception("Error in biological assembly generation")
            return False

    def _get_all_chains(self, pdb_path: str) -> List[str]:
        """Get all chain IDs from a PDB file."""
        chains = set()
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    if len(line) > 21:
                        chain_id = line[21]
                        if chain_id.strip():
                            chains.add(chain_id)
        return sorted(chains)
