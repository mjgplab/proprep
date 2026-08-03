"""
Pre-frcmod Generator
Generates preliminary frcmod files with NON/YES markers following MCPB workflow.

© 2024 ProPrep Developer. All rights reserved.

This module replicates MCPB's gene_pre_frcmod_file.py functionality, creating
an intermediate frcmod file that:
1. Marks metal-related parameters as "NON" (need QM calculation)
2. Marks non-metal parameters as "YES" (from force field)
3. Provides educational checkpoint before expensive QM calculations
4. Enables workflow resumption if QM calculations fail

The pre-frcmod file serves as:
- Template for final frcmod generation
- Documentation of which parameters need QM
- Checkpoint for workflow management

PROPRIETARY SOFTWARE: This file contains proprietary code belonging to the ProPrep project.
Unauthorized copying, distribution, or modification is strictly prohibited.
"""

from pathlib import Path
from typing import Dict, List, Tuple, Set, Optional, TYPE_CHECKING
from dataclasses import dataclass
from rich.console import Console

from proprep.forcefield_prep.ff_types import NonbondedParameter

if TYPE_CHECKING:
    from proprep.forcefield_prep.forcefield_data import ForceFieldData


@dataclass
class ParameterNeed:
    """Tracks whether a parameter needs QM calculation or is in force field."""
    atom_types: Tuple[str, ...]
    status: str  # "NON" or "YES"
    value: Optional[str] = None  # Parameter values if from FF
    comment: str = ""


class PreFrcmodGenerator:
    """
    Generate pre-frcmod files with NON/YES markers.

    Follows MCPB gene_pre_frcmod_file.py logic:
    - MASS: All atoms (YES markers)
    - BOND: NON for metal-ligand, YES for others
    - ANGL: NON for metal-containing, YES for others
    - DIHE: All metal dihedrals set to zero (YES markers)
    - NONB: All atoms (YES markers)

    The "NON" markers indicate parameters that require QM calculations.
    The "YES" markers indicate parameters available from the force field.
    """

    def __init__(self, param_provider, console: Optional[Console] = None):
        """
        Initialize pre-frcmod generator.

        Args:
            param_provider: Parameter provider with get_mass_parameter, get_bond_parameter,
                          get_angle_parameter, get_nonbonded_parameter methods.
                          Can be ForceFieldData or PrmtopParameterProvider.
            console: Rich console for output (optional)
        """
        self.ff_data = param_provider  # Keep attribute name for compatibility
        self.console = console or Console()
        self.warnings = []  # Collect warnings for missing parameters

    def generate_pre_frcmod(
        self,
        type_assignments: Dict,
        bonds: List[Tuple[int, int]],
        angles: List[Tuple[int, int, int]],
        output_file: Path,
        idx_to_coords: Dict[int, Tuple[float, float, float]]
    ) -> Dict:
        """
        Generate pre-frcmod file with NON/YES markers.

        Args:
            type_assignments: Dict mapping coords to atom type data
            bonds: List of (atom1_idx, atom2_idx) tuples
            angles: List of (atom1_idx, atom2_idx, atom3_idx) tuples
            output_file: Path for pre-frcmod file
            idx_to_coords: Mapping from atom indices to coordinates

        Returns:
            Dict with statistics about NON/YES parameters
        """
        self.console.print("[bold]Generating pre-frcmod file...[/bold]")
        self.warnings = []  # Reset warnings

        # Identify metal centers and ligands
        metal_coords = {
            coords for coords, data in type_assignments.items()
            if data.get('is_center', False)
        }

        ligand_coords = {
            coords for coords, data in type_assignments.items()
            if data.get('is_metal_ligand', False)
        }

        # Track statistics
        stats = {
            'non_bonds': [],
            'yes_bonds': [],
            'non_angles': [],
            'yes_angles': [],
            'mass_entries': 0,
            'nonb_entries': 0,
            'dihedral_entries': 0
        }

        # Open output file
        with open(output_file, 'w') as f:
            self._write_header(f)

            # MASS section
            self._write_mass_section(f, type_assignments, stats)

            # BOND section
            self._write_bond_section(
                f, bonds, type_assignments, idx_to_coords,
                metal_coords, ligand_coords, stats
            )

            # ANGL section
            self._write_angle_section(
                f, angles, type_assignments, idx_to_coords,
                metal_coords, stats
            )

            # DIHE section
            self._write_dihedral_section(
                f, bonds, angles, type_assignments, idx_to_coords,
                metal_coords, stats
            )

            # IMPR section (empty for MCPB)
            f.write("IMPR\n\n")

            # NONB section
            self._write_nonbonded_section(f, type_assignments, stats)

            # End with blank lines
            f.write("\n\n")

        # Print summary
        self.console.print(f"\n[green]✅ Pre-frcmod file generated: {output_file.name}[/green]")
        self.console.print(f"[cyan]Parameters requiring QM calculations:[/cyan]")
        self.console.print(f"  • NON bonds: {len(stats['non_bonds'])}")
        self.console.print(f"  • NON angles: {len(stats['non_angles'])}")
        self.console.print(f"[cyan]Parameters from force field:[/cyan]")
        self.console.print(f"  • YES bonds: {len(stats['yes_bonds'])}")
        self.console.print(f"  • YES angles: {len(stats['yes_angles'])}")

        # Display warnings if any
        if self.warnings:
            self.console.print(f"\n[yellow]⚠️  {len(self.warnings)} warning(s):[/yellow]")
            for warning in self.warnings:
                self.console.print(f"  [yellow]• {warning}[/yellow]")
            self.console.print("\n[yellow]Please review the pre-frcmod file and verify fallback values.[/yellow]")

        stats['warnings'] = self.warnings
        return stats

    def _write_header(self, f):
        """Write frcmod header."""
        f.write("REMARK GOES HERE, THIS FILE IS GENERATED BY PROPREP MCPB\n")
        f.write("REMARK Pre-frcmod: NON markers indicate parameters needing QM calculations\n")

    def _write_mass_section(self, f, type_assignments: Dict, stats: Dict):
        """Write MASS section with all atom types."""
        f.write("MASS\n")

        # Collect unique renamed atom types
        renamed_types = {}
        for coords, data in type_assignments.items():
            if data.get('renamed', False):
                atom_type = data['renamed_type']
                if atom_type not in renamed_types:
                    original_type = data['original_type']
                    element = data.get('element', 'UNK')

                    # Priority 1: Use pre-stored mass (metals from MetalIonDatabase)
                    if data.get('mass') is not None:
                        mass = data['mass']
                        comment = f"{element} from MetalIonDatabase"
                    else:
                        # Priority 2: Get mass from prmtop for original type
                        mass_param = self.ff_data.get_mass_parameter(original_type)
                        if mass_param:
                            mass = mass_param.mass
                            comment = mass_param.comment if mass_param.comment else original_type
                        else:
                            # Priority 3: Use pymsmt element table (same as MCPB.py)
                            mass = self._get_element_mass(element)
                            comment = f"{element} - VERIFY MASS (not in prmtop)"
                            self.warnings.append(
                                f"Mass for {atom_type} (original: {original_type}) not found. "
                                f"Using element mass {mass:.2f} for {element}."
                            )

                    renamed_types[atom_type] = (mass, comment)

        # Write mass entries
        for atom_type in sorted(renamed_types.keys()):
            mass, comment = renamed_types[atom_type]
            comment_str = f"  {comment}" if comment else ""
            f.write(f"YES {atom_type:2s}  {mass:8.2f}{comment_str}\n")
            stats['mass_entries'] += 1

        f.write("\n")

    def _write_bond_section(
        self, f, bonds: List[Tuple[int, int]],
        type_assignments: Dict, idx_to_coords: Dict,
        metal_coords: Set, ligand_coords: Set, stats: Dict
    ):
        """Write BOND section with NON/YES markers."""
        f.write("BOND\n")

        # Collect bond types
        bond_types = {}

        for atom1_idx, atom2_idx in bonds:
            coord1 = idx_to_coords[atom1_idx]
            coord2 = idx_to_coords[atom2_idx]

            data1 = type_assignments.get(coord1)
            data2 = type_assignments.get(coord2)

            if not data1 or not data2:
                continue

            type1 = data1['renamed_type']
            type2 = data2['renamed_type']

            # Create canonical bond type (alphabetical order)
            if type1 <= type2:
                bond_type = (type1, type2)
            else:
                bond_type = (type2, type1)

            # Skip if already processed
            if bond_type in bond_types:
                continue

            # Check if metal-ligand bond
            is_metal_bond = (coord1 in metal_coords or coord2 in metal_coords)

            if is_metal_bond:
                # Metal-ligand bond needs QM calculation
                bond_types[bond_type] = ('NON', None, '')
                stats['non_bonds'].append(f"{type1}-{type2}")
            else:
                # Try to find in force field (unified ForceFieldData)
                orig_type1 = data1['original_type']
                orig_type2 = data2['original_type']

                bond_param = self.ff_data.get_bond_parameter(orig_type1, orig_type2)

                if bond_param:
                    value = f"  {bond_param.force_constant:6.1f}   {bond_param.eq_length:7.4f}"
                    comment = f"    {bond_param.source}" if bond_param.source else ""
                    bond_types[bond_type] = ('YES', value, comment)
                    stats['yes_bonds'].append(f"{type1}-{type2}")
                else:
                    # Not in FF and not metal - mark as NON (will need QM or manual)
                    bond_types[bond_type] = ('NON', None, '')
                    stats['non_bonds'].append(f"{type1}-{type2}")
                    self.warnings.append(
                        f"Bond {type1}-{type2} (original: {orig_type1}-{orig_type2}) "
                        f"not found in force field."
                    )

        # Write bond parameters (NON first, then YES)
        for bond_type, (marker, value, comment) in sorted(bond_types.items()):
            if marker == 'NON':
                f.write(f"NON {bond_type[0]:2s}-{bond_type[1]:2s}\n")

        for bond_type, (marker, value, comment) in sorted(bond_types.items()):
            if marker == 'YES':
                f.write(f"YES {bond_type[0]:2s}-{bond_type[1]:2s}{value}{comment}\n")

        f.write("\n")

    def _write_angle_section(
        self, f, angles: List[Tuple[int, int, int]],
        type_assignments: Dict, idx_to_coords: Dict,
        metal_coords: Set, stats: Dict
    ):
        """Write ANGL section with NON/YES markers."""
        f.write("ANGL\n")

        # Collect angle types
        angle_types = {}

        for atom1_idx, atom2_idx, atom3_idx in angles:
            coord1 = idx_to_coords[atom1_idx]
            coord2 = idx_to_coords[atom2_idx]
            coord3 = idx_to_coords[atom3_idx]

            data1 = type_assignments.get(coord1)
            data2 = type_assignments.get(coord2)
            data3 = type_assignments.get(coord3)

            if not data1 or not data2 or not data3:
                continue

            type1 = data1['renamed_type']
            type2 = data2['renamed_type']
            type3 = data3['renamed_type']

            # Create canonical angle type (can reverse 1-2-3 to 3-2-1)
            if type1 <= type3:
                angle_type = (type1, type2, type3)
            else:
                angle_type = (type3, type2, type1)

            # Skip if already processed
            if angle_type in angle_types:
                continue

            # Check if metal-containing angle
            is_metal_angle = (
                coord1 in metal_coords or
                coord2 in metal_coords or
                coord3 in metal_coords
            )

            if is_metal_angle:
                # Metal-containing angle needs QM calculation
                angle_types[angle_type] = ('NON', None, '')
                stats['non_angles'].append(f"{type1}-{type2}-{type3}")
            else:
                # Try to find in force field (unified ForceFieldData)
                orig_type1 = data1['original_type']
                orig_type2 = data2['original_type']
                orig_type3 = data3['original_type']

                angle_param = self.ff_data.get_angle_parameter(orig_type1, orig_type2, orig_type3)

                if angle_param:
                    value = f"  {angle_param.force_constant:7.2f}   {angle_param.eq_angle:7.2f}"
                    comment = f"    {angle_param.source}" if angle_param.source else ""
                    angle_types[angle_type] = ('YES', value, comment)
                    stats['yes_angles'].append(f"{type1}-{type2}-{type3}")
                else:
                    # Not in FF and not metal - mark as NON
                    angle_types[angle_type] = ('NON', None, '')
                    stats['non_angles'].append(f"{type1}-{type2}-{type3}")
                    self.warnings.append(
                        f"Angle {type1}-{type2}-{type3} (original: {orig_type1}-{orig_type2}-{orig_type3}) "
                        f"not found in force field."
                    )

        # Write angle parameters (NON first, then YES)
        for angle_type, (marker, value, comment) in sorted(angle_types.items()):
            if marker == 'NON':
                f.write(f"NON {angle_type[0]:2s}-{angle_type[1]:2s}-{angle_type[2]:2s}\n")

        for angle_type, (marker, value, comment) in sorted(angle_types.items()):
            if marker == 'YES':
                f.write(f"YES {angle_type[0]:2s}-{angle_type[1]:2s}-{angle_type[2]:2s}{value}{comment}\n")

        f.write("\n")

    def _write_dihedral_section(
        self, f, bonds: List[Tuple[int, int]],
        angles: List[Tuple[int, int, int]],
        type_assignments: Dict, idx_to_coords: Dict,
        metal_coords: Set, stats: Dict
    ):
        """
        Write DIHE section.

        Following MCPB convention: All dihedrals involving metals are set to zero.
        """
        f.write("DIHE\n")

        # Build connectivity for dihedral generation
        connectivity = {}
        for atom1_idx, atom2_idx in bonds:
            if atom1_idx not in connectivity:
                connectivity[atom1_idx] = []
            if atom2_idx not in connectivity:
                connectivity[atom2_idx] = []
            connectivity[atom1_idx].append(atom2_idx)
            connectivity[atom2_idx].append(atom1_idx)

        # Generate dihedrals from angles (extend by one bond on each end)
        dihedral_types = set()

        for atom1_idx, atom2_idx, atom3_idx in angles:
            # Find atoms bonded to atom1 (excluding atom2)
            if atom1_idx in connectivity:
                for atom0_idx in connectivity[atom1_idx]:
                    if atom0_idx != atom2_idx:
                        # Dihedral: atom0_idx - atom1_idx - atom2_idx - atom3_idx
                        dihedral_indices = (atom0_idx, atom1_idx, atom2_idx, atom3_idx)
                        self._process_dihedral(
                            dihedral_indices, type_assignments, idx_to_coords,
                            metal_coords, dihedral_types
                        )

            # Find atoms bonded to atom3 (excluding atom2)
            if atom3_idx in connectivity:
                for atom4_idx in connectivity[atom3_idx]:
                    if atom4_idx != atom2_idx:
                        # Dihedral: atom1_idx - atom2_idx - atom3_idx - atom4_idx
                        dihedral_indices = (atom1_idx, atom2_idx, atom3_idx, atom4_idx)
                        self._process_dihedral(
                            dihedral_indices, type_assignments, idx_to_coords,
                            metal_coords, dihedral_types
                        )

        # Write dihedral parameters (all zero for metal dihedrals)
        for dihedral_type in sorted(dihedral_types):
            type_str = f"{dihedral_type[0]:2s}-{dihedral_type[1]:2s}-{dihedral_type[2]:2s}-{dihedral_type[3]:2s}"
            f.write(f"YES {type_str}    3       0.00       0.00   3.0    Treat as zero by MCPB.py\n")
            stats['dihedral_entries'] += 1

        f.write("\n")

    def _process_dihedral(
        self, dihedral_indices: Tuple[int, int, int, int],
        type_assignments: Dict, idx_to_coords: Dict,
        metal_coords: Set, dihedral_types: Set
    ):
        """Process a single dihedral and add to set if metal-containing."""
        coords = [idx_to_coords[idx] for idx in dihedral_indices]
        data = [type_assignments.get(c) for c in coords]

        if any(d is None for d in data):
            return

        # Check if any atom is a metal
        has_metal = any(c in metal_coords for c in coords)

        if has_metal:
            # Extract atom types
            types = tuple(d['renamed_type'] for d in data)

            # Add canonical form (also reverse)
            dihedral_types.add(types)
            dihedral_types.add(types[::-1])

    def _write_nonbonded_section(self, f, type_assignments: Dict, stats: Dict):
        """Write NONB section with VDW parameters."""
        f.write("NONB\n")

        # Collect unique renamed atom types
        renamed_types = {}
        for coords, data in type_assignments.items():
            if data.get('renamed', False):
                atom_type = data['renamed_type']
                if atom_type not in renamed_types:
                    # Check if VDW params are pre-stored (metals from MetalIonDatabase)
                    vdw_params = data.get('vdw_params')
                    if vdw_params:
                        renamed_types[atom_type] = NonbondedParameter(
                            atom_type=atom_type,
                            radius=vdw_params.get('radius', 0.0),
                            well_depth=vdw_params.get('epsilon', 0.0),
                            comment=vdw_params.get('source', 'Metal ion database')
                        )
                    else:
                        # Get nonbonded param from force field for original type
                        original_type = data['original_type']
                        nonb_param = self.ff_data.get_nonbonded_parameter(original_type)

                        if nonb_param:
                            renamed_types[atom_type] = nonb_param
                        else:
                            # FALLBACK: Use default values with warning
                            element = data.get('element', 'UNK')
                            radius, epsilon = self._get_element_vdw(element)
                            renamed_types[atom_type] = NonbondedParameter(
                                atom_type=atom_type,
                                radius=radius,
                                well_depth=epsilon,
                                comment=f"{element} - VERIFY VDW PARAMETERS (not in force field)"
                            )
                            self.warnings.append(
                                f"VDW parameters for {atom_type} (original: {original_type}) not found in force field. "
                                f"Using fallback values for element {element}."
                            )

        # Write nonbonded entries
        for atom_type in sorted(renamed_types.keys()):
            param = renamed_types[atom_type]
            comment_str = f"       {param.comment}" if param.comment else ""
            f.write(f"YES   {atom_type:2s}  {param.radius:8.4f}  {param.well_depth:12.6f}{comment_str}\n")
            stats['nonb_entries'] += 1

        f.write("\n")

    def _get_element_mass(self, element: str) -> float:
        """
        Get atomic mass for element.

        Uses pymsmt.mol.element.Mass (same as MCPB.py) for consistency.
        """
        from pymsmt.mol.element import Mass
        # Normalize: 'ZN' -> 'Zn', 'FE' -> 'Fe' (pymsmt uses title case)
        element_key = element.strip().capitalize()
        return Mass.get(element_key, 1.0)

    def _get_element_vdw(self, element: str) -> Tuple[float, float]:
        """
        Get VDW parameters (radius, epsilon) for element.

        Warning: This is a fallback. Ideally all VDW params should come from force field.
        Returns: (radius in Å, epsilon in kcal/mol)
        """
        # Rough estimates - user should verify these!
        vdw_table = {
            'H': (1.20, 0.0157),
            'C': (1.70, 0.0860),
            'N': (1.55, 0.1700),
            'O': (1.52, 0.2100),
            'S': (1.80, 0.2500),
            'P': (1.80, 0.2000),
            'Zn': (1.39, 0.0150),  # Common for Zn2+
            'Fe': (1.30, 0.0100),
            'Cu': (1.40, 0.0200),
            'Mg': (1.40, 0.0100),
            'Ca': (1.60, 0.0150),
        }
        return vdw_table.get(element, (1.50, 0.0100))
