"""
Seminario Method Refinement for Small Molecule Parameters

This module provides QM-based refinement of force field parameters using the
Seminario method, which derives force constants directly from the Hessian matrix
of a frequency calculation.

The Seminario method is an alternative to paramfit that:
1. Uses the Hessian (second derivatives of energy) from a single optimized structure
2. Does NOT require conformer generation or multiple QM calculations
3. Directly computes force constants by projecting Hessian eigenvalues onto bond/angle directions

Advantages over paramfit:
- Uses data already computed during opt+freq calculation
- No need for CREST conformer generation
- No iterative fitting - direct analytical result
- Well-suited for molecules where opt+freq already performed

Reference:
    J. M. Seminario, "Calculation of intramolecular force fields from
    second-derivative tensors", International Journal of Quantum Chemistry,
    1996, 60(7), 1271-1277.

Author: ProPrep Development Team
"""

import os
import subprocess
import shutil
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import re

import numpy as np
from rich.console import Console
from rich.panel import Panel
from rich.table import Table


def convert_chk_to_fchk(chk_file: str, console: Console) -> Optional[str]:
    """
    Convert Gaussian binary checkpoint (.chk) to formatted checkpoint (.fchk).

    Args:
        chk_file: Path to the .chk file
        console: Rich console for output

    Returns:
        Path to the generated .fchk file, or None if conversion failed
    """
    if not os.path.exists(chk_file):
        console.print(f"[red]Checkpoint file not found: {chk_file}[/red]")
        return None

    # Generate fchk filename
    base = os.path.splitext(chk_file)[0]
    fchk_file = f"{base}.fchk"

    # Check if fchk already exists and is newer than chk
    if os.path.exists(fchk_file):
        chk_mtime = os.path.getmtime(chk_file)
        fchk_mtime = os.path.getmtime(fchk_file)
        if fchk_mtime >= chk_mtime:
            console.print(f"[cyan]Using existing formatted checkpoint: {fchk_file}[/cyan]")
            return fchk_file

    console.print(f"[cyan]Converting checkpoint to formatted checkpoint...[/cyan]")
    console.print(f"[grey50]formchk {chk_file} {fchk_file}[/grey50]")

    try:
        result = subprocess.run(
            ["formchk", chk_file, fchk_file],
            capture_output=True,
            text=True,
            timeout=60
        )

        if result.returncode != 0:
            console.print(f"[red]formchk failed: {result.stderr}[/red]")
            return None

        console.print(f"[green]✓ Generated formatted checkpoint: {fchk_file}[/green]")
        return fchk_file

    except FileNotFoundError:
        console.print("[red]formchk not found. Ensure Gaussian is installed and in PATH.[/red]")
        return None
    except subprocess.TimeoutExpired:
        console.print("[red]formchk timed out[/red]")
        return None
    except Exception as e:
        console.print(f"[red]Error running formchk: {e}[/red]")
        return None


def parse_mol2_connectivity(mol2_file: str) -> Dict:
    """
    Parse MOL2 file to extract atom indices, types, and connectivity.

    Args:
        mol2_file: Path to MOL2 file

    Returns:
        Dictionary containing:
        - atom_types: List of GAFF2 atom types by index
        - atom_names: List of atom names by index
        - bonds: List of (idx1, idx2) tuples (0-indexed)
        - elements: List of element symbols
    """
    atom_types = []
    atom_names = []
    elements = []
    bonds = []

    with open(mol2_file, 'r') as f:
        lines = f.readlines()

    section = None
    for line in lines:
        stripped = line.strip()

        # Section headers
        if stripped.startswith('@<TRIPOS>'):
            section = stripped.replace('@<TRIPOS>', '')
            continue

        if section == 'ATOM' and stripped:
            parts = line.split()
            if len(parts) >= 6:
                # MOL2 format: idx name x y z type [subst] [subst_name] [charge]
                atom_name = parts[1]
                atom_type = parts[5]  # GAFF2 type

                # Extract element from atom name (first 1-2 uppercase letters)
                element_match = re.match(r'([A-Z][a-z]?)', atom_name)
                element = element_match.group(1) if element_match else atom_name[0]

                atom_names.append(atom_name)
                atom_types.append(atom_type)
                elements.append(element)

        elif section == 'BOND' and stripped:
            parts = line.split()
            if len(parts) >= 3:
                try:
                    idx1 = int(parts[1]) - 1  # Convert to 0-indexed
                    idx2 = int(parts[2]) - 1
                    bonds.append((idx1, idx2))
                except ValueError:
                    continue

    return {
        'atom_types': atom_types,
        'atom_names': atom_names,
        'elements': elements,
        'bonds': bonds,
        'n_atoms': len(atom_types)
    }


def find_angles_from_bonds(bonds: List[Tuple[int, int]]) -> List[Tuple[int, int, int]]:
    """
    Derive all angles from bond connectivity.

    An angle A-B-C exists when atoms A and C are both bonded to B.

    Args:
        bonds: List of (idx1, idx2) bond tuples

    Returns:
        List of (A, B, C) angle tuples where B is the central atom
    """
    # Build adjacency list
    adjacency = {}
    for idx1, idx2 in bonds:
        if idx1 not in adjacency:
            adjacency[idx1] = []
        if idx2 not in adjacency:
            adjacency[idx2] = []
        adjacency[idx1].append(idx2)
        adjacency[idx2].append(idx1)

    # Find all angles
    angles = []
    for center in adjacency:
        neighbors = adjacency[center]
        if len(neighbors) >= 2:
            # All pairs of neighbors form angles through center
            for i in range(len(neighbors)):
                for j in range(i + 1, len(neighbors)):
                    # Order so smaller index comes first
                    a, c = sorted([neighbors[i], neighbors[j]])
                    angles.append((a, center, c))

    return angles


def map_penalty_to_indices(
    selected_params: List[Tuple[str, float, str, str]],
    mol2_connectivity: Dict,
    console: Console
) -> Tuple[List[Tuple[int, int, str, str]], List[Tuple[int, int, int, str, str, str]]]:
    """
    Map penalty parameter atom types to actual atom indices.

    Args:
        selected_params: List of (param_name, score, status, section) tuples
        mol2_connectivity: Dictionary from parse_mol2_connectivity
        console: Rich console for output

    Returns:
        Tuple of (bonds_to_refine, angles_to_refine):
        - bonds_to_refine: List of (idx1, idx2, type1, type2)
        - angles_to_refine: List of (idx1, idx2, idx3, type1, type2, type3)
    """
    atom_types = mol2_connectivity['atom_types']
    bonds = mol2_connectivity['bonds']
    angles = find_angles_from_bonds(bonds)

    bonds_to_refine = []
    angles_to_refine = []

    for param_name, score, status, section in selected_params:
        # Parse atom types from parameter name
        types = param_name.replace(' ', '').split('-')

        if section == 'BOND' and len(types) == 2:
            t1, t2 = types
            # Find matching bonds
            for idx1, idx2 in bonds:
                at1, at2 = atom_types[idx1], atom_types[idx2]
                if (at1 == t1 and at2 == t2) or (at1 == t2 and at2 == t1):
                    bonds_to_refine.append((idx1, idx2, at1, at2))
                    console.print(f"[grey50]  Mapped BOND {param_name} → atoms {idx1+1}-{idx2+1}[/grey50]")
                    break
            else:
                console.print(f"[yellow]  Could not find bond for {param_name}[/yellow]")

        elif section == 'ANGLE' and len(types) == 3:
            t1, t2, t3 = types
            # Find matching angles (t2 is central atom)
            for idx1, idx2, idx3 in angles:
                at1, at2, at3 = atom_types[idx1], atom_types[idx2], atom_types[idx3]
                # Check both orderings (A-B-C or C-B-A)
                if at2 == t2:
                    if (at1 == t1 and at3 == t3) or (at1 == t3 and at3 == t1):
                        angles_to_refine.append((idx1, idx2, idx3, at1, at2, at3))
                        console.print(f"[grey50]  Mapped ANGLE {param_name} → atoms {idx1+1}-{idx2+1}-{idx3+1}[/grey50]")
                        break
            else:
                console.print(f"[yellow]  Could not find angle for {param_name}[/yellow]")

        elif section == 'DIHE':
            console.print(f"[grey50]  Skipping DIHE {param_name} (Seminario does not compute dihedrals)[/grey50]")

    return bonds_to_refine, angles_to_refine


def run_seminario_refinement(
    fchk_file: str,
    mol2_file: str,
    frcmod_file: str,
    selected_params: List[Tuple[str, float, str, str]],
    console: Console,
    scale_factor: float = 1.0
) -> Dict[str, Any]:
    """
    Apply Seminario method to refine force constants.

    Args:
        fchk_file: Path to Gaussian formatted checkpoint file
        mol2_file: Path to MOL2 file with atom types
        frcmod_file: Path to original frcmod file
        selected_params: List of (param_name, score, status, section) tuples
        console: Rich console for output
        scale_factor: Frequency scaling factor (default 1.0)

    Returns:
        Dictionary with:
        - success: True if refinement succeeded
        - fitted_frcmod: Path to refined frcmod file
        - bond_params: List of refined bond parameters
        - angle_params: List of refined angle parameters
        - message: Status message
    """
    results = {
        'success': False,
        'fitted_frcmod': None,
        'bond_params': [],
        'angle_params': [],
        'message': ''
    }

    try:
        # Import Seminario components
        from proprep.forcefield_prep.mcpb.hessian_parser import HessianParser, HessianParseError
        from proprep.forcefield_prep.mcpb.seminario import SeminarioMethod

        console.print(f"\n[bold cyan]Seminario Method Parameter Refinement[/bold cyan]")
        console.print(f"[grey50]Reference: Seminario, Int. J. Quantum Chem. 1996, 60, 1271-1277[/grey50]")

        # Check that fchk_file is valid
        if not fchk_file:
            results['message'] = "No formatted checkpoint file available. Ensure formchk is in PATH."
            console.print(f"[red]{results['message']}[/red]")
            return results

        # Parse Hessian from fchk file
        console.print(f"\n[cyan]Parsing Hessian from {os.path.basename(fchk_file)}...[/cyan]")
        try:
            hessian, coords = HessianParser.parse_gaussian_fchk(Path(fchk_file))
            console.print(f"[green]✓ Hessian extracted ({hessian.shape[0]}×{hessian.shape[1]})[/green]")
        except HessianParseError as e:
            results['message'] = f"Failed to parse Hessian: {e}"
            console.print(f"[red]{results['message']}[/red]")
            return results

        # Validate Hessian
        is_valid, msg = HessianParser.validate_hessian(hessian, coords)
        if not is_valid:
            results['message'] = f"Invalid Hessian: {msg}"
            console.print(f"[red]{results['message']}[/red]")
            return results

        # Parse MOL2 for connectivity and atom types
        console.print(f"\n[cyan]Parsing connectivity from {os.path.basename(mol2_file)}...[/cyan]")
        mol2_data = parse_mol2_connectivity(mol2_file)
        console.print(f"[green]✓ Found {mol2_data['n_atoms']} atoms, {len(mol2_data['bonds'])} bonds[/green]")

        # Verify atom counts match
        n_atoms_hessian = hessian.shape[0] // 3
        if n_atoms_hessian != mol2_data['n_atoms']:
            results['message'] = f"Atom count mismatch: Hessian has {n_atoms_hessian}, MOL2 has {mol2_data['n_atoms']}"
            console.print(f"[red]{results['message']}[/red]")
            return results

        # Map selected parameters to atom indices
        console.print(f"\n[cyan]Mapping parameters to atom indices...[/cyan]")
        bonds_to_refine, angles_to_refine = map_penalty_to_indices(
            selected_params, mol2_data, console
        )

        if not bonds_to_refine and not angles_to_refine:
            console.print("[yellow]No bonds or angles to refine with Seminario method[/yellow]")
            # Check if all were dihedrals
            dihe_count = sum(1 for _, _, _, section in selected_params if section == 'DIHE')
            if dihe_count > 0:
                console.print(f"[grey50]({dihe_count} dihedral(s) skipped - Seminario only handles bonds/angles)[/grey50]")
            results['message'] = "No bonds or angles selected for refinement"
            return results

        # Initialize Seminario method
        console.print(f"\n[cyan]Computing force constants with Seminario method...[/cyan]")
        seminario = SeminarioMethod(hessian, coords, scale_factor=scale_factor, console=console)

        # Compute bond force constants
        bond_results = []
        for idx1, idx2, type1, type2 in bonds_to_refine:
            param = seminario.compute_bond_force_constant(
                idx1, idx2, type1, type2
            )
            bond_results.append(param)
            console.print(
                f"  [green]BOND {type1}-{type2}:[/green] "
                f"Kr = {param.force_constant:.2f} kcal/(mol·Å²), "
                f"Req = {param.eq_length:.4f} Å"
            )

        # Compute angle force constants
        angle_results = []
        for idx1, idx2, idx3, type1, type2, type3 in angles_to_refine:
            param = seminario.compute_angle_force_constant(
                idx1, idx2, idx3, type1, type2, type3
            )
            angle_results.append(param)
            console.print(
                f"  [green]ANGLE {type1}-{type2}-{type3}:[/green] "
                f"Kθ = {param.force_constant:.2f} kcal/(mol·rad²), "
                f"θeq = {param.eq_angle:.2f}°"
            )

        # Generate refined frcmod file
        output_frcmod = generate_seminario_frcmod(
            frcmod_file, bond_results, angle_results, console
        )

        if output_frcmod:
            results['success'] = True
            results['fitted_frcmod'] = output_frcmod
            results['bond_params'] = bond_results
            results['angle_params'] = angle_results
            results['message'] = f"Successfully refined {len(bond_results)} bond(s) and {len(angle_results)} angle(s)"
            console.print(f"\n[green]✓ {results['message']}[/green]")
        else:
            results['message'] = "Failed to generate refined frcmod file"

        return results

    except ImportError as e:
        results['message'] = f"Import error: {e}"
        console.print(f"[red]{results['message']}[/red]")
        return results
    except Exception as e:
        results['message'] = f"Error during Seminario refinement: {e}"
        console.print(f"[red]{results['message']}[/red]")
        import traceback
        console.print(f"[grey50]{traceback.format_exc()}[/grey50]")
        return results


def generate_seminario_frcmod(
    original_frcmod: str,
    bond_params: List,
    angle_params: List,
    console: Console
) -> Optional[str]:
    """
    Generate updated frcmod file with Seminario-derived parameters.

    Args:
        original_frcmod: Path to original frcmod file
        bond_params: List of BondParameter objects
        angle_params: List of AngleParameter objects
        console: Rich console for output

    Returns:
        Path to generated frcmod file, or None on failure
    """
    if not os.path.exists(original_frcmod):
        console.print(f"[red]Original frcmod not found: {original_frcmod}[/red]")
        return None

    # Read original frcmod
    with open(original_frcmod, 'r') as f:
        lines = f.readlines()

    # Create mapping of parameters to refine
    bond_map = {}
    for param in bond_params:
        # Normalize key (alphabetical order for consistency)
        types = sorted([param.atom1_type, param.atom2_type])
        key = f"{types[0]}-{types[1]}"
        bond_map[key] = param

    angle_map = {}
    for param in angle_params:
        # Key format: type1-type2-type3 where type2 is central
        key = f"{param.atom1_type}-{param.atom2_type}-{param.atom3_type}"
        angle_map[key] = param
        # Also add reverse
        rev_key = f"{param.atom3_type}-{param.atom2_type}-{param.atom1_type}"
        angle_map[rev_key] = param

    # Process frcmod and update parameters
    output_lines = []
    current_section = None
    updated_bonds = 0
    updated_angles = 0

    for line in lines:
        stripped = line.strip()

        # Track sections
        if stripped in ['MASS', 'BOND', 'ANGLE', 'DIHE', 'IMPROPER', 'NONBON', 'END']:
            current_section = stripped
            output_lines.append(line)
            continue

        # Process BOND section
        if current_section == 'BOND' and stripped and not stripped.startswith('Remark'):
            # Parse bond line (fixed width: 2+1+2 = 5 chars for atom types)
            if len(line) >= 5:
                atom1 = line[0:2].strip()
                atom2 = line[3:5].strip()
                types = sorted([atom1, atom2])
                key = f"{types[0]}-{types[1]}"

                if key in bond_map:
                    param = bond_map[key]
                    # Generate new line with Seminario values
                    # Format: A2-A2  force  length   comment
                    new_line = f"{atom1}-{atom2}  {param.force_constant:8.2f}   {param.eq_length:.4f}       Seminario method\n"
                    output_lines.append(new_line)
                    updated_bonds += 1
                    continue

        # Process ANGLE section
        if current_section == 'ANGLE' and stripped and not stripped.startswith('Remark'):
            # Parse angle line (fixed width: 2+1+2+1+2 = 8 chars for atom types)
            if len(line) >= 8:
                atom1 = line[0:2].strip()
                atom2 = line[3:5].strip()
                atom3 = line[6:8].strip()
                key = f"{atom1}-{atom2}-{atom3}"

                if key in angle_map:
                    param = angle_map[key]
                    # Generate new line with Seminario values
                    # Format: A2-A2-A2  force  angle   comment
                    new_line = f"{atom1}-{atom2}-{atom3}   {param.force_constant:8.3f}   {param.eq_angle:8.3f}     Seminario method\n"
                    output_lines.append(new_line)
                    updated_angles += 1
                    continue

        # Keep original line
        output_lines.append(line)

    # Write output file
    base_name = os.path.splitext(original_frcmod)[0]
    output_file = f"{base_name}_seminario.frcmod"

    with open(output_file, 'w') as f:
        f.writelines(output_lines)

    console.print(f"\n[cyan]Generated refined frcmod: {output_file}[/cyan]")
    console.print(f"  Updated {updated_bonds} bond parameter(s)")
    console.print(f"  Updated {updated_angles} angle parameter(s)")

    return output_file


def run_seminario_refinement_workflow(
    mol_name: str,
    mol2_file: str,
    frcmod_file: str,
    chk_file: str,
    selected_params: List[Tuple[str, float, str, str]],
    console: Console,
    interactive: bool = True,
    scale_factor: float = 1.0
) -> Dict[str, Any]:
    """
    Complete Seminario refinement workflow.

    This is the main entry point for Seminario-based parameter refinement.

    Args:
        mol_name: Molecule/residue name
        mol2_file: Path to MOL2 file with atom types
        frcmod_file: Path to original frcmod file
        chk_file: Path to Gaussian checkpoint file
        selected_params: List of (param_name, score, status, section) tuples
        console: Rich console for output
        interactive: Whether to show prompts
        scale_factor: Frequency scaling factor

    Returns:
        Dictionary with refinement results
    """
    results = {
        'refinement_success': False,
        'fitted_frcmod': None,
        'method': 'seminario',
        'message': ''
    }

    console.print(Panel(
        "[bold]Seminario Method Parameter Refinement[/bold]\n\n"
        "The Seminario method derives force constants directly from the Hessian matrix\n"
        "(second derivatives of energy) computed during the frequency calculation.\n\n"
        "[cyan]Advantages:[/cyan]\n"
        "• Uses existing opt+freq data (no additional QM needed)\n"
        "• Direct analytical result (no iterative fitting)\n"
        "• Well-established method for force constant extraction\n\n"
        "[yellow]Limitations:[/yellow]\n"
        "• Only refines bonds and angles (not dihedrals)\n"
        "• Based on single optimized structure",
        title="Seminario Refinement",
        expand=False
    ))

    # Check for dihedral parameters
    dihe_params = [p for p in selected_params if p[3] == 'DIHE']
    if dihe_params:
        console.print(f"\n[yellow]Note: {len(dihe_params)} dihedral parameter(s) cannot be refined with Seminario method[/yellow]")
        console.print("[grey50]Dihedrals require sampling multiple conformations (use paramfit instead)[/grey50]")

    # Convert checkpoint to formatted checkpoint
    fchk_file = convert_chk_to_fchk(chk_file, console)
    if not fchk_file:
        results['message'] = "Failed to convert checkpoint file"
        return results

    # Run Seminario refinement
    seminario_results = run_seminario_refinement(
        fchk_file=fchk_file,
        mol2_file=mol2_file,
        frcmod_file=frcmod_file,
        selected_params=selected_params,
        console=console,
        scale_factor=scale_factor
    )

    if seminario_results['success']:
        results['refinement_success'] = True
        results['fitted_frcmod'] = seminario_results['fitted_frcmod']
        results['bond_params'] = seminario_results.get('bond_params', [])
        results['angle_params'] = seminario_results.get('angle_params', [])
        results['message'] = seminario_results['message']

        # Display summary table
        display_seminario_summary(seminario_results, console)
    else:
        results['message'] = seminario_results['message']

    return results


def display_seminario_summary(results: Dict[str, Any], console: Console) -> None:
    """Display summary table of Seminario-derived parameters."""
    bond_params = results.get('bond_params', [])
    angle_params = results.get('angle_params', [])

    if not bond_params and not angle_params:
        return

    table = Table(title="Seminario-Derived Force Constants")
    table.add_column("Type", style="cyan")
    table.add_column("Atoms", style="white")
    table.add_column("Force Constant", style="green")
    table.add_column("Equilibrium", style="green")
    table.add_column("Std Dev", style="grey50")

    for param in bond_params:
        table.add_row(
            "BOND",
            f"{param.atom1_type}-{param.atom2_type}",
            f"{param.force_constant:.2f} kcal/(mol·Å²)",
            f"{param.eq_length:.4f} Å",
            f"±{param.std_dev:.2f}" if param.std_dev > 0 else "-"
        )

    for param in angle_params:
        table.add_row(
            "ANGLE",
            f"{param.atom1_type}-{param.atom2_type}-{param.atom3_type}",
            f"{param.force_constant:.2f} kcal/(mol·rad²)",
            f"{param.eq_angle:.2f}°",
            f"±{param.std_dev:.2f}" if param.std_dev > 0 else "-"
        )

    console.print(table)

    # Quality checks
    console.print("\n[bold]Parameter Quality Check:[/bold]")

    issues = []
    for param in bond_params:
        if param.force_constant < 50:
            issues.append(f"  [yellow]⚠ {param.atom1_type}-{param.atom2_type}: Low bond force constant ({param.force_constant:.1f})[/yellow]")
        elif param.force_constant > 1000:
            issues.append(f"  [yellow]⚠ {param.atom1_type}-{param.atom2_type}: High bond force constant ({param.force_constant:.1f})[/yellow]")

    for param in angle_params:
        if param.force_constant < 10:
            issues.append(f"  [yellow]⚠ {param.atom1_type}-{param.atom2_type}-{param.atom3_type}: Low angle force constant ({param.force_constant:.1f})[/yellow]")
        elif param.force_constant > 500:
            issues.append(f"  [yellow]⚠ {param.atom1_type}-{param.atom2_type}-{param.atom3_type}: High angle force constant ({param.force_constant:.1f})[/yellow]")

    if issues:
        for issue in issues:
            console.print(issue)
        console.print("[grey50]Note: Unusual values may indicate problematic QM geometry or basis set issues[/grey50]")
    else:
        console.print("[green]✓ All parameters are within typical ranges[/green]")
