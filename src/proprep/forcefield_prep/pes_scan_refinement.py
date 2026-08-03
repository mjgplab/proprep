"""
PES Scan Refinement Module for Dihedral Parameters

This module provides dihedral parameter refinement using relaxed potential
energy surface (PES) scans. Each dihedral is scanned systematically to
map the torsional potential, then paramfit is used to fit the dihedral
parameters (V, phase, periodicity).

Workflow:
1. For each selected dihedral:
   a. Generate Gaussian input with opt=modredundant and dihedral scan
   b. User runs Gaussian (or wait for completion)
   c. Parse scan energies and geometries from Gaussian log
2. Convert geometries to AMBER mdcrd trajectory format
3. Run paramfit to fit dihedral parameters (K first, then parameters)
4. Generate refined frcmod file

Author: ProPrep Development Team
"""

import os
import re
import subprocess
import shutil
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any

import numpy as np
from rich.console import Console
from rich.panel import Panel
from proprep.utils.prompts import prompt_with_context, confirm_with_context, int_prompt_with_context
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn

# Import working paramfit functions from the CREST module
from proprep.forcefield_prep.paramfit_refinement import (
    run_paramfit_k_fitting,
    run_paramfit_parameter_fitting,
    run_paramfit_set_params_automated,
    plot_energy_residuals,
    check_fitted_parameters,
    check_paramfit_availability,
)
from proprep.utils.prompts import confirm_with_context


# Conversion factors
HARTREE_TO_KCAL = 627.5095


def map_dihedral_to_atoms(
    dihedral_param: Tuple[str, float, str, str],
    mol2_file: str,
    console: Console
) -> Optional[Tuple[int, int, int, int]]:
    """
    Map dihedral atom types to atom indices.

    Args:
        dihedral_param: (param_name, score, status, section) tuple
        mol2_file: Path to MOL2 file
        console: Rich console for output

    Returns:
        Tuple of (atom1, atom2, atom3, atom4) indices (1-indexed for Gaussian),
        or None if not found
    """
    param_name = dihedral_param[0]
    types = param_name.replace(' ', '').split('-')

    if len(types) != 4:
        console.print(f"[yellow]Invalid dihedral format: {param_name}[/yellow]")
        return None

    # Parse MOL2 to get atom types and connectivity
    from proprep.forcefield_prep.seminario_refinement import parse_mol2_connectivity

    mol2_data = parse_mol2_connectivity(mol2_file)
    atom_types = mol2_data['atom_types']
    bonds = mol2_data['bonds']

    # Build adjacency for finding dihedral paths
    adjacency = {}
    for idx1, idx2 in bonds:
        if idx1 not in adjacency:
            adjacency[idx1] = []
        if idx2 not in adjacency:
            adjacency[idx2] = []
        adjacency[idx1].append(idx2)
        adjacency[idx2].append(idx1)

    # Search for matching dihedral
    t1, t2, t3, t4 = types

    # For each bond (potential central bond of dihedral)
    for idx2, idx3 in bonds:
        if not ((atom_types[idx2] == t2 and atom_types[idx3] == t3) or
                (atom_types[idx2] == t3 and atom_types[idx3] == t2)):
            continue

        # Swap if needed to match order
        if atom_types[idx2] == t3:
            idx2, idx3 = idx3, idx2

        # Find atom1 (bonded to idx2, type t1)
        for idx1 in adjacency.get(idx2, []):
            if idx1 == idx3:
                continue
            if atom_types[idx1] != t1:
                continue

            # Find atom4 (bonded to idx3, type t4)
            for idx4 in adjacency.get(idx3, []):
                if idx4 == idx2:
                    continue
                if atom_types[idx4] != t4:
                    continue

                # Found matching dihedral (convert to 1-indexed for Gaussian)
                console.print(f"[grey50]  Mapped {param_name} → atoms {idx1+1}-{idx2+1}-{idx3+1}-{idx4+1}[/grey50]")
                return (idx1 + 1, idx2 + 1, idx3 + 1, idx4 + 1)

    console.print(f"[yellow]Could not find atoms for dihedral: {param_name}[/yellow]")
    return None


def generate_pes_scan_input(
    mol_name: str,
    pdb_file: str,
    dihedral_atoms: Tuple[int, int, int, int],
    dihedral_name: str,
    charge: int,
    multiplicity: int,
    settings: Dict[str, Any],
    console: Console,
    n_steps: int = 24,
    step_size: float = 15.0,
    output_dir: str = None
) -> str:
    """
    Generate Gaussian input for relaxed PES scan.

    Args:
        mol_name: Molecule name
        pdb_file: Path to PDB file with coordinates
        dihedral_atoms: (atom1, atom2, atom3, atom4) 1-indexed
        dihedral_name: Name for output files
        charge: Net molecular charge
        multiplicity: Spin multiplicity
        settings: Gaussian settings dict
        console: Rich console
        n_steps: Number of scan steps (default 24 for 360° with 15° steps)
        step_size: Step size in degrees
        output_dir: Directory to write output files (default: current directory)

    Returns:
        Path to generated Gaussian input file
    """
    # Read coordinates from PDB
    coords = []
    elements = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                element = line[76:78].strip()
                if not element:
                    element = line[12:16].strip()[0]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append((x, y, z))
                elements.append(element)

    # Generate scan input
    base_name = f"{mol_name}_{dihedral_name}_scan"
    if output_dir:
        output_file = os.path.join(output_dir, f"{base_name}.gjf")
    else:
        output_file = f"{base_name}.gjf"
    a1, a2, a3, a4 = dihedral_atoms

    with open(output_file, 'w') as f:
        # Checkpoint file is relative (in same directory as input)
        f.write(f"%chk={base_name}.chk\n")
        f.write(f"%nprocshared={settings.get('processors', 4)}\n")
        f.write(f"%mem={settings.get('memory', '10GB')}\n")
        # Use opt=modredundant for relaxed scan
        keywords = settings.get('keywords', 'b3lyp/6-31G(d)')
        # Remove 'freq' and 'IOp' if present, add modredundant
        keywords = re.sub(r'freq\s*', '', keywords)
        keywords = re.sub(r'IOp\([^)]+\)\s*', '', keywords)
        keywords = keywords.replace('  ', ' ').strip()
        if 'opt' not in keywords.lower():
            keywords = f"opt=modredundant {keywords}"
        elif 'modredundant' not in keywords.lower():
            keywords = keywords.replace('opt', 'opt=modredundant')
        f.write(f"#p {keywords}\n")
        f.write("\n")
        f.write(f"{mol_name} dihedral scan: {dihedral_name}\n")
        f.write("\n")
        f.write(f"{charge} {multiplicity}\n")

        for elem, (x, y, z) in zip(elements, coords):
            f.write(f" {elem:2s}    {x:14.8f}  {y:14.8f}  {z:14.8f}\n")

        f.write("\n")
        # Add dihedral scan specification
        # Format: D atom1 atom2 atom3 atom4 S n_steps step_size
        f.write(f"D {a1} {a2} {a3} {a4} S {n_steps} {step_size:.1f}\n")
        f.write("\n")

    console.print(f"[green]✓ Generated PES scan input: {output_file}[/green]")
    console.print(f"[grey50]  Scanning dihedral {a1}-{a2}-{a3}-{a4} in {n_steps} steps of {step_size}°[/grey50]")

    return output_file


def parse_pes_scan_log(log_file: str, console: Console) -> Dict[str, Any]:
    """
    Parse PES scan results from Gaussian log file.

    Gaussian PES scans with opt=modredundant output:
    1. Multiple optimization cycles, each ending with "Stationary point found"
    2. A "Summary of Optimized Potential Surface Scan" section at the end
       containing the scan coordinate values and eigenvalues (energies)
    3. Archive entries with "HF=..." containing final energies

    This parser uses multiple strategies:
    - Primary: Parse the summary section for energies and scan coordinates
    - Fallback: Extract geometries at each optimization step
    - Additional: Parse archive entries for energy validation

    Args:
        log_file: Path to Gaussian log file
        console: Rich console

    Returns:
        Dictionary with:
        - 'energies': List of energies in Hartree
        - 'geometries': List of coordinate arrays
        - 'angles': List of dihedral angles
        - 'elements': List of element symbols
        - 'n_atoms': Number of atoms
        - 'success': True if parsing succeeded
    """
    result = {
        'energies': [],
        'geometries': [],
        'angles': [],
        'elements': [],
        'n_atoms': 0,
        'success': False
    }

    if not os.path.exists(log_file):
        console.print(f"[red]Log file not found: {log_file}[/red]")
        return result

    with open(log_file, 'r') as f:
        content = f.read()

    # Check for normal termination
    if "Normal termination" not in content:
        console.print(f"[yellow]Warning: Gaussian job may not have completed normally[/yellow]")

    lines = content.split('\n')

    # Parse number of atoms
    for line in lines:
        if "NAtoms=" in line:
            match = re.search(r'NAtoms=\s*(\d+)', line)
            if match:
                result['n_atoms'] = int(match.group(1))
                break

    if result['n_atoms'] == 0:
        console.print(f"[red]Could not determine number of atoms[/red]")
        return result

    # Atomic number to symbol mapping (elements 1-86)
    atomic_symbols = {
        1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O',
        9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',
        16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti',
        23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu',
        30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr',
        37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc',
        44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
        51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La',
        58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd',
        65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu',
        72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt',
        79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At',
        86: 'Rn'
    }

    # Strategy 1: Parse "Summary of Optimized Potential Surface Scan" section
    # This section appears at the end and contains the definitive scan results
    #
    # Gaussian 16+ format with base energy offset:
    #     Summary of Optimized Potential Surface Scan (add -382.0 to eigenvalues):
    #   D(1,2,3,4)               Eigenvalues --   -.012345
    #     0.0000                      -.012345
    #    15.0000                      -.011234
    #
    # Older format without offset:
    #     Summary of Optimized Potential Surface Scan
    #   D(1,2,3,4)               Eigenvalues -- -382.012345
    #     0.0000                  -382.012345
    #    15.0000                  -382.011234
    summary_energies = []
    summary_angles = []
    base_energy_offset = 0.0

    # Look for the summary section with optional base energy offset (Gaussian 16+)
    summary_match = re.search(
        r'Summary of Optimized Potential Surface Scan'
        r'(?:\s*\(add\s*(-?\d+\.?\d*)\s*to\s*eigenvalues\))?'
        r'[^\n]*\n(.*?)(?:GradGrad|Normal termination|$)',
        content, re.DOTALL | re.IGNORECASE
    )

    if summary_match:
        # Extract base energy offset if present (Gaussian 16+ format)
        if summary_match.group(1):
            base_energy_offset = float(summary_match.group(1))
            console.print(f"[grey50]  Found energy offset: {base_energy_offset}[/grey50]")

        summary_text = summary_match.group(2)

        # Parse the summary table
        # cclib uses fixed-width column parsing (10 chars each for energies)
        # Each line has: param_value(s)   energy_value(s)
        for line in summary_text.split('\n'):
            # Skip header lines
            if 'Eigenvalues' in line or 'D(' in line or 'A(' in line or not line.strip():
                continue

            # Extract values from line - can have multiple columns
            parts = line.split()
            if len(parts) >= 2:
                try:
                    # First value is scan coordinate (angle)
                    angle = float(parts[0])
                    # Last value is energy (may have multiple columns)
                    energy = float(parts[-1])
                    # Apply base energy offset
                    energy = energy + base_energy_offset
                    summary_angles.append(angle)
                    summary_energies.append(energy)
                except ValueError:
                    continue
            elif len(parts) == 1:
                # Some formats put angle and energy on separate lines
                try:
                    val = float(parts[0])
                    # Apply offset and add to appropriate list
                    if abs(val) < 360:
                        summary_angles.append(val)
                    else:
                        summary_energies.append(val + base_energy_offset)
                except ValueError:
                    continue

    # Also extract scan length if available
    scan_length_match = re.search(r'Number of optimizations in scan\s*=\s*(\d+)', content)
    expected_scan_length = int(scan_length_match.group(1)) if scan_length_match else None

    # Strategy 2: Parse archive entries for HF= values
    # Format: \\HF=-382.012345\\
    archive_energies = []
    archive_matches = re.findall(r'\\HF=(-?\d+\.\d+)\\', content)
    archive_energies = [float(e) for e in archive_matches]

    # Strategy 3: Parse step-by-step for geometries and energies
    step_geometries = []
    step_energies = []
    current_geometry = []
    current_energy = None
    in_geometry = False
    geometry_line_count = 0
    skip_lines = 0

    for i, line in enumerate(lines):
        # Look for SCF energy (take the last one before stationary point)
        if "SCF Done:" in line:
            match = re.search(r'E\([^)]+\)\s*=\s*(-?\d+\.\d+)', line)
            if match:
                current_energy = float(match.group(1))

        # Look for optimized geometry. Match BOTH orientations: a scan written
        # with NoSymm (as the from-structure route does, to keep the atom order
        # stable for the ESP job) prints "Input orientation:" and never
        # "Standard orientation:", so keying only on the latter finds zero
        # geometries and the whole scan looks empty downstream.
        if "Standard orientation:" in line or "Input orientation:" in line:
            in_geometry = True
            current_geometry = []
            geometry_line_count = 0
            skip_lines = 4  # Skip header lines
            continue

        if in_geometry:
            if skip_lines > 0:
                skip_lines -= 1
                continue

            # End of geometry section (dashed line after coordinates)
            if re.match(r'\s*-{50,}', line):
                in_geometry = False
                continue

            # Parse coordinate line
            # Format: Center  Atomic   Atomic     Coordinates (Angstroms)
            #         Number  Number   Type       X           Y           Z
            parts = line.split()
            if len(parts) >= 6:
                try:
                    atom_num = int(parts[1])
                    x = float(parts[3])
                    y = float(parts[4])
                    z = float(parts[5])
                    current_geometry.append((atom_num, x, y, z))
                except (ValueError, IndexError):
                    pass

        # Look for scan step completion
        # "-- Stationary point found" indicates successful optimization at a scan point
        if "-- Stationary point found" in line or "Optimization completed" in line:
            if current_energy is not None and len(current_geometry) == result['n_atoms']:
                step_energies.append(current_energy)
                step_geometries.append(current_geometry.copy())
            current_geometry = []

    # Determine which energy source to use.
    #
    # Prefer step-by-step: those energies are paired 1:1 with the geometries we
    # extracted (each is the SCF energy at that scan point's stationary point),
    # which is exactly what the ESP job at every point needs. It is also immune
    # to the column-batched summary layout that the summary parser above
    # mis-reads — Gaussian prints the summary as blocks of five points with an
    # index header row ("1 2 3 4 5"), and that header parses as a spurious
    # (angle=1, energy=5) pair, yielding one garbage point per block. Fall back
    # to archive HF= entries, then the summary, only when step-by-step is
    # unavailable (and geometries are then absent, so downstream can't use it
    # anyway — but the energies are still reported for diagnostics).
    if step_energies and step_geometries and len(step_energies) == len(step_geometries):
        result['energies'] = step_energies
        result['angles'] = []  # regenerated from the point count below
        console.print(f"[grey50]  Using energies from step-by-step parsing "
                      f"(paired with geometries)[/grey50]")
    elif archive_energies and len(archive_energies) > 1:
        result['energies'] = archive_energies
        console.print(f"[grey50]  Using energies from archive entries[/grey50]")
    elif summary_energies and len(summary_energies) > 1:
        result['energies'] = summary_energies
        result['angles'] = summary_angles
        console.print(f"[grey50]  Using energies from scan summary section[/grey50]")
    elif step_energies:
        result['energies'] = step_energies
        console.print(f"[grey50]  Using energies from step-by-step parsing[/grey50]")

    # Use geometries from step-by-step parsing (only source)
    result['geometries'] = step_geometries

    # Extract elements from first geometry
    if result['geometries']:
        first_geom = result['geometries'][0]
        result['elements'] = [atomic_symbols.get(atom[0], 'X') for atom in first_geom]

    # Generate angles if not found in summary
    if not result['angles'] and result['energies']:
        n_points = len(result['energies'])
        # Assume 15° steps starting from 0
        result['angles'] = [i * 15.0 for i in range(n_points)]
        console.print(f"[grey50]  Generated angles (15° steps assumed)[/grey50]")

    # Validate and report
    if result['energies']:
        n_energies = len(result['energies'])
        n_geometries = len(result['geometries'])

        # Check against expected scan length if available
        if expected_scan_length and n_energies != expected_scan_length:
            console.print(f"[yellow]Warning: Expected {expected_scan_length} scan points, found {n_energies}[/yellow]")

        # Warn if mismatch between energies and geometries
        if n_geometries > 0 and n_energies != n_geometries:
            console.print(f"[yellow]Warning: {n_energies} energies but {n_geometries} geometries[/yellow]")
            # Try to align - take the shorter set
            if n_geometries < n_energies:
                result['energies'] = result['energies'][:n_geometries]
                result['angles'] = result['angles'][:n_geometries] if result['angles'] else []
            else:
                result['geometries'] = result['geometries'][:n_energies]

        result['success'] = True
        console.print(f"[green]✓ Parsed {len(result['energies'])} scan points[/green]")

        # Show energy range
        e_min = min(result['energies'])
        e_max = max(result['energies'])
        e_range_kcal = (e_max - e_min) * HARTREE_TO_KCAL
        console.print(f"[grey50]  Energy range: {e_range_kcal:.2f} kcal/mol[/grey50]")

        # Show angle range if available
        if result['angles']:
            a_min, a_max = min(result['angles']), max(result['angles'])
            console.print(f"[grey50]  Scan range: {a_min:.1f}° to {a_max:.1f}°[/grey50]")
    else:
        console.print(f"[yellow]No scan data found in log file[/yellow]")

    return result


def write_scan_mdcrd(
    scan_data: Dict[str, Any],
    mol_name: str,
    output_dir: str,
    console: Console
) -> str:
    """
    Write AMBER mdcrd trajectory file containing all scan conformations.

    Paramfit requires a single mdcrd file with all conformations, not individual
    coordinate files. The mdcrd format is:
    - Title line
    - For each frame: coordinates in 10F8.3 format (10 values per line, 8.3f)

    Args:
        scan_data: Parsed scan data from parse_pes_scan_log
        mol_name: Molecule name
        output_dir: Directory for output files
        console: Rich console

    Returns:
        Path to mdcrd trajectory file
    """
    os.makedirs(output_dir, exist_ok=True)
    mdcrd_file = os.path.join(output_dir, f"{mol_name}_scan.mdcrd")

    n_atoms = scan_data['n_atoms']
    n_frames = len(scan_data['geometries'])

    with open(mdcrd_file, 'w') as f:
        # Title line
        f.write(f"{mol_name} PES scan trajectory with {n_frames} conformations\n")

        # Write each frame
        for frame_idx, geometry in enumerate(scan_data['geometries']):
            # Collect all coordinates for this frame
            coords = []
            for atom in geometry:
                coords.extend([atom[1], atom[2], atom[3]])

            # Write coordinates in mdcrd format: 10F8.3 (10 values per line)
            for j in range(0, len(coords), 10):
                line_coords = coords[j:j+10]
                f.write(''.join(f"{c:8.3f}" for c in line_coords) + '\n')

    console.print(f"[green]✓ Wrote mdcrd trajectory: {mdcrd_file} ({n_frames} frames)[/green]")
    return mdcrd_file


def write_scan_energy_file(
    scan_data: Dict[str, Any],
    mol_name: str,
    output_dir: str,
    console: Console
) -> str:
    """
    Write QM energies file for paramfit.

    Paramfit expects energies in Hartree (one per line) when using
    QM_ENERGY_UNITS=HARTREE in the job control file.

    Args:
        scan_data: Parsed scan data
        mol_name: Molecule name
        output_dir: Directory for output
        console: Rich console

    Returns:
        Path to energy file
    """
    os.makedirs(output_dir, exist_ok=True)
    energy_file = os.path.join(output_dir, f"{mol_name}_scan_qm.dat")

    # Write energies in Hartree (paramfit expects this format)
    energies_hartree = scan_data['energies']

    with open(energy_file, 'w') as f:
        for e in energies_hartree:
            f.write(f"{e}\n")

    console.print(f"[green]✓ Wrote QM energies: {energy_file} ({len(energies_hartree)} values in Hartree)[/green]")

    # Show energy range in kcal/mol for user info
    e_min = min(energies_hartree)
    e_max = max(energies_hartree)
    e_range_kcal = (e_max - e_min) * HARTREE_TO_KCAL
    console.print(f"[grey50]  Energy range: {e_range_kcal:.2f} kcal/mol[/grey50]")

    return energy_file


def merge_dihedral_into_frcmod(
    original_frcmod: str,
    fitted_frcmod: str,
    dihedral_types: str,
    output_file: str,
    console: Console
) -> str:
    """
    Merge fitted dihedral parameters into original frcmod.

    Args:
        original_frcmod: Path to original frcmod
        fitted_frcmod: Path to fitted frcmod with new dihedral
        dihedral_types: Dihedral type string
        output_file: Output frcmod path
        console: Rich console

    Returns:
        Path to merged frcmod
    """
    # Read original frcmod
    with open(original_frcmod, 'r') as f:
        original_lines = f.readlines()

    # Read fitted frcmod to extract dihedral parameters
    fitted_dihedrals = []
    if os.path.exists(fitted_frcmod):
        with open(fitted_frcmod, 'r') as f:
            in_dihe = False
            for line in f:
                if line.strip() == 'DIHE':
                    in_dihe = True
                    continue
                elif line.strip() in ['IMPROPER', 'NONBON', 'END', '']:
                    in_dihe = False
                elif in_dihe and line.strip():
                    fitted_dihedrals.append(line)

    if not fitted_dihedrals:
        console.print(f"[yellow]No fitted dihedrals found in {fitted_frcmod}[/yellow]")
        shutil.copy(original_frcmod, output_file)
        return output_file

    # Process original frcmod and replace matching dihedral
    output_lines = []
    in_dihe = False
    normalized_target = dihedral_types.replace(' ', '').replace('-', '')

    for line in original_lines:
        if line.strip() == 'DIHE':
            in_dihe = True
            output_lines.append(line)
            continue
        elif line.strip() in ['IMPROPER', 'NONBON', 'END']:
            if in_dihe:
                # Add fitted dihedrals before leaving section
                for fitted_line in fitted_dihedrals:
                    output_lines.append(fitted_line)
            in_dihe = False
            output_lines.append(line)
            continue

        if in_dihe and line.strip():
            # Check if this is the dihedral we're replacing
            # AMBER format: first 11 chars are atom types
            if len(line) >= 11:
                line_types = line[:11].replace(' ', '').replace('-', '')
                if line_types == normalized_target or line_types == normalized_target[::-1]:
                    # Skip original line, we'll add fitted version
                    continue

        output_lines.append(line)

    # Write output
    with open(output_file, 'w') as f:
        f.writelines(output_lines)

    console.print(f"[green]✓ Merged dihedral into: {output_file}[/green]")
    return output_file


def run_pes_scan_workflow(
    mol_name: str,
    mol2_file: str,
    prmtop_file: str,
    frcmod_file: str,
    selected_dihedrals: List[Tuple[str, float, str, str]],
    charge: int,
    multiplicity: int,
    console: Console,
    interactive: bool = True,
    processor=None,
    gaussian_settings: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    Run PES scan workflow for dihedral refinement.

    This workflow processes each selected dihedral:
    1. Check for existing Gaussian PES scan log file
    2. If no log: generate Gaussian input with opt=modredundant
    3. If log exists: parse scan, fit K, select params, fit dihedral, merge to frcmod

    The user can run multiple dihedrals - each gets its own scan and fit,
    with results merged into the final frcmod file.

    Args:
        mol_name: Molecule name
        mol2_file: Path to MOL2 file
        prmtop_file: Path to prmtop file
        frcmod_file: Path to frcmod file
        selected_dihedrals: List of (param_name, score, status, section) tuples
        charge: Net molecular charge
        multiplicity: Spin multiplicity
        console: Rich console
        interactive: Whether to prompt user
        processor: Optional processor for session recording
        gaussian_settings: Optional settings from earlier Gaussian step (carries SCRF, memory, etc.)

    Returns:
        Dictionary with refinement results
    """
    results = {
        'refinement_success': False,
        'fitted_frcmod': None,
        'method': 'pes_scan',
        'message': '',
        'scan_inputs': [],
        'completed_scans': []
    }

    console.print(Panel(
        "[bold]PES Scan Dihedral Refinement[/bold]\n\n"
        "This workflow:\n"
        "1. Generates Gaussian input for relaxed PES scan of each dihedral\n"
        "2. Parses scan energies and geometries from Gaussian output\n"
        "3. Fits K (energy offset) first, then dihedral parameters\n"
        "4. Uses paramfit genetic algorithm to optimize V, phase, periodicity\n"
        "5. Merges fitted parameters into final frcmod file\n\n"
        "[grey50]Each scan systematically rotates the dihedral in 15° steps (24 points)[/grey50]",
        title="PES Scan Workflow",
        expand=False
    ))

    # Check paramfit availability
    paramfit_ok, paramfit_msg = check_paramfit_availability()
    if not paramfit_ok:
        results['message'] = f"paramfit not available: {paramfit_msg}"
        console.print(f"[red]{results['message']}[/red]")
        return results

    # Find PDB file for coordinates - prefer hydrogenated version
    pdb_file = f"{mol_name}_H.pdb"
    if not os.path.exists(pdb_file):
        pdb_file = f"{mol_name}.pdb"
        if not os.path.exists(pdb_file):
            results['message'] = f"PDB file not found: {mol_name}.pdb or {mol_name}_H.pdb"
            console.print(f"[red]{results['message']}[/red]")
            return results
        console.print(f"[yellow]Note: Using {pdb_file} (no hydrogenated version found)[/yellow]")

    # Create output directory for scan files
    scan_dir = f"{mol_name}_pes_scans"
    os.makedirs(scan_dir, exist_ok=True)
    original_dir = os.getcwd()

    # Extract SCRF keyword from upstream Gaussian settings if present
    upstream_scrf = ""
    default_memory = "10GB"
    default_nproc = 4
    if gaussian_settings and gaussian_settings.get("opt", {}).get("keywords"):
        opt_kw = gaussian_settings["opt"]["keywords"]
        scrf_match = re.search(r'SCRF=\([^)]+\)', opt_kw, re.IGNORECASE)
        if scrf_match:
            upstream_scrf = scrf_match.group(0)
        default_memory = gaussian_settings["opt"].get("memory", "10GB")
        default_nproc = gaussian_settings["opt"].get("processors", 4)

    # Get Gaussian settings from user if interactive
    if interactive:
        console.print("\n[cyan]Gaussian PES scan settings:[/cyan]")

        default_method = "B3LYP/6-31G(d)"
        if upstream_scrf:
            default_method = f"B3LYP/6-31G(d) {upstream_scrf}"
            console.print(f"[grey50]Carrying forward implicit solvation from optimization step[/grey50]")

        keywords = prompt_with_context(
            processor,
            "DFT method/basis for PES scan",
            default=default_method,
            module="PES Scan Refinement",
            description="DFT method/basis for PES scan",
        )

        memory = prompt_with_context(
            processor,
            "Memory allocation",
            default=default_memory,
            module="PES Scan Refinement",
            description="Memory allocation",
        )

        nproc = int_prompt_with_context(
            processor,
            "Number of processors",
            default=default_nproc,
            module="PES Scan Refinement",
            description="Number of processors",
        )

        n_steps = int_prompt_with_context(
            processor,
            "Number of scan steps (24 = 15° increments for 360°)",
            default=24,
            module="PES Scan Refinement",
            description="Number of PES scan steps",
        )

        step_size = 360.0 / n_steps
    else:
        keywords = "B3LYP/6-31G(d)"
        if upstream_scrf:
            keywords = f"{keywords} {upstream_scrf}"
        memory = default_memory
        nproc = default_nproc
        n_steps = 24
        step_size = 15.0

    settings = {
        'memory': memory,
        'processors': nproc,
        'keywords': f'opt=modredundant {keywords} nosym'
    }

    # Process each dihedral
    current_frcmod = os.path.abspath(frcmod_file)
    prmtop_path = os.path.abspath(prmtop_file)

    for dihe_idx, dihe_param in enumerate(selected_dihedrals):
        param_name = dihe_param[0]
        safe_name = param_name.replace(' ', '').replace('-', '_')

        console.print(f"\n[bold cyan]═══ Dihedral {dihe_idx + 1}/{len(selected_dihedrals)}: {param_name} ═══[/bold cyan]")

        # Map to atom indices
        dihedral_atoms = map_dihedral_to_atoms(dihe_param, mol2_file, console)
        if not dihedral_atoms:
            console.print(f"[yellow]Could not map dihedral atoms, skipping[/yellow]")
            continue

        # Check for existing log file (in scan directory)
        base_name = f"{mol_name}_{safe_name}_scan"
        log_file = os.path.join(scan_dir, f"{base_name}.log")

        if os.path.exists(log_file):
            console.print(f"[green]Found existing scan log: {log_file}[/green]")

            # Parse scan data
            scan_data = parse_pes_scan_log(log_file, console)

            if not scan_data['success']:
                console.print(f"[yellow]Could not parse scan data, skipping this dihedral[/yellow]")
                continue

            if len(scan_data['geometries']) < 5:
                console.print(f"[yellow]Only {len(scan_data['geometries'])} scan points - need at least 5 for fitting[/yellow]")
                continue

            # Create working directory for this dihedral
            dihe_dir = os.path.join(scan_dir, safe_name)
            os.makedirs(dihe_dir, exist_ok=True)
            os.chdir(dihe_dir)

            try:
                # Write mdcrd trajectory file (all scan conformations)
                mdcrd_file = write_scan_mdcrd(scan_data, mol_name, ".", console)

                # Write energies file (Hartree, one per line)
                energy_file = write_scan_energy_file(scan_data, mol_name, ".", console)

                n_structures = len(scan_data['geometries'])

                # Step 1: Fit K (energy offset constant)
                console.print(f"\n[bold]Fitting K (energy offset)...[/bold]")
                k_value = run_paramfit_k_fitting(
                    prmtop_file=prmtop_path,
                    mdcrd_file=mdcrd_file,
                    energies_file=energy_file,
                    n_structures=n_structures,
                    console=console
                )

                if k_value is None:
                    console.print(f"[yellow]Failed to fit K, skipping this dihedral[/yellow]")
                    os.chdir(original_dir)
                    continue

                # Step 2: Select dihedral parameters using automated SET_PARAMS
                console.print(f"\n[bold]Selecting dihedral parameter for fitting...[/bold]")
                params_file = f"{mol_name}_params.in"

                # Create a param tuple list for just this dihedral
                dihe_param_list = [dihe_param]

                success = run_paramfit_set_params_automated(
                    prmtop_file=prmtop_path,
                    params_file=params_file,
                    selected_params=dihe_param_list,
                    console=console,
                    force_constants_only=False  # Fit KP, NP, and PHASE for dihedrals
                )

                if not success:
                    console.print(f"[yellow]Failed to select parameters, skipping this dihedral[/yellow]")
                    os.chdir(original_dir)
                    continue

                # Step 3: Run paramfit to fit the dihedral
                console.print(f"\n[bold]Fitting dihedral parameters...[/bold]")
                output_frcmod = f"{mol_name}_{safe_name}_fitted.frcmod"

                success = run_paramfit_parameter_fitting(
                    prmtop_file=prmtop_path,
                    mdcrd_file=mdcrd_file,
                    energies_file=energy_file,
                    params_file=params_file,
                    n_structures=n_structures,
                    k_value=k_value,
                    output_frcmod=output_frcmod,
                    console=console
                )

                if success and os.path.exists(output_frcmod):
                    # Check fitted parameters for issues
                    check_fitted_parameters(output_frcmod, console)

                    # Merge into cumulative frcmod
                    merged_frcmod = os.path.join(original_dir, scan_dir, f"{mol_name}_merged_{dihe_idx}.frcmod")
                    current_frcmod = merge_dihedral_into_frcmod(
                        current_frcmod,
                        os.path.abspath(output_frcmod),
                        param_name,
                        merged_frcmod,
                        console
                    )
                    results['completed_scans'].append(param_name)
                    console.print(f"[green]✓ Dihedral {param_name} fitted successfully[/green]")
                else:
                    console.print(f"[yellow]Parameter fitting failed for {param_name}[/yellow]")

            finally:
                os.chdir(original_dir)

        else:
            # Generate scan input file in scan directory
            console.print(f"[cyan]No scan log found. Generating Gaussian input...[/cyan]")

            input_file = generate_pes_scan_input(
                mol_name, pdb_file, dihedral_atoms, safe_name,
                charge, multiplicity, settings, console,
                n_steps=n_steps, step_size=step_size,
                output_dir=scan_dir
            )
            results['scan_inputs'].append({
                'param_name': param_name,
                'input_file': input_file,
                'log_file': log_file
            })

    # Final Summary
    console.print(f"\n[bold]═══ PES Scan Workflow Summary ═══[/bold]")

    if results['completed_scans']:
        results['refinement_success'] = True
        results['fitted_frcmod'] = current_frcmod
        results['message'] = f"Refined {len(results['completed_scans'])} dihedral(s)"

        console.print(Panel(
            f"[bold green]Successfully refined {len(results['completed_scans'])} dihedral(s):[/bold green]\n\n" +
            '\n'.join(f"  [green]✓[/green] {name}" for name in results['completed_scans']) +
            f"\n\n[bold]Output frcmod:[/bold]\n  {current_frcmod}",
            title="Refinement Complete",
            border_style="green",
            expand=False
        ))

    if results['scan_inputs']:
        console.print(Panel(
            f"[bold yellow]Gaussian calculations needed for {len(results['scan_inputs'])} dihedral(s):[/bold yellow]\n\n" +
            '\n'.join(f"  • {s['param_name']}: {s['input_file']}" for s in results['scan_inputs']) +
            "\n\n[bold]To run the scans:[/bold]\n" +
            '\n'.join(f"  g16 {s['input_file']} > {s['log_file']}" for s in results['scan_inputs']) +
            "\n\n[grey50]After Gaussian completes, re-run this step to fit parameters.[/grey50]",
            title="Action Required",
            border_style="yellow",
            expand=False
        ))

        # Ask user if they've completed the Gaussian scans
        if interactive and not results['completed_scans']:
            console.print(f"\n[grey50]To resume later: Run Gaussian, then re-run this workflow.[/grey50]")
            proceed = confirm_with_context(
                processor,
                "Have you completed the PES scan calculations?",
                default=False,
                module="PES Scan Refinement",
                description="Confirm PES scan calculations completed"
            )

            if not proceed:
                results['message'] = "PES scan calculations required"
                return results

            # Re-process dihedrals now that logs should exist
            console.print(f"\n[cyan]Re-checking for completed scan logs...[/cyan]")
            for scan_info in results['scan_inputs']:
                log_file = scan_info['log_file']
                param_name = scan_info['param_name']

                if not os.path.exists(log_file):
                    console.print(f"[yellow]Log file not found: {log_file}[/yellow]")
                    continue

                console.print(f"[green]Found scan log: {log_file}[/green]")
                scan_data = parse_pes_scan_log(log_file, console)

                if not scan_data['success']:
                    console.print(f"[yellow]Could not parse scan data for {param_name}[/yellow]")
                    continue

                if len(scan_data['geometries']) < 5:
                    console.print(f"[yellow]Only {len(scan_data['geometries'])} scan points for {param_name} - need at least 5[/yellow]")
                    continue

                safe_name = param_name.replace(' ', '').replace('-', '_')
                dihe_idx = next(
                    (i for i, s in enumerate(results['scan_inputs']) if s['param_name'] == param_name), 0
                )

                dihe_dir = os.path.join(scan_dir, safe_name)
                os.makedirs(dihe_dir, exist_ok=True)
                os.chdir(dihe_dir)

                try:
                    mdcrd_file = write_scan_mdcrd(scan_data, mol_name, ".", console)
                    energy_file = write_scan_energy_file(scan_data, mol_name, ".", console)
                    n_structures = len(scan_data['geometries'])

                    console.print(f"\n[bold]Fitting K (energy offset)...[/bold]")
                    k_value = run_paramfit_k_fitting(
                        prmtop_file=prmtop_path,
                        mdcrd_file=mdcrd_file,
                        energies_file=energy_file,
                        n_structures=n_structures,
                        console=console
                    )

                    if k_value is None:
                        console.print(f"[yellow]Failed to fit K for {param_name}[/yellow]")
                        continue

                    console.print(f"\n[bold]Selecting dihedral parameter for fitting...[/bold]")
                    params_file = f"{mol_name}_params.in"

                    # Find the original dihedral param tuple
                    dihe_param = next(
                        (d for d in selected_dihedrals if d[0] == param_name), None
                    )
                    if not dihe_param:
                        console.print(f"[yellow]Could not find param info for {param_name}[/yellow]")
                        continue

                    success = run_paramfit_set_params_automated(
                        prmtop_file=prmtop_path,
                        params_file=params_file,
                        selected_params=[dihe_param],
                        console=console,
                        force_constants_only=False
                    )

                    if not success:
                        console.print(f"[yellow]Failed to select parameters for {param_name}[/yellow]")
                        continue

                    console.print(f"\n[bold]Running paramfit for {param_name}...[/bold]")
                    output_frcmod = f"{mol_name}_{safe_name}_fitted.frcmod"

                    fit_success = run_paramfit_parameter_fitting(
                        prmtop_file=prmtop_path,
                        mdcrd_file=mdcrd_file,
                        energies_file=energy_file,
                        params_file=params_file,
                        n_structures=n_structures,
                        k_value=k_value,
                        output_frcmod=output_frcmod,
                        console=console
                    )

                    if fit_success and os.path.exists(output_frcmod):
                        check_fitted_parameters(output_frcmod, console)

                        merged_frcmod = os.path.join(original_dir, scan_dir, f"{mol_name}_merged_{dihe_idx}.frcmod")
                        current_frcmod = merge_dihedral_into_frcmod(
                            current_frcmod,
                            os.path.abspath(output_frcmod),
                            param_name,
                            merged_frcmod,
                            console
                        )
                        results['completed_scans'].append(param_name)
                        console.print(f"[green]✓ Dihedral {param_name} fitted successfully[/green]")
                    else:
                        console.print(f"[yellow]Parameter fitting failed for {param_name}[/yellow]")

                finally:
                    os.chdir(original_dir)

            # Remove successfully fitted dihedrals from scan_inputs
            results['scan_inputs'] = [
                s for s in results['scan_inputs']
                if s['param_name'] not in results['completed_scans']
            ]

            # Update summary if we completed scans
            if results['completed_scans']:
                results['refinement_success'] = True
                results['fitted_frcmod'] = current_frcmod
                results['message'] = f"Refined {len(results['completed_scans'])} dihedral(s)"

                console.print(Panel(
                    f"[bold green]Successfully refined {len(results['completed_scans'])} dihedral(s):[/bold green]\n\n" +
                    '\n'.join(f"  [green]✓[/green] {name}" for name in results['completed_scans']) +
                    f"\n\n[bold]Output frcmod:[/bold]\n  {current_frcmod}",
                    title="Refinement Complete",
                    border_style="green",
                    expand=False
                ))

        if not results['completed_scans'] and results['scan_inputs']:
            results['message'] = f"Generated {len(results['scan_inputs'])} scan input(s); awaiting Gaussian"

    if not results['completed_scans'] and not results['scan_inputs']:
        results['message'] = "No dihedrals could be processed"
        console.print(f"[yellow]{results['message']}[/yellow]")

    return results
