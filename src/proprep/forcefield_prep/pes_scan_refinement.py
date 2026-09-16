"""
Relaxed PES scan parsing for dihedral refinement.

What remains here is the Gaussian side that every route shares:
``map_dihedral_to_atoms`` (a type quad → one atom quadruple in a MOL2) and
``parse_pes_scan_log`` (a relaxed ``opt=modredundant`` scan log → geometries
and energies). The per-dihedral fitting workflow that used to live here was
replaced by the joint fit in :mod:`proprep.forcefield_prep.torsion_refinement`,
which the small-molecule step sm-7 and the modified-amino-acid step 9 both
call.
"""

import os
import re
from typing import Dict, List, Tuple, Optional, Any

from rich.console import Console

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


def parse_pes_scan_log(log_file: str, console: Console, step_size: float = 15.0) -> Dict[str, Any]:
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
        # Regenerate from the scan's own step size (the caller reads it from the
        # input's ``D ... S n step`` line); 15° is only the historical default.
        result['angles'] = [i * float(step_size) for i in range(n_points)]
        console.print(f"[grey50]  Generated angles ({step_size:g}° steps)[/grey50]")

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
