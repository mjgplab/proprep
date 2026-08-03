"""
Paramfit Refinement Module for Small Molecule Parameters

This module provides QM-based refinement of force field parameters using:
1. CREST (Conformer-Rotamer Ensemble Sampling Tool) for conformer generation
2. Gaussian for single-point QM energies
3. AmberTools paramfit for parameter fitting

The workflow addresses high-penalty parameters identified by parmchk2 by fitting
them to QM energies across multiple molecular conformations.

Scientific rationale:
- parmchk2 estimates parameters by analogy to similar atom types
- High penalties indicate unreliable estimates
- Fitting to QM energies across conformations improves accuracy
- CREST uses xTB (semiempirical QM) for unbiased conformer generation

Author: ProPrep Development Team
"""

import os
import subprocess
import shutil
import pty
import re
import select
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any, Set

from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn

from proprep.utils.prompts import (
    prompt_with_context,
    confirm_with_context,
    int_prompt_with_context,
)


def _add_session_context(processor, module: str, description: str):
    """Add context to the last session recorder interaction.

    Args:
        processor: The processor object with session_manager
        module: Module name for context (e.g., "Paramfit Refinement")
        description: Description of what the prompt is asking
    """
    if processor and hasattr(processor, 'session_manager'):
        recorder = processor.session_manager.recorder
        if recorder and recorder.recording and recorder.session_data["interactions"]:
            last_interaction = recorder.session_data["interactions"][-1]
            last_interaction["context"]["module"] = module
            last_interaction["context"]["description"] = description


def normalize_param_atoms(param_string: str) -> str:
    """
    Normalize parameter atom type string for comparison.

    Removes all whitespace so that:
    - "nc-c -ns" (paramfit format) → "nc-c-ns"
    - "nc-c-ns" (ProPrep stored format) → "nc-c-ns"

    This allows matching between the two formats.

    Args:
        param_string: Atom type string like "nc-c -ns" or "nc-c-ns"

    Returns:
        Normalized string with whitespace removed
    """
    return param_string.replace(" ", "")


def plot_energy_residuals(energies_file: str, console: Console) -> None:
    """
    Create an ASCII plot of energy residuals (Amber+K - Quantum).

    Args:
        energies_file: Path to the paramfit energies output file
        console: Rich console for output
    """
    if not os.path.exists(energies_file):
        return

    try:
        # Parse energies file
        residuals = []
        with open(energies_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        idx = int(parts[0])
                        amber_k = float(parts[1])
                        quantum = float(parts[2])
                        residual = amber_k - quantum
                        residuals.append((idx, residual))
                    except (ValueError, IndexError):
                        continue

        if not residuals:
            return

        # Calculate statistics
        res_values = [r[1] for r in residuals]
        mean_res = sum(res_values) / len(res_values)
        max_res = max(res_values)
        min_res = min(res_values)
        rmse = (sum(r**2 for r in res_values) / len(res_values)) ** 0.5

        # Create ASCII plot
        plot_width = 50
        plot_height = 12

        # Determine scale
        abs_max = max(abs(max_res), abs(min_res))
        if abs_max < 0.001:
            abs_max = 1.0  # Avoid division by zero

        # Build plot grid
        grid = [[' ' for _ in range(plot_width + 10)] for _ in range(plot_height)]

        # Add zero line
        zero_row = plot_height // 2
        for col in range(5, plot_width + 5):
            grid[zero_row][col] = '-'

        # Add y-axis
        for row in range(plot_height):
            grid[row][4] = '|'

        # Add axis labels
        grid[0][0:4] = list(f"{abs_max:+.1f}"[:4].rjust(4))
        grid[zero_row][0:4] = list("0.0 ".rjust(4))
        grid[plot_height-1][0:4] = list(f"{-abs_max:+.1f}"[:4].rjust(4))

        # Plot data points
        n_points = len(residuals)
        for i, (idx, res) in enumerate(residuals):
            # X position: spread across plot width
            col = 5 + int((i / max(n_points - 1, 1)) * (plot_width - 1))
            # Y position: scale residual to plot height
            row = zero_row - int((res / abs_max) * (plot_height // 2 - 1))
            row = max(0, min(plot_height - 1, row))

            # Use different symbols based on magnitude
            if abs(res) < rmse * 0.5:
                symbol = '·'
            elif abs(res) < rmse:
                symbol = 'o'
            else:
                symbol = '*'

            if 0 <= col < len(grid[row]):
                grid[row][col] = symbol

        # Convert grid to string
        plot_lines = [''.join(row).rstrip() for row in grid]

        # Create output
        console.print(f"\n[bold]Energy Residuals (Amber+K - Quantum):[/bold]")
        console.print("[grey50]" + "─" * 60 + "[/grey50]")
        for line in plot_lines:
            console.print(f"[grey50]{line}[/grey50]")
        console.print(f"[grey50]     {'─' * plot_width}[/grey50]")
        console.print(f"[grey50]     Conformer (0 to {n_points-1})[/grey50]")
        console.print("[grey50]" + "─" * 60 + "[/grey50]")
        console.print(f"[cyan]  RMSE: {rmse:.2f} kcal/mol[/cyan]")
        console.print(f"[cyan]  Mean: {mean_res:.2f} kcal/mol[/cyan]")
        console.print(f"[cyan]  Range: {min_res:.2f} to {max_res:.2f} kcal/mol[/cyan]")

        # Add interpretation
        if rmse < 1.0:
            console.print(f"[green]  ✓ Excellent fit (RMSE < 1 kcal/mol)[/green]")
        elif rmse < 2.0:
            console.print(f"[green]  ✓ Good fit (RMSE < 2 kcal/mol)[/green]")
        elif rmse < 5.0:
            console.print(f"[yellow]  ⚠ Moderate fit (RMSE < 5 kcal/mol)[/yellow]")
        else:
            console.print(f"[red]  ✗ Poor fit (RMSE ≥ 5 kcal/mol) - parameters may be unreliable[/red]")

    except Exception as e:
        console.print(f"[grey50]Could not plot residuals: {e}[/grey50]")


def check_fitted_parameters(frcmod_file: str, console: Console) -> bool:
    """
    Check fitted parameters for unphysical values.

    Warns about:
    - Negative force constants
    - Extreme equilibrium values
    - Negative dihedral periodicities

    Args:
        frcmod_file: Path to the fitted frcmod file
        console: Rich console for output

    Returns:
        True if parameters look reasonable, False if problems detected
    """
    if not os.path.exists(frcmod_file):
        return True

    problems = []

    try:
        with open(frcmod_file, 'r') as f:
            lines = f.readlines()

        current_section = None

        for line in lines:
            stripped = line.strip()
            if not stripped or stripped.startswith('#') or stripped.startswith('Generated'):
                continue

            # Track sections
            if stripped in ['BOND', 'ANGL', 'ANGLE', 'DIHE', 'DIHEDRAL', 'IMPROPER']:
                current_section = stripped
                continue

            # Parse parameter values based on section
            parts = line.split()
            if len(parts) < 2:
                continue

            if current_section == 'BOND':
                # BOND format: "c -ns      4199.7513         1.2649"
                # First 5 chars are atom types, rest are values
                try:
                    param_name = line[:5].strip()
                    values = line[5:].split()
                    if len(values) >= 2:
                        force_const = float(values[0])
                        eq_dist = float(values[1])

                        if force_const < 0:
                            problems.append(f"Bond {param_name}: negative force constant ({force_const:.1f})")
                        elif force_const > 2000:
                            problems.append(f"Bond {param_name}: extreme force constant ({force_const:.1f}, typical: 200-600)")
                        if eq_dist < 0.5 or eq_dist > 3.0:
                            problems.append(f"Bond {param_name}: unusual equilibrium distance ({eq_dist:.3f} Å)")
                except (ValueError, IndexError):
                    pass

            elif current_section in ['ANGL', 'ANGLE']:
                # ANGLE format: "nc-c -ns       -47.2097       246.8317"
                # First 8 chars are atom types, rest are values
                try:
                    param_name = line[:8].strip()
                    values = line[8:].split()
                    if len(values) >= 2:
                        force_const = float(values[0])
                        eq_angle = float(values[1])

                        if force_const < 0:
                            problems.append(f"Angle {param_name}: negative force constant ({force_const:.1f})")
                        if eq_angle < 60 or eq_angle > 180:
                            problems.append(f"Angle {param_name}: unusual equilibrium angle ({eq_angle:.1f}°)")
                except (ValueError, IndexError):
                    pass

            elif current_section in ['DIHE', 'DIHEDRAL']:
                # DIHE format: "nc-c -ns-c     1       -11.0086       113.6335        -4.2852"
                # First 11 chars are atom types, rest are values
                try:
                    param_name = line[:11].strip()
                    values = line[11:].split()
                    if len(values) >= 4:
                        barrier = float(values[1])  # n, barrier, phase, periodicity
                        periodicity = float(values[3])

                        if barrier < -20 or barrier > 50:
                            problems.append(f"Dihedral {param_name}: extreme barrier ({barrier:.1f} kcal/mol)")
                        if periodicity < 0:
                            problems.append(f"Dihedral {param_name}: negative periodicity ({periodicity:.2f})")
                except (ValueError, IndexError):
                    pass

        # Report problems
        if problems:
            console.print(f"\n[bold red]⚠️  Warning: Fitted parameters may be unphysical![/bold red]")
            for problem in problems:
                console.print(f"[red]  • {problem}[/red]")
            console.print(f"\n[yellow]This often indicates:[/yellow]")
            console.print(f"[yellow]  • Too many parameters fitted relative to conformer diversity[/yellow]")
            console.print(f"[yellow]  • Poor conformer coverage of the parameter space[/yellow]")
            console.print(f"[yellow]  • Insufficient QM/MM energy correlation[/yellow]")
            console.print(f"\n[yellow]Consider:[/yellow]")
            console.print(f"[yellow]  • Fitting fewer parameters (focus on highest-penalty ones)[/yellow]")
            console.print(f"[yellow]  • Generating more diverse conformers[/yellow]")
            console.print(f"[yellow]  • Using the original parmchk2 parameters instead[/yellow]")
            return False

        return True

    except Exception as e:
        console.print(f"[grey50]Could not check parameters: {e}[/grey50]")
        return True


# =============================================================================
# Dependency Checking
# =============================================================================

def check_crest_availability() -> Tuple[bool, str]:
    """
    Check if CREST is available in the system PATH.

    Returns:
        Tuple of (available, message)
    """
    try:
        result = subprocess.run(
            ["crest", "--version"],
            capture_output=True,
            text=True,
            timeout=10
        )
        if result.returncode == 0:
            # Look for line containing "Version" in CREST output
            # Example: "Version 2.12,   Thu 19. Mai 16:32:32 CEST 2022"
            lines = result.stdout.strip().split('\n') if result.stdout else []
            version = "available"
            for line in lines:
                line = line.strip()
                if 'version' in line.lower() and not line.startswith('Using'):
                    # Extract just the version part
                    version = line
                    break
            return True, version
        else:
            return False, "CREST found but returned error"
    except FileNotFoundError:
        return False, "CREST not found in PATH"
    except subprocess.TimeoutExpired:
        return False, "CREST check timed out"
    except Exception as e:
        return False, f"Error checking CREST: {str(e)}"


def check_xtb_availability() -> Tuple[bool, str]:
    """
    Check if xTB is available in the system PATH.

    xTB is required by CREST for the GFN2-xTB calculations.

    Returns:
        Tuple of (available, message)
    """
    try:
        result = subprocess.run(
            ["xtb", "--version"],
            capture_output=True,
            text=True,
            timeout=10
        )
        if result.returncode == 0:
            # Look for line containing "xtb version" in output
            # Example: "* xtb version 6.5.1 (conda-forge) compiled by..."
            lines = result.stdout.strip().split('\n') if result.stdout else []
            version = "available"
            for line in lines:
                line = line.strip()
                if 'xtb version' in line.lower():
                    # Clean up the line (remove leading *)
                    version = line.lstrip('* ').strip()
                    break
            return True, version
        else:
            return False, "xTB found but returned error"
    except FileNotFoundError:
        return False, "xTB not found in PATH"
    except subprocess.TimeoutExpired:
        return False, "xTB check timed out"
    except Exception as e:
        return False, f"Error checking xTB: {str(e)}"


def check_paramfit_availability() -> Tuple[bool, str]:
    """
    Check if paramfit is available in the system PATH.

    Returns:
        Tuple of (available, message)
    """
    try:
        result = subprocess.run(
            ["paramfit", "--help"],
            capture_output=True,
            text=True,
            timeout=10
        )
        # paramfit --help returns 0 and prints usage
        if result.returncode == 0 or "Usage" in result.stdout or "usage" in result.stderr:
            return True, "paramfit available"
        else:
            return False, "paramfit found but returned unexpected output"
    except FileNotFoundError:
        return False, "paramfit not found in PATH"
    except subprocess.TimeoutExpired:
        return False, "paramfit check timed out"
    except Exception as e:
        return False, f"Error checking paramfit: {str(e)}"


def check_all_dependencies(console: Console) -> Dict[str, bool]:
    """
    Check all dependencies required for paramfit refinement.

    Args:
        console: Rich console for output

    Returns:
        Dict with dependency names as keys and availability as values
    """
    console.print("\n[bold]Checking dependencies for parameter refinement...[/bold]")

    dependencies = {}

    # Check CREST
    crest_ok, crest_msg = check_crest_availability()
    dependencies['crest'] = crest_ok
    if crest_ok:
        console.print(f"  [green]✓[/green] {crest_msg}")
    else:
        console.print(f"  [red]✗[/red] {crest_msg}")

    # Check xTB
    xtb_ok, xtb_msg = check_xtb_availability()
    dependencies['xtb'] = xtb_ok
    if xtb_ok:
        console.print(f"  [green]✓[/green] {xtb_msg}")
    else:
        console.print(f"  [red]✗[/red] {xtb_msg}")

    # Check paramfit
    paramfit_ok, paramfit_msg = check_paramfit_availability()
    dependencies['paramfit'] = paramfit_ok
    if paramfit_ok:
        console.print(f"  [green]✓[/green] {paramfit_msg}")
    else:
        console.print(f"  [red]✗[/red] {paramfit_msg}")

    # Summary
    if all(dependencies.values()):
        console.print("\n[green]All dependencies available for parameter refinement.[/green]")
    else:
        missing = [k for k, v in dependencies.items() if not v]
        console.print(f"\n[yellow]Missing dependencies: {', '.join(missing)}[/yellow]")
        console.print("[grey50]Install with: conda install conda-forge::crest conda-forge::xtb[/grey50]")

    return dependencies


# =============================================================================
# File Format Conversions
# =============================================================================

def mol2_to_xyz(mol2_file: str, xyz_file: str, console: Console) -> bool:
    """
    Convert MOL2 file to XYZ format for CREST.

    Args:
        mol2_file: Input MOL2 file path
        xyz_file: Output XYZ file path
        console: Rich console for output

    Returns:
        True if conversion successful
    """
    console.print(f"[cyan]Converting {mol2_file} to XYZ format...[/cyan]")

    try:
        atoms = []
        with open(mol2_file, 'r') as f:
            in_atom_section = False
            for line in f:
                if "@<TRIPOS>ATOM" in line:
                    in_atom_section = True
                    continue
                if line.startswith("@<TRIPOS>"):
                    in_atom_section = False
                    continue
                if in_atom_section and line.strip():
                    parts = line.split()
                    if len(parts) >= 6:
                        # parts: atom_id, atom_name, x, y, z, atom_type, ...
                        # Atom name format: element + index (e.g., N1, C12, Cl1, Fe1)
                        atom_name = parts[1]
                        x = float(parts[2])
                        y = float(parts[3])
                        z = float(parts[4])

                        # Extract element from atom name by taking leading letters
                        element = ""
                        for char in atom_name:
                            if char.isalpha():
                                element += char
                            else:
                                break
                        # Capitalize properly (first upper, rest lower)
                        element = element.capitalize()

                        atoms.append((element, x, y, z))

        if not atoms:
            console.print(f"[red]No atoms found in {mol2_file}[/red]")
            return False

        with open(xyz_file, 'w') as f:
            f.write(f"{len(atoms)}\n")
            f.write(f"Converted from {mol2_file}\n")
            for elem, x, y, z in atoms:
                f.write(f"{elem:2s}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")

        console.print(f"[green]✓ Created {xyz_file} with {len(atoms)} atoms[/green]")
        return True

    except Exception as e:
        console.print(f"[red]Error converting MOL2 to XYZ: {str(e)}[/red]")
        return False


def parse_crest_ensemble(ensemble_file: str, console: Console) -> List[Dict]:
    """
    Parse CREST conformer ensemble XYZ file.

    The CREST output format has multiple structures concatenated:
    - Line 1: Number of atoms
    - Line 2: Energy in Hartree (or comment)
    - Lines 3+: Atom coordinates
    - Repeat for each conformer

    Args:
        ensemble_file: Path to crest_conformers.xyz
        console: Rich console for output

    Returns:
        List of conformer dictionaries with 'id', 'energy_xtb', 'atoms'
    """
    console.print(f"[cyan]Parsing CREST ensemble from {ensemble_file}...[/cyan]")

    conformers = []

    try:
        with open(ensemble_file, 'r') as f:
            while True:
                # Read number of atoms
                natoms_line = f.readline()
                if not natoms_line:
                    break  # End of file

                try:
                    natoms = int(natoms_line.strip())
                except ValueError:
                    break  # Not a valid structure start

                # Read energy/comment line
                energy_line = f.readline().strip()
                try:
                    energy = float(energy_line)
                except ValueError:
                    energy = 0.0  # Comment line, no energy

                # Read atom coordinates
                atoms = []
                for _ in range(natoms):
                    atom_line = f.readline()
                    if not atom_line:
                        break
                    parts = atom_line.split()
                    if len(parts) >= 4:
                        element = parts[0]
                        x = float(parts[1])
                        y = float(parts[2])
                        z = float(parts[3])
                        atoms.append((element, x, y, z))

                if len(atoms) == natoms:
                    conformers.append({
                        'id': len(conformers),
                        'energy_xtb': energy,
                        'atoms': atoms
                    })

        console.print(f"[green]✓ Parsed {len(conformers)} conformers from ensemble[/green]")

        # Show energy range
        if conformers:
            energies = [c['energy_xtb'] for c in conformers if c['energy_xtb'] != 0.0]
            if energies:
                e_min = min(energies)
                e_max = max(energies)
                e_range_kcal = (e_max - e_min) * 627.509  # Hartree to kcal/mol
                console.print(f"[grey50]Energy range: {e_range_kcal:.2f} kcal/mol[/grey50]")

        return conformers

    except Exception as e:
        console.print(f"[red]Error parsing CREST ensemble: {str(e)}[/red]")
        return []


def write_amber_mdcrd(conformers: List[Dict], output_file: str, console: Console) -> bool:
    """
    Write conformers to AMBER trajectory (mdcrd) format.

    AMBER mdcrd format:
    - Line 1: Title
    - Following lines: 10 coordinates per line, F8.3 format
    - Each structure's coordinates are written sequentially

    Args:
        conformers: List of conformer dictionaries
        output_file: Output mdcrd file path
        console: Rich console for output

    Returns:
        True if successful
    """
    console.print(f"[cyan]Writing {len(conformers)} conformers to AMBER trajectory...[/cyan]")

    try:
        with open(output_file, 'w') as f:
            f.write("CREST conformers for paramfit\n")

            for conf in conformers:
                # Flatten coordinates
                coords = []
                for elem, x, y, z in conf['atoms']:
                    coords.extend([x, y, z])

                # Write 10 values per line, F8.3 format
                for i in range(0, len(coords), 10):
                    chunk = coords[i:i+10]
                    line = "".join(f"{c:8.3f}" for c in chunk)
                    f.write(line + "\n")

        console.print(f"[green]✓ Created {output_file}[/green]")
        return True

    except Exception as e:
        console.print(f"[red]Error writing mdcrd: {str(e)}[/red]")
        return False


# =============================================================================
# CREST Conformer Generation
# =============================================================================

def run_crest_conformer_search(
    xyz_file: str,
    charge: int,
    multiplicity: int,
    energy_window: float,
    nproc: int,
    console: Console,
    working_dir: str = None,
    use_solvation: bool = False,
    gfn_method: str = "gfn2"
) -> Optional[str]:
    """
    Run CREST conformer search using GFN-xTB.

    CREST generates diverse conformers by metadynamics sampling on the
    xTB potential energy surface. This is parameter-independent and
    provides unbiased conformational sampling.

    Args:
        xyz_file: Input XYZ file with starting structure
        charge: Molecular charge
        multiplicity: Spin multiplicity (1=singlet, 2=doublet, etc.)
        energy_window: Energy window in kcal/mol for conformer selection
        nproc: Number of CPU threads
        console: Rich console for output
        working_dir: Working directory (default: current)
        use_solvation: Whether to use ALPB implicit solvation (water)
        gfn_method: xTB method to use ("gfn1" or "gfn2")

    Returns:
        Path to crest_conformers.xyz if successful, None otherwise
    """
    # Build settings summary
    method_name = "GFN2-xTB" if gfn_method == "gfn2" else "GFN1-xTB"
    solvation_str = "ALPB water" if use_solvation else "gas phase"

    # Display explanation
    console.print(Panel(
        "[bold cyan]CREST Conformer Search[/bold cyan]\n\n"
        "[bold]What this does:[/bold]\n"
        f"  CREST uses metadynamics with {method_name} (semiempirical QM) to\n"
        "  explore conformational space and find low-energy conformers.\n\n"
        "[bold]Why xTB instead of MM?[/bold]\n"
        "  Using xTB avoids bias from the MM parameters we're trying to\n"
        "  improve. The conformers are generated on a QM-level surface,\n"
        "  independent of GAFF2 parameters.\n\n"
        f"[bold]Settings:[/bold]\n"
        f"  Charge: {charge}\n"
        f"  Multiplicity: {multiplicity}\n"
        f"  Energy window: {energy_window} kcal/mol\n"
        f"  Method: {method_name}\n"
        f"  Solvation: {solvation_str}\n"
        f"  Threads: {nproc}",
        title="Conformer Generation",
        border_style="cyan",
        expand=False
    ))

    # Convert multiplicity to unpaired electrons for xTB
    uhf = multiplicity - 1

    # Build command
    cmd = [
        "crest", xyz_file,
        f"--{gfn_method}",
        "--chrg", str(charge),
        "--uhf", str(uhf),
        "--ewin", str(energy_window),
        "-T", str(nproc)
    ]

    # Add solvation if requested (helps convergence for charged molecules)
    if use_solvation:
        cmd.extend(["--alpb", "water"])

    console.print(f"\n[yellow]Running CREST command:[/yellow]")
    console.print(f"[grey50]{' '.join(cmd)}[/grey50]")
    console.print(f"\n[grey50]This may take several minutes...[/grey50]")

    try:
        # Run CREST
        cwd = working_dir if working_dir else os.getcwd()

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
            transient=True
        ) as progress:
            task = progress.add_task("Running CREST conformer search...", total=None)

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=cwd
                # No timeout - CREST can take hours for complex molecules
            )

        # Check for output file
        ensemble_file = os.path.join(cwd, "crest_conformers.xyz")

        if result.returncode == 0 and os.path.exists(ensemble_file):
            console.print(f"[green]✓ CREST completed successfully[/green]")

            # Parse output for summary
            if "number of unique conformers" in result.stdout:
                for line in result.stdout.split('\n'):
                    if "number of unique conformers" in line.lower():
                        console.print(f"[grey50]{line.strip()}[/grey50]")
                        break

            return ensemble_file
        else:
            console.print(f"[red]✗ CREST failed[/red]")
            if result.stderr:
                console.print(f"[red]Error: {result.stderr[:500]}[/red]")
            return None

    except Exception as e:
        console.print(f"[red]✗ Error running CREST: {str(e)}[/red]")
        return None


# =============================================================================
# Gaussian Single-Point Calculations
# =============================================================================

def write_gaussian_sp_inputs(
    conformers: List[Dict],
    output_dir: str,
    mol_name: str,
    charge: int,
    multiplicity: int,
    method: str,
    memory: str,
    nproc: int,
    console: Console
) -> List[str]:
    """
    Write Gaussian single-point input files for each conformer.

    Args:
        conformers: List of conformer dictionaries
        output_dir: Directory for output files
        mol_name: Molecule name for file naming
        charge: Molecular charge
        multiplicity: Spin multiplicity
        method: Gaussian route line keywords (e.g., "B3LYP/6-31G*" or
                "wB97XD/def2-TZVP EmpiricalDispersion=GD3BJ")
        memory: Memory allocation (e.g., "4GB")
        nproc: Number of processors
        console: Rich console for output

    Returns:
        List of generated input file paths (in order!)
    """
    console.print(Panel(
        "[bold cyan]Gaussian Single-Point Inputs[/bold cyan]\n\n"
        "[bold]What this does:[/bold]\n"
        "  Creates Gaussian input files for single-point energy calculations\n"
        "  on each conformer. These energies will be used to fit parameters.\n\n"
        "[bold]Why single-point (no optimization)?[/bold]\n"
        "  We want energies at the exact CREST geometries. Optimization would\n"
        "  change the structures and break the correspondence with the mdcrd.\n\n"
        f"[bold]Settings:[/bold]\n"
        f"  Route keywords: {method} SP\n"
        f"  Memory: {memory}\n"
        f"  Processors: {nproc}\n"
        f"  Conformers: {len(conformers)}",
        title="QM Energy Calculations",
        border_style="cyan",
        expand=False
    ))

    os.makedirs(output_dir, exist_ok=True)

    input_files = []

    for conf in conformers:
        conf_id = conf['id']
        # Use zero-padded numbering for proper sorting
        filename = f"{mol_name}_conf_{conf_id:03d}.gjf"
        filepath = os.path.join(output_dir, filename)

        with open(filepath, 'w') as f:
            f.write(f"%chk={mol_name}_conf_{conf_id:03d}.chk\n")
            f.write(f"%nproc={nproc}\n")
            f.write(f"%mem={memory}\n")
            f.write(f"# {method} SP\n")  # SP = Single Point
            f.write("\n")
            f.write(f"{mol_name} conformer {conf_id} single-point energy\n")
            f.write("\n")
            f.write(f"{charge} {multiplicity}\n")

            for elem, x, y, z in conf['atoms']:
                f.write(f"{elem:2s}  {x:14.8f}  {y:14.8f}  {z:14.8f}\n")

            f.write("\n")

        input_files.append(filepath)

    console.print(f"[green]✓ Created {len(input_files)} Gaussian input files in {output_dir}/[/green]")
    console.print(f"[grey50]Files: {mol_name}_conf_000.gjf to {mol_name}_conf_{len(conformers)-1:03d}.gjf[/grey50]")

    return input_files


def write_gaussian_batch_script(
    output_dir: str,
    mol_name: str,
    n_conformers: int,
    parallel_jobs: int,
    console: Console
) -> str:
    """
    Generate a bash script to run Gaussian jobs in parallel batches.

    Args:
        output_dir: Directory containing Gaussian input files
        mol_name: Molecule name used in file naming
        n_conformers: Number of conformer input files
        parallel_jobs: Number of jobs to run simultaneously
        console: Rich console for output

    Returns:
        Path to the generated run script
    """
    script_path = os.path.join(output_dir, "run_gaussian.sh")

    script_content = f'''#!/bin/bash
# Gaussian batch run script for {mol_name} conformers
# Generated by ProPrep
#
# This script runs {n_conformers} Gaussian jobs in batches of {parallel_jobs}
#
# Usage: cd {output_dir} && bash run_gaussian.sh
#        or: bash {script_path}

# Configuration
PARALLEL_JOBS={parallel_jobs}
GAUSSIAN_CMD="g16"  # Change to g09 if using Gaussian 09

# Change to script directory
cd "$(dirname "$0")"

echo "=============================================="
echo "Gaussian Batch Runner - {mol_name}"
echo "=============================================="
echo "Total jobs: {n_conformers}"
echo "Parallel jobs: $PARALLEL_JOBS"
echo ""

# Count total and completed jobs
TOTAL={n_conformers}
COMPLETED=0
RUNNING=0

# Function to check if a job is done (log file exists and contains "Normal termination")
is_complete() {{
    local logfile="$1"
    if [ -f "$logfile" ] && grep -q "Normal termination" "$logfile" 2>/dev/null; then
        return 0
    fi
    return 1
}}

# Get list of all input files
INPUT_FILES=({mol_name}_conf_*.gjf)

# Count already completed
for f in "${{INPUT_FILES[@]}}"; do
    logfile="${{f%.gjf}}.log"
    if is_complete "$logfile"; then
        ((COMPLETED++))
    fi
done

if [ $COMPLETED -gt 0 ]; then
    echo "Found $COMPLETED already completed jobs"
fi

if [ $COMPLETED -eq $TOTAL ]; then
    echo "All jobs already complete!"
    exit 0
fi

echo "Starting calculations..."
echo ""

# Run jobs in parallel batches
PIDS=()
for f in "${{INPUT_FILES[@]}}"; do
    logfile="${{f%.gjf}}.log"

    # Skip if already complete
    if is_complete "$logfile"; then
        continue
    fi

    # Wait if we have too many running jobs
    while [ ${{#PIDS[@]}} -ge $PARALLEL_JOBS ]; do
        # Check which jobs have finished
        NEW_PIDS=()
        for pid in "${{PIDS[@]}}"; do
            if kill -0 "$pid" 2>/dev/null; then
                NEW_PIDS+=("$pid")
            fi
        done
        PIDS=("${{NEW_PIDS[@]}}")

        if [ ${{#PIDS[@]}} -ge $PARALLEL_JOBS ]; then
            sleep 5
        fi
    done

    # Start new job
    echo "Starting: $f"
    $GAUSSIAN_CMD "$f" > "$logfile" 2>&1 &
    PIDS+=($!)
done

# Wait for remaining jobs
echo ""
echo "Waiting for remaining jobs to complete..."
for pid in "${{PIDS[@]}}"; do
    wait "$pid"
done

# Final status check
echo ""
echo "=============================================="
echo "Final Status"
echo "=============================================="
COMPLETED=0
FAILED=0
for f in "${{INPUT_FILES[@]}}"; do
    logfile="${{f%.gjf}}.log"
    if is_complete "$logfile"; then
        ((COMPLETED++))
        echo "✓ $logfile"
    elif [ -f "$logfile" ]; then
        ((FAILED++))
        echo "✗ $logfile (check for errors)"
    else
        ((FAILED++))
        echo "✗ $logfile (not created)"
    fi
done

echo ""
echo "Completed: $COMPLETED / $TOTAL"
if [ $FAILED -gt 0 ]; then
    echo "Failed/Missing: $FAILED"
    echo ""
    echo "Check failed log files for errors."
    exit 1
fi

echo ""
echo "All Gaussian calculations completed successfully!"
'''

    with open(script_path, 'w') as f:
        f.write(script_content)

    # Make executable
    os.chmod(script_path, 0o755)

    console.print(f"[green]✓ Created batch run script: {script_path}[/green]")

    return script_path


def extract_gaussian_energies(
    log_dir: str,
    mol_name: str,
    n_conformers: int,
    output_file: str,
    console: Console
) -> Optional[List[float]]:
    """
    Extract SCF energies from Gaussian log files.

    CRITICAL: Files must be processed in order (conf_000, conf_001, etc.)
    to maintain correspondence with the mdcrd trajectory.

    Args:
        log_dir: Directory containing Gaussian log files
        mol_name: Molecule name used in file naming
        n_conformers: Expected number of conformers
        output_file: Output file for energies (one per line)
        console: Rich console for output

    Returns:
        List of energies in Hartree, or None if extraction failed
    """
    console.print(f"[cyan]Extracting energies from Gaussian log files...[/cyan]")

    energies = []
    missing_files = []
    failed_extractions = []

    for conf_id in range(n_conformers):
        log_file = os.path.join(log_dir, f"{mol_name}_conf_{conf_id:03d}.log")

        if not os.path.exists(log_file):
            missing_files.append(log_file)
            continue

        energy = None
        try:
            with open(log_file, 'r') as f:
                for line in f:
                    if "SCF Done" in line:
                        # Format: SCF Done:  E(RB3LYP) =  -XXX.XXXXXXXX     A.U.
                        parts = line.split()
                        for i, part in enumerate(parts):
                            if part == "=":
                                energy = float(parts[i + 1])
                                break
        except Exception as e:
            failed_extractions.append((log_file, str(e)))
            continue

        if energy is None:
            failed_extractions.append((log_file, "No SCF Done found"))
        else:
            energies.append(energy)

    # Report issues
    if missing_files:
        console.print(f"[red]Missing {len(missing_files)} log files:[/red]")
        for f in missing_files[:5]:
            console.print(f"[red]  - {f}[/red]")
        if len(missing_files) > 5:
            console.print(f"[red]  ... and {len(missing_files) - 5} more[/red]")
        return None

    if failed_extractions:
        console.print(f"[red]Failed to extract energy from {len(failed_extractions)} files:[/red]")
        for f, err in failed_extractions[:5]:
            console.print(f"[red]  - {f}: {err}[/red]")
        return None

    if len(energies) != n_conformers:
        console.print(f"[red]Expected {n_conformers} energies, got {len(energies)}[/red]")
        return None

    # Write energies file
    with open(output_file, 'w') as f:
        for e in energies:
            f.write(f"{e}\n")

    console.print(f"[green]✓ Extracted {len(energies)} energies to {output_file}[/green]")

    # Show energy range
    e_min = min(energies)
    e_max = max(energies)
    e_range_kcal = (e_max - e_min) * 627.509
    console.print(f"[grey50]Energy range: {e_range_kcal:.2f} kcal/mol[/grey50]")

    return energies


# =============================================================================
# Paramfit Integration
# =============================================================================

def run_paramfit_k_fitting(
    prmtop_file: str,
    mdcrd_file: str,
    energies_file: str,
    n_structures: int,
    console: Console
) -> Optional[float]:
    """
    Run paramfit to determine K (QM-MM energy offset).

    K must be fit first before fitting any other parameters.
    It represents the intrinsic difference between QM and MM
    energy origins.

    Args:
        prmtop_file: AMBER topology file
        mdcrd_file: AMBER trajectory with conformers
        energies_file: QM energies (one per line, Hartree)
        n_structures: Number of structures
        console: Rich console for output

    Returns:
        Fitted K value, or None if fitting failed
    """
    console.print(Panel(
        "[bold cyan]Fitting K (Energy Offset Constant)[/bold cyan]\n\n"
        "[bold]What is K?[/bold]\n"
        "  K represents the intrinsic difference between QM and MM energy\n"
        "  scales. QM energies include all electrons, while MM only models\n"
        "  classical interactions, so they have different 'zero points'.\n\n"
        "[bold]Why fit K first?[/bold]\n"
        "  Paramfit minimizes (E_MM - E_QM + K)^2 over all conformers.\n"
        "  If K is wrong, the algorithm will try to adjust K instead of\n"
        "  the force field parameters, giving poor results.",
        title="K Fitting",
        border_style="cyan",
        expand=False
    ))

    # Create job control file
    job_file = "fit_K.in"
    with open(job_file, 'w') as f:
        f.write(f"""# Paramfit job control file - K fitting only
RUNTYPE=FIT
PARAMETERS_TO_FIT=K_ONLY
FUNC_TO_FIT=SUM_SQUARES_AMBER_STANDARD
QM_ENERGY_UNITS=HARTREE
NSTRUCTURES={n_structures}
COORDINATE_FORMAT=TRAJECTORY
ALGORITHM=SIMPLEX
""")

    console.print(f"\n[yellow]Running paramfit to fit K...[/yellow]")

    cmd = [
        "paramfit",
        "-i", job_file,
        "-p", prmtop_file,
        "-c", mdcrd_file,
        "-q", energies_file
    ]

    console.print(f"[grey50]{' '.join(cmd)}[/grey50]")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
        )

        # Parse output for K value
        k_value = None
        r_squared = None

        for line in result.stdout.split('\n'):
            if 'K =' in line or 'K=' in line:
                try:
                    parts = line.split('=')
                    k_value = float(parts[-1].strip().split()[0])
                except (ValueError, IndexError):
                    pass
            if 'R^2' in line or 'R2' in line:
                try:
                    parts = line.replace('R^2', 'R2').split('=')
                    r_squared = float(parts[-1].strip().split()[0])
                except (ValueError, IndexError):
                    pass

        if k_value is not None:
            console.print(f"[green]✓ K fitting complete[/green]")
            console.print(f"[cyan]  K = {k_value:.6f} kcal/mol[/cyan]")
            if r_squared is not None:
                console.print(f"[cyan]  R² = {r_squared:.4f}[/cyan]")
            return k_value
        else:
            console.print(f"[red]✗ Could not parse K value from paramfit output[/red]")
            if result.stdout:
                console.print(f"[grey50]Output: {result.stdout[:500]}[/grey50]")
            return None

    except Exception as e:
        console.print(f"[red]✗ Error running paramfit: {str(e)}[/red]")
        return None


def run_paramfit_parameter_fitting(
    prmtop_file: str,
    mdcrd_file: str,
    energies_file: str,
    params_file: str,
    n_structures: int,
    k_value: float,
    output_frcmod: str,
    console: Console
) -> bool:
    """
    Run paramfit to fit force field parameters.

    Args:
        prmtop_file: AMBER topology file
        mdcrd_file: AMBER trajectory with conformers
        energies_file: QM energies file
        params_file: Parameter specification file from SET_PARAMS
        n_structures: Number of structures
        k_value: Pre-determined K value
        output_frcmod: Output frcmod file path
        console: Rich console for output

    Returns:
        True if fitting succeeded
    """
    console.print(Panel(
        "[bold cyan]Fitting Force Field Parameters[/bold cyan]\n\n"
        "[bold]What this does:[/bold]\n"
        "  Paramfit adjusts the selected parameters (bonds, angles, dihedrals)\n"
        "  to minimize the difference between MM and QM energies across all\n"
        "  conformers using a genetic algorithm.\n\n"
        "[bold]Output:[/bold]\n"
        f"  The fitted parameters will be saved to: {output_frcmod}",
        title="Parameter Fitting",
        border_style="cyan",
        expand=False
    ))

    # Create job control file
    job_file = "fit_params.in"
    energy_output = output_frcmod.replace('.frcmod', '_energies.dat')

    with open(job_file, 'w') as f:
        f.write(f"""# Paramfit job control file - parameter fitting
RUNTYPE=FIT
PARAMETERS_TO_FIT=LOAD
PARAMETER_FILE_NAME={params_file}
FUNC_TO_FIT=SUM_SQUARES_AMBER_STANDARD
QM_ENERGY_UNITS=HARTREE
NSTRUCTURES={n_structures}
COORDINATE_FORMAT=TRAJECTORY
K={k_value}

# Genetic algorithm settings
ALGORITHM=GENETIC
OPTIMIZATIONS=500
GENERATIONS_TO_CONV=20
GENERATIONS_TO_SIMPLEX=5
GENERATIONS_WITHOUT_SIMPLEX=5

# Output
WRITE_FRCMOD={output_frcmod}
WRITE_ENERGY={energy_output}
""")

    console.print(f"\n[yellow]Running paramfit parameter fitting...[/yellow]")
    console.print(f"[grey50]This may take several minutes...[/grey50]")

    cmd = [
        "paramfit",
        "-i", job_file,
        "-p", prmtop_file,
        "-c", mdcrd_file,
        "-q", energies_file,
        "-v", "MEDIUM"
    ]

    console.print(f"[grey50]{' '.join(cmd)}[/grey50]")

    try:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
            transient=True
        ) as progress:
            task = progress.add_task("Fitting parameters...", total=None)

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True
            )

        # Check for output file
        if os.path.exists(output_frcmod):
            console.print(f"[green]✓ Parameter fitting complete[/green]")
            console.print(f"[green]  Output: {output_frcmod}[/green]")

            # Parse for R² value
            for line in result.stdout.split('\n'):
                if 'R^2' in line or 'R2' in line:
                    console.print(f"[cyan]  {line.strip()}[/cyan]")
                    break

            # Plot energy residuals if energies file exists
            if os.path.exists(energy_output):
                plot_energy_residuals(energy_output, console)

            # Check for unphysical parameters
            check_fitted_parameters(output_frcmod, console)

            return True
        else:
            console.print(f"[red]✗ Paramfit did not create output file[/red]")
            if result.stderr:
                console.print(f"[red]Error: {result.stderr[:500]}[/red]")
            return False

    except Exception as e:
        console.print(f"[red]✗ Error running paramfit: {str(e)}[/red]")
        return False


def run_paramfit_set_params_wizard(prmtop_file: str, params_file: str, console: Console) -> bool:
    """
    Run paramfit in SET_PARAMS mode (interactive wizard).

    This launches paramfit's interactive mode where the user can
    select which parameters to fit.

    Args:
        prmtop_file: AMBER topology file
        params_file: Output file for parameter specification
        console: Rich console for output

    Returns:
        True if parameter file was created
    """
    console.print(Panel(
        "[bold cyan]Interactive Parameter Selection[/bold cyan]\n\n"
        "Paramfit will now launch in interactive mode.\n"
        "You will be asked which bonds, angles, and dihedrals to fit.\n\n"
        "[bold]Tips:[/bold]\n"
        "  - Focus on parameters with high penalty scores\n"
        "  - Fitting fewer parameters gives more reliable results\n"
        "  - Dihedrals are most commonly problematic\n\n"
        f"[grey50]Parameter specification will be saved to: {params_file}[/grey50]",
        title="Parameter Selection",
        border_style="cyan",
        expand=False
    ))

    # Create minimal job file for SET_PARAMS mode
    job_file = "set_params.in"
    with open(job_file, 'w') as f:
        f.write(f"""RUNTYPE=SET_PARAMS
PARAMETER_FILE_NAME={params_file}
""")

    console.print(f"\n[yellow]Launching paramfit interactive mode...[/yellow]")
    console.print(f"[grey50]Follow the prompts to select parameters to fit.[/grey50]\n")

    try:
        # Run interactively (not capturing output)
        result = subprocess.run(
            ["paramfit", "-i", job_file, "-p", prmtop_file]
        )

        if os.path.exists(params_file):
            console.print(f"\n[green]✓ Parameter file created: {params_file}[/green]")
            return True
        else:
            console.print(f"\n[red]✗ Parameter file was not created[/red]")
            return False

    except Exception as e:
        console.print(f"[red]✗ Error: {str(e)}[/red]")
        return False


def run_paramfit_set_params_automated(
    prmtop_file: str,
    params_file: str,
    selected_params: List[Tuple[str, float, str, str]],
    console: Console,
    force_constants_only: bool = True
) -> bool:
    """
    Run paramfit SET_PARAMS mode with automated responses based on pre-selected parameters.

    This uses pty to interact with paramfit's prompts, answering based on the
    parameters the user selected from the penalty table earlier.

    Args:
        prmtop_file: AMBER topology file
        params_file: Output file for parameter specification
        selected_params: List of (param_name, score, status, section) tuples
        console: Rich console for output
        force_constants_only: If True, only fit force constants (KR, KT, KP).
                             If False, also fit equilibrium values (REQ, THEQ, NP, PHASE).

    Returns:
        True if parameter file was created successfully
    """
    # Categorize selected parameters by section
    selected_by_section: Dict[str, Set[str]] = {
        'BOND': set(),
        'ANGLE': set(),
        'DIHE': set(),
    }

    for param_name, score, status, section in selected_params:
        normalized = normalize_param_atoms(param_name)
        if section == 'BOND':
            selected_by_section['BOND'].add(normalized)
        elif section == 'ANGLE':
            selected_by_section['ANGLE'].add(normalized)
        elif section in ['DIHE', 'DIHEDRAL', 'IMPROPER']:
            selected_by_section['DIHE'].add(normalized)

    has_bonds = len(selected_by_section['BOND']) > 0
    has_angles = len(selected_by_section['ANGLE']) > 0
    has_dihedrals = len(selected_by_section['DIHE']) > 0

    fit_mode_desc = "Force constants only (KR, KT, KP)" if force_constants_only else "All properties (force constants + equilibrium values)"

    console.print(Panel(
        "[bold cyan]Automated Parameter Selection[/bold cyan]\n\n"
        "ProPrep will automatically configure paramfit based on your earlier\n"
        "parameter selections from the penalty table.\n\n"
        f"[bold]Parameters to fit:[/bold]\n"
        f"  Bonds: {len(selected_by_section['BOND'])} selected\n"
        f"  Angles: {len(selected_by_section['ANGLE'])} selected\n"
        f"  Dihedrals: {len(selected_by_section['DIHE'])} selected\n\n"
        f"[bold]Fit mode:[/bold] {fit_mode_desc}",
        title="Parameter Selection",
        border_style="cyan",
        expand=False
    ))

    # Create job file for SET_PARAMS mode
    job_file = "set_params.in"
    with open(job_file, 'w') as f:
        f.write(f"""RUNTYPE=SET_PARAMS
PARAMETER_FILE_NAME={params_file}
""")

    console.print(f"\n[yellow]Running paramfit SET_PARAMS with automated responses...[/yellow]")

    # Build command
    cmd = ["paramfit", "-i", job_file, "-p", prmtop_file]

    try:
        # Use pty to interact with paramfit
        master_fd, slave_fd = pty.openpty()

        process = subprocess.Popen(
            cmd,
            stdin=slave_fd,
            stdout=slave_fd,
            stderr=slave_fd,
            close_fds=True
        )

        os.close(slave_fd)

        # Buffer for accumulating output
        output_buffer = ""
        params_selected = 0
        params_skipped = 0

        def send_response(response: str):
            """Send a response to paramfit."""
            os.write(master_fd, (response + "\n").encode())

        def read_output(timeout: float = 2.0) -> str:
            """Read available output from paramfit."""
            nonlocal output_buffer
            ready, _, _ = select.select([master_fd], [], [], timeout)
            if ready:
                try:
                    data = os.read(master_fd, 4096).decode('utf-8', errors='replace')
                    output_buffer += data
                    return data
                except OSError:
                    return ""
            return ""

        def wait_for_prompt(timeout: float = 5.0) -> str:
            """Wait for and return prompt text."""
            nonlocal output_buffer
            start_buffer = output_buffer
            elapsed = 0.0
            while elapsed < timeout:
                read_output(0.5)
                elapsed += 0.5
                # Check if we have a y/n prompt
                if "(y/n):" in output_buffer[len(start_buffer):]:
                    return output_buffer
                # Check if done
                if "Filename to be written" in output_buffer:
                    return output_buffer
            return output_buffer

        # Main interaction loop
        max_iterations = 500  # Safety limit
        iteration = 0

        while iteration < max_iterations:
            iteration += 1
            wait_for_prompt()

            # Check for completion
            if "Filename to be written" in output_buffer or process.poll() is not None:
                break

            # Parse the latest prompt from the buffer
            # Look for the most recent (y/n): prompt
            lines = output_buffer.split('\n')

            # Find prompts we haven't answered yet
            prompt_line = None
            for line in reversed(lines):
                if "(y/n):" in line:
                    prompt_line = line.strip()
                    break

            if not prompt_line:
                # No prompt found, wait a bit more
                read_output(0.5)
                continue

            # Determine response based on prompt type
            response = "n"  # Default to no

            # Global type questions
            if "Would you like to optimise any bond parameters?" in prompt_line:
                response = "y" if has_bonds else "n"
            elif "Would you like to optimise any angle parameters?" in prompt_line:
                response = "y" if has_angles else "n"
            elif "Would you like to optimise any dihedral parameters?" in prompt_line:
                response = "y" if has_dihedrals else "n"

            # Sub-property questions - handle based on force_constants_only setting
            # Force constants: KR (bond), KT (angle), KP (dihedral)
            # Equilibrium values: REQ (bond), THEQ (angle), NP/PHASE (dihedral)
            elif "Any bond KR?" in prompt_line:
                response = "y" if has_bonds else "n"
            elif "Any bond REQ?" in prompt_line:
                response = "y" if (has_bonds and not force_constants_only) else "n"
            elif "Any angle KT?" in prompt_line:
                response = "y" if has_angles else "n"
            elif "Any angle THEQ?" in prompt_line:
                response = "y" if (has_angles and not force_constants_only) else "n"
            elif "Any dihedral KP?" in prompt_line:
                response = "y" if has_dihedrals else "n"
            elif "Any dihedral NP?" in prompt_line or "Any dihedral PHASE?" in prompt_line:
                response = "y" if (has_dihedrals and not force_constants_only) else "n"

            # Individual parameter prompts
            elif "Fit Parameter:" in prompt_line:
                # Parse: "Fit Parameter: ANGLE (nc-c -ns) KT? (y/n):"
                # or: "Fit Parameter: DIHEDRAL (o -c -ns-hn) term 1 KP? (y/n):"
                match = re.search(r'Fit Parameter:\s*(BOND|ANGLE|DIHEDRAL)\s*\(([^)]+)\)', prompt_line)
                if match:
                    param_type = match.group(1)
                    atom_types = match.group(2)
                    normalized_atoms = normalize_param_atoms(atom_types)

                    # Check if this parameter was selected
                    if param_type == 'BOND' and normalized_atoms in selected_by_section['BOND']:
                        response = "y"
                        params_selected += 1
                    elif param_type == 'ANGLE' and normalized_atoms in selected_by_section['ANGLE']:
                        response = "y"
                        params_selected += 1
                    elif param_type == 'DIHEDRAL' and normalized_atoms in selected_by_section['DIHE']:
                        response = "y"
                        params_selected += 1
                    else:
                        response = "n"
                        params_skipped += 1

            # Send response and clear relevant part of buffer
            send_response(response)

            # Clear the processed prompt from consideration
            # by marking how far we've processed
            output_buffer += f"\n[RESPONDED: {response}]\n"

        # Close pty and wait for process
        os.close(master_fd)
        process.wait(timeout=5)

        # Check if parameter file was created
        if os.path.exists(params_file):
            console.print(f"\n[green]✓ Parameter file created: {params_file}[/green]")
            console.print(f"[green]  Parameters selected: {params_selected}[/green]")
            console.print(f"[grey50]  Parameters skipped: {params_skipped}[/grey50]")
            return True
        else:
            console.print(f"\n[red]✗ Parameter file was not created[/red]")
            console.print(f"[grey50]Paramfit output may contain errors[/grey50]")
            return False

    except Exception as e:
        console.print(f"[red]✗ Error in automated parameter selection: {str(e)}[/red]")
        return False


# =============================================================================
# High-Level Workflow Functions
# =============================================================================

def offer_paramfit_refinement(
    penalties: List[Tuple[str, float, str, str]],
    mol_name: str,
    mol2_file: str,
    prmtop_file: str,
    frcmod_file: str,
    charge: int,
    multiplicity: int,
    console: Console,
    interactive: bool = True,
    processor=None
) -> Dict[str, Any]:
    """
    Main entry point for paramfit refinement workflow.

    Called by small_molecule_parameterizer after penalty analysis.

    Args:
        penalties: List of (param_name, score, status, section) from analyze_frcmod_penalties
        mol_name: Molecule/residue name
        mol2_file: MOL2 file with charges and atom types
        prmtop_file: AMBER topology file
        frcmod_file: Original frcmod file from parmchk2
        charge: Molecular charge
        multiplicity: Spin multiplicity
        console: Rich console for output
        interactive: Whether to run in interactive mode
        processor: Optional processor for session recording context

    Returns:
        Dict with results of refinement workflow
    """
    results = {
        'refinement_attempted': False,
        'refinement_success': False,
        'fitted_frcmod': None,
        'n_conformers': 0,
        'message': ''
    }

    if not interactive:
        results['message'] = "Paramfit refinement requires interactive mode"
        return results

    # Categorize penalties
    attn_params = [p for p in penalties if p[2] == 'ATTN']
    high_penalty = [p for p in penalties if p[1] > 10 and p[2] != 'ATTN']
    medium_penalty = [p for p in penalties if 1 < p[1] <= 10]
    low_penalty = [p for p in penalties if p[1] <= 1 and p[2] != 'ATTN']

    # Build penalty summary table
    console.print("\n[bold cyan]Parameter Penalty Summary[/bold cyan]")

    table = Table(expand=False, show_header=True)
    table.add_column("Category", style="bold")
    table.add_column("Count", justify="right")
    table.add_column("Parameters")

    if attn_params:
        param_list = ", ".join([p[0] for p in attn_params[:5]])
        if len(attn_params) > 5:
            param_list += f" (+{len(attn_params) - 5} more)"
        table.add_row(
            "[red bold]⚠️  ATTN (needs attention)[/red bold]",
            f"[red]{len(attn_params)}[/red]",
            f"[red]{param_list}[/red]"
        )

    if high_penalty:
        param_list = ", ".join([p[0] for p in high_penalty[:5]])
        if len(high_penalty) > 5:
            param_list += f" (+{len(high_penalty) - 5} more)"
        table.add_row(
            "[yellow]High penalty (>10)[/yellow]",
            f"[yellow]{len(high_penalty)}[/yellow]",
            param_list
        )

    if medium_penalty:
        param_list = ", ".join([p[0] for p in medium_penalty[:5]])
        if len(medium_penalty) > 5:
            param_list += f" (+{len(medium_penalty) - 5} more)"
        table.add_row(
            "Medium penalty (1-10)",
            str(len(medium_penalty)),
            param_list
        )

    if low_penalty:
        table.add_row(
            "[green]Low penalty (<1)[/green]",
            f"[green]{len(low_penalty)}[/green]",
            "[grey50]acceptable[/grey50]"
        )

    console.print(table)

    # Highlight ATTN warnings specifically
    if attn_params:
        console.print(Panel(
            "[red bold]ATTN parameters detected![/red bold]\n\n"
            "ATTN (Attention Required) means parmchk2 could not find any similar\n"
            "atom type to estimate these parameters. The values are essentially\n"
            "guesses and [bold]should be refined[/bold] for reliable simulations.\n\n"
            f"Parameters requiring attention: [yellow]{', '.join([p[0] for p in attn_params])}[/yellow]",
            border_style="red",
            expand=False
        ))

    # Ask user if they want to refine
    console.print(Panel(
        "[bold cyan]What is paramfit refinement?[/bold cyan]\n\n"
        "An advanced workflow that fits force field parameters to QM energies\n"
        "across multiple molecular conformations. This improves accuracy for\n"
        "parameters that parmchk2 could only estimate by analogy.\n\n"
        "[bold]Workflow:[/bold]\n"
        "  1. CREST generates diverse conformers (using xTB - no parameters needed)\n"
        "  2. Gaussian calculates QM energies for each conformer\n"
        "  3. paramfit optimizes parameters to reproduce QM energy differences\n\n"
        "[bold]Requirements:[/bold]\n"
        "  • CREST and xTB (for conformer generation)\n"
        "  • Gaussian (for QM single-point energies)\n"
        "  • Additional computation time\n\n"
        "[grey50]This is optional - the current parameters may work for your application.\n"
        "Consider refinement especially if ATTN warnings are present.[/grey50]",
        title="Advanced Parameter Refinement",
        border_style="cyan",
        expand=False
    ))

    proceed = confirm_with_context(
        processor,
        "Would you like to refine parameters using paramfit?",
        default=False,
        module="Paramfit Refinement",
        description="Proceed with paramfit refinement",
    )

    if not proceed:
        results['message'] = "User declined paramfit refinement"
        return results

    results['refinement_attempted'] = True

    # Check dependencies
    deps = check_all_dependencies(console)

    if not all(deps.values()):
        console.print("\n[red]Cannot proceed without required dependencies.[/red]")
        results['message'] = "Missing dependencies"
        return results

    # Run the refinement workflow
    return run_paramfit_refinement_workflow(
        mol_name=mol_name,
        mol2_file=mol2_file,
        prmtop_file=prmtop_file,
        frcmod_file=frcmod_file,
        charge=charge,
        multiplicity=multiplicity,
        console=console,
        processor=processor
    )


def run_paramfit_refinement_workflow(
    mol_name: str,
    mol2_file: str,
    prmtop_file: str,
    frcmod_file: str,
    charge: int,
    multiplicity: int,
    console: Console,
    selected_params: Optional[List[Tuple[str, float, str, str]]] = None,
    processor=None
) -> Dict[str, Any]:
    """
    Execute the full paramfit refinement workflow.

    Steps:
    1. Convert MOL2 to XYZ for CREST
    2. Run CREST conformer search
    3. Select conformers and prepare for QM
    4. Generate Gaussian single-point inputs
    5. [User runs Gaussian externally]
    6. Extract QM energies
    7. Fit K (energy offset)
    8. Select parameters to fit (using paramfit wizard)
    9. Run paramfit
    10. Output refined frcmod

    Args:
        mol_name: Molecule name
        mol2_file: Input MOL2 file
        prmtop_file: AMBER topology file
        frcmod_file: Original frcmod file
        charge: Molecular charge
        multiplicity: Spin multiplicity
        console: Rich console for output
        selected_params: List of (param_name, score, status, section) tuples
                        that the user selected for refinement
        processor: Optional processor for session recording context

    Returns:
        Dict with workflow results
    """
    results = {
        'refinement_attempted': True,
        'refinement_success': False,
        'fitted_frcmod': None,
        'n_conformers': 0,
        'message': ''
    }

    # Create working directory
    work_dir = f"{mol_name}_paramfit"
    os.makedirs(work_dir, exist_ok=True)
    original_dir = os.getcwd()
    os.chdir(work_dir)

    console.print(f"\n[cyan]Working directory: {os.getcwd()}[/cyan]")

    # Display selected parameters for refinement
    if selected_params:
        console.print(f"\n[bold]Parameters selected for refinement:[/bold]")
        for param_name, score, status, section in selected_params:
            score_str = f"{score:.1f}" if score != float('inf') else "ATTN"
            console.print(f"  • {section}: [cyan]{param_name}[/cyan] (penalty: {score_str})")

    # Check dependencies before starting
    deps = check_all_dependencies(console)
    if not all(deps.values()):
        results['message'] = "Missing dependencies for paramfit refinement"
        os.chdir(original_dir)
        return results

    try:
        # =====================================================================
        # Step 1: Convert MOL2 to XYZ
        # =====================================================================
        console.print(f"\n[bold]Step 1: Prepare input structure[/bold]")

        xyz_file = f"{mol_name}.xyz"
        mol2_path = os.path.join(original_dir, mol2_file) if not os.path.isabs(mol2_file) else mol2_file

        if not mol2_to_xyz(mol2_path, xyz_file, console):
            results['message'] = "Failed to convert MOL2 to XYZ"
            return results

        # =====================================================================
        # Step 2: Configure and run CREST
        # =====================================================================
        console.print(f"\n[bold]Step 2: CREST conformer search[/bold]")

        # Check for existing CREST output
        existing_ensemble = "crest_conformers.xyz"
        ensemble_file = None

        if os.path.exists(existing_ensemble):
            # Count conformers in existing file
            try:
                with open(existing_ensemble, 'r') as f:
                    first_line = f.readline().strip()
                    n_atoms = int(first_line)
                # Count structures (each structure is n_atoms + 2 lines)
                with open(existing_ensemble, 'r') as f:
                    content = f.read()
                    lines = content.strip().split('\n')
                    n_existing = len(lines) // (n_atoms + 2)
                console.print(f"[green]Found existing CREST output: {existing_ensemble}[/green]")
                console.print(f"[green]  Contains {n_existing} conformers[/green]")

                use_existing = confirm_with_context(
                    processor,
                    "Use existing CREST conformers?",
                    default=True,
                    module="Paramfit Refinement",
                    description="Use existing CREST conformers",
                )

                if use_existing:
                    ensemble_file = existing_ensemble
                    console.print(f"[cyan]Skipping CREST calculation, using existing conformers[/cyan]")
            except Exception as e:
                console.print(f"[yellow]Could not read existing file: {e}[/yellow]")
                console.print("[yellow]Will run new CREST calculation[/yellow]")

        # Run CREST if no existing file or user chose to re-run
        if ensemble_file is None:
            console.print("\n[cyan]CREST settings:[/cyan]")

            energy_window = prompt_with_context(
                processor,
                "Energy window for conformers (kcal/mol)",
                default="6.0",
                module="Paramfit Refinement",
                description="CREST energy window (kcal/mol)",
            )
            energy_window = float(energy_window)

            n_conformers_target = int_prompt_with_context(
                processor,
                "Maximum conformers to use for fitting",
                default=30,
                module="Paramfit Refinement",
                description="Maximum conformers for fitting",
            )

            nproc_crest = int_prompt_with_context(
                processor,
                "Number of CPU threads for CREST",
                default=4,
                module="Paramfit Refinement",
                description="CREST CPU threads",
            )

            # Solvation option - recommended for charged molecules
            use_solvation = False
            if abs(charge) >= 1:
                console.print(f"\n[yellow]Note: Molecule has charge {charge}. Implicit solvation can help xTB convergence.[/yellow]")
                use_solvation = confirm_with_context(
                    processor,
                    "Use implicit solvation (ALPB water)?",
                    default=True,
                    module="Paramfit Refinement",
                    description="Use ALPB water implicit solvation for CREST",
                )

            # GFN method selection
            console.print("\n[cyan]xTB method:[/cyan]")
            console.print("  1: GFN2-xTB (more accurate, default)")
            console.print("  2: GFN1-xTB (more robust for difficult cases)")
            gfn_choice = prompt_with_context(
                processor,
                "Method",
                choices=["1", "2"],
                default="1",
                module="Paramfit Refinement",
                description="xTB method selection",
                options_map={"1": "GFN2-xTB (more accurate)", "2": "GFN1-xTB (more robust)"},
            )
            gfn_method = "gfn2" if gfn_choice == "1" else "gfn1"

            # Run CREST
            ensemble_file = run_crest_conformer_search(
                xyz_file=xyz_file,
                charge=charge,
                multiplicity=multiplicity,
                energy_window=energy_window,
                nproc=nproc_crest,
                console=console,
                use_solvation=use_solvation,
                gfn_method=gfn_method
            )

            if not ensemble_file:
                results['message'] = "CREST conformer search failed"
                return results
        else:
            # When using existing file, still need n_conformers_target
            n_conformers_target = int_prompt_with_context(
                processor,
                "Maximum conformers to use for fitting",
                default=30,
                module="Paramfit Refinement",
                description="Maximum conformers for fitting",
            )

        # =====================================================================
        # Step 3: Parse and select conformers
        # =====================================================================
        console.print(f"\n[bold]Step 3: Select conformers[/bold]")

        all_conformers = parse_crest_ensemble(ensemble_file, console)

        if not all_conformers:
            results['message'] = "No conformers found in CREST output"
            return results

        # Select top N conformers
        conformers = all_conformers[:n_conformers_target]
        results['n_conformers'] = len(conformers)

        console.print(f"[green]Selected {len(conformers)} conformers for QM calculations[/green]")

        # Write mdcrd for paramfit (do this now to ensure order matches)
        mdcrd_file = f"{mol_name}_conformers.mdcrd"
        if not write_amber_mdcrd(conformers, mdcrd_file, console):
            results['message'] = "Failed to write mdcrd file"
            return results

        # =====================================================================
        # Step 4: Generate Gaussian inputs
        # =====================================================================
        console.print(f"\n[bold]Step 4: Configure Gaussian single-point calculations[/bold]")

        console.print(Panel(
            "[bold cyan]Gaussian Route Line Keywords[/bold cyan]\n\n"
            "Enter any Gaussian route line keywords for the single-point calculations.\n"
            "This appears after '#' in the input file (e.g., '# B3LYP/6-31G* SP').\n\n"
            "[bold]Examples:[/bold]\n"
            "  • B3LYP/6-31G*              (simple DFT)\n"
            "  • wB97XD/def2-TZVP          (dispersion-corrected DFT)\n"
            "  • MP2/aug-cc-pVDZ           (post-HF method)\n"
            "  • B3LYP/6-31G* EmpiricalDispersion=GD3BJ   (with dispersion)\n\n"
            "[grey50]The 'SP' keyword for single-point is added automatically.[/grey50]",
            title="QM Method Settings",
            border_style="cyan",
            expand=False
        ))

        console.print("\n[cyan]Gaussian settings:[/cyan]")

        method = prompt_with_context(
            processor,
            "Route line keywords (method/basis + options)",
            default="B3LYP/6-31G*",
            module="Paramfit Refinement",
            description="Gaussian route line keywords",
        )

        memory = prompt_with_context(
            processor,
            "Memory per job",
            default="4GB",
            module="Paramfit Refinement",
            description="Gaussian memory allocation",
        )

        nproc_gaussian = int_prompt_with_context(
            processor,
            "Processors per job",
            default=4,
            module="Paramfit Refinement",
            description="Gaussian processors per job",
        )

        # Ask about parallel batch execution
        parallel_jobs = int_prompt_with_context(
            processor,
            "Simultaneous jobs in batch script",
            default=2,
            module="Paramfit Refinement",
            description="Simultaneous Gaussian jobs in batch script",
        )

        gaussian_dir = "gaussian_inputs"
        input_files = write_gaussian_sp_inputs(
            conformers=conformers,
            output_dir=gaussian_dir,
            mol_name=mol_name,
            charge=charge,
            multiplicity=multiplicity,
            method=method,
            memory=memory,
            nproc=nproc_gaussian,
            console=console
        )

        # Generate batch run script
        batch_script = write_gaussian_batch_script(
            output_dir=gaussian_dir,
            mol_name=mol_name,
            n_conformers=len(conformers),
            parallel_jobs=parallel_jobs,
            console=console
        )

        # =====================================================================
        # Pause point: User runs Gaussian
        # =====================================================================
        console.print(Panel(
            "[bold yellow]Gaussian Calculations Required[/bold yellow]\n\n"
            f"[bold]{len(input_files)}[/bold] Gaussian input files have been created in:\n"
            f"  [cyan]{os.path.join(os.getcwd(), gaussian_dir)}/[/cyan]\n\n"
            "[bold]Run script (recommended):[/bold]\n"
            f"  [green]cd {os.path.join(os.getcwd(), gaussian_dir)} && bash run_gaussian.sh[/green]\n\n"
            f"  This runs {parallel_jobs} jobs simultaneously and tracks progress.\n"
            f"  Edit PARALLEL_JOBS in the script to change concurrency.\n\n"
            "[bold]Alternative (sequential):[/bold]\n"
            f"  cd {os.path.join(os.getcwd(), gaussian_dir)}\n"
            f"  for f in *.gjf; do g16 $f > ${{f%.gjf}}.log; done\n\n"
            "[grey50]The .log files must have the same naming pattern:[/grey50]\n"
            f"[grey50]  {mol_name}_conf_000.log, {mol_name}_conf_001.log, etc.[/grey50]",
            title="Action Required",
            border_style="yellow",
            expand=False
        ))

        console.print(f"\n[grey50]To resume later: Re-run paramfit refinement for {mol_name}[/grey50]")
        console.print(f"[grey50]The workflow will detect completed Gaussian calculations.[/grey50]")

        proceed = confirm_with_context(
            processor,
            "\nHave you completed all Gaussian calculations?",
            default=False,
            module="Paramfit Refinement",
            description="Confirm Gaussian calculations completed",
        )

        if not proceed:
            results['message'] = "Waiting for Gaussian calculations"
            console.print(f"\n[cyan]Files saved in: {os.getcwd()}[/cyan]")
            console.print("[cyan]Re-run this workflow after completing Gaussian calculations.[/cyan]")
            return results

        # =====================================================================
        # Step 5: Extract Gaussian energies
        # =====================================================================
        console.print(f"\n[bold]Step 5: Extract QM energies[/bold]")

        energies_file = f"{mol_name}_energies.dat"
        energies = extract_gaussian_energies(
            log_dir=gaussian_dir,
            mol_name=mol_name,
            n_conformers=len(conformers),
            output_file=energies_file,
            console=console
        )

        if energies is None:
            results['message'] = "Failed to extract Gaussian energies"
            return results

        # =====================================================================
        # Step 6: Fit K
        # =====================================================================
        console.print(f"\n[bold]Step 6: Fit energy offset constant (K)[/bold]")

        # Use absolute path for prmtop
        prmtop_path = os.path.join(original_dir, prmtop_file) if not os.path.isabs(prmtop_file) else prmtop_file

        k_value = run_paramfit_k_fitting(
            prmtop_file=prmtop_path,
            mdcrd_file=mdcrd_file,
            energies_file=energies_file,
            n_structures=len(conformers),
            console=console
        )

        if k_value is None:
            results['message'] = "Failed to fit K value"
            return results

        # =====================================================================
        # Step 7: Select parameters to fit
        # =====================================================================
        console.print(f"\n[bold]Step 7: Select parameters to fit[/bold]")

        params_file = f"{mol_name}_params.in"

        success = False

        # Try automated selection if we have pre-selected parameters
        if selected_params:
            console.print("\n[cyan]Using your earlier parameter selections from the penalty table...[/cyan]")

            # Ask user what properties to fit
            console.print("\n[bold]What properties should be fitted?[/bold]")
            console.print("  1) Force constants only (KR, KT, KP) - [green]Recommended[/green]")
            console.print("     Fewer parameters = more stable fitting")
            console.print("  2) All properties (force constants + equilibrium values)")
            console.print("     More parameters = risk of overfitting")

            fit_choice = prompt_with_context(
                processor,
                "Fit mode",
                choices=["1", "2"],
                default="1",
                module="Paramfit Refinement",
                description="Fit mode selection",
                options_map={
                    "1": "Force constants only (KR, KT, KP)",
                    "2": "All properties (force constants + equilibrium values)",
                },
            )

            force_constants_only = (fit_choice == "1")

            success = run_paramfit_set_params_automated(
                prmtop_file=prmtop_path,
                params_file=params_file,
                selected_params=selected_params,
                console=console,
                force_constants_only=force_constants_only
            )

            if not success:
                console.print("\n[yellow]Automated parameter selection failed.[/yellow]")
                console.print("[yellow]Falling back to interactive mode...[/yellow]\n")

        # Fall back to interactive wizard if no pre-selection or automation failed
        if not success:
            console.print("\n[cyan]You will now use paramfit's interactive wizard to select which[/cyan]")
            console.print("[cyan]parameters to optimize. This gives you full control over the fitting.[/cyan]")

            success = run_paramfit_set_params_wizard(prmtop_path, params_file, console)

        if not success:
            results['message'] = "Parameter selection cancelled"
            return results

        # =====================================================================
        # Step 8: Run parameter fitting
        # =====================================================================
        console.print(f"\n[bold]Step 8: Fit parameters[/bold]")

        output_frcmod = f"{mol_name}_fitted.frcmod"

        success = run_paramfit_parameter_fitting(
            prmtop_file=prmtop_path,
            mdcrd_file=mdcrd_file,
            energies_file=energies_file,
            params_file=params_file,
            n_structures=len(conformers),
            k_value=k_value,
            output_frcmod=output_frcmod,
            console=console
        )

        if not success:
            results['message'] = "Parameter fitting failed"
            return results

        # =====================================================================
        # Step 9: Summary
        # =====================================================================
        results['refinement_success'] = True
        results['fitted_frcmod'] = os.path.abspath(output_frcmod)
        results['message'] = "Parameter refinement completed successfully"

        console.print(Panel(
            f"[bold green]Parameter Refinement Complete![/bold green]\n\n"
            f"[bold]Results:[/bold]\n"
            f"  Conformers used: {len(conformers)}\n"
            f"  K value: {k_value:.2f} kcal/mol\n"
            f"  Fitted frcmod: {output_frcmod}\n\n"
            f"[bold]Next steps:[/bold]\n"
            f"  1. Review the fitted parameters in {output_frcmod}\n"
            f"  2. Compare with original parmchk2 parameters\n"
            f"  3. Use the fitted frcmod in your tleap workflow\n\n"
            f"[grey50]Files are in: {os.getcwd()}[/grey50]",
            title="Refinement Complete",
            border_style="green",
            expand=False
        ))

        return results

    except Exception as e:
        console.print(f"[red]Error in paramfit workflow: {str(e)}[/red]")
        results['message'] = f"Error: {str(e)}"
        return results

    finally:
        os.chdir(original_dir)
