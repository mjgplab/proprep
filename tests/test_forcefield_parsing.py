#!/usr/bin/env python3
"""
Test script for force field file parsing.

Verifies that ForcefieldDataCollector correctly:
1. Parses leaprc files and follows source directives
2. Loads lib files and extracts atom definitions
3. Loads frcmod files and extracts bond/angle parameters
4. Handles residue mapping (HEM → HEH)
5. Provides correct lookup chains

Run with: python tests/test_forcefield_parsing.py
"""

import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import pytest

from rich.console import Console
from rich.table import Table
from rich.panel import Panel

console = Console()

# These parse real leaprc/lib/frcmod files out of an Amber installation, so
# they cannot run without one. Skip rather than fail: an absent AMBERHOME says
# nothing about whether the parsing code is correct, and failing on it hides
# real regressions in the rest of the suite behind noise.
pytestmark = pytest.mark.skipif(
    not os.environ.get("AMBERHOME"),
    reason="requires an Amber installation (AMBERHOME not set)",
)


def check_amberhome():
    """Verify AMBERHOME is set."""
    amberhome = os.environ.get('AMBERHOME')
    if not amberhome:
        console.print("[red]ERROR: AMBERHOME environment variable not set[/red]")
        console.print("Please set AMBERHOME to your AMBER installation directory")
        return None

    amberhome = Path(amberhome)
    if not amberhome.exists():
        console.print(f"[red]ERROR: AMBERHOME directory does not exist: {amberhome}[/red]")
        return None

    console.print(f"[green]AMBERHOME: {amberhome}[/green]")
    return amberhome


def test_leaprc_parsing():
    """Test 1: Verify leaprc file parsing."""
    console.print(Panel("[bold]Test 1: leaprc Parsing[/bold]"))

    from proprep.forcefield_prep.forcefield_data_collector import ForcefieldDataCollector

    # Test with ff14SB
    console.print("\n[cyan]Testing leaprc.protein.ff14SB...[/cyan]")
    ff = ForcefieldDataCollector(console=console)

    amberhome = ff.amberhome
    leaprc_path = amberhome / "dat/leap/cmd/leaprc.protein.ff14SB"

    if not leaprc_path.exists():
        console.print(f"[yellow]Skipping: {leaprc_path} not found[/yellow]")
        return False

    ff._collect_from_leaprc(leaprc_path)

    # Check results
    n_atoms = len(ff.atom_definitions)
    n_bonds = len(ff.bond_parameters) // 2  # Stored both directions
    n_angles = len(ff.angle_parameters) // 2

    console.print(f"  Atom definitions loaded: {n_atoms}")
    console.print(f"  Bond parameters loaded: {n_bonds}")
    console.print(f"  Angle parameters loaded: {n_angles}")
    console.print(f"  Leaprcs processed: {len(ff.selected_leaprcs)}")
    console.print(f"  Frcmod files: {len(ff.selected_frcmod_files)}")

    # Verify some standard residues
    # Note: ff14SB uses HID/HIE/HIP instead of HIS
    test_residues = ['ALA', 'GLY', 'HID', 'HIE', 'CYS', 'ACE', 'NME']
    missing = []
    for res in test_residues:
        if not ff.has_residue(res):
            missing.append(res)

    if missing:
        console.print(f"[red]  Missing expected residues: {missing}[/red]")
        return False
    else:
        console.print(f"[green]  All standard residues found: {test_residues}[/green]")
        console.print(f"[dim]  (Note: HIS is represented as HID/HIE/HIP in AMBER)[/dim]")

    # Verify some atom definitions
    ala_ca = ff.get_atom_type('ALA', 'CA')
    if ala_ca:
        console.print(f"  ALA-CA: type={ala_ca.atom_type}, charge={ala_ca.atom_charge:.4f}")
    else:
        console.print("[red]  ALA-CA not found![/red]")
        return False

    return True


def test_conste_parsing():
    """Test 2: Verify leaprc.conste parsing (special case with heme)."""
    console.print(Panel("[bold]Test 2: leaprc.conste Parsing (Heme)[/bold]"))

    from proprep.forcefield_prep.forcefield_data_collector import ForcefieldDataCollector

    ff = ForcefieldDataCollector(console=console)
    amberhome = ff.amberhome

    leaprc_path = amberhome / "dat/leap/cmd/leaprc.conste"

    if not leaprc_path.exists():
        console.print(f"[yellow]Skipping: {leaprc_path} not found[/yellow]")
        console.print("[yellow]leaprc.conste may not be installed in your AMBER version[/yellow]")
        return None  # Skip, not failure

    console.print(f"\n[cyan]Testing {leaprc_path.name}...[/cyan]")
    ff._collect_from_leaprc(leaprc_path)

    # Check for heme residues (HEH = reduced, HEO = oxidized)
    heme_residues = ['HEH', 'HEO', 'PRN', 'HIO', 'CYO']
    found = []
    missing = []

    for res in heme_residues:
        if ff.has_residue(res):
            found.append(res)
        else:
            missing.append(res)

    console.print(f"  Found heme residues: {found}")
    if missing:
        console.print(f"[yellow]  Missing heme residues: {missing}[/yellow]")

    # Check HEH-FE if available
    if 'HEH' in found:
        fe_atom = ff.get_atom_type('HEH', 'FE')
        if fe_atom:
            console.print(f"  HEH-FE: type={fe_atom.atom_type}, charge={fe_atom.atom_charge:.4f}")
        else:
            console.print("[yellow]  HEH-FE not found (atom name may differ)[/yellow]")
            # Try alternate names
            for name in ['FE', 'FE2', 'FE1']:
                fe_atom = ff.get_atom_type('HEH', name)
                if fe_atom:
                    console.print(f"  Found as HEH-{name}: type={fe_atom.atom_type}, charge={fe_atom.atom_charge:.4f}")
                    break

    # Check frcmod parameters for Fe-N bonds
    console.print(f"\n  Bond parameters loaded: {len(ff.bond_parameters) // 2}")
    console.print(f"  Angle parameters loaded: {len(ff.angle_parameters) // 2}")

    # Look for Fe-N type bonds
    fe_bonds = [k for k in ff.bond_parameters.keys() if 'FE' in k.upper() or 'M1' in k.upper()]
    if fe_bonds:
        console.print(f"  Fe-related bond types: {fe_bonds[:10]}...")  # Show first 10
    else:
        console.print("[yellow]  No Fe-related bond types found[/yellow]")

    return len(found) > 0


def test_frcmod_parsing():
    """Test 3: Verify frcmod file parsing."""
    console.print(Panel("[bold]Test 3: frcmod File Parsing[/bold]"))

    from proprep.forcefield_prep.forcefield_data_collector import ForcefieldDataCollector

    ff = ForcefieldDataCollector(console=console)
    amberhome = ff.amberhome

    # Find an frcmod file to test
    frcmod_candidates = [
        amberhome / "dat/leap/parm/frcmod.ff14SB",
        amberhome / "dat/leap/parm/frcmod.ff19SB",
        amberhome / "dat/leap/parm/frcmod.conste",
    ]

    frcmod_path = None
    for candidate in frcmod_candidates:
        if candidate.exists():
            frcmod_path = candidate
            break

    if not frcmod_path:
        console.print("[yellow]No frcmod files found to test[/yellow]")
        return None

    console.print(f"\n[cyan]Testing {frcmod_path.name}...[/cyan]")

    # Clear any existing parameters
    ff.bond_parameters = {}
    ff.angle_parameters = {}

    ff._load_frcmod_file(frcmod_path, source=frcmod_path.name)

    n_bonds = len(ff.bond_parameters) // 2
    n_angles = len(ff.angle_parameters) // 2

    console.print(f"  Bond parameters: {n_bonds}")
    console.print(f"  Angle parameters: {n_angles}")

    # Show some examples
    if ff.bond_parameters:
        console.print("\n  Sample bond parameters:")
        for i, (key, param) in enumerate(ff.bond_parameters.items()):
            if i >= 5:
                break
            console.print(f"    {key}: k={param.force_constant:.1f}, eq={param.eq_length:.3f}")

    if ff.angle_parameters:
        console.print("\n  Sample angle parameters:")
        for i, (key, param) in enumerate(ff.angle_parameters.items()):
            if i >= 5:
                break
            console.print(f"    {key}: k={param.force_constant:.2f}, eq={param.eq_angle:.1f}")

    return n_bonds > 0 or n_angles > 0


def test_residue_mapping():
    """Test 4: Verify residue mapping (HEM → HEH)."""
    console.print(Panel("[bold]Test 4: Residue Mapping[/bold]"))

    from proprep.forcefield_prep.forcefield_data_collector import ForcefieldDataCollector

    ff = ForcefieldDataCollector(console=console)
    amberhome = ff.amberhome

    # First load conste to get HEH
    leaprc_path = amberhome / "dat/leap/cmd/leaprc.conste"
    if not leaprc_path.exists():
        console.print("[yellow]Skipping: leaprc.conste not available[/yellow]")
        return None

    ff._collect_from_leaprc(leaprc_path)

    if not ff.has_residue('HEH'):
        console.print("[yellow]HEH residue not found, cannot test mapping[/yellow]")
        return None

    # Get HEH atoms before mapping
    heh_atoms = ff.get_residue_atoms('HEH')
    console.print(f"  HEH atoms before mapping: {len(heh_atoms)}")

    # Verify HEM is NOT present
    if ff.has_residue('HEM'):
        console.print("[yellow]  HEM already exists (unexpected)[/yellow]")
    else:
        console.print("  HEM not present (as expected)")

    # Simulate mapping HEM → HEH
    console.print("\n[cyan]  Mapping HEM → HEH...[/cyan]")
    ff.residue_aliases['HEM'] = 'HEH'
    ff._copy_residue_definitions('HEH', 'HEM')

    # Verify HEM now exists
    if ff.has_residue('HEM'):
        console.print("[green]  HEM now recognized via mapping[/green]")
    else:
        console.print("[red]  HEM still not found after mapping![/red]")
        return False

    # Get HEM atoms after mapping
    hem_atoms = ff.get_residue_atoms('HEM')
    console.print(f"  HEM atoms after mapping: {len(hem_atoms)}")

    if len(hem_atoms) == len(heh_atoms):
        console.print(f"[green]  Atom count matches: {len(hem_atoms)}[/green]")
    else:
        console.print(f"[red]  Atom count mismatch! HEH={len(heh_atoms)}, HEM={len(hem_atoms)}[/red]")
        return False

    # Verify charge lookup works for HEM
    hem_fe = ff.get_atom_type('HEM', 'FE')
    heh_fe = ff.get_atom_type('HEH', 'FE')

    if hem_fe and heh_fe:
        if hem_fe.atom_charge == heh_fe.atom_charge:
            console.print(f"[green]  Charges match: HEM-FE={hem_fe.atom_charge:.4f}, HEH-FE={heh_fe.atom_charge:.4f}[/green]")
            return True
        else:
            console.print(f"[red]  Charge mismatch: HEM-FE={hem_fe.atom_charge:.4f}, HEH-FE={heh_fe.atom_charge:.4f}[/red]")
            return False
    elif hem_fe is None and heh_fe is None:
        console.print("[yellow]  FE atom not found in either (name may differ)[/yellow]")
        # Still a pass if mapping worked
        return True
    else:
        console.print("[red]  Inconsistent FE lookup[/red]")
        return False


def test_bond_parameter_lookup():
    """Test 5: Verify bond parameter lookup in ForcefieldDataCollector.

    NOTE: ForcefieldDataCollector only loads frcmod files (modifications).
    Base parameters (like C-N) are in parmXX.dat and require FFParameterReader.
    """
    console.print(Panel("[bold]Test 5: Bond Parameter Lookup (frcmod)[/bold]"))

    from proprep.forcefield_prep.forcefield_data_collector import ForcefieldDataCollector

    ff = ForcefieldDataCollector(console=console)
    amberhome = ff.amberhome

    # Load ff14SB
    leaprc_path = amberhome / "dat/leap/cmd/leaprc.protein.ff14SB"
    if leaprc_path.exists():
        ff._collect_from_leaprc(leaprc_path)

    console.print("\n[cyan]Testing bond lookups in ForcefieldDataCollector...[/cyan]")
    console.print("[dim]Note: Only frcmod parameters are here, base params need FFParameterReader[/dim]\n")

    # Test frcmod-specific bonds (from frcmod.ff14SB)
    test_bonds = [
        ('C', '2C'),   # ff14SB specific (new 2C type)
        ('C8', 'C8'),  # ff14SB specific (new C8 type)
        ('C*', '2C'),  # ff14SB specific
        ('C', 'N'),    # Base param - should NOT be found here (in parm10.dat)
    ]

    results = []
    for type1, type2 in test_bonds:
        param = ff.get_bond_parameter(type1, type2)
        if param:
            results.append((f"{type1}-{type2}", f"k={param.force_constant:.1f}, eq={param.eq_length:.3f}", "Found", True))
        else:
            # Check if this is expected to be missing
            expected_missing = (type1, type2) == ('C', 'N')
            status = "Not found (expected)" if expected_missing else "Not found"
            results.append((f"{type1}-{type2}", "-", status, expected_missing))

    table = Table(title="Bond Parameter Lookup (frcmod only)")
    table.add_column("Bond", style="cyan")
    table.add_column("Parameters", style="white")
    table.add_column("Status", style="green")

    for bond, params, status, _ in results:
        if "expected" in status.lower():
            style = "dim"
        elif status == "Found":
            style = "green"
        else:
            style = "yellow"
        table.add_row(bond, params, f"[{style}]{status}[/{style}]")

    console.print(table)

    # C-2C should be found (it's in frcmod.ff14SB)
    c_2c_param = ff.get_bond_parameter('C', '2C')
    return c_2c_param is not None


def test_angle_parameter_lookup():
    """Test 6: Verify angle parameter lookup in ForcefieldDataCollector.

    NOTE: ForcefieldDataCollector only loads frcmod files (modifications).
    Base parameters are in parmXX.dat and require FFParameterReader.
    """
    console.print(Panel("[bold]Test 6: Angle Parameter Lookup (frcmod)[/bold]"))

    from proprep.forcefield_prep.forcefield_data_collector import ForcefieldDataCollector

    ff = ForcefieldDataCollector(console=console)
    amberhome = ff.amberhome

    # Load ff14SB
    leaprc_path = amberhome / "dat/leap/cmd/leaprc.protein.ff14SB"
    if leaprc_path.exists():
        ff._collect_from_leaprc(leaprc_path)

    console.print("\n[cyan]Testing angle lookups in ForcefieldDataCollector...[/cyan]")
    console.print("[dim]Note: Only frcmod parameters are here, base params need FFParameterReader[/dim]\n")

    # Test frcmod-specific angles (from frcmod.ff14SB)
    test_angles = [
        ('N', 'C', '2C'),   # ff14SB specific (new 2C type)
        ('O', 'C', '2C'),   # ff14SB specific
        ('C8', 'C8', 'C8'), # ff14SB specific (new C8 type)
        ('N', 'CA', 'C'),   # Base param - should NOT be found here
    ]

    results = []
    for type1, type2, type3 in test_angles:
        param = ff.get_angle_parameter(type1, type2, type3)
        if param:
            results.append((f"{type1}-{type2}-{type3}", f"k={param.force_constant:.2f}, eq={param.eq_angle:.1f}", "Found", True))
        else:
            expected_missing = (type1, type2, type3) == ('N', 'CA', 'C')
            status = "Not found (expected)" if expected_missing else "Not found"
            results.append((f"{type1}-{type2}-{type3}", "-", status, expected_missing))

    table = Table(title="Angle Parameter Lookup (frcmod only)")
    table.add_column("Angle", style="cyan")
    table.add_column("Parameters", style="white")
    table.add_column("Status", style="green")

    for angle, params, status, _ in results:
        if "expected" in status.lower():
            style = "dim"
        elif status == "Found":
            style = "green"
        else:
            style = "yellow"
        table.add_row(angle, params, f"[{style}]{status}[/{style}]")

    console.print(table)

    # N-C-2C should be found (it's in frcmod.ff14SB)
    n_c_2c_param = ff.get_angle_parameter('N', 'C', '2C')
    return n_c_2c_param is not None


def test_ff_library():
    """Test 7: Verify ForceFieldLibrary (legacy mol2 parser)."""
    console.print(Panel("[bold]Test 7: ForceFieldLibrary (mol2)[/bold]"))

    try:
        from proprep.forcefield_prep.mcpb.ff_library import ForceFieldLibrary
    except ImportError as e:
        console.print(f"[yellow]Could not import ForceFieldLibrary: {e}[/yellow]")
        return None

    console.print("\n[cyan]Testing ForceFieldLibrary...[/cyan]")

    try:
        ff_lib = ForceFieldLibrary("ff14SB")

        # Test some lookups
        test_atoms = ['ALA-CA', 'GLY-N', 'ACE-C', 'NME-CH3']

        results = []
        for atom_key in test_atoms:
            info = ff_lib.get_atom_type(atom_key)
            if info:
                results.append((atom_key, f"type={info.get('type', '?')}, charge={info.get('charge', 0):.4f}", "Found"))
            else:
                results.append((atom_key, "-", "Not found"))

        table = Table(title="ForceFieldLibrary Lookup")
        table.add_column("Atom", style="cyan")
        table.add_column("Info", style="white")
        table.add_column("Status", style="green")

        for atom, info, status in results:
            style = "green" if status == "Found" else "yellow"
            table.add_row(atom, info, f"[{style}]{status}[/{style}]")

        console.print(table)

        return any(s == "Found" for _, _, s in results)

    except Exception as e:
        console.print(f"[red]Error: {e}[/red]")
        return False


def test_ff_parameter_reader():
    """Test 8: Verify FFParameterReader (parmXX.dat parser)."""
    console.print(Panel("[bold]Test 8: FFParameterReader (parmXX.dat)[/bold]"))

    try:
        from proprep.forcefield_prep.mcpb.ff_parameter_reader import FFParameterReader
    except ImportError as e:
        console.print(f"[yellow]Could not import FFParameterReader: {e}[/yellow]")
        return None

    amberhome = os.environ.get('AMBERHOME')
    if not amberhome:
        console.print("[yellow]AMBERHOME not set[/yellow]")
        return None

    console.print("\n[cyan]Testing FFParameterReader...[/cyan]")

    try:
        ff_reader = FFParameterReader("ff14SB", amberhome)

        stats = ff_reader.get_statistics()
        console.print(f"  Bonds loaded: {stats.get('n_bonds', 0)}")
        console.print(f"  Angles loaded: {stats.get('n_angles', 0)}")
        console.print(f"  Dihedrals loaded: {stats.get('n_dihedrals', 0)}")

        # Test some lookups
        cn_bond = ff_reader.get_bond_parameter('C', 'N')
        if cn_bond:
            console.print(f"  C-N bond: k={cn_bond.force_constant:.1f}, eq={cn_bond.eq_length:.3f}")

        nca_c_angle = ff_reader.get_angle_parameter('N', 'CA', 'C')
        if nca_c_angle:
            console.print(f"  N-CA-C angle: k={nca_c_angle.force_constant:.2f}, eq={nca_c_angle.eq_angle:.1f}")

        return stats.get('n_bonds', 0) > 0

    except Exception as e:
        console.print(f"[red]Error: {e}[/red]")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all tests."""
    console.print(Panel("[bold magenta]Force Field Parsing Tests[/bold magenta]",
                       subtitle="Verifying AMBER file handling"))

    # Check prerequisites
    if not check_amberhome():
        return 1

    console.print()

    # Run tests
    tests = [
        ("leaprc Parsing", test_leaprc_parsing),
        ("leaprc.conste Parsing", test_conste_parsing),
        ("frcmod Parsing", test_frcmod_parsing),
        ("Residue Mapping", test_residue_mapping),
        ("Bond Parameter Lookup", test_bond_parameter_lookup),
        ("Angle Parameter Lookup", test_angle_parameter_lookup),
        ("ForceFieldLibrary", test_ff_library),
        ("FFParameterReader", test_ff_parameter_reader),
    ]

    results = []
    for name, test_func in tests:
        try:
            result = test_func()
            results.append((name, result))
        except Exception as e:
            console.print(f"[red]Test '{name}' raised exception: {e}[/red]")
            import traceback
            traceback.print_exc()
            results.append((name, False))
        console.print()

    # Summary
    console.print(Panel("[bold]Test Summary[/bold]"))

    table = Table()
    table.add_column("Test", style="cyan")
    table.add_column("Result", style="white")

    passed = 0
    failed = 0
    skipped = 0

    for name, result in results:
        if result is True:
            table.add_row(name, "[green]PASSED[/green]")
            passed += 1
        elif result is False:
            table.add_row(name, "[red]FAILED[/red]")
            failed += 1
        else:
            table.add_row(name, "[yellow]SKIPPED[/yellow]")
            skipped += 1

    console.print(table)
    console.print(f"\n[bold]Passed: {passed}, Failed: {failed}, Skipped: {skipped}[/bold]")

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
