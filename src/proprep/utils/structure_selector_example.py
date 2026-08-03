#!/usr/bin/env python3
"""
Example demonstrating the StructureSelector module.

Run this script to see how structure selection works.
"""

from proprep.utils.structure_selector import (
    StructureSelector,
    register_structure_type,
    get_priority_pdb_file,
    get_interactive_pdb_file,
    StructureRegistry
)
from rich.console import Console


def example_1_basic_priority():
    """Example 1: Basic priority-based selection"""
    console = Console()
    console.print("\n[bold cyan]Example 1: Basic Priority Selection[/bold cyan]\n")

    # Mock workspace
    class MockWorkspace:
        def __init__(self):
            self.data = {
                "filtered_pdb_file": "/tmp/protein_filtered.pdb",
                "original_pdb_file": "/tmp/protein_original.pdb",
                "repaired_pdb_file": "/tmp/protein_repaired.pdb"
            }

        def get(self, key, default=None):
            return self.data.get(key, default)

    workspace = MockWorkspace()

    # Create some fake files
    import os
    import tempfile
    os.makedirs("/tmp", exist_ok=True)
    for filepath in workspace.data.values():
        with open(filepath, 'w') as f:
            f.write("HEADER TEST PDB\n")

    # Use simple helper function
    selected = get_priority_pdb_file(workspace, console)
    console.print(f"\n[bold]Result:[/bold] {selected}\n")

    # Cleanup
    for filepath in workspace.data.values():
        if os.path.exists(filepath):
            os.remove(filepath)


def example_2_display_available():
    """Example 2: Display all available structures"""
    console = Console()
    console.print("\n[bold cyan]Example 2: Display Available Structures[/bold cyan]\n")

    # Mock workspace with more structures
    class MockWorkspace:
        def __init__(self):
            self.data = {
                "structure_with_prot_resnames": "/data/prot_updated.pdb",
                "transformed_pdb_file": "/data/transformed.pdb",
                "repaired_pdb_file": "/data/repaired.pdb",
                "filtered_pdb_file": "/data/filtered.pdb",
                "original_pdb_file": "/data/original.pdb"
            }

        def get(self, key, default=None):
            return self.data.get(key, default)

    workspace = MockWorkspace()

    # Create files
    import os
    os.makedirs("/data", exist_ok=True)
    for filepath in workspace.data.values():
        with open(filepath, 'w') as f:
            f.write("HEADER TEST PDB\n" * 100)  # Make files different sizes

    selector = StructureSelector(workspace, console)
    selector.display_available_structures()

    # Cleanup
    for filepath in workspace.data.values():
        if os.path.exists(filepath):
            os.remove(filepath)
    if os.path.exists("/data"):
        os.rmdir("/data")


def example_3_custom_structure_type():
    """Example 3: Register and use custom structure type"""
    console = Console()
    console.print("\n[bold cyan]Example 3: Custom Structure Type[/bold cyan]\n")

    # Register a new structure type for QM/MM optimization
    console.print("[yellow]Registering custom structure type...[/yellow]")
    register_structure_type(
        workspace_key="qmmm_optimized_structure",
        display_name="QM/MM Optimized",
        priority=2,  # Higher priority than transformed (3)
        description="After QM/MM geometry optimization"
    )

    # Show it's registered
    console.print("\n[green]All registered structure types:[/green]")
    for st in StructureRegistry.get_all():
        console.print(f"  {st.priority}. {st.display_name} ({st.workspace_key})")

    # Mock workspace with our custom structure
    class MockWorkspace:
        def __init__(self):
            self.data = {
                "qmmm_optimized_structure": "/results/qmmm_opt.pdb",
                "repaired_pdb_file": "/results/repaired.pdb",
                "original_pdb_file": "/results/original.pdb"
            }

        def get(self, key, default=None):
            return self.data.get(key, default)

    workspace = MockWorkspace()

    # Create files
    import os
    os.makedirs("/results", exist_ok=True)
    for filepath in workspace.data.values():
        with open(filepath, 'w') as f:
            f.write("HEADER TEST PDB\n")

    # Selection should pick our custom type (priority 2) over repaired (priority 3)
    selected = get_priority_pdb_file(workspace, console)
    console.print(f"\n[bold]Selected:[/bold] {selected}")
    console.print("[grey50]Note: QM/MM Optimized (priority 2) was chosen over repaired (priority 3)[/grey50]\n")

    # Cleanup
    for filepath in workspace.data.values():
        if os.path.exists(filepath):
            os.remove(filepath)
    if os.path.exists("/results"):
        os.rmdir("/results")

    # Reset registry for other examples
    StructureRegistry.clear()


def example_4_filtering():
    """Example 4: Filtering and custom selection"""
    console = Console()
    console.print("\n[bold cyan]Example 4: Filtering Options[/bold cyan]\n")

    class MockWorkspace:
        def __init__(self):
            self.data = {
                "structure_with_prot_resnames": "/work/protonated.pdb",
                "transformed_pdb_file": "/work/transformed.pdb",
                "repaired_pdb_file": "/work/repaired.pdb",
                "filtered_pdb_file": "/work/filtered.pdb",
                "original_pdb_file": "/work/original.pdb",
                "pdb_file": "/work/legacy.pdb"  # Legacy key
            }

        def get(self, key, default=None):
            return self.data.get(key, default)

    workspace = MockWorkspace()

    # Create files
    import os
    os.makedirs("/work", exist_ok=True)
    for filepath in workspace.data.values():
        with open(filepath, 'w') as f:
            f.write("HEADER TEST PDB\n")

    selector = StructureSelector(workspace, console)

    # Test 1: Exclude legacy key
    console.print("[yellow]Test 1: Exclude legacy key[/yellow]")
    result = selector.get_structure(include_legacy=False, silent=False)

    # Test 2: Only consider structures from priority 3-5
    console.print("\n[yellow]Test 2: Only priority 3-5 (repaired, filtered, original)[/yellow]")
    result = selector.get_structure(min_priority=3, max_priority=5, silent=False)

    # Test 3: Exclude specific keys
    console.print("\n[yellow]Test 3: Exclude protonated and transformed[/yellow]")
    result = selector.get_structure(
        exclude_keys=["structure_with_prot_resnames", "transformed_pdb_file"],
        silent=False
    )

    console.print()

    # Cleanup
    for filepath in workspace.data.values():
        if os.path.exists(filepath):
            os.remove(filepath)
    if os.path.exists("/work"):
        os.rmdir("/work")


def example_5_interactive_mode():
    """Example 5: Interactive selection (commented out - requires user input)"""
    console = Console()
    console.print("\n[bold cyan]Example 5: Interactive Selection[/bold cyan]\n")

    console.print("[yellow]This example is commented out because it requires user input.[/yellow]")
    console.print("[grey50]To try it, uncomment the code in this function.[/grey50]\n")

    # Uncomment to test interactively:
    """
    class MockWorkspace:
        def __init__(self):
            self.data = {
                "structure_with_prot_resnames": "/interactive/protonated.pdb",
                "repaired_pdb_file": "/interactive/repaired.pdb",
                "filtered_pdb_file": "/interactive/filtered.pdb"
            }
        def get(self, key, default=None):
            return self.data.get(key, default)

    workspace = MockWorkspace()

    # Create files
    import os
    os.makedirs("/interactive", exist_ok=True)
    for filepath in workspace.data.values():
        with open(filepath, 'w') as f:
            f.write("HEADER TEST PDB\n" * 50)

    # This will prompt user to select
    selected = get_interactive_pdb_file(workspace, console)
    console.print(f"\n[bold]You selected:[/bold] {selected}\n")

    # Cleanup
    for filepath in workspace.data.values():
        if os.path.exists(filepath):
            os.remove(filepath)
    if os.path.exists("/interactive"):
        os.rmdir("/interactive")
    """


def run_all_examples():
    """Run all examples"""
    console = Console()
    console.print("\n" + "="*70)
    console.print("[bold green]Structure Selector Examples[/bold green]")
    console.print("="*70)

    example_1_basic_priority()
    input("\nPress Enter to continue to Example 2...")

    example_2_display_available()
    input("\nPress Enter to continue to Example 3...")

    example_3_custom_structure_type()
    input("\nPress Enter to continue to Example 4...")

    example_4_filtering()
    input("\nPress Enter to continue to Example 5...")

    example_5_interactive_mode()

    console.print("\n" + "="*70)
    console.print("[bold green]Examples Complete![/bold green]")
    console.print("="*70 + "\n")


if __name__ == "__main__":
    run_all_examples()
