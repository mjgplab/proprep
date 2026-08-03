"""
Prototype: Step 0 with New Layout System
This demonstrates the simplified, cleaner UX for structure file selection.
"""

from pathlib import Path
from typing import List, Dict, Optional 
from rich.console import Console

from proprep.utils.prompts import prompt_with_context
from layout_helpers import WorkflowLayoutFormatter, format_file_path


class Step0Prototype:
    """Prototype of refactored Step 0: Structure File Selection."""

    def __init__(self, console: Optional[Console] = None, processor=None):
        self.console = console or Console()
        self.processor = processor
        self.layout = WorkflowLayoutFormatter(self.console)

    def execute_step0(self, structure_pairs: List[Dict] = None) -> Optional[List[Dict]]:
        """
        Execute Step 0: Structure File Selection.

        Args:
            structure_pairs: Existing structure pairs (for demo)

        Returns:
            List of structure pairs, or None if user exits
        """
        if structure_pairs is None:
            structure_pairs = []

        # For demo purposes, simulate detected pairs
        suggested_pairs = self._simulate_detection()

        while True:
            # Build status lines
            status_lines = self._build_status_lines(structure_pairs, suggested_pairs)

            # Build actions menu
            actions, option_map = self._build_actions_menu(structure_pairs, suggested_pairs)

            # Display using new layout system
            self.layout.simple_prompt(
                step_num=0,
                status_lines=status_lines,
                actions=actions
            )

            # Get user choice
            valid_choices = list(option_map.keys())
            choice = prompt_with_context(
                self.processor,
                "Enter choice",
                choices=valid_choices,
                module="MD Manager - Step 0",
                description="Structure file selection action",
            )

            # Handle choice
            result = self._handle_choice(choice, option_map, structure_pairs, suggested_pairs)

            if result == "exit":
                return None
            elif result == "next":
                return structure_pairs
            elif result == "continue":
                continue  # Loop again

    def _build_status_lines(self, structure_pairs: List[Dict], suggested_pairs: List[Dict]) -> List[str]:
        """Build status lines for display."""
        status_lines = []

        # Current pairs
        if structure_pairs:
            status_lines.append(f"Structure Pairs: {len(structure_pairs)} selected")
            for i, pair in enumerate(structure_pairs, 1):
                prmtop_name = Path(pair['prmtop']).name
                rst7_name = Path(pair['rst7']).name
                status_lines.append(f"  {i}. {prmtop_name} + {rst7_name}")
        else:
            status_lines.append("Structure Pairs: None selected")

        # Suggested pairs
        if suggested_pairs:
            status_lines.append("")
            status_lines.append(f"Auto-detected: {len(suggested_pairs)} pairs available")

        return status_lines

    def _build_actions_menu(
        self,
        structure_pairs: List[Dict],
        suggested_pairs: List[Dict]
    ) -> tuple[List[tuple[str, str]], Dict[str, any]]:
        """Build actions menu with option mapping."""
        actions = []
        option_map = {}

        # Suggested pair options
        for i, pair in enumerate(suggested_pairs):
            key = str(len(actions) + 1)
            match_type = pair.get('match_type', 'exact match')
            actions.append((key, f"Use {pair['name']} ({match_type})"))
            option_map[key] = ("use_suggested", i)

        # Use all suggested
        if len(suggested_pairs) > 1:
            key = str(len(actions) + 1)
            actions.append((key, f"Use all suggested pairs ({len(suggested_pairs)})"))
            option_map[key] = "use_all_suggested"

        # Browse for files
        key = str(len(actions) + 1)
        actions.append((key, "Browse for files"))
        option_map[key] = "browse_files"

        # Remove pair option
        if structure_pairs:
            key = str(len(actions) + 1)
            actions.append((key, "Remove structure pair"))
            option_map[key] = "remove"

        # Navigation
        if structure_pairs:
            actions.append(("n", "Next step (Workflow Configuration)"))
            option_map["n"] = "next"
        else:
            # Show dimmed next option
            actions.append(("n", "[grey50]Next step (requires structure pairs)[/grey50]"))

        actions.append(("h", "Help"))
        option_map["h"] = "help"

        actions.append(("x", "Exit workflow"))
        option_map["x"] = "exit"

        return actions, option_map

    def _handle_choice(
        self,
        choice: str,
        option_map: Dict,
        structure_pairs: List[Dict],
        suggested_pairs: List[Dict]
    ) -> str:
        """Handle user choice. Returns 'exit', 'next', or 'continue'."""

        action = option_map.get(choice)

        if action == "exit":
            return "exit"

        elif action == "next":
            if structure_pairs:
                return "next"
            else:
                self.console.print("[yellow]⚠ Please select at least one structure pair first[/yellow]")
                return "continue"

        elif action == "help":
            self._show_help()
            return "continue"

        elif isinstance(action, tuple) and action[0] == "use_suggested":
            # Add suggested pair
            idx = action[1]
            pair = suggested_pairs[idx]
            structure_pairs.append(pair)
            self.console.print(f"[green]✓ Added pair: {pair['name']}[/green]")
            return "continue"

        elif action == "use_all_suggested":
            # Add all suggested pairs
            for pair in suggested_pairs:
                if pair not in structure_pairs:
                    structure_pairs.append(pair)
            self.console.print(f"[green]✓ Added {len(suggested_pairs)} pairs[/green]")
            return "continue"

        elif action == "browse_files":
            self.console.print("[cyan]→ Opening file browser...[/cyan]")
            # In real implementation, would call file browser
            return "continue"

        elif action == "remove":
            # Remove pair
            if structure_pairs:
                self.console.print("\nSelect pair to remove:")
                for i, pair in enumerate(structure_pairs, 1):
                    self.console.print(f"  {i}. {pair['name']}")

                remove_options = {str(i + 1): f"Pair {i + 1}" for i in range(len(structure_pairs))}
                remove_options["c"] = "Cancel"
                remove_choice = prompt_with_context(
                    self.processor,
                    "Enter number to remove",
                    choices=[str(i) for i in range(1, len(structure_pairs) + 1)] + ["c"],
                    default="c",
                    module="MD Manager - Step 0",
                    description="Select structure pair to remove",
                    options_map=remove_options,
                )

                if remove_choice != "c":
                    removed = structure_pairs.pop(int(remove_choice) - 1)
                    self.console.print(f"[yellow]✓ Removed: {removed['name']}[/yellow]")

            return "continue"

        return "continue"

    def _show_help(self):
        """Show contextual help."""
        from rich.panel import Panel

        help_text = """
[bold cyan]Structure File Selection Help[/bold cyan]

[bold]What are structure files?[/bold]
• Topology (.prmtop, .parm7): Contains system parameters
• Coordinates (.rst7, .inpcrd): Contains atom positions

[bold]Auto-detection:[/bold]
• ProPrep automatically finds matching pairs
• Looks in workspace and working directory
• Match types:
  - [green]exact match[/green]: Same basename (system.prmtop + system.rst7)
  - [yellow]modified[/yellow]: Topology modified for const. pH (system_modified.prmtop + system.rst7)

[bold]Tips:[/bold]
• Use suggested pairs for quick setup
• Browse files if your files are elsewhere
• You can select multiple pairs for batch simulations

[grey50]Press Enter to continue...[/grey50]
"""
        self.console.print(Panel(help_text, border_style="cyan", padding=(1, 2)))
        input()

    def _simulate_detection(self) -> List[Dict]:
        """Simulate auto-detected pairs for demo."""
        return [
            {
                'name': 'WT_2LL6_4Ca',
                'prmtop': '/Users/mjgp/demo/new4/WT_2LL6_4Ca.prmtop',
                'rst7': '/Users/mjgp/demo/new4/WT_2LL6_4Ca.rst7',
                'match_type': 'exact match',
                'source': 'directory'
            },
            {
                'name': 'OxRed_HS_LS',
                'prmtop': '/Users/mjgp/demo/new4/OxRed_HS_LS.prmtop',
                'rst7': '/Users/mjgp/demo/new4/OxRed_HS_LS.rst7',
                'match_type': 'exact match',
                'source': 'directory'
            }
        ]


def demo_step0():
    """Run demo of Step 0 with new layout."""
    console = Console()

    console.print("\n[bold cyan]═══ PROTOTYPE: Step 0 - New Layout ═══[/bold cyan]\n")
    console.print("[grey50]This demonstrates the simplified structure file selection UX.[/grey50]\n")

    step0 = Step0Prototype(console)
    result = step0.execute_step0()

    if result:
        console.print(f"\n[green]✓ Step 0 completed with {len(result)} pairs[/green]")
    else:
        console.print("\n[yellow]Step 0 cancelled[/yellow]")


if __name__ == "__main__":
    demo_step0()
