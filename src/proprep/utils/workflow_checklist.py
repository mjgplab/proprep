"""
Workflow Checklist Framework

A reusable framework for multi-step workflows with:
- Visual checklist display with progress indicators
- State persistence for resume capability
- Step dependencies and re-run support
- Full workspace snapshot for seamless resumption

Usage:
    steps = [
        WorkflowStep(id="1", name="First Step", ...),
        WorkflowStep(id="2", name="Second Step", dependencies=["1"], ...),
    ]
    checklist = WorkflowChecklist(steps, executor=my_object, processor=processor)
    checklist.run()

Author: ProPrep Development Team
"""

import json
import os
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Callable, TYPE_CHECKING

from rich.console import Console, Group
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from proprep.utils.prompts import prompt_with_context, confirm_with_context

if TYPE_CHECKING:
    from proprep.core.processor import Processor


# =============================================================================
# Data Structures
# =============================================================================

@dataclass
class WorkflowStep:
    """
    Definition of a single workflow step.

    Steps are defined declaratively and can be easily added/removed/rearranged.
    """
    id: str                      # Unique identifier (e.g., "0a", "1")
    name: str                    # Display name (e.g., "Structure Filtering")
    description: str             # One-line description shown in checklist
    handler: str                 # Method name on executor (e.g., "_step_0a_...")
    section: str = ""            # Section header (e.g., "Structure Preparation")
    dependencies: List[str] = field(default_factory=list)  # Step IDs required first
    optional: bool = False       # Can be skipped?
    checkpoint: bool = False     # Is this a checkpoint (e.g., requires external calculation)?
    checkpoint_message: str = "" # Message to show at checkpoint
    info: str = ""               # Long-form educational text shown via [#i] keybind


@dataclass
class StepStatus:
    """Runtime status of a workflow step."""
    status: str = "pending"      # 'pending', 'in_progress', 'completed', 'skipped', 'failed'
    started_at: Optional[str] = None
    completed_at: Optional[str] = None
    output_summary: Optional[str] = None  # Brief description of outputs
    error_message: Optional[str] = None

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> 'StepStatus':
        return cls(**d)


@dataclass
class WorkflowState:
    """
    Complete state of a workflow run.

    Includes step statuses and full workspace snapshot for resumption.
    """
    workflow_id: str
    workflow_name: str
    working_directory: str
    created_at: str
    updated_at: str
    step_statuses: Dict[str, dict]  # step_id -> StepStatus as dict
    workspace_snapshot: Dict[str, Any]
    current_step: Optional[str] = None
    pdb_file: Optional[str] = None  # Original input file

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> 'WorkflowState':
        return cls(**d)

    def get_step_status(self, step_id: str) -> StepStatus:
        if step_id in self.step_statuses:
            return StepStatus.from_dict(self.step_statuses[step_id])
        return StepStatus()

    def set_step_status(self, step_id: str, status: StepStatus):
        self.step_statuses[step_id] = status.to_dict()
        self.updated_at = datetime.now().isoformat()


# =============================================================================
# Workspace Serialization Helpers
# =============================================================================

def _serialize_workspace(workspace) -> Dict[str, Any]:
    """
    Serialize workspace to JSON-compatible dict.

    Handles special types like Path, RedoxSite, etc.
    """
    if workspace is None:
        return {}

    snapshot = {}

    # Get all workspace keys
    if hasattr(workspace, '_data'):
        data = workspace._data
    elif hasattr(workspace, 'get_all'):
        data = workspace.get_all()
    else:
        return {}

    for key, value in data.items():
        try:
            snapshot[key] = _serialize_value(value)
        except Exception:
            # Skip non-serializable values
            pass

    return snapshot


def _serialize_value(value) -> Any:
    """Serialize a single value to JSON-compatible format."""
    if value is None:
        return None
    if isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Path):
        return {"__type__": "Path", "value": str(value)}
    if isinstance(value, dict):
        # JSON keys must be strings - convert tuple keys to string representation
        result = {}
        for k, v in value.items():
            if isinstance(k, tuple):
                # Convert tuple key to string: (1.0, 2.0, 3.0) -> "tuple:(1.0, 2.0, 3.0)"
                key_str = f"tuple:{k}"
            elif isinstance(k, (str, int, float, bool)):
                key_str = str(k) if not isinstance(k, str) else k
            else:
                key_str = str(k)
            result[key_str] = _serialize_value(v)
        return result
    if isinstance(value, (list, tuple)):
        serialized = [_serialize_value(v) for v in value]
        if isinstance(value, tuple):
            return {"__type__": "tuple", "value": serialized}
        return serialized
    if hasattr(value, 'to_dict'):
        return {"__type__": type(value).__name__, "value": value.to_dict()}
    if hasattr(value, '__dict__'):
        return {"__type__": type(value).__name__, "value": _serialize_value(vars(value))}
    # Fallback: convert to string
    return {"__type__": "str", "value": str(value)}


def _deserialize_workspace(snapshot: Dict[str, Any], workspace) -> None:
    """
    Restore workspace from serialized snapshot.

    Note: Complex objects (RedoxSite, etc.) are restored as dicts.
    The workflow should handle re-creating them if needed.
    """
    if workspace is None or not snapshot:
        return

    for key, value in snapshot.items():
        try:
            restored = _deserialize_value(value)
            workspace.set(key, restored)
        except Exception:
            pass


def _deserialize_value(value) -> Any:
    """Deserialize a single value from JSON format."""
    if value is None:
        return None
    if isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, dict):
        if "__type__" in value:
            type_name = value["__type__"]
            if type_name == "Path":
                return Path(value["value"])
            if type_name == "tuple":
                return tuple(_deserialize_value(v) for v in value["value"])
            # For other types, return the value dict (caller can reconstruct)
            return value
        # Handle dict with potentially tuple keys
        result = {}
        for k, v in value.items():
            # Check if key was a serialized tuple
            if isinstance(k, str) and k.startswith("tuple:"):
                # Parse tuple string back to tuple: "tuple:(1.0, 2.0, 3.0)" -> (1.0, 2.0, 3.0)
                try:
                    import ast
                    tuple_str = k[6:]  # Remove "tuple:" prefix
                    key = ast.literal_eval(tuple_str)
                except (ValueError, SyntaxError):
                    key = k  # Keep as string if parsing fails
            else:
                key = k
            result[key] = _deserialize_value(v)
        return result
    if isinstance(value, list):
        return [_deserialize_value(v) for v in value]
    return value


# =============================================================================
# Main Checklist Class
# =============================================================================

class WorkflowChecklist:
    """
    Interactive checklist-driven workflow runner.

    Displays a visual checklist, tracks progress, handles state persistence,
    and allows re-running steps.
    """

    STATE_FILENAME = "workflow_state.json"

    def __init__(
        self,
        steps: List[WorkflowStep],
        executor: Any,
        processor: Optional['Processor'] = None,
        workflow_name: str = "Workflow",
        console: Optional[Console] = None,
        state_dir: Optional[Path] = None,
        info_text: str = "",
    ):
        """
        Initialize the workflow checklist.

        Args:
            steps: List of WorkflowStep definitions
            executor: Object containing the step handler methods
            processor: ProPrep processor for workspace access
            workflow_name: Display name for the workflow
            console: Rich console for output
            state_dir: Directory for state file (defaults to current dir)
            info_text: Long-form module overview shown via the [i] keybind
        """
        self.steps = steps
        self.executor = executor
        self.processor = processor
        self.workflow_name = workflow_name
        self.console = console or Console()
        self.state_dir = Path(state_dir) if state_dir else None
        self.info_text = info_text

        self.state: Optional[WorkflowState] = None
        self.workspace = None
        if processor:
            self.workspace = processor._get_workspace()

    # =========================================================================
    # Main Entry Points
    # =========================================================================

    def run(self, pdb_file: Optional[str] = None) -> bool:
        """
        Run the workflow with interactive checklist.

        Args:
            pdb_file: Optional PDB file to process

        Returns:
            True if workflow completed successfully
        """
        # Check for existing state
        state_file = self._get_state_file_path()
        if state_file.exists():
            if self._offer_resume(state_file):
                # State loaded, continue from where we left off
                pass
            else:
                # Start fresh
                self._initialize_state(pdb_file)
        else:
            self._initialize_state(pdb_file)

        # Main loop
        return self._run_loop()

    def _run_loop(self) -> bool:
        """Main interactive loop."""
        while True:
            self._display_checklist()
            action = self._get_user_action()

            if action == 'quit':
                if confirm_with_context(
                    self.processor,
                    "Quit without saving? (Progress will be lost)",
                    default=False,
                    module=f"{self.workflow_name} Checklist",
                    description="Confirm quit workflow without saving",
                ):
                    return False
                continue

            elif action == 'save':
                self._save_state()
                self.console.print("[green]State saved. You can resume later.[/green]")
                return False

            elif action == 'view':
                self._launch_structure_viewer()
                continue

            elif action == 'info':
                self._show_info()
                continue

            elif action.startswith('info:'):
                self._show_info(step_id=action[5:])
                continue

            elif action == 'next':
                next_step = self._get_next_pending_step()
                if next_step:
                    self._run_step(next_step.id)
                else:
                    self.console.print("[green]All steps completed![/green]")
                    self._cleanup_state_file()
                    return True

            elif action.startswith('run:'):
                step_id = action[4:]
                self._run_step(step_id)

            elif action.startswith('skip:'):
                step_id = action[5:]
                step = self._get_step_by_id(step_id)
                if self._skip_step(step_id):
                    self.console.print(
                        f"[grey50]Skipped {step.name if step else step_id}.[/grey50]")

            elif action == 'done':
                self.console.print("[green]Workflow completed![/green]")
                self._cleanup_state_file()
                return True

    # =========================================================================
    # Display Methods
    # =========================================================================

    def _build_step_number_map(self) -> Dict[int, str]:
        """Build mapping from consecutive display numbers to step IDs."""
        return {i + 1: step.id for i, step in enumerate(self.steps)}

    def _display_checklist(self):
        """Display the checklist with current status.

        Layout: one row per step in an aligned grid (indicator | number |
        name | description) so each step occupies a single line and the
        descriptions line up in their own column. The name column is pinned
        to the longest step name across ALL sections so the description edge
        stays straight down the whole panel even though each section is its
        own grid. Uses terminal width (expand=True) to trade height for
        width -- descriptions live in horizontal space rather than on a
        second wrapped line per step.
        """
        self.console.print()

        # Pin the name column to the widest name (incl. "(optional)") across
        # every section so descriptions align globally, not just per-section.
        name_width = max(
            (len(step.name) + (11 if step.optional else 0)) for step in self.steps
        ) if self.steps else 0
        # Pin the number column too (e.g. "13.") so single- vs double-digit
        # sections start their name column at the same x.
        num_width = len(f"{len(self.steps)}.") if self.steps else 0

        def _new_grid() -> Table:
            g = Table.grid(padding=(0, 2, 0, 0), expand=True)
            g.add_column(no_wrap=True)                              # indicator
            g.add_column(no_wrap=True, justify="right",
                         min_width=num_width)                      # number
            g.add_column(no_wrap=True, min_width=name_width)       # name
            g.add_column(overflow="fold", ratio=1)                 # description
            return g

        renderables = []
        grid = None
        current_section = ""
        step_num = 0

        for step in self.steps:
            step_num += 1

            # Section header: close the running grid, emit a spacer + heading,
            # then start a fresh grid for the new section.
            if step.section and step.section != current_section:
                if grid is not None:
                    renderables.append(grid)
                    renderables.append(Text(""))
                renderables.append(
                    Text.from_markup(f"[bold blue]{step.section.upper()}[/bold blue]")
                )
                grid = _new_grid()
                current_section = step.section
            elif grid is None:
                # No section labels at all -- single unnamed grid.
                grid = _new_grid()

            # Step status indicator (orange, not yellow, for white-bg legibility)
            status = self.state.get_step_status(step.id) if self.state else StepStatus()

            if status.status == "completed":
                indicator = "[green]\\[✓][/green]"
            elif status.status == "in_progress":
                indicator = "[dark_orange3]\\[→][/dark_orange3]"
            elif status.status == "skipped":
                indicator = "\\[○]"
            elif status.status == "failed":
                indicator = "[red]\\[✗][/red]"
            else:
                indicator = "\\[ ]"

            optional_tag = " (optional)" if step.optional else ""
            # No explicit style on the number: bold-on-default-foreground emits a
            # bare SGR 1, which the "use bright colors for bold text" terminal
            # setting remaps to bright white -> invisible on a white background
            # (same trap as the repr.brace bolding handled in _get_user_action).
            # Default foreground is theme-adaptive, so it flips to black on white.
            grid.add_row(
                Text.from_markup(indicator),
                Text(f"{step_num}."),
                Text(f"{step.name}{optional_tag}"),
                Text.from_markup(f"[grey50]{step.description}[/grey50]"),
            )

            # Output/error summary on its own row, aligned under the description.
            if status.status == "completed" and status.output_summary:
                grid.add_row(
                    "", "", "",
                    Text.from_markup(f"[green]→ {status.output_summary}[/green]"),
                )
            elif status.status == "failed" and status.error_message:
                grid.add_row(
                    "", "", "",
                    Text.from_markup(f"[red]→ {status.error_message}[/red]"),
                )

        if grid is not None:
            renderables.append(grid)

        # Legend inside the panel
        renderables.append(Text(""))
        renderables.append(Text.from_markup(
            "[green]\\[✓][/green] Completed  "
            "[dark_orange3]\\[→][/dark_orange3] In Progress  "
            "\\[ ] Pending  \\[○] Skipped  [red]\\[✗][/red] Failed"
        ))

        panel = Panel(
            Group(*renderables),
            title=f"[bold]{self.workflow_name}[/bold]",
            border_style="blue",
            padding=(1, 2),
            expand=True,
        )
        self.console.print(panel)

    def _get_user_action(self) -> str:
        """Get user's next action."""
        # Build number-to-step mapping
        num_to_step = self._build_step_number_map()

        # Find next step and its display number
        next_step = self._get_next_pending_step()
        next_step_num = None
        if next_step:
            for num, step_id in num_to_step.items():
                if step_id == next_step.id:
                    next_step_num = num
                    break

        completed_count = sum(
            1 for s in self.steps
            if self.state and self.state.get_step_status(s.id).status == "completed"
        )

        line1 = []
        if next_step:
            line1.append(f"\\[n]ext ({next_step_num}. {next_step.name})")
        if completed_count == len(self.steps):
            line1.append("\\[d]one")
        line1.append("\\[v]iew structure")
        line1.append("\\[s]ave & quit")

        line2 = ["\\[q]uit", "\\[#] to run/re-run step"]
        if any(s.optional for s in self.steps):
            line2.append("\\[#s] skip optional step")
        if self.info_text or any(s.info for s in self.steps):
            line2.append("\\[i]nfo (or \\[#i] for step info)")

        # highlight=False: Rich's ReprHighlighter otherwise colors the
        # next-step number bold-cyan (repr.number) and styles the [n]/[v]/...
        # brackets bold-default (repr.brace), which the "bright bold text"
        # terminal setting remaps to invisible on a white background.
        self.console.print("  " + "  |  ".join(line1), highlight=False)
        self.console.print("  " + "  |  ".join(line2), highlight=False)
        self.console.print()

        action_options_map = {
            "n": "Next step",
            "v": "View structure",
            "s": "Save and quit",
            "q": "Quit",
            "d": "Done",
            "i": "Info / theory panel",
        }
        response = prompt_with_context(
            self.processor,
            "Select action",
            default="n" if next_step else "d",
            module=f"{self.workflow_name} Checklist",
            description="Select next checklist action",
            options_map=action_options_map,
        ).strip().lower()

        if response in ('n', 'next', ''):
            return 'next'
        elif response in ('v', 'view'):
            return 'view'
        elif response in ('s', 'save'):
            return 'save'
        elif response in ('q', 'quit'):
            return 'quit'
        elif response in ('d', 'done'):
            return 'done'
        elif response in ('i', 'info'):
            return 'info'
        else:
            # Per-step info: "<num>i" or "i<num>"
            stripped = response.replace(" ", "")
            num_part = None
            if stripped.endswith("i") and stripped[:-1].isdigit():
                num_part = stripped[:-1]
            elif stripped.startswith("i") and stripped[1:].isdigit():
                num_part = stripped[1:]
            if num_part is not None:
                step_num = int(num_part)
                if step_num in num_to_step:
                    return f'info:{num_to_step[step_num]}'
                self.console.print(f"[yellow]Invalid step number: {step_num}[/yellow]")
                return self._get_user_action()

            # Per-step skip: "<num>s" or "s<num>" (only marks optional steps).
            # Bare "s" is handled above as save, so it never reaches here.
            skip_part = None
            if stripped.endswith("s") and stripped[:-1].isdigit():
                skip_part = stripped[:-1]
            elif stripped.startswith("s") and stripped[1:].isdigit():
                skip_part = stripped[1:]
            if skip_part is not None:
                step_num = int(skip_part)
                if step_num in num_to_step:
                    return f'skip:{num_to_step[step_num]}'
                self.console.print(f"[yellow]Invalid step number: {step_num}[/yellow]")
                return self._get_user_action()

            # Check if it's a number
            try:
                step_num = int(response)
                if step_num in num_to_step:
                    return f'run:{num_to_step[step_num]}'
                else:
                    self.console.print(f"[yellow]Invalid step number: {step_num}[/yellow]")
                    return self._get_user_action()
            except ValueError:
                pass

            # Check if it's a valid step ID (for backwards compatibility)
            step_ids = [s.id for s in self.steps]
            if response in step_ids:
                return f'run:{response}'

            self.console.print(f"[yellow]Unknown action: {response}[/yellow]")
            return self._get_user_action()

    def _show_info(self, step_id: Optional[str] = None):
        """
        Display the educational info panel.

        With no step_id, shows the workflow-level overview (`info_text`).
        With a step_id, shows that step's `info` field.
        """
        if step_id is None:
            body = self.info_text.strip() if self.info_text else ""
            if not body:
                self.console.print("[yellow]No overview info defined for this workflow.[/yellow]")
                return
            num_to_step = self._build_step_number_map()
            steps_with_info = [
                (num, self._get_step_by_id(sid))
                for num, sid in num_to_step.items()
                if self._get_step_by_id(sid) and self._get_step_by_id(sid).info
            ]
            if steps_with_info:
                hint = "\n".join(
                    f"  • {num}i — {step.name}" for num, step in steps_with_info
                )
                body = f"{body}\n\n[bold blue]Per-step info:[/bold blue]\n{hint}"
            self.console.print()
            self.console.print(Panel(
                body,
                title=f"[bold]{self.workflow_name} — Overview & Theory[/bold]",
                border_style="cyan",
                padding=(1, 2),
                expand=False,
            ))
            self.console.print("[grey50]Press Enter to return to the checklist...[/grey50]")
            input()
            return

        step = self._get_step_by_id(step_id)
        if step is None:
            self.console.print(f"[yellow]Unknown step: {step_id}[/yellow]")
            return
        if not step.info:
            self.console.print(f"[yellow]No info defined for step '{step.name}'.[/yellow]")
            return
        # Build numbered title
        num_to_step = self._build_step_number_map()
        step_num = next((n for n, sid in num_to_step.items() if sid == step.id), None)
        title = f"[bold]Step {step_num}: {step.name}[/bold]" if step_num else f"[bold]{step.name}[/bold]"
        self.console.print()
        self.console.print(Panel(
            step.info.strip(),
            title=title,
            border_style="cyan",
            padding=(1, 2),
            expand=False,
        ))
        self.console.print("[grey50]Press Enter to return to the checklist...[/grey50]")
        input()

    def _launch_structure_viewer(self):
        """
        Launch the interactive structure viewer.

        Shows the current structure with redox sites highlighted if available.
        This allows users to inspect the metal coordination at any point in
        the workflow. Routed through the ``ViewerCoordinator`` singleton so
        the launch shares state with all other coordinator-driven viewer
        hooks across the session (annotations from earlier tools remain
        visible, the same viewer server is reused). Previously this method
        instantiated its own ``InteractiveStructureViewer``, which created
        a phantom viewer instance that didn't share annotation state with
        the rest of ProPrep.
        """
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer
        except ImportError:
            self.console.print("[yellow]Interactive viewer not available[/yellow]")
            return

        # Get PDB file - try prepared first, then original, then workspace
        pdb_file = None
        if hasattr(self.executor, '_final_pdb') and self.executor._final_pdb:
            pdb_file = str(self.executor._final_pdb)
        elif hasattr(self.executor, '_pdb_file') and self.executor._pdb_file:
            pdb_file = str(self.executor._pdb_file)
        elif self.state and getattr(self.state, 'pdb_file', None):
            pdb_file = self.state.pdb_file

        # Also try workspace keys
        if not pdb_file and self.processor:
            workspace = self.processor._get_workspace() if hasattr(self.processor, '_get_workspace') else {}
            for key in ['prepared_pdb', 'processed_pdb_file', 'local_pdb_file', 'rcsb_pdb_file']:
                pdb_file = workspace.get(key)
                if pdb_file and Path(pdb_file).exists():
                    break

        if not pdb_file or not Path(pdb_file).exists():
            self.console.print("[yellow]No PDB file available to view[/yellow]")
            return

        # Get redox sites if available
        redox_sites = []
        if hasattr(self.executor, 'redox_sites') and self.executor.redox_sites:
            redox_sites = self.executor.redox_sites
        elif self.processor:
            workspace = self.processor._get_workspace() if hasattr(self.processor, '_get_workspace') else {}
            # Check both keys - 'detected_redox_sites' is used for imported sites
            workspace_sites = workspace.get('detected_redox_sites', []) or workspace.get('redox_sites', [])
            if workspace_sites:
                redox_sites = workspace_sites

        # Wire processor for session recording, then route the launch through
        # the coordinator. force=True is correct: the user reached this method
        # by typing 'view' in the checklist menu, so this is a user-initiated
        # view.
        _viewer.set_processor(self.processor)
        _viewer.show_structure(pdb_file, force=True)

        # Build NGL selection for redox sites and overlay them as a single
        # focused-mode highlight (focused=True hides default protein/ligand/
        # ion/water reps, matching the original "show redox site only"
        # behavior). When no sites are detected, fall back to the default
        # structure view rather than an empty stage.
        if redox_sites:
            ngl_selections = set()
            for site in redox_sites:
                for atom in site.atoms:
                    chain = atom.chain if atom.chain else ""
                    if chain:
                        ngl_selections.add(f"({atom.resid}:{chain})")
                    else:
                        ngl_selections.add(f"({atom.resid})")

            if ngl_selections:
                ngl_selection = " or ".join(sorted(ngl_selections))
                _viewer.highlight(
                    ngl_selection,
                    style="ball+stick",
                    color="element",
                    label="redox_site",
                    focused=True,
                    force=True,
                )
                self.console.print(
                    f"\n[cyan]Launching viewer with {len(redox_sites)} redox site(s)[/cyan]"
                )
                self.console.print(
                    f"[grey50]Showing {len(ngl_selections)} residues (redox site only)[/grey50]"
                )
        else:
            self.console.print("\n[cyan]Launching viewer (no redox sites detected yet)[/cyan]")

        self.console.print("[grey50]Viewer launched. Press Enter to continue...[/grey50]")
        input()

    # =========================================================================
    # Step Execution
    # =========================================================================

    def _run_step(self, step_id: str) -> bool:
        """
        Run a specific workflow step.

        Args:
            step_id: ID of the step to run

        Returns:
            True if step completed successfully
        """
        step = self._get_step_by_id(step_id)
        if not step:
            self.console.print(f"[red]Unknown step: {step_id}[/red]")
            return False

        # Check dependencies
        missing_deps = self._check_dependencies(step)
        if missing_deps:
            self.console.print(f"[yellow]Step {step_id} requires these steps first: {', '.join(missing_deps)}[/yellow]")
            if not confirm_with_context(
                self.processor,
                "Run anyway? (may fail)",
                default=False,
                module=f"{self.workflow_name} Checklist",
                description=f"Confirm run step {step_id} despite missing dependencies",
            ):
                return False

        # Check if re-running a completed step
        status = self.state.get_step_status(step_id)
        if status.status == "completed":
            if not confirm_with_context(
                self.processor,
                f"Step {step_id} is already completed. Re-run it?",
                default=False,
                module=f"{self.workflow_name} Checklist",
                description=f"Confirm re-run completed step {step_id}",
            ):
                return False
            # Clear downstream steps if re-running
            self._clear_downstream_steps(step_id)

        # Update status to in_progress
        status = StepStatus(
            status="in_progress",
            started_at=datetime.now().isoformat()
        )
        self.state.set_step_status(step_id, status)
        self._auto_save()

        # Run the step.
        # Display the user-facing checklist position (1..N), not the internal
        # step id. Internal ids (e.g. "param-4"/"param-6"/"param-5") are ordered
        # historically and can read out of sequence; the number map follows the
        # displayed checklist order. Purely cosmetic — ids/dependencies/state
        # keys are unchanged.
        num_to_step = self._build_step_number_map()
        step_num = next((n for n, sid in num_to_step.items() if sid == step.id), None)
        step_label = f"Step {step_num}" if step_num else f"Step {step.id}"
        self.console.print()
        self.console.print(f"[bold cyan]═══ {step_label}: {step.name} ═══[/bold cyan]")
        self.console.print(f"[grey50]{step.description}[/grey50]")
        self.console.print()

        try:
            # Get the handler method
            handler = getattr(self.executor, step.handler, None)
            if handler is None:
                raise AttributeError(f"Handler method not found: {step.handler}")

            # Call the handler
            result = handler()

            # Handle checkpoint steps - only show checkpoint message if handler requests it
            # Handler returns {"checkpoint": True} when external calculation is needed
            show_checkpoint = (
                step.checkpoint and
                isinstance(result, dict) and
                result.get('checkpoint', False)
            )

            if show_checkpoint:
                self.console.print()
                self.console.print(Panel(
                    step.checkpoint_message or "External calculation required. Save and resume later.",
                    title="[yellow]Checkpoint[/yellow]",
                    border_style="yellow",
                    expand=False,
                ))
                status = StepStatus(
                    status="completed",
                    started_at=status.started_at,
                    completed_at=datetime.now().isoformat(),
                    output_summary="Checkpoint reached - external calculation required"
                )
            else:
                # Normal completion
                output_summary = None
                if isinstance(result, dict) and 'summary' in result:
                    output_summary = result['summary']
                elif isinstance(result, str):
                    output_summary = result

                status = StepStatus(
                    status="completed",
                    started_at=status.started_at,
                    completed_at=datetime.now().isoformat(),
                    output_summary=output_summary
                )

            self.state.set_step_status(step_id, status)
            self._snapshot_workspace()
            self._auto_save()

            self.console.print()
            self.console.print(f"[green]✓ {step_label} completed[/green]")
            return True

        except Exception as e:
            status = StepStatus(
                status="failed",
                started_at=status.started_at,
                completed_at=datetime.now().isoformat(),
                error_message=str(e)
            )
            self.state.set_step_status(step_id, status)
            self._auto_save()

            self.console.print()
            self.console.print(f"[red]✗ {step_label} failed: {e}[/red]")
            import traceback
            self.console.print(f"[grey50]{traceback.format_exc()}[/grey50]")
            return False

    def _skip_step(self, step_id: str) -> bool:
        """Mark a step as skipped."""
        step = self._get_step_by_id(step_id)
        if not step or not step.optional:
            self.console.print(f"[red]Step {step_id} cannot be skipped[/red]")
            return False

        status = StepStatus(
            status="skipped",
            completed_at=datetime.now().isoformat()
        )
        self.state.set_step_status(step_id, status)
        self._auto_save()
        return True

    # =========================================================================
    # State Management
    # =========================================================================

    def _initialize_state(self, pdb_file: Optional[str] = None):
        """Initialize fresh workflow state."""
        import uuid

        self.state = WorkflowState(
            workflow_id=str(uuid.uuid4())[:8],
            workflow_name=self.workflow_name,
            working_directory=os.getcwd(),
            created_at=datetime.now().isoformat(),
            updated_at=datetime.now().isoformat(),
            step_statuses={},
            workspace_snapshot={},
            pdb_file=pdb_file,
        )

        # Initialize all steps as pending
        for step in self.steps:
            self.state.set_step_status(step.id, StepStatus())

    def _offer_resume(self, state_file: Path) -> bool:
        """Offer to resume from existing state."""
        try:
            with open(state_file, 'r') as f:
                data = json.load(f)

            state = WorkflowState.from_dict(data)

            # Show resume info
            self.console.print()
            self.console.print(Panel(
                f"[bold blue]Found saved workflow state[/bold blue]\n\n"
                f"Workflow: {state.workflow_name}\n"
                f"Started: {state.created_at[:19].replace('T', ' ')}\n"
                f"Last updated: {state.updated_at[:19].replace('T', ' ')}\n"
                f"Working directory: {state.working_directory}",
                title="Resume Available",
                border_style="cyan",
                expand=False,
            ))

            # Count progress
            completed = sum(
                1 for sid, s in state.step_statuses.items()
                if s.get('status') == 'completed'
            )
            total = len(self.steps)

            self.console.print(f"  Progress: {completed}/{total} steps completed")
            self.console.print()

            if confirm_with_context(
                self.processor,
                "Resume from saved state?",
                default=True,
                module=f"{self.workflow_name} Checklist",
                description="Resume workflow from saved state",
            ):
                self.state = state
                # Restore workspace
                if self.workspace and state.workspace_snapshot:
                    _deserialize_workspace(state.workspace_snapshot, self.workspace)
                return True
            else:
                # Start fresh - remove old state
                state_file.unlink()
                return False

        except Exception as e:
            self.console.print(f"[yellow]Could not load saved state: {e}[/yellow]")
            return False

    def _save_state(self):
        """Save current state to file."""
        if not self.state:
            return

        self._snapshot_workspace()

        state_file = self._get_state_file_path()
        with open(state_file, 'w') as f:
            json.dump(self.state.to_dict(), f, indent=2)

        self.console.print(f"[grey50]State saved to {state_file}[/grey50]")

    def _auto_save(self):
        """Auto-save after each step."""
        self._save_state()

    def _snapshot_workspace(self):
        """Take a snapshot of current workspace."""
        if self.workspace and self.state:
            self.state.workspace_snapshot = _serialize_workspace(self.workspace)

    def _cleanup_state_file(self):
        """Remove state file after successful completion."""
        state_file = self._get_state_file_path()
        if state_file.exists():
            state_file.unlink()

    # =========================================================================
    # Helper Methods
    # =========================================================================

    def _get_state_file_path(self) -> Path:
        """Get the path to the state file, using state_dir if set."""
        if self.state_dir:
            return self.state_dir / self.STATE_FILENAME
        return Path(self.STATE_FILENAME)

    def _get_step_by_id(self, step_id: str) -> Optional[WorkflowStep]:
        """Get step definition by ID."""
        for step in self.steps:
            if step.id == step_id:
                return step
        return None

    def _get_next_pending_step(self) -> Optional[WorkflowStep]:
        """Get the next step that should be run."""
        for step in self.steps:
            status = self.state.get_step_status(step.id) if self.state else StepStatus()
            if status.status == "pending":
                return step
        return None

    def _check_dependencies(self, step: WorkflowStep) -> List[str]:
        """Check if step's dependencies are satisfied. Returns list of missing deps."""
        missing = []
        for dep_id in step.dependencies:
            dep_status = self.state.get_step_status(dep_id) if self.state else StepStatus()
            if dep_status.status not in ("completed", "skipped"):
                missing.append(dep_id)
        return missing

    def _clear_downstream_steps(self, step_id: str):
        """Clear status of steps that depend on this one (when re-running)."""
        # Find all steps that come after this one
        found = False
        for step in self.steps:
            if step.id == step_id:
                found = True
                continue
            if found:
                # Clear this step's status
                self.state.set_step_status(step.id, StepStatus())

        self.console.print("[grey50]Cleared downstream step statuses[/grey50]")
