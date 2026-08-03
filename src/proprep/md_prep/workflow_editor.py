"""
Workflow Editor

Interactive editor for custom workflows with add/edit/remove/reorder capabilities.
Integrates with template creation wizard and maintains dependency tracking.
"""

import copy
import os
import sys
import tempfile
import subprocess
from typing import Dict, List, Optional, Any, Tuple
from pathlib import Path
from datetime import datetime
import json

from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from proprep.utils.prompts import (
    prompt_with_context,
    confirm_with_context,
    int_prompt_with_context,
)

from .user_data_manager import UserDataManager
from .amber_wizard import AmberWizard


class WorkflowEditor:
    """
    Interactive editor for simulation protocols.

    Provides add/edit/remove/reorder functionality with dependency management
    and integration with template creation systems. Edits are stored inline
    on protocol steps (as mdin_content_override), not as separate template files.
    """

    def __init__(self, console: Console, user_data_manager: UserDataManager, processor=None):
        self.console = console
        self.user_data_manager = user_data_manager
        self.processor = processor

    def _resolve_step_mdin_content(self, step: Dict) -> Optional[str]:
        """Resolve the MDIN content for a protocol step.

        Resolution order:
        1. mdin_content_override — user edited via nano
        2. template_ref (builtin path) — read from disk, apply parameter_overrides
        3. parameter_overrides alone — generate from config dict (wizard-created)

        Returns:
            MDIN content string, or None if nothing can be resolved
        """
        # 1. Direct override takes priority
        if step.get('mdin_content_override'):
            return step['mdin_content_override']

        # 2. Builtin template reference
        template_ref = step.get('template_ref', '')
        if template_ref and template_ref.startswith('builtin/'):
            try:
                content, _ = self.user_data_manager.get_template_content(template_ref)
                overrides = step.get('parameter_overrides', {})
                if overrides:
                    content = self._apply_parameter_overrides(content, overrides)
                return content
            except Exception:
                pass

        # 3. Generate from parameter_overrides (wizard-created steps)
        overrides = step.get('parameter_overrides', {})
        if overrides:
            return self._generate_mdin_from_config(
                overrides,
                step.get('description', 'Custom template'),
                step.get('nmr_section', '')
            )

        # Legacy fallback: try custom_template_id if present
        custom_id = step.get('custom_template_id')
        if custom_id:
            template_data = self.user_data_manager.load_custom_template(custom_id)
            if template_data:
                ref = template_data.get('template', '')
                if ref and ref.startswith('builtin/'):
                    try:
                        content, _ = self.user_data_manager.get_template_content(ref)
                        return content
                    except Exception:
                        pass
                return self._generate_mdin_content(template_data)

        return None

    def _apply_parameter_overrides(self, mdin_content: str, overrides: Dict) -> str:
        """Apply parameter overrides to MDIN content by modifying values in &cntrl."""
        if not overrides:
            return mdin_content

        lines = mdin_content.split('\n')
        result_lines = []
        in_cntrl = False
        applied = set()

        for line in lines:
            stripped = line.strip()

            if stripped.startswith('&cntrl'):
                in_cntrl = True
                result_lines.append(line)
                continue
            elif stripped == '/' and in_cntrl:
                # Insert any overrides not yet applied
                for key in sorted(overrides.keys()):
                    if key not in applied:
                        value = overrides[key]
                        if isinstance(value, str) and ' ' in value:
                            result_lines.append(f"  {key}='{value}',")
                        else:
                            result_lines.append(f"  {key}={value},")
                in_cntrl = False
                result_lines.append(line)
                continue

            if in_cntrl and '=' in stripped:
                # Check if this parameter has an override
                param_name = stripped.split('=')[0].strip()
                if param_name in overrides:
                    value = overrides[param_name]
                    # Preserve comment if present
                    comment = ''
                    if '!' in stripped:
                        comment_idx = stripped.index('!')
                        comment = '  ' + stripped[comment_idx:]
                    if isinstance(value, str) and ' ' in value:
                        result_lines.append(f"  {param_name}='{value}',{comment}")
                    else:
                        result_lines.append(f"  {param_name}={value},{comment}")
                    applied.add(param_name)
                    continue

            result_lines.append(line)

        return '\n'.join(result_lines)

    def _generate_mdin_from_config(self, config: Dict, description: str = 'Custom template',
                                    nmr_section: str = '') -> str:
        """Generate MDIN content from a config dictionary."""
        lines = [description, "&cntrl"]
        for key in sorted(config.keys()):
            value = config[key]
            if isinstance(value, str) and ' ' in value:
                lines.append(f"  {key}='{value}',")
            else:
                lines.append(f"  {key}={value},")
        lines.append("/")
        if nmr_section:
            lines.append("")
            lines.extend(nmr_section.split('\n'))
        return '\n'.join(lines) + '\n'
        
    def edit_workflow(self, workflow: Dict) -> Optional[Dict]:
        """
        Main workflow editing interface with direct command support.

        Args:
            workflow: Workflow dictionary to edit

        Returns:
            Modified workflow dictionary or None if cancelled
        """
        edited_workflow = copy.deepcopy(workflow)

        while True:
            self._display_workflow_with_commands(edited_workflow)

            command = prompt_with_context(
                self.processor,
                "\n[bold]Your command[/bold]",
                module="Workflow Editor",
                description="Workflow editor command",
            ).strip().lower()

            # Parse and execute command
            try:
                if command == "save":
                    return edited_workflow
                elif command == "cancel":
                    confirm = prompt_with_context(
                        self.processor,
                        "Discard all changes?",
                        choices=["y", "n"],
                        default="n",
                        module="Workflow Editor",
                        description="Discard workflow edits",
                        options_map={"y": "Yes, discard", "n": "No, continue editing"},
                    )
                    if confirm == "y":
                        return None
                elif command == "preview":
                    self._preview_workflow_output(edited_workflow)
                elif command == "info":
                    self._edit_workflow_metadata(edited_workflow)
                elif command.startswith("edit "):
                    step_num = int(command.split()[1])
                    self._edit_step_direct(edited_workflow, step_num)
                elif command.startswith("name "):
                    step_num = int(command.split()[1])
                    self._rename_step(edited_workflow, step_num)
                elif command.startswith("files "):
                    step_num = int(command.split()[1])
                    self._edit_step_files(edited_workflow, step_num)
                elif command.startswith("remove "):
                    step_num = int(command.split()[1])
                    self._remove_step_direct(edited_workflow, step_num)
                elif command.startswith("add before "):
                    step_num = int(command.split()[2])
                    self._add_step_at_position(edited_workflow, step_num, before=True)
                elif command.startswith("add after "):
                    step_num = int(command.split()[2])
                    self._add_step_at_position(edited_workflow, step_num, before=False)
                elif command.startswith("move "):
                    # Parse "move 3 to 1"
                    parts = command.split()
                    from_pos = int(parts[1])
                    to_pos = int(parts[3])
                    self._move_step_direct(edited_workflow, from_pos, to_pos)
                elif command == "help" or command == "h":
                    self._show_command_help()
                else:
                    self.console.print("[yellow]Unknown command. Type 'help' for available commands.[/yellow]")
            except (IndexError, ValueError) as e:
                self.console.print(f"[red]Invalid command format. Type 'help' for examples.[/red]")
                
    def _display_workflow_with_commands(self, workflow: Dict):
        """Display workflow with command-based interface."""
        self.console.print(f"\n[bold cyan]Editing Workflow: {workflow.get('name', 'Unnamed')}[/bold cyan]")

        description = workflow.get('description', '')
        if description:
            self.console.print(f"[grey50]{description}[/grey50]\n")

        steps = workflow.get('steps', [])
        if not steps:
            self.console.print("[yellow]No steps in workflow[/yellow]\n")
        else:
            self.console.print("[bold]Workflow Steps:[/bold]")
            for i, step in enumerate(steps, 1):
                step_name = step.get('name', f"Step {i}")
                step_type = step.get('type', 'unknown')
                self.console.print(f"  {i}. {step_name} ({step_type})")

        self.console.print("\n[bold]Commands:[/bold]")
        self.console.print("  [cyan]edit <step>[/cyan]              Edit MDIN parameters (e.g., 'edit 2')")
        self.console.print("  [cyan]name <step>[/cyan]              Rename step (e.g., 'name 1')")
        self.console.print("  [cyan]files <step>[/cyan]             Change input/output files (e.g., 'files 3')")
        self.console.print("  [cyan]add before <step>[/cyan]        Add new step before (e.g., 'add before 2')")
        self.console.print("  [cyan]add after <step>[/cyan]         Add new step after (e.g., 'add after 3')")
        self.console.print("  [cyan]remove <step>[/cyan]            Remove step (e.g., 'remove 5')")
        self.console.print("  [cyan]move <from> to <pos>[/cyan]     Reorder (e.g., 'move 3 to 1')")
        self.console.print("  [cyan]preview[/cyan]                  Preview generated files")
        self.console.print("  [cyan]info[/cyan]                     Show workflow info")
        self.console.print("  [green]save[/green]                     Save and finish")
        self.console.print("  [yellow]cancel[/yellow]                   Cancel changes")
        self.console.print("  [grey50]help[/grey50]                     Show command examples")

    def _display_workflow_overview(self, workflow: Dict):
        """Display current workflow structure (legacy table view)."""
        self.console.print(f"\n[bold cyan]Editing Workflow: {workflow.get('name', 'Unnamed')}[/bold cyan]")

        description = workflow.get('description', '')
        if description:
            self.console.print(f"[grey50]{description}[/grey50]")

        steps = workflow.get('steps', [])
        if not steps:
            self.console.print("\n[yellow]No steps in workflow[/yellow]")
            return

        # Create workflow steps table
        table = Table(title="Workflow Steps")
        table.add_column("Step", style="cyan", width=6)
        table.add_column("Name", style="yellow")
        table.add_column("Type", style="green")
        table.add_column("Dependencies", style="blue")
        table.add_column("Input", style="magenta")
        table.add_column("Output", style="magenta")

        for i, step in enumerate(steps, 1):
            deps = ", ".join(step.get('dependencies', []))
            table.add_row(
                str(i),
                step.get('name', f"Step {i}"),
                step.get('type', 'unknown'),
                deps or "none",
                step.get('input_coord', ''),
                step.get('output_coord', '')
            )

        self.console.print(table)

    def _show_command_help(self):
        """Show detailed command examples."""
        self.console.print("\n[bold cyan]Command Examples:[/bold cyan]\n")
        self.console.print("[bold]Editing steps:[/bold]")
        self.console.print("  edit 2           - Edit MDIN parameters for step 2")
        self.console.print("  name 1           - Rename step 1")
        self.console.print("  files 3          - Change input/output files for step 3\n")
        self.console.print("[bold]Adding steps:[/bold]")
        self.console.print("  add before 2     - Add new step before step 2")
        self.console.print("  add after 3      - Add new step after step 3\n")
        self.console.print("[bold]Removing/Moving:[/bold]")
        self.console.print("  remove 5         - Remove step 5")
        self.console.print("  move 3 to 1      - Move step 3 to position 1\n")
        self.console.print("[bold]Other commands:[/bold]")
        self.console.print("  preview          - Preview generated MDIN files")
        self.console.print("  info             - Edit workflow name/description")
        self.console.print("  save             - Save changes and exit")
        self.console.print("  cancel           - Discard changes and exit")
        input("\nPress Enter to continue...")

    # Direct command methods
    def _edit_step_direct(self, workflow: Dict, step_num: int):
        """Edit MDIN parameters for a step directly."""
        steps = workflow.get('steps', [])
        if step_num < 1 or step_num > len(steps):
            self.console.print(f"[red]Invalid step number. Must be between 1 and {len(steps)}[/red]")
            return

        step = steps[step_num - 1]
        self._show_step_detail_and_edit(step, step_num)

    def _show_step_detail_and_edit(self, step: Dict, step_num: int):
        """Show step details with explanations and edit options. Loops until user chooses 'back'."""
        from pathlib import Path

        while True:
            self.console.print(f"\n[bold cyan]Editing: {step.get('name', 'Unnamed')} (Step {step_num})[/bold cyan]\n")

            # Show current configuration
            self.console.print("[bold]Current Configuration:[/bold]")
            self.console.print(f"  Name: {step.get('name', 'Unnamed')}")
            self.console.print(f"  Type: {step.get('type', 'unknown')}")
            desc = step.get('description', '')
            if desc:
                self.console.print(f"  Description: {desc}\n")

            self.console.print("  [bold]Files:[/bold]")
            self.console.print(f"    Input:  {step.get('input_coord', 'inpcrd')}")
            self.console.print(f"    Output: {step.get('output_coord', 'unknown')}\n")

            # Show MDIN parameter summary
            mdin_content = self._resolve_step_mdin_content(step)
            if mdin_content:
                self._show_mdin_parameters_from_content(mdin_content)

            self.console.print("\n[bold]What would you like to edit?[/bold]")
            self.console.print("  [cyan]params[/cyan]    Edit MDIN parameters (opens text editor)")
            self.console.print("  [cyan]name[/cyan]      Change step name")
            self.console.print("  [cyan]desc[/cyan]      Change description")
            self.console.print("  [cyan]files[/cyan]     Change input/output coordinate files")
            self.console.print("  [cyan]back[/cyan]      Back to workflow\n")

            choice = prompt_with_context(
                self.processor,
                "Enter choice",
                choices=["params", "name", "desc", "files", "back"],
                default="params",
                module="Workflow Editor",
                description="Step edit action",
                options_map={
                    "params": "Edit parameters",
                    "name": "Edit step name",
                    "desc": "Edit description",
                    "files": "Edit files",
                    "back": "Back",
                },
            )

            if choice == "params":
                self._edit_step_template_parameters(step)
                # Loop back to show updated step details
            elif choice == "name":
                new_name = prompt_with_context(
                    self.processor,
                    "New step name",
                    default=step.get('name', ''),
                    module="Workflow Editor",
                    description="New step name",
                )
                step['name'] = new_name
                self.console.print(f"[green]✅ Step renamed to '{new_name}'[/green]")
                # Loop back to show updated step details
            elif choice == "desc":
                new_desc = prompt_with_context(
                    self.processor,
                    "New description",
                    default=step.get('description', ''),
                    module="Workflow Editor",
                    description="New step description",
                )
                step['description'] = new_desc
                self.console.print(f"[green]✅ Description updated[/green]")
                # Loop back to show updated step details
            elif choice == "files":
                self._edit_step_files_interactive(step)
                # Loop back to show updated step details
            elif choice == "back":
                break  # Exit loop and return to main workflow menu

    def _show_mdin_parameters_from_content(self, content: str):
        """Show all &cntrl parameters extracted from MDIN content."""
        self.console.print("  [bold]MDIN Parameters:[/bold]")

        import re
        # Collect parameters in order of appearance
        params = []
        in_cntrl = False
        for line in content.split('\n'):
            stripped = line.strip()
            if stripped.startswith('&cntrl'):
                in_cntrl = True
                continue
            if in_cntrl and stripped == '/':
                break
            if in_cntrl:
                # Strip AMBER comments (! to end of line), but not ! inside quotes
                in_quote = False
                comment_pos = -1
                for ci, ch in enumerate(stripped):
                    if ch == "'":
                        in_quote = not in_quote
                    elif ch == '!' and not in_quote:
                        comment_pos = ci
                        break
                code_part = (stripped[:comment_pos] if comment_pos >= 0 else stripped).strip()
                if not code_part:
                    continue
                # Handle multiple params on one line, respecting quoted values
                # Match: key = 'quoted value', or key = unquoted_value,
                for match in re.finditer(r"(\w+)\s*=\s*('[^']*'|\"[^\"]*\"|[^,=]+?)(?:,|$)", code_part):
                    params.append((match.group(1).strip(), match.group(2).strip()))

        if params:
            max_key_len = max(len(k) for k, _ in params)
            for key, value in params:
                self.console.print(f"    {key:<{max_key_len}} = {value}")
        else:
            self.console.print("    [grey50]No &cntrl parameters found[/grey50]")

    def _rename_step(self, workflow: Dict, step_num: int):
        """Rename a step directly."""
        steps = workflow.get('steps', [])
        if step_num < 1 or step_num > len(steps):
            self.console.print(f"[red]Invalid step number. Must be between 1 and {len(steps)}[/red]")
            return

        step = steps[step_num - 1]
        old_name = step.get('name', f'Step {step_num}')
        new_name = prompt_with_context(
            self.processor,
            f"New name for '{old_name}'",
            default=old_name,
            module="Workflow Editor",
            description=f"Rename step from {old_name}",
        )
        step['name'] = new_name
        self.console.print(f"[green]✅ Renamed to '{new_name}'[/green]")

    def _edit_step_files(self, workflow: Dict, step_num: int):
        """Edit step input/output files directly."""
        steps = workflow.get('steps', [])
        if step_num < 1 or step_num > len(steps):
            self.console.print(f"[red]Invalid step number. Must be between 1 and {len(steps)}[/red]")
            return

        step = steps[step_num - 1]
        self._edit_step_files_interactive(step)

    def _edit_step_files_interactive(self, step: Dict):
        """Interactive file editing."""
        self.console.print(f"\n[bold]File Configuration for {step.get('name', 'Step')}:[/bold]")
        input_coord = prompt_with_context(
            self.processor,
            "Input coordinate file",
            default=step.get('input_coord', 'inpcrd'),
            module="Workflow Editor",
            description="Step input coordinate file",
        )
        output_coord = prompt_with_context(
            self.processor,
            "Output coordinate file",
            default=step.get('output_coord', 'output.rst'),
            module="Workflow Editor",
            description="Step output coordinate file",
        )
        step['input_coord'] = input_coord
        step['output_coord'] = output_coord
        self.console.print(f"[green]✅ Files updated[/green]")

    def _remove_step_direct(self, workflow: Dict, step_num: int):
        """Remove a step directly."""
        steps = workflow.get('steps', [])
        if step_num < 1 or step_num > len(steps):
            self.console.print(f"[red]Invalid step number. Must be between 1 and {len(steps)}[/red]")
            return

        step = steps[step_num - 1]
        step_name = step.get('name', f'Step {step_num}')
        confirm = prompt_with_context(
            self.processor,
            f"Remove '{step_name}'?",
            choices=["y", "n"],
            default="n",
            module="Workflow Editor",
            description=f"Confirm remove step {step_name}",
            options_map={"y": "Yes, remove", "n": "No, keep"},
        )
        if confirm == "y":
            steps.pop(step_num - 1)
            self._update_workflow_dependencies(workflow)
            self.console.print(f"[green]✅ Removed step '{step_name}'[/green]")

    def _add_step_at_position(self, workflow: Dict, step_num: int, before: bool = True):
        """Add a step before or after a position."""
        steps = workflow.get('steps', [])
        if step_num < 1 or step_num > len(steps):
            self.console.print(f"[red]Invalid step number. Must be between 1 and {len(steps)}[/red]")
            return

        position = step_num - 1 if before else step_num
        self.console.print(f"\n[bold cyan]Add New Step {'Before' if before else 'After'} Step {step_num}[/bold cyan]")

        # Prompt user for step creation method
        steps = workflow.get('steps', [])

        self.console.print("\n[bold]New Step Options:[/bold]")
        self.console.print("  1. Create new template via wizard")
        self.console.print("  2. Use existing template from library")
        self.console.print("  3. Copy existing step in workflow")
        self.console.print("\n[bold]Navigation:[/bold]")
        self.console.print("  c. Cancel")

        choice = prompt_with_context(
            self.processor,
            "Select option",
            choices=["1", "2", "3", "c"],
            default="1",
            module="Workflow Editor",
            description="Add step source (scratch/template/copy or cancel)",
            options_map={
                "1": "Create from scratch",
                "2": "From template",
                "3": "Copy existing step",
                "c": "Cancel",
            },
        )

        if choice == "c":
            return

        new_step = None
        if choice == "1":
            new_step = self._create_step_from_wizard()
        elif choice == "2":
            new_step = self._create_step_from_existing_template()
        elif choice == "3":
            new_step = self._create_step_from_copy(steps)

        if new_step:
            steps.insert(position, new_step)
            self._update_workflow_dependencies(workflow)
            self.console.print(f"[green]Added step '{new_step['name']}' at position {position + 1}[/green]")

    def _move_step_direct(self, workflow: Dict, from_pos: int, to_pos: int):
        """Move a step from one position to another."""
        steps = workflow.get('steps', [])
        if from_pos < 1 or from_pos > len(steps) or to_pos < 1 or to_pos > len(steps):
            self.console.print(f"[red]Invalid positions. Must be between 1 and {len(steps)}[/red]")
            return

        step = steps.pop(from_pos - 1)
        steps.insert(to_pos - 1, step)
        self._update_workflow_dependencies(workflow)
        self.console.print(f"[green]✅ Moved step from position {from_pos} to {to_pos}[/green]")

    def _show_editing_menu(self) -> str:
        """Display editing options menu."""
        self.console.print(f"\n[bold]Workflow Editing Options:[/bold]")
        self.console.print("  1. ➕ Add new step")
        self.console.print("  2. ✏️  Edit existing step") 
        self.console.print("  3. ❌ Remove step")
        self.console.print("  4. 🔄 Reorder steps")
        self.console.print("  5. 📝 Edit workflow info")
        self.console.print("  6. 👁️  Preview workflow output")
        
        self.console.print("\n[bold]Navigation:[/bold]")
        self.console.print("  s. ✅ Save and finish")
        self.console.print("  c. ❌ Cancel changes")
        
        return prompt_with_context(
            self.processor,
            "Select option",
            choices=["1", "2", "3", "4", "5", "6", "s", "c"],
            default="1",
            module="Workflow Editor",
            description="Top-level workflow edit action",
        )
                         
    def _add_workflow_step(self, workflow: Dict):
        """Add a new step to the workflow."""
        steps = workflow.get('steps', [])
        
        self.console.print(f"\n[bold cyan]Add New Workflow Step[/bold cyan]")
        
        # Step addition options
        self.console.print("\n[bold]New Step Options:[/bold]")
        self.console.print("  1. Create new template via wizard")
        self.console.print("  2. Use existing template from library")
        self.console.print("  3. Copy existing step in workflow")
        
        self.console.print("\n[bold]Navigation:[/bold]")
        self.console.print("  c. ← Cancel")
        
        choice = prompt_with_context(
            self.processor,
            "Select option",
            choices=["1", "2", "3", "c"],
            default="1",
            module="Workflow Editor",
            description="Add step source (scratch/template/copy or cancel)",
            options_map={
                "1": "Create from scratch",
                "2": "From template",
                "3": "Copy existing step",
                "c": "Cancel",
            },
        )
        
        if choice == "c":
            return
        elif choice == "1":
            new_step = self._create_step_from_wizard()
        elif choice == "2":
            new_step = self._create_step_from_existing_template()
        elif choice == "3":
            new_step = self._create_step_from_copy(steps)
        else:
            return
            
        if new_step:
            # Determine insertion position
            if not steps:
                position = 0
            else:
                position = self._choose_insertion_position(steps)
                
            # Insert step and update dependencies
            steps.insert(position, new_step)
            self._update_workflow_dependencies(workflow)
            
            self.console.print(f"[green]✅ Added step '{new_step['name']}' at position {position + 1}[/green]")
            
    def _create_step_from_wizard(self) -> Optional[Dict]:
        """Create new step using AmberWizard."""
        self.console.print("\n[bold]Template Creation Wizard[/bold]")

        try:
            wizard_config = AmberWizard.configure(console=self.console, processor=self.processor)

            if not wizard_config:
                return None

            # Get simulation type from wizard config
            selected_type = wizard_config.pop("_simulation_type", "custom")

            # Get step details
            step_name = prompt_with_context(
                self.processor,
                "Step name",
                default=f"New {selected_type}",
                module="Workflow Editor",
                description="New step name",
            )
            step_description = prompt_with_context(
                self.processor,
                "Step description",
                default="",
                module="Workflow Editor",
                description="New step description",
            )

            # Create protocol step with config stored inline
            step_id = f"step_new"
            return {
                'id': step_id,
                'name': step_name,
                'type': selected_type,
                'template_ref': '',  # No builtin reference — wizard-created
                'description': step_description,
                'dependencies': [],
                'input_coord': 'inpcrd',
                'output_coord': f"{step_id}.rst",
                'parameter_overrides': wizard_config,
                'mdin_content_override': None,
                'nmr_section': ''
            }

        except Exception as e:
            self.console.print(f"[red]Wizard failed: {e}[/red]")
            return None
            
    def _create_step_from_existing_template(self) -> Optional[Dict]:
        """Create step from a builtin template."""
        templates = self.user_data_manager.list_templates(show_builtin=True, show_user=False)
        if not templates:
            self.console.print("[yellow]No templates available[/yellow]")
            return None

        self.console.print("\n[bold]Available Templates:[/bold]")
        template_list = list(templates.items())

        for i, (template_id, template_data) in enumerate(template_list, 1):
            name = template_data.get('name', template_id)
            t_type = template_data.get('simulation_type') or template_data.get('type', 'unknown')
            description = template_data.get('description', '')

            self.console.print(f"  {i}. [cyan]{name}[/cyan] ({t_type})")
            if description:
                self.console.print(f"     [grey50]{description}[/grey50]")

        try:
            choice = int_prompt_with_context(
                self.processor,
                f"Select template (1-{len(template_list)})",
                module="Workflow Editor",
                description="Select template from library",
            )
            if 1 <= choice <= len(template_list):
                template_id, template_data = template_list[choice - 1]

                step_name = prompt_with_context(
                    self.processor,
                    "Step name",
                    default=template_data.get('name', 'Step'),
                    module="Workflow Editor",
                    description="Step name (from template)",
                )
                sim_type = template_data.get('simulation_type') or template_data.get('type', 'unknown')

                step_id = f"step_new"
                return {
                    'id': step_id,
                    'name': step_name,
                    'type': sim_type,
                    'template_ref': template_data.get('template_path', ''),
                    'description': template_data.get('description', ''),
                    'dependencies': [],
                    'input_coord': 'inpcrd',
                    'output_coord': f"{step_id}.rst",
                    'parameter_overrides': {},
                    'mdin_content_override': None,
                    'nmr_section': ''
                }
            return None
        except Exception:
            return None
            
    def _create_step_from_copy(self, existing_steps: List[Dict]) -> Optional[Dict]:
        """Create step by copying existing step in workflow."""
        if not existing_steps:
            self.console.print("[yellow]No existing steps to copy[/yellow]")
            return None
            
        self.console.print("\n[bold]Copy Existing Step:[/bold]")
        for i, step in enumerate(existing_steps, 1):
            self.console.print(f"  {i}. {step.get('name', f'Step {i}')} ({step.get('type', 'unknown')})")
            
        try:
            choice = int_prompt_with_context(
                self.processor,
                f"Select step to copy (1-{len(existing_steps)})",
                module="Workflow Editor",
                description="Select existing step to copy",
            )
            if 1 <= choice <= len(existing_steps):
                source_step = existing_steps[choice - 1]
                
                # Create copy
                copied_step = copy.deepcopy(source_step)
                
                # Update identifiers
                step_id = f"step_{len(existing_steps) + 1}"
                copied_step['id'] = step_id
                copied_step['name'] = prompt_with_context(
                    self.processor,
                    "Step name",
                    default=f"{copied_step['name']}_copy",
                    module="Workflow Editor",
                    description="Copied step name",
                )
                copied_step['output_coord'] = f"{step_id}.rst"
                copied_step['dependencies'] = []  # Will be recalculated
                
                return copied_step
            return None
        except:
            return None
            
    def _choose_insertion_position(self, steps: List[Dict]) -> int:
        """Choose where to insert new step in workflow."""
        self.console.print(f"\n[bold]Insert Position:[/bold]")
        self.console.print("  0. At beginning")
        
        for i, step in enumerate(steps, 1):
            self.console.print(f"  {i}. After '{step.get('name', f'Step {i}')}'")
            
        try:
            position = int_prompt_with_context(
                self.processor,
                f"Insert position (0-{len(steps)})",
                default=len(steps),
                module="Workflow Editor",
                description="Step insert position",
            )
            return max(0, min(position, len(steps)))
        except:
            return len(steps)
            
    def _edit_workflow_step(self, workflow: Dict):
        """Edit an existing workflow step."""
        steps = workflow.get('steps', [])
        if not steps:
            self.console.print("[yellow]No steps to edit[/yellow]")
            return
            
        self.console.print(f"\n[bold cyan]Edit Workflow Step[/bold cyan]")
        for i, step in enumerate(steps, 1):
            self.console.print(f"  {i}. {step.get('name', f'Step {i}')} ({step.get('type', 'unknown')})")
            
        try:
            choice = int_prompt_with_context(
                self.processor,
                f"Select step to edit (1-{len(steps)})",
                module="Workflow Editor",
                description="Select step to edit",
            )
            if 1 <= choice <= len(steps):
                step = steps[choice - 1]
                self._edit_single_step(step)
        except:
            pass
            
    def _edit_single_step(self, step: Dict):
        """Edit a single workflow step."""
        while True:
            self.console.print(f"\n[bold]Editing Step: {step.get('name', 'Unnamed')}[/bold]")
            
            # Display current step info
            table = Table(title="Current Step Configuration")
            table.add_column("Property", style="cyan")
            table.add_column("Value", style="yellow")
            
            table.add_row("Name", step.get('name', ''))
            table.add_row("Type", step.get('type', ''))
            table.add_row("Description", step.get('description', ''))
            table.add_row("Input Coord", step.get('input_coord', ''))
            table.add_row("Output Coord", step.get('output_coord', ''))
            
            self.console.print(table)
            
            self.console.print(f"\n[bold]Edit Options:[/bold]")
            self.console.print("  1. Edit step metadata")
            self.console.print("  2. Edit template parameters")
            self.console.print("  3. Edit coordinate files")
            
            self.console.print(f"\n[bold]Navigation:[/bold]")
            self.console.print("  b. ← Back to workflow")
            
            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=["1", "2", "3", "b"],
                default="1",
                module="Workflow Editor",
                description="Single-step edit action",
                options_map={
                    "1": "Edit step metadata",
                    "2": "Edit template parameters",
                    "3": "Edit coordinate files",
                    "b": "Back",
                },
            )

            if choice == "1":
                step['name'] = prompt_with_context(
                    self.processor,
                    "Step name",
                    default=step.get('name', ''),
                    module="Workflow Editor",
                    description="Step name (single-step edit)",
                )
                step['description'] = prompt_with_context(
                    self.processor,
                    "Description",
                    default=step.get('description', ''),
                    module="Workflow Editor",
                    description="Step description (single-step edit)",
                )

            elif choice == "2":
                self._edit_step_template_parameters(step)

            elif choice == "3":
                step['input_coord'] = prompt_with_context(
                    self.processor,
                    "Input coordinate file",
                    default=step.get('input_coord', ''),
                    module="Workflow Editor",
                    description="Step input coord (single-step edit)",
                )
                step['output_coord'] = prompt_with_context(
                    self.processor,
                    "Output coordinate file",
                    default=step.get('output_coord', ''),
                    module="Workflow Editor",
                    description="Step output coord (single-step edit)",
                )
                
            elif choice == "b":
                break
                
    def _edit_step_template_parameters(self, step: Dict):
        """Edit template parameters for a step using nano editor."""
        # Resolve current MDIN content from the protocol step
        mdin_content = self._resolve_step_mdin_content(step)
        if not mdin_content:
            self.console.print("[yellow]No template content available for this step[/yellow]")
            return

        self._direct_edit_step_mdin(step, mdin_content)

    def _direct_edit_step_mdin(self, step: Dict, mdin_content: str):
        """Open MDIN content in the user's preferred editor and store result as mdin_content_override.

        If $EDITOR is not set we try nano, then vi, then vim. If the chosen editor
        is killed by a signal (e.g. SIGSEGV from a broken ncurses/terminfo on a
        remote host) we also fall through to the next candidate.
        """
        try:
            # Create temporary file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.mdin', delete=False) as temp_file:
                temp_file.write(mdin_content)
                temp_file_path = temp_file.name

            user_editor = os.environ.get('EDITOR')
            if user_editor:
                candidates = [user_editor]
            else:
                candidates = ['nano', 'vi', 'vim']

            last_error = None
            edited = False

            for editor in candidates:
                self.console.print(f"\n[bold cyan]Opening MDIN file in {editor}...[/bold cyan]")
                if editor == 'nano':
                    self.console.print("[grey50]Make your changes, then press Ctrl+X to save and exit[/grey50]")
                else:
                    self.console.print("[grey50]Make your changes, save, and exit the editor[/grey50]")

                try:
                    result = subprocess.run(
                        [editor, temp_file_path],
                        stdin=sys.stdin,
                        stdout=sys.stdout,
                        stderr=sys.stderr,
                    )
                except FileNotFoundError:
                    last_error = f"Editor '{editor}' not found"
                    self.console.print(f"[yellow]{last_error}; trying next candidate...[/yellow]")
                    continue

                # Negative returncode on POSIX means killed by signal N (e.g. -11 = SIGSEGV)
                if result.returncode < 0:
                    signal_num = -result.returncode
                    last_error = f"Editor '{editor}' killed by signal {signal_num} (likely ncurses/terminfo issue)"
                    self.console.print(f"[yellow]{last_error}[/yellow]")
                    # Only auto-fall-through if the user did not explicitly choose this editor
                    if user_editor:
                        break
                    self.console.print("[grey50]Trying next editor...[/grey50]")
                    continue

                edited = True
                break

            # Read back the modified content (may be unchanged if editor failed)
            with open(temp_file_path, 'r') as f:
                modified_content = f.read()

            if not edited and modified_content.strip() == mdin_content.strip():
                self.console.print(
                    f"[red]Could not launch an editor successfully.[/red] "
                    f"[yellow]Last error: {last_error or 'unknown'}[/yellow]"
                )
                self.console.print(
                    "[grey50]Tip: set $EDITOR to a working editor before starting proprep "
                    "(e.g. 'export EDITOR=vim'), or check TERM and terminfo on this host.[/grey50]"
                )
            elif modified_content.strip() != mdin_content.strip():
                step['mdin_content_override'] = modified_content
                self.console.print("[green]MDIN parameters updated[/green]")
            elif result.returncode != 0:
                self.console.print(
                    f"[yellow]Editor exited with status {result.returncode}; no changes saved[/yellow]"
                )
            else:
                self.console.print("[grey50]No changes made[/grey50]")

            Path(temp_file_path).unlink()

        except Exception as e:
            self.console.print(f"[red]Edit failed: {e}[/red]")
            
    def _generate_mdin_content(self, template_data: Dict) -> str:
        """Generate mdin file content from template data."""
        config = template_data.get('config', {})
        description = template_data.get('description', 'Custom template')
        nmr_section = template_data.get('nmr_section', '')
        
        lines = [
            description,
            "&cntrl"
        ]
        
        # Add parameters in sorted order for consistency
        for key in sorted(config.keys()):
            value = config[key]
            if isinstance(value, str) and " " in value:
                lines.append(f"  {key}='{value}',")
            else:
                lines.append(f"  {key}={value},")
                
        lines.append("/")
        
        # Add NMR section if present
        if nmr_section:
            lines.append("")
            lines.extend(nmr_section.split('\n'))
            
        return "\n".join(lines) + "\n"
        
    def _parse_mdin_content(self, content: str) -> Optional[Dict[str, Any]]:
        """Parse mdin file content back to config dictionary."""
        config = {}
        
        lines = content.split('\n')
        in_cntrl_section = False
        
        for line in lines:
            line = line.strip()
            
            # Skip empty lines and comments
            if not line or line.startswith('#') or line.startswith('!'):
                continue
                
            # Check for control section
            if line.startswith('&cntrl'):
                in_cntrl_section = True
                continue
            elif line.startswith('/'):
                in_cntrl_section = False
                continue
                
            # Parse parameters in control section
            if in_cntrl_section and '=' in line:
                # Remove trailing comma and whitespace
                line = line.rstrip(',').strip()
                
                try:
                    key, value = line.split('=', 1)
                    key = key.strip()
                    value = value.strip()
                    
                    # Remove comments after value
                    if '!' in value:
                        value = value.split('!')[0].strip()
                    
                    # Remove trailing comma from value
                    value = value.rstrip(',')
                    
                    # Remove quotes if present
                    if (value.startswith("'") and value.endswith("'")) or \
                       (value.startswith('"') and value.endswith('"')):
                        value = value[1:-1]
                        config[key] = value
                    else:
                        # Try to convert to appropriate type
                        try:
                            if '.' in value:
                                config[key] = float(value)
                            else:
                                config[key] = int(value)
                        except ValueError:
                            config[key] = value
                            
                except ValueError:
                    continue  # Skip malformed lines
                    
        return config if config else None
            
    def _edit_individual_parameters(self, config: Dict, template_data: Dict, 
                                  template_id: str, step: Dict):
        """Edit individual template parameters."""
        while True:
            param_list = list(config.keys())
            if not param_list:
                self.console.print("[yellow]No parameters to edit[/yellow]")
                break
                
            self.console.print(f"\n[bold]Current Parameters:[/bold]")
            for i, param in enumerate(param_list, 1):
                value = config[param]
                self.console.print(f"  {i}. {param} = {value}")
                
            self.console.print(f"  {len(param_list) + 1}. Add new parameter")
            self.console.print(f"  {len(param_list) + 2}. ← Finish editing")
            
            try:
                choice = int_prompt_with_context(
                    self.processor,
                    f"Select parameter (1-{len(param_list) + 2})",
                    module="Workflow Editor",
                    description="Select template parameter to edit",
                )
                
                if choice <= len(param_list):
                    # Edit existing parameter
                    param = param_list[choice - 1]
                    current_value = config[param]
                    
                    new_value = prompt_with_context(
                        self.processor,
                        f"New value for {param}",
                        default=str(current_value),
                        module="Workflow Editor",
                        description=f"New value for parameter {param}",
                    )
                    
                    # Try to preserve type
                    if isinstance(current_value, int):
                        try:
                            config[param] = int(new_value)
                        except:
                            config[param] = new_value
                    elif isinstance(current_value, float):
                        try:
                            config[param] = float(new_value)
                        except:
                            config[param] = new_value
                    else:
                        config[param] = new_value
                        
                elif choice == len(param_list) + 1:
                    # Add new parameter
                    param_name = prompt_with_context(
                        self.processor,
                        "Parameter name",
                        module="Workflow Editor",
                        description="New parameter name",
                    )
                    param_value = prompt_with_context(
                        self.processor,
                        "Parameter value",
                        module="Workflow Editor",
                        description="New parameter value",
                    )
                    config[param_name] = param_value
                    
                else:
                    # Finish editing
                    break
                    
            except:
                break
                
        # Store updated config as parameter overrides on the step
        step['parameter_overrides'] = config
        self.console.print("[green]Parameters saved[/green]")
        
    def _remove_workflow_step(self, workflow: Dict):
        """Remove a step from the workflow."""
        steps = workflow.get('steps', [])
        if not steps:
            self.console.print("[yellow]No steps to remove[/yellow]")
            return
            
        self.console.print(f"\n[bold cyan]Remove Workflow Step[/bold cyan]")
        for i, step in enumerate(steps, 1):
            self.console.print(f"  {i}. {step.get('name', f'Step {i}')} ({step.get('type', 'unknown')})")
            
        try:
            choice = int_prompt_with_context(
                self.processor,
                f"Select step to remove (1-{len(steps)})",
                module="Workflow Editor",
                description="Select step to remove",
            )
            if 1 <= choice <= len(steps):
                removed_step = steps.pop(choice - 1)
                
                # Update dependencies
                self._update_workflow_dependencies(workflow)
                
                self.console.print(f"[green]✅ Removed step '{removed_step.get('name', 'step')}'[/green]")
        except:
            pass
            
    def _reorder_workflow_steps(self, workflow: Dict):
        """Reorder workflow steps using syntax like '1,3,2,4'."""
        steps = workflow.get('steps', [])
        if len(steps) < 2:
            self.console.print("[yellow]Need at least 2 steps to reorder[/yellow]")
            return
            
        self.console.print(f"\n[bold cyan]Reorder Workflow Steps[/bold cyan]")
        self.console.print("Current order:")
        for i, step in enumerate(steps, 1):
            self.console.print(f"  {i}. {step.get('name', f'Step {i}')} ({step.get('type', 'unknown')})")
            
        self.console.print(f"\n[bold]Reorder Instructions:[/bold]")
        self.console.print("Enter new order as comma-separated numbers")
        self.console.print(f"Example: '2,1,3' to swap first two steps")
        
        current_order = ",".join(str(i) for i in range(1, len(steps) + 1))
        new_order = prompt_with_context(
            self.processor,
            "New order",
            default=current_order,
            module="Workflow Editor",
            description="New step order (comma-separated 1-based indices)",
        )
        
        try:
            indices = [int(x.strip()) - 1 for x in new_order.split(',')]
            
            if len(indices) != len(steps) or set(indices) != set(range(len(steps))):
                self.console.print("[red]Invalid order - must include all steps exactly once[/red]")
                return
                
            # Reorder steps
            reordered_steps = [steps[i] for i in indices]
            workflow['steps'] = reordered_steps
            
            # Update dependencies
            self._update_workflow_dependencies(workflow)
            
            self.console.print("[green]✅ Steps reordered successfully[/green]")
            
        except ValueError:
            self.console.print("[red]Invalid format - use comma-separated numbers[/red]")
            
    def _edit_workflow_metadata(self, workflow: Dict):
        """Edit workflow name and description."""
        self.console.print(f"\n[bold cyan]Edit Workflow Metadata[/bold cyan]")
        
        workflow['name'] = prompt_with_context(
            self.processor,
            "Workflow name",
            default=workflow.get('name', ''),
            module="Workflow Editor",
            description="Workflow name",
        )
        workflow['description'] = prompt_with_context(
            self.processor,
            "Description",
            default=workflow.get('description', ''),
            module="Workflow Editor",
            description="Workflow description",
        )
        
        self.console.print("[green]✅ Workflow metadata updated[/green]")
        
    def _update_workflow_dependencies(self, workflow: Dict):
        """Update workflow step dependencies based on current order."""
        steps = workflow.get('steps', [])
        
        for i, step in enumerate(steps):
            # Update step IDs
            step['id'] = f"step_{i+1}"
            step['output_coord'] = f"step_{i+1}.rst"
            
            if i == 0:
                # First step uses the initial coordinate file
                step['input_coord'] = workflow.get('initial_coords', 'inpcrd')
                step['dependencies'] = []
            else:
                # Subsequent steps depend on previous step
                step['input_coord'] = f"step_{i}.rst"
                step['dependencies'] = [f"step_{i}"]
                
    def _preview_workflow_output(self, workflow: Dict):
        """Preview the workflow's generated output files."""
        steps = workflow.get('steps', [])
        if not steps:
            self.console.print("[yellow]No steps to preview[/yellow]")
            return
            
        self.console.print(f"\n[bold cyan]Workflow Output Preview[/bold cyan]")
        
        # Show file flow
        table = Table(title="File Flow")
        table.add_column("Step", style="cyan")
        table.add_column("Input", style="green")
        table.add_column("Output", style="yellow")
        table.add_column("Command Preview", style="white")
        
        for i, step in enumerate(steps, 1):
            input_coord = step.get('input_coord', '')
            output_coord = step.get('output_coord', '')
            step_type = step.get('type', 'unknown')
            
            # Generate sample command
            if step_type == 'minimization':
                cmd = f"pmemd -i step_{i}.in -p system.prmtop -c {input_coord} -o step_{i}.out -r {output_coord}"
            else:
                cmd = f"pmemd -i step_{i}.in -p system.prmtop -c {input_coord} -o step_{i}.out -r {output_coord} -x step_{i}.nc"
                
            table.add_row(f"Step {i}", input_coord, output_coord, cmd[:60] + "..." if len(cmd) > 60 else cmd)
            
        self.console.print(table)
        
        # Show template preview for selected step
        if confirm_with_context(
            self.processor,
            "Preview template content for a step?",
            default=False,
            module="Workflow Editor",
            description="Preview template content for a step",
        ):
            try:
                choice = int_prompt_with_context(
                    self.processor,
                    f"Select step (1-{len(steps)})",
                    default=1,
                    module="Workflow Editor",
                    description="Select step to preview template",
                )
                if 1 <= choice <= len(steps):
                    step = steps[choice - 1]
                    self._preview_step_template(step)
            except:
                pass
                
    def _preview_step_template(self, step: Dict):
        """Preview the MDIN content for a protocol step."""
        mdin_content = self._resolve_step_mdin_content(step)
        if not mdin_content:
            self.console.print("[yellow]No template content available for this step[/yellow]")
            return

        self.console.print(f"\n[bold]Template Preview: {step.get('name', 'Step')}[/bold]")
        self.console.print(Panel(mdin_content, title="MDIN content", border_style="green", expand=False))