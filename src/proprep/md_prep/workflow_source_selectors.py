"""
Protocol Source Selectors

Three paths for protocol source selection in Step 1:
1. Predefined Protocols - Browse standard multi-step protocols by category
2. Template Catalog Browser - Browse builtin templates with protocol provenance
3. Template Creation Wizard - Create steps from scratch via AmberWizard
"""

import json
from typing import Dict, List, Optional, Any, Tuple
from pathlib import Path
from collections import defaultdict
from datetime import datetime

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
from .workflow_loader import WorkflowLoader, WorkflowPreset
from .amber_wizard import AmberWizard


class PredefinedWorkflowSelector:
    """Handles selection from predefined protocols with category grouping."""

    def __init__(self, console: Console, user_data_manager: UserDataManager, processor=None):
        self.console = console
        self.user_data_manager = user_data_manager
        self.processor = processor
        self.workflow_loader = WorkflowLoader(console=console)

    def select_workflow(self, structure_pair: Dict) -> Optional[Dict]:
        """
        Select from predefined protocols.

        Returns:
            Configuration dict with selected workflow and steps
        """
        structure_name = structure_pair['name']

        workflows = self.workflow_loader.get_available_workflows()
        if not workflows:
            self.console.print("[yellow]No predefined protocols available[/yellow]")
            return None

        categories = self._group_workflows_by_category(workflows)

        while True:
            self.console.print(f"\n[bold cyan]Select Predefined Protocol[/bold cyan]")
            self.console.print(f"[grey50]Structure: {structure_name}[/grey50]\n")

            self.console.print("[bold]Available Protocols:[/bold]\n")

            all_workflows = []
            current_index = 1

            # Show only builtin protocols (user-saved are in option 4)
            builtin_categories = {k: v for k, v in categories.items() if k != 'custom'}

            if builtin_categories:
                for category in sorted(builtin_categories.keys()):
                    for workflow_id, workflow in builtin_categories[category]:
                        all_workflows.append((workflow_id, workflow))
                        step_count = len(workflow.steps) if hasattr(workflow, 'steps') else 0
                        self.console.print(f"  {current_index}. {workflow.name} ({step_count} steps)")
                        if hasattr(workflow, 'description'):
                            desc = workflow.description[:80] + "..." if len(workflow.description) > 80 else workflow.description
                            self.console.print(f"     [grey50]{desc}[/grey50]")
                        current_index += 1
                self.console.print()

            self.console.print("   Enter number to select, or 'b' to go back")
            self.console.print()

            valid_choices = [str(i) for i in range(1, len(all_workflows) + 1)] + ["b"]
            choice = prompt_with_context(
                self.processor,
                "Your selection",
                choices=valid_choices,
                module="MD Manager - Workflow Source",
                description="Select predefined workflow",
            ).strip()

            if choice == "b":
                return None

            workflow_idx = int(choice) - 1
            workflow_id, selected_workflow = all_workflows[workflow_idx]

            result = self._show_workflow_and_get_action(workflow_id, selected_workflow, structure_name)
            if result:
                return result

    def _show_workflow_and_get_action(self, workflow_id: str, workflow, structure_name: str) -> Optional[Dict]:
        """Show protocol details and get user action."""
        from .layout_helpers import TemplatePreviewFormatter

        if not hasattr(workflow, 'steps') or not workflow.steps:
            self.console.print("[red]No steps found in protocol[/red]")
            return None

        previewer = TemplatePreviewFormatter(self.console)

        while True:
            self.console.print(f"\n")
            header = f"{workflow.name}\n{workflow.description if hasattr(workflow, 'description') else ''}"
            self.console.print(Panel(header, border_style="cyan", padding=(0, 1), expand=False))

            self.console.print("\n[bold]Protocol Steps:[/bold]")
            import textwrap
            term_width = self.console.width or 100
            indent = "     "
            wrap_width = max(40, term_width - len(indent))

            for i, step in enumerate(workflow.steps, 1):
                step_name = step.name if hasattr(step, 'name') else f"Step {i}"
                step_type = step.type if hasattr(step, 'type') else 'unknown'
                self.console.print(f"  {i}. {step_name} ({step_type})")
                step_desc = getattr(step, 'description', '')
                if step_desc:
                    wrapped = textwrap.fill(step_desc, width=wrap_width,
                                            initial_indent=indent, subsequent_indent=indent)
                    self.console.print(f"[grey50]{wrapped}[/grey50]")

            total_time = self._compute_workflow_total_time(workflow.steps)
            if total_time:
                self.console.print(f"\n[bold]Total simulation time:[/bold] {total_time}")

            self.console.print("\n[bold]Commands:[/bold]")
            self.console.print("   p <n>  Preview step template (e.g., 'p 3')")
            self.console.print("   u      Use as-is")
            self.console.print("   s      Select specific steps only")
            self.console.print("   e      Edit protocol (customize steps/parameters)")
            self.console.print("   b      Back to protocol list")
            self.console.print()

            choice = prompt_with_context(
                self.processor,
                "Enter choice",
                module="MD Manager - Workflow Source",
                description="Workflow browser action (select, preview, or back)",
            ).strip().lower()

            if choice.startswith('p '):
                try:
                    step_num = int(choice.split()[1])
                    if 1 <= step_num <= len(workflow.steps):
                        step = workflow.steps[step_num - 1]
                        template_ref = getattr(step, 'template', None)
                        tdata = {'name': step.name, 'description': getattr(step, 'description', '')}

                        if template_ref and template_ref.startswith('builtin/'):
                            tdata['template_path'] = template_ref

                        if 'template_path' in tdata:
                            self.console.print()
                            previewer.preview_template(tdata)
                            input("\nPress Enter to continue...")
                        else:
                            self.console.print("[yellow]No template content available for this step[/yellow]")
                    else:
                        self.console.print(f"[red]Invalid step number. Must be between 1 and {len(workflow.steps)}[/red]")
                except (ValueError, IndexError):
                    self.console.print("[red]Invalid preview command. Use 'p <number>'[/red]")
                continue

            elif choice == "b":
                return None

            elif choice == "u":
                return {
                    'source_type': 'predefined',
                    'parent_workflow_id': workflow_id,
                    'workflow_name': workflow.name,
                    'selected_steps': list(workflow.steps),
                    'original_workflow': workflow,
                    'edit_requested': False
                }

            elif choice == "s":
                selected_steps = self._select_specific_steps(workflow.steps)
                if not selected_steps:
                    continue
                return {
                    'source_type': 'predefined',
                    'parent_workflow_id': workflow_id,
                    'workflow_name': workflow.name,
                    'selected_steps': selected_steps,
                    'original_workflow': workflow,
                    'edit_requested': False
                }

            elif choice == "e":
                return {
                    'source_type': 'predefined',
                    'parent_workflow_id': workflow_id,
                    'workflow_name': workflow.name,
                    'selected_steps': list(workflow.steps),
                    'original_workflow': workflow,
                    'edit_requested': True
                }

            else:
                self.console.print("[yellow]Unknown command. Enter u, s, e, b, or 'p <number>'[/yellow]")

    def _compute_workflow_total_time(self, steps) -> Optional[str]:
        """Compute total simulation time across all protocol steps."""
        from .layout_helpers import TemplatePreviewFormatter
        total_ps = 0.0
        has_time = False
        previewer = TemplatePreviewFormatter()

        for step in steps:
            template_ref = getattr(step, 'template', None)
            if not template_ref or not template_ref.startswith('builtin/'):
                continue

            try:
                tdata = {'template_path': template_ref}
                content = previewer._load_mdin_content(tdata)
                if not content:
                    continue

                params = previewer._parse_mdin_parameters(content)
                imin = int(params.get('imin', '0'))
                if imin == 1:
                    continue

                nstlim = params.get('nstlim')
                dt = params.get('dt')
                if nstlim and dt:
                    total_ps += int(nstlim) * float(dt)
                    has_time = True
            except Exception:
                continue

        if not has_time:
            return None

        if total_ps >= 1000:
            return f"{total_ps / 1000:.1f} ns"
        else:
            return f"{total_ps:.1f} ps"

    def _group_workflows_by_category(self, workflows: Dict) -> Dict[str, List[Tuple[str, Any]]]:
        """Group workflows by their category field."""
        categories = defaultdict(list)
        for workflow_id, workflow in workflows.items():
            category = getattr(workflow, 'category', 'General')
            categories[category].append((workflow_id, workflow))
        return dict(sorted(categories.items()))

    def _select_specific_steps(self, steps: List) -> Optional[List]:
        """Allow user to select specific steps from protocol."""
        self.console.print(f"\n[bold]Select Steps (comma-separated numbers):[/bold]")
        self.console.print("Example: 1,3,4 to select steps 1, 3, and 4")

        for i, step in enumerate(steps, 1):
            self.console.print(f"  {i}. {step.name} ({step.type})")

        selection = prompt_with_context(
            self.processor,
            "Enter step numbers",
            default="1",
            module="MD Manager - Workflow Source",
            description="Comma-separated step numbers to select",
        )

        try:
            step_indices = [int(x.strip()) - 1 for x in selection.split(',')]
            selected_steps = []
            for idx in step_indices:
                if 0 <= idx < len(steps):
                    selected_steps.append(steps[idx])
                else:
                    self.console.print(f"[red]Invalid step number: {idx + 1}[/red]")
                    return None
            return selected_steps
        except ValueError:
            self.console.print("[red]Invalid format. Please use comma-separated numbers.[/red]")
            return None


class TemplateCatalogBrowser:
    """
    Browse builtin templates with protocol provenance labels.

    Shows each template with its parent protocol name (e.g.,
    "Protein Equil: Initial Heating") so the user knows which
    protocol each template comes from.
    """

    def __init__(self, console: Console, user_data_manager: UserDataManager, processor=None):
        self.console = console
        self.user_data_manager = user_data_manager
        self.processor = processor
        self.workflow_loader = WorkflowLoader(console=console)

    def build_protocol(self, structure_pair: Dict) -> Optional[Dict]:
        """
        Build a custom protocol by selecting templates from the catalog.

        Returns:
            Configuration dict with selected templates and source_type='catalog'
        """
        structure_name = structure_pair['name']

        catalog = self._build_template_catalog()
        if not catalog:
            self.console.print("[yellow]No templates available[/yellow]")
            return None

        # Group catalog entries by protocol
        by_protocol = defaultdict(list)
        for path, entry in catalog.items():
            by_protocol[entry['protocol_name']].append(entry)

        # Sort protocols by their discovery order, then steps within each by step_number
        sorted_protocols = sorted(
            by_protocol.keys(),
            key=lambda name: min(e['protocol_order'] for e in by_protocol[name])
        )
        for entries in by_protocol.values():
            entries.sort(key=lambda e: e.get('step_number', ''))

        selected = []

        self.console.print(f"\n[bold cyan]Template Catalog[/bold cyan]")
        self.console.print(f"[grey50]Structure: {structure_name}[/grey50]")
        self.console.print("\nBuild your protocol by selecting steps from the catalog.\n")

        while True:
            # Show current selection
            if selected:
                self.console.print(f"[bold green]Current Protocol ({len(selected)} step{'s' if len(selected) != 1 else ''}):[/bold green]")
                for i, entry in enumerate(selected, 1):
                    provenance = f"{entry['protocol_name']}: " if entry.get('protocol_name') else ""
                    self.console.print(f"  {i}. {provenance}{entry['step_name']}")
                self.console.print()
            else:
                self.console.print("[grey50]No steps selected yet[/grey50]\n")

            # Build flat numbered list with protocol headers
            self.console.print("[bold]Template Catalog:[/bold]\n")
            flat_list = []
            display_idx = 1

            for protocol_name in sorted_protocols:
                entries = by_protocol[protocol_name]
                step_word = "step" if len(entries) == 1 else "steps"
                self.console.print(f"[bold cyan]{protocol_name} ({len(entries)} {step_word}):[/bold cyan]")

                for entry in entries:
                    flat_list.append(entry)
                    step_num = entry.get('step_number', '')
                    type_tag = entry['step_type']
                    num_str = f"{display_idx}.".rjust(4)
                    self.console.print(f"  {num_str}  {step_num}  {entry['step_name']}  [grey50][{type_tag}][/grey50]")
                    if entry.get('description'):
                        desc = entry['description'][:70] + "..." if len(entry['description']) > 70 else entry['description']
                        self.console.print(f"          [grey50]{desc}[/grey50]")
                    display_idx += 1
                self.console.print()

            self.console.print("[bold]Commands:[/bold]")
            self.console.print("  <n>        Add step (e.g., '3')")
            self.console.print("  <n>,<n>    Add multiple steps (e.g., '1,3,5')")
            self.console.print("  p <n>      Preview template")
            self.console.print("  c <n> <n>  Compare two templates")
            if selected:
                self.console.print("  r          Remove last added step")
                self.console.print("  d          Done building protocol")
            self.console.print("  b          Back to protocol source")
            self.console.print()

            choice = prompt_with_context(
                self.processor,
                "Your choice",
                module="MD Manager - Workflow Source",
                description="Template catalog action",
            ).strip().lower()

            if choice == 'b':
                return None

            if choice == 'd' and selected:
                # Convert catalog entries to the template_info format expected by ProtocolCreator
                selected_templates = []
                for entry in selected:
                    selected_templates.append({
                        'id': entry['template_path'],
                        'data': {
                            'name': entry['step_name'],
                            'description': entry.get('description', ''),
                            'template_path': entry['template_path'],
                        },
                        'type': entry['step_type']
                    })
                return {
                    'source_type': 'catalog',
                    'selected_templates': selected_templates,
                    'template_types': [t['type'] for t in selected_templates]
                }

            if choice == 'r' and selected:
                removed = selected.pop()
                self.console.print(f"[yellow]Removed: {removed['step_name']}[/yellow]")
                continue

            if choice.startswith('p '):
                try:
                    idx = int(choice.split()[1]) - 1
                    if 0 <= idx < len(flat_list):
                        from .layout_helpers import TemplatePreviewFormatter
                        previewer = TemplatePreviewFormatter(self.console)
                        tdata = {
                            'name': flat_list[idx]['step_name'],
                            'description': flat_list[idx].get('description', ''),
                            'template_path': flat_list[idx]['template_path']
                        }
                        self.console.print()
                        previewer.preview_template(tdata)
                        input("\nPress Enter to continue...")
                    else:
                        self.console.print("[red]Invalid template number[/red]")
                except (ValueError, IndexError):
                    self.console.print("[red]Invalid preview command. Use 'p <number>'[/red]")
                continue

            if choice.startswith('c '):
                try:
                    parts = choice.split()
                    idx1 = int(parts[1]) - 1
                    idx2 = int(parts[2]) - 1
                    if 0 <= idx1 < len(flat_list) and 0 <= idx2 < len(flat_list):
                        from .layout_helpers import TemplatePreviewFormatter
                        previewer = TemplatePreviewFormatter(self.console)
                        tdata1 = {'name': flat_list[idx1]['step_name'], 'template_path': flat_list[idx1]['template_path']}
                        tdata2 = {'name': flat_list[idx2]['step_name'], 'template_path': flat_list[idx2]['template_path']}
                        previewer.compare_templates(tdata1, tdata2)
                        input("\nPress Enter to continue...")
                    else:
                        self.console.print("[red]Invalid template numbers[/red]")
                except (ValueError, IndexError):
                    self.console.print("[red]Invalid compare command. Use 'c <n1> <n2>'[/red]")
                continue

            # Number selection (single or comma-separated)
            try:
                if ',' in choice:
                    indices = [int(x.strip()) - 1 for x in choice.split(',')]
                else:
                    indices = [int(choice) - 1]

                for idx in indices:
                    if 0 <= idx < len(flat_list):
                        selected.append(flat_list[idx])
                        name = flat_list[idx]['step_name']
                        self.console.print(f"[green]Added: {name}[/green]")
                    else:
                        self.console.print(f"[red]Invalid number: {idx + 1}[/red]")
            except ValueError:
                self.console.print("[red]Invalid input[/red]")

    def _build_template_catalog(self) -> Dict[str, Dict]:
        """Build catalog mapping each builtin template to its parent protocol."""
        catalog = {}

        # Map templates to their parent workflows (protocol order preserved)
        workflows = self.workflow_loader.get_available_workflows()
        protocol_idx = 0
        for wf_id, workflow in workflows.items():
            if not hasattr(workflow, 'steps'):
                continue
            for step in workflow.steps:
                template_path = getattr(step, 'template', '')
                if template_path and template_path.startswith('builtin/'):
                    if template_path not in catalog:
                        # Extract step number from filename (e.g., "00" from "00_initial_minimization.mdin")
                        filename = template_path.rsplit('/', 1)[-1]
                        step_number = filename.split('_', 1)[0] if filename[:1].isdigit() else ""

                        catalog[template_path] = {
                            'protocol_name': workflow.name,
                            'protocol_order': protocol_idx,
                            'step_name': step.name,
                            'step_type': step.type,
                            'step_number': step_number,
                            'description': getattr(step, 'description', ''),
                            'template_path': template_path
                        }
            protocol_idx += 1

        # Include any builtin templates not referenced by any workflow
        all_builtins = self.user_data_manager.list_templates(show_builtin=True, show_user=False)
        for tid, tdata in all_builtins.items():
            path = tdata.get('template_path', '')
            if path and path not in catalog:
                filename = path.rsplit('/', 1)[-1]
                step_number = filename.split('_', 1)[0] if filename[:1].isdigit() else ""
                # Extract protocol dir name from path
                parts = path.split('/')
                protocol_dir = parts[1] if len(parts) >= 3 else "unknown"

                catalog[path] = {
                    'protocol_name': protocol_dir.replace('_', ' ').title(),
                    'protocol_order': 999,
                    'step_name': tdata.get('name', ''),
                    'step_type': tdata.get('simulation_type') or tdata.get('type', 'unknown'),
                    'step_number': step_number,
                    'description': tdata.get('description', ''),
                    'template_path': path
                }

        return catalog


# Keep old name as alias for backward compatibility
StandaloneTemplateSelector = TemplateCatalogBrowser
CustomWorkflowBuilder = TemplateCatalogBrowser


class TemplateCreationWizard:
    """Handles new step creation via the AmberWizard."""

    def __init__(self, console: Console, user_data_manager: UserDataManager, processor=None):
        self.console = console
        self.user_data_manager = user_data_manager
        self.processor = processor

    def create_template(self, structure_pair: Dict) -> Optional[Dict]:
        """
        Create a new step via AmberWizard.

        Returns config with parameter_overrides stored inline — no custom
        template file is created.
        """
        structure_name = structure_pair['name']

        self.console.print(f"\n[bold cyan]Create New Step[/bold cyan]")
        self.console.print(f"[grey50]Structure: {structure_name}[/grey50]")

        try:
            wizard_config = AmberWizard.configure(console=self.console, processor=self.processor)

            if not wizard_config:
                return None

            # Get simulation type from wizard config
            selected_type = wizard_config.pop("_simulation_type", "custom")

            template_name = prompt_with_context(
                self.processor,
                "Step name",
                default=f"Custom_{selected_type}_{structure_name}",
                module="MD Manager - Workflow Source",
                description="Custom template step name",
            )
            template_description = prompt_with_context(
                self.processor,
                "Step description",
                default=f"Custom {selected_type} step",
                module="MD Manager - Workflow Source",
                description="Custom template step description",
            )

            # Return config with inline data — no custom template file created
            return {
                'source_type': 'new',
                'new_template': {
                    'id': '',
                    'data': {
                        'name': template_name,
                        'description': template_description,
                        'type': selected_type,
                        'config': wizard_config,
                    },
                    'type': selected_type
                }
            }

        except Exception as e:
            self.console.print(f"[red]Step creation failed: {e}[/red]")
            return None


class SavedProtocolSelector:
    """Load a user-saved protocol or import one from a file."""

    def __init__(self, console: Console, user_data_manager: UserDataManager, processor=None):
        self.console = console
        self.user_data_manager = user_data_manager
        self.processor = processor
        self.workflow_loader = WorkflowLoader(console=console)

    def select_protocol(self, structure_pair: Dict) -> Optional[Dict]:
        """
        Select from user-saved protocols or import from file.

        Returns:
            Configuration dict with selected workflow and steps, or None
        """
        structure_name = structure_pair['name']

        while True:
            self.console.print(f"\n[bold cyan]Load Saved Protocol[/bold cyan]")
            self.console.print(f"[grey50]Structure: {structure_name}[/grey50]\n")

            # Get user-saved protocols only
            saved = self._get_saved_protocols()

            if saved:
                self.console.print("[bold]Your Saved Protocols:[/bold]\n")
                all_protocols = []
                for i, (wf_id, workflow, metadata) in enumerate(saved, 1):
                    all_protocols.append((wf_id, workflow))
                    step_count = len(workflow.steps) if hasattr(workflow, 'steps') else 0
                    created = metadata.get('created', '')
                    if created:
                        try:
                            from datetime import datetime as dt
                            date_str = dt.fromisoformat(created).strftime('%Y-%m-%d')
                        except (ValueError, TypeError):
                            date_str = created[:10]
                    else:
                        date_str = ''
                    source = metadata.get('source', '')
                    date_info = f"  [grey50]{date_str}[/grey50]" if date_str else ""
                    source_info = f"  [grey50]({source})[/grey50]" if source else ""
                    self.console.print(
                        f"  {i}. {workflow.name} ({step_count} steps){date_info}{source_info}"
                    )
                    if hasattr(workflow, 'description') and workflow.description:
                        desc = workflow.description[:80] + "..." if len(workflow.description) > 80 else workflow.description
                        self.console.print(f"     [grey50]{desc}[/grey50]")
                self.console.print()
            else:
                all_protocols = []
                self.console.print("[grey50]No saved protocols found.[/grey50]\n")

            # Always show import option
            import_num = len(all_protocols) + 1
            self.console.print(f"  {import_num}. [cyan]Import protocol from file[/cyan]")
            self.console.print(f"     [grey50]Load a protocol JSON exported from another project[/grey50]\n")

            self.console.print("   Enter number to select, or 'b' to go back")
            self.console.print()

            valid_choices = [str(i) for i in range(1, import_num + 1)] + ["b"]
            choice = prompt_with_context(
                self.processor,
                "Your selection",
                choices=valid_choices,
                module="MD Manager - Workflow Source",
                description="Select predefined workflow",
            ).strip()

            if choice == "b":
                return None

            idx = int(choice) - 1

            if idx == len(all_protocols):
                # Import from file
                result = self._import_from_file(structure_name)
                if result:
                    return result
                continue

            workflow_id, selected_workflow = all_protocols[idx]
            result = self._show_protocol_and_get_action(workflow_id, selected_workflow, structure_name)
            if result:
                return result

    def _get_saved_protocols(self) -> List[Tuple[str, Any, Dict]]:
        """Get user-saved protocols with metadata."""
        all_workflows = self.workflow_loader.get_available_workflows()
        if not all_workflows:
            return []

        result = []
        for wf_id, workflow in all_workflows.items():
            category = getattr(workflow, 'category', '')
            if category == 'custom':
                metadata = self.user_data_manager.workflow_metadata.get(wf_id, {})
                result.append((wf_id, workflow, metadata))
        return result

    def _show_protocol_and_get_action(self, workflow_id: str, workflow, structure_name: str) -> Optional[Dict]:
        """Show protocol details and get user action (reuses PredefinedWorkflowSelector pattern)."""
        from .layout_helpers import TemplatePreviewFormatter
        import textwrap

        if not hasattr(workflow, 'steps') or not workflow.steps:
            self.console.print("[red]No steps found in protocol[/red]")
            return None

        previewer = TemplatePreviewFormatter(self.console)

        while True:
            self.console.print(f"\n")
            header = f"{workflow.name}\n{workflow.description if hasattr(workflow, 'description') else ''}"
            self.console.print(Panel(header, border_style="yellow", padding=(0, 1), expand=False))

            self.console.print("\n[bold]Protocol Steps:[/bold]")
            term_width = self.console.width or 100
            indent = "     "
            wrap_width = max(40, term_width - len(indent))

            for i, step in enumerate(workflow.steps, 1):
                step_name = step.name if hasattr(step, 'name') else f"Step {i}"
                step_type = step.type if hasattr(step, 'type') else 'unknown'
                self.console.print(f"  {i}. {step_name} ({step_type})")
                step_desc = getattr(step, 'description', '')
                if step_desc:
                    wrapped = textwrap.fill(step_desc, width=wrap_width,
                                            initial_indent=indent, subsequent_indent=indent)
                    self.console.print(f"[grey50]{wrapped}[/grey50]")

            self.console.print("\n[bold]Commands:[/bold]")
            self.console.print("   p <n>  Preview step template (e.g., 'p 3')")
            self.console.print("   u      Use as-is")
            self.console.print("   s      Select specific steps only")
            self.console.print("   e      Edit protocol (customize steps/parameters)")
            self.console.print("   b      Back to protocol list")
            self.console.print()

            choice = prompt_with_context(
                self.processor,
                "Enter choice",
                module="MD Manager - Workflow Source",
                description="Workflow browser action (select, preview, or back)",
            ).strip().lower()

            if choice.startswith('p '):
                try:
                    step_num = int(choice.split()[1])
                    if 1 <= step_num <= len(workflow.steps):
                        step = workflow.steps[step_num - 1]
                        template_ref = getattr(step, 'template', None)
                        tdata = {'name': step.name, 'description': getattr(step, 'description', '')}
                        if template_ref and template_ref.startswith('builtin/'):
                            tdata['template_path'] = template_ref
                        if 'template_path' in tdata:
                            self.console.print()
                            previewer.preview_template(tdata)
                            input("\nPress Enter to continue...")
                        else:
                            self.console.print("[yellow]No template content available for this step[/yellow]")
                    else:
                        self.console.print(f"[red]Invalid step number. Must be between 1 and {len(workflow.steps)}[/red]")
                except (ValueError, IndexError):
                    self.console.print("[red]Invalid preview command. Use 'p <number>'[/red]")
                continue

            elif choice == "b":
                return None

            elif choice == "u":
                return {
                    'source_type': 'saved',
                    'parent_workflow_id': workflow_id,
                    'workflow_name': workflow.name,
                    'selected_steps': list(workflow.steps),
                    'original_workflow': workflow,
                    'edit_requested': False
                }

            elif choice == "s":
                selected_steps = self._select_specific_steps(workflow.steps)
                if not selected_steps:
                    continue
                return {
                    'source_type': 'saved',
                    'parent_workflow_id': workflow_id,
                    'workflow_name': workflow.name,
                    'selected_steps': selected_steps,
                    'original_workflow': workflow,
                    'edit_requested': False
                }

            elif choice == "e":
                return {
                    'source_type': 'saved',
                    'parent_workflow_id': workflow_id,
                    'workflow_name': workflow.name,
                    'selected_steps': list(workflow.steps),
                    'original_workflow': workflow,
                    'edit_requested': True
                }

            else:
                self.console.print("[yellow]Unknown command. Enter u, s, e, b, or 'p <number>'[/yellow]")

    def _select_specific_steps(self, steps: List) -> Optional[List]:
        """Allow user to select specific steps from protocol."""
        self.console.print(f"\n[bold]Select Steps (comma-separated numbers):[/bold]")
        self.console.print("Example: 1,3,4 to select steps 1, 3, and 4")

        for i, step in enumerate(steps, 1):
            self.console.print(f"  {i}. {step.name} ({step.type})")

        selection = prompt_with_context(
            self.processor,
            "Enter step numbers",
            default="1",
            module="MD Manager - Workflow Source",
            description="Comma-separated step numbers to select",
        )

        try:
            step_indices = [int(x.strip()) - 1 for x in selection.split(',')]
            selected_steps = []
            for idx in step_indices:
                if 0 <= idx < len(steps):
                    selected_steps.append(steps[idx])
                else:
                    self.console.print(f"[red]Invalid step number: {idx + 1}[/red]")
                    return None
            return selected_steps
        except ValueError:
            self.console.print("[red]Invalid format. Please use comma-separated numbers.[/red]")
            return None

    def _import_from_file(self, structure_name: str) -> Optional[Dict]:
        """Import a protocol from a JSON export file."""
        self.console.print("\n[bold cyan]Import Protocol from File[/bold cyan]")
        self.console.print("[grey50]Enter the path to a protocol JSON file exported from another project.[/grey50]\n")

        file_path = prompt_with_context(
            self.processor,
            "File path (or 'b' to go back)",
            module="MD Manager - Workflow Source",
            description="Saved protocol file path",
        )
        if file_path.strip().lower() == 'b':
            return None

        path = Path(file_path.strip())
        if not path.exists():
            self.console.print(f"[red]File not found: {path}[/red]")
            return None

        try:
            with open(path, 'r') as f:
                import_data = json.load(f)
        except (json.JSONDecodeError, IOError) as e:
            self.console.print(f"[red]Could not read file: {e}[/red]")
            return None

        # Handle both export format (type/metadata/content) and raw protocol format (steps)
        if 'type' in import_data and import_data['type'] == 'workflow':
            # Export format from user_data_manager.export_content()
            content_id = self.user_data_manager.import_content(path)
            if not content_id:
                return None
            self.console.print(f"[green]✓ Protocol imported successfully[/green]")
            # Reload and present to user
            workflow = self.user_data_manager.load_custom_workflow(content_id)
            if workflow:
                # Wrap in WorkflowPreset if it's a raw dict
                if isinstance(workflow, dict):
                    workflow = WorkflowPreset.from_dict(content_id, workflow)
                return self._show_protocol_and_get_action(content_id, workflow, structure_name)
        elif 'steps' in import_data:
            # Raw protocol dict — import it
            name = import_data.get('name', path.stem)
            content_id = self.user_data_manager.save_custom_workflow(import_data)
            if not content_id:
                self.console.print("[red]Failed to import protocol[/red]")
                return None
            self.console.print(f"[green]✓ Protocol '{name}' imported successfully[/green]")
            workflow = self.user_data_manager.load_custom_workflow(content_id)
            if workflow:
                if isinstance(workflow, dict):
                    workflow = WorkflowPreset.from_dict(content_id, workflow)
                return self._show_protocol_and_get_action(content_id, workflow, structure_name)
        else:
            self.console.print("[red]Unrecognized file format. Expected a protocol with 'steps' or an exported workflow.[/red]")

        return None
