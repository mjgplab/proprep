"""
Protocol-Centric Step 1 Manager

Step 1 interface that handles multiple structure pairs with a protocol-first approach.
Protocols reference builtin templates directly and carry edits inline.
"""

import json
import uuid
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, asdict
import copy

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
from .workflow_loader import WorkflowLoader


@dataclass
class StructureWorkflowAssignment:
    """Assignment of a protocol to a structure pair."""
    structure_name: str
    structure_pair: Dict[str, Any]
    custom_workflow_id: str
    custom_workflow_name: str
    source_type: str  # "predefined", "catalog", "wizard", "saved"
    parent_workflow_id: Optional[str] = None


class WorkflowCentricStep1Manager:
    """
    Protocol-centric Step 1 manager.

    Handles multiple structure pairs sequentially, creating simulation protocols
    that reference builtin templates directly. No custom template copies are
    created — edits are stored inline on protocol steps.
    """

    def __init__(self, console: Console, user_data_manager: UserDataManager, workspace=None, processor=None):
        self.console = console
        self.user_data_manager = user_data_manager
        self.workspace = workspace
        self.processor = processor
        self.workflow_loader = WorkflowLoader(console=console)

        # Layout formatters for consistent UX
        from .layout_helpers import WorkflowLayoutFormatter, TemplatePreviewFormatter
        self.layout = WorkflowLayoutFormatter(console)
        self.template_preview = TemplatePreviewFormatter(console)

        # Multi-structure tracking
        self.current_structure_index = 0
        self.total_structures = 0
        self.structure_assignments: List[StructureWorkflowAssignment] = []

        # Source selectors
        from .workflow_source_selectors import (
            PredefinedWorkflowSelector,
            TemplateCatalogBrowser,
            TemplateCreationWizard,
            SavedProtocolSelector,
        )
        self.predefined_selector = PredefinedWorkflowSelector(console, user_data_manager)
        self.catalog_browser = TemplateCatalogBrowser(console, user_data_manager)
        self.template_wizard = TemplateCreationWizard(console, user_data_manager)
        self.saved_protocol_selector = SavedProtocolSelector(console, user_data_manager)

        # Protocol creator
        from .custom_workflow_creator import ProtocolCreator
        self.protocol_creator = ProtocolCreator(console, user_data_manager)

    def execute_step1(
        self,
        structure_pairs: List[Dict],
        multi_structure_mode: Optional[str] = None,
    ) -> List[StructureWorkflowAssignment]:
        """
        Main entry point for protocol-centric Step 1.

        Args:
            structure_pairs: List of structure pair dictionaries from Step 0
            multi_structure_mode: The wizard-level choice ("all", "separate",
                or None). When set, the apply-to-all prompt is skipped because
                the user already answered the equivalent question upstream.

        Returns:
            List of structure-workflow assignments for Step 2
        """
        self.total_structures = len(structure_pairs)
        self.structure_assignments = []

        for i, structure_pair in enumerate(structure_pairs):
            self.current_structure_index = i
            assignment = self._process_structure_pair(structure_pair)
            if assignment:
                self.structure_assignments.append(assignment)

            if i == 0 and len(structure_pairs) > 1 and assignment:
                if multi_structure_mode == "all":
                    self._apply_assignment_to_remaining(structure_pairs[1:], assignment)
                    break
                if multi_structure_mode == "separate":
                    continue
                if self._offer_apply_to_all(structure_pairs[1:], assignment):
                    break

        return self._finalize_step1()

    def display_step_summary(self):
        """Display protocol selection overview."""
        self.console.print(f"\n[bold cyan]Protocol Selection Overview[/bold cyan]")
        self.console.print(f"[grey50]Configure the simulation protocol for your structures[/grey50]\n")

        self.console.print("[bold cyan]Protocol Sources:[/bold cyan]")
        self.console.print("  1. Use a predefined protocol  - Standard multi-step equilibration protocols")
        self.console.print("  2. Load a saved protocol      - Reuse a protocol from a previous project")
        self.console.print("  3. Build a custom protocol    - Select steps from the template catalog")
        self.console.print("  4. Create a step from scratch - Configure parameters via guided wizard\n")

        self.console.print("[bold cyan]Protocol Customization:[/bold cyan]")
        self.console.print("  - Preview templates before selecting")
        self.console.print("  - Compare templates side-by-side")
        self.console.print("  - Select specific steps from predefined protocols")
        self.console.print("  - Edit MDIN parameters directly (opens text editor)")
        self.console.print("  - Add, remove, or reorder steps\n")

        self.console.print("[bold cyan]Multi-Structure Options:[/bold cyan]")
        self.console.print("  - Apply to All: Use same protocol for remaining structures")
        self.console.print("  - Individual: Configure each structure separately")
        self.console.print("  - Skip: Reuse previous structure's protocol\n")

    def _process_structure_pair(self, structure_pair: Dict) -> Optional[StructureWorkflowAssignment]:
        """Process a single structure pair through protocol selection."""
        structure_name = structure_pair['name']

        while True:
            self.console.print(self.layout.step_header(step_num=1))
            self.console.print("[grey50]  A typical MD protocol proceeds through phases: minimization (relax clashes),[/grey50]")
            self.console.print("[grey50]  heating (bring to target temperature), equilibration (stabilize density/pressure),[/grey50]")
            self.console.print("[grey50]  and production (collect data for analysis).[/grey50]")
            self.console.print(f"\n[grey50]Structure: {structure_name} ({self.current_structure_index + 1} of {self.total_structures})[/grey50]\n")

            self.console.print("[bold]Protocol Source:[/bold]\n")

            self.console.print("   1. Use a predefined protocol     [grey50]— standard multi-step protocols[/grey50]")
            self.console.print("   2. Load a saved protocol         [grey50]— reuse a protocol from a previous project[/grey50]")
            self.console.print("   3. Build a custom protocol       [grey50]— assemble steps from the template catalog[/grey50]")
            self.console.print("   4. Create a step from scratch    [grey50]— launch the parameter wizard[/grey50]")

            if self.current_structure_index > 0:
                self.console.print("\n   5. Skip (use previous protocol)")

            self.console.print(f"\n[bold cyan]Navigation:[/bold cyan]")
            self.console.print("   [cyan]h[/cyan]. Show overview")
            self.console.print("   [yellow]b[/yellow]. Back to Structure Files")
            self.console.print("   [grey50]x[/grey50]. Exit setup")
            self.console.print()

            valid_choices = ["1", "2", "3", "4", "h", "b", "x"]
            if self.current_structure_index > 0:
                valid_choices.append("5")

            choice = prompt_with_context(
                self.processor,
                "Enter choice",
                choices=valid_choices,
                default="1",
                module="MD Manager - Workflow Assignment",
                description="Protocol source for structure",
            )

            if choice == "1":
                config = self._handle_predefined_protocol(structure_pair)
                if config:
                    return self._create_protocol_assignment(structure_pair, config)

            elif choice == "2":
                config = self._handle_load_saved_protocol(structure_pair)
                if config:
                    return self._create_protocol_assignment(structure_pair, config)

            elif choice == "3":
                config = self._handle_custom_protocol(structure_pair)
                if config:
                    return self._create_protocol_assignment(structure_pair, config)

            elif choice == "4":
                config = self._handle_wizard_creation(structure_pair)
                if config:
                    return self._create_protocol_assignment(structure_pair, config)

            elif choice == "5" and self.current_structure_index > 0:
                if self.structure_assignments:
                    previous = self.structure_assignments[-1]
                    return self._clone_assignment_for_structure(previous, structure_pair)

            elif choice == "h":
                self.display_step_summary()

            elif choice == "b":
                return None

            elif choice == "x":
                raise KeyboardInterrupt("User requested exit")

    def _handle_predefined_protocol(self, structure_pair: Dict) -> Optional[Dict]:
        """Handle predefined protocol selection."""
        return self.predefined_selector.select_workflow(structure_pair)

    def _handle_custom_protocol(self, structure_pair: Dict) -> Optional[Dict]:
        """Handle custom protocol building from template catalog."""
        return self.catalog_browser.build_protocol(structure_pair)

    def _handle_wizard_creation(self, structure_pair: Dict) -> Optional[Dict]:
        """Handle new step creation via wizard."""
        return self.template_wizard.create_template(structure_pair)

    def _handle_load_saved_protocol(self, structure_pair: Dict) -> Optional[Dict]:
        """Handle loading a saved protocol."""
        return self.saved_protocol_selector.select_protocol(structure_pair)

    def _create_protocol_assignment(self, structure_pair: Dict, config: Dict) -> StructureWorkflowAssignment:
        """Create protocol and assignment from configuration."""
        structure_name = structure_pair['name']

        default_name = f"Protocol_{structure_name}_{datetime.now().strftime('%Y%m%d_%H%M')}"
        protocol_name = prompt_with_context(
            self.processor,
            "Protocol name",
            default=default_name,
            module="MD Manager - Workflow Assignment",
            description="Name for new protocol",
        )

        # Pass the actual coordinate filename so first step shows it instead of 'inpcrd'
        rst7_path = structure_pair.get('rst7')
        if rst7_path:
            from pathlib import Path
            config['initial_coords'] = Path(str(rst7_path)).name

        edit_requested = config.get('edit_requested', False)
        protocol = self.protocol_creator.create_protocol(
            config, protocol_name, structure_name, edit_requested=edit_requested
        )

        # Check for DISANG integration
        disang_config = self._check_disang_integration(protocol)
        if disang_config:
            protocol = self._apply_disang_integration(protocol, disang_config)

        # Save protocol
        protocol_id = self.user_data_manager.save_custom_workflow(protocol)

        return StructureWorkflowAssignment(
            structure_name=structure_name,
            structure_pair=structure_pair,
            custom_workflow_id=protocol_id,
            custom_workflow_name=protocol_name,
            source_type=config['source_type'],
            parent_workflow_id=config.get('parent_workflow_id')
        )

    def _clone_assignment_for_structure(self, source_assignment: StructureWorkflowAssignment,
                                        new_structure_pair: Dict) -> StructureWorkflowAssignment:
        """Clone a protocol assignment for a new structure pair."""
        structure_name = new_structure_pair['name']

        source_protocol = self.user_data_manager.load_custom_workflow(source_assignment.custom_workflow_id)

        new_name = f"{source_assignment.custom_workflow_name}_{structure_name}"

        # Deep copy — no external template files to clone
        cloned = copy.deepcopy(source_protocol)
        cloned['id'] = str(uuid.uuid4())
        cloned['name'] = new_name
        cloned['structure_name'] = structure_name
        cloned['created_date'] = datetime.now().isoformat()

        cloned_id = self.user_data_manager.save_custom_workflow(cloned)

        self.console.print(f"[green]Cloned protocol for {structure_name}[/green]")

        return StructureWorkflowAssignment(
            structure_name=structure_name,
            structure_pair=new_structure_pair,
            custom_workflow_id=cloned_id,
            custom_workflow_name=new_name,
            source_type=source_assignment.source_type,
            parent_workflow_id=source_assignment.parent_workflow_id
        )

    def _apply_assignment_to_remaining(self, remaining_structures: List[Dict],
                                       first_assignment: StructureWorkflowAssignment) -> None:
        """Clone the first structure's protocol onto every remaining structure."""
        for structure_pair in remaining_structures:
            cloned = self._clone_assignment_for_structure(first_assignment, structure_pair)
            self.structure_assignments.append(cloned)

    def _offer_apply_to_all(self, remaining_structures: List[Dict],
                            first_assignment: StructureWorkflowAssignment) -> bool:
        """Offer to apply first structure's protocol to remaining structures."""
        if not remaining_structures:
            return False

        self.console.print(f"\n[bold cyan]Apply to Remaining Structures?[/bold cyan]")
        self.console.print(f"You have {len(remaining_structures)} more structures to configure.")
        self.console.print(f"Current protocol: [bold]{first_assignment.custom_workflow_name}[/bold]")

        apply_all = confirm_with_context(
            self.processor,
            "Apply this protocol to all remaining structures?",
            default=True,
            module="MD Manager - Workflow Assignment",
            description="Apply protocol to all remaining structures",
        )

        if apply_all:
            self._apply_assignment_to_remaining(remaining_structures, first_assignment)
            return True

        return False

    def _check_disang_integration(self, protocol: Dict) -> Optional[Dict]:
        """Check for DISANG restraints and configure integration."""
        disang_files = []
        if self.workspace:
            workspace_files = self.workspace.get('files', [])
            disang_files = [f for f in workspace_files if f.lower().endswith('.disang')]

        if not disang_files:
            return None

        self.console.print(f"\n[yellow]DISANG restraint files detected: {', '.join(disang_files)}[/yellow]")

        use_disang = confirm_with_context(
            self.processor,
            "Integrate DISANG restraints into protocol?",
            default=True,
            module="MD Manager - Workflow Assignment",
            description="Integrate DISANG restraints into protocol",
        )
        if not use_disang:
            return None

        protocol_steps = protocol.get('steps', [])
        if not protocol_steps:
            return None

        self.console.print("\n[bold]Select protocol steps for DISANG integration:[/bold]")

        step_selections = {}
        for i, step in enumerate(protocol_steps):
            step_name = step.get('name', f"Step {i+1}")
            step_type = step.get('type', 'unknown')
            include = confirm_with_context(
                self.processor,
                f"Include DISANG in {step_name} ({step_type})?",
                default=False,
                module="MD Manager - Workflow Assignment",
                description=f"Include DISANG restraints in step {step_name}",
            )
            step_selections[i] = include

        return {
            'disang_files': disang_files,
            'step_selections': step_selections
        }

    def _apply_disang_integration(self, protocol: Dict, disang_config: Dict) -> Dict:
        """Apply DISANG integration to selected protocol steps."""
        protocol_steps = protocol.get('steps', [])

        for step_idx, include_disang in disang_config['step_selections'].items():
            if include_disang and step_idx < len(protocol_steps):
                step = protocol_steps[step_idx]

                # Store nmropt override
                overrides = step.get('parameter_overrides', {})
                overrides['nmropt'] = 1
                step['parameter_overrides'] = overrides

                step['nmr_section'] = step.get('nmr_section', '') + "\n&wt type='END' /"

                self.console.print(f"[green]DISANG integration applied to {step.get('name', 'step')}[/green]")

        return protocol

    def _finalize_step1(self) -> List[StructureWorkflowAssignment]:
        """Finalize Step 1 and prepare for Step 2."""
        if not self.structure_assignments:
            self.console.print("[red]No protocol assignments created[/red]")
            return []

        self.console.print(f"\n[bold green]Step 2 Complete - {len(self.structure_assignments)} protocol(s) configured[/bold green]")

        table = Table(title="Protocol Assignments")
        table.add_column("Structure", style="cyan")
        table.add_column("Protocol", style="yellow")
        table.add_column("Source", style="green")

        for assignment in self.structure_assignments:
            table.add_row(
                assignment.structure_name,
                assignment.custom_workflow_name,
                assignment.source_type
            )

        self.console.print(table)

        if self.workspace:
            workspace_assignments = []
            for assignment in self.structure_assignments:
                workspace_assignments.append({
                    'structure_name': assignment.structure_name,
                    'structure_pair': assignment.structure_pair,
                    'custom_workflow_id': assignment.custom_workflow_id,
                    'custom_workflow_name': assignment.custom_workflow_name,
                    'source_type': assignment.source_type,
                    'parent_workflow_id': assignment.parent_workflow_id
                })
            self.workspace.set('workflow_assignments', workspace_assignments)

        return self.structure_assignments
