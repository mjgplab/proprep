"""
Protocol Creator

Creates simulation protocols from various sources (predefined workflows,
template catalog selections, or wizard-created steps). Protocols reference
builtin templates directly and carry edits inline — no custom template copies
are created during selection.
"""

import uuid
from datetime import datetime
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, asdict, field
import copy

from rich.console import Console

from proprep.utils.prompts import (
    confirm_with_context,
)

from .user_data_manager import UserDataManager
from .workflow_editor import WorkflowEditor


@dataclass
class ProtocolStep:
    """Represents a step in a simulation protocol."""
    id: str
    name: str
    type: str  # minimization, heating, equilibration, production
    template_ref: str  # builtin path (e.g., "builtin/york_lab_tutorial/00_initial_minimization.mdin")
    description: str
    dependencies: List[str]
    input_coord: str
    output_coord: str
    parameter_overrides: Dict[str, Any] = field(default_factory=dict)
    mdin_content_override: Optional[str] = None
    nmr_section: str = ""


@dataclass
class Protocol:
    """Represents a complete simulation protocol."""
    id: str
    name: str
    description: str
    structure_name: str
    source_type: str  # "predefined", "catalog", "wizard"
    parent_workflow_id: Optional[str]
    steps: List[ProtocolStep]
    created_date: str
    modified_date: str
    version: str = "1.0"
    author: str = "user"
    initial_coords: str = "inpcrd"  # Actual coordinate filename for first step display


class ProtocolCreator:
    """
    Creates simulation protocols from various sources.

    Protocols reference builtin templates directly via template_ref and carry
    parameter edits inline. No custom template files are created during
    protocol assembly — MDIN files are materialized at write time (Step 5).
    """

    def __init__(self, console: Console, user_data_manager: UserDataManager, processor=None):
        self.console = console
        self.processor = processor
        self.user_data_manager = user_data_manager
        self.workflow_editor = WorkflowEditor(console, user_data_manager, processor=processor)

    def create_protocol(self, config: Dict, protocol_name: str, structure_name: str,
                        edit_requested: bool = False) -> Dict:
        """
        Create protocol from configuration returned by a source selector.

        Args:
            config: Configuration from workflow source selectors
            protocol_name: Name for the protocol
            structure_name: Target structure name
            edit_requested: Whether user requested to edit the protocol

        Returns:
            Protocol dictionary ready for saving
        """
        source_type = config['source_type']

        if source_type == 'predefined':
            return self._create_from_predefined(config, protocol_name, structure_name, edit_requested)
        elif source_type in ('standalone', 'custom_workflow', 'catalog'):
            return self._create_from_catalog_selection(config, protocol_name, structure_name, edit_requested)
        elif source_type == 'new':
            return self._create_from_wizard(config, protocol_name, structure_name)
        else:
            raise ValueError(f"Unknown source type: {source_type}")

    def _create_from_predefined(self, config: Dict, protocol_name: str,
                                structure_name: str, edit_requested: bool = False) -> Dict:
        """Create protocol from a predefined workflow."""
        original_workflow = config['original_workflow']
        selected_steps = config['selected_steps']
        parent_workflow_id = config['parent_workflow_id']
        initial_coords = config.get('initial_coords', 'inpcrd')

        self.console.print(f"\n[bold cyan]Creating protocol from '{original_workflow.name}'[/bold cyan]")

        protocol_steps = []
        for i, original_step in enumerate(selected_steps):
            default_input = initial_coords if i == 0 else getattr(original_step, 'input_coord', f"step_{i}.rst")
            step = ProtocolStep(
                id=f"step_{i+1}",
                name=original_step.name,
                type=original_step.type,
                template_ref=getattr(original_step, 'template', ''),
                description=getattr(original_step, 'description', ''),
                dependencies=self._calculate_dependencies(i),
                input_coord=default_input,
                output_coord=getattr(original_step, 'output_coord', f"step_{i+1}.rst"),
                parameter_overrides=getattr(original_step, 'parameter_overrides', {}),
                nmr_section=getattr(original_step, 'nmr_section', '')
            )
            protocol_steps.append(step)

        protocol = Protocol(
            id=str(uuid.uuid4()),
            name=protocol_name,
            description=f"Based on {original_workflow.name}",
            structure_name=structure_name,
            source_type='predefined',
            parent_workflow_id=parent_workflow_id,
            steps=protocol_steps,
            created_date=datetime.now().isoformat(),
            modified_date=datetime.now().isoformat(),
            initial_coords=initial_coords,
        )

        if edit_requested:
            protocol = self._edit_protocol(protocol)

        return asdict(protocol)

    def _create_from_catalog_selection(self, config: Dict, protocol_name: str,
                                       structure_name: str, edit_requested: bool = False) -> Dict:
        """Create protocol from template catalog selections."""
        selected_templates = config['selected_templates']
        initial_coords = config.get('initial_coords', 'inpcrd')

        self.console.print(f"\n[bold cyan]Creating protocol from {len(selected_templates)} step(s)[/bold cyan]")

        protocol_steps = []
        for i, template_info in enumerate(selected_templates):
            template_data = template_info['data']
            template_type = template_info['type']
            # Use template_ref (builtin path) directly — no custom copy
            template_ref = template_data.get('template_path', template_info.get('id', ''))

            step = ProtocolStep(
                id=f"step_{i+1}",
                name=template_data.get('name', f"Step {i+1}"),
                type=template_type,
                template_ref=template_ref,
                description=template_data.get('description', ''),
                dependencies=self._calculate_dependencies(i),
                input_coord=initial_coords if i == 0 else f"step_{i}.rst",
                output_coord=f"step_{i+1}.rst",
            )
            protocol_steps.append(step)

        protocol = Protocol(
            id=str(uuid.uuid4()),
            name=protocol_name,
            description="Custom protocol from template catalog",
            structure_name=structure_name,
            source_type='catalog',
            parent_workflow_id=None,
            steps=protocol_steps,
            created_date=datetime.now().isoformat(),
            modified_date=datetime.now().isoformat(),
            initial_coords=initial_coords,
        )

        if edit_requested:
            protocol = self._edit_protocol(protocol)

        return asdict(protocol)

    def _create_from_wizard(self, config: Dict, protocol_name: str,
                            structure_name: str) -> Dict:
        """Create protocol from a wizard-created step."""
        new_template = config['new_template']
        template_data = new_template['data']
        template_type = new_template['type']
        initial_coords = config.get('initial_coords', 'inpcrd')

        self.console.print(f"\n[bold cyan]Creating protocol from new {template_type} step[/bold cyan]")

        step = ProtocolStep(
            id="step_1",
            name=template_data.get('name', f"Custom {template_type}"),
            type=template_type,
            template_ref="",  # No builtin reference — content generated from overrides
            description=template_data.get('description', ''),
            dependencies=[],
            input_coord=initial_coords,
            output_coord='step_1.rst',
            parameter_overrides=template_data.get('config', {}),
            nmr_section=template_data.get('nmr_section', '')
        )

        protocol = Protocol(
            id=str(uuid.uuid4()),
            name=protocol_name,
            description=f"Protocol from new {template_type} step",
            structure_name=structure_name,
            source_type='wizard',
            parent_workflow_id=None,
            steps=[step],
            created_date=datetime.now().isoformat(),
            modified_date=datetime.now().isoformat(),
            initial_coords=initial_coords,
        )

        if confirm_with_context(
            self.processor,
            "Add more steps to this protocol?",
            default=True,
            module="MD Manager - Protocol Creator",
            description="Add more steps to protocol",
        ):
            protocol = self._edit_protocol(protocol)

        return asdict(protocol)

    def _calculate_dependencies(self, step_index: int) -> List[str]:
        """Calculate dependencies for a protocol step (sequential)."""
        if step_index == 0:
            return []
        return [f"step_{step_index}"]

    def _edit_protocol(self, protocol: Protocol) -> Protocol:
        """Edit protocol using WorkflowEditor."""
        self.console.print(f"\n[bold green]Editing Protocol: {protocol.name}[/bold green]")

        workflow_dict = asdict(protocol)
        edited = self.workflow_editor.edit_workflow(workflow_dict)

        if edited:
            edited['modified_date'] = datetime.now().isoformat()
            steps = []
            for step_dict in edited.get('steps', []):
                step = ProtocolStep(**step_dict)
                steps.append(step)

            protocol.steps = steps
            protocol.modified_date = edited['modified_date']
            protocol.name = edited.get('name', protocol.name)
            protocol.description = edited.get('description', protocol.description)

        return protocol

    def _edit_protocol_dict(self, protocol_dict: Dict) -> Dict:
        """Edit a protocol dict (for saved protocols that aren't Protocol objects)."""
        edited = self.workflow_editor.edit_workflow(protocol_dict)
        if edited:
            edited['modified_date'] = datetime.now().isoformat()
            return edited
        return protocol_dict


class ProtocolCloner:
    """Handles cloning protocols for multiple structures."""

    def __init__(self, console: Console, user_data_manager: UserDataManager):
        self.console = console
        self.user_data_manager = user_data_manager

    def clone_protocol_for_structure(self, source_protocol: Dict,
                                     target_structure_name: str) -> Dict:
        """Clone a protocol for a different structure.

        Since protocols carry all edits inline (no external template files),
        cloning is just a deep copy with updated metadata.
        """
        cloned = copy.deepcopy(source_protocol)

        original_name = source_protocol.get('name', 'Protocol')
        cloned['id'] = str(uuid.uuid4())
        cloned['name'] = f"{original_name}_{target_structure_name}"
        cloned['structure_name'] = target_structure_name
        cloned['created_date'] = datetime.now().isoformat()
        cloned['modified_date'] = datetime.now().isoformat()

        return cloned


# Backward-compatible aliases
CustomWorkflowStep = ProtocolStep
CustomWorkflow = Protocol
CustomWorkflowCreator = ProtocolCreator
WorkflowCloner = ProtocolCloner
