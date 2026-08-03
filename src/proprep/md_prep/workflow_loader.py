"""
JSON-Based Workflow Preset System

Loads workflow presets from JSON files and templates from organized directory structure.
Replaces hardcoded workflow definitions with data-driven approach.
"""

import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, asdict
from rich.console import Console
from rich.table import Table

from proprep.utils.paths import get_package_dir
from proprep.utils.prompts import prompt_with_context


@dataclass
class WorkflowStep:
    """Represents a single step in a workflow."""
    id: str
    name: str
    type: str  # minimization, heating, equilibration, production
    template: str  # relative path to template file
    description: str
    parameter_overrides: Dict[str, Any]
    dependencies: List[str]
    input_coord: str
    output_coord: str


@dataclass 
class WorkflowPreset:
    """Represents a complete workflow preset."""
    name: str
    description: str
    version: str
    author: str
    category: str
    steps: List[WorkflowStep]
    
    
class WorkflowLoader:
    """Loads workflow presets from JSON files and templates from organized structure."""
    
    def __init__(self, console: Optional[Console] = None):
        self.console = console or Console()

        # Initialize user data manager
        from .user_data_manager import UserDataManager
        self.user_data_manager = UserDataManager(console=self.console)

        self.package_dir = get_package_dir()

        self.template_dir = self.package_dir / "md_templates"
        self.workflow_dir = self.package_dir / "md_workflows"

        # Ensure directories exist
        if not self.template_dir.exists():
            raise FileNotFoundError(f"Template directory not found: {self.template_dir}")
        if not self.workflow_dir.exists():
            raise FileNotFoundError(f"Workflow directory not found: {self.workflow_dir}")
            
    def get_available_workflows(self) -> Dict[str, WorkflowPreset]:
        """Load all available workflow presets using user data manager."""
        workflows = {}

        # Get workflows from user data manager
        workflow_list = self.user_data_manager.list_workflows()

        for workflow_id, workflow_meta in workflow_list.items():
            try:
                if workflow_meta["source"] == "builtin":
                    # Load builtin workflow
                    workflow_path = workflow_meta["workflow_path"]
                    workflow = self._load_workflow_from_json(self.workflow_dir / workflow_path)
                else:
                    # Load user workflow
                    workflow_content, _ = self.user_data_manager.get_workflow_content(workflow_id)
                    workflow = self._load_workflow_from_dict(workflow_content, workflow_id)

                workflows[workflow_id] = workflow

            except Exception as e:
                # Silently skip workflows that can't be loaded (e.g., files deleted but still in metadata)
                # This is expected behavior when workflows are cleaned up manually
                continue

        return workflows
        
    def _load_workflow_from_json(self, json_path: Path) -> WorkflowPreset:
        """Load a workflow preset from a JSON file."""
        with open(json_path, 'r') as f:
            data = json.load(f)
            
        # Parse steps
        steps = []
        for step_data in data.get("steps", []):
            # Load template content
            template_path = self.template_dir / step_data["template"]
            if not template_path.exists():
                raise FileNotFoundError(f"Template not found: {template_path}")
                
            with open(template_path, 'r') as f:
                template_content = f.read()
                
            step = WorkflowStep(
                id=step_data["id"],
                name=step_data["name"], 
                type=step_data["type"],
                template=step_data["template"],
                description=step_data["description"],
                parameter_overrides=step_data.get("parameter_overrides", {}),
                dependencies=step_data.get("dependencies", []),
                input_coord=step_data["input_coord"],
                output_coord=step_data["output_coord"]
            )
            steps.append(step)
            
        workflow = WorkflowPreset(
            name=data["name"],
            description=data["description"], 
            version=data["version"],
            author=data["author"],
            category=data["category"],
            steps=steps
        )
        
        return workflow
        
    def _load_workflow_from_dict(self, workflow_data: Dict, workflow_id: str) -> WorkflowPreset:
        """Load a workflow preset from a dictionary (for user workflows)."""
        # Parse steps
        steps = []
        for i, step_data in enumerate(workflow_data.get("steps", []), 1):
            try:
                # Use template_ref (new-style) with fallback to custom_template_id (legacy)
                template = step_data.get("template_ref", step_data.get("custom_template_id", ""))

                step = WorkflowStep(
                    id=step_data["id"],
                    name=step_data["name"],
                    type=step_data["type"],
                    template=template,
                    description=step_data.get("description", ""),
                    parameter_overrides=step_data.get("parameter_overrides", {}),
                    dependencies=step_data.get("dependencies", []),
                    input_coord=step_data.get("input_coord", "inpcrd"),
                    output_coord=step_data.get("output_coord", f"step_{i}.rst")
                )
                steps.append(step)
            except KeyError as e:
                raise KeyError(
                    f"Step {i} in workflow '{workflow_data.get('name', workflow_id)}' is missing required field: {e}. "
                    f"Required fields: id, name, type, template_ref (or custom_template_id)"
                )
            
        workflow = WorkflowPreset(
            name=workflow_data["name"],
            description=workflow_data["description"],
            version=workflow_data["version"],
            author=workflow_data["author"],
            category=workflow_data.get("category", "custom"),
            steps=steps
        )
        
        return workflow
        
    def get_template_content(self, template_path: str) -> str:
        """Get the content of a template file."""
        full_path = self.template_dir / template_path
        if not full_path.exists():
            raise FileNotFoundError(f"Template not found: {full_path}")
            
        with open(full_path, 'r') as f:
            return f.read()
            
    def apply_parameter_overrides(self, template_content: str, overrides: Dict[str, Any]) -> str:
        """Apply parameter overrides to template content."""
        if not overrides:
            return template_content
            
        # Simple parameter replacement (could be enhanced with more sophisticated parsing)
        content = template_content
        for param, value in overrides.items():
            # Look for parameter in format "param = old_value,"
            import re
            pattern = rf"(\s*{param}\s*=\s*)[^,\n]+([,\n])"
            replacement = rf"\g<1>{value}\g<2>"
            content = re.sub(pattern, replacement, content)
            
        return content
        
    def display_workflow_menu(self) -> Optional[str]:
        """Display workflow selection menu and return selected workflow key."""
        # Use user data manager to get workflow list with proper categorization
        workflow_list = self.user_data_manager.list_workflows()
        
        if not workflow_list:
            self.console.print("[yellow]No workflow presets available[/yellow]")
            return None
            
        self.console.print("\n[bold cyan]===== Available Workflow Presets =====[/bold cyan]")
        
        choices = {}
        choice_num = 1
        
        # Built-in workflows
        builtin_workflows = {k: v for k, v in workflow_list.items() if v["source"] == "builtin"}
        if builtin_workflows:
            self.console.print("\n[bold]Built-in Workflows:[/bold]")
            for workflow_id, metadata in builtin_workflows.items():
                self.console.print(f"  {choice_num}. {metadata['name']} - {metadata['description']}")
                choices[str(choice_num)] = workflow_id
                choice_num += 1
                
        # User custom workflows  
        custom_workflows = {k: v for k, v in workflow_list.items() if v["source"] == "custom"}
        if custom_workflows:
            self.console.print("\n[bold]Custom Workflows:[/bold]")
            for workflow_id, metadata in custom_workflows.items():
                self.console.print(f"  {choice_num}. {metadata['name']} - {metadata['description']}")
                choices[str(choice_num)] = workflow_id
                choice_num += 1
                
        # Modified workflows
        modified_workflows = {k: v for k, v in workflow_list.items() if v["source"].startswith("modified")}
        if modified_workflows:
            self.console.print("\n[bold]Modified Workflows:[/bold]")
            for workflow_id, metadata in modified_workflows.items():
                self.console.print(f"  {choice_num}. {metadata['name']} - {metadata['description']}")
                choices[str(choice_num)] = workflow_id
                choice_num += 1
                
        # Add management options
        self.console.print(f"\n[bold]Management:[/bold]")
        self.console.print(f"  {choice_num}. + Create new workflow")
        choices[str(choice_num)] = "create_new"
        choice_num += 1
        
        self.console.print(f"  {choice_num}. ✏ Manage templates & workflows")
        choices[str(choice_num)] = "manage"
        choice_num += 1
                
        # Add return option
        self.console.print(f"  {choice_num}. ← Return to AMBER Input Generator")
        choices[str(choice_num)] = None  # Return option
        
        choice = prompt_with_context(None, "Select option", choices=list(choices.keys()), default="1")
        return choices.get(choice)
        
    def display_workflow_overview(self, workflow: WorkflowPreset):
        """Display detailed overview of a workflow."""
        self.console.print(f"\n[bold cyan]===== {workflow.name} =====[/bold cyan]")
        self.console.print(f"[bold]Description:[/bold] {workflow.description}")
        self.console.print(f"[bold]Category:[/bold] {workflow.category}")
        self.console.print(f"[bold]Version:[/bold] {workflow.version}")
        self.console.print(f"[bold]Author:[/bold] {workflow.author}")
        self.console.print(f"[bold]Steps:[/bold] {len(workflow.steps)}")
        
        # Create overview table
        table = Table(title="Workflow Steps")
        table.add_column("Step", style="cyan", width=15)
        table.add_column("Type", style="green", width=15) 
        table.add_column("Description", style="white", width=50)
        
        for step in workflow.steps:
            table.add_row(step.id, step.type.title(), step.name)
            
        self.console.print(table)
        
    def save_workflow_to_json(self, workflow: WorkflowPreset, output_path: Path):
        """Save workflow preset to JSON file."""
        # Convert to serializable format
        workflow_dict = asdict(workflow)
        
        with open(output_path, 'w') as f:
            json.dump(workflow_dict, f, indent=2)
            
    def generate_workflow_sequence_json(self, workflow: WorkflowPreset, output_dir: Path) -> Path:
        """Generate workflow sequence JSON for Workflow Manager integration."""
        sequence_data = {
            "workflow_name": workflow.name,
            "description": workflow.description,
            "version": workflow.version,
            "steps": []
        }
        
        for step in workflow.steps:
            step_data = {
                "id": step.id,
                "name": step.name,
                "type": step.type,
                "input_file": f"{step.id}.mdin",
                "input_coord": step.input_coord,
                "output_coord": step.output_coord,
                "dependencies": step.dependencies
            }
            sequence_data["steps"].append(step_data)
            
        sequence_file = output_dir / f"{workflow.name.lower().replace(' ', '_')}_sequence.json"
        with open(sequence_file, 'w') as f:
            json.dump(sequence_data, f, indent=2)
            
        return sequence_file