"""
Template and Workflow Management Interface

Provides user-friendly interface for creating, modifying, and managing templates and workflows.
"""

import tempfile
import subprocess
import os
from pathlib import Path
from typing import Dict, List, Optional, Any
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.text import Text

from proprep.utils.prompts import (
    prompt_with_context,
    confirm_with_context,
    int_prompt_with_context,
)

from .user_data_manager import UserDataManager
from .workflow_loader import WorkflowLoader


class TemplateWorkflowManager:
    """Interactive management interface for templates and workflows."""
    
    def __init__(self, console: Optional[Console] = None, processor=None):
        self.console = console or Console()
        self.processor = processor
        self.user_data_manager = UserDataManager(console=self.console)
        self.workflow_loader = WorkflowLoader(console=self.console)
        
    def show_management_menu(self):
        """Display main management menu."""
        while True:
            self.console.print("\n[bold cyan]===== Template & Workflow Management =====[/bold cyan]")
            
            choices = {
                "1": "template_menu",
                "2": "workflow_menu", 
                "3": "import_export",
                "4": "back"
            }
            
            self.console.print("\n[bold]Options:[/bold]")
            self.console.print("  1. Manage Templates")
            self.console.print("  2. Manage Workflows")
            self.console.print("  3. Import/Export Content")
            self.console.print("  4. ← Return to workflow selection")
            
            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=list(choices.keys()),
                default="1",
                module="Template & Workflow Manager",
                description="Select menu option",
            )
            action = choices.get(choice)
            
            if action == "back":
                break
            elif action == "template_menu":
                self._show_template_menu()
            elif action == "workflow_menu":
                self._show_workflow_menu()
            elif action == "import_export":
                self._show_import_export_menu()
                
    def _show_template_menu(self):
        """Display template management menu."""
        while True:
            self.console.print("\n[bold cyan]===== Template Management =====[/bold cyan]")
            
            # Show template statistics
            templates = self.user_data_manager.list_templates()
            builtin_count = len([t for t in templates.values() if t["source"] == "builtin"])
            user_count = len([t for t in templates.values() if t["source"] != "builtin"])
            
            self.console.print(f"[grey50]Built-in templates: {builtin_count} | User templates: {user_count}[/grey50]\n")
            
            choices = {
                "1": "list_templates",
                "2": "create_template",
                "3": "edit_template",
                "4": "delete_template",
                "5": "back"
            }
            
            self.console.print("[bold]Options:[/bold]")
            self.console.print("  1. List all templates")
            self.console.print("  2. Create new template")  
            self.console.print("  3. Edit existing template")
            self.console.print("  4. Delete user template")
            self.console.print("  5. ← Back to main menu")
            
            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=list(choices.keys()),
                default="1",
                module="Template & Workflow Manager",
                description="Select menu option",
            )
            action = choices.get(choice)
            
            if action == "back":
                break
            elif action == "list_templates":
                self._list_templates()
            elif action == "create_template":
                self._create_template_wizard()
            elif action == "edit_template":
                self._edit_template_wizard()
            elif action == "delete_template":
                self._delete_template_wizard()
                
    def _show_workflow_menu(self):
        """Display workflow management menu."""
        while True:
            self.console.print("\n[bold cyan]===== Workflow Management =====[/bold cyan]")
            
            # Show workflow statistics
            workflows = self.user_data_manager.list_workflows()
            builtin_count = len([w for w in workflows.values() if w["source"] == "builtin"])
            user_count = len([w for w in workflows.values() if w["source"] != "builtin"])
            
            self.console.print(f"[grey50]Built-in workflows: {builtin_count} | User workflows: {user_count}[/grey50]\n")
            
            choices = {
                "1": "list_workflows",
                "2": "create_workflow",
                "3": "edit_workflow", 
                "4": "delete_workflow",
                "5": "back"
            }
            
            self.console.print("[bold]Options:[/bold]")
            self.console.print("  1. List all workflows")
            self.console.print("  2. Create new workflow")
            self.console.print("  3. Edit existing workflow")
            self.console.print("  4. Delete user workflow")
            self.console.print("  5. ← Back to main menu")
            
            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=list(choices.keys()),
                default="1",
                module="Template & Workflow Manager",
                description="Select menu option",
            )
            action = choices.get(choice)
            
            if action == "back":
                break
            elif action == "list_workflows":
                self._list_workflows()
            elif action == "create_workflow":
                self._create_workflow_wizard()
            elif action == "edit_workflow":
                self._edit_workflow_wizard()
            elif action == "delete_workflow":
                self._delete_workflow_wizard()
                
    def _show_import_export_menu(self):
        """Display import/export menu."""
        while True:
            self.console.print("\n[bold cyan]===== Import/Export Content =====[/bold cyan]")
            
            choices = {
                "1": "export_template",
                "2": "export_workflow",
                "3": "import_content",
                "4": "back"
            }
            
            self.console.print("[bold]Options:[/bold]")
            self.console.print("  1. Export template")
            self.console.print("  2. Export workflow")
            self.console.print("  3. Import content")
            self.console.print("  4. ← Back to main menu")
            
            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=list(choices.keys()),
                default="1",
                module="Template & Workflow Manager",
                description="Select menu option",
            )
            action = choices.get(choice)
            
            if action == "back":
                break
            elif action == "export_template":
                self._export_content_wizard("template")
            elif action == "export_workflow":
                self._export_content_wizard("workflow") 
            elif action == "import_content":
                self._import_content_wizard()
                
    def _list_templates(self):
        """List all templates with metadata."""
        templates = self.user_data_manager.list_templates()
        
        if not templates:
            self.console.print("[yellow]No templates found[/yellow]")
            return
            
        # Group by source
        builtin_templates = {k: v for k, v in templates.items() if v["source"] == "builtin"}
        user_templates = {k: v for k, v in templates.items() if v["source"] != "builtin"}
        
        # Display builtin templates
        if builtin_templates:
            table = Table(title="Built-in Templates")
            table.add_column("Name", style="cyan")
            table.add_column("Type", style="green")
            table.add_column("Description", style="white")
            
            for template_id, metadata in builtin_templates.items():
                table.add_row(
                    metadata["name"],
                    metadata.get("simulation_type", "unknown").title(),
                    metadata["description"][:50] + "..." if len(metadata["description"]) > 50 else metadata["description"]
                )
                
            self.console.print(table)
            
        # Display user templates
        if user_templates:
            table = Table(title="User Templates")
            table.add_column("ID", style="grey50")
            table.add_column("Name", style="cyan")
            table.add_column("Type", style="green")
            table.add_column("Source", style="yellow")
            table.add_column("Author", style="blue")
            table.add_column("Created", style="grey50")
            
            for template_id, metadata in user_templates.items():
                created_date = metadata.get("created_date", "")[:10]  # Just the date part
                table.add_row(
                    template_id,
                    metadata["name"],
                    metadata.get("simulation_type", "unknown").title(),
                    metadata["source"].replace("_", " ").title(),
                    metadata.get("author", "Unknown"),
                    created_date
                )
                
            self.console.print(table)
            
        input("[grey50]Press Enter to continue...[/grey50]")
        
    def _list_workflows(self):
        """List all workflows with metadata."""
        workflows = self.user_data_manager.list_workflows()
        
        if not workflows:
            self.console.print("[yellow]No workflows found[/yellow]")
            return
            
        # Group by source
        builtin_workflows = {k: v for k, v in workflows.items() if v["source"] == "builtin"}
        user_workflows = {k: v for k, v in workflows.items() if v["source"] != "builtin"}
        
        # Display builtin workflows
        if builtin_workflows:
            table = Table(title="Built-in Workflows")
            table.add_column("Name", style="cyan")
            table.add_column("Category", style="green")
            table.add_column("Steps", style="yellow")
            table.add_column("Description", style="white")
            
            for workflow_id, metadata in builtin_workflows.items():
                table.add_row(
                    metadata["name"],
                    metadata.get("category", "unknown").title(),
                    str(metadata.get("step_count", 0)),
                    metadata["description"][:50] + "..." if len(metadata["description"]) > 50 else metadata["description"]
                )
                
            self.console.print(table)
            
        # Display user workflows
        if user_workflows:
            table = Table(title="User Workflows")
            table.add_column("ID", style="grey50")
            table.add_column("Name", style="cyan")
            table.add_column("Category", style="green")
            table.add_column("Steps", style="yellow")
            table.add_column("Source", style="yellow")
            table.add_column("Author", style="blue")
            table.add_column("Created", style="grey50")
            
            for workflow_id, metadata in user_workflows.items():
                created_date = metadata.get("created_date", "")[:10]  # Just the date part
                table.add_row(
                    workflow_id,
                    metadata["name"],
                    metadata.get("category", "unknown").title(),
                    str(metadata.get("step_count", 0)),
                    metadata["source"].replace("_", " ").title(),
                    metadata.get("author", "Unknown"),
                    created_date
                )
                
            self.console.print(table)
            
        input("[grey50]Press Enter to continue...[/grey50]")
        
    def _create_template_wizard(self):
        """Wizard for creating new templates."""
        self.console.print("\n[bold]Create New Template[/bold]")
        
        # Get basic info
        name = prompt_with_context(
            self.processor,
            "Template name",
            module="Template & Workflow Manager",
            description="New template name",
        )
        description = prompt_with_context(
            self.processor,
            "Description",
            module="Template & Workflow Manager",
            description="New template description",
        )
        
        # Get simulation type
        sim_types = ["minimization", "heating", "equilibration", "production"]
        self.console.print("\\nSimulation types:")
        for i, sim_type in enumerate(sim_types, 1):
            self.console.print(f"  {i}. {sim_type.title()}")
            
        sim_type_options = {str(i + 1): t for i, t in enumerate(sim_types)}
        type_choice = int_prompt_with_context(
            self.processor,
            "Select simulation type",
            module="Template & Workflow Manager",
            description="Simulation type",
            options_map=sim_type_options,
        )
        simulation_type = sim_types[type_choice - 1]
        
        # Get content
        self.console.print("\\n[bold]Template Content Options:[/bold]")
        self.console.print("  1. Create from scratch")
        self.console.print("  2. Copy from existing template")
        
        content_choice = prompt_with_context(
            self.processor,
            "Select option",
            choices=["1", "2"],
            default="1",
            module="Template & Workflow Manager",
            description="Template content source",
            options_map={"1": "Create from scratch", "2": "Copy existing template"},
        )
        
        if content_choice == "1":
            content = self._create_template_from_scratch(simulation_type)
        else:
            content = self._copy_existing_template()
            
        if not content:
            self.console.print("[yellow]Template creation cancelled[/yellow]")
            return
            
        # Get author
        author = prompt_with_context(
            self.processor,
            "Author name",
            default="User",
            module="Template & Workflow Manager",
            description="Template author name",
        )

        # Get tags
        tags_input = prompt_with_context(
            self.processor,
            "Tags (comma-separated, optional)",
            default="",
            module="Template & Workflow Manager",
            description="Template tags (comma-separated, optional)",
        )
        tags = [tag.strip() for tag in tags_input.split(",") if tag.strip()]
        
        # Create template
        try:
            template_id = self.user_data_manager.create_template(
                name=name,
                description=description,
                simulation_type=simulation_type,
                content=content,
                author=author,
                tags=tags
            )
            
            self.console.print(f"[green]✓ Template created successfully![/green]")
            self.console.print(f"Template ID: {template_id}")
            
        except Exception as e:
            self.console.print(f"[red]Error creating template: {e}[/red]")
            
    def _create_template_from_scratch(self, simulation_type: str) -> Optional[str]:
        """Create template content from scratch using editor."""
        # Create basic template structure
        basic_content = f"""# {simulation_type.title()} Template
# Created with ProPrep Template Manager

&cntrl
! Run Control
  imin = 0,                 ! Flag to run minimization (0=MD, 1=minimize)
  nstlim = 100000,          ! Number of MD steps
  dt = 0.001,               ! Time step (ps)

! Input/Output Control
  ntx = 1,                  ! Read coordinates: 1=start, 5=restart
  irest = 0,                ! Start (0) or restart (1) simulation
  ntpr = 1000,              ! Print energy every ntpr steps
  ntwx = 1000,              ! Write coordinates every ntwx steps
  ntwr = 1000,              ! Write restart every ntwr steps
  ioutfm = 1,               ! Binary NetCDF trajectory format

! Temperature Control
  ntt = 3,                  ! Langevin thermostat
  gamma_ln = 5.0,           ! Collision frequency for Langevin
  temp0 = 300.0,            ! Reference temperature (K)
  ig = -1,                  ! Random seed based on wallclock

! Miscellaneous
  cut = 12.0,               ! Non-bonded cutoff (Angstroms)

! Periodic Boundary Conditions
  ntb = 1,                  ! Constant volume PBC
  
! SHAKE constraints
  ntc = 2,                  ! SHAKE to constrain bonds involving hydrogen
  ntf = 2,                  ! Omit force evaluation for bonds constrained by SHAKE
/
"""
        
        # Write to temporary file for editing
        with tempfile.NamedTemporaryFile(mode='w', suffix='.mdin', delete=False) as tmp:
            tmp.write(basic_content)
            tmp_path = tmp.name
            
        try:
            # Get editor
            editor = os.environ.get('EDITOR', 'nano')
            
            self.console.print(f"[yellow]Opening template in {editor}...[/yellow]")
            self.console.print("[grey50]Save and close the editor when done.[/grey50]")
            input("Press Enter to continue...")
            
            # Open editor
            result = subprocess.run([editor, tmp_path])
            
            if result.returncode == 0:
                # Read modified content
                with open(tmp_path, 'r') as f:
                    content = f.read()
                    
                if content.strip() != basic_content.strip():
                    return content
                else:
                    self.console.print("[yellow]No changes made[/yellow]")
                    return None
            else:
                self.console.print("[red]Editor exited with error[/red]")
                return None
                
        finally:
            # Clean up
            os.unlink(tmp_path)
            
    def _copy_existing_template(self) -> Optional[str]:
        """Copy content from existing template."""
        templates = self.user_data_manager.list_templates()
        
        if not templates:
            self.console.print("[yellow]No templates available to copy[/yellow]")
            return None
            
        # Display templates for selection
        self.console.print("\\n[bold]Available Templates:[/bold]")
        choices = {}
        choice_num = 1
        
        for template_id, metadata in templates.items():
            self.console.print(f"  {choice_num}. {metadata['name']} ({metadata.get('simulation_type', 'unknown')})")
            choices[str(choice_num)] = template_id
            choice_num += 1
            
        template_choice = prompt_with_context(
            self.processor,
            "Select template to copy",
            choices=list(choices.keys()),
            module="Template & Workflow Manager",
            description="Select template to copy",
        )
        selected_id = choices.get(template_choice)
        
        try:
            content, metadata = self.user_data_manager.get_template_content(selected_id)
            self.console.print(f"[green]✓ Copied content from '{metadata['name']}'[/green]")
            return content
        except Exception as e:
            self.console.print(f"[red]Error copying template: {e}[/red]")
            return None
            
    def _edit_template_wizard(self):
        """Wizard for editing existing templates."""
        # Get user templates only (can't edit builtin)
        templates = self.user_data_manager.list_templates(show_builtin=False, show_user=True)
        user_templates = {k: v for k, v in templates.items() if v["source"] != "builtin"}
        
        if not user_templates:
            self.console.print("[yellow]No user templates available to edit[/yellow]")
            self.console.print("[grey50]You can only edit templates you created or modified[/grey50]")
            return
            
        # Display templates for selection
        self.console.print("\\n[bold]User Templates:[/bold]")
        choices = {}
        choice_num = 1
        
        for template_id, metadata in user_templates.items():
            self.console.print(f"  {choice_num}. {metadata['name']} ({metadata.get('simulation_type', 'unknown')})")
            choices[str(choice_num)] = template_id
            choice_num += 1
            
        template_choice = prompt_with_context(
            self.processor,
            "Select template to edit",
            choices=list(choices.keys()),
            module="Template & Workflow Manager",
            description="Select template to edit",
        )
        selected_id = choices.get(template_choice)
        
        try:
            content, metadata = self.user_data_manager.get_template_content(selected_id)
            
            # Edit in external editor
            edited_content = self._edit_content_in_editor(content, f"{metadata['name']}.mdin")
            
            if edited_content and edited_content != content:
                # Update existing template (this creates a new version)
                new_name = prompt_with_context(
                    self.processor,
                    "New template name",
                    default=metadata['name'],
                    module="Template & Workflow Manager",
                    description="New template name for edit",
                )
                new_description = prompt_with_context(
                    self.processor,
                    "New description",
                    default=metadata['description'],
                    module="Template & Workflow Manager",
                    description="New template description for edit",
                )
                
                new_id = self.user_data_manager.modify_template(
                    source_path=metadata.get('template_path', ''),
                    name=new_name,
                    description=new_description,
                    content=edited_content
                )
                
                self.console.print(f"[green]✓ Template updated! New ID: {new_id}[/green]")
            else:
                self.console.print("[yellow]No changes made[/yellow]")
                
        except Exception as e:
            self.console.print(f"[red]Error editing template: {e}[/red]")
            
    def _edit_content_in_editor(self, content: str, filename: str) -> Optional[str]:
        """Edit content in external editor."""
        with tempfile.NamedTemporaryFile(mode='w', suffix=f'_{filename}', delete=False) as tmp:
            tmp.write(content)
            tmp_path = tmp.name
            
        try:
            editor = os.environ.get('EDITOR', 'nano')
            
            self.console.print(f"[yellow]Opening {filename} in {editor}...[/yellow]")
            input("Press Enter to continue...")
            
            result = subprocess.run([editor, tmp_path])
            
            if result.returncode == 0:
                with open(tmp_path, 'r') as f:
                    return f.read()
            else:
                self.console.print("[red]Editor exited with error[/red]")
                return None
                
        finally:
            os.unlink(tmp_path)
            
    def _delete_template_wizard(self):
        """Wizard for deleting user templates."""
        templates = self.user_data_manager.list_templates(show_builtin=False, show_user=True)
        user_templates = {k: v for k, v in templates.items() if v["source"] != "builtin"}
        
        if not user_templates:
            self.console.print("[yellow]No user templates available to delete[/yellow]")
            return
            
        # Display templates for selection
        self.console.print("\\n[bold]User Templates:[/bold]")
        choices = {}
        choice_num = 1
        
        for template_id, metadata in user_templates.items():
            self.console.print(f"  {choice_num}. {metadata['name']} ({metadata.get('simulation_type', 'unknown')})")
            choices[str(choice_num)] = template_id
            choice_num += 1
            
        template_choice = prompt_with_context(
            self.processor,
            "Select template to delete",
            choices=list(choices.keys()),
            module="Template & Workflow Manager",
            description="Select template to delete",
        )
        selected_id = choices.get(template_choice)
        selected_metadata = user_templates[selected_id]
        
        # Confirm deletion
        if confirm_with_context(
            self.processor,
            f"Delete template '{selected_metadata['name']}'?",
            default=False,
            module="Template & Workflow Manager",
            description=f"Confirm delete template {selected_metadata['name']}",
        ):
            if self.user_data_manager.delete_user_content(selected_id):
                self.console.print("[green]✓ Template deleted successfully[/green]")
            else:
                self.console.print("[red]Failed to delete template[/red]")
                
    def _create_workflow_wizard(self):
        """Wizard for creating new workflows."""
        self.console.print("[yellow]Workflow creation wizard - Coming soon![/yellow]")
        self.console.print("[grey50]This feature will allow you to interactively build workflows by selecting templates and configuring steps.[/grey50]")
        
    def _edit_workflow_wizard(self):
        """Wizard for editing existing workflows."""
        self.console.print("[yellow]Workflow editing wizard - Coming soon![/yellow]")
        self.console.print("[grey50]This feature will allow you to modify existing workflows by adding, removing, or reordering steps.[/grey50]")
        
    def _delete_workflow_wizard(self):
        """Wizard for deleting user workflows."""
        workflows = self.user_data_manager.list_workflows(show_builtin=False, show_user=True)
        user_workflows = {k: v for k, v in workflows.items() if v["source"] != "builtin"}
        
        if not user_workflows:
            self.console.print("[yellow]No user workflows available to delete[/yellow]")
            return
            
        # Display workflows for selection
        self.console.print("\\n[bold]User Workflows:[/bold]")
        choices = {}
        choice_num = 1
        
        for workflow_id, metadata in user_workflows.items():
            self.console.print(f"  {choice_num}. {metadata['name']} ({metadata.get('category', 'unknown')})")
            choices[str(choice_num)] = workflow_id
            choice_num += 1
            
        workflow_choice = prompt_with_context(
            self.processor,
            "Select workflow to delete",
            choices=list(choices.keys()),
            module="Template & Workflow Manager",
            description="Select workflow to delete",
        )
        selected_id = choices.get(workflow_choice)
        selected_metadata = user_workflows[selected_id]
        
        # Confirm deletion
        if confirm_with_context(
            self.processor,
            f"Delete workflow '{selected_metadata['name']}'?",
            default=False,
            module="Template & Workflow Manager",
            description=f"Confirm delete workflow {selected_metadata['name']}",
        ):
            if self.user_data_manager.delete_user_content(selected_id):
                self.console.print("[green]✓ Workflow deleted successfully[/green]")
            else:
                self.console.print("[red]Failed to delete workflow[/red]")
                
    def _export_content_wizard(self, content_type: str):
        """Wizard for exporting templates or workflows."""
        if content_type == "template":
            content_dict = self.user_data_manager.list_templates(show_builtin=False, show_user=True)
            content_dict = {k: v for k, v in content_dict.items() if v["source"] != "builtin"}
        else:
            content_dict = self.user_data_manager.list_workflows(show_builtin=False, show_user=True)
            content_dict = {k: v for k, v in content_dict.items() if v["source"] != "builtin"}
            
        if not content_dict:
            self.console.print(f"[yellow]No user {content_type}s available to export[/yellow]")
            return
            
        # Display content for selection
        self.console.print(f"\\n[bold]User {content_type.title()}s:[/bold]")
        choices = {}
        choice_num = 1
        
        for content_id, metadata in content_dict.items():
            self.console.print(f"  {choice_num}. {metadata['name']}")
            choices[str(choice_num)] = content_id
            choice_num += 1
            
        content_choice = prompt_with_context(
            self.processor,
            f"Select {content_type} to export",
            choices=list(choices.keys()),
            module="Template & Workflow Manager",
            description=f"Select {content_type} to export",
        )
        selected_id = choices.get(content_choice)
        selected_metadata = content_dict[selected_id]
        
        # Get export path
        default_filename = f"{selected_metadata['name'].lower().replace(' ', '_')}.proprep"
        export_filename = prompt_with_context(
            self.processor,
            "Export filename",
            default=default_filename,
            module="Template & Workflow Manager",
            description="Export file name",
        )
        export_path = Path.cwd() / export_filename
        
        # Export
        if self.user_data_manager.export_content(selected_id, export_path):
            self.console.print(f"[green]✓ {content_type.title()} exported to {export_path}[/green]")
        else:
            self.console.print(f"[red]Failed to export {content_type}[/red]")
            
    def _import_content_wizard(self):
        """Wizard for importing templates or workflows."""
        import_path = prompt_with_context(
            self.processor,
            "Import file path",
            module="Template & Workflow Manager",
            description="Import file path",
        )
        import_file = Path(import_path)
        
        if not import_file.exists():
            self.console.print(f"[red]File not found: {import_path}[/red]")
            return
            
        author = prompt_with_context(
            self.processor,
            "Author name for imported content",
            default="User",
            module="Template & Workflow Manager",
            description="Author name for imported content",
        )
        
        content_id = self.user_data_manager.import_content(import_file, author=author)
        if content_id:
            self.console.print(f"[green]✓ Content imported successfully! ID: {content_id}[/green]")
        else:
            self.console.print("[red]Failed to import content[/red]")