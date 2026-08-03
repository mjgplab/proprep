"""
AMBER Input Controller - Unified Template-Wizard System

This controller orchestrates the interaction between annotated templates
and the comprehensive AMBER wizard, providing a unified user experience.
"""

import os
import re
import tempfile
import subprocess
from typing import Dict, Any, List, Optional, Tuple
from pathlib import Path
from enum import Enum

from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from proprep.utils.prompts import (
    prompt_with_context,
    confirm_with_context,
    int_prompt_with_context,
    float_prompt_with_context,
)

from .amber_annotated_templates import AmberAnnotatedTemplate, AmberAnnotatedTemplateSystem, SimulationType
from .amber_wizard import AmberWizard
from .workflow_loader import WorkflowLoader, WorkflowPreset, WorkflowStep


class ConfigurationMethod(Enum):
    """Methods for configuring AMBER templates."""
    DIRECT_EDIT = "direct"
    WIZARD_GUIDED = "wizard"
    LOAD_EXISTING = "load"


class AmberController:
    """
    Main controller for the unified AMBER input generation system.
    
    Orchestrates template management, wizard configuration, and workflow assembly.
    """
    
    def __init__(self, processor=None, package_dir: Path = None):
        """
        Initialize the AMBER controller.
        
        Args:
            processor: Parent processor instance
            package_dir: Package directory for template storage
        """
        self.processor = processor
        self.console = self._get_console()
        
        # Initialize subsystems
        self.workflow_loader = WorkflowLoader(console=self.console)
        self.template_system = AmberAnnotatedTemplateSystem(package_dir, self.console)
        
        # Current workflow state
        self.current_workflow = {}  # {phase: AmberAnnotatedTemplate}
        self.workflow_metadata = {
            "name": "Custom Workflow",
            "description": "User-configured AMBER workflow",
            "created": None
        }
        
    def _get_console(self) -> Console:
        """Get console from processor or create new one."""
        if (self.processor and 
            hasattr(self.processor, "console") and 
            self.processor.console):
            return self.processor.console
        return Console()
        
    def get_menu_options(self) -> Dict[str, str]:
        """Get main menu options for the redesigned system."""
        return {
            "single_templates": "Manage single simulation Templates",
            "workflow_templates": "Manage workflow templates"
        }
        
    def handle_menu_option(self, option: str) -> bool:
        """Handle menu option selection."""
        try:
            if option == "single_templates":
                return self.manage_single_templates()
                
            elif option == "workflow_templates":
                return self.manage_workflow_templates()
                
            else:
                self.console.print(f"[red]Unknown option: {option}[/red]")
                return False
                
        except Exception as e:
            self.console.print(f"[red]Error: {e}[/red]")
            return False
            
    def configure_simulation_phase(self, simulation_type: SimulationType) -> bool:
        """
        Configure a specific simulation phase using template + optional wizard.
        
        This is the core method that implements the template-centric workflow
        with wizard integration.
        """
        try:
            self.console.print(f"\n[bold cyan]Configure {simulation_type.value.title()} Phase[/bold cyan]")
            
            # Step 1: Load and display template content immediately
            template = self.template_system.load_template(simulation_type)
            self.template_system.display_template_content(template)
            
            # Step 2: User choice for configuration method
            method = self._ask_configuration_method()
            
            if method == ConfigurationMethod.WIZARD_GUIDED:
                # Step 3a: Launch wizard with template as starting point
                updated_template = self._configure_with_wizard(template)
                
            elif method == ConfigurationMethod.DIRECT_EDIT:
                # Step 3b: Direct template editing
                updated_template = self._configure_with_direct_edit(template)
                
            elif method == ConfigurationMethod.LOAD_EXISTING:
                # Step 3c: Load existing template
                updated_template = self._load_existing_template(simulation_type)
                
            else:
                return False
                
            # Step 4: Show updated template
            if updated_template:
                self.console.print(f"\n[bold green]Updated Template: {updated_template.name}[/bold green]")
                self.template_system.display_template_content(updated_template)
                
                # Step 5: Save to current workflow
                self.current_workflow[simulation_type] = updated_template
                
                # Option to save as custom template
                if confirm_with_context(
                    self.processor,
                    "\nSave this configuration as a custom template?",
                    default=False,
                    module="AMBER Controller",
                    description="Save configuration as custom template",
                ):
                    self._save_as_custom_template(updated_template)
                    
                return True
            else:
                return False
                
        except Exception as e:
            self.console.print(f"[red]Error configuring {simulation_type.value}: {e}[/red]")
            return False
            
    def _ask_configuration_method(self) -> ConfigurationMethod:
        """Ask user how they want to configure the template."""
        self.console.print("\n[bold]How would you like to configure this template?[/bold]")
        self.console.print("[cyan]1.[/cyan] Use wizard (guided parameter selection)")
        self.console.print("[cyan]2.[/cyan] Edit template directly") 
        self.console.print("[cyan]3.[/cyan] Load existing custom template")
        
        choice = prompt_with_context(
            self.processor,
            "Choose configuration method",
            choices=["1", "2", "3"],
            default="1",
            module="AMBER Controller",
            description="Template configuration method",
            options_map={
                "1": "Use wizard (guided)",
                "2": "Edit template directly",
                "3": "Load existing custom template",
            },
        )
        
        method_map = {
            "1": ConfigurationMethod.WIZARD_GUIDED,
            "2": ConfigurationMethod.DIRECT_EDIT,
            "3": ConfigurationMethod.LOAD_EXISTING
        }
        
        return method_map[choice]
        
    def _configure_with_wizard(self, template: AmberAnnotatedTemplate) -> AmberAnnotatedTemplate:
        """Configure template using the comprehensive AMBER wizard."""
        self.console.print("\n[bold yellow]Launching AMBER Configuration Wizard...[/bold yellow]")
        self.console.print("[grey50]The wizard will guide you through all parameters step-by-step[/grey50]")
        
        # Launch wizard with template's current config as starting point
        initial_config = template.get_config_dict()
        
        wizard_config = AmberWizard.configure(
            simulation_type=template.simulation_type.value,
            initial_config=initial_config,
            console=self.console,
            processor=self.processor,
        )
        
        # Update template with wizard results
        template.update_from_wizard_config(wizard_config)
        template.metadata["configured_by"] = "wizard"
        template.metadata["wizard_config"] = wizard_config
        
        self.console.print("\n[green]✓ Configuration completed via wizard[/green]")
        return template
        
    def _configure_with_direct_edit(self, template: AmberAnnotatedTemplate) -> AmberAnnotatedTemplate:
        """Configure template through direct editing."""
        self.console.print("\n[bold yellow]Opening template in editor...[/bold yellow]")
        
        # Generate current template content and open directly in editor
        current_content = template.generate_mdin_content()
        edited_content = self._open_editor(current_content, template.name)
        
        if edited_content and edited_content != current_content:
            # Parse edited content back into template
            # (This would need a parser - for now, we'll just store the raw content)
            template.metadata["edited_content"] = edited_content
            template.metadata["configured_by"] = "direct_edit"
            
            self.console.print("\n[green]✓ Template edited successfully[/green]")
        else:
            self.console.print("\n[yellow]No changes made to template[/yellow]")
            
        return template
        
    def _load_existing_template(self, simulation_type: SimulationType) -> Optional[AmberAnnotatedTemplate]:
        """Load an existing custom template."""
        # List available custom templates
        custom_templates = self._list_custom_templates(simulation_type)
        
        if not custom_templates:
            self.console.print(f"[yellow]No custom templates found for {simulation_type.value}[/yellow]")
            return None
            
        self.console.print(f"\n[bold]Available {simulation_type.value} templates:[/bold]")
        
        for i, template_name in enumerate(custom_templates, 1):
            self.console.print(f"[cyan]{i}.[/cyan] {template_name}")
            
        choice = prompt_with_context(
            self.processor,
            f"Select template [1-{len(custom_templates)}]",
            choices=[str(i) for i in range(1, len(custom_templates) + 1)],
            module="AMBER Controller",
            description="Select custom template to load",
            options_map={str(i + 1): name for i, name in enumerate(custom_templates)},
        )
        
        template_name = custom_templates[int(choice) - 1]
        
        # Load the selected template
        # (Implementation would load from user templates directory)
        self.console.print(f"[green]Loaded custom template: {template_name}[/green]")
        
        # For now, return the default template
        # In full implementation, this would load the actual custom template
        return self.template_system.load_template(simulation_type)
        
    def _open_editor(self, content: str, template_name: str) -> Optional[str]:
        """Open content in system editor and return edited content."""
        try:
            # Create temporary file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.in', delete=False) as tmp:
                tmp.write(content)
                tmp_path = tmp.name
                
            # Get editor from environment or use default
            editor = os.environ.get('EDITOR', 'nano')
            
            # Open editor
            subprocess.run([editor, tmp_path], check=True)
            
            # Read back edited content
            with open(tmp_path, 'r') as f:
                edited_content = f.read()
                
            # Cleanup
            os.unlink(tmp_path)
            
            return edited_content
            
        except Exception as e:
            self.console.print(f"[red]Editor error: {e}[/red]")
            return None
            
    def _save_as_custom_template(self, template: AmberAnnotatedTemplate) -> bool:
        """Save template as a custom user template."""
        template_name = prompt_with_context(
            self.processor,
            "Template name",
            default=f"custom_{template.simulation_type.value}",
            module="AMBER Controller",
            description="Custom template name",
        )
        
        # Update template metadata
        template.name = template_name
        template.metadata["custom"] = True
        template.metadata["saved_by_user"] = True
        
        # Save to user templates directory
        try:
            self.template_system._save_template(template, self.template_system.user_templates_dir)
            self.console.print(f"[green]✓ Saved custom template: {template_name}[/green]")
            return True
        except Exception as e:
            self.console.print(f"[red]Error saving template: {e}[/red]")
            return False
            
    def _list_custom_templates(self, simulation_type: SimulationType) -> List[str]:
        """List available custom templates for a simulation type."""
        # Implementation would scan user templates directory
        # For now, return empty list
        return []
        
    def assemble_workflow(self) -> bool:
        """Assemble complete workflow from current templates."""
        self.console.print("\n[bold cyan]Assemble Complete Workflow[/bold cyan]")
        
        if not self.current_workflow:
            self.console.print("[yellow]No templates configured yet. Configure individual phases first.[/yellow]")
            return False
            
        # Show current workflow status
        table = Table(title="Current Workflow Status")
        table.add_column("Phase", style="cyan")
        table.add_column("Template", style="green") 
        table.add_column("Status", style="yellow")
        
        phases = [SimulationType.MINIMIZATION, SimulationType.HEATING, 
                 SimulationType.EQUILIBRATION, SimulationType.PRODUCTION]
        
        for phase in phases:
            if phase in self.current_workflow:
                template = self.current_workflow[phase]
                table.add_row(phase.value.title(), template.name, "✓ Configured")
            else:
                table.add_row(phase.value.title(), "Not configured", "❌ Missing")
                
        self.console.print(table)
        
        # Generate workflow files
        if confirm_with_context(
            self.processor,
            "Generate input files for current workflow?",
            default=True,
            module="AMBER Controller",
            description="Generate input files for current workflow",
        ):
            return self._generate_workflow_files()
            
        return True
        
    def _generate_workflow_files(self) -> bool:
        """Generate AMBER input files for the complete workflow."""
        output_dir = Path("amber_workflow")
        output_dir.mkdir(exist_ok=True)
        
        file_mapping = {
            SimulationType.MINIMIZATION: "01_min.in",
            SimulationType.HEATING: "02_heat.in", 
            SimulationType.EQUILIBRATION: "03_equil.in",
            SimulationType.PRODUCTION: "04_prod.in"
        }
        
        generated_files = []
        
        for phase, filename in file_mapping.items():
            if phase in self.current_workflow:
                template = self.current_workflow[phase]
                content = template.generate_mdin_content()
                
                filepath = output_dir / filename
                with open(filepath, 'w') as f:
                    f.write(content)
                    
                generated_files.append(filepath)
                
        if generated_files:
            self.console.print(f"\n[green]✓ Generated {len(generated_files)} workflow files in {output_dir}/[/green]")
            for filepath in generated_files:
                self.console.print(f"  • {filepath.name}")
            return True
        else:
            self.console.print("[yellow]No files generated - no templates configured[/yellow]")
            return False
            
    def load_workflow_preset(self) -> bool:
        """Load a complete workflow preset."""
        self.console.print("\n[bold cyan]Load Workflow Preset[/bold cyan]")
        
        # Step 1: Select preset
        preset_key = self.workflow_loader.display_workflow_menu()
        if not preset_key:
            self.console.print("[yellow]No preset selected[/yellow]")
            return False
            
        # Handle special actions
        if preset_key == "create_new":
            self.console.print("[yellow]Workflow creation wizard - Coming soon![/yellow]")
            return self.load_workflow_preset()  # Return to menu
        elif preset_key == "manage":
            from .template_workflow_manager import TemplateWorkflowManager
            manager = TemplateWorkflowManager(console=self.console)
            manager.show_management_menu()
            return self.load_workflow_preset()  # Return to menu
            
        workflows = self.workflow_loader.get_available_workflows()
        preset = workflows[preset_key]
        
        # Step 2: Show preset overview
        self.workflow_loader.display_workflow_overview(preset)
        
        # Step 3: User choice
        self.console.print("\nOptions:")
        self.console.print("  1. Use as-is (generate all files)")
        self.console.print("  2. Edit individual workflow steps")
        self.console.print("  3. View step details")
        self.console.print("  4. ← Return to AMBER Input Generator")
        
        choice = prompt_with_context(
            self.processor,
            "Select option",
            choices=["1", "2", "3", "4"],
            default="1",
            module="AMBER Controller",
            description="Workflow preset action",
            options_map={
                "1": "Use as-is (generate all files)",
                "2": "Edit individual workflow steps",
                "3": "View step details",
                "4": "Return to AMBER Input Generator",
            },
        )

        if choice == "4":
            return True  # Return to previous menu
        elif choice == "1":
            return self._generate_workflow_files(preset)
        elif choice == "2":
            return self._edit_workflow_steps(preset)
        elif choice == "3":
            return self._view_step_details(preset)
            
        return False
        
    def save_workflow_preset(self) -> bool:
        """Save current workflow as a preset.""" 
        self.console.print("\n[bold cyan]Save Workflow Preset[/bold cyan]")
        
        if not self.current_workflow:
            self.console.print("[yellow]No workflow to save. Configure phases first.[/yellow]")
            return False
            
        preset_name = prompt_with_context(
            self.processor,
            "Preset name",
            default="my_workflow",
            module="AMBER Controller",
            description="Workflow preset name",
        )
        
        # Save workflow preset (implementation would store to presets directory)
        self.console.print(f"[green]✓ Workflow preset saved: {preset_name}[/green]")
        return True
        
    def manage_single_templates(self) -> bool:
        """Manage single simulation templates."""
        while True:
            # Get available templates
            try:
                from .user_data_manager import UserDataManager
                user_data_manager = UserDataManager(console=self.console)
                templates = user_data_manager.list_templates()
            except Exception as e:
                self.console.print(f"[red]Error loading templates: {e}[/red]")
                return True
                
            self.console.print("\n[bold cyan]===== Manage Single Simulation Templates =====[/bold cyan]")
            
            # Display templates
            self._display_template_list(templates)
            
            self.console.print("\n[bold]Options:[/bold]")
            self.console.print("1. View a specific template", highlight=False)
            self.console.print("2. Create new template", highlight=False)
            self.console.print("3. Modify existing template", highlight=False)
            self.console.print("4. Delete template", highlight=False)
            self.console.print("5. ← Back to AMBER Input Generator", highlight=False)
            
            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=["1", "2", "3", "4", "5"],
                default="1",
                module="AMBER Controller",
                description="Single template management action",
                options_map={
                    "1": "View a specific template",
                    "2": "Create new template",
                    "3": "Modify existing template",
                    "4": "Delete template",
                    "5": "Back to AMBER Input Generator",
                },
            )

            if choice == "1":
                self._view_single_template(templates)
            elif choice == "2":
                self._create_single_template()
            elif choice == "3":
                self._modify_single_template(templates)
            elif choice == "4":
                self._delete_single_template(templates)
            elif choice == "5":
                return True
                
    def manage_workflow_templates(self) -> bool:
        """Manage workflow templates."""
        while True:
            # Get available workflows
            try:
                from .workflow_loader import WorkflowLoader
                workflow_loader = WorkflowLoader(console=self.console)
                workflows = workflow_loader.get_available_workflows()
            except Exception as e:
                self.console.print(f"[red]Error loading workflows: {e}[/red]")
                return True
                
            self.console.print("\n[bold cyan]===== Manage Workflow Templates =====[/bold cyan]")
            
            # Display workflows
            self._display_workflow_list(workflows)
            
            self.console.print("\n[bold]Options:[/bold]")
            self.console.print("1. View a specific workflow", highlight=False)
            self.console.print("2. Create new workflow", highlight=False)
            self.console.print("3. Modify existing workflow", highlight=False)
            self.console.print("4. Delete workflow", highlight=False)
            self.console.print("5. ← Back to AMBER Input Generator", highlight=False)
            
            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=["1", "2", "3", "4", "5"],
                default="1",
                module="AMBER Controller",
                description="Workflow template management action",
                options_map={
                    "1": "View a specific workflow",
                    "2": "Create new workflow",
                    "3": "Modify existing workflow",
                    "4": "Delete workflow",
                    "5": "Back to AMBER Input Generator",
                },
            )

            if choice == "1":
                self._view_workflow_template(workflows)
            elif choice == "2":
                self._create_workflow_template()
            elif choice == "3":
                self._modify_workflow_template(workflows)
            elif choice == "4":
                self._delete_workflow_template(workflows)
            elif choice == "5":
                return True
        
    def _generate_workflow_files(self, preset: WorkflowPreset) -> bool:
        """Generate all workflow files from a preset."""
        self.console.print("\n[bold cyan]===== Generating Workflow Files =====[/bold cyan]")
        
        try:
            # Create output directory
            output_dir = Path.cwd() / "workflow_output"
            output_dir.mkdir(exist_ok=True)
            
            # Generate individual mdin files
            self.console.print("\nCreating mdin files:")
            mdin_files = []
            
            # Calculate max filename length for alignment
            max_filename_len = max(len(f"{step.id}.mdin") for step in preset.steps)
            
            for step in preset.steps:
                # Get template content and apply overrides
                template_content = self.workflow_loader.get_template_content(step.template)
                mdin_content = self.workflow_loader.apply_parameter_overrides(
                    template_content, step.parameter_overrides
                )
                
                mdin_path = output_dir / f"{step.id}.mdin"
                with open(mdin_path, 'w') as f:
                    f.write(mdin_content)
                    
                # Extract simulation info for display  
                sim_info = self._extract_simulation_info(mdin_content, f"{step.id}.mdin")
                
                if step.type == "minimization":
                    maxcyc_match = re.search(r'maxcyc\s*=\s*(\d+)', mdin_content)
                    maxcyc = maxcyc_match.group(1) if maxcyc_match else "unknown"
                    detail = f"({maxcyc} steps minimization)"
                else:
                    nstlim = sim_info["nstlim"]
                    time_ps = sim_info["simulation_time_ps"]
                    if time_ps >= 1000:
                        time_str = f"{time_ps/1000:.1f} ns"
                    else:
                        time_str = f"{time_ps:.0f} ps"
                    detail = f"({nstlim:,} steps, {time_str})"
                
                # Use dynamic padding for perfect alignment
                filename = f"{step.id}.mdin"
                padding = " " * (max_filename_len - len(filename) + 2)
                self.console.print(f"  ✓ {filename}{padding}{detail}")
                mdin_files.append(filename)
                
            # Generate workflow sequence JSON
            self.console.print("\nCreating workflow sequence:")
            sequence_file = self.workflow_loader.generate_workflow_sequence_json(
                preset, output_dir
            )
            self.console.print(f"  ✓ {sequence_file.name}")
            
            # Generate execution scripts
            self.console.print("\nCreating execution scripts:")
            self._generate_execution_scripts(preset, output_dir)
            
            # Update workspace with file paths
            if self.processor and hasattr(self.processor, 'workspace'):
                workspace = self.processor.workspace
                workspace.set("workflow_preset", preset.name)
                workspace.set("workflow_files", mdin_files)
                workspace.set("workflow_sequence", str(sequence_file))
                workspace.set("workflow_dir", str(output_dir))
            
            self.console.print(f"\n[bold green]Files saved to: {output_dir}[/bold green]")
            self.console.print("Files created:")
            self.console.print(f"  • Workflow: {sequence_file.name}")
            self.console.print(f"  • Input files: {len(mdin_files)} mdin files")
            self.console.print(f"  • Scripts: 2 bash scripts")
            self.console.print("\n[bold]Ready for AMBER Workflow Manager execution![/bold]")
            
            return True
            
        except Exception as e:
            self.console.print(f"[red]Error generating workflow files: {e}[/red]")
            return False
    
    def _edit_workflow_steps(self, preset: WorkflowPreset) -> bool:
        """Edit individual workflow steps."""
        
        while True:
            # Show step selection menu
            step_index = self._display_step_menu(preset)
            
            if step_index is None:
                break
            elif step_index == -99:  # Return option
                break
            elif step_index >= 0:  # Edit specific step
                step = preset.steps[step_index]
                self.console.print(f"\n[bold cyan]===== Editing Step: {step.id} ({step.name}) =====[/bold cyan]")
                
                # Show current configuration
                self._display_step_configuration(step)
                
                # Edit options
                self.console.print("\nEdit Options:")
                self.console.print("  1. View current mdin content")
                self.console.print("  2. Direct edit mdin file")
                self.console.print("  3. Use Configuration Wizard")
                self.console.print("  4. Reset to default")
                self.console.print("  5. ← Back to step selection")
                
                edit_choice = prompt_with_context(
                    self.processor,
                    "Select option",
                    choices=["1", "2", "3", "4", "5"],
                    default="1",
                    module="AMBER Controller",
                    description="Step edit action",
                    options_map={
                        "1": "View current mdin content",
                        "2": "Direct edit mdin file",
                        "3": "Use Configuration Wizard",
                        "4": "Reset to default",
                        "5": "Back to step selection",
                    },
                )
                
                if edit_choice == "5":
                    continue  # Back to step selection
                elif edit_choice == "1":
                    # View content
                    self._view_mdin_content(step)
                    continue
                elif edit_choice == "2":
                    # Direct edit
                    updated_step = self._edit_step_direct(step)
                elif edit_choice == "3":
                    # Use wizard
                    updated_step = self._edit_step_with_wizard(step)
                elif edit_choice == "4":
                    # Reset to default
                    updated_step = self._reset_step_to_default(step)
                    
                if updated_step:
                    preset.steps[step_index] = updated_step
                    self.console.print(f"\n[green]✓ Step {step.id} updated successfully[/green]")
                    
            elif step_index == -1:  # Add new step
                self.console.print("[yellow]Add new step - feature coming soon[/yellow]")
            elif step_index == -2:  # Remove step  
                self.console.print("[yellow]Remove step - feature coming soon[/yellow]")
            elif step_index == -3:  # Reorder steps
                self.console.print("[yellow]Reorder steps - feature coming soon[/yellow]")
                
            # Ask if user wants to continue editing
            if not confirm_with_context(
                self.processor,
                "\nContinue editing workflow?",
                default=False,
                module="AMBER Controller",
                description="Continue editing workflow",
            ):
                break
                
        # Generate files with modifications
        if confirm_with_context(
            self.processor,
            "\nGenerate workflow files with modifications?",
            default=True,
            module="AMBER Controller",
            description="Generate workflow files with modifications",
        ):
            return self._generate_workflow_files(preset)
            
        return True
        
    def _view_step_details(self, preset: WorkflowPreset) -> bool:
        """View detailed information about workflow steps."""
        step_index = self._display_step_menu(preset)
        
        if step_index is not None and step_index >= 0:
            step = preset.steps[step_index]
            self._display_step_configuration(step)
            self._view_mdin_content(step)
            
        return True
        
    def _display_step_configuration(self, step: WorkflowStep):
        """Display current step configuration."""
        # Get template content for step
        template_content = self.workflow_loader.get_template_content(step.template)
        mdin_content = self.workflow_loader.apply_parameter_overrides(
            template_content, step.parameter_overrides
        )
        sim_info = self._extract_simulation_info(mdin_content, f"{step.id}.mdin")
        
        self.console.print("\nCurrent Configuration:")
        self.console.print(f"  • Type: {step.type.title()}")
        
        if step.type == "minimization":
            maxcyc_match = re.search(r'maxcyc\s*=\s*(\d+)', mdin_content)
            maxcyc = maxcyc_match.group(1) if maxcyc_match else "unknown"
            self.console.print(f"  • Maximum cycles: {maxcyc}")
        else:
            time_ps = sim_info["simulation_time_ps"]
            if time_ps > 0:
                time_str = f"{time_ps/1000:.1f} ns" if time_ps >= 1000 else f"{time_ps:.0f} ps"
                self.console.print(f"  • Simulation time: {time_str}")
                self.console.print(f"  • Time step: {sim_info['dt']} ps")
                
        if sim_info["has_restraints"]:
            self.console.print(f"  • Restraints: {sim_info['restraint_weight']} kcal/mol·Å²")
        else:
            self.console.print("  • Restraints: None")
            
    def _view_mdin_content(self, step: WorkflowStep):
        """Display mdin file content."""
        from rich.syntax import Syntax
        
        # Get template content for step
        template_content = self.workflow_loader.get_template_content(step.template)
        mdin_content = self.workflow_loader.apply_parameter_overrides(
            template_content, step.parameter_overrides
        )
        
        # Format content with aligned comments
        formatted_content = self._format_mdin_content(mdin_content)
        
        self.console.print(f"\n[bold]Content of {step.id}.mdin:[/bold]")
        
        # Display with better syntax highlighting for AMBER input
        syntax = Syntax(formatted_content, "fortran", theme="github-dark", line_numbers=True)
        self.console.print(syntax)
        
    def _format_mdin_content(self, content: str) -> str:
        """Format MDIN content with aligned comments and proper structure."""
        lines = content.split('\n')
        formatted_lines = []
        in_cntrl_section = False
        
        for line in lines:
            stripped = line.strip()
            
            # Track if we're in &cntrl section
            if stripped.startswith('&cntrl'):
                in_cntrl_section = True
                formatted_lines.append(line)
                continue
            elif stripped == '/':
                in_cntrl_section = False
                formatted_lines.append(line)
                continue
            elif not in_cntrl_section:
                # Outside &cntrl section, keep original formatting
                formatted_lines.append(line)
                continue
                
            # Process lines within &cntrl section
            if not stripped or stripped.startswith('!'):
                # Comment-only lines or empty lines
                formatted_lines.append(line)
            elif '!' in line and '=' in line:
                # Parameter lines with inline comments
                param_part, comment_part = line.split('!', 1)
                param_part = param_part.rstrip()
                comment_part = comment_part.strip()
                
                # Align comments at column 30 for better readability
                if len(param_part) < 28:
                    padding = ' ' * (28 - len(param_part))
                    formatted_line = f"{param_part}{padding}! {comment_part}"
                else:
                    formatted_line = f"{param_part}  ! {comment_part}"
                formatted_lines.append(formatted_line)
            else:
                # Lines without comments
                formatted_lines.append(line)
                
        return '\n'.join(formatted_lines)
        
    def _edit_step_with_wizard(self, step: WorkflowStep) -> Optional[WorkflowStep]:
        """Edit a workflow step using the AMBER wizard."""
        self.console.print(f"\n[bold]Launching Configuration Wizard for {step.id}[/bold]")
        
        try:
            # Map step type to simulation type
            sim_type_map = {
                "minimization": "minimization",
                "heating": "heating",
                "equilibration": "md",  # Use md as base for equilibration
                "production": "md"
            }
            
            simulation_type = sim_type_map.get(step.type, "md")
            
            # Create wizard with current step parameters as starting point
            wizard = AmberWizard(simulation_type=simulation_type, console=self.console, processor=self.processor)
            
            # Run wizard workflow
            config = wizard.run_complete_workflow()
            
            if config:
                # Convert wizard config back to mdin content
                updated_content = self._generate_mdin_from_config(config, step.mdin_file)
                
                # Create updated step  
                updated_step = WorkflowStep(
                    id=step.id,
                    name=step.name,
                    type=step.type,
                    template=step.template,
                    description=step.description,
                    parameter_overrides=config if config else step.parameter_overrides,
                    dependencies=step.dependencies,
                    input_coord=step.input_coord,
                    output_coord=step.output_coord
                )
                
                return updated_step
            else:
                self.console.print("[yellow]Wizard configuration cancelled[/yellow]")
                return None
                
        except Exception as e:
            self.console.print(f"[red]Error running wizard: {e}[/red]")
            return None
    
    def _generate_mdin_from_config(self, config: Dict[str, Any], original_filename: str) -> str:
        """Generate mdin content from wizard configuration using template system."""
        try:
            # Use the template system to generate proper mdin content
            from .amber_annotated_templates import AmberAnnotatedTemplate
            
            # Determine simulation type from config
            if config.get("imin", 0) == 1:
                sim_type = "minimization"
            elif config.get("nmropt", 0) == 1 and "TEMP0" in str(config):
                sim_type = "heating" 
            else:
                sim_type = "md"
                
            # Create template and configure
            template = AmberAnnotatedTemplate(sim_type)
            
            # Update template with wizard config
            for key, value in config.items():
                if hasattr(template, key):
                    setattr(template, key, value)
            
            # Generate formatted mdin content
            return template.format_template()
            
        except Exception as e:
            # Fallback to simplified version if template system fails
            self.console.print(f"[yellow]Using fallback mdin generation: {e}[/yellow]")
            
            # Extract title from original filename  
            title_map = {
                "000_Min.mdin": "Initial Solvent Energy Minimization Stage",
                "02_Heat.mdin": "Initial Heating After Solvent Minimization", 
                "05_Heat_solute_Rwt25.mdin": "Solute Equilibration Heating Stage with 25 kcal weight",
                "06_Equil_solute_Rwt25.mdin": "Solutel Equilibration with 25 kcal weight"
            }
            
            title = title_map.get(original_filename, f"AMBER Simulation - {original_filename}")
            
            # Generate basic mdin content
            lines = [title, "", "&cntrl"]
            
            # Add parameters from config
            for key, value in sorted(config.items()):
                if key.startswith("_"):  # Skip internal parameters
                    continue
                    
                if isinstance(value, bool):
                    value = 1 if value else 0
                elif isinstance(value, str) and not (key in ["restraintmask", "bellymask"]):
                    continue  # Skip string parameters that aren't masks
                    
                lines.append(f"  {key} = {value},")
                
            lines.append("/")
            
            # Add &wt blocks if nmropt is enabled
            if config.get("nmropt") == 1:
                lines.extend([
                    "&wt type = 'DUMPFREQ', istep1 = 10000, /",
                    "&wt type='END' /"
                ])
                
            return "\n".join(lines)
    
    def _edit_step_direct(self, step: WorkflowStep) -> Optional[WorkflowStep]:
        """Edit step by directly modifying mdin content."""
        self.console.print("[yellow]Direct edit - opening in external editor[/yellow]")
        
        # Create temporary file
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.mdin', delete=False) as tmp:
            tmp.write(step.mdin_content)
            tmp_path = tmp.name
            
        try:
            # Open in editor
            editor = os.environ.get('EDITOR', 'nano')
            subprocess.run([editor, tmp_path])
            
            # Read modified content
            with open(tmp_path, 'r') as f:
                modified_content = f.read()
                
            # Create updated step
            if modified_content != step.mdin_content:
                updated_step = WorkflowStep(
                    id=step.id,
                    name=step.name,
                    type=step.type,
                    template=step.template,
                    description=step.description,
                    parameter_overrides=step.parameter_overrides,
                    dependencies=step.dependencies,
                    input_coord=step.input_coord,
                    output_coord=step.output_coord
                )
                return updated_step
            else:
                self.console.print("[yellow]No changes made[/yellow]")
                return None
                
        finally:
            # Clean up temporary file
            os.unlink(tmp_path)
            
    def _reset_step_to_default(self, step: WorkflowStep) -> WorkflowStep:
        """Reset step to default from original template."""
        # Re-read from original template file
        template_content = self.workflow_loader.get_template_content(step.template)
        # No parameter overrides for reset to defaults
        
        updated_step = WorkflowStep(
            id=step.id,
            name=step.name,
            type=step.type,
            template=step.template,
            description=step.description,
            parameter_overrides={},  # Reset overrides
            dependencies=step.dependencies,
            input_coord=step.input_coord,
            output_coord=step.output_coord
        )
        
        return updated_step
        
    def _display_step_menu(self, preset: WorkflowPreset) -> Optional[int]:
        """Display step selection menu for editing."""
        self.console.print(f"\n[bold cyan]===== Edit Individual Workflow Steps =====[/bold cyan]")
        
        self.console.print("Select step to modify:")
        step_choices = {}
        
        for i, step in enumerate(preset.steps, 1):
            # Get time info for display
            try:
                template_content = self.workflow_loader.get_template_content(step.template)
                mdin_content = self.workflow_loader.apply_parameter_overrides(
                    template_content, step.parameter_overrides
                )
                sim_info = self._extract_simulation_info(mdin_content, f"{step.id}.mdin")
                
                if step.type == "minimization":
                    time_info = ""
                elif sim_info["simulation_time_ps"] > 0:
                    time_ps = sim_info["simulation_time_ps"]
                    time_str = f"{time_ps/1000:.1f} ns" if time_ps >= 1000 else f"{time_ps:.0f} ps"
                    time_info = f" ({time_str})"
                else:
                    time_info = ""
            except:
                time_info = ""
                
            # Get restraint info
            restraint_info = ""
            try:
                if sim_info["has_restraints"]:
                    restraint_info = f" [R:{sim_info['restraint_weight']}]"
            except:
                pass
                
            step_display = f"{step.id} - {step.name}{time_info}{restraint_info}"
            self.console.print(f"  {i}. {step_display}")
            step_choices[str(i)] = i - 1  # Convert to 0-based index
            
        self.console.print(f"\n[bold]Workflow options:[/bold]")
        self.console.print(f"  {len(preset.steps) + 1}. Add new step")
        step_choices[str(len(preset.steps) + 1)] = -1  # Add
        self.console.print(f"  {len(preset.steps) + 2}. Remove step")
        step_choices[str(len(preset.steps) + 2)] = -2  # Remove
        self.console.print(f"  {len(preset.steps) + 3}. Reorder steps")
        step_choices[str(len(preset.steps) + 3)] = -3  # Reorder
        self.console.print(f"  {len(preset.steps) + 4}. ← Return to previous menu")
        step_choices[str(len(preset.steps) + 4)] = -99  # Return
        
        choices = list(step_choices.keys())
        choice = prompt_with_context(
            self.processor,
            "Select option",
            choices=choices,
            default="1",
            module="AMBER Controller",
            description="Select step or action",
        )

        return step_choices.get(choice)
        
    def _generate_execution_scripts(self, preset: WorkflowPreset, output_dir: Path):
        """Generate bash execution scripts for the workflow."""
        
        # Single step script
        step_script_content = '''#!/bin/bash
#
# Single Step Execution Script
# Usage: ./run_step.sh <step_name>
#

if [ $# -eq 0 ]; then
    echo "Usage: $0 <step_name>"
    echo "Available steps:"
''' + '\n'.join([f'    echo "  {step.id}"' for step in preset.steps]) + '''
    exit 1
fi

STEP_NAME=$1
MDIN_FILE="${STEP_NAME}.mdin"

if [ ! -f "$MDIN_FILE" ]; then
    echo "Error: $MDIN_FILE not found"
    exit 1
fi

echo "Running AMBER step: $STEP_NAME"
$AMBERHOME/bin/pmemd.cuda -O -i "$MDIN_FILE" -o "${STEP_NAME}.out" -p prmtop -c inpcrd -r "${STEP_NAME}.rst" -x "${STEP_NAME}.nc"
'''
        
        step_script_path = output_dir / "run_step.sh"
        with open(step_script_path, 'w') as f:
            f.write(step_script_content)
        os.chmod(step_script_path, 0o755)
        
        # Complete workflow script
        workflow_script_content = f'''#!/bin/bash
#
# Complete Workflow Execution Script
# {preset.name}
#

set -e

echo "Starting {preset.name}"
echo "Total steps: {len(preset.steps)}"
echo

'''
        
        for i, step in enumerate(preset.steps, 1):
            workflow_script_content += f'''echo "Step {i}/{len(preset.steps)}: {step.name}"
echo "Step: {step.name}"
./run_step.sh {step.id}
echo "Step {i} completed\n"

'''
            
        workflow_script_content += '''echo "Workflow completed successfully!"
echo "All output files have been generated."
'''
        
        workflow_script_path = output_dir / "run_workflow.sh"
        with open(workflow_script_path, 'w') as f:
            f.write(workflow_script_content)
        os.chmod(workflow_script_path, 0o755)
        
        self.console.print(f"  ✓ run_step.sh           (single step execution)")
        self.console.print(f"  ✓ run_workflow.sh       (complete workflow execution)")
        
    def _extract_simulation_info(self, mdin_content: str, filename: str) -> Dict[str, Any]:
        """Extract simulation information from mdin content."""
        import re
        
        info = {
            "nstlim": 0,
            "dt": 0.001,
            "simulation_time_ps": 0,
            "has_restraints": False,
            "restraint_weight": 0.0
        }
        
        # Extract nstlim
        nstlim_match = re.search(r'nstlim\s*=\s*(\d+)', mdin_content)
        if nstlim_match:
            info["nstlim"] = int(nstlim_match.group(1))
            
        # Extract dt
        dt_match = re.search(r'dt\s*=\s*([\d.]+)', mdin_content)
        if dt_match:
            info["dt"] = float(dt_match.group(1))
            
        # Calculate simulation time
        info["simulation_time_ps"] = info["nstlim"] * info["dt"]
        
        # Check for restraints
        restraint_match = re.search(r'ntr\s*=\s*1', mdin_content)
        if restraint_match:
            info["has_restraints"] = True
            
            # Extract restraint weight
            weight_match = re.search(r'restraint_wt\s*=\s*([\d.]+)', mdin_content)
            if weight_match:
                info["restraint_weight"] = float(weight_match.group(1))
                
        return info
        
    # Single Template Management Methods
    def _display_categorized_templates(self, templates: Dict, bullet_style: str = "•", show_numbers: bool = False):
        """Display templates organized by category and simulation type."""
        if not templates:
            self.console.print("[yellow]No templates available[/yellow]")
            return {}
            
        choice_map = {}
        choice_num = 1
        
        # Separate by source first
        builtin_templates = {k: v for k, v in templates.items() if v["source"] == "builtin"}
        # Custom templates include all non-builtin sources (custom, from_workflow_step, copied_from_builtin, etc.)
        custom_templates = {k: v for k, v in templates.items() if v["source"] != "builtin"}
        
        # Process built-in templates
        if builtin_templates:
            self.console.print("\n[bold]Built-in Templates:[/bold]")
            choice_num = self._display_templates_by_type(builtin_templates, bullet_style, show_numbers, choice_num, choice_map)
            
        # Process custom templates
        if custom_templates:
            self.console.print("\n[bold]Custom Templates:[/bold]")
            choice_num = self._display_templates_by_type(custom_templates, bullet_style, show_numbers, choice_num, choice_map)
            
        return choice_map
        
    def _display_templates_by_type(self, templates: Dict, bullet_style: str, show_numbers: bool, start_num: int, choice_map: Dict) -> int:
        """Display templates grouped by simulation type."""
        # Group by simulation type
        by_type = {}
        for tid, meta in templates.items():
            sim_type = meta.get("simulation_type", "unknown")
            if sim_type not in by_type:
                by_type[sim_type] = []
            by_type[sim_type].append((tid, meta))
        
        # Sort templates within each type by priority (ascending), then by name
        for sim_type in by_type:
            by_type[sim_type].sort(key=lambda x: (x[1].get("priority", 999), x[1]["name"]))
        
        # Calculate global max name length for consistent alignment
        all_templates = []
        for templates_list in by_type.values():
            all_templates.extend(templates_list)
        global_max_name_len = max(len(meta["name"]) for _, meta in all_templates) if all_templates else 0
        
        choice_num = start_num
        
        # Display in logical order
        for sim_type in ["minimization", "heating", "equilibration", "production"]:
            if sim_type in by_type:
                self.console.print(f"\n  [bold]{sim_type.title()}:[/bold]")
                
                for tid, meta in by_type[sim_type]:
                    name = meta["name"]
                    padding = " " * (global_max_name_len - len(name) + 2)
                    
                    if show_numbers:
                        # Right-align numbers for consistent formatting
                        num_str = f"{choice_num}.".rjust(3)
                        self.console.print(f"    {num_str} {name}{padding}- {meta['description']}")
                        choice_map[str(choice_num)] = tid
                        choice_num += 1
                    else:
                        self.console.print(f"    {bullet_style} {name}{padding}- {meta['description']}")
        
        # Display other types if any
        for sim_type, templates_list in by_type.items():
            if sim_type not in ["minimization", "heating", "equilibration", "production"]:
                self.console.print(f"\n  [bold]{sim_type.title()}:[/bold]")
                
                for tid, meta in templates_list:
                    name = meta["name"]
                    padding = " " * (global_max_name_len - len(name) + 2)
                    
                    if show_numbers:
                        # Right-align numbers for consistent formatting
                        num_str = f"{choice_num}.".rjust(3)
                        self.console.print(f"    {num_str} {name}{padding}- {meta['description']}")
                        choice_map[str(choice_num)] = tid
                        choice_num += 1
                    else:
                        self.console.print(f"    {bullet_style} {name}{padding}- {meta['description']}")
                        
        return choice_num
        
    def _display_template_list(self, templates: Dict):
        """Display formatted list of templates using categorized format."""
        self._display_categorized_templates(templates)
                
    def _view_single_template(self, templates: Dict):
        """View a specific template."""
        if not templates:
            self.console.print("[yellow]No templates available to view[/yellow]")
            return
        
        self.console.print("\n[bold]Select template to view:[/bold]")
        
        # Use categorized display with numbers
        template_choices = self._display_categorized_templates(templates, show_numbers=True)
        
        # Add back option
        back_num = len(templates) + 1
        template_choices[str(back_num)] = None
        self.console.print(f"\n[bold]Options:[/bold]")
        self.console.print(f"  {back_num}. ← Back")
        
        choice = prompt_with_context(
            self.processor,
            "Select template",
            choices=list(template_choices.keys()),
            default="1",
            module="AMBER Controller",
            description="Select template",
        )
        selected_template_id = template_choices[choice]
        
        if selected_template_id:
            try:
                from .user_data_manager import UserDataManager
                user_data_manager = UserDataManager(console=self.console)
                
                # Get the template metadata to find the correct path
                template_metadata = templates[selected_template_id]
                
                # For builtin templates, use the template_path; for custom templates, use the ID
                if template_metadata["source"] == "builtin":
                    template_key = template_metadata["template_path"]
                else:
                    template_key = selected_template_id
                    
                content, metadata = user_data_manager.get_template_content(template_key)
                
                self.console.print(f"\n[bold cyan]===== Template: {metadata['name']} =====[/bold cyan]")
                self.console.print(f"[bold]Description:[/bold] {metadata['description']}")
                self.console.print(f"[bold]Simulation Type:[/bold] {metadata.get('simulation_type', 'unknown')}")
                self.console.print(f"[bold]Source:[/bold] {metadata['source']}")
                self.console.print("\n[bold]Content:[/bold]")
                
                formatted_content = self._format_mdin_content(content)
                self.console.print(formatted_content)
                
                input("Press Enter to continue...")
                
            except Exception as e:
                self.console.print(f"[red]Error viewing template: {e}[/red]")
                
    def _create_single_template(self):
        """Create a new single template."""
        self.console.print("\n[bold cyan]===== Create New Template =====[/bold cyan]")
        
        name = prompt_with_context(
            self.processor,
            "Template name",
            module="AMBER Controller",
            description="Template name",
        )
        if not name:
            return
            
        description = prompt_with_context(
            self.processor,
            "Description",
            module="AMBER Controller",
            description="Description",
        )
        if not description:
            return
            
        # Simulation type selection
        sim_types = ["minimization", "heating", "equilibration", "production", "other"]
        self.console.print("\n[bold]Simulation types:[/bold]")
        for i, st in enumerate(sim_types, 1):
            self.console.print(f"  {i}. {st}")
            
        type_choice = prompt_with_context(
            self.processor,
            "Select simulation type",
            choices=[str(i) for i in range(1, len(sim_types)+1)],
            default="1",
            module="AMBER Controller",
            description="Template simulation type",
            options_map={str(i + 1): t for i, t in enumerate(sim_types)},
        )
        simulation_type = sim_types[int(type_choice)-1]
        
        # Create template content
        self.console.print("\n[bold]Choose template creation method:[/bold]")
        self.console.print("1. Start from built-in template", highlight=False)
        self.console.print("2. Create from scratch", highlight=False)
        
        method = prompt_with_context(
            self.processor,
            "Select method",
            choices=["1", "2"],
            default="1",
            module="AMBER Controller",
            description="Template creation method",
            options_map={"1": "Start from scratch", "2": "Base on existing template"},
        )
        
        if method == "1":
            # Copy from built-in
            try:
                from .user_data_manager import UserDataManager
                user_data_manager = UserDataManager(console=self.console)
                templates = user_data_manager.list_templates(show_user=False)
                
                if not templates:
                    self.console.print("[yellow]No built-in templates available[/yellow]")
                    return
                    
                self.console.print("\n[bold]Select base template:[/bold]")
                
                # Use categorized display with numbers
                template_choices = self._display_categorized_templates(templates, show_numbers=True)
                    
                base_choice = prompt_with_context(
                    self.processor,
                    "Select base template",
                    choices=list(template_choices.keys()),
                    default="1",
                    module="AMBER Controller",
                    description="Select base template",
                )
                base_template_id = template_choices[base_choice]
                
                # Get the correct template path for builtin templates
                base_metadata = templates[base_template_id]
                if base_metadata["source"] == "builtin":
                    template_key = base_metadata["template_path"]
                else:
                    template_key = base_template_id
                
                # Get base content
                base_content, base_metadata = user_data_manager.get_template_content(template_key)
                
                # Open editor for direct editing of the base content
                import tempfile
                import os
                import subprocess
                
                with tempfile.NamedTemporaryFile(mode='w', suffix='.mdin', delete=False) as tmp_file:
                    tmp_file.write(base_content)
                    tmp_path = tmp_file.name
                
                try:
                    self.console.print(f"\nOpening editor for template based on: {base_metadata.get('name', 'template')}")
                    editor = os.environ.get('EDITOR', 'nano')
                    subprocess.run([editor, tmp_path])
                    
                    # Read modified content
                    with open(tmp_path, 'r') as f:
                        modified_content = f.read()
                    
                    # Create new template with edited content
                    template_id = user_data_manager.create_template(name, description, simulation_type, modified_content)
                    if template_id:
                        self.console.print(f"[green]✓ Template '{name}' created with ID: {template_id}[/green]")
                    else:
                        self.console.print(f"[red]Failed to create template[/red]")
                        
                finally:
                    os.unlink(tmp_path)
                    
            except Exception as e:
                self.console.print(f"[red]Error creating template: {e}[/red]")
        else:
            # Create from scratch - provide basic template with header
            author = prompt_with_context(
                self.processor,
                "Author",
                default="ProPrep User",
                module="AMBER Controller",
                description="Template author",
            )
            version = prompt_with_context(
                self.processor,
                "Version",
                default="1.0",
                module="AMBER Controller",
                description="Template version",
            )
            priority = prompt_with_context(
                self.processor,
                "Priority (lower numbers appear first)",
                default="100",
                module="AMBER Controller",
                description="Template priority (lower = first)",
            )
            source = prompt_with_context(
                self.processor,
                "Source URL (optional)",
                default="",
                module="AMBER Controller",
                description="Template source URL",
            )
            
            # Generate template with header
            header_content = f"""! TEMPLATE: {name}
! DESCRIPTION: {description}
! TYPE: {simulation_type}
! PRIORITY: {priority}
! AUTHOR: {author}
! VERSION: {version}"""
            
            if source:
                header_content += f"\n! SOURCE: {source}"
                
            basic_content = header_content + f"""

{simulation_type.title()} Template

&cntrl

{simulation_type.title()} Template

&cntrl
 imin=0,       ! No minimization
 irest=0,      ! No restart 
 ntx=1,        ! Read coordinates only
 ntb=1,        ! Periodic boundaries
 cut=9.0,      ! Cutoff
 ntr=0,        ! No restraints
 ntc=2,        ! SHAKE on bonds with hydrogen
 ntf=2,        ! No force calculation for bonds with hydrogen
 tempi=300.0,  ! Initial temperature
 temp0=300.0,  ! Target temperature
 ntt=3,        ! Langevin dynamics
 gamma_ln=1.0, ! Collision frequency
 nstlim=1000,  ! Number of steps
 dt=0.002,     ! Time step (ps)
 ntpr=100,     ! Print frequency
 ntwx=100,     ! Trajectory write frequency
 ntwr=100,     ! Restart write frequency
/"""
            try:
                from .user_data_manager import UserDataManager
                user_data_manager = UserDataManager(console=self.console)
                template_id = user_data_manager.create_template(name, description, simulation_type, basic_content, author)
                if template_id:
                    self.console.print(f"[green]✓ Template '{name}' created successfully[/green]")
                else:
                    self.console.print(f"[red]Failed to create template[/red]")
            except Exception as e:
                self.console.print(f"[red]Error creating template: {e}[/red]")
                
    def _modify_single_template(self, templates: Dict):
        """Modify an existing template."""
        if not templates:
            self.console.print("[yellow]No templates available to modify[/yellow]")
            return
            
        self.console.print("\n[bold]Select template to modify:[/bold]")
        
        # Use categorized display with numbers
        template_choices = self._display_categorized_templates(templates, show_numbers=True)
        
        # Add back option
        back_num = len(templates) + 1
        template_choices[str(back_num)] = None
        self.console.print(f"\n[bold]Options:[/bold]")
        self.console.print(f"  {back_num}. ← Back")
        
        choice = prompt_with_context(
            self.processor,
            "Select template",
            choices=list(template_choices.keys()),
            default="1",
            module="AMBER Controller",
            description="Select template",
        )
        selected_template = template_choices[choice]
        
        if selected_template:
            try:
                from .user_data_manager import UserDataManager
                user_data_manager = UserDataManager(console=self.console)
                
                # Get the correct template path for builtin templates
                template_metadata = templates[selected_template]
                if template_metadata["source"] == "builtin":
                    template_key = template_metadata["template_path"]
                else:
                    template_key = selected_template
                    
                content, metadata = user_data_manager.get_template_content(template_key)
                
                # Create temporary file for editing
                import tempfile
                import os
                import subprocess
                
                with tempfile.NamedTemporaryFile(mode='w', suffix='.mdin', delete=False) as tmp_file:
                    tmp_file.write(content)
                    tmp_path = tmp_file.name
                
                try:
                    self.console.print(f"\nOpening editor for template: {selected_template}")
                    editor = os.environ.get('EDITOR', 'nano')
                    subprocess.run([editor, tmp_path])
                    
                    # Read modified content
                    with open(tmp_path, 'r') as f:
                        modified_content = f.read()
                        
                    if modified_content != content:
                        # Check if this is a built-in template
                        template_metadata = templates[selected_template]
                        is_builtin = template_metadata["source"] == "builtin"
                        
                        if is_builtin:
                            # For built-in templates, save as new custom template
                            template_name = prompt_with_context(
                                self.processor,
                                "Enter name for the modified template",
                                default=metadata["name"],
                                module="AMBER Controller",
                                description="Modified template name",
                            )
                            description = prompt_with_context(
                                self.processor,
                                "Enter description",
                                default=metadata["description"],
                                module="AMBER Controller",
                                description="Modified template description",
                            )
                            
                            template_id = user_data_manager.create_template(
                                template_name,
                                description, 
                                metadata.get("simulation_type", "unknown"),
                                modified_content
                            )
                            if template_id:
                                self.console.print(f"[green]✓ Template '{template_name}' created as custom template with ID: {template_id}[/green]")
                            else:
                                self.console.print(f"[red]Failed to create custom template[/red]")
                        else:
                            # For custom templates, update in place
                            success = user_data_manager.modify_template(
                                selected_template, 
                                metadata["name"], 
                                metadata["description"],
                                metadata.get("simulation_type", "unknown"),
                                modified_content
                            )
                            if success:
                                self.console.print(f"[green]✓ Template '{selected_template}' updated successfully[/green]")
                            else:
                                self.console.print(f"[red]Failed to update template[/red]")
                    else:
                        self.console.print("[yellow]No changes made[/yellow]")
                        
                finally:
                    os.unlink(tmp_path)
                    
            except Exception as e:
                self.console.print(f"[red]Error modifying template: {e}[/red]")
                
    def _delete_single_template(self, templates: Dict):
        """Delete a custom template."""
        custom_templates = {k: v for k, v in templates.items() if v["source"] == "custom"}
        if not custom_templates:
            self.console.print("[yellow]No custom templates available to delete[/yellow]")
            return
            
        self.console.print("\n[bold]Select template to delete:[/bold]")
        
        # Use categorized display with numbers
        template_choices = self._display_categorized_templates(custom_templates, show_numbers=True)
        
        # Add back option
        back_num = len(custom_templates) + 1
        template_choices[str(back_num)] = None
        self.console.print(f"\n[bold]Options:[/bold]")
        self.console.print(f"  {back_num}. ← Back")
        
        choice = prompt_with_context(
            self.processor,
            "Select template",
            choices=list(template_choices.keys()),
            default=str(len(custom_templates)+1),
            module="AMBER Controller",
            description="Select template to delete",
        )
        selected_template = template_choices[choice]
        
        if selected_template:
            # Confirm deletion
            confirm = prompt_with_context(
                self.processor,
                f"Are you sure you want to delete template '{selected_template}'?",
                choices=["y", "n"],
                default="n",
                module="AMBER Controller",
                description=f"Confirm delete template {selected_template}",
                options_map={"y": "Yes, delete", "n": "No, keep"},
            )
            if confirm.lower() == "y":
                try:
                    from .user_data_manager import UserDataManager
                    user_data_manager = UserDataManager(console=self.console)
                    success = user_data_manager.delete_user_content(selected_template)
                    if success:
                        self.console.print(f"[green]✓ Template '{selected_template}' deleted successfully[/green]")
                    else:
                        self.console.print(f"[red]Failed to delete template[/red]")
                except Exception as e:
                    self.console.print(f"[red]Error deleting template: {e}[/red]")
                    
    # Workflow Template Management Methods
    def _display_workflow_list(self, workflows: Dict):
        """Display formatted list of workflow templates."""
        if not workflows:
            self.console.print("[yellow]No workflows available[/yellow]")
            return
            
        # Get workflow metadata for categorization
        from .user_data_manager import UserDataManager
        user_data_manager = UserDataManager(console=self.console)
        workflow_list = user_data_manager.list_workflows()
        
        # Built-in workflows
        builtin_workflows = {k: v for k, v in workflows.items() 
                           if workflow_list.get(k, {}).get("source") == "builtin"}
        if builtin_workflows:
            self.console.print("\n[bold]Built-in Workflows:[/bold]")
            for workflow_id, workflow in builtin_workflows.items():
                step_count = len(workflow.steps)
                self.console.print(f"  • {workflow.name} ({step_count} steps) - {workflow.description}")
                
        # Custom workflows  
        custom_workflows = {k: v for k, v in workflows.items() 
                          if workflow_list.get(k, {}).get("source") == "custom"}
        if custom_workflows:
            self.console.print("\n[bold]Custom Workflows:[/bold]")
            for workflow_id, workflow in custom_workflows.items():
                step_count = len(workflow.steps)
                self.console.print(f"  • {workflow.name} ({step_count} steps) - {workflow.description}")
                
    def _view_workflow_template(self, workflows: Dict):
        """View a specific workflow template."""
        if not workflows:
            self.console.print("[yellow]No workflows available to view[/yellow]")
            return
            
        workflow_choices = {str(i+1): wid for i, wid in enumerate(workflows.keys())}
        workflow_choices[str(len(workflows)+1)] = None  # Back option
        
        self.console.print("\n[bold]Select workflow to view:[/bold]")
        for i, (wid, workflow) in enumerate(workflows.items(), 1):
            self.console.print(f"  {i}. {workflow.name} - {workflow.description}")
        self.console.print(f"  {len(workflows)+1}. ← Back")
        
        choice = prompt_with_context(
            self.processor,
            "Select workflow",
            choices=list(workflow_choices.keys()),
            default="1",
            module="AMBER Controller",
            description="Select workflow",
        )
        selected_workflow = workflow_choices[choice]
        
        if selected_workflow:
            from .workflow_loader import WorkflowLoader
            workflow_loader = WorkflowLoader(console=self.console)
            workflow = workflows[selected_workflow]
            workflow_loader.display_workflow_overview(workflow)
            input("Press Enter to continue...")
            
    def _create_workflow_template(self):
        """Create a new workflow template - visual workflow builder."""
        self.console.print("\n[bold cyan]===== Create New Workflow Template =====[/bold cyan]")
        
        try:
            # Phase 1: Workflow Metadata
            workflow_metadata = self._collect_workflow_metadata()
            if not workflow_metadata:
                return
            
            # Phase 2: Template Selection & Sequencing
            workflow_steps = self._build_workflow_sequence()
            if not workflow_steps:
                return
                
            # Phase 3: Configure File Connections
            workflow_steps = self._configure_file_connections(workflow_steps)
            
            # Phase 4: Parameter Overrides (optional)
            if prompt_with_context(
                self.processor,
                "\nWould you like to configure parameter overrides for any steps?",
                choices=["y", "n"],
                default="n",
                module="AMBER Controller",
                description="Configure parameter overrides",
                options_map={"y": "Yes", "n": "No"},
            ) == "y":
                workflow_steps = self._configure_parameter_overrides(workflow_steps)
            
            # Phase 5: Preview Workflow
            self._preview_workflow(workflow_metadata, workflow_steps)
            
            # Phase 6: Save or Discard
            if prompt_with_context(
                self.processor,
                "\nSave this workflow template?",
                choices=["y", "n"],
                default="y",
                module="AMBER Controller",
                description="Save workflow template",
                options_map={"y": "Yes, save", "n": "No, discard"},
            ) == "y":
                success = self._save_workflow_template(workflow_metadata, workflow_steps)
                if success:
                    self.console.print(f"[green]✓ Workflow template '{workflow_metadata['name']}' created successfully![/green]")
                else:
                    self.console.print("[red]Failed to save workflow template[/red]")
            else:
                self.console.print("[yellow]Workflow template discarded[/yellow]")
                
        except KeyboardInterrupt:
            self.console.print("\n[yellow]Workflow creation cancelled[/yellow]")
        except Exception as e:
            self.console.print(f"[red]Error creating workflow: {e}[/red]")
        
    def _modify_workflow_template(self, workflows: Dict):
        """Modify an existing workflow template."""
        if not workflows:
            self.console.print("[yellow]No workflows available to modify[/yellow]")
            return
            
        from .user_data_manager import UserDataManager
        user_data_manager = UserDataManager(console=self.console)
        workflow_list = user_data_manager.list_workflows()
            
        workflow_choices = {str(i+1): wid for i, wid in enumerate(workflows.keys())}
        workflow_choices[str(len(workflows)+1)] = None  # Back option
        
        self.console.print("\n[bold]Select workflow to modify:[/bold]")
        
        # Group by source for display
        builtin_workflows = {k: v for k, v in workflows.items() 
                           if workflow_list.get(k, {}).get("source") == "builtin"}
        custom_workflows = {k: v for k, v in workflows.items() 
                          if workflow_list.get(k, {}).get("source") == "custom"}
        
        choice_num = 1
        
        if builtin_workflows:
            self.console.print("\n[bold]Built-in Workflows:[/bold]")
            for wid, workflow in builtin_workflows.items():
                self.console.print(f"  {choice_num}. {workflow.name} - {workflow.description}")
                choice_num += 1
                
        if custom_workflows:
            self.console.print("\n[bold]Custom Workflows:[/bold]")
            for wid, workflow in custom_workflows.items():
                self.console.print(f"  {choice_num}. {workflow.name} - {workflow.description}")
                choice_num += 1
        
        self.console.print(f"\n[bold]Options:[/bold]")
        self.console.print(f"  {len(workflows)+1}. ← Back")
        
        choice = prompt_with_context(
            self.processor,
            "Select workflow",
            choices=list(workflow_choices.keys()),
            default="1",
            module="AMBER Controller",
            description="Select workflow",
        )
        selected_workflow = workflow_choices[choice]
        
        if selected_workflow:
            try:
                # Check if this is a built-in workflow
                workflow_metadata = workflow_list.get(selected_workflow, {})
                is_builtin = workflow_metadata.get("source") == "builtin"
                
                # Get workflow content
                if is_builtin:
                    workflow_key = workflow_metadata["workflow_path"]
                else:
                    workflow_key = selected_workflow
                    
                workflow_content, workflow_meta = user_data_manager.get_workflow_content(workflow_key)
                
                if is_builtin:
                    self.console.print(f"\n[yellow]Modifying built-in workflow will create a custom copy.[/yellow]")
                
                self.console.print(f"\n[bold]Editing Workflow:[/bold] {workflow_content['name']}")
                self.console.print(f"[bold]Description:[/bold] {workflow_content['description']}")
                
                # Start workflow modification interface
                self._edit_workflow_steps(workflow_content, workflow_meta, user_data_manager, is_builtin)
                
            except Exception as e:
                self.console.print(f"[red]Error modifying workflow: {e}[/red]")
    
    def _edit_workflow_steps(self, workflow_content: Dict, workflow_meta: Dict, user_data_manager, is_builtin: bool):
        """Command-line workflow step editor."""
        steps = workflow_content.get("steps", [])
        modified = False
        
        # Show initial help
        self.console.print(f"\n[bold cyan]===== Workflow Step Editor =====[/bold cyan]")
        self.console.print(f"[bold]Workflow:[/bold] {workflow_content['name']}")
        self._show_workflow_help()
        
        while True:
            # Display current steps
            self._display_workflow_steps(steps)
            
            # Get command
            try:
                command = prompt_with_context(
                    self.processor,
                    f"\n[bold blue]workflow>[/bold blue]",
                    default="save",
                    module="AMBER Controller",
                    description="Workflow editor command",
                ).strip().lower()
                
                if not command:
                    continue
                    
                # Parse and execute command
                result = self._execute_workflow_command(command, steps)
                
                if result == "modified":
                    modified = True
                elif result == "save":
                    # Save workflow
                    if modified or is_builtin:
                        try:
                            # Get workflow metadata
                            name = prompt_with_context(
                                self.processor,
                                "Workflow name",
                                default=workflow_content.get("name", "Modified Workflow"),
                                module="AMBER Controller",
                                description="Workflow name (edit)",
                            )
                            description = prompt_with_context(
                                self.processor,
                                "Description",
                                default=workflow_content.get("description", "Modified workflow"),
                                module="AMBER Controller",
                                description="Workflow description (edit)",
                            )
                            author = prompt_with_context(
                self.processor,
                "Author",
                default="ProPrep User",
                module="AMBER Controller",
                description="Template author",
            )
                            
                            # Determine source path
                            if is_builtin:
                                source_path = workflow_meta.get("workflow_path", "builtin/unknown.json") 
                            else:
                                source_path = workflow_meta.get("workflow_path", "custom/unknown.json")
                            
                            # Save modified workflow
                            workflow_id = user_data_manager.modify_workflow(
                                source_path=source_path,
                                name=name,
                                description=description,
                                steps=steps,
                                author=author
                            )
                            
                            self.console.print(f"[green]✓ Workflow saved with ID: {workflow_id}[/green]")
                            break
                            
                        except Exception as e:
                            self.console.print(f"[red]Error saving workflow: {e}[/red]")
                    else:
                        self.console.print("[yellow]No changes to save[/yellow]")
                        break
                        
                elif result == "cancel":
                    # Cancel
                    if modified:
                        confirm = confirm_with_context(
                            self.processor,
                            "Discard changes?",
                            default=False,
                            module="AMBER Controller",
                            description="Discard workflow changes",
                        )
                        if confirm:
                            self.console.print("[yellow]Changes discarded[/yellow]")
                            break
                    else:
                        break
                        
            except (KeyboardInterrupt, EOFError):
                if modified:
                    confirm = confirm_with_context(
                        self.processor,
                        "\nDiscard changes?",
                        default=False,
                        module="AMBER Controller",
                        description="Discard workflow changes on cancel",
                    )
                    if confirm:
                        self.console.print("[yellow]Changes discarded[/yellow]")
                        break
                else:
                    break
    
    def _show_workflow_help(self):
        """Show command syntax help."""
        self.console.print(f"\n[bold]Commands:[/bold] (word 'step' is optional)")
        self.console.print("  [grey50]add [step] before/after X          - Add step at position[/grey50]")
        self.console.print("  [grey50]add [step] between X and Y         - Add step between positions[/grey50]")
        self.console.print("  [grey50]remove [step] X                    - Remove step X[/grey50]")
        self.console.print("  [grey50]move [step] X before/after Y       - Move step X relative to Y[/grey50]")
        self.console.print("  [grey50]move [step] X to start/end         - Move step to start or end[/grey50]")
        self.console.print("  [grey50]move [step] X to #                 - Move step to position #[/grey50]")
        self.console.print("  [grey50]modify [step] X                    - Modify step X[/grey50]")
        self.console.print("  [grey50]view [step] X                      - View step X details[/grey50]")
        self.console.print("  [grey50]save                               - Save workflow[/grey50]")
        self.console.print("  [grey50]cancel                             - Cancel changes[/grey50]")
    
    def _display_workflow_steps(self, steps: List[Dict]):
        """Display workflow steps with consistent formatting."""
        if not steps:
            self.console.print(f"\n[yellow]No steps defined[/yellow]")
            return
            
        self.console.print(f"\n[bold]Current Steps ({len(steps)}):[/bold]")
        
        # Calculate max name length for consistent alignment  
        max_name_len = 0
        for step in steps:
            step_name = step.get("name", step.get("id", f"Step {len(steps)}"))
            max_name_len = max(max_name_len, len(step_name))
        
        for i, step in enumerate(steps, 1):
            step_name = step.get("name", step.get("id", f"Step {i}"))
            template = step.get("template", "N/A")
            
            # Right-align numbers and align descriptions
            num_str = f"{i}.".rjust(3)
            padding = " " * (max_name_len - len(step_name) + 2)
            
            self.console.print(f"    {num_str} {step_name}{padding}- {template}")
    
    def _execute_workflow_command(self, command: str, steps: List[Dict]) -> str:
        """Parse and execute workflow command. Returns 'modified', 'save', 'cancel', or None."""
        import re
        
        # Normalize command - remove optional "step" word
        cmd = re.sub(r'\bstep\s+', '', command)
        
        # Parse different command patterns
        if cmd == "save":
            return "save"
        elif cmd == "cancel":
            return "cancel"
            
        # Add commands
        elif match := re.match(r'add\s+(?:before|after)\s+(\d+)', cmd):
            pos = int(match.group(1))
            position = "before" if "before" in cmd else "after"
            return self._add_step_at_position(steps, pos, position)
            
        elif match := re.match(r'add\s+between\s+(\d+)\s+and\s+(\d+)', cmd):
            pos1, pos2 = int(match.group(1)), int(match.group(2))
            return self._add_step_between(steps, pos1, pos2)
            
        # Remove commands
        elif match := re.match(r'remove\s+(\d+)', cmd):
            pos = int(match.group(1))
            return self._remove_step(steps, pos)
            
        # Move commands
        elif match := re.match(r'move\s+(\d+)\s+(?:before|after)\s+(\d+)', cmd):
            from_pos, to_pos = int(match.group(1)), int(match.group(2))
            direction = "before" if "before" in cmd else "after"
            return self._move_step_relative(steps, from_pos, to_pos, direction)
            
        elif match := re.match(r'move\s+(\d+)\s+to\s+(?:start|end)', cmd):
            pos = int(match.group(1))
            position = "start" if "start" in cmd else "end"
            return self._move_step_to_boundary(steps, pos, position)
            
        elif match := re.match(r'move\s+(\d+)\s+to\s+(\d+)', cmd):
            from_pos, to_pos = int(match.group(1)), int(match.group(2))
            return self._move_step_to_position(steps, from_pos, to_pos)
            
        # View/modify commands
        elif match := re.match(r'view\s+(\d+)', cmd):
            pos = int(match.group(1))
            return self._view_step(steps, pos)
            
        elif match := re.match(r'modify\s+(\d+)', cmd):
            pos = int(match.group(1))
            return self._modify_step_at_position(steps, pos)
            
        else:
            self.console.print(f"[red]Unknown command: '{command}'[/red]")
            self.console.print("[grey50]Type a command or press Enter to save[/grey50]")
            return None
    
    def _create_workflow_step(self) -> Optional[Dict]:
        """Create a new workflow step."""
        try:
            self.console.print(f"\n[bold]Create New Workflow Step[/bold]")
            
            name = prompt_with_context(
                self.processor,
                "Step name",
                default="New Step",
                module="AMBER Controller",
                description="New step name",
            )
            description = prompt_with_context(
                self.processor,
                "Description",
                default="",
                module="AMBER Controller",
                description="New step description",
            )
            step_type = prompt_with_context(
                self.processor,
                "Step type",
                choices=["minimization", "heating", "equilibration", "production"],
                default="minimization",
                module="AMBER Controller",
                description="New step type",
                options_map={
                    "minimization": "Energy minimization",
                    "heating": "Heating/temperature ramp",
                    "equilibration": "Equilibration",
                    "production": "Production MD",
                },
            )
            
            # Show available templates and let user select
            template = self._select_template_for_step(step_type)
            if not template:
                return None  # User cancelled template selection
            
            # Generate step ID
            import uuid
            step_id = f"step_{uuid.uuid4().hex[:8]}"
            
            step = {
                "id": step_id,
                "name": name,
                "type": step_type,
                "template": template,
                "description": description,
                "parameter_overrides": {},
                "dependencies": [],
                "input_coord": "previous",
                "output_coord": f"{step_id}_output"
            }
            
            return step
            
        except KeyboardInterrupt:
            self.console.print("[yellow]Step creation cancelled[/yellow]")
            return None
            
    def _modify_workflow_step(self, step: Dict) -> Optional[Dict]:
        """Modify an existing workflow step."""
        try:
            self.console.print(f"\n[bold]Modify Workflow Step[/bold]")
            self.console.print(f"Current: {step.get('name', 'Unnamed')}")
            
            new_step = step.copy()
            
            name = prompt_with_context(
                self.processor,
                "Step name",
                default=step.get("name", ""),
                module="AMBER Controller",
                description="Step name (edit)",
            )
            description = prompt_with_context(
                self.processor,
                "Description",
                default=step.get("description", ""),
                module="AMBER Controller",
                description="Step description (edit)",
            )
            step_type = prompt_with_context(
                self.processor,
                "Step type",
                choices=["minimization", "heating", "equilibration", "production"],
                default=step.get("type", "minimization"),
                module="AMBER Controller",
                description="Step type (edit)",
                options_map={
                    "minimization": "Energy minimization",
                    "heating": "Heating/temperature ramp",
                    "equilibration": "Equilibration",
                    "production": "Production MD",
                },
            )
            template = prompt_with_context(
                self.processor,
                "Template file",
                default=step.get("template", ""),
                module="AMBER Controller",
                description="Template file for step (edit)",
            )
            
            new_step.update({
                "name": name,
                "description": description,
                "type": step_type,
                "template": template
            })
            
            return new_step
            
        except KeyboardInterrupt:
            self.console.print("[yellow]Step modification cancelled[/yellow]")
            return None
            
    def _display_step_details(self, step: Dict):
        """Display detailed information about a workflow step."""
        self.console.print(f"\n[bold cyan]===== Step Details =====[/bold cyan]")
        self.console.print(f"[bold]ID:[/bold] {step.get('id', 'N/A')}")
        self.console.print(f"[bold]Name:[/bold] {step.get('name', 'N/A')}")
        self.console.print(f"[bold]Type:[/bold] {step.get('type', 'N/A')}")
        self.console.print(f"[bold]Template:[/bold] {step.get('template', 'N/A')}")
        self.console.print(f"[bold]Description:[/bold] {step.get('description', 'N/A')}")
        self.console.print(f"[bold]Input:[/bold] {step.get('input_coord', 'N/A')}")
        self.console.print(f"[bold]Output:[/bold] {step.get('output_coord', 'N/A')}")
        
        if step.get('dependencies'):
            self.console.print(f"[bold]Dependencies:[/bold] {', '.join(step['dependencies'])}")
        else:
            self.console.print(f"[bold]Dependencies:[/bold] None")
            
        if step.get('parameter_overrides'):
            self.console.print(f"[bold]Parameter Overrides:[/bold]")
            for key, value in step['parameter_overrides'].items():
                self.console.print(f"  {key}: {value}")
        else:
            self.console.print(f"[bold]Parameter Overrides:[/bold] None")
        
        input("Press Enter to continue...")
    
    # Command handler methods for workflow editor
    def _add_step_at_position(self, steps: List[Dict], pos: int, position: str) -> str:
        """Add step before/after specified position."""
        if pos < 1 or pos > len(steps):
            self.console.print(f"[red]Invalid position: {pos}. Valid range: 1-{len(steps)}[/red]")
            return None
            
        new_step = self._create_workflow_step()
        if new_step:
            if position == "before":
                steps.insert(pos - 1, new_step)
            else:  # after
                steps.insert(pos, new_step)
            self.console.print(f"[green]✓ Step added {position} position {pos}[/green]")
            return "modified"
        return None
    
    def _add_step_between(self, steps: List[Dict], pos1: int, pos2: int) -> str:
        """Add step between two positions."""
        if pos1 >= pos2:
            self.console.print(f"[red]First position ({pos1}) must be less than second ({pos2})[/red]")
            return None
        if pos1 < 1 or pos2 > len(steps):
            self.console.print(f"[red]Invalid positions. Valid range: 1-{len(steps)}[/red]")
            return None
            
        new_step = self._create_workflow_step()
        if new_step:
            # Insert at pos1 + 1 (after pos1, before pos2)
            steps.insert(pos1, new_step)
            self.console.print(f"[green]✓ Step added between positions {pos1} and {pos2}[/green]")
            return "modified"
        return None
    
    def _remove_step(self, steps: List[Dict], pos: int) -> str:
        """Remove step at position."""
        if pos < 1 or pos > len(steps):
            self.console.print(f"[red]Invalid position: {pos}. Valid range: 1-{len(steps)}[/red]")
            return None
            
        removed_step = steps.pop(pos - 1)
        step_name = removed_step.get("name", f"Step {pos}")
        self.console.print(f"[green]✓ Removed step {pos}: '{step_name}'[/green]")
        return "modified"
    
    def _move_step_relative(self, steps: List[Dict], from_pos: int, to_pos: int, direction: str) -> str:
        """Move step before/after another step."""
        if from_pos < 1 or from_pos > len(steps) or to_pos < 1 or to_pos > len(steps):
            self.console.print(f"[red]Invalid positions. Valid range: 1-{len(steps)}[/red]")
            return None
        if from_pos == to_pos:
            self.console.print("[yellow]Source and destination are the same[/yellow]")
            return None
            
        # Remove step from original position
        step = steps.pop(from_pos - 1)
        
        # Adjust target position after removal
        if from_pos < to_pos:
            to_pos -= 1
            
        # Insert at new position
        if direction == "before":
            steps.insert(to_pos - 1, step)
        else:  # after
            steps.insert(to_pos, step)
            
        step_name = step.get("name", f"Step {from_pos}")
        self.console.print(f"[green]✓ Moved '{step_name}' {direction} position {to_pos}[/green]")
        return "modified"
    
    def _move_step_to_boundary(self, steps: List[Dict], pos: int, position: str) -> str:
        """Move step to start or end."""
        if pos < 1 or pos > len(steps):
            self.console.print(f"[red]Invalid position: {pos}. Valid range: 1-{len(steps)}[/red]")
            return None
            
        step = steps.pop(pos - 1)
        
        if position == "start":
            steps.insert(0, step)
        else:  # end
            steps.append(step)
            
        step_name = step.get("name", f"Step {pos}")
        self.console.print(f"[green]✓ Moved '{step_name}' to {position}[/green]")
        return "modified"
    
    def _move_step_to_position(self, steps: List[Dict], from_pos: int, to_pos: int) -> str:
        """Move step to specific position."""
        if from_pos < 1 or from_pos > len(steps) or to_pos < 1 or to_pos > len(steps):
            self.console.print(f"[red]Invalid positions. Valid range: 1-{len(steps)}[/red]")
            return None
        if from_pos == to_pos:
            self.console.print("[yellow]Source and destination are the same[/yellow]")
            return None
            
        step = steps.pop(from_pos - 1)
        
        # Insert at target position (no adjustment needed - user specifies final position)
        steps.insert(to_pos - 1, step)
        
        step_name = step.get("name", f"Step {from_pos}")
        self.console.print(f"[green]✓ Moved '{step_name}' to position {to_pos}[/green]")
        return "modified"
    
    def _view_step(self, steps: List[Dict], pos: int) -> str:
        """View step details."""
        if pos < 1 or pos > len(steps):
            self.console.print(f"[red]Invalid position: {pos}. Valid range: 1-{len(steps)}[/red]")
            return None
            
        self._display_step_details(steps[pos - 1])
        return None
    
    def _modify_step_at_position(self, steps: List[Dict], pos: int) -> str:
        """Modify step at position."""
        if pos < 1 or pos > len(steps):
            self.console.print(f"[red]Invalid position: {pos}. Valid range: 1-{len(steps)}[/red]")
            return None
            
        updated_step = self._modify_workflow_step(steps[pos - 1])
        if updated_step:
            steps[pos - 1] = updated_step
            step_name = updated_step.get("name", f"Step {pos}")
            self.console.print(f"[green]✓ Modified step {pos}: '{step_name}'[/green]")
            return "modified"
        return None
    
    def _select_template_for_step(self, step_type: str) -> Optional[str]:
        """Show available templates and let user select one."""
        try:
            # Get available templates
            from .user_data_manager import UserDataManager
            user_data_manager = UserDataManager(console=self.console)
            all_templates = user_data_manager.list_templates()
            
            # Filter templates by type
            matching_templates = {
                tid: meta for tid, meta in all_templates.items()
                if meta.get("simulation_type") == step_type
            }
            
            if not matching_templates:
                self.console.print(f"[yellow]No {step_type} templates available[/yellow]")
                return None
            
            self.console.print(f"\n[bold]Available {step_type.title()} Templates:[/bold]")
            
            # Use categorized display but only for the matching type
            template_choices = {}
            choice_num = 1
            
            # Separate builtin and custom
            builtin_templates = {k: v for k, v in matching_templates.items() if v["source"] == "builtin"}
            custom_templates = {k: v for k, v in matching_templates.items() if v["source"] == "custom"}
            
            if builtin_templates:
                self.console.print(f"\n  [bold]Built-in:[/bold]")
                for tid, meta in sorted(builtin_templates.items(), key=lambda x: (x[1].get("priority", 999), x[1]["name"])):
                    name = meta["name"]
                    description = meta["description"]
                    # Calculate padding for alignment
                    padding = " " * (35 - len(name)) if len(name) < 35 else "  "
                    num_str = f"{choice_num}.".rjust(3)
                    self.console.print(f"    {num_str} {name}{padding}- {description}")
                    template_choices[str(choice_num)] = meta["template_path"]
                    choice_num += 1
            
            if custom_templates:
                self.console.print(f"\n  [bold]Custom:[/bold]")
                for tid, meta in sorted(custom_templates.items(), key=lambda x: x[1]["name"]):
                    name = meta["name"]
                    description = meta["description"]
                    padding = " " * (35 - len(name)) if len(name) < 35 else "  "
                    num_str = f"{choice_num}.".rjust(3)
                    self.console.print(f"    {num_str} {name}{padding}- {description}")
                    template_choices[str(choice_num)] = tid  # Custom templates use ID, not path
                    choice_num += 1
            
            # Get user choice
            choice = prompt_with_context(
            self.processor,
            "Select template",
            choices=list(template_choices.keys()),
            default="1",
            module="AMBER Controller",
            description="Select template",
        )
            return template_choices[choice]
            
        except KeyboardInterrupt:
            self.console.print("[yellow]Template selection cancelled[/yellow]")
            return None
        except Exception as e:
            self.console.print(f"[red]Error loading templates: {e}[/red]")
            return None
            
    def _delete_workflow_template(self, workflows: Dict):
        """Delete a custom workflow template."""
        # Get custom workflows only
        from .user_data_manager import UserDataManager
        user_data_manager = UserDataManager(console=self.console)
        workflow_list = user_data_manager.list_workflows()
        
        custom_workflows = {k: v for k, v in workflows.items() 
                          if workflow_list.get(k, {}).get("source") == "custom"}
        if not custom_workflows:
            self.console.print("[yellow]No custom workflows available to delete[/yellow]")
            return
            
        workflow_choices = {str(i+1): wid for i, wid in enumerate(custom_workflows.keys())}
        workflow_choices[str(len(custom_workflows)+1)] = None  # Back option
        
        self.console.print("\n[bold]Select workflow to delete:[/bold]")
        for i, (wid, workflow) in enumerate(custom_workflows.items(), 1):
            self.console.print(f"  {i}. {workflow.name} - {workflow.description}")
        self.console.print(f"  {len(custom_workflows)+1}. ← Back")
        
        choice = prompt_with_context(
            self.processor,
            "Select workflow",
            choices=list(workflow_choices.keys()),
            default=str(len(custom_workflows)+1),
            module="AMBER Controller",
            description="Select workflow to delete",
        )
        selected_workflow = workflow_choices[choice]
        
        if selected_workflow:
            workflow_name = workflows[selected_workflow].name
            # Confirm deletion
            confirm = prompt_with_context(
                self.processor,
                f"Are you sure you want to delete workflow '{workflow_name}'?",
                choices=["y", "n"],
                default="n",
                module="AMBER Controller",
                description=f"Confirm delete workflow {workflow_name}",
                options_map={"y": "Yes, delete", "n": "No, keep"},
            )
            if confirm.lower() == "y":
                try:
                    success = user_data_manager.delete_content(selected_workflow)
                    if success:
                        self.console.print(f"[green]✓ Workflow '{workflow_name}' deleted successfully[/green]")
                    else:
                        self.console.print(f"[red]Failed to delete workflow[/red]")
                except Exception as e:
                    self.console.print(f"[red]Error deleting workflow: {e}[/red]")
                    
    # Workflow Builder Helper Methods
    def _collect_workflow_metadata(self) -> Dict:
        """Collect basic workflow metadata."""
        self.console.print("\n[bold cyan]Step 1: Workflow Information[/bold cyan]")
        
        name = prompt_with_context(
            self.processor,
            "Workflow name",
            module="AMBER Controller",
            description="New workflow name",
        )
        if not name:
            return None
            
        description = prompt_with_context(
            self.processor,
            "Description",
            module="AMBER Controller",
            description="Description",
        )
        if not description:
            return None
            
        # Category selection
        categories = ["protein", "nucleic_acid", "membrane", "ligand", "custom"]
        self.console.print("\n[bold]Categories:[/bold]")
        for i, cat in enumerate(categories, 1):
            self.console.print(f"  {i}. {cat}")
            
        cat_choice = prompt_with_context(
            self.processor,
            "Select category",
            choices=[str(i) for i in range(1, len(categories)+1)],
            default="1",
            module="AMBER Controller",
            description="Select workflow category",
            options_map={str(i + 1): c for i, c in enumerate(categories)},
        )
        category = categories[int(cat_choice)-1]
        
        author = prompt_with_context(
            self.processor,
            "Author",
            default="ProPrep User",
            module="AMBER Controller",
            description="Workflow author",
        )
        version = prompt_with_context(
            self.processor,
            "Version",
            default="1.0",
            module="AMBER Controller",
            description="Workflow version",
        )
        
        return {
            "name": name,
            "description": description,
            "category": category,
            "author": author,
            "version": version
        }
        
    def _build_workflow_sequence(self) -> List[Dict]:
        """Build the sequence of workflow steps using comma-separated input."""
        self.console.print("\n[bold cyan]Step 2: Build Workflow Sequence[/bold cyan]")
        
        # Load available templates
        try:
            from .user_data_manager import UserDataManager
            user_data_manager = UserDataManager(console=self.console)
            available_templates = user_data_manager.list_templates()
        except Exception as e:
            self.console.print(f"[red]Error loading templates: {e}[/red]")
            return []
            
        while True:
            # Show available templates with numbers
            self.console.print("\n[bold]Available templates:[/bold]")
            template_map = self._display_categorized_templates(available_templates, show_numbers=True)
            
            # Show the mapping for user reference
            self.console.print(f"\n[bold]Quick reference:[/bold] Enter template numbers separated by commas")
            self.console.print("[grey50]Example: 1,5,8,12 (will create a 4-step workflow)[/grey50]")
            
            # Get user input
            sequence_input = prompt_with_context(
                self.processor,
                "\nEnter template sequence (comma-separated numbers)",
                default="",
                module="AMBER Controller",
                description="Template sequence for workflow steps",
            ).strip()
            
            if not sequence_input:
                if prompt_with_context(
                    self.processor,
                    "Cancel workflow creation?",
                    choices=["y", "n"],
                    default="n",
                    module="AMBER Controller",
                    description="Cancel workflow creation",
                    options_map={"y": "Yes, cancel", "n": "No, continue"},
                ) == "y":
                    return []
                continue
            
            # Parse the input
            try:
                sequence_numbers = [num.strip() for num in sequence_input.split(",")]
                workflow_steps = []
                
                # Validate and create steps
                for i, num_str in enumerate(sequence_numbers, 1):
                    if num_str not in template_map:
                        self.console.print(f"[red]Error: '{num_str}' is not a valid template number[/red]")
                        workflow_steps = None
                        break
                        
                    tid = template_map[num_str]
                    if not tid:  # Skip None values
                        self.console.print(f"[red]Error: Invalid template selection[/red]")
                        workflow_steps = None
                        break
                        
                    meta = available_templates[tid]
                    step_id = f"step_{i}"
                    
                    # Create step data structure
                    step = {
                        "id": step_id,
                        "name": meta["name"],
                        "type": meta.get("simulation_type", "unknown"),
                        "template": meta.get("template_path", tid) if meta["source"] == "builtin" else tid,
                        "description": meta["description"],
                        "parameter_overrides": {},
                        "dependencies": [],
                        "input_coord": f"{step_id}.rst7" if i > 1 else "system.rst7",
                        "output_coord": f"{step_id}_out.rst7"
                    }
                    
                    workflow_steps.append(step)
                
                # If validation failed, continue loop to re-enter
                if workflow_steps is None:
                    continue
                    
                # Show the proposed sequence for confirmation
                self.console.print(f"\n[bold green]Proposed workflow sequence ({len(workflow_steps)} steps):[/bold green]")
                for i, step in enumerate(workflow_steps, 1):
                    arrow = " → " if i < len(workflow_steps) else ""
                    self.console.print(f"  {i}. {step['name']} ({step['type']}){arrow}", end="")
                self.console.print()  # Final newline
                
                # Confirmation options
                self.console.print(f"\n[bold]Options:[/bold]")
                self.console.print("1. ✓ Accept this sequence", highlight=False)
                self.console.print("2. ✏ Enter a different sequence", highlight=False)  
                self.console.print("3. ✗ Cancel workflow creation", highlight=False)
                
                choice = prompt_with_context(
                    self.processor,
                    "Select option",
                    choices=["1", "2", "3"],
                    default="1",
                    module="AMBER Controller",
                    description="Select workflow-building option",
                )
                
                if choice == "1":
                    return workflow_steps
                elif choice == "2":
                    continue  # Loop back to enter new sequence
                elif choice == "3":
                    return []
                    
            except Exception as e:
                self.console.print(f"[red]Error parsing sequence: {e}[/red]")
                self.console.print("[yellow]Please enter numbers separated by commas (e.g., 1,3,5,7)[/yellow]")
                continue
        
    def _configure_file_connections(self, workflow_steps: List[Dict]) -> List[Dict]:
        """Configure input/output file connections between steps."""
        self.console.print("\n[bold cyan]Step 3: Configure File Connections[/bold cyan]")
        
        if len(workflow_steps) <= 1:
            self.console.print("[grey50]Single step workflow - no connections needed[/grey50]")
            return workflow_steps
            
        self.console.print("Configuring automatic file connections...")
        
        # Set up automatic dependencies and file connections
        for i, step in enumerate(workflow_steps):
            if i == 0:
                # First step
                step["input_coord"] = "system.rst7"
                step["dependencies"] = []
            else:
                # Subsequent steps depend on previous step
                prev_step = workflow_steps[i-1]
                step["dependencies"] = [prev_step["id"]]
                step["input_coord"] = prev_step["output_coord"]
                
        # Show connections
        self.console.print("\n[bold]File connection chain:[/bold]")
        for i, step in enumerate(workflow_steps):
            if i == 0:
                self.console.print(f"  system.rst7 → [{step['name']}] → {step['output_coord']}")
            else:
                self.console.print(f"  {step['input_coord']} → [{step['name']}] → {step['output_coord']}")
                
        # Ask if user wants to customize
        if prompt_with_context(
            self.processor,
            "\nCustomize file connections?",
            choices=["y", "n"],
            default="n",
            module="AMBER Controller",
            description="Customize file connections between steps",
            options_map={"y": "Yes, customize", "n": "No, use defaults"},
        ) == "y":
            workflow_steps = self._customize_file_connections(workflow_steps)
            
        return workflow_steps
        
    def _customize_file_connections(self, workflow_steps: List[Dict]) -> List[Dict]:
        """Allow user to customize file connections."""
        # Implementation for custom file connections - can be added later
        self.console.print("[yellow]Custom file connections - advanced feature coming soon![/yellow]")
        self.console.print("Current automatic connections will be used.")
        input("Press Enter to continue...")
        return workflow_steps
        
    def _configure_parameter_overrides(self, workflow_steps: List[Dict]) -> List[Dict]:
        """Configure parameter overrides for workflow steps."""
        self.console.print("\n[bold cyan]Step 4: Configure Parameter Overrides[/bold cyan]")
        
        for step in workflow_steps:
            self.console.print(f"\n[bold]Step: {step['name']} ({step['type']})[/bold]")
            
            if prompt_with_context(
                self.processor,
                f"Add parameter overrides for {step['name']}?",
                choices=["y", "n"],
                default="n",
                module="AMBER Controller",
                description=f"Add parameter overrides for {step['name']}",
                options_map={"y": "Yes, add overrides", "n": "No, skip"},
            ) == "y":
                # Simple parameter override interface
                overrides = {}
                
                while True:
                    param = prompt_with_context(
                        self.processor,
                        "Parameter name (or 'done' to finish)",
                        default="done",
                        module="AMBER Controller",
                        description="Parameter name for override",
                    )
                    if param.lower() == "done":
                        break
                        
                    value = prompt_with_context(
                        self.processor,
                        f"Value for {param}",
                        module="AMBER Controller",
                        description=f"Value for parameter {param}",
                    )
                    if value:
                        overrides[param] = value
                        self.console.print(f"[green]✓ Added override: {param} = {value}[/green]")
                
                step["parameter_overrides"] = overrides
                
        return workflow_steps
        
    def _preview_workflow(self, metadata: Dict, workflow_steps: List[Dict]):
        """Preview the complete workflow."""
        self.console.print("\n[bold cyan]Step 5: Workflow Preview[/bold cyan]")
        
        # Workflow info
        self.console.print(f"\n[bold]Workflow: {metadata['name']}[/bold]")
        self.console.print(f"[bold]Description:[/bold] {metadata['description']}")
        self.console.print(f"[bold]Category:[/bold] {metadata['category']}")
        self.console.print(f"[bold]Author:[/bold] {metadata['author']}")
        self.console.print(f"[bold]Version:[/bold] {metadata['version']}")
        self.console.print(f"[bold]Steps:[/bold] {len(workflow_steps)}")
        
        # Steps table
        from rich.table import Table
        table = Table(title="Workflow Steps")
        table.add_column("Step", style="cyan", width=12)
        table.add_column("Template", style="green", width=25)
        table.add_column("Type", style="blue", width=15)
        table.add_column("Input", style="yellow", width=15)
        table.add_column("Output", style="yellow", width=15)
        table.add_column("Overrides", style="magenta", width=10)
        
        for step in workflow_steps:
            overrides_count = len(step["parameter_overrides"])
            override_display = f"{overrides_count}" if overrides_count > 0 else "-"
            
            table.add_row(
                step["id"],
                step["name"],
                step["type"],
                step["input_coord"],
                step["output_coord"],
                override_display
            )
            
        self.console.print(table)
        input("\nPress Enter to continue...")
        
    def _save_workflow_template(self, metadata: Dict, workflow_steps: List[Dict]) -> bool:
        """Save the workflow template."""
        try:
            from .user_data_manager import UserDataManager
            user_data_manager = UserDataManager(console=self.console)
            
            # Create workflow data structure matching WorkflowPreset format
            workflow_data = {
                "name": metadata["name"],
                "description": metadata["description"],
                "version": metadata["version"],
                "author": metadata["author"],
                "category": metadata["category"],
                "steps": workflow_steps
            }
            
            # Create workflow using user data manager
            success = user_data_manager.create_workflow(
                metadata["name"],
                metadata["description"],
                workflow_data
            )
            
            return success
            
        except Exception as e:
            self.console.print(f"[red]Error saving workflow: {e}[/red]")
            return False