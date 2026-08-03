"""
AMBER Input Generator Module - Redesigned Unified Architecture

Generates AMBER input files using a template-centric approach with wizard integration.
This redesigned system provides a clean separation between comprehensive parameter
configuration (wizard) and template management, with a unified user experience.

Version 2.0 - Complete architectural redesign
"""

import os
from typing import Dict, Any, Optional, List
from pathlib import Path

from proprep.utils.paths import get_package_dir
from proprep.utils.module_registry import ProcessingModule, register_module
from proprep.md_prep.amber_controller import AmberController


@register_module
class AmberInputGenerator(ProcessingModule):
    """
    Generate AMBER MD input files using unified template-wizard system.

    This redesigned system provides:
    - Template-centric workflow (always start/end with annotated templates)
    - Choice of configuration method (wizard vs direct edit)
    - Individual simulation phase configuration
    - Complete workflow assembly
    - Custom template and preset management

    Architecture:
    - AmberController: Orchestrates template-wizard interactions
    - AmberTemplateSystem: Manages rich, annotated templates
    - AmberWizard: Comprehensive 15-step parameter configuration
    """

    NAME = "AMBER Input Generator"
    CATEGORY = "MD Preparation"
    DESCRIPTION = "Create optimized mdin files using annotated templates and comprehensive wizard"
    VERSION = "2.0.0"

    def __init__(self):
        super().__init__()

        # Package directory for template storage
        self.package_dir = get_package_dir() / "md_prep"
        
        # Initialize the unified controller (lazy initialization)
        self.controller = None
        
    def set_processor(self, processor):
        """Set the processor reference and initialize controller."""
        self.processor = processor
        
        # Initialize controller with processor context
        self.controller = AmberController(processor, self.package_dir)

    @property
    def console(self):
        """Get console from processor if available."""
        if (
            hasattr(self, "processor")
            and self.processor
            and hasattr(self.processor, "console")
        ):
            return self.processor.console
        else:
            from rich.console import Console
            return Console()

    def get_menu_options(self) -> Dict[str, str]:
        """Get redesigned menu options for simulation-phase-based workflow."""
        if self.controller:
            return self.controller.get_menu_options()
        else:
            # Fallback menu before controller is initialized
            return {
                "minimization": "Configure Minimization Settings",
                "heating": "Configure Heating/Thermal Equilibration", 
                "equilibration": "Configure Pressure Equilibration",
                "production": "Configure Production MD Settings",
                "assemble_workflow": "Assemble Complete Workflow",
                "load_workflow": "Load Workflow Preset",
                "save_workflow": "Save Current Workflow as Preset",
                "manage_templates": "Manage Custom Templates"
            }

    def get_workspace_requirements(self) -> List[str]:
        """No specific workspace requirements - can work standalone."""
        return []

    def can_process(self, workspace) -> bool:
        """Can always process."""
        return True

    def handle_menu_option(self, option: str) -> bool:
        """Handle menu option selection using the unified controller."""
        try:
            if not self.controller:
                self.console.print("[red]Controller not initialized. Please set processor first.[/red]")
                return False
                
            # Show redesign announcement on first use
            self._show_redesign_announcement()
            
            return self.controller.handle_menu_option(option)
            
        except Exception as e:
            self.console.print(f"[red]Error executing option '{option}': {e}[/red]")
            self.console.print("[yellow]If you encounter issues, please report at https://github.com/anthropics/claude-code/issues[/yellow]")
            return False
            
    def _show_redesign_announcement(self):
        """Show announcement about the redesigned system (once per session)."""
        if not hasattr(self, '_announcement_shown'):
            from rich.panel import Panel
            
            announcement = """[bold cyan]AMBER Input Generator 2.0[/bold cyan]

[green]✨ Redesigned System Features:[/green]
• [bold]Template-centric workflow[/bold] - Always start with rich, annotated templates
• [bold]Flexible configuration[/bold] - Choose wizard guidance OR direct editing  
• [bold]Phase-based organization[/bold] - Configure individual simulation phases
• [bold]Complete workflow assembly[/bold] - Build full protocols from components
• [bold]Custom template management[/bold] - Save and reuse your configurations

[cyan]Each simulation phase shows an annotated template with parameter explanations,
then you choose how to configure it - guided wizard or direct editing.[/cyan]"""

            self.console.print(Panel(announcement, 
                                   title="🚀 New Architecture", 
                                   border_style="cyan",
                                   padding=(1, 2)))
            self.console.print()
            
            self._announcement_shown = True
            
    # Backward compatibility methods (if needed)
    def process(self, workspace) -> Dict[str, Any]:
        """Process method for module registry compatibility."""
        return workspace
        
    def cleanup(self):
        """Clean up resources."""
        if self.controller:
            # Controller cleanup if needed
            pass
            
    def initialize(self):
        """Initialize the module."""
        # Module initialization if needed
        pass