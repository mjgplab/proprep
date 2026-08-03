"""
Metal Module Integration Layer for MPSA

This module provides a thin integration layer between the MPSA main program
and the specialized metal-related modules for detecting and processing
metal binding sites in protein structures.
"""

import logging
from collections import OrderedDict
from typing import Any, Dict, List, Optional
from pathlib import Path

from proprep.utils.prompts import prompt_with_context, confirm_with_context

# Setup logging
logger = logging.getLogger(__name__)

# Import required components for registration
from proprep.utils.module_registry import ProcessingModule, register_module

# Import worker and commands
from .metallo_worker import MetalWorker
from .metallo_commands import (
    DetectMetalSitesCommand,
    ImportMetalSitesCommand,
    ViewResultsCommand,
    ExportMetalSitesCommand,
    CreateCustomTransformerCommand,
    GenerateSingleMicrostateCommand,
    GenerateMicrostatesCommand, 
    GenerateTleapBondsCommand,
    ExtractClustersCommand,
    GenerateRestraintMasksCommand,
)

from .microstate_generation.redox_microstate_generator import RedoxMicrostateGenerator

# Import new RedoxSite transformation system
try:
    from proprep.redoxsite_prep.transformation.redox_transformation_manager import RedoxTransformationManager
    REDOX_TRANSFORMATION_AVAILABLE = True
except ImportError:
    REDOX_TRANSFORMATION_AVAILABLE = False
    logger.warning("RedoxSite transformation system not available")

# @register_module  # DISABLED: Replaced by redoxsite_prep.RedoxSitePreparationModule
class MetalSiteModule(ProcessingModule):
    """
    Metal Site module for MPSA

    This module serves as an integration layer between MPSA and the specialized
    metal-related modules. It coordinates metal site detection, transformation,
    bond generation, and cluster model extraction using the command pattern.
    """

    NAME = "Redox Site Preparation (Legacy - DISABLED)"
    DESCRIPTION = (
        "Detects, defines, and prepares redox-active sites (metal and non-metallic) for MD simulations "
    )
    VERSION = "1.0.0"
    CATEGORY = "analysis"
    REQUIRES = ["PDB Loader", "PDB Filter"]
    PRIORITY = 5  # After structure completeness analysis

    def initialize(self):
        """Initialize module resources"""
        self.worker = MetalWorker(self.processor, self)
        self._record_local_command = self._record_command_history

    def _record_command_history(self, command):
        """Record command in local history if available"""
        if hasattr(self, "_command_history") and self._command_history:
            self._command_history.append(command)

    def get_workspace_requirements(self) -> List[str]:
        """Get workspace requirements"""
        return ["structure"]  # Require structure from PDB Loader

    def get_workspace_outputs(self) -> List[str]:
        """Get workspace outputs"""
        return ["metal_sites", "excluded_residues", "cross_chain_coordination"]

    def can_process(self, workspace) -> bool:
        """Check if the module can process the current workspace"""
        return self.worker.can_process(workspace)

    def get_menu_options(self) -> Dict[str, str]:
        """Get module menu options"""
        # Use an ordered dictionary to maintain option order
        menu = OrderedDict()

        menu["import"] = "Import redox site information from file"
        menu["detect"] = "Detect redox sites (metal and non-metallic)"
        menu["extract_clusters"] = "Extract active site cluster models"
        menu["view_results"] = "View redox site analysis results"
        menu["export"] = "Export redox site information"
        menu["create_transformer"] = "Create custom transformer for new redox site type"
        menu["generate_single"] = "Generate single redox microstate"
        menu["generate_microstates"] = "Generate all redox microstates" 
        menu["generate_tleap"] = "Generate TLEaP bond commands for redox sites"
        menu["generate_restraints"] = "Generate restraint masks for MD minimization"

        return menu

    def handle_menu_option(self, option: str) -> bool:
        """Handle menu option selections using command pattern"""
        try:
            if option == "detect":
                command = DetectMetalSitesCommand(self.processor, interactive=True)
                return command.execute_with_error_handling()

            elif option == "import":
                command = ImportMetalSitesCommand(self.processor)
                return command.execute_with_error_handling()

            elif option == "view_results":
                command = ViewResultsCommand(self.processor)
                return command.execute_with_error_handling()

            elif option == "export":
                command = ExportMetalSitesCommand(self.processor)
                return command.execute_with_error_handling()

            elif option == "create_transformer":
                command = CreateCustomTransformerCommand(self.processor)
                return command.execute_with_error_handling()

            elif option == "generate_single":
                command = GenerateSingleMicrostateCommand(self.processor)
                return command.execute_with_error_handling()

            elif option == "generate_microstates":
                command = GenerateMicrostatesCommand(self.processor)
                return command.execute_with_error_handling()

            elif option == "generate_tleap":
                command = GenerateTleapBondsCommand(self.processor)
                return command.execute_with_error_handling()

            elif option == "extract_clusters":
                command = ExtractClustersCommand(self.processor)
                return command.execute_with_error_handling()

            elif option == "generate_restraints":
                command = GenerateRestraintMasksCommand(self.processor)
                return command.execute_with_error_handling()

            else:
                logger.warning(f"Unknown menu option: {option}")
                return False

        except Exception as e:
            logger.error(f"Error handling menu option '{option}': {str(e)}")
            logger.exception(e)
            if hasattr(self.processor, "console"):
                self.processor.console.print(
                    f"[red]Error executing {option}: {str(e)}[/red]"
                )
            return False
        
    def process(self, workspace: Dict[str, Any]) -> Dict[str, Any]:
        """Process the workspace using metal site detection"""
        return self.worker.process(workspace)

    def get_workspace_display(
        self, workspace_key: str, value: Any, console: Any
    ) -> bool:
        """
        Custom display method for metal site data in workspace
        """
        # Delegate to worker if it has display method
        if hasattr(self.worker, "get_workspace_display"):
            return self.worker.get_workspace_display(workspace_key, value, console)
        return False

    # Command action methods - these are called by ModuleActionCommand

    def detect_metal_sites(self, interactive: bool = True) -> Dict[str, Any]:
        """Command action: Detect metal sites in structure."""
        workspace = self.processor._get_workspace()
        
        # Run the detection
        result = self.worker.detect_metal_sites(workspace, interactive)
        
        # Check if cross-coordinating residues were identified and included
        metal_sites = workspace.get("metal_sites", [])
        has_coordinating_residues = False
        
        for site in metal_sites:
            if isinstance(site, dict) and site.get("coordinating_residues"):
                has_coordinating_residues = True
                break
        
        # If cross-coordinating residues were found, offer to save filtered structure
        if has_coordinating_residues and interactive:
            if confirm_with_context(
                processor=self.processor,
                prompt="Save filtered structure with coordinating residues?",
                default=True,
                module="Metal Site Detection",
                description="Save filtered structure with coordinating residues"
            ):
                save_success = self.worker.save_structure(workspace)
                if save_success:
                    self.processor.console.print("[green]✓ Filtered structure with coordinating residues saved[/green]")
        
        return result

    def import_metal_sites(self, file_path: Optional[str] = None) -> Dict[str, Any]:
        """Command action: Import metal site information from file."""
        workspace = self.processor._get_workspace()
        return self.worker.import_metal_sites(workspace, file_path)

    def view_combined_results(self) -> bool:
        """Command action: View combined results (summary + detailed if requested)."""
        workspace = self.processor._get_workspace()
        
        # Show summary first
        summary_success = self.worker.show_metal_sites_summary(workspace)
        if not summary_success:
            return False
        
        # Ask if user wants detailed view
        if confirm_with_context(
            processor=self.processor,
            prompt="Show detailed reports?",
            module="Metal Site Detection",
            description="Show detailed reports"
        ):
            return self.worker.view_metal_sites(workspace)
        
        return True

    def export_metal_sites(self) -> bool:
        """Command action: Export detection results."""
        workspace = self.processor._get_workspace()
        return self.worker.export_metal_sites(workspace)

    def save_structure(self) -> bool:
        """Command action: Save filtered structure with coordinating residues."""
        workspace = self.processor._get_workspace()
        return self.worker.save_structure(workspace)

    def transform_metal_sites(self) -> Dict[str, Any]:
        """Command action: Transform metal sites for MD simulation."""
        workspace = self.processor._get_workspace()
        return self.worker.transform_metal_sites(workspace)

    def save_transformed_structure(self) -> bool:
        """Command action: Save the transformed structure."""
        workspace = self.processor._get_workspace()
        return self.worker.save_transformed_structure(workspace)

    def create_custom_transformer(self) -> bool:
        """Launch the transformer creator."""
        try:
            from proprep.metallo_prep.transformation.transformer_creator import TransformerCreator
            
            # Get the path to the transformers directory in the MPSA installation
            import proprep.metallo_prep.transformation.transformers as transformers_module
            transformers_dir = Path(transformers_module.__file__).parent
            
            # Create transformer creator instance pointing to the actual transformers directory
            creator = TransformerCreator(output_dir=transformers_dir)
            creator.console = self.processor.console
            
            # Run the creation process
            creator.create_transformer()
            
            return True
            
        except Exception as e:
            self.processor.console.print(f"[red]Error creating transformer: {e}[/red]")
            return False
        
    def generate_single_microstate(self) -> bool:
        """Command action: Generate a single redox microstate using new RedoxSite transformation system."""
        workspace = self.processor._get_workspace()
        
        # Check for RedoxSite data from comprehensive detector
        redox_sites = workspace.get("detected_redox_sites", [])
        if redox_sites:
            # Use new RedoxSite transformation system
            return self._transform_redox_sites_with_new_system(workspace)
        else:
            # Fallback to legacy MetalSite transformation for backward compatibility
            self.processor.console.print("[yellow]No RedoxSites detected, using legacy MetalSite transformation[/yellow]")
            transform_success = self.worker.transform_metal_sites(workspace)
            if not transform_success:
                return False
            return self.worker.save_transformed_structure(workspace)

    def _transform_redox_sites_with_new_system(self, workspace) -> bool:
        """Transform RedoxSites using the new redoxsite_prep transformation system."""
        if not REDOX_TRANSFORMATION_AVAILABLE:
            self.processor.console.print("[red]RedoxSite transformation system not available[/red]")
            return False
        
        try:
            # Initialize the new transformation manager
            transformation_manager = RedoxTransformationManager(
                processor=self.processor,
                console=self.processor.console
            )
            
            # Run the transformation process
            self.processor.console.print("\n[bold blue]🔄 Starting RedoxSite Transformation Process[/bold blue]")
            success = transformation_manager.transform_redox_sites()
            
            if success:
                self.processor.console.print("[bold green]✅ RedoxSite transformation completed successfully[/bold green]")
                
                # Update workspace with any results from the transformation
                if hasattr(transformation_manager, 'redox_sites'):
                    workspace["detected_redox_sites"] = transformation_manager.redox_sites
                
                # The transformation manager handles structure saving internally
                return True
            else:
                self.processor.console.print("[red]❌ RedoxSite transformation failed[/red]")
                return False
                
        except Exception as e:
            logger.error(f"RedoxSite transformation failed: {e}")
            self.processor.console.print(f"[red]RedoxSite transformation error: {e}[/red]")
            return False

    def generate_microstates(self, generate_all: bool = False) -> bool:
        """Command action: Generate all redox microstates."""
        workspace = self.processor._get_workspace()
        return self.worker.generate_microstates(workspace, generate_all=generate_all)

    def setup_transformer_assignments(self) -> bool:
        """Command action: Setup transformer assignments for redox sites."""
        workspace = self.processor._get_workspace()
        return self.worker.setup_transformer_assignments(workspace)

    def view_available_microstates(self) -> bool:
        """Command action: View available redox microstates."""
        workspace = self.processor._get_workspace()
        return self.worker.view_available_microstates(workspace)

    def select_combinations_interactive(self) -> bool:
        """Command action: Select redox/spin state combinations for each site."""
        workspace = self.processor._get_workspace()
        return self.worker.select_combinations_interactive(workspace)

    def select_microstates_interactive(self) -> bool:
        """Command action: Select microstates for generation."""
        workspace = self.processor._get_workspace()
        return self.worker.select_microstates_interactive(workspace)

    def show_microstate_menu(self) -> bool:
        """Command action: Show redox microstate menu."""
        workspace = self.processor._get_workspace()
        return self.worker.show_microstate_menu(workspace) 


    def _update_transformer_registration(self, transformers_dir, transformer_name, class_name, category):
        """Update the register_transformers.py file to include the new transformer."""
        try:
            register_file = transformers_dir / 'register_transformers.py'
            
            if register_file.exists():
                # Read the current content
                content = register_file.read_text()
                
                # Check if already registered
                if transformer_name in content and class_name in content:
                    self.processor.console.print(f"[yellow]Transformer '{transformer_name}' already registered[/yellow]")
                    return
                
                # Find where to insert the new transformer registration
                lines = content.split('\n')
                
                # Look for the end of the existing transformer imports
                # We'll insert the new transformer after the last "except ImportError" block
                insert_index = -1
                for i in range(len(lines) - 1, -1, -1):
                    if "except ImportError as e:" in lines[i]:
                        # Find the end of this except block
                        j = i + 1
                        while j < len(lines) and lines[j].strip() and not lines[j].startswith('try:'):
                            j += 1
                        insert_index = j
                        break
                
                if insert_index == -1:
                    # If we can't find a good spot, insert before the final lines
                    for i, line in enumerate(lines):
                        if "available_transformers = MetalSiteTransformer.get_available_transformers()" in line:
                            insert_index = i
                            break
                
                if insert_index > 0:
                    # Create the new registration block
                    new_registration = [
                        "",
                        f"# Register {transformer_name}",
                        "try:",
                        f"    from proprep.metallo_prep.transformation.transformers.{category}.{transformer_name} import {class_name}",
                        f'    MetalSiteTransformer.register_transformer("{transformer_name}", {class_name})',
                        f'    logger.debug(f"Successfully registered {class_name} as \'{transformer_name}\'")',
                        "except ImportError as e:",
                        f'    logger.error(f"Failed to import {class_name}: {{e}}")',
                    ]
                    
                    # Insert the new registration
                    lines[insert_index:insert_index] = new_registration
                    content = '\n'.join(lines)
                    
                    # Write back
                    register_file.write_text(content)
                    self.processor.console.print(f"[green]Updated registration for '{transformer_name}'[/green]")
                else:
                    raise ValueError("Could not find appropriate location to insert registration")
                    
        except Exception as e:
            self.processor.console.print(f"[yellow]Could not auto-update registration: {e}[/yellow]")
            self.processor.console.print("[yellow]You may need to manually add the following to register_transformers.py:[/yellow]")
            self.processor.console.print(f"\ntry:")
            self.processor.console.print(f"    from proprep.metallo_prep.transformation.transformers.{category}.{transformer_name} import {class_name}")
            self.processor.console.print(f'    MetalSiteTransformer.register_transformer("{transformer_name}", {class_name})')
            self.processor.console.print(f"except ImportError as e:")
            self.processor.console.print(f'    logger.error(f"Failed to import {class_name}: {{e}}")')
            
    def generate_tleap_bonds(self) -> bool:
        """Command action: Generate TLEaP bond commands for metal sites."""
        workspace = self.processor._get_workspace()
        return self.worker.generate_tleap_bonds(workspace)

    def extract_clusters(self) -> bool:
        """Command action: Extract active site cluster models."""
        workspace = self.processor._get_workspace()
        return self.worker.extract_clusters(workspace)

    def generate_restraint_masks(self) -> bool:
        """Command action: Generate restraint masks for MD minimization."""
        workspace = self.processor._get_workspace()
        return self.worker.generate_restraint_masks(workspace)

    # Helper methods for backward compatibility

    def get_metal_sites(self):
        """Get the current metal sites"""
        return self.worker.get_metal_sites()

    def get_excluded_residues(self):
        """Get residues to exclude from protonation"""
        return self.worker.get_excluded_residues()

    def cleanup(self):
        """Clean up module resources"""
        if self.worker:
            self.worker.cleanup()

