#!/usr/bin/env python3
"""
RedoxSite Preparation Integration Layer

Modern integration layer for the RedoxSite preparation system.
Provides a clean interface between the main application and redox transformation functionality.

Author: Claude Code Implementation
"""

import logging
from collections import OrderedDict
from typing import Any, Dict, List, Optional
from pathlib import Path

# Setup logging
logger = logging.getLogger(__name__)

# Import required components for registration
from proprep.utils.module_registry import ProcessingModule, register_module
from proprep.utils.prompts import prompt_with_context

# Import new RedoxSite transformation system
from .transformation.redox_transformation_manager import RedoxTransformationManager

# Import commands
from .redox_commands import (
    GenerateSingleRedoxMicrostateCommand,
    GenerateAllRedoxMicrostatesCommand,
    CreateCustomTransformerCommand,
)


@register_module
class RedoxSitePreparationModule(ProcessingModule):
    """
    Modern RedoxSite Preparation module
    
    Handles redox site transformation and microstate generation using
    the new RedoxSite system that replaced the legacy metallo_prep.
    """
    
    NAME = "Redox Site Preparer"
    DESCRIPTION = "Apply FF-compatible names and atom partitioning to redox sites"
    VERSION = "2.0.0"
    CATEGORY = "structure_prep"
    PRIORITY = 60  # After PDB Filter (50) but before tLEaP (70)
    REQUIRES = ["detected_redox_sites"]  # Requires redox sites from PDB Filter
    
    def __init__(self):
        super().__init__()
        
        # Initialize transformation manager
        self.transformation_manager = None
        self.processor = None  # Will be set when methods are called
        
        # Define available commands (reordered: custom transformer first)
        self.commands = OrderedDict([
            ("create_custom_transformer", {
                "description": "Create transformer (interactive PDB editor)",
                "category": "advanced",
                "enabled": True
            }),
            ("generate_single_microstate", {
                "description": "Generate single redox microstate",
                "category": "transformation"
            }),
            ("generate_all_microstates", {
                "description": "Generate all redox microstates (batch processing)",
                "category": "transformation"
            }),
        ])
    
    def get_menu_options(self) -> Dict[str, str]:
        """Get module menu options"""
        menu = OrderedDict()

        for key, info in self.commands.items():
            if info.get("enabled", True):  # Skip disabled commands
                menu[key] = info["description"]

        return menu

    def get_enhanced_menu_options(self, workspace):
        """
        Get menu options with enhanced status information.

        Args:
            workspace: Current workspace

        Returns:
            List of MenuOption objects with status
        """
        from proprep.utils.enhanced_menu import MenuOption, OptionStatus

        options = []

        # Check workspace state
        has_transformed = workspace.get("transformed_pdb_file") is not None
        has_batch = workspace.get("generated_microstate_pdbs") is not None
        has_generated = has_transformed or has_batch
        has_sites = bool(workspace.get("detected_redox_sites"))

        # Option 1: Table-driven transformer creator (the recommended path).
        # Gated on a loaded structure with detected redox sites -- the editor
        # shows the atoms of a detected site, so it needs one.
        if has_generated:
            status = OptionStatus.BLOCKED
            dep_text = "[Already generated microstate(s)] ○"
        elif not has_sites:
            status = OptionStatus.BLOCKED
            dep_text = "[Load a structure & detect redox sites first] ○"
        else:
            status = OptionStatus.AVAILABLE
            dep_text = ""

        options.append(MenuOption(
            key="1",
            description="Create transformer (interactive PDB editor)",
            status=status,
            dependency_text=dep_text
        ))

        # Option 2: Generate single microstate - blocked if already generated
        if has_generated:
            status = OptionStatus.BLOCKED
            dep_text = "[Already generated microstate(s)] ○"
        else:
            status = OptionStatus.AVAILABLE
            dep_text = ""

        options.append(MenuOption(
            key="2",
            description="Generate single redox microstate",
            status=status,
            dependency_text=dep_text
        ))

        # Option 3: Generate all microstates - blocked if already generated
        if has_generated:
            status = OptionStatus.BLOCKED
            dep_text = "[Already generated microstate(s)] ○"
        else:
            status = OptionStatus.AVAILABLE
            dep_text = ""

        options.append(MenuOption(
            key="3",
            description="Generate all redox microstates (batch processing)",
            status=status,
            dependency_text=dep_text
        ))

        return options

    def get_menu_suggestion(self, workspace):
        """
        Get a suggestion for the next recommended action.

        Args:
            workspace: Current workspace

        Returns:
            Suggestion text or None
        """
        has_transformed = workspace.get("transformed_pdb_file") is not None
        has_batch = workspace.get("generated_microstate_pdbs") is not None
        redox_sites = workspace.get("detected_redox_sites", [])
        num_sites = len(redox_sites)

        if has_batch:
            batch_files = workspace.get("generated_microstate_pdbs", [])
            return f"✓ Generated {len(batch_files)} microstate PDB files. Press [m] to return to the main menu"
        elif has_transformed:
            return "✓ Generated single microstate PDB. Press [m] to return to the main menu"
        else:
            if num_sites == 0:
                return "No redox sites detected. Create custom transformer first (option 1) if needed, or return to the Redox Site Detector to detect sites"
            elif num_sites == 1:
                return f"Found {num_sites} redox site. Create custom transformer first (option 1) if needed, then generate microstate (option 2)"
            else:
                return f"Found {num_sites} redox sites. Create custom transformer first (option 1) if needed, then generate single (option 2) or batch (option 3)"

    def get_workspace_requirements(self) -> List[str]:
        """Get workspace requirements for this module"""
        return self.REQUIRES

    def get_workspace_outputs(self) -> List[str]:
        """Get workspace outputs"""
        return [
            "untransformed_structure",
            "transformed_pdb_file",
            "transformed_structure",
            "generated_microstate_pdbs",
            "metal_sites",
            "excluded_residues",
        ]

    def handle_menu_option(self, option: str) -> bool:
        """Handle a menu option selection"""
        try:
            if option == "create_custom_transformer":
                return self.create_custom_transformer()
            elif option == "generate_single_microstate":
                return self.generate_single_microstate()
            elif option == "generate_all_microstates":
                return self.generate_all_microstates()
        except Exception as e:
            logger.error(f"Error executing menu option '{option}': {e}")
            if self.processor:
                self.processor.console.print(f"[red]Error: {e}[/red]")
        return False
    
    # Action implementations
    
    def generate_single_microstate(self) -> bool:
        """Generate a single redox microstate using RedoxSite transformation system."""
        try:
            # Get workspace from processor
            workspace = self.processor._get_workspace()
            
            # Check for RedoxSite data from comprehensive detector
            redox_sites = workspace.get("detected_redox_sites", [])
            if not redox_sites:
                self.processor.console.print("[red]No redox sites detected. Please run PDB Filter first.[/red]")
                return False
            
            # Initialize transformation manager
            if not self.transformation_manager:
                self.transformation_manager = RedoxTransformationManager(self.processor)
            
            # Run single microstate generation
            return self.transformation_manager.transform_redox_sites()
            
        except Exception as e:
            logger.error(f"Failed to generate single redox microstate: {e}")
            self.processor.console.print(f"[red]Error: {e}[/red]")
            return False
    
    def generate_all_microstates(self) -> bool:
        """Generate all possible redox microstates in batch mode."""
        try:
            # Get workspace from processor
            workspace = self.processor._get_workspace()
            
            # Check for RedoxSite data
            redox_sites = workspace.get("detected_redox_sites", [])
            if not redox_sites:
                self.processor.console.print("[red]No redox sites detected. Please run PDB Filter first.[/red]")
                return False
            
            # Initialize transformation manager  
            if not self.transformation_manager:
                self.transformation_manager = RedoxTransformationManager(self.processor)
            
            # Run batch microstate generation
            return self.transformation_manager.batch_transformation()
            
        except Exception as e:
            logger.error(f"Failed to generate all redox microstates: {e}")
            self.processor.console.print(f"[red]Error: {e}[/red]")
            return False
    
    def create_custom_transformer(self) -> bool:
        """Create a transformer with the interactive, table-driven PDB editor.

        Shows the atoms of a detected redox site grouped by residue and applies
        one edit at a time, saving a reusable JSON transformer spec. This is the
        transformer-creation path exposed by the menu; it replaces the earlier
        guided v3 creator in place.
        """
        try:
            from .transformation.table_transformer_creator import TableTransformerCreator

            workspace = self.processor._get_workspace()
            redox_sites = workspace.get("detected_redox_sites", []) or []
            if not redox_sites:
                self.processor.console.print(
                    "[red]No redox sites detected.[/red] Load a structure and run the "
                    "Redox Site Detector (or PDB Filter) first, then create a transformer."
                )
                return False

            # Pick a site whose atoms the editor will show.
            redox_site = redox_sites[0]
            if len(redox_sites) > 1:
                self.processor.console.print(
                    f"\n[bold]{len(redox_sites)} detected redox site(s):[/bold]"
                )
                site_map = {}
                for i, site in enumerate(redox_sites, 1):
                    site_id = getattr(site, "site_id", f"site_{i}")
                    n_atoms = len(getattr(site, "atoms", []))
                    self.processor.console.print(f"  [cyan]{i}[/cyan]. {site_id} ({n_atoms} atoms)")
                    site_map[str(i)] = site_id
                choice = prompt_with_context(
                    self.processor, "Select site to build a transformer for",
                    module="Table Transformer Creator",
                    description="Select detected redox site", options_map=site_map,
                ).strip()
                try:
                    idx = int(choice) - 1
                    if not (0 <= idx < len(redox_sites)):
                        raise ValueError
                    redox_site = redox_sites[idx]
                except ValueError:
                    self.processor.console.print("[red]Invalid selection.[/red]")
                    return False

            TableTransformerCreator(self.processor, redox_site).create()
            return True

        except Exception as e:
            logger.error(f"Failed to create transformer: {e}", exc_info=True)
            self.processor.console.print(f"[red]Error creating transformer: {e}[/red]")
            return False

    def get_status_info(self, workspace) -> Dict[str, Any]:
        """Get status information for the workspace display."""
        redox_sites = workspace.get("detected_redox_sites", [])
        transformed_pdb = workspace.get("transformed_pdb_file")
        
        status = {
            "redox_sites_detected": len(redox_sites),
            "transformed_structure": bool(transformed_pdb),
        }
        
        if transformed_pdb:
            status["transformed_file"] = Path(transformed_pdb).name
        
        return status