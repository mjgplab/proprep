"""
Metal Module Worker

Core functionality for metal site detection, transformation, bond generation, and cluster extraction.
This worker provides the backend functionality for the metal site processing pipeline.
"""

import logging
from typing import Any, Dict, List, Optional

from proprep.metallo_prep.detection.detection_manager import MetalSiteDetection
from proprep.metallo_prep.transformation.transformer_manager import TransformerManager
from proprep.metallo_prep.bond_generation.metal_bond_generator import MetalBondGenerator
from proprep.metallo_prep.extraction.metal_cluster_generator import MetalClusterGenerator

logger = logging.getLogger(__name__)


class MetalWorker:
    """Handles metal site processing operations."""

    def __init__(self, processor, metal_module=None):
        """
        Initialize Metal Worker.

        Args:
            processor: The main processor instance
            metal_module: The parent MetalSiteModule instance
        """
        self.processor = processor
        self.metal_module = metal_module
        self.detector = None
        self.transformer = None
        self.bond_generator = None
        self.cluster_generator = None
        self.microstate_generator = None
        self.metal_sites = []
        self.excluded_residues = set()

        # Set logging level for all metal detection modules
        import logging
        logging.getLogger('proprep.metallo_prep.detection').setLevel(logging.WARNING)

    def _initialize_detector(self):
        """Initialize the metal site detector on-demand."""
        if self.detector is None:
            # Require metal_module to be provided
            if not self.metal_module:
                logger.error("Metal module is required for detector initialization")
                return False

            self.detector = MetalSiteDetection(self.metal_module)
            logger.debug("Metal site detector initialized")
        return True

    def _initialize_transformer(self):
        """Initialize the transformer manager on-demand."""
        if self.transformer is None:
            self.transformer = TransformerManager(
                self.processor, verbose=True, debug=True
            )
            logger.debug("Transformer manager initialized")
        return True

    def _initialize_bond_generator(self):
        """Initialize the bond generator on-demand."""
        if self.bond_generator is None:
            self.bond_generator = MetalBondGenerator(self.processor)
            logger.debug("Metal bond generator initialized")
        return True

    def _initialize_cluster_generator(self):
        """Initialize the cluster generator on-demand."""
        if self.cluster_generator is None:
            self.cluster_generator = MetalClusterGenerator(self.processor)
            logger.debug("Metal cluster generator initialized")
        return True
    
    def _initialize_microstate_generator(self):
        """Initialize the microstate generator on-demand."""
        if self.microstate_generator is None:
            from .microstate_generation.redox_microstate_generator import RedoxMicrostateGenerator
            self.microstate_generator = RedoxMicrostateGenerator()
            self.microstate_generator.set_processor(self.processor)
            logger.debug("Redox microstate generator initialized")
        return True

    # Detection operations
    def detect_metal_sites(
        self, workspace: Dict[str, Any], interactive: bool = True
    ) -> Dict[str, Any]:
        """
        Detect metal sites in the structure.

        Args:
            workspace: Current workspace
            interactive: Whether to run interactively

        Returns:
            Updated workspace
        """
        if not self._initialize_detector():
            return workspace

        try:
            # Use the detector's process method or interactive methods
            if interactive and hasattr(self.detector, "handle_menu_option"):
                result = self.detector.handle_menu_option("detect")
                if result and hasattr(self.detector, "metal_sites"):
                    self.metal_sites = self.detector.metal_sites
                    self.excluded_residues = getattr(
                        self.detector, "excluded_residues", set()
                    )
                    workspace.set("metal_sites", self.metal_sites)
                    workspace.set("excluded_residues", self.excluded_residues)
            else:
                # Non-interactive processing
                workspace = self.detector.process(workspace)
                if hasattr(self.detector, "metal_sites"):
                    self.metal_sites = self.detector.metal_sites
                if hasattr(self.detector, "excluded_residues"):
                    self.excluded_residues = self.detector.excluded_residues

        except Exception as e:
            logger.error(f"Error during metal site detection: {str(e)}")
            raise

        return workspace

    def import_metal_sites(
        self, workspace: Dict[str, Any], file_path: Optional[str] = None
    ) -> Dict[str, Any]:
        """Import metal site information from file."""
        if not self._initialize_detector():
            return workspace

        try:
            if hasattr(self.detector, "handle_menu_option"):
                result = self.detector.handle_menu_option("import")
                if result:
                    workspace = self._sync_detector_state(workspace)
        except Exception as e:
            logger.error(f"Error importing metal sites: {str(e)}")
            logger.exception(e)
            raise

        return workspace

    def view_metal_sites(self, workspace: Dict[str, Any]) -> bool:
        """View metal site analysis results."""
        if not self._initialize_detector():
            return False

        try:
            if hasattr(self.detector, "handle_menu_option"):
                return self.detector.handle_menu_option("view")
        except Exception as e:
            logger.error(f"Error viewing metal sites: {str(e)}")
            return False

        return True

    def show_metal_sites_summary(self, workspace: Dict[str, Any]) -> bool:
        """Show summary of detected metal sites."""
        if not self._initialize_detector():
            return False

        try:
            if hasattr(self.detector, "handle_menu_option"):
                return self.detector.handle_menu_option("summary")
        except Exception as e:
            logger.error(f"Error showing metal sites summary: {str(e)}")
            return False

        return True

    def export_metal_sites(self, workspace: Dict[str, Any]) -> bool:
        """Export metal site information."""
        if not self._initialize_detector():
            return False

        try:
            if hasattr(self.detector, "handle_menu_option"):
                return self.detector.handle_menu_option("export")
        except Exception as e:
            logger.error(f"Error exporting metal sites: {str(e)}")
            return False

        return True

    def save_structure(self, workspace: Dict[str, Any]) -> bool:
        """Save the filtered structure with coordinating residues."""
        if not self._initialize_detector():
            return False

        try:
            if hasattr(self.detector, "handle_menu_option"):
                return self.detector.handle_menu_option("save_structure")
        except Exception as e:
            logger.error(f"Error saving structure: {str(e)}")
            return False

        return True

    # Transformation operations

    def transform_metal_sites(self, workspace: Dict[str, Any]) -> Dict[str, Any]:
        """Transform metal sites for MD simulation."""
        
        # CRITICAL: Capture the current structure BEFORE any transformations
        if not workspace.get("untransformed_structure"):
            # Get the current structure using StructureSelector (prioritize repaired over filtered)
            from proprep.utils.structure_selector import StructureSelector
            selector = StructureSelector(workspace, self.processor.console, processor=self.processor)
            current_structure = selector.get_structure_by_key(
                "repaired_structure",
                fallback_keys=["filtered_structure", "structure"],
                require_exists=False
            )
            
            if current_structure:
                # Store this as the untransformed structure
                workspace.set("untransformed_structure", current_structure)
                logger.debug("Stored untransformed structure for mapping")
        
        if not self._initialize_transformer():
            return workspace

        try:
            if hasattr(self.transformer, "handle_menu_option"):
                result = self.transformer.handle_menu_option("transform")
                if result:
                    workspace = self._sync_transformer_state(workspace)
        except Exception as e:
            logger.error(f"Error transforming metal sites: {str(e)}")
            raise

        return workspace

    def save_transformed_structure(self, workspace: Dict[str, Any]) -> bool:
        """Save the transformed structure."""
        if not self._initialize_transformer():
            return False

        try:
            if hasattr(self.transformer, "handle_menu_option"):
                return self.transformer.handle_menu_option("save_transformed")
        except Exception as e:
            logger.error(f"Error saving transformed structure: {str(e)}")
            return False

        return True

    #Redox microstate generation
    def show_microstate_menu(self, workspace: Dict[str, Any]) -> bool:
        """Generate redox microstates for metal sites."""
        if not self._initialize_microstate_generator():
            return False
        
        # The microstate generator handles its own menu system
        # (just like detection_manager does)
        return self.microstate_generator.run_microstate_menu()

    def setup_transformer_assignments(self, workspace: Dict[str, Any]) -> bool:
        """Setup transformer assignments for redox sites."""
        if not self._initialize_microstate_generator():
            return False
        
        return self.microstate_generator.setup_transformer_assignments()

    def view_available_microstates(self, workspace: Dict[str, Any]) -> bool:
        """View available redox microstates."""
        if not self._initialize_microstate_generator():
            return False
        
        return self.microstate_generator.view_available_microstates()

    def select_combinations_interactive(self, workspace: Dict[str, Any]) -> bool:
        """Select redox/spin state combinations for each site."""
        if not self._initialize_microstate_generator():
            return False
        
        return self.microstate_generator.select_combinations_interactive()

    def select_microstates_interactive(self, workspace: Dict[str, Any]) -> bool:
        """Select microstates for generation."""
        if not self._initialize_microstate_generator():
            return False
        
        return self.microstate_generator.select_microstates_interactive()

    def generate_microstates(self, workspace: Dict[str, Any], generate_all: bool = False) -> bool:
        """Generate microstate PDB files."""
        if not self._initialize_microstate_generator():
            return False
        
        return self.microstate_generator.generate_microstates(generate_all)

    # Bond generation operations
    def generate_tleap_bonds(self, workspace: Dict[str, Any]) -> bool:
        """Generate TLEaP bond commands for metal sites."""
        if not self._initialize_bond_generator():
            return False

        try:
            if hasattr(self.bond_generator, "handle_menu_option"):
                return self.bond_generator.handle_menu_option("generate_tleap")
        except Exception as e:
            logger.error(f"Error generating TLEaP bonds: {str(e)}")
            return False

        return True

    # Cluster extraction operations
    def extract_clusters(self, workspace: Dict[str, Any]) -> bool:
        """Extract active site cluster models."""
        if not self._initialize_cluster_generator():
            return False

        try:
            if hasattr(self.cluster_generator, "handle_menu_option"):
                return self.cluster_generator.handle_menu_option("extract_clusters")
        except Exception as e:
            logger.error(f"Error extracting clusters: {str(e)}")
            return False

        return True

    # Utility methods
    def _sync_detector_state(self, workspace: Dict[str, Any]) -> Dict[str, Any]:
        """Synchronize detector state with workspace."""
        if self.detector:
            if hasattr(self.detector, "metal_sites"):
                self.metal_sites = self.detector.metal_sites
                workspace.set("metal_sites", self.metal_sites)
            if hasattr(self.detector, "excluded_residues"):
                self.excluded_residues = self.detector.excluded_residues
                workspace.set("excluded_residues", self.excluded_residues)
        return workspace

    def _sync_transformer_state(self, workspace: Dict[str, Any]) -> Dict[str, Any]:
        """Synchronize transformer state with workspace."""
        if self.transformer:
            # Update workspace with any transformed structures or data
            # This depends on what the transformer produces
            pass
        return workspace

    def get_metal_sites(self) -> List[Any]:
        """Get the current metal sites."""
        if self.detector and hasattr(self.detector, "metal_sites"):
            return self.detector.metal_sites
        return self.metal_sites

    def get_excluded_residues(self) -> set:
        """Get residues to exclude from protonation."""
        if self.detector and hasattr(self.detector, "excluded_residues"):
            return self.detector.excluded_residues
        return self.excluded_residues

    def generate_restraint_masks(self, workspace: Dict[str, Any]) -> bool:
        """Generate restraint masks for MD minimization using transformed metal sites."""
        try:
            import os
            from pathlib import Path
            # Check prerequisites
            metal_sites = workspace.get("metal_sites", [])
            if not metal_sites:
                self.processor.console.print("[yellow]No metal sites found in workspace.[/yellow]")
                self.processor.console.print("[yellow]Please run 'Detect redox sites' first.[/yellow]")
                return False

            # Check for transformed structure using StructureSelector (prioritize file paths)
            from proprep.utils.structure_selector import StructureSelector
            selector = StructureSelector(workspace, self.processor.console, processor=self.processor)

            structure_source = None
            pdb_content = None

            # Try to get structure file path with priority
            result = selector.get_structure(
                priority_override=[
                    "repaired_pdb_file",
                    "filtered_pdb_file",
                    "hstripped_pdb_file",
                    "transformed_pdb_file",
                ],
                return_key=True,
                silent=True,
            )

            if result:
                structure_file, structure_key = result
                if os.path.exists(str(structure_file)):
                    structure_source = structure_key.replace("_pdb_file", "_file")
                    with open(str(structure_file), 'r') as f:
                        pdb_content = f.read()
                    self.processor.console.print(f"[green]Using transformed structure from: {structure_file}[/green]")

            if not pdb_content:
                # Check for PDB content directly
                pdb_content = workspace.get("filtered_structure_pdb_content")
                if pdb_content:
                    structure_source = "pdb_content"
                    self.processor.console.print("[green]Using transformed structure from workspace PDB content[/green]")
                
            if not pdb_content:
                self.processor.console.print("[yellow]No transformed structure found.[/yellow]")
                self.processor.console.print("[yellow]Please run 'Generate single redox microstate' first to transform the structure.[/yellow]")
                return False

            # Initialize restraint mask generator
            from .restraint_mask_generator import RestraintMaskGenerator
            generator = RestraintMaskGenerator(self.processor.console)

            # Generate restraint masks
            result = generator.generate_masks(metal_sites, pdb_content)
            
            if result:
                # Save the restraint mask to workspace
                workspace["redox_restraint_mask"] = result["restraint_mask"]
                workspace["redox_restraint_info"] = result["info"]
                
                self.processor.console.print(f"\n[green]✓ Restraint mask generated successfully[/green]")
                self.processor.console.print(f"[cyan]Restraint mask:[/cyan] {result['restraint_mask']}")
                
                # For residue-based restraints, don't show atom count as it's unreliable when residues are excluded
                if result['info']['restraint_by_residues']:
                    self.processor.console.print(f"[cyan]Method:[/cyan] Residue-based restraints")
                else:
                    self.processor.console.print(f"[cyan]Atoms included:[/cyan] {len(result['info']['included_atoms'])} redox site atoms")
                    
                self.processor.console.print(f"[cyan]Restraint mask saved to workspace as 'redox_restraint_mask'[/cyan]")
                
                return True
            else:
                self.processor.console.print("[red]Failed to generate restraint mask[/red]")
                return False
                
        except Exception as e:
            logger.error(f"Error generating restraint masks: {str(e)}")
            self.processor.console.print(f"[red]Error generating restraint masks: {str(e)}[/red]")
            return False

    def can_process(self, workspace: Dict[str, Any]) -> bool:
        """Check if the worker can process the current workspace."""
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(workspace, self.processor.console)
        status = selector.get_structure_status()
        return status.get("has_any", False)

    def process(self, workspace: Dict[str, Any]) -> Dict[str, Any]:
        """Process the workspace using metal site detection."""
        if self.can_process(workspace):
            # Check if we should auto-detect metal sites
            auto_detect = workspace.get("auto_detect_metals", False)
            if auto_detect:
                workspace = self.detect_metal_sites(workspace, interactive=False)

        return workspace

    def cleanup(self):
        """Clean up worker resources."""
        self.detector = None
        self.transformer = None
        self.bond_generator = None
        self.cluster_generator = None
        self.metal_sites = []
        self.excluded_residues = set()
