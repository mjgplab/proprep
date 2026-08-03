"""
Unified QM/MM Preparator Module

Single entry point for Gaussian (ONIOM) and ORCA QM/MM input generation.
Handles shared topology loading and trajectory frame extraction, then
dispatches to the selected interface.
"""

import logging
import os
from typing import Any, Dict, List, Optional

from proprep.utils.module_registry import ProcessingModule, register_module
from proprep.utils.prompts import prompt_with_context, confirm_with_context

from .frame_extractor import FrameExtractor, FrameExtractionResult

logger = logging.getLogger(__name__)


@register_module
class QMMMPreparator(ProcessingModule):
    """Unified QM/MM preparation module for Gaussian (ONIOM) and ORCA."""

    NAME = "QM/MM Preparator"
    DESCRIPTION = "Generate QM/MM input files (Gaussian ONIOM or ORCA)"
    CATEGORY = "preparation"
    VERSION = "1.0.0"

    def __init__(self):
        self.processor = None
        self._oniom = None
        self._orca = None
        self._frame_result: Optional[FrameExtractionResult] = None
        self._active_interface: str = ""

    def get_workspace_requirements(self) -> List[str]:
        return ["detected_redox_sites"]

    def get_workspace_outputs(self) -> List[str]:
        return ["oniom_input_file", "orca_input_file", "orca_ff_file"]

    def get_menu_options(self) -> Dict[str, str]:
        return {
            "run_workflow": "Run QM/MM preparation workflow",
            "view_config": "View current QM/MM configuration",
            "reset": "Reset configuration to defaults",
        }

    def get_menu_suggestion(self, workspace):
        redox_sites = workspace.get("detected_redox_sites", [])
        if not redox_sites:
            return "No redox sites detected. Run RedoxSite Detector first."

        prmtop = workspace.get("parm7_file")
        if not prmtop:
            return "No prmtop file found. Run Topology Generator first."

        # Check if either interface has written output
        if workspace.get("oniom_input_file") or workspace.get("orca_input_file"):
            return "QM/MM input file written. Setup complete."

        return "Run the QM/MM workflow (option 1) to prepare input files."

    def handle_menu_option(self, option: str) -> bool:
        if option == "run_workflow":
            return self._run_workflow()
        elif option == "view_config":
            return self._view_config()
        elif option == "reset":
            return self._reset()
        return False

    # ------------------------------------------------------------------
    # Main workflow
    # ------------------------------------------------------------------

    def _run_workflow(self) -> bool:
        """Shared frame extraction followed by interface-specific workflow."""
        console = self.processor.console
        workspace = self._get_workspace_obj()

        # Phase 1: Shared frame extraction
        console.print("\n[bold cyan]===== QM/MM Preparation =====[/bold cyan]\n")

        extractor = FrameExtractor(self.processor, workspace)
        try:
            self._frame_result = extractor.run()
        except (FileNotFoundError, ValueError, RuntimeError) as e:
            console.print(f"\n[red]Frame extraction failed: {e}[/red]")
            return False

        # Phase 2: Interface selection
        console.print("\n[bold]Select QM/MM interface:[/bold]")
        console.print("  [cyan]1.[/cyan] Gaussian (ONIOM)")
        console.print("  [cyan]2.[/cyan] ORCA")

        choice = prompt_with_context(
            self.processor,
            "Select interface",
            choices=["1", "2"],
            default="1",
            module="QM/MM Preparator",
            description="QM/MM software interface",
            options_map={
                "1": "Gaussian (ONIOM)",
                "2": "ORCA",
            },
        )

        # Phase 3: Dispatch
        if choice == "1":
            return self._dispatch_oniom()
        else:
            return self._dispatch_orca()

    def _dispatch_oniom(self) -> bool:
        """Set up and run the Gaussian (ONIOM) workflow."""
        import parmed
        from proprep.oniom_prep.oniom_qmmm_preparator import ONIOMQMMMPreparator

        if not self._oniom:
            self._oniom = ONIOMQMMMPreparator()
            self._oniom.set_processor(self.processor)

        result = self._frame_result

        # Load parmed once and pass it down
        parm = parmed.load_file(result.prmtop_path, result.frame_pdb_paths[0])

        self._oniom.set_frame_data(
            result.prmtop_path,
            result.frame_pdb_paths,
            result.frame_labels,
            parm=parm,
        )
        self._active_interface = "oniom"
        return self._oniom._run_oniom_workflow()

    def _dispatch_orca(self) -> bool:
        """Set up and run the ORCA QM/MM workflow."""
        from proprep.orca_prep.orca_qmmm_preparator import ORCAQMMMPreparator

        if not self._orca:
            self._orca = ORCAQMMMPreparator()
            self._orca.set_processor(self.processor)

        result = self._frame_result
        self._orca.set_frame_data(
            result.prmtop_path,
            result.frame_rst7_paths,
            result.frame_labels,
            pdb_paths=result.frame_pdb_paths,
        )
        self._active_interface = "orca"
        return self._orca._run_orca_workflow()

    # ------------------------------------------------------------------
    # View / Reset
    # ------------------------------------------------------------------

    def _view_config(self) -> bool:
        """Delegate to the active interface's view command."""
        console = self.processor.console
        if self._active_interface == "oniom" and self._oniom:
            return self._oniom.view_current_setup()
        elif self._active_interface == "orca" and self._orca:
            return self._orca._view_current_config()
        else:
            console.print("[grey50]No QM/MM workflow has been run yet.[/grey50]")
            return True

    def _reset(self) -> bool:
        """Reset both interfaces and frame data."""
        self._oniom = None
        self._orca = None
        self._frame_result = None
        self._active_interface = ""
        self.processor.console.print("[green]QM/MM configuration reset.[/green]")
        return True

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _get_workspace_obj(self):
        """Get the workspace object."""
        if hasattr(self.processor, 'workspace'):
            return self.processor.workspace
        return self.processor._get_workspace()
