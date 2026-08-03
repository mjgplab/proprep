"""
Disulfide Bond Commands

Commands for the MPSA processor command pattern implementation for disulfide bond operations.
"""

from typing import Any, Dict, Optional

from rich.prompt import Confirm, Prompt
from rich.status import Status
from rich.table import Table

from proprep.application.processor_command import ModuleActionCommand


class DetectDisulfideBondsCommand(ModuleActionCommand):
    """Command to detect disulfide bonds in the structure."""

    def __init__(self, processor, interactive=True, use_ssbond_records=True):
        super().__init__(
            processor=processor,
            module_name="Disulfide Bond Detector",
            action_name="detect_disulfide_bonds",
            args={"interactive": interactive, "use_ssbond_records": use_ssbond_records},
        )


class ViewDisulfideBondsCommand(ModuleActionCommand):
    """Command to view detected disulfide bonds."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Disulfide Bond Detector",
            action_name="view_disulfide_bonds",
            args={},
        )


class UpdateCysResiduesCommand(ModuleActionCommand):
    """Command to update CYS residues to CYX for disulfide bonds."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Disulfide Bond Detector",
            action_name="update_cys_residues",
            args={},
        )


class GenerateTLeapCommandsCommand(ModuleActionCommand):
    """Command to generate tLEaP commands for disulfide bonds."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Disulfide Bond Detector",
            action_name="generate_tleap_commands",
            args={},
        )
