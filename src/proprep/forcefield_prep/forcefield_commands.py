"""
Forcefield Commands

Commands for the MPSA processor command pattern implementation for forcefield operations.
"""

from typing import Any, Dict, Optional

from proprep.application.processor_command import ModuleActionCommand


class AnalyzeNonstandardResiduesCommand(ModuleActionCommand):
    """Command to analyze non-standard residues in structure."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Force Field Parameterizer",
            action_name="analyze_nonstandard_residues",
            args={},
        )

class ConfigureClassificationSettingsCommand(ModuleActionCommand):
    """Command to configure classification settings."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Force Field Parameterizer",
            action_name="configure_classification_settings",
            args={},
        )



class ParameterizeAACommand(ModuleActionCommand):
    """Command to parameterize a specific non-standard amino acid."""

    def __init__(self, processor, residue_name: Optional[str] = None):
        super().__init__(
            processor=processor,
            module_name="Force Field Parameterizer",
            action_name="parameterize_residue",  
            args={"residue_name": residue_name} if residue_name else {},
        )

class ManageMappingsCommand(ModuleActionCommand):
    """Command to manage amino acid and residue classification mappings."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Force Field Parameterizer",
            action_name="manage_mappings",
            args={},
        )

class ReclassifyResidueCommand(ModuleActionCommand):
    """Command to reclassify a residue to a different type."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Force Field Parameterizer",
            action_name="reclassify_residue",
            args={},
        )

class DisplayHelpCommand(ModuleActionCommand):
    """Command to display help information."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Force Field Parameterizer",
            action_name="display_help",
            args={},
        )

    def can_execute(self) -> bool:
        # Help is informational and never requires a loaded structure. The
        # base ModuleActionCommand gates on module.can_process() (which needs a
        # structure); override so Help runs from the menu even with nothing
        # loaded — matching its always-available (✓) status.
        return True
