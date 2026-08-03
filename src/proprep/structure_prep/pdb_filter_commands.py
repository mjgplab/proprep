"""
PDB Filter Commands

Commands for the MPSA processor command pattern implementation for PDB filtering operations.
"""

from typing import Any, Dict, Optional

from proprep.application.processor_command import ModuleActionCommand


class FilterPDBStructureCommand(ModuleActionCommand):
    """Command to filter PDB structure."""

    def __init__(self, processor, interactive=True):
        super().__init__(
            processor=processor,
            module_name="PDB Filter",
            action_name="filter_pdb_structure",
            args={"interactive": interactive},
        )


class ShowFilterStatusCommand(ModuleActionCommand):
    """Command to show filtering status."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="PDB Filter",
            action_name="show_filter_status",
            args={},
        )


class ExportFilterStatisticsCommand(ModuleActionCommand):
    """Command to export filter statistics."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="PDB Filter",
            action_name="export_filter_statistics",
            args={},
        )

