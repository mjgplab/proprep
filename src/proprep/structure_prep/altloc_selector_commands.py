"""
AltLoc Selector Commands

Commands for the MPSA processor command pattern implementation for alternate location selection.
"""

from typing import Any, Dict, Optional

from proprep.application.processor_command import ModuleActionCommand


class IdentifyAlternateLocationsCommand(ModuleActionCommand):
    """Command to identify alternate locations in the structure."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="AltLoc Selector",
            action_name="identify_alternate_locations",
            args={},
        )


class SelectConformationsCommand(ModuleActionCommand):
    """Command to select conformations for alternate locations."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="AltLoc Selector",
            action_name="select_conformations",
            args={},
        )


class ProcessStructureCommand(ModuleActionCommand):
    """Command to process structure with selected conformations."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="AltLoc Selector",
            action_name="process_structure",
            args={},
        )


class ShowSelectionResultsCommand(ModuleActionCommand):
    """Command to show selection results."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="AltLoc Selector",
            action_name="show_selection_results",
            args={},
        )
