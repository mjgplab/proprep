"""
Structure Alignment Commands

Commands for the ProPrep processor command pattern implementation.
"""

from proprep.application.processor_command import ModuleActionCommand


class AlignStructuresCommand(ModuleActionCommand):
    """Command to align multiple PDB structures interactively."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Structure Aligner",
            action_name="_align_structures_interactive",
            args={},
        )


class AlignOnRedoxSitesCommand(ModuleActionCommand):
    """Command to align structures using redox site residues."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Structure Aligner",
            action_name="_align_on_redox_sites_interactive",
            args={},
        )


class ViewAlignmentResultsCommand(ModuleActionCommand):
    """Command to view alignment results."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Structure Aligner",
            action_name="_view_alignment_results",
            args={},
        )
