"""
BLAST Commands

Commands for the MPSA processor command pattern implementation for BLAST operations.
"""

from proprep.application.processor_command import ModuleActionCommand


class RunBLASTCommand(ModuleActionCommand):
    """Command to run BLAST homology search."""

    def __init__(self, processor, interactive=True):
        super().__init__(
            processor=processor,
            module_name="Homology Searcher",
            action_name="run_blast_search",
            args={"interactive": interactive},
        )


class ViewBLASTResultsCommand(ModuleActionCommand):
    """Command to view BLAST results."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Homology Searcher",
            action_name="view_blast_results",
            args={},
        )


class ExportBLASTResultsCommand(ModuleActionCommand):
    """Command to export BLAST alignments and annotations."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Homology Searcher",
            action_name="export_blast_results",
            args={},
        )


class BuildHomologyModelCommand(ModuleActionCommand):
    """Command to build a homology model with MODELLER."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Homology Searcher",
            action_name="build_homology_model",
            args={},
        )
