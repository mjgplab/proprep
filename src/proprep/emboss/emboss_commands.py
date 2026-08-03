"""
EMBOSS Commands

Command pattern implementations for EMBOSS module actions.
Each command encapsulates a specific EMBOSS operation following ProPrep's command pattern.
"""

from typing import Any, Dict, Optional

from proprep.application.processor_command import ModuleActionCommand


class SequenceAnalysisCommand(ModuleActionCommand):
    """Command to run EMBOSS sequence analysis tools (pepstats, pepinfo, charge, iep)"""

    def __init__(self, processor, interactive=True):
        super().__init__(
            processor=processor,
            module_name="EMBOSS Analysis",
            action_name="_sequence_analysis",
            args={"interactive": interactive},
        )


class PairwiseAlignmentCommand(ModuleActionCommand):
    """Command to run EMBOSS pairwise alignment tools (needle, water, stretcher)"""

    def __init__(self, processor, interactive=True):
        super().__init__(
            processor=processor,
            module_name="EMBOSS Analysis",
            action_name="_pairwise_alignment",
            args={"interactive": interactive},
        )


class PatternSearchCommand(ModuleActionCommand):
    """Command to run EMBOSS pattern/motif search tools (patmatmotifs, fuzzpro, sigcleave)"""

    def __init__(self, processor, interactive=True):
        super().__init__(
            processor=processor,
            module_name="EMBOSS Analysis",
            action_name="_pattern_search",
            args={"interactive": interactive},
        )


class BatchAnalysisCommand(ModuleActionCommand):
    """Command to run batch analysis on multiple sequences"""

    def __init__(self, processor, interactive=True):
        super().__init__(
            processor=processor,
            module_name="EMBOSS Analysis",
            action_name="_batch_analysis",
            args={"interactive": interactive},
        )


class ViewEMBOSSResultsCommand(ModuleActionCommand):
    """Command to view previous EMBOSS results"""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="EMBOSS Analysis",
            action_name="_view_results",
            args={},
        )


class ExportEMBOSSResultsCommand(ModuleActionCommand):
    """Command to export EMBOSS results to files"""

    def __init__(self, processor, format="text"):
        super().__init__(
            processor=processor,
            module_name="EMBOSS Analysis",
            action_name="_export_results",
            args={"format": format},
        )


# Integration Commands for other modules

class EMBOSSSequencePropertiesCommand(ModuleActionCommand):
    """Command to add EMBOSS sequence analysis to PDB Filter"""

    def __init__(self, processor, sequence_source="filtered"):
        super().__init__(
            processor=processor,
            module_name="EMBOSS Analysis",
            action_name="_analyze_sequence_properties",
            args={"sequence_source": sequence_source},
        )


class EMBOSSAlignmentAnalysisCommand(ModuleActionCommand):
    """Command to add EMBOSS alignment to Homology Searcher"""

    def __init__(self, processor, comparison_target="blast_hits"):
        super().__init__(
            processor=processor,
            module_name="EMBOSS Analysis",
            action_name="_alignment_analysis",
            args={"comparison_target": comparison_target},
        )


class EMBOSSChargeValidationCommand(ModuleActionCommand):
    """Command to add EMBOSS charge analysis to Protonation State Analyzer"""

    def __init__(self, processor, ph_range=(0, 14)):
        super().__init__(
            processor=processor,
            module_name="EMBOSS Analysis",
            action_name="_charge_validation",
            args={"ph_range": ph_range},
        )
