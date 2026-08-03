"""
Protonation Commands

Commands for the MPSA processor command pattern implementation for protonation state analysis.
"""

from typing import Any, Dict, Optional

from proprep.application.processor_command import ModuleActionCommand


class AnalyzeProtonationStatesCommand(ModuleActionCommand):
    """Command to analyze protonation states of titratable residues."""

    def __init__(self, processor, interactive=True):
        super().__init__(
            processor=processor,
            module_name="Protonation State Analyzer",
            action_name="analyze_protonation_states",
            args={"interactive": interactive},
        )


class RunPropkaCommand(ModuleActionCommand):
    """Command to run PROPKA for improved pKa prediction."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Protonation State Analyzer",
            action_name="run_propka",
            args={},
        )


class ViewProtonationResultsCommand(ModuleActionCommand):
    """Command to view protonation analysis results."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Protonation State Analyzer",
            action_name="view_results",
            args={},
        )


class SetMDResidueNamesCommand(ModuleActionCommand):
    """Command to set residue names for MD simulation."""

    def __init__(self, processor, constant_ph=False):
        super().__init__(
            processor=processor,
            module_name="Protonation State Analyzer",
            action_name="set_residue_names",
            args={"constant_ph": constant_ph},
        )


class UpdateStructureCommand(ModuleActionCommand):
    """Command to update structure with MD residue names."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Protonation State Analyzer",
            action_name="update_structure",
            args={},
        )


class SaveProtonationResultsCommand(ModuleActionCommand):
    """Command to save analysis results to file."""

    def __init__(self, processor, output_format="text"):
        super().__init__(
            processor=processor,
            module_name="Protonation State Analyzer",
            action_name="save_results",
            args={"output_format": output_format},
        )


class GeneratePHProfileCommand(ModuleActionCommand):
    """Command to generate pH titration profile."""

    def __init__(self, processor, ph_range=(0, 14), step=0.5):
        super().__init__(
            processor=processor,
            module_name="Protonation State Analyzer",
            action_name="generate_ph_profile",
            args={"ph_range": ph_range, "step": step},
        )


class CalculateIdealCapacitanceCommand(ModuleActionCommand):
    """Command to calculate ideal charge capacitance."""
    
    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Protonation State Analyzer",
            action_name="calculate_capacitance",
            args={},  # Empty args, will prompt in the action method
        )


class GenerateCapacitanceProfileCommand(ModuleActionCommand):
    """Command to generate capacitance profile over pH range."""
    
    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Protonation State Analyzer",
            action_name="generate_capacitance_profile",
            args={},  # Empty args, will prompt in the action method
        )