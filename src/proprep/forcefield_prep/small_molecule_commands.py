"""
Small Molecule Parameterizer Commands

Commands for the ProPrep processor command pattern implementation.
Handles the parameterization of non-covalently bound small molecules.
"""

from proprep.application.processor_command import ModuleActionCommand


class AnalyzeSmallMoleculesCommand(ModuleActionCommand):
    """Command to analyze structure for small molecules."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Small Molecule Parameterizer",
            action_name="_analyze_small_molecules",
            args={},
        )


class GenerateGaussianInputCommand(ModuleActionCommand):
    """Command to generate Gaussian input files for RESP charge calculation."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Small Molecule Parameterizer",
            action_name="_generate_gaussian_input",
            args={},
        )


class ProcessGaussianOutputCommand(ModuleActionCommand):
    """Command to process Gaussian output and generate MOL2 files."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Small Molecule Parameterizer",
            action_name="_process_gaussian_output",
            args={},
        )


class GenerateParametersCommand(ModuleActionCommand):
    """Command to generate force field parameters using parmchk2."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Small Molecule Parameterizer",
            action_name="_generate_parameters",
            args={},
        )


class CreateTLeapInputCommand(ModuleActionCommand):
    """Command to create tLEaP input files for gas and aqueous phases."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Small Molecule Parameterizer",
            action_name="_create_tleap_input",
            args={},
        )


class DisplayHelpCommand(ModuleActionCommand):
    """Command to display detailed help information."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Small Molecule Parameterizer",
            action_name="_display_help",
            args={},
        )


class RunWorkflowCommand(ModuleActionCommand):
    """Command to run the complete parameterization workflow."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Small Molecule Parameterizer",
            action_name="_run_complete_workflow",
            args={},
        )