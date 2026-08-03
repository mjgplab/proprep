"""
AlphaFold Commands

Commands for the AlphaFold predictor module using the command pattern.
"""

from proprep.application.processor_command import ModuleActionCommand


class PredictStructureCommand(ModuleActionCommand):
    """Command to run AlphaFold structure prediction."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="AlphaFold Predictor",
            action_name="predict_structure_interactive",
            args={},
        )


class ConfigureParametersCommand(ModuleActionCommand):
    """Command to configure AlphaFold prediction parameters."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="AlphaFold Predictor",
            action_name="configure_parameters",
            args={},
        )


class ViewResultsCommand(ModuleActionCommand):
    """Command to view AlphaFold prediction results."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="AlphaFold Predictor",
            action_name="view_results",
            args={},
        )


class ExportBestModelCommand(ModuleActionCommand):
    """Command to export the best AlphaFold model."""

    def __init__(self, processor, output_file=None):
        super().__init__(
            processor=processor,
            module_name="AlphaFold Predictor",
            action_name="export_best_model",
            args={"output_file": output_file},
        )


class CheckInstallationCommand(ModuleActionCommand):
    """Command to check AlphaFold/ColabFold installation status."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="AlphaFold Predictor",
            action_name="check_installation",
            args={},
        )
