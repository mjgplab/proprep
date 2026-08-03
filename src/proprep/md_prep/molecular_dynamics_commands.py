"""
Molecular Dynamics Manager Commands

Commands for the ProPrep processor command pattern implementation.
"""

from proprep.application.processor_command import ModuleActionCommand


class SetupSingleSimulationCommand(ModuleActionCommand):
    """Command to set up single simulations."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Molecular Dynamics Manager",
            action_name="_setup_single_simulations",
            args={},
        )


class PerformSimulationsCommand(ModuleActionCommand):
    """Command to perform/execute simulations."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Molecular Dynamics Manager",
            action_name="_perform_simulations",
            args={},
        )


class ConfigureHardwareCommand(ModuleActionCommand):
    """Command to configure hardware settings."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Molecular Dynamics Manager",
            action_name="_configure_hardware",
            args={},
        )


class ImportTemplateCommand(ModuleActionCommand):
    """Command to import .mdin template files."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Molecular Dynamics Manager",
            action_name="_import_template",
            args={},
        )


class ImportWorkflowCommand(ModuleActionCommand):
    """Command to import .json workflow files."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Molecular Dynamics Manager",
            action_name="_import_workflow",
            args={},
        )


class AnalyzeSimulationsCommand(ModuleActionCommand):
    """Command to analyze completed simulations."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Molecular Dynamics Manager",
            action_name="_analyze_simulations",
            args={},
        )