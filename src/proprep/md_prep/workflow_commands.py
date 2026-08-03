"""
AMBER Workflow Manager Commands

Commands for the ProPrep processor command pattern implementation.
"""

from proprep.application.processor_command import ModuleActionCommand


class StartWorkflowCommand(ModuleActionCommand):
    """Command to start new workflow execution."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="AMBER Workflow Manager",
            action_name="_start_new_workflow",
            args={},
        )


class ResumeWorkflowCommand(ModuleActionCommand):
    """Command to resume workflow from checkpoint."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="AMBER Workflow Manager",
            action_name="_resume_workflow",
            args={},
        )


class MonitorSimulationCommand(ModuleActionCommand):
    """Command to monitor running simulation."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="AMBER Workflow Manager",
            action_name="_monitor_current_simulation",
            args={},
        )


class AnalyzeHistoryCommand(ModuleActionCommand):
    """Command to analyze historical simulation data."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="AMBER Workflow Manager",
            action_name="_analyze_historical_data",
            args={},
        )


class ConfigureSettingsCommand(ModuleActionCommand):
    """Command to configure execution settings."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="AMBER Workflow Manager",
            action_name="_configure_settings",
            args={},
        )


class ManageCheckpointsCommand(ModuleActionCommand):
    """Command to manage checkpoints."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="AMBER Workflow Manager",
            action_name="_manage_checkpoints",
            args={},
        )