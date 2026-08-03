"""
AMBER MD Input Generator Commands

Commands for the MPSA processor command pattern implementation for AMBER input file generation operations.
"""

from proprep.application.processor_command import ModuleActionCommand


class GenerateInputCommand(ModuleActionCommand):
    """Command to generate MD input file."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="AMBER Input Generator",
            action_name="_generate_input_workflow",
            args={},
        )


class GenerateWorkflowInputsCommand(ModuleActionCommand):
    """Command to generate complete workflow input files."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="AMBER Input Generator",
            action_name="_generate_workflow_inputs",
            args={},
        )


class LoadTemplateCommand(ModuleActionCommand):
    """Command to load configuration from template."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="AMBER Input Generator",
            action_name="_load_template",
            args={},
        )


class SaveTemplateCommand(ModuleActionCommand):
    """Command to save current configuration as template."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="AMBER Input Generator",
            action_name="_save_template",
            args={},
        )


class ReviewConfigurationCommand(ModuleActionCommand):
    """Command to review current configuration."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="AMBER Input Generator",
            action_name="_review_configuration",
            args={},
        )


class ClearConfigurationCommand(ModuleActionCommand):
    """Command to clear configuration."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="AMBER Input Generator",
            action_name="_clear_configuration",
            args={},
        )

