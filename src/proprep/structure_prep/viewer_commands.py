"""
Interactive Structure Viewer Commands

Commands for the Interactive Structure Viewer module following ProPrep's
command pattern implementation.

Author: ProPrep Development Team
Date: 2025-11-14
"""

from proprep.application.processor_command import ModuleActionCommand


class LaunchViewerCommand(ModuleActionCommand):
    """Command to launch the interactive structure viewer."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Structure Viewer",
            action_name="_launch_viewer_workflow",
            args={},
        )


class ShowAnnotationInfoCommand(ModuleActionCommand):
    """Command to show available annotation information."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Structure Viewer",
            action_name="_show_annotation_info",
            args={},
        )


class SelectStructuresCommand(ModuleActionCommand):
    """Command to select structures for viewing."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Structure Viewer",
            action_name="_select_structures_for_viewing",
            args={},
        )


class ConfigureAnnotationsCommand(ModuleActionCommand):
    """Command to configure annotation display settings."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Structure Viewer",
            action_name="_configure_annotation_display",
            args={},
        )
