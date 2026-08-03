"""
MD Restraint Commands

Commands for the ProPrep processor command pattern implementation for MD restraint management.
"""

from proprep.application.processor_command import ModuleActionCommand


class ConfigureRestraintsCommand(ModuleActionCommand):
    """Command to configure MD restraints interactively."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="MD Restraint Manager",
            action_name="_configure_restraints",
            args={},
        )


class DisplayRestraintsCommand(ModuleActionCommand):
    """Command to display current MD restraints."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="MD Restraint Manager",
            action_name="_display_restraints",
            args={},
        )


class ImportRestraintsCommand(ModuleActionCommand):
    """Command to import restraints from file."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="MD Restraint Manager",
            action_name="_import_restraints",
            args={},
        )


class ExportDisangCommand(ModuleActionCommand):
    """Command to export restraints to DISANG file."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="MD Restraint Manager",
            action_name="_export_disang",
            args={},
        )


class GenerateRestraintMaskCommand(ModuleActionCommand):
    """Command to generate restraint mask for MD minimization."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="MD Restraint Manager",
            action_name="_generate_restraint_mask",
            args={},
        )