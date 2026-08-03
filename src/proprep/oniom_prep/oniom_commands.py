"""
ONIOM QM/MM Commands

Commands for the ProPrep processor command pattern implementation for ONIOM input generation operations.
"""

from proprep.application.processor_command import ModuleActionCommand


class SelectRedoxSiteCommand(ModuleActionCommand):
    """Command to select a RedoxSite for ONIOM calculation."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="ONIOM QM/MM Preparator",
            action_name="select_redox_site",
            args={},
        )


class ConfigureLayersCommand(ModuleActionCommand):
    """Command to configure ONIOM layer assignments."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="ONIOM QM/MM Preparator",
            action_name="configure_layers",
            args={},
        )


class ConfigureQMSettingsCommand(ModuleActionCommand):
    """Command to configure QM calculation settings."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="ONIOM QM/MM Preparator",
            action_name="configure_qm_settings",
            args={},
        )


class ConfigureMMSettingsCommand(ModuleActionCommand):
    """Command to configure MM force field settings."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="ONIOM QM/MM Preparator",
            action_name="configure_mm_settings",
            args={},
        )


class ConfigureJobSettingsCommand(ModuleActionCommand):
    """Command to configure job execution settings."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="ONIOM QM/MM Preparator",
            action_name="configure_job_settings",
            args={},
        )


class PrepareONIOMSetupCommand(ModuleActionCommand):
    """Command to prepare complete ONIOM setup."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="ONIOM QM/MM Preparator",
            action_name="prepare_oniom_setup",
            args={},
        )


class ViewValidationReportCommand(ModuleActionCommand):
    """Command to view detailed validation report."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="ONIOM QM/MM Preparator",
            action_name="view_validation_report",
            args={},
        )


class WriteONIOMInputCommand(ModuleActionCommand):
    """Command to write Gaussian ONIOM input file."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="ONIOM QM/MM Preparator",
            action_name="write_oniom_input",
            args={},
        )


class WriteDiagnosticCommand(ModuleActionCommand):
    """Command to write diagnostic ONIOM input file (InputFiles + Test)."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="ONIOM QM/MM Preparator",
            action_name="write_diagnostic_input",
            args={},
        )


class WriteSummaryReportCommand(ModuleActionCommand):
    """Command to write ONIOM setup summary report."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="ONIOM QM/MM Preparator",
            action_name="write_summary_report",
            args={},
        )


class ViewCurrentSetupCommand(ModuleActionCommand):
    """Command to view current ONIOM configuration."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="ONIOM QM/MM Preparator",
            action_name="view_current_setup",
            args={},
        )


class ResetConfigurationCommand(ModuleActionCommand):
    """Command to reset ONIOM configuration to defaults."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="ONIOM QM/MM Preparator",
            action_name="reset_configuration",
            args={},
        )
