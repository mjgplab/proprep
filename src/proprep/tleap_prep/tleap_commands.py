"""
tLEaP Commands

Commands for the MPSA processor command pattern implementation for tLEaP input generation operations.
"""

from proprep.application.processor_command import ModuleActionCommand


class GatherBondDefinitionsCommand(ModuleActionCommand):
    """Command to gather bond definitions from all modules."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Topology Generator",
            action_name="gather_bond_definitions",
            args={},
        )


class EditCombinedBondsCommand(ModuleActionCommand):
    """Command to edit combined bond definitions."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Topology Generator",
            action_name="edit_combined_bonds",
            args={},
        )


class ConfigureTLeapParametersCommand(ModuleActionCommand):
    """Command to configure tLEaP parameters."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Topology Generator",
            action_name="configure_tleap_parameters",
            args={},
        )


class WriteTLeapInputFileCommand(ModuleActionCommand):
    """Command to write tLEaP input file."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Topology Generator",
            action_name="write_tleap_input_file",
            args={},
        )


class ViewBondSummaryCommand(ModuleActionCommand):
    """Command to display bond definitions summary."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Topology Generator",
            action_name="display_bond_summary",
            args={},
        )


class ValidateTLeapParametersCommand(ModuleActionCommand):
    """Command to validate tLEaP parameters."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Topology Generator",
            action_name="validate_tleap_parameters",
            args={},
        )


class ExportTLeapConfigCommand(ModuleActionCommand):
    """Command to export tLEaP configuration to file."""

    def __init__(self, processor, filename=None):
        super().__init__(
            processor=processor,
            module_name="Topology Generator",
            action_name="export_tleap_config",
            args={"filename": filename},
        )


class ImportTLeapConfigCommand(ModuleActionCommand):
    """Command to import tLEaP configuration from file."""

    def __init__(self, processor, filename=None):
        super().__init__(
            processor=processor,
            module_name="Topology Generator",
            action_name="import_tleap_config",
            args={"filename": filename},
        )


class GenerateSingleStateTLeapCommand(ModuleActionCommand):
    """Command to generate tLEaP input file for a single state."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Topology Generator",
            action_name="generate_single_state_tleap",
            args={},
        )


class GenerateMicrostateInputsCommand(ModuleActionCommand):
    """Command to generate tLEaP input files for all redox microstates."""

    def __init__(self, processor, metadata_file=None):
        super().__init__(
            processor=processor,
            module_name="Topology Generator",
            action_name="generate_microstate_inputs",
            args={"metadata_file": metadata_file},
        )


class GenerateTopologyFilesCommand(ModuleActionCommand):
    """Command to generate topology files."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Topology Generator",
            action_name="generate_topology_files",
            args={},
        )


class GenerateCpinCommand(ModuleActionCommand):
    """Command to generate the titration input file for constant pH / redox MD.

    Dispatches to a single flow that runs cpinutil.py, ceinutil.py or
    cpeinutil.py depending on the titration mode detected from the topology.
    The command name is unchanged so recorded sessions keep replaying.
    """

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Topology Generator",
            action_name="generate_cpin_file",
            args={},
        )


class RunPBTitrateCommand(ModuleActionCommand):
    """Command to refine protonation states via PB (PBSA)."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Topology Generator",
            action_name="run_pb_titrate",
            args={},
        )
