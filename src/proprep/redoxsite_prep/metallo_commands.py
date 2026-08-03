"""
Metal Module Commands

Commands for the MPSA processor command pattern implementation for metal site operations.
"""

from typing import Any, Dict, Optional

from proprep.application.processor_command import ModuleActionCommand

class DetectMetalSitesCommand(ModuleActionCommand):
    """Command to detect metal sites in structure."""

    def __init__(self, processor, interactive=True):
        super().__init__(
            processor=processor,
            module_name="Redox Site Preparer",
            action_name="detect_metal_sites",
            args={"interactive": interactive},
        )


class ImportMetalSitesCommand(ModuleActionCommand):
    """Command to import metal site information from file."""

    def __init__(self, processor, file_path: Optional[str] = None):
        super().__init__(
            processor=processor,
            module_name="Redox Site Preparer",
            action_name="import_metal_sites",
            args={"file_path": file_path} if file_path else {},
        )


class ViewMetalSitesCommand(ModuleActionCommand):
    """Command to view metal site analysis results."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Redox Site Preparer",
            action_name="view_metal_sites",
            args={},
        )


class ShowMetalSitesSummaryCommand(ModuleActionCommand):
    """Command to show summary of detected metal sites."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Redox Site Preparer",
            action_name="show_metal_sites_summary",
            args={},
        )


class ViewResultsCommand(ModuleActionCommand):
    """Command to view combined results (summary + detailed if requested)."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Redox Site Preparer",
            action_name="view_combined_results",
            args={},
        )


class ExportMetalSitesCommand(ModuleActionCommand):
    """Command to export detection results."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Redox Site Preparer",
            action_name="export_metal_sites",
            args={},
        )


class CreateCustomTransformerCommand(ModuleActionCommand):
    """Command to create a custom transformer for new metal site types."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Redox Site Preparer",
            action_name="create_custom_transformer",
            args={},
        )


class TransformMetalSitesCommand(ModuleActionCommand):
    """Command to transform metal sites for MD simulation."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Redox Site Preparer",
            action_name="transform_metal_sites",
            args={},
        )


class SaveTransformedStructureCommand(ModuleActionCommand):
    """Command to save the transformed structure."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Redox Site Preparer",
            action_name="save_transformed_structure",
            args={},
        )


class SaveFilteredStructureCommand(ModuleActionCommand):
    """Command to save filtered structure with coordinating residues."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Redox Site Preparer",
            action_name="save_structure",  # Note: this maps to the worker method
            args={},
        )


class GenerateSingleMicrostateCommand(ModuleActionCommand):
    """Command to generate a single redox microstate (transform + save)."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Redox Site Preparer",
            action_name="generate_single_microstate",
            args={},
        )


class GenerateMicrostatesCommand(ModuleActionCommand):
    """Command to show redox microstate menu for all microstates."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Redox Site Preparer",
            action_name="show_microstate_menu",
            args={},
        )
                

class GenerateTleapBondsCommand(ModuleActionCommand):
    """Command to generate TLEaP bond commands for metal sites."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Redox Site Preparer",
            action_name="generate_tleap_bonds",
            args={},
        )


class ExtractClustersCommand(ModuleActionCommand):
    """Command to extract active site cluster models."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Redox Site Preparer",
            action_name="extract_clusters",
            args={},
        )


class GenerateRestraintMasksCommand(ModuleActionCommand):
    """Command to generate restraint masks for MD minimization."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Redox Site Preparer",
            action_name="generate_restraint_masks",
            args={},
        )
