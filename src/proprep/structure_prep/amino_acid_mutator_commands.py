"""
Amino Acid Mutator Commands

Commands for the MPSA processor command pattern implementation for amino acid mutations.
"""

from typing import Any, Dict, Optional

from proprep.application.processor_command import ModuleActionCommand


class DisplaySequenceCommand(ModuleActionCommand):
    """Command to display protein sequences with numbering."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Amino Acid Mutator",
            action_name="_display_sequences",
            args={},
        )


class SpecifyMutationsCommand(ModuleActionCommand):
    """Command to specify amino acid mutations interactively."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Amino Acid Mutator",
            action_name="_specify_mutations",
            args={},
        )


class ApplyMutationsCommand(ModuleActionCommand):
    """Command to store mutations in workspace for Structure Completeness."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Amino Acid Mutator",
            action_name="_apply_mutations",
            args={},
        )


class ClearMutationsCommand(ModuleActionCommand):
    """Command to clear all pending mutations."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Amino Acid Mutator",
            action_name="_clear_mutations",
            args={},
        )

