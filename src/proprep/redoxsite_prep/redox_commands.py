#!/usr/bin/env python3
"""
RedoxSite Preparation Commands

Modern command implementations for the RedoxSite preparation system.
Replaces the legacy metallo_prep command structure.

Author: Claude Code Implementation  
"""

from typing import Optional
from proprep.application.processor_command import ProcessorCommand


class GenerateSingleRedoxMicrostateCommand(ProcessorCommand):
    """Command to generate a single redox microstate using the new RedoxSite transformation system."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Redox Site Preparer",
            action_name="generate_single_microstate",
            args={},
        )


class GenerateAllRedoxMicrostatesCommand(ProcessorCommand):
    """Command to generate all possible redox microstates in batch."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Redox Site Preparer", 
            action_name="generate_all_microstates",
            args={},
        )


class CreateCustomTransformerCommand(ProcessorCommand):
    """Command to create a custom transformer for new redox site types."""

    def __init__(self, processor):
        super().__init__(
            processor=processor,
            module_name="Redox Site Preparer",
            action_name="create_custom_transformer",
            args={},
        )