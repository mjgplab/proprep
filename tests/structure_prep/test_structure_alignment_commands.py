"""
Test module for Structure Aligner command classes.

Comprehensive tests for all command pattern implementations.
"""

import unittest
from unittest.mock import Mock
from typing import Dict, Any

from proprep.structure_prep.structure_alignment_commands import (
    AlignStructuresCommand,
    AlignOnRedoxSitesCommand,
    ViewAlignmentResultsCommand,
)
from proprep.application.processor_command import ModuleActionCommand


class MockProcessor:
    """Mock processor for testing."""

    def __init__(self):
        self.current_module = Mock()
        self.current_module._align_structures_interactive = Mock()
        self.current_module._align_on_redox_sites_interactive = Mock()
        self.current_module._view_alignment_results = Mock()
        # ModuleActionCommand.execute iterates this to snapshot workspace state
        # before running, so a bare Mock (not iterable) fails the call rather
        # than the assertion under test.
        self.current_module.get_workspace_requirements = Mock(return_value=[])

        self._workspace = {}
        self.console = Mock()
        self.command_history = Mock()

    def get_module_instance(self, module_name):
        return self.current_module

    def _get_workspace(self):
        return self._workspace


class TestAlignStructuresCommand(unittest.TestCase):
    """Test AlignStructuresCommand class."""

    def setUp(self):
        self.processor = MockProcessor()
        self.command = AlignStructuresCommand(self.processor)

    def test_initialization(self):
        """Test command initialization."""
        self.assertIsInstance(self.command, ModuleActionCommand)
        self.assertEqual(self.command.processor, self.processor)
        self.assertEqual(self.command.module_name, "Structure Aligner")
        self.assertEqual(self.command._action_to_run, "_align_structures_interactive")
        self.assertEqual(self.command.args, {})

    def test_execute_success(self):
        """Test successful command execution."""
        self.processor.current_module._align_structures_interactive.return_value = True
        result = self.command.execute()
        self.assertTrue(result)
        self.processor.current_module._align_structures_interactive.assert_called_once()

    def test_execute_failure(self):
        """Test command execution failure."""
        self.processor.current_module._align_structures_interactive.return_value = False
        result = self.command.execute()
        self.assertFalse(result)

    def test_breadcrumb_generation(self):
        """Test breadcrumb generation."""
        expected = "Structure Aligner._align_structures_interactive"
        self.assertEqual(self.command.get_breadcrumb(), expected)


class TestAlignOnRedoxSitesCommand(unittest.TestCase):
    """Test AlignOnRedoxSitesCommand class."""

    def setUp(self):
        self.processor = MockProcessor()
        self.command = AlignOnRedoxSitesCommand(self.processor)

    def test_initialization(self):
        """Test command initialization."""
        self.assertIsInstance(self.command, ModuleActionCommand)
        self.assertEqual(self.command.processor, self.processor)
        self.assertEqual(self.command.module_name, "Structure Aligner")
        self.assertEqual(self.command._action_to_run, "_align_on_redox_sites_interactive")
        self.assertEqual(self.command.args, {})

    def test_execute_success(self):
        """Test successful command execution."""
        self.processor.current_module._align_on_redox_sites_interactive.return_value = True
        result = self.command.execute()
        self.assertTrue(result)
        self.processor.current_module._align_on_redox_sites_interactive.assert_called_once()

    def test_execute_failure(self):
        """Test command execution failure."""
        self.processor.current_module._align_on_redox_sites_interactive.return_value = False
        result = self.command.execute()
        self.assertFalse(result)

    def test_breadcrumb_generation(self):
        """Test breadcrumb generation."""
        expected = "Structure Aligner._align_on_redox_sites_interactive"
        self.assertEqual(self.command.get_breadcrumb(), expected)


class TestViewAlignmentResultsCommand(unittest.TestCase):
    """Test ViewAlignmentResultsCommand class."""

    def setUp(self):
        self.processor = MockProcessor()
        self.command = ViewAlignmentResultsCommand(self.processor)

    def test_initialization(self):
        """Test command initialization."""
        self.assertIsInstance(self.command, ModuleActionCommand)
        self.assertEqual(self.command.processor, self.processor)
        self.assertEqual(self.command.module_name, "Structure Aligner")
        self.assertEqual(self.command._action_to_run, "_view_alignment_results")
        self.assertEqual(self.command.args, {})

    def test_execute_success(self):
        """Test successful command execution."""
        self.processor.current_module._view_alignment_results.return_value = True
        result = self.command.execute()
        self.assertTrue(result)
        self.processor.current_module._view_alignment_results.assert_called_once()

    def test_execute_failure(self):
        """Test command execution failure."""
        self.processor.current_module._view_alignment_results.return_value = False
        result = self.command.execute()
        self.assertFalse(result)

    def test_breadcrumb_generation(self):
        """Test breadcrumb generation."""
        expected = "Structure Aligner._view_alignment_results"
        self.assertEqual(self.command.get_breadcrumb(), expected)


class TestStructureAlignmentIntegration(unittest.TestCase):
    """Integration tests for all structure alignment commands."""

    def setUp(self):
        self.processor = MockProcessor()

    def test_all_commands_inherit_from_base(self):
        """Test all commands inherit from ModuleActionCommand."""
        commands = [
            AlignStructuresCommand(self.processor),
            AlignOnRedoxSitesCommand(self.processor),
            ViewAlignmentResultsCommand(self.processor),
        ]

        for command in commands:
            self.assertIsInstance(command, ModuleActionCommand)

    def test_all_commands_have_correct_module_name(self):
        """Test all commands have correct module name."""
        commands = [
            AlignStructuresCommand(self.processor),
            AlignOnRedoxSitesCommand(self.processor),
            ViewAlignmentResultsCommand(self.processor),
        ]

        for command in commands:
            self.assertEqual(command.module_name, "Structure Aligner")

    def test_all_commands_have_unique_actions(self):
        """Test all commands have unique action names."""
        commands = [
            AlignStructuresCommand(self.processor),
            AlignOnRedoxSitesCommand(self.processor),
            ViewAlignmentResultsCommand(self.processor),
        ]

        action_names = [cmd._action_to_run for cmd in commands]
        self.assertEqual(len(action_names), len(set(action_names)),
                        "Action names should be unique")


if __name__ == "__main__":
    unittest.main()
