"""
Test module for AlphaFold Predictor command classes.

Comprehensive tests for all command pattern implementations.
"""

import unittest
from unittest.mock import Mock, MagicMock
from typing import Dict, Any

from proprep.structure_prep.alphafold_commands import (
    PredictStructureCommand,
    ConfigureParametersCommand,
    ViewResultsCommand,
    ExportBestModelCommand,
    CheckInstallationCommand,
)
from proprep.application.processor_command import ModuleActionCommand


class MockProcessor:
    """Mock processor for testing."""

    def __init__(self):
        self.current_module = Mock()
        self.current_module.predict_structure_interactive = Mock()
        self.current_module.configure_parameters = Mock()
        self.current_module.view_results = Mock()
        self.current_module.export_best_model = Mock()
        self.current_module.check_installation = Mock()
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


class TestPredictStructureCommand(unittest.TestCase):
    """Test PredictStructureCommand class."""

    def setUp(self):
        self.processor = MockProcessor()
        self.command = PredictStructureCommand(self.processor)

    def test_initialization(self):
        """Test command initialization."""
        self.assertIsInstance(self.command, ModuleActionCommand)
        self.assertEqual(self.command.processor, self.processor)
        self.assertEqual(self.command.module_name, "AlphaFold Predictor")
        self.assertEqual(self.command._action_to_run, "predict_structure_interactive")
        self.assertEqual(self.command.args, {})

    def test_execute_success(self):
        """Test successful command execution."""
        self.processor.current_module.predict_structure_interactive.return_value = True
        result = self.command.execute()
        self.assertTrue(result)
        self.processor.current_module.predict_structure_interactive.assert_called_once()

    def test_execute_failure(self):
        """Test command execution failure."""
        self.processor.current_module.predict_structure_interactive.return_value = False
        result = self.command.execute()
        self.assertFalse(result)

    def test_breadcrumb_generation(self):
        """Test breadcrumb generation."""
        expected = "AlphaFold Predictor.predict_structure_interactive"
        self.assertEqual(self.command.get_breadcrumb(), expected)


class TestConfigureParametersCommand(unittest.TestCase):
    """Test ConfigureParametersCommand class."""

    def setUp(self):
        self.processor = MockProcessor()
        self.command = ConfigureParametersCommand(self.processor)

    def test_initialization(self):
        """Test command initialization."""
        self.assertIsInstance(self.command, ModuleActionCommand)
        self.assertEqual(self.command.processor, self.processor)
        self.assertEqual(self.command.module_name, "AlphaFold Predictor")
        self.assertEqual(self.command._action_to_run, "configure_parameters")
        self.assertEqual(self.command.args, {})

    def test_execute_success(self):
        """Test successful command execution."""
        self.processor.current_module.configure_parameters.return_value = True
        result = self.command.execute()
        self.assertTrue(result)
        self.processor.current_module.configure_parameters.assert_called_once()

    def test_breadcrumb_generation(self):
        """Test breadcrumb generation."""
        expected = "AlphaFold Predictor.configure_parameters"
        self.assertEqual(self.command.get_breadcrumb(), expected)


class TestViewResultsCommand(unittest.TestCase):
    """Test ViewResultsCommand class."""

    def setUp(self):
        self.processor = MockProcessor()
        self.command = ViewResultsCommand(self.processor)

    def test_initialization(self):
        """Test command initialization."""
        self.assertIsInstance(self.command, ModuleActionCommand)
        self.assertEqual(self.command.processor, self.processor)
        self.assertEqual(self.command.module_name, "AlphaFold Predictor")
        self.assertEqual(self.command._action_to_run, "view_results")
        self.assertEqual(self.command.args, {})

    def test_execute_success(self):
        """Test successful command execution."""
        self.processor.current_module.view_results.return_value = True
        result = self.command.execute()
        self.assertTrue(result)
        self.processor.current_module.view_results.assert_called_once()

    def test_breadcrumb_generation(self):
        """Test breadcrumb generation."""
        expected = "AlphaFold Predictor.view_results"
        self.assertEqual(self.command.get_breadcrumb(), expected)


class TestExportBestModelCommand(unittest.TestCase):
    """Test ExportBestModelCommand class."""

    def setUp(self):
        self.processor = MockProcessor()

    def test_initialization_no_args(self):
        """Test command initialization without arguments."""
        command = ExportBestModelCommand(self.processor)
        self.assertIsInstance(command, ModuleActionCommand)
        self.assertEqual(command.processor, self.processor)
        self.assertEqual(command.module_name, "AlphaFold Predictor")
        self.assertEqual(command._action_to_run, "export_best_model")
        self.assertEqual(command.args, {"output_file": None})

    def test_initialization_with_args(self):
        """Test command initialization with arguments."""
        output_file = "/path/to/output.pdb"
        command = ExportBestModelCommand(self.processor, output_file=output_file)
        self.assertEqual(command.args, {"output_file": output_file})

    def test_execute_success(self):
        """Test successful command execution."""
        command = ExportBestModelCommand(self.processor, output_file="test.pdb")
        self.processor.current_module.export_best_model.return_value = True
        result = command.execute()
        self.assertTrue(result)
        self.processor.current_module.export_best_model.assert_called_once()

    def test_breadcrumb_generation(self):
        """Test breadcrumb generation."""
        command = ExportBestModelCommand(self.processor)
        expected = "AlphaFold Predictor.export_best_model"
        self.assertEqual(command.get_breadcrumb(), expected)


class TestCheckInstallationCommand(unittest.TestCase):
    """Test CheckInstallationCommand class."""

    def setUp(self):
        self.processor = MockProcessor()
        self.command = CheckInstallationCommand(self.processor)

    def test_initialization(self):
        """Test command initialization."""
        self.assertIsInstance(self.command, ModuleActionCommand)
        self.assertEqual(self.command.processor, self.processor)
        self.assertEqual(self.command.module_name, "AlphaFold Predictor")
        self.assertEqual(self.command._action_to_run, "check_installation")
        self.assertEqual(self.command.args, {})

    def test_execute_success(self):
        """Test successful command execution."""
        self.processor.current_module.check_installation.return_value = True
        result = self.command.execute()
        self.assertTrue(result)
        self.processor.current_module.check_installation.assert_called_once()

    def test_breadcrumb_generation(self):
        """Test breadcrumb generation."""
        expected = "AlphaFold Predictor.check_installation"
        self.assertEqual(self.command.get_breadcrumb(), expected)


class TestAlphaFoldCommandsIntegration(unittest.TestCase):
    """Integration tests for all AlphaFold commands."""

    def setUp(self):
        self.processor = MockProcessor()

    def test_all_commands_inherit_from_base(self):
        """Test all commands inherit from ModuleActionCommand."""
        commands = [
            PredictStructureCommand(self.processor),
            ConfigureParametersCommand(self.processor),
            ViewResultsCommand(self.processor),
            ExportBestModelCommand(self.processor),
            CheckInstallationCommand(self.processor),
        ]

        for command in commands:
            self.assertIsInstance(command, ModuleActionCommand)

    def test_all_commands_have_correct_module_name(self):
        """Test all commands have correct module name."""
        commands = [
            PredictStructureCommand(self.processor),
            ConfigureParametersCommand(self.processor),
            ViewResultsCommand(self.processor),
            ExportBestModelCommand(self.processor),
            CheckInstallationCommand(self.processor),
        ]

        for command in commands:
            self.assertEqual(command.module_name, "AlphaFold Predictor")

    def test_all_commands_have_unique_actions(self):
        """Test all commands have unique action names."""
        commands = [
            PredictStructureCommand(self.processor),
            ConfigureParametersCommand(self.processor),
            ViewResultsCommand(self.processor),
            ExportBestModelCommand(self.processor),
            CheckInstallationCommand(self.processor),
        ]

        action_names = [cmd._action_to_run for cmd in commands]
        self.assertEqual(len(action_names), len(set(action_names)),
                        "Action names should be unique")

    def test_command_execution_flow(self):
        """Test typical command execution flow."""
        # Configure parameters
        config_cmd = ConfigureParametersCommand(self.processor)
        self.processor.current_module.configure_parameters.return_value = True
        self.assertTrue(config_cmd.execute())

        # Predict structure
        predict_cmd = PredictStructureCommand(self.processor)
        self.processor.current_module.predict_structure_interactive.return_value = True
        self.assertTrue(predict_cmd.execute())

        # View results
        view_cmd = ViewResultsCommand(self.processor)
        self.processor.current_module.view_results.return_value = True
        self.assertTrue(view_cmd.execute())

        # Export model
        export_cmd = ExportBestModelCommand(self.processor, output_file="test.pdb")
        self.processor.current_module.export_best_model.return_value = True
        self.assertTrue(export_cmd.execute())


if __name__ == "__main__":
    unittest.main()
