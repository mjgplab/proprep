"""
Integration Test for Workflow-Centric Step 1

Simple test to validate the new workflow-centric Step 1 implementation.
"""

import json
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

from rich.console import Console

from .workflow_centric_step1 import WorkflowCentricStep1Manager
from .user_data_manager import UserDataManager


def test_workflow_centric_step1_basic():
    """Basic integration test for the new workflow-centric Step 1."""
    
    # Create temporary directory for testing
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # Create mock console
        console = Console(file=open('/dev/null', 'w'))
        
        # Create user data manager with temp directory
        user_data_manager = UserDataManager(package_dir=temp_path, console=console)
        
        # Create mock workspace
        workspace = Mock()
        workspace.get.return_value = []
        workspace.set.return_value = None
        
        # Create test structure pairs
        test_structure_pairs = [
            {
                'name': 'protein_1abc',
                'prmtop_path': '/path/to/1abc.prmtop',
                'rst7_path': '/path/to/1abc.rst7'
            },
            {
                'name': 'protein_2def',
                'prmtop_path': '/path/to/2def.prmtop', 
                'rst7_path': '/path/to/2def.rst7'
            }
        ]
        
        # Create some basic workflow presets for testing
        _create_test_workflow_presets(user_data_manager)
        
        # Create Step 1 manager
        step1_manager = WorkflowCentricStep1Manager(
            console=console,
            user_data_manager=user_data_manager,
            workspace=workspace
        )
        
        print("✓ WorkflowCentricStep1Manager created successfully")
        print("✓ Basic integration test passed")
        
        return True


def _create_test_workflow_presets(user_data_manager: UserDataManager):
    """Create basic test workflow presets."""
    
    # Create test template
    test_template = {
        'name': 'Test Minimization',
        'description': 'Basic minimization template',
        'type': 'minimization',
        'config': {
            'imin': 1,
            'ntx': 1,
            'maxcyc': 1000,
            'ncyc': 500,
            'cut': 10.0,
            'ntpr': 100
        },
        'created_date': '2025-01-15T10:00:00'
    }
    
    template_id = user_data_manager.save_custom_template(test_template)
    
    # Create test workflow
    test_workflow = {
        'name': 'Test Protein Workflow',
        'description': 'Basic protein workflow for testing',
        'category': 'Testing',
        'version': '1.0',
        'author': 'Test',
        'steps': [
            {
                'id': 'step_1',
                'name': 'Minimization',
                'type': 'minimization',
                'template_ref': template_id,
                'description': 'Energy minimization step',
                'dependencies': [],
                'input_coord': 'inpcrd',
                'output_coord': 'step_1.rst',
                'parameter_overrides': test_template.get('config', {})
            }
        ],
        'created_date': '2025-01-15T10:00:00'
    }
    
    workflow_id = user_data_manager.save_custom_workflow(test_workflow)
    
    print(f"✓ Created test template: {template_id}")
    print(f"✓ Created test workflow: {workflow_id}")


def run_integration_tests():
    """Run all integration tests."""
    print("Running Workflow-Centric Step 1 Integration Tests...")
    print("=" * 50)
    
    try:
        test_workflow_centric_step1_basic()
        print("\n✅ All integration tests passed!")
        return True
    except Exception as e:
        print(f"\n❌ Integration tests failed: {e}")
        return False


if __name__ == "__main__":
    run_integration_tests()