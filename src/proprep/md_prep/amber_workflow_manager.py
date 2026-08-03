"""
AMBER Workflow Manager Module

Complete AMBER workflow orchestration and monitoring for ProPrep.
Integrates all functionality from the original amber_workflow.py script.
"""

import os
import time
import glob
from pathlib import Path
from typing import Dict, Any, List, Optional

from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from proprep.utils.prompts import (
    prompt_with_context,
    confirm_with_context,
    int_prompt_with_context,
)
from proprep.utils.file_browser import remap_recorded_index, annotate_selected_path

from proprep.utils.module_registry import ProcessingModule, register_module
from proprep.md_prep.workflow_commands import (
    StartWorkflowCommand,
    MonitorSimulationCommand,
    ResumeWorkflowCommand,
    AnalyzeHistoryCommand,
    ConfigureSettingsCommand,
    ManageCheckpointsCommand,
)

from .amber_workflow_components import (
    AMBERMonitor,
    CheckpointManager,
    WorkflowExecutor,
    AMBERWorkflowCore
)


# Temporarily disabled - replaced by MolecularDynamicsManager
# @register_module 
class AmberWorkflowManager(ProcessingModule):
    """Complete AMBER workflow orchestration and monitoring."""

    NAME = "AMBER Workflow Manager"
    CATEGORY = "MD Execution"
    DESCRIPTION = "Execute and monitor complete AMBER MD workflows"
    VERSION = "2.0.0"
    REQUIRES = ["AMBER Input Generator"]
    PRIORITY = 10  # Run after input generation

    def __init__(self):
        super().__init__()
        
        # Core workflow functionality from original amber_workflow.py
        self.workflow_core = AMBERWorkflowCore()
        
    def set_processor(self, processor):
        """Set the processor reference."""
        self.processor = processor
        # Propagate to the core so its interactive prompts record context
        if hasattr(self, "workflow_core") and self.workflow_core is not None:
            self.workflow_core.processor = processor

    @property
    def console(self):
        """Get console from processor if available."""
        if (
            hasattr(self, "processor")
            and self.processor
            and hasattr(self.processor, "console")
        ):
            return self.processor.console
        else:
            from rich.console import Console
            return Console()

    def get_menu_options(self) -> Dict[str, str]:
        """Get available menu options."""
        options = {
            "start_workflow": "Start new workflow execution",
            "monitor": "Monitor running simulation",
            "analyze_history": "Analyze completed simulations", 
            "configure": "Configure execution settings",
            "help": "Show workflow manager help",
        }
        
        # Add resume option if checkpoint exists
        if self.workflow_core.checkpoint_manager.checkpoint_exists():
            options["resume"] = "Resume from checkpoint"
            
        return options

    def get_workspace_requirements(self) -> List[str]:
        """List what this module needs from workspace."""
        return []  # Can work independently

    def can_process(self, workspace) -> bool:
        """Check if module can process current workspace."""
        return True  # Always available

    def handle_menu_option(self, option: str) -> bool:
        """Handle menu option selection."""
        if option == "start_workflow":
            command = StartWorkflowCommand(self.processor)
            return command.execute()
        elif option == "resume":
            command = ResumeWorkflowCommand(self.processor)
            return command.execute()
        elif option == "monitor":
            command = MonitorSimulationCommand(self.processor)
            return command.execute()
        elif option == "analyze_history":
            command = AnalyzeHistoryCommand(self.processor)
            return command.execute()
        elif option == "configure":
            command = ConfigureSettingsCommand(self.processor)
            return command.execute()
        elif option == "help":
            return self._show_help()
        return False

    def _start_new_workflow(self) -> bool:
        """Start a new workflow execution."""
        self.console.print("[bold]🚀 Start New AMBER Workflow[/bold]")
        
        # Check for existing checkpoint
        if self.workflow_core.checkpoint_manager.checkpoint_exists():
            overwrite = confirm_with_context(
                self.processor,
                "⚠️  Existing checkpoint found. Overwrite?",
                default=False,
                module="AMBER Workflow Manager",
                description="Overwrite existing checkpoint",
            )
            if not overwrite:
                self.console.print("Use 'resume' option to continue existing workflow.")
                return False
            else:
                self.workflow_core.checkpoint_manager.delete_checkpoint()
        
        # Verify required files
        if not self._verify_system_files():
            return False
        
        # Verify input files
        if not self._verify_input_files():
            return False
        
        # Configure execution settings
        if not self._quick_configure():
            return False
        
        # Choose execution mode
        self.console.print("\n⚙️  Execution Mode")
        self.console.print("-" * 20)
        self.console.print("1. Interactive (step-by-step with monitoring)", highlight=False)
        self.console.print("2. Continuous (run all steps automatically)", highlight=False)
        
        mode_choice = prompt_with_context(
            self.processor,
            "Choose mode [1-2]",
            choices=["1", "2"],
            default="1",
            module="AMBER Workflow Manager",
            description="Workflow mode selection",
        )
        
        if mode_choice == "1":
            return self._run_interactive_workflow()
        else:
            return self._run_continuous_workflow()

    def _verify_system_files(self) -> bool:
        """Verify required system files exist."""
        # Look for topology files
        prmtop_files = list(Path(".").glob("*.prmtop")) + list(Path(".").glob("*.parm7"))
        if not prmtop_files:
            self.console.print("❌ No topology file (.prmtop/.parm7) found!")
            return False
        
        # Look for coordinate files  
        coord_files = list(Path(".").glob("*.inpcrd")) + list(Path(".").glob("*.rst7"))
        if not coord_files:
            self.console.print("❌ No coordinate file (.inpcrd/.rst7) found!")
            return False
        
        self.console.print("✅ System files verified")
        return True

    def _verify_input_files(self) -> bool:
        """Verify input files for workflow steps exist."""
        missing_files = []
        
        input_files = {
            "minimization": "min.in",
            "heating": "heat.in", 
            "equilibration_1": "equil_1.in",
            "equilibration_2": "equil_2.in",
            "production": "prod_1.in"
        }
        
        for step, filename in input_files.items():
            if not Path(filename).exists():
                missing_files.append(filename)
        
        if missing_files:
            self.console.print(f"❌ Missing input files: {', '.join(missing_files)}")
            self.console.print("💡 Use 'AMBER Input Generator' to create workflow inputs")
            return False
        
        self.console.print("✅ Input files verified")
        return True

    def _quick_configure(self) -> bool:
        """Quick configuration of execution settings."""
        self.console.print("\n🖥️  Quick Configuration")
        self.console.print("-" * 25)
        
        # MPI settings
        self.workflow_core.use_mpi = confirm_with_context(
            self.processor,
            "Use MPI for CPU steps?",
            default=True,
            module="AMBER Workflow Manager",
            description="Use MPI for CPU steps",
        )
        
        if self.workflow_core.use_mpi:
            self.workflow_core.num_cores = int_prompt_with_context(
                self.processor,
                "Number of MPI cores",
                default=32,
                module="AMBER Workflow Manager",
                description="Number of MPI cores",
            )
        
        # GPU settings
        self.workflow_core.gpu_device = int_prompt_with_context(
            self.processor,
            "GPU device ID",
            default=0,
            module="AMBER Workflow Manager",
            description="GPU device ID",
        )
        
        self.console.print(f"✅ Configuration: MPI={self.workflow_core.use_mpi}")
        if self.workflow_core.use_mpi:
            self.console.print(f"   Cores: {self.workflow_core.num_cores}")
        self.console.print(f"   GPU: {self.workflow_core.gpu_device}")
        
        return True

    def _run_interactive_workflow(self) -> bool:
        """Run workflow interactively with step-by-step monitoring."""
        self.console.print("\n🎮 Interactive Workflow Mode")
        self.console.print("=" * 50)
        
        for i, step in enumerate(self.workflow_core.workflow_steps):
            if not self._execute_step_interactive(step, i):
                return False
        
        self.console.print("\n🎉 Workflow completed successfully!")
        self._offer_final_analysis()
        return True

    def _execute_step_interactive(self, step: str, step_index: int) -> bool:
        """Execute a single step interactively."""
        self.console.print(f"\n{'='*70}")
        self.console.print(f"Step {step_index+1}/{len(self.workflow_core.workflow_steps)}: {step.replace('_', ' ').title()}")
        self.console.print(f"{'='*70}")
        
        # Build command
        command = self.workflow_core.workflow_executor.build_amber_command(
            step, self.workflow_core.use_mpi, self.workflow_core.num_cores, self.workflow_core.gpu_device
        )
        self.console.print(f"🔧 Command: {command}")
        
        # Save checkpoint
        self.workflow_core.checkpoint_manager.save_checkpoint(
            step_index, step, self.workflow_core.completed_steps, 
            "system.prmtop", "system.inpcrd", f"{step}.rst7", 1
        )
        
        # Ask user if ready to proceed
        proceed = confirm_with_context(
            self.processor,
            f"🚀 Execute {step}?",
            default=True,
            module="AMBER Workflow Manager",
            description=f"Execute step {step}",
        )
        if not proceed:
            self.console.print("⏭️  Skipped step")
            return True
        
        # Start simulation
        output_file = f"{step}.out"
        if not self.workflow_core.workflow_executor.start_simulation_background(command, step, output_file):
            return False
        
        # Setup monitor
        self.workflow_core.current_monitor = AMBERMonitor(output_file)
        
        # Interactive monitoring using original functionality
        result = self.workflow_core.interactive_monitor(step)
        if result:
            self.workflow_core.completed_steps.append(step)
        
        return result

    def _run_continuous_workflow(self) -> bool:
        """Run workflow continuously without interaction."""
        self.console.print("\n🔄 Continuous Workflow Mode")
        self.console.print("=" * 50)
        
        for i, step in enumerate(self.workflow_core.workflow_steps):
            self.console.print(f"\n📍 Step {i+1}/{len(self.workflow_core.workflow_steps)}: {step.replace('_', ' ').title()}")
            
            # Build command
            command = self.workflow_core.workflow_executor.build_amber_command(
                step, self.workflow_core.use_mpi, self.workflow_core.num_cores, self.workflow_core.gpu_device
            )
            
            # Save checkpoint
            self.workflow_core.checkpoint_manager.save_checkpoint(
                i, step, self.workflow_core.completed_steps, 
                "system.prmtop", "system.inpcrd", f"{step}.rst7", 2
            )
            
            # Start simulation
            output_file = f"{step}.out"
            self.console.print(f"🚀 Starting {step}...")
            
            if not self.workflow_core.workflow_executor.start_simulation_background(command, step, output_file):
                return False
            
            # Wait for completion
            while self.workflow_core.workflow_executor.is_simulation_running():
                time.sleep(5)
                self.console.print(".", end="", flush=True)
            
            # Check result
            is_running, return_code = self.workflow_core.workflow_executor.get_simulation_status()
            if return_code == 0:
                self.console.print(f"\n✅ {step} completed!")
                self.workflow_core.completed_steps.append(step)
            else:
                self.console.print(f"\n❌ {step} failed with return code {return_code}!")
                return False
        
        self.console.print("\n🎉 Continuous workflow completed successfully!")
        self._offer_final_analysis()
        return True

    def _resume_workflow(self) -> bool:
        """Resume workflow from checkpoint."""
        checkpoint = self.workflow_core.checkpoint_manager.load_checkpoint()
        if not checkpoint:
            self.console.print("❌ No checkpoint found to resume from")
            return False
        
        self.console.print(f"📄 Resuming from step: {checkpoint['step_name']}")
        self.workflow_core.completed_steps = checkpoint.get('completed_steps', [])
        
        # Offer to view historical plots
        if self.workflow_core.completed_steps:
            self.console.print(f"📊 Found {len(self.workflow_core.completed_steps)} completed step(s): {', '.join(self.workflow_core.completed_steps)}")
            view_history = confirm_with_context(
                self.processor,
                "📈 Review plots from previous steps?",
                default=False,
                module="AMBER Workflow Manager",
                description="Review plots from previous steps",
            )
            
            if view_history:
                self._show_historical_plots()
        
        # Check if the last step in checkpoint was actually completed
        last_step = checkpoint['step_name']
        last_step_output = f"{last_step}.out"
        
        if last_step in self.workflow_core.completed_steps:
            # Last step was completed, move to next step
            self.console.print(f"✅ {last_step} was completed, proceeding to next step")
            resume_step_index = checkpoint['step_index'] + 1
        else:
            # Last step was not completed, ask user what to do
            self.console.print(f"⚠️  {last_step} was started but may not have completed")
            
            if os.path.exists(last_step_output):
                # Check file size to estimate if it ran
                file_size = os.path.getsize(last_step_output)
                self.console.print(f"   📄 Output file exists ({file_size} bytes)")
                
                if file_size > 1000:  # Reasonable threshold
                    self.console.print("   💡 Appears simulation ran, but may have been interrupted")
                    action = prompt_with_context(
                        self.processor,
                        "   🤔 [r]estart this step, [c]ontinue to next, or [a]nalyze output?",
                        choices=["r", "c", "a"],
                        default="r",
                        module="AMBER Workflow Manager",
                        description="Resume action for partially complete step",
                        options_map={"r": "Restart this step", "c": "Continue to next", "a": "Analyze output"},
                    )
                    
                    if action == "a":
                        # Let user analyze the potentially incomplete output
                        self.workflow_core.analyze_single_step(last_step)
                        action = prompt_with_context(
                            self.processor,
                            "   🤔 After analysis, [r]estart this step or [c]ontinue to next?",
                            choices=["r", "c"],
                            default="r",
                            module="AMBER Workflow Manager",
                            description="Post-analysis resume action",
                            options_map={"r": "Restart this step", "c": "Continue to next"},
                        )
                    
                    if action == "c":
                        # Add to completed steps and move on
                        self.workflow_core.completed_steps.append(last_step)
                        resume_step_index = checkpoint['step_index'] + 1
                        self.console.print(f"   ✅ Marking {last_step} as completed, proceeding to next step")
                    else:
                        # Restart the step
                        resume_step_index = checkpoint['step_index']
                        self.console.print(f"   🔄 Restarting {last_step}")
                else:
                    self.console.print("   ⚠️  Output file is very small, likely simulation failed")
                    resume_step_index = checkpoint['step_index']
                    self.console.print(f"   🔄 Restarting {last_step}")
            else:
                self.console.print("   ❌ No output file found")
                resume_step_index = checkpoint['step_index']
                self.console.print(f"   🔄 Restarting {last_step}")
        
        # Continue from determined step
        return self._run_workflow_from_step(resume_step_index)

    def _run_workflow_from_step(self, start_step_index: int) -> bool:
        """Run workflow starting from a specific step index."""
        remaining_steps = self.workflow_core.workflow_steps[start_step_index:]
        
        self.console.print(f"\n🎮 Resuming Interactive Workflow")
        self.console.print(f"Remaining steps: {', '.join(remaining_steps)}")
        self.console.print("=" * 50)
        
        for i, step in enumerate(remaining_steps):
            step_index = start_step_index + i
            if not self._execute_step_interactive(step, step_index):
                return False
        
        self.console.print("\n🎉 Workflow completed successfully!")
        self._offer_final_analysis()
        return True

    def _monitor_current_simulation(self) -> bool:
        """Monitor currently running simulation."""
        # Look for running simulation by checking for .info files
        info_files = list(Path(".").glob("*.info"))
        
        if not info_files:
            self.console.print("❌ No active simulation found (no .info files)")
            return False
        
        # Find the most recently modified .info file
        latest_info = max(info_files, key=lambda f: f.stat().st_mtime)
        output_file = str(latest_info).replace('.info', '.out')
        
        if not os.path.exists(output_file):
            self.console.print(f"❌ Output file {output_file} not found")
            return False
        
        step_name = latest_info.stem
        self.console.print(f"📊 Monitoring {step_name} simulation")
        
        # Create monitor
        monitor = AMBERMonitor(output_file)
        
        # Simple monitoring loop
        self.console.print("💡 Commands: plot <metric>, timing, status, quit")
        
        while True:
            try:
                user_input = prompt_with_context(
                    self.processor, f"[{step_name}] Monitor>",
                    module="MD Manager - Monitor",
                    description="Interactive monitor command",
                ).strip().lower()
                
                if user_input == 'quit' or user_input == 'exit':
                    break
                elif user_input == 'status':
                    monitor.parse_amber_output()
                    if 'step' in monitor.data and len(monitor.data['step']) > 0:
                        last_step = monitor.data['step'][-1]
                        last_time = monitor.data['time'][-1] if 'time' in monitor.data else 0
                        self.console.print(f"📈 Last step: {last_step}, Time: {last_time:.1f} ps")
                    else:
                        self.console.print("⏳ No data available yet")
                elif user_input == 'timing':
                    self.console.print(monitor.get_timing_info())
                elif user_input.startswith('plot '):
                    metric = user_input.split(' ', 1)[1]
                    self._show_monitor_plot(monitor, metric)
                else:
                    self.console.print("❓ Available commands: plot <metric>, timing, status, quit")
                    
            except KeyboardInterrupt:
                break
        
        return True

    def _show_monitor_plot(self, monitor: AMBERMonitor, metric: str):
        """Show plot for monitoring."""
        monitor.parse_amber_output()
        
        plot_map = {
            'temp': ('temperature', 'Temperature', '(K)'),
            'energy': ('total_energy', 'Total Energy', '(kcal/mol)'),
            'pressure': ('pressure', 'Pressure', '(bar)'),
            'rms': ('rms_gradient', 'RMS Gradient', '(kcal/mol/Å)'),
            'kinetic': ('kinetic_energy', 'Kinetic Energy', '(kcal/mol)'),
            'potential': ('potential_energy', 'Potential Energy', '(kcal/mol)')
        }
        
        if metric in plot_map:
            data_key, title, units = plot_map[metric]
            self.console.print(monitor.ascii_plot(data_key, title, units))
        else:
            available = ', '.join(plot_map.keys())
            self.console.print(f"❌ Unknown metric. Available: {available}")

    def _analyze_historical_data(self) -> bool:
        """Analyze historical simulation data."""
        # Find completed output files
        output_files = list(Path(".").glob("*.out"))
        
        if not output_files:
            self.console.print("❌ No output files found for analysis")
            return False
        
        # Filter to workflow step files
        workflow_outputs = []
        for output_file in output_files:
            stem = output_file.stem
            if any(step in stem for step in self.workflow_core.workflow_steps):
                workflow_outputs.append(str(output_file))
        
        if not workflow_outputs:
            self.console.print("❌ No workflow output files found")
            return False
        
        self.console.print("📊 Available output files for analysis:")
        
        table = Table()
        table.add_column("Option", style="cyan")
        table.add_column("File", style="green")
        table.add_column("Size", style="white")
        
        for i, output_file in enumerate(workflow_outputs, 1):
            file_size = Path(output_file).stat().st_size
            size_str = f"{file_size:,} bytes"
            table.add_row(str(i), output_file, size_str)
        
        table.add_row(str(len(workflow_outputs) + 1), "Analyze all files", "")
        
        self.console.print(table)
        
        choices = [str(i) for i in range(1, len(workflow_outputs) + 2)]
        file_options = {str(i + 1): f for i, f in enumerate(workflow_outputs)}
        file_options[str(len(workflow_outputs) + 1)] = "Analyze all files"
        choice = prompt_with_context(
            self.processor,
            f"Select file to analyze [1-{len(workflow_outputs) + 1}]",
            choices=choices,
            default="1",
            module="AMBER Workflow Manager",
            description="Select workflow output file to analyze",
            options_map=file_options,
        )
        
        # Remap a recorded single-file pick by basename; the 'all' option carries
        # no basename and passes through to the len+1 branch.
        choice = remap_recorded_index(self.processor, workflow_outputs, str(choice))
        if choice == str(len(workflow_outputs) + 1):
            # Analyze all files
            for output_file in workflow_outputs:
                step_name = Path(output_file).stem
                self.console.print(f"\n📍 Analyzing {step_name}")
                self.console.print("-" * 40)
                self.workflow_core.analyze_single_step(step_name)
        else:
            # Analyze single file
            selected_file = workflow_outputs[int(choice) - 1]
            annotate_selected_path(self.processor, selected_file)
            step_name = Path(selected_file).stem
            self.workflow_core.analyze_single_step(step_name)
        
        return True

    def _show_historical_plots(self):
        """Show plots for completed simulation steps."""
        if not self.workflow_core.completed_steps:
            self.console.print("📊 No completed steps to analyze yet.")
            return
        
        self.console.print(f"\n📈 Historical Plot Viewer")
        self.console.print("=" * 50)
        self.console.print("Available completed steps:")
        
        for i, step in enumerate(self.workflow_core.completed_steps, 1):
            output_file = f"{step}.out"
            if os.path.exists(output_file):
                self.console.print(f"  {i}. {step.replace('_', ' ').title()}")
            else:
                self.console.print(f"  {i}. {step.replace('_', ' ').title()} (output file missing)")
        
        self.console.print(f"  {len(self.workflow_core.completed_steps) + 1}. View all steps")
        self.console.print(f"  {len(self.workflow_core.completed_steps) + 2}. Skip")
        
        choices = [str(i) for i in range(1, len(self.workflow_core.completed_steps) + 3)]
        step_options = {str(i + 1): s for i, s in enumerate(self.workflow_core.completed_steps)}
        step_options[str(len(self.workflow_core.completed_steps) + 1)] = "View all steps"
        step_options[str(len(self.workflow_core.completed_steps) + 2)] = "Skip"
        choice = prompt_with_context(
            self.processor,
            f"Select step to analyze [1-{len(self.workflow_core.completed_steps) + 2}]",
            choices=choices,
            default=str(len(self.workflow_core.completed_steps) + 2),
            module="AMBER Workflow Manager",
            description="Select completed step to re-analyze",
            options_map=step_options,
        )
        
        if choice == str(len(self.workflow_core.completed_steps) + 2):  # Skip
            return
        elif choice == str(len(self.workflow_core.completed_steps) + 1):  # View all
            self._analyze_complete_workflow()
        else:
            idx = int(choice) - 1
            if 0 <= idx < len(self.workflow_core.completed_steps):
                self.workflow_core.analyze_single_step(self.workflow_core.completed_steps[idx])

    def _analyze_complete_workflow(self):
        """Analyze the complete workflow."""
        self.console.print("\n📊 Complete Workflow Analysis")
        self.console.print("=" * 50)
        
        for step in self.workflow_core.completed_steps:
            output_file = f"{step}.out"
            if not os.path.exists(output_file):
                self.console.print(f"⚠️  Skipping {step} - output file not found")
                continue
            
            self.console.print(f"\n📍 {step.replace('_', ' ').title()}:")
            self.console.print("-" * 30)
            
            monitor = AMBERMonitor.create_historical_monitor(output_file, max_points=100)
            self.workflow_core.show_step_summary(monitor, step)

    def _configure_settings(self) -> bool:
        """Configure execution settings."""
        self.console.print("[bold]⚙️  Configure Workflow Settings[/bold]")
        
        # Current settings display
        current_table = Table(title="Current Settings")
        current_table.add_column("Setting", style="cyan")
        current_table.add_column("Value", style="white")
        
        current_table.add_row("Use MPI", str(self.workflow_core.use_mpi))
        current_table.add_row("MPI Cores", str(self.workflow_core.num_cores) if self.workflow_core.use_mpi else "N/A")
        current_table.add_row("GPU Device", str(self.workflow_core.gpu_device))
        
        self.console.print(current_table)
        
        # Allow changes
        change = confirm_with_context(
            self.processor,
            "\nModify settings?",
            default=False,
            module="AMBER Workflow Manager",
            description="Modify compute settings",
        )
        if not change:
            return True
        
        # MPI Configuration
        self.console.print("\n🖥️  CPU/MPI Configuration")
        self.workflow_core.use_mpi = confirm_with_context(
            self.processor,
            "Use MPI for CPU steps?",
            default=self.workflow_core.use_mpi,
            module="AMBER Workflow Manager",
            description="Use MPI for CPU steps (modify settings)",
        )
        
        if self.workflow_core.use_mpi:
            self.workflow_core.num_cores = int_prompt_with_context(
                self.processor,
                "Number of MPI cores",
                default=self.workflow_core.num_cores,
                module="AMBER Workflow Manager",
                description="Number of MPI cores (modify settings)",
            )
        
        # GPU Configuration
        self.console.print("\n🚀 GPU Configuration")
        self.workflow_core.gpu_device = int_prompt_with_context(
            self.processor,
            "GPU device ID",
            default=self.workflow_core.gpu_device,
            module="AMBER Workflow Manager",
            description="GPU device ID (modify settings)",
        )
        
        # Show updated settings
        self.console.print("\n✅ Updated Settings:")
        updated_table = Table()
        updated_table.add_column("Setting", style="cyan")
        updated_table.add_column("Value", style="green")
        
        updated_table.add_row("Use MPI", str(self.workflow_core.use_mpi))
        updated_table.add_row("MPI Cores", str(self.workflow_core.num_cores) if self.workflow_core.use_mpi else "N/A")
        updated_table.add_row("GPU Device", str(self.workflow_core.gpu_device))
        
        self.console.print(updated_table)
        
        return True

    def _manage_checkpoints(self) -> bool:
        """Manage checkpoint files."""
        self.console.print("[bold]💾 Checkpoint Management[/bold]")
        
        if not self.workflow_core.checkpoint_manager.checkpoint_exists():
            self.console.print("❌ No checkpoint file found")
            return True
        
        checkpoint = self.workflow_core.checkpoint_manager.load_checkpoint()
        if not checkpoint:
            self.console.print("❌ Failed to load checkpoint")
            return False
        
        # Display checkpoint info
        info_table = Table(title="Current Checkpoint")
        info_table.add_column("Property", style="cyan")
        info_table.add_column("Value", style="white")
        
        info_table.add_row("Timestamp", checkpoint['timestamp'])
        info_table.add_row("Current Step", checkpoint['step_name'])
        info_table.add_row("Completed Steps", ', '.join(checkpoint.get('completed_steps', [])))
        info_table.add_row("Mode", str(checkpoint['mode']))
        
        self.console.print(info_table)
        
        # Management options
        self.console.print("\nOptions:")
        self.console.print("1. View checkpoint details", highlight=False)
        self.console.print("2. Delete checkpoint", highlight=False)
        self.console.print("3. Export checkpoint", highlight=False)
        self.console.print("4. Return to menu", highlight=False)
        
        choice = prompt_with_context(
            self.processor,
            "Choose option [1-4]",
            choices=["1", "2", "3", "4"],
            default="4",
            module="AMBER Workflow Manager",
            description="Checkpoint management action",
        )
        
        if choice == "1":
            # Show detailed checkpoint info
            self.console.print("\n📋 Detailed Checkpoint Information:")
            for key, value in checkpoint.items():
                self.console.print(f"  {key}: {value}")
        
        elif choice == "2":
            # Delete checkpoint
            confirm = confirm_with_context(
                self.processor,
                "⚠️  Really delete checkpoint?",
                default=False,
                module="AMBER Workflow Manager",
                description="Confirm delete checkpoint",
            )
            if confirm:
                self.workflow_core.checkpoint_manager.delete_checkpoint()
                self.console.print("✅ Checkpoint deleted")
        
        elif choice == "3":
            # Export checkpoint
            export_name = prompt_with_context(
                self.processor,
                "Export filename",
                default="exported_checkpoint.json",
                module="AMBER Workflow Manager",
                description="Export filename for checkpoint",
            )
            try:
                import shutil
                shutil.copy2(self.workflow_core.checkpoint_manager.checkpoint_file, export_name)
                self.console.print(f"✅ Checkpoint exported to {export_name}")
            except Exception as e:
                self.console.print(f"❌ Export failed: {e}")
        
        return True

    def _offer_final_analysis(self):
        """Offer final workflow analysis."""
        analyze = confirm_with_context(
            self.processor,
            "📊 View complete workflow analysis?",
            default=True,
            module="AMBER Workflow Manager",
            description="View complete workflow analysis",
        )
        if analyze:
            self._analyze_complete_workflow()

    def _show_help(self) -> bool:
        """Show comprehensive help for workflow manager."""
        help_text = """
[bold]🚀 AMBER Workflow Manager Help[/bold]

The Workflow Manager orchestrates complete AMBER MD simulations with:

[bold]🎮 Execution Modes:[/bold]
• Interactive: Step-by-step with real-time monitoring
• Continuous: Automated execution of all steps

[bold]📊 Real-time Monitoring:[/bold]
• ASCII plots of energy, temperature, pressure
• Performance metrics and ETA calculations
• Live simulation status updates

[bold]💾 Checkpoint System:[/bold]
• Automatic checkpointing between steps
• Resume interrupted workflows
• Manage multiple checkpoint states

[bold]📈 Analysis Tools:[/bold]
• Historical data analysis
• Multi-step workflow comparison
• Performance optimization insights

[bold]⚙️ Configuration:[/bold]
• CPU/GPU execution settings
• MPI core allocation
• Custom simulation parameters

[bold]🔗 Integration:[/bold]
• Works with AMBER Input Generator presets
• Automatic file detection and validation
• Export results for further analysis

[bold]Available Commands in Monitoring:[/bold]
• plot temp/energy/pressure/rms - Show real-time plots
• timing - Show performance metrics
• status - Check simulation status
• stop - Terminate simulation
• continue - Proceed to next step

[bold]Original amber_workflow.py Features:[/bold]
All functionality from the original standalone script is preserved:
• Complete interactive workflow management
• Real-time ASCII plotting and monitoring
• Checkpoint/resume system
• Performance tracking with ETA calculations
• Educational command explanations
• Background process management
"""
        
        self.console.print(Panel(help_text, title="AMBER Workflow Manager"))
        return True